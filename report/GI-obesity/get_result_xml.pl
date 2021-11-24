#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(cwd);
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use POSIX;
use Template;
use Data::Dumper;
use File::Slurp;
use Scalar::Util 'refaddr';
use XML::Simple qw(XMLin XMLout);
use Data::Printer colored => 1;
use utf8;
use Encode;

use Statistics::R;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my $version = "1.0.0";
my $R = Statistics::R->new();

###
my ($datadir,$indir,$referdir,$dirout,$language,$barcode,$historydir,$confdir,$level_set,$help);
GetOptions (
	"data=s"	=> \$datadir,
	"in=s"		=> \$indir,
	"referdir=s"	=> \$referdir,
	"out=s"		=> \$dirout,
	"language=s"	=> \$language,
	"barcode=s"	=> \$barcode,
	"historydir=s"	=> \$historydir,
	"confdir=s"	=> \$confdir,
	"l=s"		=> \$level_set,
	"h"	=> \$help,
) or &USAGE;
&USAGE if ((!defined $datadir || !defined $barcode) || defined $help) ;

########## AbsolutePath and default ##############
$dirout||="./";
$indir||="$Bin/db2xml/current";
$referdir||="$Bin/../../reference_data/current";
$confdir||="$Bin/../conf/current";
$historydir||="/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
mkdir $dirout if (! -d $dirout);
$indir=&AbsolutePath("dir",$indir);
$datadir=&AbsolutePath("dir",$datadir);
$referdir=&AbsolutePath("dir",$referdir);
$dirout=&AbsolutePath("dir",$dirout);
$historydir=&AbsolutePath("dir",$historydir);
$confdir=&AbsolutePath("dir",$confdir);

$language||="CN";

$level_set||="N";

################# main process ###################
my %vars=();
## check file and dir
&check_file($datadir,$indir,$referdir,$confdir,$barcode);

## load sample information file
my ($samInfo)=&load_sample_info($barcode,$datadir,\%vars);

## load $type.list.order$sex
my ($listH,$term2oriid,$balancecomposition,$compositionweigth)=&load_list($indir,$samInfo,$confdir);

## calculate meta level and per by sample data and ref data
my ($Data)=&calculate_meta_level_per($datadir,$indir,$referdir,$barcode,$term2oriid,$compositionweigth,$level_set);
&calculate_balance_level_per($Data,$balancecomposition,$compositionweigth,$level_set);

## main process
#&main_process($Data,$listH,\%{$diseaseOrder{question}},\%{$diseaseOrder{predict}},\%vars,$datadir,$indir,$language,$barcode,$samInfo);
&main_process($Data,$listH,\%vars,$datadir,$indir,$language,$barcode,$samInfo,$historydir);

## sample result csv for erp import
&sample_result_csv(\%vars,$listH->{pathogens},$barcode,$dirout);

############################# sub ###############################################

#&main_process($Data,$listH,\%vars,$datadir,$indir,$language,$barcode,$samInfo);
sub main_process {
	my ($Data,$listH,$vars,$datadir,$indir,$language,$barcode,$samInfo,$historydir)=@_;

	## 菌群分布(菌属分布 致病菌)
	# $listH->{"distribution"}; $listH->{"pathogens"}
	&Intestinal_Flora($Data,$listH->{distribution},$vars,"$indir/XML/distribution/$language",$language);
	&Intestinal_Pathogen($Data,$listH->{pathogens},$vars,"$indir/XML/pathogens/$language",$language);

	## 肠道微生态平衡能力(多样性、有益菌、有害菌)
	# $listH->{"overview"}
	&Micro_Overvwiew($Data,$listH->{overview},$vars,"$indir/XML/overview/$language",$language);

	## 菌群总览指标(平衡能力、胖菌、瘦菌)
	# $listH->{"metasummary"}
	&Micro_Metasummary($Data,$listH->{metasummary},$vars,"$indir/XML/metasummary/$language",$language);

	## 疾病风险(T2D)
	# $listH->{"disease"};
	&Intestinal_Diseaserisk($Data,$listH->{disease},$vars,"$indir/XML/disease/$language",$language);

	## 肠道调节方案
	# $listH->{"advise"};
	&Program($vars,"$indir/XML/advise",$language,$samInfo,$listH->{advise});

	## 多次检测结果比较
	# $samInfo->{batch}; $samInfo->{'历次检测编号'}
	#&Compare_history($vars,$samInfo,$historydir) if ($samInfo->{batch} > 1);

	## 报告日期
	#$vars->{reportdate}=strftime("%Y-%m-%d",localtime());
	
	## output xml
	&output_xml($vars,$barcode,$dirout);
	
	return;
}

#&Intestinal_Diseaserisk($Data,$listH->{disease},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Diseaserisk {
	my ($Data,$list1,$vars,$dir1,$language)=@_;
	
	##
	$vars->{diseaserisk}{name}="潜在疾病风险";
	$vars->{diseaserisk}{ename}="DISEASERISK";

	##
	&item_Intestinal_Diseaserisk('desc',$Data,$list1,$vars,$dir1,$language);

	return;
}

#&item_Intestinal_Diseaserisk('fenbu',$Data,$list,$vars,$dir,$language);
sub item_Intestinal_Diseaserisk {
	my ($prefix,$Data,$list,$vars,$dir,$language)=@_;
	
	##
	my (@good,@bad);
	##
	$vars->{diseaserisk}{$prefix}{abnormalN}=0;
	$vars->{diseaserisk}{$prefix}{goodN}=0;
	$vars->{diseaserisk}{$prefix}{badN}=0;
	##

	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		push @{$vars->{diseaserisk}{$prefix}{item}},$item;


		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{diseaserisk}{$prefix}{abnormalN}+=$add;
		$vars->{diseaserisk}{$prefix}{goodN}+=$addG;
		$vars->{diseaserisk}{$prefix}{badN}+=$addN;
	}

	## good,bad
	($vars->{diseaserisk}{$prefix}{good},$vars->{diseaserisk}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	
}

#&Compare_history($vars,$samInfo,$historydir);
sub Compare_history {
	my ($vars,$samInfo,$historydir) = @_;

	my $report_version=$samInfo->{version};
	my @historys=split/,/,$samInfo->{'历次检测编号'};
	my $batchlast=pop(@historys);
	&extract_compare_info($vars,$batchlast,$historydir,$report_version);

	if (@historys) {
		for (1..@historys) {
			my $batch_str = "batch"."$_";
			my $xmlfile="$historydir/$report_version/$historys[$_-1].xml";
			if (!-f $xmlfile) {
				print STDERR "$historys[$_] $historydir/$report_version/$historys[$_].xml result unexists, please Check!\n\n";
				die;
			}
			my $xml=&read_xml($xmlfile);
			$vars->{$batch_str}{gaikuang}{desc}=\%{$xml->{gaikuang}{desc}};
			$vars->{$batch_str}{changdaojunqun}{fenbu}{goodN}=$xml->{changdaojunqun}{fenbu}{goodN};
			$vars->{$batch_str}{changdaojunqun}{fenbu}{badN}=$xml->{changdaojunqun}{fenbu}{badN};
			$vars->{$batch_str}{changdaojunqun}{fenbu}{good}=$xml->{changdaojunqun}{fenbu}{good};
			$vars->{$batch_str}{changdaojunqun}{fenbu}{bad}=$xml->{changdaojunqun}{fenbu}{bad};
		}
	}
}

#&extract_compare_info($vars,$batchlast,$historydir);
sub extract_compare_info {
	my ($vars,$batchlast,$historydir,$report_version) = @_;

	my $xmlfile="$historydir/$report_version/$batchlast.xml";
	chomp $xmlfile;
	if (!-f $xmlfile) {
		print STDERR "$batchlast $historydir/$report_version/$batchlast.xml result unexists, please Check!\n\n";
		die;
	}
	my $xml=&read_xml($xmlfile);

	my $batch_str = "batchlast";

	$vars->{$batch_str}{gaikuang}{desc}=\%{$xml->{gaikuang}{desc}};
	&compare_gaikuang($vars, "youyijun", $vars->{gaikuang}{desc}{youyijun}{per100}, $xml->{gaikuang}{desc}{youyijun}{per100});

	$vars->{$batch_str}{changdaojunqun}{fenbu}{goodN}=$xml->{changdaojunqun}{fenbu}{goodN};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{badN}=$xml->{changdaojunqun}{fenbu}{badN};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{good}=$xml->{changdaojunqun}{fenbu}{good};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{bad}=$xml->{changdaojunqun}{fenbu}{bad};

}

#&compare_gaikuang($vars, "duoyangxing", $vars->{gaikuang}{desc}{duoyangxing}{per}, $xml->{gaikuang}{desc}{duoyangxing}{per});
sub compare_gaikuang {
	my ($vars, $key, $valnow, $vallast) = @_;

	if (($valnow - $vallast) >= 5) {
		$vars->{comparelast}{$key} = "升高";
	}
	elsif (($valnow - $vallast) <= -5) {
		$vars->{comparelast}{$key} = "降低";
	}
	else {
		$vars->{comparelast}{$key} = "保持不变";
	}
}

#&sample_result_csv(\%vars,$listH,$barcode,$dirout);
sub sample_result_csv {
	my ($vars,$list,$barcode,$dirout) = @_;


	my $outfile="$dirout/$barcode".".result.csv";
	my $outfile1="$dirout/$barcode".".result_stat.csv";

	open my $out, '>:encoding(utf-8)', $outfile || die $!;
	open my $out1, '>:encoding(utf-8)', $outfile1 || die $!;

	my $stat_head = "检测编号\t姓名\t";
	my $stat_result = "$vars->{barcode}\t$vars->{name}\t";
	my $stat_bad = 0;
	my $stat_good = 0;
	my $bad_per = 0;
	my $good_per = 0;

	print $out "检测者信息\n检测编号\t姓名\n$vars->{barcode}\t$vars->{name}\n";

#	my $head_str = "";
	my $term_str = "";
	my $result_str = "";
	my $value_str = "";
#	my $range_str = "\t";
	my $per_str = "";

	my $level;

	print $out "检测结果总体评价\n胖菌\t瘦菌\n";
#	$level = &transfer_item_level($vars->{metasummary}{desc}{pinghengnengli}{level});
#	$result_str ="$level($vars->{metasummary}{desc}{pinghengnengli}{risk})\t";
#	$value_str .= "\t";
#	$range_str .= "\t$vars->{metasummary}{desc}{shoujun}{range}";
#	$per_str .= $vars->{metasummary}{desc}{pinghengnengli}{per}*100;
#	$per_str .= "(%)\t";

$stat_head = "检测编号\t姓名\t";                   
$stat_result = "$vars->{barcode}\t$vars->{name}\t";

	$level = &transfer_item_level($vars->{metasummary}{desc}{pangjun}{level});
	$result_str .="$level($vars->{metasummary}{desc}{pangjun}{risk})\t";
	$value_str .= "$vars->{metasummary}{desc}{pangjun}{value}\t";
#	$range_str .= "\t$vars->{metasummary}{desc}{shoujun}{range}";
	$per_str .= $vars->{metasummary}{desc}{pangjun}{per}*100;
	$per_str .= "(%)\t";

	$stat_head .= "胖菌BAD\t胖菌BAD(%)\t胖菌GOOD\t胖菌GOOD(%)\t";
	$stat_bad = $vars->{changdaojunqun}{pangjun}{badN} if (defined $vars->{changdaojunqun}{pangjun}{badN});
	$bad_per = sprintf "%.2f", $stat_bad/@{$vars{changdaojunqun}{pangjun}{item}}*100;
	$stat_good = $vars->{changdaojunqun}{pangjun}{goodN} if (defined $vars->{changdaojunqun}{pangjun}{goodN});
	$good_per = sprintf "%.2f", $stat_good/@{$vars{changdaojunqun}{pangjun}{item}}*100;
	$stat_result .= "$stat_bad\t$bad_per\t$stat_good\t$good_per\t";


	$level = &transfer_item_level($vars->{metasummary}{desc}{shoujun}{level});
	$result_str .="$level($vars->{metasummary}{desc}{shoujun}{risk})";
	$value_str .= "$vars->{metasummary}{desc}{shoujun}{value}";
#	$range_str .= "\t$vars->{metasummary}{desc}{shoujun}{range}";
	$per_str .= $vars->{metasummary}{desc}{shoujun}{per}*100;
	$per_str .= "(%)";
	$result_str =~s/\s+$//;

	$stat_head .= "瘦菌BAD\t瘦菌BAD(%)\t瘦菌GOOD\t瘦菌GOOD(%)\t";
	$stat_bad = $vars->{changdaojunqun}{shoujun}{badN} if (defined $vars->{changdaojunqun}{shoujun}{badN});
	$bad_per = sprintf "%.2f", $stat_bad/@{$vars{changdaojunqun}{shoujun}{item}}*100;
	$stat_good = $vars->{changdaojunqun}{shoujun}{goodN} if (defined $vars->{changdaojunqun}{shoujun}{goodN});
	$good_per = sprintf "%.2f", $stat_good/@{$vars{changdaojunqun}{shoujun}{item}}*100;
	$stat_result .= "$stat_bad\t$bad_per\t$stat_good\t$good_per\t";

	print $out "$result_str\n";
	print $out "$value_str\n$per_str\n";

#	$head_str = "";
	$term_str = "";
	$result_str = "";
	$value_str = "";
#	$range_str = "\t";
	$per_str = "";

	print $out "微生态平衡能力\n肠道菌群多样性\t有益菌\t有害菌\n";
	$level = &transfer_item_level($vars->{gaikuang}{desc}{duoyangxing}{level});
	$result_str ="$level($vars->{gaikuang}{desc}{duoyangxing}{risk})\t";
	$value_str .= "$vars->{gaikuang}{desc}{duoyangxing}{value}\t";
#	$range_str .= "\t$vars->{gaikuang}{desc}{youhaijun}{range}";
	$per_str .= $vars->{gaikuang}{desc}{duoyangxing}{per}*100;
	$per_str .= "(%)\t";

	$stat_head .= "肠道菌群多样性值\t肠道菌群多样性百分比\t肠道菌群多样性风险";
	my $xx = sprintf "%.2f", $vars->{gaikuang}{desc}{duoyangxing}{per}*100;
	$stat_result .= "$vars->{gaikuang}{desc}{duoyangxing}{value}\t$xx\t$level($vars->{gaikuang}{desc}{duoyangxing}{risk})";

	print $out1 "$stat_head\n$stat_result\n";
	close $out1;

	$level = &transfer_item_level($vars->{gaikuang}{desc}{youyijun}{level});
	$result_str .="$level($vars->{gaikuang}{desc}{youyijun}{risk})\t";
	$value_str .= "$vars->{gaikuang}{desc}{youyijun}{value}\t";
#	$range_str .= "\t$vars->{gaikuang}{desc}{youhaijun}{range}";
	$per_str .= $vars->{gaikuang}{desc}{youyijun}{per}*100;
	$per_str .= "(%)\t";

	$level = &transfer_item_level($vars->{gaikuang}{desc}{youhaijun}{level});
	$result_str .="$level($vars->{gaikuang}{desc}{youhaijun}{risk})";
	$value_str .= "$vars->{gaikuang}{desc}{youhaijun}{value}";
#	$range_str .= "\t$vars->{gaikuang}{desc}{youhaijun}{range}";
	$per_str .= $vars->{gaikuang}{desc}{youhaijun}{per}*100;
	$per_str .= "(%)";
	$result_str =~s/\s+$//;
	print $out "$result_str\n";
	print $out "$value_str\n$per_str\n";


#	$head_str = "";
	$term_str = "";
	$result_str = "";
	$value_str = "";
#	$range_str = "\t";
	$per_str = "";

#	$head_str .= "\t菌群分布";
	print $out "胖菌菌属分布\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{changdaojunqun}{pangjun}{item}}) {
#		$head_str .= "\t";
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
		$value_str .= "$item->{value}\t";
#		$range_str .= "\t$item->{range}";
		$per_str .= $item->{per}*100;
		$per_str .= "(%)\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";
#	print $out "$value_str\n$per_str\n";

#	$head_str = "";
	$term_str = "";
	$result_str = "";
	$value_str = "";
#	$range_str = "\t";
	$per_str = "";

#	$head_str .= "\t菌群分布";
	print $out "瘦菌菌属分布\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{changdaojunqun}{shoujun}{item}}) {
#		$head_str .= "\t";
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
		$value_str .= "$item->{value}\t";
#		$range_str .= "\t$item->{range}";
		$per_str .= $item->{per}*100;
		$per_str .= "(%)\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";
#	print $out "$value_str\n$per_str\n";

#	$head_str .= "致病菌含量";
	print $out "致病菌\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{changdaojunqun}{pathogen}{item}}) {
#		$head_str .= "\t";
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
#		$value_str .= "\t$item->{value}";
#		$range_str .= "\t$item->{range}";
#		$per_str .= "\t$item->{per}";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	close $out;

	return;
}

#&Program($vars,"$indir/XML",$language,$samInfo,$listH->{diseaseList});
sub Program {
	my ($vars,$dir,$language,$samInfo,$list)=@_;
	#my ($foodKey,$regulateKey,$sportKey)=('-','-','-');
	my $ageKey=&get_program_age($samInfo->{age});

	## process food and others
	&process_food($vars,$dir,$ageKey,$samInfo->{sex},$language);
	## process regulate
	&process_regulate($vars,$dir,$ageKey,$samInfo->{sex},$language);
	
	return;
}

#&process_regulate($vars,$dir,$regulateKey,$language);
sub process_regulate {
	my ($vars,$dir,$ageKey,$sex,$language)=@_;
	
	## 补充益生元和益生菌
	my $vd=XMLin("$dir/$language/tzgltiaojiefangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{regulate}{name} = $vd->{$language}{title};
	$vars->{fangan}{regulate}{desc} = '-';
	my $bmicondition = &get_bmicondition($vars->{BMI});
	my @actdescs = split/\n/,$vd->{$language}{bmisuggestion}{$bmicondition}{$ageKey}{action};
	my $note = $vd->{$language}{bmisuggestion}{$bmicondition}{$ageKey}{note};
	foreach my $actdesc (@actdescs) {
		$actdesc =~ /(^\d+):(.*)/;
		if ($vars->{metasummary}{desc}{badN} == $1) {
			$vars->{fangan}{regulate}{desc} = $2;
			$vars->{fangan}{regulate}{desc} =~s/\[|\]//g;
			$vars->{fangan}{regulate}{note} = $note;
		}
	}

	return;
}

sub get_bmicondition {
	my ($bmi) = @_;
	my $bmikey = "bmi1";

	if ($bmi >= 32 or ($bmi >= 22 && $bmi < 24)) {
		$bmikey = "bmi1";
	}
	elsif ($bmi >= 30 && $bmi < 32) {
		$bmikey = "bmi2";
	}
	elsif ($bmi >= 28 && $bmi < 30) {
		$bmikey = "bmi3";
	}
	elsif ($bmi >= 24 && $bmi < 28) {
		$bmikey = "bmi4";
	}
	elsif ($bmi < 22) {
		$bmikey = "bmi5";
	}
}

#&process_food_sport($vars,$dir,$diseaseF,$ageF,$badDiseaseF,$language);
sub process_food {
	my ($vars,$dir,$ageKey,$sex,$language)=@_;
	
	## food
	if ($vars->{sex} eq "M" && ($vars->{BMI} >= 24 && $vars->{BMI} < 28)) {
		$vars->{fangan}{food}{desc} = "gaodanbaijianzhongyinshi.pdf";
	}
	else {
		if ($vars->{BMI} < 24) {
			$vars->{fangan}{food}{desc} = "dizhonghaiyinshi.pdf";
		}
		elsif ($vars->{BMI} >= 24 && $vars->{BMI} < 32) {
			$vars->{fangan}{food}{desc} = "ditanshuihuahewujianzhongyinshi.pdf";
		}
		else {
			$vars->{fangan}{food}{desc} = "chaodikalulijianzhongyinshi.pdf";
		}
	}

	$vars->{fangan}{sport}{desc} = "ranzhiyundongfangan.pdf";

#	my $vd_food=XMLin("$dir/$language/qpcryinshijianyi.xml",NoAttr=>1,SuppressEmpty => "");
#	$vars->{fangan}{food}{name} = $vd_food->{$language}{title};
#	#$vars->{fangan}{food}{desc} = $vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	my @food_notes = split/\n/,$vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	$vars->{fangan}{food}{head} = shift @food_notes;
#	my $food_note_index = 0;
#	foreach my $note (@food_notes) {
#		$note=~s/^\s+|\s+$//;
#		next if ($note=~/^$/);
#		$food_note_index ++;
#		my ($index_name,$index_desc) = split/：/,$note;
#		#$vars->{fangan}{others}{"desc$note_index"} = $note;
#		push @{$vars->{fangan}{food}{desc}{item}}, {name => $index_name, desc => $index_desc, num => $food_note_index};
#	}
#
#
#	my $vd_others=XMLin("$dir/$language/qpcrrichangdiandi.xml",NoAttr=>1,SuppressEmpty => "");
#	$vars->{fangan}{others}{name} = $vd_others->{$language}{title};
#	my @others_notes = split/\n/,$vd_others->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	$vars->{fangan}{others}{head} = shift @others_notes;
#	my $others_note_index = 0;
#	foreach my $note (@others_notes) {
#		$note=~s/^\s+|\s+$//;
#		next if ($note=~/^$/);
#		$others_note_index ++;
#		#$vars->{fangan}{others}{"desc$note_index"} = $note;
#		push @{$vars->{fangan}{others}{desc}{item}}, {desc => $note, num => $others_note_index};
#	}

	return;
}

#($age)=&get_program_age($samInfo->{age});
sub get_program_age {
	my ($num)=@_;
	my $age="age5";
	
	##
	if ($num <= 1) {
		$age="age0";
	}
	elsif ($num < 3) {
		$age="age1";
	}
	elsif ($num < 7) {
		$age="age2";
	}
	elsif ($num < 13) {
		$age="age3";
	}
	elsif ($num < 19) {
		$age="age4";
	}
	elsif ($num < 45) {
		$age="age5";
	}
	elsif ($num < 65) {
		$age="age6";
	}
	else {
		$age="age7";
	}
	
	return ($age);
}

#($vars->{changdaojunqun}{$prefix}{goodSource},$vars->{changdaojunqun}{$prefix}{badSource})=&process_good_badSource(\@goodSource,\@badSource);
sub process_good_badSource {
	my ($goodSource,$badSource)=@_;
	my ($good_source,$bad_source)=('-','-');

	##
	if (@{$goodSource}) {
		$good_source=join "",@{$goodSource};
	}

	if (@{$badSource}) {
		$bad_source=join "",@{$badSource};
	}

	return ($good_source,$bad_source);
}

#($vars->{changdaojunqun}{desc}{$prefix}{goodDesc},$vars->{changdaojunqun}{desc}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);
sub process_good_badDesc {
	my ($goodDesc,$badDesc)=@_;
	my ($good_desc,$bad_desc)=('-','-');

	##
	my ($good)=&uniq_array($goodDesc);
	my ($bad)=&uniq_array($badDesc);
	##
	if (@{$good}) {
		if (scalar @{$good} <= 5) {
			$good_desc=join "；",@{$good};
		}
		else {
			$good_desc=join("；",$good->[0],$good->[1],$good->[2],$good->[3],$good->[4]);
			$good_desc.="等";
		}
	}

	if (@{$bad}) {
		if (scalar @{$bad} <= 5) {
			$bad_desc=join "；",@{$bad};
		}
		else {
			$bad_desc=join("；",$bad->[0],$bad->[1],$bad->[2],$bad->[3],$bad->[4]);
			$bad_desc.="等";
		}
	}

	return ($good_desc,$bad_desc);
}

#my ($dis)=&uniq_array($Disease);
sub uniq_array {
	my ($array)=@_;
	my @dis=();
	
	##
	my $n=1;
	my %hash;
	foreach my $unit (@{$array}) {
		$hash{$unit}=$n++ if (! exists $hash{$unit});
	}
	##
	@dis=sort {$hash{$a} <=> $hash{$b}} keys %hash;
	
	return (\@dis);
}

#($vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{badEffect})=&process_goodbad_DiseaseEffect(\@goodDisease,\@badDisease,\@goodEffect,\@badEffect);
sub process_goodbad_DiseaseEffect {
	my ($good_Disease,$bad_Disease,$good_Effect,$bad_Effect)=@_;
	my ($goodDisease,$badDisease,$goodEffect,$badEffect)=('-','-','-','-');

	##
	my ($gooddis)=&split_uniq_array($good_Disease);
	my ($baddis)=&split_uniq_array($bad_Disease);
	my ($goodeff)=&split_uniq_array($good_Effect);
	my ($badeff)=&split_uniq_array($bad_Effect);
	##
	if (@{$gooddis}) {
		if (scalar @{$gooddis} <= 3) {
			$goodDisease=join "、",@{$gooddis};
		}
		else {
			$goodDisease=join("、",$gooddis->[0],$gooddis->[1],$gooddis->[2]);
			$goodDisease.="等";
		}
	}
	if (@{$baddis}) {
		if (scalar @{$baddis} <= 3) {
			$badDisease=join "、",@{$baddis};
		}
		else {
			$badDisease=join("、",$baddis->[0],$baddis->[1],$baddis->[2]);
			$badDisease.="等";
		}
	}
	if (@{$goodeff}) {
		if (scalar @{$goodeff} <= 3) {
			$goodEffect=join "、",@{$goodeff};
		}
		else {
			$goodEffect=join("、",$goodeff->[0],$goodeff->[1],$goodeff->[2]);
			$goodEffect.="等";
		}
	}
	if (@{$badeff}) {
		if (scalar @{$badeff} <= 3) {
			$badEffect=join "、",@{$badeff};
		}
		else {
			$badEffect=join("、",$badeff->[0],$badeff->[1],$badeff->[2]);
			$badEffect.="等";
		}
	}

	$goodDisease='-' if ($goodDisease eq "");
	$badDisease='-' if ($badDisease eq "");
	$goodEffect='-' if ($goodEffect eq "");
	$badEffect='-' if ($badEffect eq "");

	return ($goodDisease,$badDisease,$goodEffect,$badEffect);
}

#&split_uniq_hash($Disease,$badDiseaseTotal);
sub split_uniq_hash {
	my ($array,$hash)=@_;
	
	##
	foreach my $unit (@{$array}) {
		my @d=split /\；/,$unit;
		foreach my $dd (@d) {
			$hash->{$dd}=1;
		}
	}
	
	return;
}

#my ($dis)=&split_uniq_array($Disease);
sub split_uniq_array {
	my ($array)=@_;
	my @dis=();
	
	##
	my $n=1;
	my %hash;
	foreach my $unit (@{$array}) {
		my @d=split /\；/,$unit;
		foreach my $dd (@d) {
			$hash{$dd}=$n++ if (! exists $hash{$dd});
		}
	}
	##
	@dis=sort {$hash{$a} <=> $hash{$b}} keys %hash;
	
	return (\@dis);
}

#&Micro_Overvwiew($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Micro_Overvwiew {
	my ($Data,$list1,$vars,$dir1,$language)=@_;
	
	##
	$vars->{gaikuang}{name}="肠道菌群概况";
	$vars->{gaikuang}{ename}="OVERVIEW";

	##
	&item_Micro_Overvwiew('desc',$Data,$list1,$vars,$dir1,$language);

	return;
}

#$risk=&transfer_item_level($risk);
sub transfer_item_level {
	my ($risk)=@_;
	
	##
	if ($risk eq 'l1') {
		$risk="低";
	}
	elsif ($risk eq 'l2') {
		$risk="偏低";
	}
	elsif ($risk eq 'l3') {
		$risk="正常";
	}
	elsif ($risk eq 'l4') {
		$risk="偏高";
	}
	else {
		$risk="高";
	}
	
	return ($risk);
}

#$risk=&transfer_pathogen_level($risk);
sub transfer_pathogen_level {
	my ($risk)=@_;
	
	##
	if ($risk eq 'l0') {
		$risk="未检出";
	}
	elsif ($risk eq 'l1' || $risk eq 'l2' || $risk eq 'l3') {
		$risk="含量未超标";
	}
	elsif ($risk eq 'l4') {
		$risk="偏高";
	}
	else {
		$risk="高";
	}
	
	return ($risk);
}

#$level=&get_summary_level($risk,$protect);
sub get_summary_level {
	my ($risk,$protect)=@_;
	my $level='-';
	
	##
	if ($risk eq "正常" && $protect eq "正常") {
		$level='l1';
	}
	elsif ($risk eq "正常" && $protect eq "低") {
		$level='l2';
	}
	elsif ($risk eq "正常" && $protect eq "高") {
		$level='l3';
	}
	elsif ($risk eq "低" && $protect eq "正常") {
		$level='l4';
	}
	elsif ($risk eq "低" && $protect eq "低") {
		$level='l5';
	}
	elsif ($risk eq "低" && $protect eq "高") {
		$level='l6';
	}
	elsif ($risk eq "高" && $protect eq "正常") {
		$level='l7';
	}
	elsif ($risk eq "高" && $protect eq "低") {
		$level='l8';
	}
	elsif ($risk eq "高" && $protect eq "高") {
		$level='l9';
	}

	return ($level);
}

#&item_Micro_Overvwiew('fenbu',$Data,$list,$vars,$dir,$language);
sub item_Micro_Overvwiew {
	my ($prefix,$Data,$list,$vars,$dir,$language)=@_;
	
	##
	my (@good,@bad);
	##
	$vars->{gaikuang}{$prefix}{abnormalN}=0;
	$vars->{gaikuang}{$prefix}{goodN}=0;
	$vars->{gaikuang}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		$vars->{gaikuang}{$prefix}{duoyangxing}=$item if ($meta =~ /junqunduoyangxing/);
		if ($meta =~ /youyijun/) {
			$vars->{gaikuang}{$prefix}{youyijun}=$item;
		}
		if ($meta =~ /youhaijun/) {
			$vars->{gaikuang}{$prefix}{youhaijun}=$item;
		}

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{gaikuang}{$prefix}{abnormalN}+=$add;
		$vars->{gaikuang}{$prefix}{goodN}+=$addG;
		$vars->{gaikuang}{$prefix}{badN}+=$addN;
	}

	## good,bad
	($vars->{gaikuang}{$prefix}{good},$vars->{gaikuang}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	
}

#&Micro_Metasummary($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Micro_Metasummary {
	my ($Data,$list1,$vars,$dir1,$language)=@_;
	
	##
	$vars->{metasummary}{name}="肠道菌群总览";
	$vars->{metasummary}{ename}="METASUMMARY";

	##
	&item_Micro_Metasummary('desc',$Data,$list1,$vars,$dir1,$language);
	&get_Metasummary_pic('desc',$vars);

	return;
}

sub get_Metasummary_pic {
	my ($prefix, $vars) = @_;

	if ($vars->{BMI} >= 32) {
		$vars->{metasummary}{desc}{difficulty} = "difficulty-l5.pdf";
		$vars->{metasummary}{desc}{regulatetime} = "regulatetime-l4.pdf";
	}
	elsif ($vars->{BMI} >= 30 || $vars->{BMI} < 32) {
		$vars->{metasummary}{desc}{difficulty} = "difficulty-l4.pdf";
		$vars->{metasummary}{desc}{regulatetime} = "regulatetime-l3.pdf";
	}
#	elsif (($vars->{BMI} >= 22 && $vars->{BMI} < 24) || $vars->{BMI} < 22) {
	elsif ($vars->{BMI} < 24) {
		$vars->{metasummary}{desc}{difficulty} = "difficulty-l3.pdf";
		$vars->{metasummary}{desc}{regulatetime} = "regulatetime-l3.pdf";
	}
	elsif ($vars->{BMI} >= 28 || $vars->{BMI} < 30) {
		$vars->{metasummary}{desc}{difficulty} = "difficulty-l2.pdf";
		$vars->{metasummary}{desc}{regulatetime} = "regulatetime-l2.pdf";
	}
	else {
		$vars->{metasummary}{desc}{difficulty} = "difficulty-l1.pdf";
		$vars->{metasummary}{desc}{regulatetime} = "regulatetime-l1.pdf";
	}

	if ($vars->{BMI} < 18.5) {
		$vars->{metasummary}{desc}{weightstatus} = "偏瘦";
	}
	elsif ($vars->{BMI} >= 18.5 && $vars->{BMI} < 24) {
		$vars->{metasummary}{desc}{weightstatus} = "正常";
	}
	elsif ($vars->{BMI} >= 24 && $vars->{BMI} < 28) {
		$vars->{metasummary}{desc}{weightstatus} = "偏胖";
	}
	else {
		$vars->{metasummary}{desc}{weightstatus} = "肥胖";
	}

}

#&item_Micro_Metasummary('fenbu',$Data,$list,$vars,$dir,$language);
sub item_Micro_Metasummary {
	my ($prefix,$Data,$list,$vars,$dir,$language)=@_;
	
	##
	my (@good,@bad,@normal,@goodDesc,@badDesc,@normalDesc,@goodEffect,@badEffect,@normalEffect,@goodActrule,@badActrule,@normalActrule);
	##
	$vars->{metasummary}{$prefix}{abnormalN}=0;
	$vars->{metasummary}{$prefix}{goodN}=0;
	$vars->{metasummary}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");

		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		$vars->{metasummary}{$prefix}{pinghengnengli}=$item if ($meta =~ /weishengtaipinghengnengli/);
		$vars->{metasummary}{$prefix}{pangjun}=$item if ($meta =~ /zengpangjun/);
		$vars->{metasummary}{$prefix}{shoujun}=$item if ($meta =~ /shoushenjun/);

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{metasummary}{$prefix}{abnormalN}+=$add;
		$vars->{metasummary}{$prefix}{goodN}+=$addG;
		$vars->{metasummary}{$prefix}{badN}+=$addN;

		## get goodDesc,badDesc
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			### Get risk disease info
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badActrule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
		}
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @goodActrule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
		}
		if (defined $item->{risk} && $item->{risk} eq "正常") {
			push @normal, $item->{name};
			### Get risk disease info
			push @normalDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @normalEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @normalActrule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
		}
	}

	## good,bad
	($vars->{metasummary}{$prefix}{good},$vars->{metasummary}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	## badDesc,goodDesc
	($vars->{metasummary}{$prefix}{goodDesc},$vars->{metasummary}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);
	## badActrule,goodActrule
	($vars->{metasummary}{$prefix}{goodActrule},$vars->{metasummary}{$prefix}{badActrule})=&process_good_badDesc(\@goodActrule,\@badActrule);

	## process metasummary normal item
	($vars->{metasummary}{$prefix}{normalitem})=&process_normal_metasummary(\@normal);

}

sub process_normal_metasummary {
	my ($normal)=@_;
	my ($normalitem)=('-','-');
	
	##
	if (@{$normal}) {
		if (scalar @{$normal} <= 3) {
			$normalitem=join "、",@{$normal};
		}
		else {
			$normalitem=join("、",$normal->[0],$normal->[1],$normal->[2]);
			#$good.="等";
		}
	}
	return ($normalitem);
}

#&Intestinal_Flora($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Flora {
	my ($Data,$list,$vars,$dir,$language)=@_;
	
	##
	$vars->{changdaojunqun}{name}="肠道菌群";
	$vars->{changdaojunqun}{ename}="MICROBIOTA";

	## 分布
	&item_Intestinal_Flora('pangjun',$Data,$list->{pangjun},$vars,$dir,$language);
	&item_Intestinal_Flora('shoujun',$Data,$list->{shoujun},$vars,$dir,$language);


	return;
}

#&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language);
sub item_Intestinal_Flora {
	my ($prefix,$Data,$list,$vars,$dir,$language)=@_;
	
	##
	my (@good,@bad,@gooditem,@baditem);
	my (@goodDesc,@badDesc,@goodDisease,@goodSource,@badDisease,@goodEffect,@badEffect,@badSource);
	my (@biotics_good,@biotics_bad,@probiotics_sup,@prebiotics_sup,@food_sup);
	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;

	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## Filt if the bacteria genus percent is zero
		# next if ($Data->{$meta}{sample} == 0);

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		$item=&get_per_count($item);
		push @{$vars->{changdaojunqun}{$prefix}{item}},$item;

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
		$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
		$vars->{changdaojunqun}{$prefix}{badN}+=$addN;

		## get goodDesc,badDesc
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			push @baditem,$item;
			### Get risk disease info
			push @badSource,$vd->{$language}{summary}{source} if (defined $vd->{$language}{summary}{source} && $vd->{$language}{summary}{source} ne "");
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			#$item->{name}=~/(.*)\（.*/;
			#unshift @biotics_bad,$1;
			unshift @biotics_bad,$item->{name};
			&process_PreProbiotics_and_food($vars->{age},\@prebiotics_sup,\@probiotics_sup,\@food_sup,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc}) if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
		}
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			push @gooditem,$item;
			### Get risk disease info
			push @goodSource,$vd->{$language}{summary}{source} if (defined $vd->{$language}{summary}{source} && $vd->{$language}{summary}{source} ne "");
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			#$item->{name}=~/(.*)\（.*/;
			#unshift @biotics_good,$1;
			unshift @biotics_good,$item->{name};
		}
	}

	my $prebiotics_sup_uniq=&uniq_array(\@prebiotics_sup);
	my $probiotics_sup_uniq=&uniq_array(\@probiotics_sup);
	my $food_sup_uniq=&uniq_array(\@food_sup);

	##
	# push @{$vars->{changdaojunqun}{$prefix}{abnormalitem}},@baditem,@gooditem;
	## good,bad
	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	## goodSource,badSource
	($vars->{changdaojunqun}{$prefix}{goodSource},$vars->{changdaojunqun}{$prefix}{badSource})=&process_good_badSource(\@goodSource,\@badSource);

	($vars->{changdaojunqun}{$prefix}{biotics_good},$vars->{changdaojunqun}{$prefix}{biotics_bad})=&process_good_bad(\@biotics_good,\@biotics_bad);
	($vars->{changdaojunqun}{$prefix}{prebiotics_sup},$vars->{changdaojunqun}{$prefix}{probiotics_sup},$vars->{changdaojunqun}{$prefix}{food_sup})=&process_bad_sups($prebiotics_sup_uniq,$probiotics_sup_uniq,$food_sup_uniq);
	## goodDisease,badDisease,goodEffect,badEffect
	($vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{badEffect})=&process_goodbad_DiseaseEffect(\@goodDisease,\@badDisease,\@goodEffect,\@badEffect);
	## badDesc,goodDesc
	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);

	return;
}

#&Intestinal_Pathogen($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Pathogen {
	my ($Data,$list,$vars,$dir,$language,$species,$summary)=@_;
	
	##
	$vars->{changdaojunqun}{pathogen}{name}="致病菌";
	$vars->{changdaojunqun}{pathogen}{ename}="PATHOGEN";
	### num
	$vars->{changdaojunqun}{pathogen}{num}=0;
	
	## 分布
	&item_Intestinal_Pathogen('pathogen',$Data,$list,$vars,$dir,$language,$species);

	return;
}

#&item_Intestinal_Pathogen($prefix,$Data,$list,$vars,$dir,$language,$species);
sub item_Intestinal_Pathogen {
	my ($prefix,$Data,$list,$vars,$dir,$language,$species)=@_;
	
	##
	my (@good,@bad);
	my (@goodDesc,@badDesc,@goodDisease,@badDisease,@goodEffect,@badEffect);
	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;
	##
	my @patho_names=();
	my @norm_pathogen=();
	my @bad_pathogen=();
	my %pathogen2oriid=();
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## Filt if the bacteria genus percent is zero
		# next if ($Data->{$meta}{sample} == 0 && $prefix eq "pathogen");
		next if ($Data->{$meta}{sample} == 0);

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);

		## classify pathogen exceed standard and normal interval
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			push @bad_pathogen, $item;
			$pathogen2oriid{$item->{name}} = $meta;
		}
		else {
			push @norm_pathogen, $item;
			$pathogen2oriid{$item->{name}} = $meta;
		}
	}

	my $pathogen_max = 3;
	my $pathogen_number = 1;
	if (@bad_pathogen) {
		foreach my $item (@bad_pathogen) {
			if ($pathogen_number <= $pathogen_max) {
				push @patho_names, $item->{name};
				$vars->{changdaojunqun}{pathogen}{num} += 1;
				push @{$vars->{changdaojunqun}{$prefix}{item}},$item;
				my $vd=XMLin("$dir/$pathogen2oriid{$item->{name}}.xml",NoAttr=>1,SuppressEmpty => "");
				## get abnormalN,badN,good,bad
				my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
				$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
				$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
				$vars->{changdaojunqun}{$prefix}{badN}+=$addN;

				## Get risk disease & badDesc info
				push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
				push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
				push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
				$pathogen_number ++;
			}
		}
	}

#	if (@norm_pathogen && $pathogen_number <= $pathogen_max) {
#		foreach my $item (@norm_pathogen) {
#			if ($pathogen_number <= $pathogen_max) {
#				push @patho_names, $item->{name};
#				$vars->{changdaojunqun}{pathogen}{num} += 1;
#				push @{$vars->{changdaojunqun}{$prefix}{item}},$item;
#				my $vd=XMLin("$dir/$pathogen2oriid{$item->{name}}.xml",NoAttr=>1,SuppressEmpty => "");
#				## get abnormalN,badN,good,bad
#				my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
#				$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
#				$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
#				$vars->{changdaojunqun}{$prefix}{badN}+=$addN;
#
#				## Get risk disease & badDesc info
#				push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
#				push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
#				push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
#				$pathogen_number ++;
#			}
#		}
#	}

	#&undected_pathogen($prefix, $vars) if ($vars->{changdaojunqun}{pathogen}{num} == 0);

	($vars->{changdaojunqun}{$prefix}{nameDesc})=&process_pathogen_names(\@patho_names);

	## good,bad
	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	## goodDisease,badDisease,goodEffect,badEffect
	($vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{badEffect})=&process_goodbad_DiseaseEffect(\@goodDisease,\@badDisease,\@goodEffect,\@badEffect);
	## badDesc,goodDesc
	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);
	
	return;
}

# &undected_pathogen($vars);
#sub undected_pathogen {
#	my ($prefix, $vars) = @_;
#
#	my %undected_item1 = (
#		"name" => "脆弱拟杆菌",
#		"desc" => "可能导致菌血症、腹内感染、腹膜炎",
#		"risk" => "正常",
#		"level" => "l0",
#	);
#	my %undected_item2 = (
#		"name" => "艰难梭菌",
#		"desc" => "可引起艰难梭菌感染，导致腹痛腹泻等",
#		"risk" => "正常",
#		"level" => "l0",
#	);
#	my %undected_item3 = (
#		"name" => "肠道沙门氏菌",
#		"desc" => "可能引发肠炎，导致腹泻、发烧或腹部痉挛",
#		"risk" => "正常",
#		"level" => "l0",
#	);
#	my %undected_item4 = (
#		"name" => "空肠弯曲杆菌",
#		"desc" => "可导致肠胃炎等",
#		"risk" => "正常",
#		"level" => "l0",
#	);
#	push @{$vars->{changdaojunqun}{$prefix}{item}}, (\%undected_item1, \%undected_item2, \%undected_item3, \%undected_item4);
#}

#&process_Pre_and_Probiotics(\@prebiotics_sup,\@probiotics_sup,$vd->{$language}{suggestion}{$item->{level}}{actiondesc})
sub process_PreProbiotics_and_food {
	my ($age,$probiotics_sup,$prebiotics_sup,$food_sup,$Actdesc)=@_;

	my @Actarray = split/\n/,$Actdesc;
	foreach my $Act (@Actarray) {
		$Act =~s/^\s+|\s+$//;
		next if ($Act eq "");
		if ($Act =~ /(\d+)\\textasciitilde (\d+)）：/) {
			if ($age >= $1 && $age <= $2) {
				if ($Act =~/^益生元/) {
					$Act=~/^[^:]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$probiotics_sup},@sups;
				}
				if ($Act =~/^益生菌/) {
					$Act=~/^[^:]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$prebiotics_sup},@sups;
				}
				if ($Act =~/^食物/) {
					$Act=~/^[^:]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$food_sup},@sups;
				}
			}
		}
		elsif ($Act =~ /(\d+)\+）：/) {
			if ($age >= $1) {
				if ($Act =~/^益生元/) {
					$Act=~/^[^:]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$probiotics_sup},@sups;
				}
				if ($Act =~/^益生菌/) {
					$Act=~/^[^:]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$prebiotics_sup},@sups;
				}
				if ($Act =~/^食物/) {
					$Act=~/^[^:]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$food_sup},@sups;
				}
			}
		}
		else {
			if ($Act =~/^益生元/) {
				$Act=~/^[^:]+：(.*)/;
				my @sups=split/；/,$1;
				unshift @{$probiotics_sup},@sups;
			}
			if ($Act =~/^益生菌/) {
				$Act=~/^[^:]+：(.*)/;
				my @sups=split/；/,$1;
				unshift @{$prebiotics_sup},@sups;
			}
			if ($Act =~/^食物/) {
				$Act=~/^[^:]+：(.*)/;
				my @sups=split/；/,$1;
				unshift @{$food_sup},@sups;
			}
		}
	}
}

#my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
sub get_item {
	my ($vd,$Data,$language,$key)=@_;
	my %item=();
	
	##
	$item{value}=(int($Data->{sample}*100000))/100 if (defined $Data->{sample});
	$item{Val30}=$Data->{Val30} if (defined $Data->{Val30});
	$item{Val70}=$Data->{Val70} if (defined $Data->{Val70});
	$item{level}=$Data->{level};

	#$item{per100}=sprintf "%.2f",$Data->{per}*100;
	$item{per100}= int($Data->{per}*100);
	$item{per100} = "1" if ($Data->{per}*100 < 1);
	$item{per100} = "99" if ($item{per100} >= 100);
	if ($key eq "tjsjcjunqunduoyangxing") {
		$item{value}=sprintf "%.2f",$Data->{sample};
		$item{mean}=sprintf "%.2f",$Data->{ref}{mean};

		$item{Val30}=sprintf "%.2f",$Data->{Val30}/1000;
		$item{Val70}=sprintf "%.2f",$Data->{Val70}/1000;
	}

	if ($key eq "tjsjcyouyijun" || $key eq "tjsjcyouhaijun") {
		$item{value} = "0.01" if ($item{value} eq "0.00" && $key eq "tjsjcyouyijun");
	}

	$item{per}=sprintf "%.4f",$Data->{per};
	if ($Data->{per} eq "-") {
		print "$key\n\n";
	}

	$item{range}=$Data->{range} if (exists $Data->{range});
	$item{name}=$vd->{$language}{title};
	$item{desc}=$vd->{$language}{summary}{desc};
	$item{descfull}=$vd->{$language}{summary}{descfull};
	$item{source}=$vd->{$language}{summary}{source};
	if ($Data->{level} ne '-') {
		$item{risk}=$vd->{$language}{suggestion}{healthy}{$Data->{level}}{effect};
		$item{actionrule}=$vd->{$language}{suggestion}{healthy}{$Data->{level}}{actionrule} if ($vd->{$language}{suggestion}{healthy}{$Data->{level}}{actionrule} ne "" && $key eq "tjsjcyouyijun");
		if ($vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanism} ne "") {
			$item{riskdesc}=$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanism};
			$item{riskdescfull}=$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanismfull} if ($key eq "tjsjcjunqunduoyangxing");
		}
		else {
			$item{riskdesc}='-';
		}
	}
	else {
		$item{risk}='-';
		$item{riskdesc}='-';
	}
	
	return (\%item);
}

sub get_per_count {
	my ($item) = @_;

	$item->{cnt1}=int($item->{per100}/10);

	my $remainder=sprintf "%.2f", (($item->{per100}*100)%1000)/100;
	if ($remainder == 0) {
		$item->{cnt2}= (10 - $item->{cnt1});
		$item->{cnt3}= 0;
	}
	else {
		$item->{cnt2}= (10 - $item->{cnt1} - 1);
		$item->{cnt3}=int($remainder+0.5);
		$item->{cnt3}=9 if ($item->{cnt3} == 10);
	}

	return($item);
}

#my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
sub get_good_bad_array {
	my ($item,$good,$bad)=@_;
	my ($add,$addN,$addG)=(0,0,0);

	##
	if (! defined $item->{risk}) {
		$add=1;
	}
	else {
		if ($item->{risk} ne "正常") {
			$add=1;
			if ($item->{risk} eq "有益") {
				push @{$good},$item->{name};
				$addG=1;
			}
			elsif ($item->{risk} eq "有害") {
				push @{$bad},$item->{name};
				$addN=1;
			}
		}
	}

	return ($add,$addN,$addG);
}

#($vars->{yingyanggongneng}{desc}{$prefix}{good},$vars->{yingyanggongneng}{desc}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
sub process_good_bad {
	my ($goodA,$badA)=@_;
	my ($good,$bad)=('-','-');
	
	##
	if (@{$goodA}) {
		if (scalar @{$goodA} <= 5) {
			$good=join "、",@{$goodA};
		}
		else {
			$good=join("、",$goodA->[0],$goodA->[1],$goodA->[2],$goodA->[3],$goodA->[4]);
			#$good.="等";
		}
	}
	##
	if (@{$badA}) {
		if (scalar @{$badA} <= 5) {
			$bad=join "、",@{$badA};
		}
		else {
			$bad=join("、",$badA->[0],$badA->[1],$badA->[2],$badA->[3],$badA->[4]);
			#$bad.="等";
		}
	}
	
	return ($good,$bad);
}

#($vars->{changdaojunqun}{desc}{$prefix}{nameDesc})=&process_pathogen_names(\@patho_names);
sub process_pathogen_names {
	my ($namesA)=@_;
	my ($nameDesc)=('-');
	
	##
	if (@{$namesA}) {
		if (scalar @{$namesA} <= 2) {
			$nameDesc=join "、",@{$namesA};
		}
		else {
			$nameDesc=join("、",$namesA->[0],$namesA->[1]);
			$nameDesc.="等";
		}
	}
	
	return ($nameDesc);
}

#($vars->{yingyanggongneng}{desc}{$prefix}{prebiotics_sup},$vars->{yingyanggongneng}{desc}{$prefix}{probiotics_sup},$vars->{yingyanggongneng}{desc}{$prefix}{food_sup})=&process_bad_sups(\@pre,\@pro,\@food);
sub process_bad_sups {
	my ($prebioticsA,$probioticsA,$foodA)=@_;
	my ($prebiotics_sup,$probiotics_sup,$food_sup)=('-','-','-');
	
	##
	if (@{$prebioticsA}) {
		if (scalar @{$prebioticsA} <= 5) {
			$prebiotics_sup=join "、",@{$prebioticsA};
		}
		else {
			$prebiotics_sup=join("、",$prebioticsA->[0],$prebioticsA->[1],$prebioticsA->[2],$prebioticsA->[3],$prebioticsA->[4]);
			#$prebiotics_sup.="等";
		}
	}
	##
	if (@{$probioticsA}) {
		if (scalar @{$probioticsA} <= 5) {
			$probiotics_sup=join "、",@{$probioticsA};
		}
		else {
			$probiotics_sup=join("、",$probioticsA->[0],$probioticsA->[1],$probioticsA->[2],$probioticsA->[3],$probioticsA->[4]);
			#$probiotics_sup.="等";
		}
	}
	##
	if (@{$foodA}) {
		if (scalar @{$foodA} <= 5) {
			$food_sup=join "、",@{$foodA};
		}
		else {
			$food_sup=join("、",$foodA->[0],$foodA->[1],$foodA->[2],$foodA->[3],$foodA->[4]);
			#$food_sup.="等";
		}
	}
	
	return ($prebiotics_sup,$probiotics_sup,$food_sup);
}

#&check_file($datadir,$indir,$referdir,$confdir,$barcode);
sub check_file {
	my ($datadir,$indir,$referdir,$confdir,$barcode)=@_;
	
	my $file;
	my $check=0;

	## check reference_population
	#
	$file="$referdir/refer.extract.xls";
	$check+=&checkF($file);
	#
	$file="$referdir/refer.extract.stat.xls";
	$check+=&checkF($file);
	
	## check input data file
	#
	$file="$datadir/$barcode".".extract.xls";
	$check+=&checkF($file);
	
	## check sample information file
	#
	$file="$datadir/$barcode".".csv";
	$check+=&checkF($file);
	
	## check config file
	#
	$file="$confdir/obesity_balance.conf";
	$check+=&checkF($file);

	##
	if ($check > 0) {
		print STDOUT "Files not exist \n";
		die;
	}
	
	return;
}

#$check+=&checkF($file);
sub checkF {
	my ($file)=@_;
	
	##
	my $check=0;
	if (! -f $file) {
		print STDOUT "file not exists : $file\n";
		$check=1;
	}
	
	return($check);
}

#my ($samInfo)=&load_sample_info($barcode,$datadir,\%vars);
sub load_sample_info {
	my ($barcode,$datadir,$vars)=@_;
	my %samInfo=();
	
	##
	my $fatN=-1;
	$/="\n";
	my $file="$datadir/$barcode".".csv";
	open IN,$file or die $!;
	my $head=<IN>;
	my @head=split /\s+/,$head;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @unit=split /\t/,$_;
		for (my $i=0;$i<@head;$i++) {
			if ($head[$i] eq "barcode" || $head[$i] eq "sex" || $head[$i] eq "weight" || $head[$i] eq "height" || $head[$i] eq "batch" || $head[$i] eq "receivedate" || $head[$i] eq "testdate") 
			{
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
			}
			elsif ($head[$i] eq "name") {
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
			}
			elsif ($head[$i] eq "age") {
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$vars->{'agestr'}=$unit[$i];
				if ($unit[$i] =~/未知/) {
					$samInfo{$head[$i]}=30;
					$vars->{$head[$i]}=30;
				}
			}
			elsif ($head[$i] eq "source") {
				if ($unit[$i] eq "") {
					$unit[$i]=decode("UTF-8","-");
				}
				else {
					$unit[$i]=decode("UTF-8",$unit[$i]);
				}
				$unit[$i]="人和未来" if ($unit[$i]=~/人和未来/);
				$unit[$i]="新疆医科大学第一附属医院" if ($unit[$i]=~/新疆医科大学第一附属医院/);
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
			}
			elsif ($head[$i] eq "version") {
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
			}
			else {
				$head[$i]=decode("UTF-8",$head[$i]);
				$head[$i]="2型糖尿病" if ($head[$i] =~ /糖尿病/);
				$head[$i]="炎症性肠病" if ($head[$i] =~ /肠炎/ || $head[$i] =~ /炎症性肠病/ || $head[$i] =~ /腹泻/);
				$head[$i]="粪便颜色" if ($head[$i] =~ /粪便颜色/ || $head[$i] =~ /颜色/);
				if ($head[$i] eq "粪便质地" || $head[$i] eq "粪便颜色" || $head[$i] =~ /历次检测编号/) {
					$samInfo{$head[$i]}=decode("UTF-8",$unit[$i]);
				}
				else {
					$samInfo{disease}{$head[$i]}{'risk'}=decode("UTF-8",$unit[$i]);
				}
			}
		}
	}
	close IN;
	
	## BMI
	my $BMI="NA";
	$BMI=sprintf "%.2f",$samInfo{weight}/(($samInfo{height}/100)**2) if ($samInfo{height} != 0);
	$vars->{BMI} = $BMI;
	my $fat="肥胖";
	if ($BMI ne "NA" && $BMI > 24) {
		$samInfo{disease}{$fat}{risk}='Y';
	}
	else {
		$samInfo{disease}{$fat}{risk}='N';
	}
	
	return (\%samInfo);
}

#my ($listH)=&load_list($indir,$samInfo);
sub load_list {
	my ($indir,$samInfo,$confdir)=@_;
	my %listH=();
	my %term2oriid=();
	my %balancecomposition=();
	my %compositionweigth=();
	
	## pcoress sex
	$samInfo->{sex}='A' if (! defined $samInfo->{sex} || ($samInfo->{sex} ne 'F' && $samInfo->{sex} ne 'M'));
	
	my $file;
	## advise.list.order$sex
	$file="$indir/advise.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"advise");
	}
	## distribution.list.order$sex
	$file="$indir/distribution.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_dis_list($file,\%listH,\%term2oriid);
	}
	## overview.list.order$sex
	$file="$indir/overview.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"overview");
	}
	## pathogens.list.order$sex
	$file="$indir/pathogens.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"pathogens");
	}
	## metasummary.list.order$sex
	$file="$indir/metasummary.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"metasummary");
	}
	## disease.list.order$sex
	$file="$indir/disease.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"disease");
	}
	
	## microbiota balance term composition of oriid and weigth
	$file="$confdir/obesity_balance.conf";
	if (-f $file) {
		&get_balance_composition($file, \%balancecomposition, \%compositionweigth);
	}

	return (\%listH,\%term2oriid,\%balancecomposition,\%compositionweigth);
}

#&get_balance_composition($file,\%balancecomposition);
sub get_balance_composition {
	my ($file,$balancecomposition,$compositionweigth)=@_;
	
	##
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($balance_oriid, $composition_oriid, $weigth, $effects) = (split/\t/,$_)[0,1,2,3];
		$balancecomposition->{$balance_oriid}{$composition_oriid} = $effects;
		$compositionweigth->{$composition_oriid} = $weigth;
	}
	$/="\n";
	close IN;

	return;
}

#&get_list($file,\%listH,"meta");
sub get_list {
	my ($file,$listH,$term2oriid,$type)=@_;
	
	##
	$/="\#\#\#";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my ($section1,$units)=split /\n/,$_,2;
		$section1=(split /\s+/,$section1)[0];
		$section1=decode('UTF-8',$section1);
		my @unit=split /\#\#/,$units;
		for (my $j=1;$j<@unit;$j++) {
			my @lines=split /\n/,$unit[$j];
			my $section2;
			if ($lines[0] eq "") {
				$section2="NA";
			}
			else {
				$section2=(split /\s+/,$lines[0])[0];
				$section2=decode('UTF-8',$section2);
			}
#			$listH->{$type}{$section1}{$j}{name}=$section2;
			for (my $i=1;$i<@lines;$i++) {
				my @aa=split /\t/,$lines[$i];
				my ($oriid,$en,$cn);
				$oriid=$aa[0];
				$cn=$aa[-1];
				$cn=decode('UTF-8',$cn);
				if (scalar @aa == 3) {
					$en=$aa[1];
				}
				$listH->{$type}{$oriid}{order}=$i;
				$listH->{$type}{$oriid}{cn}=$cn;
				if (defined $en) {
					$listH->{$type}{$oriid}{en}=$en;
					$term2oriid->{$en}=$oriid;
				}
			}
		}
	}
	$/="\n";
	close IN;
	
	return;
}

#&get_dis_list($file,\%listH,"meta");
sub get_dis_list {
	my ($file,$listH,$term2oriid)=@_;

	my $type = "distribution";
	##
	$/="\#\#\#";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my ($section1,$units)=split /\n/,$_,2;
		$section1=(split /\s+/,$section1)[0];
		$section1=decode('UTF-8',$section1);
		my @unit=split /\#\#/,$units;
		for (my $j=1;$j<@unit;$j++) {
			my @lines=split /\n/,$unit[$j];
			my $section2=(split /\s+/,$lines[0])[0];
			$section2=decode('UTF-8',$section2);
			$type = "pangjun" if ($section2 eq "胖菌");
			$type = "shoujun" if ($section2 eq "瘦菌");
#			$listH->{$type}{$section1}{$j}{name}=$section2;
			for (my $i=1;$i<@lines;$i++) {
				my @aa=split /\t/,$lines[$i];
				my ($oriid,$en,$cn);
				$oriid=$aa[0];
				$cn=$aa[-1];
				$cn=decode('UTF-8',$cn);
				if (scalar @aa == 3) {
					$en=$aa[1];
				}
				$listH->{distribution}{$type}{$oriid}{order}=$i;
				$listH->{distribution}{$type}{$oriid}{cn}=$cn;
				if (defined $en) {
					$listH->{distribution}{$type}{$oriid}{en}=$en;
					$term2oriid->{$en}=$oriid;
				}
			}
		}
	}
	$/="\n";
	close IN;
	
	return;
}

#&calculate_balance_level_per($Data,$balancecomposition,$level_set);
sub calculate_balance_level_per {
	my ($Data, $balancecomposition, $compositionweigth, $level_set) = @_;

	foreach my $balance_oriid (keys %{$balancecomposition}) {
		my $item_score = 0;
		my $item_num = 0;
		# my $item_val = 0;
		foreach my $oriid (sort keys %{$balancecomposition->{$balance_oriid}}) {
			if (exists $Data->{$oriid}) {
				$item_score += $Data->{$oriid}{per}*$compositionweigth->{$oriid} if ($balancecomposition->{$balance_oriid}{$oriid} eq "positive");
				$item_score += (1 - $Data->{$oriid}{per})*$compositionweigth->{$oriid} if ($balancecomposition->{$balance_oriid}{$oriid} eq "negative");
				# $item_val += $Data->{$oriid}{sample};
				$item_num ++;
			}
			else {
				print STDERR "WARN: Unexists of $oriid of overview index $balance_oriid...\n\n";
			}
		}
		$Data->{$balance_oriid}{per} = $item_score/$item_num;
		# $Data->{$balance_oriid}{sample} = $item_val;

		$Data->{$balance_oriid}{level} = &get_balance_level($Data->{$balance_oriid}{per}, $level_set);
	}
}

#my ($Data)=&calculate_meta_level_per($datadir,$indir,$barcode,$level_set);
sub calculate_meta_level_per {
	my ($datadir,$indir,$referdir,$barcode,$term2oriid,$compositionweigth,$level_set)=@_;
	my $Data=();
	
	## load sample and reference data
	($Data)=&load_sample_reference_data($datadir,$indir,$referdir,$barcode,$term2oriid);

	## get level,per,range
	&calculate_level_per($Data,$compositionweigth,$level_set);
	
	return ($Data);
}

#&calculate_level_per($Data,$level_set)
sub calculate_level_per {
	my ($Data,$compositionweigth,$level_set)=@_;
	
	##
	foreach my $meta (keys %{$Data}) {
		# $meta weigth
		my $weigth = 1;
		$weigth = $compositionweigth->{$meta} if (exists $compositionweigth->{$meta});

		# get level
		($Data->{$meta}{level},$Data->{$meta}{range},$Data->{$meta}{plusSD},$Data->{$meta}{minusSD})=&get_range($Data->{$meta},$meta,$level_set);

		## get per and level
		($Data->{$meta}{per},$Data->{$meta}{level},$Data->{$meta}{Val30},$Data->{$meta}{Val70})=&get_per($Data->{$meta},$meta,$weigth,$level_set);
	}
	
	return;
}

#$Data->{$meta}{per}=&get_per($Data->{$meta}, $key);
sub get_per {
	my ($meta, $key, $weigth, $level_set)=@_;
	my ($per,$level,$val30,$val70)=('-','-',0,0);

	my $weight_per30 = 0.3*$weigth;
	my $weight_per70 = 0.7*$weigth;

	## use norm distribution calculater
	#if ((defined $meta->{ref}{mean} && $meta->{ref}{mean} != 0) && (defined $meta->{ref}{sd} && $meta->{ref}{sd} != 0) && exists $meta->{sample}) {
	#	## use R pnorm
	#	$R->set('v', $meta->{sample});
	#	$R->set('y', $meta->{ref}{mean});
	#	$R->set('z', $meta->{ref}{sd});
	#	$R->run(q`prob <- pnorm(v, y, z)`);
	#	$per = $R->get('prob');
	#}

	## use ecdf calculate percentile
	if (exists $meta->{sample} && exists $meta->{ref}) {
		## use R ecdf
		$R->set('v', $meta->{ref}{values});
		$R->set('y', $meta->{sample});
		$R->run(q`percentile <- ecdf(v)(y)`);
		$per = $R->get('percentile');
		$R->set('l', $weight_per30);
		$R->run(q`val30 <- as.numeric(quantile(v,l))`);
		$val30 = $R->get('val30');
		$R->set('m', $weight_per70);
		$R->run(q`val70 <- as.numeric(quantile(v,m))`);
		$val70 = $R->get('val70');
	}

	## adjust some items value
	#if ($key eq "tjsjcaikemanshijunshu") {
	#	if ($meta->{sample} == 0) {
	#		$per = rand(0.25)+0.01;
	#	}
	#}

	$val30 = (int($val30*100000)+1)/100;
	$val70 = (int($val70*100000)+1)/100;

	#my $level=&get_per_level ($per*$weigth, $level_set);
	if ($meta->{sample}*1000 < $val30) {
		$level = "l2";
	}
	elsif ($meta->{sample}*1000 >= $val30 && $meta->{sample}*1000 < $val70) {
		$level = "l3";
	}
	elsif ($meta->{sample}*1000 >= $val70) {
		$level = "l4";
	}

	return ($per,$level,$val30,$val70);
}

#&get_balance_level($per,$level)
sub get_balance_level {
	my ($per, $level_set) = @_;
	my $level = '-';

	if ($per < 0.1) {
		$level = "l1";
	}
	#elsif ($per <= 0.2) {
	elsif ($per < 0.3) {
		$level = "l2";
	}
	#elsif ($per < 0.8) {
	elsif ($per < 0.7) {
		$level = "l3";
	}
	elsif ($per < 0.9) {
		$level = "l4";
	}
	else {
		$level = "l5";
	}
	$level = $level_set if ($level_set ne "N");

	return($level);
}

#&get_per_level($per,$level)
sub get_per_level {
	my ($per, $level_set) = @_;
	my $level = '-';

	if ($per < 0.1) {
		$level = "l1";
	}
	#elsif ($per <= 0.2) {
	elsif ($per < 0.3) {
		$level = "l2";
	}
	else {
		$level = "l3";
	}
	$level = $level_set if ($level_set ne "N");

	return($level);
}

#($Data->{$meta}{level},$Data->{$meta}{range})=&get_range($Data->{$meta},$meta,$level_set);
sub get_range {
	my ($meta,$k,$level_set)=@_;
	my $level="-";
	my $range="-";
	
	##
	if (exists $meta->{sample} && exists $meta->{ref}) {
		my ($mean,$sd)=($meta->{ref}{mean},$meta->{ref}{sd});
		my $lowbound = ($mean-2*$sd)*100;
		$range=(sprintf "%.2f",($mean-2*$sd)*100).'~'.(sprintf "%.2f",($mean+2*$sd)*100) if ($lowbound >= 0.01);
		$range=(sprintf "%d",0).'~'.(sprintf "%.2f",($mean+2*$sd)*100) if ($lowbound < 0.01);
		$range=(sprintf "%.2f",$mean-2*$sd).'~'.(sprintf "%.2f",$mean+2*$sd) if ($k eq "tjsjcjunqunduoyangxing");
		$range=&replace($range);
		if ($meta->{sample} <= $mean-2*$sd) {
			$level='l1';
		}
		elsif ($meta->{sample} <= $mean-$sd) {
			$level='l2';
		}
		elsif ($meta->{sample} <= $mean+$sd) {
			$level='l3';
		}
		elsif ($meta->{sample} <= $mean+2*$sd) {
			$level='l4';
		}
		else {
			$level='l5';
		}
		$level = $level_set if ($level_set ne "N");
	}

	return ($level,$range);
}

#($Data)=&load_sample_reference_data($datadir,$indir,$referdir,$barcode,$term2oriid);
sub load_sample_reference_data {
	my ($datadir,$indir,$referdir,$barcode,$term2oriid)=@_;
	my %Data=();

	my $file;
	## get sample values
	$file="$datadir/$barcode".".extract.xls";
	&get_sample_values($file,$term2oriid,\%Data);
	
	## get ref samples values && sort
	$file="$referdir/refer.extract.xls";
	&get_ref_values($file,$term2oriid,\%Data);

	## get ref mead && sd
	$file="$referdir/refer.extract.stat.xls";
	&get_ref_mean_sd($file,$term2oriid,\%Data);

	return (\%Data);
}

#&get_sample_values($file,$term2oriid,\%Data);
sub get_sample_values {
	my ($file,$term2oriid,$Data)=@_;

	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$value)=(split /\t/,$_)[0,1];
		next if (!exists $term2oriid->{$id});
		$id = $term2oriid->{$id};
		$Data->{$id}{sample}=$value;
	}
	close IN;

	return;
}

#&get_ref_mean_sd($file,$term2oriid,\%Data);
sub get_ref_mean_sd {
	my ($file,$term2oriid,$Data)=@_;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @unit=split /\t/,$_;
		next if (!exists $term2oriid->{$unit[0]});
		$unit[0] = $term2oriid->{$unit[0]};
		$Data->{$unit[0]}{ref}{mean}=$unit[3];
		$Data->{$unit[0]}{ref}{sd}=$unit[4];
	}
	close IN;
	
	return;
}

#&get_ref_values($file,$term2oriid,\%Data);
sub get_ref_values {
	my ($file,$term2oriid,$Data)=@_;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @unit=split /\t/,$_;
		next if (!exists $term2oriid->{$unit[0]});
		$unit[0] = $term2oriid->{$unit[0]};
		my $zeroN=0;
		for (my $i=1;$i<@unit;$i++) {
			push @{$Data->{$unit[0]}{ref}{values}},$unit[$i];
			$zeroN++ if ($unit[$i] == 0);
		}
		# array sort
		@{$Data->{$unit[0]}{ref}{values}}=sort {$a<=>$b} @{$Data->{$unit[0]}{ref}{values}};
		# check zero
		if ($zeroN + 1 == scalar @unit) {
			$Data->{$unit[0]}{ref}{checkZero}='Y';
		}
		else {
			$Data->{$unit[0]}{ref}{checkZero}='N';
		}
	}
	close IN;
	
	return;
}

#my $xml=&read_xml($xmlfile);
sub read_xml {
	my ($xmlfile)=@_;
	
	##
	my $xml=XMLin($xmlfile,NoAttr=>1,ForceArray => ['item','abnormalitem','disease','spot'], KeyAttr => []);
	
	return ($xml);
}

### replace
sub replace {
        my ($js)=@_;
        my $utfcode = '';
        ##
        $js=~s/{/\\{/g;
        $js=~s/}/\\}/g;
        $js=~s/%/{\\%}/g;
        $js=~s/\n/\n\n/g;
        $js=~s/\$/\\\$/g;
        $js=~s/~/\\textasciitilde /g;
        $js=~s/_/\\_/g;
        $js=~s/#/\\#/g;
        $js=~s/\//\{\/\}/g;

        return ($js);
}

#&output_xml($vars,$barcode,$dirout);
sub output_xml {
	my ($vars,$barcode,$dirout)=@_;
	
	##
	my $outfile="$dirout/$barcode".".xml";
	open my $out, '>:encoding(utf-8)', $outfile || die $!;
	XMLout($vars,  OutputFile => $out, NoAttr => 1,SuppressEmpty => "", KeyAttr => []);
	close $out;
	
	return;
}

sub AbsolutePath
{       #获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
        $/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";

    Program Function: GIhealth result to xml format;
    Version:    $version
    Contact:    fred routine <fred_routine\@163.com> 
    Program Date:   2016.11.15
    Usage:
      Options:
      -data        <dir>   sample's report result dir for input          forced
      -barcode     <str>   sample's code                                 forced

      -in         <dir>   product set database list and xml dir     optional
                          default "$Bin/db2xml/current";
      -referdir   <dir>   reference pop data dir                    optional
                          default "$Bin/../../reference_data/current";
      -out        <dir>   output file dir                           optional
                          default "./";
      -language   <str>   result language used for report           optional
                          default "CN";
      -confdir    <dir>   report conf dir                           optional
                          default "$Bin/../conf/current";
      -historydir <str>   history test xml result for report dir    optional
                          default "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
      -l          <str>   level set to check report correctness (l1~l5)          optional
                          default "N";

      -h          Help

USAGE
	print $usage;
	exit;
}
