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
$dirout=&AbsolutePath("dir",$dirout);
$datadir=&AbsolutePath("dir",$datadir);
$indir=&AbsolutePath("dir",$indir);
$referdir=&AbsolutePath("dir",$referdir);
$confdir=&AbsolutePath("dir",$confdir);
$historydir=&AbsolutePath("dir",$historydir);

$language||="CN";

$level_set||="N";

################# main process ###################
my %diseaseOrder = (
	"question" => {
			"2型糖尿病" => 3,
			"炎症性肠病" => 1,
			"心脑血管疾病" => 5,
			"结直肠癌" => 4,
			"便秘" => 2,
			"肥胖" => 6,
		},
	"predict" => {
			"2型糖尿病" => 3,
			"炎症性肠病" => 1,
			"心脑血管疾病" => 5,
			"结直肠癌" => 4,
			"便秘" => 2,
			"肥胖" => 6,
		},
	"names" => {
			"2型糖尿病" => "disease3",
			"炎症性肠病" => "disease1",
			"心脑血管疾病" => "disease5",
			"结直肠癌" => "disease4",
			"便秘" => "disease2",
			"肥胖" => "disease6",
		},
);

my %vars=();
## check file and dir
&check_file($datadir,$indir,$referdir,$confdir,$barcode);

## load sample information file
my ($samInfo)=&load_sample_info($barcode,$datadir,\%vars);

## load $type.list.order$sex
my ($listH,$term2oriid)=&load_list($indir,$samInfo);

## calculate meta level and per by sample data and ref data
my ($Data)=&calculate_meta_level_per($datadir,$indir,$referdir,$barcode,$term2oriid,$level_set);

## main process
#&main_process($Data,$listH,\%{$diseaseOrder{question}},\%{$diseaseOrder{predict}},\%vars,$datadir,$indir,$language,$barcode,$samInfo);
&main_process($Data,$listH,\%{$diseaseOrder{question}},\%{$diseaseOrder{predict}},\%{$diseaseOrder{names}},\%vars,$datadir,$indir,$confdir,$language,$barcode,$samInfo,$historydir);

## sample result csv for erp import
&sample_result_csv(\%vars,$listH->{pathogens},$barcode,$dirout);

############################# sub ###############################################

#&main_process($Data,$listH,\%vars,$datadir,$indir,$language,$barcode,$samInfo);
sub main_process {
	my ($Data,$listH,$listQ,$listP,$listN,$vars,$datadir,$indir,$confdir,$language,$barcode,$samInfo,$historydir)=@_;

	## 菌群分布(菌属分布+致病菌)
	# $listH->{"distribution"}; $listH->{"pathogens"}
	# 所有菌属
	my ($species)=&load_species("$confdir/tjsjc_genus_fenbu.conf", "$confdir/genus_latin2zh");
	# summary信息
	my ($summary)=&load_summary("$datadir/$barcode.summary.xls");
	&Intestinal_Flora($Data,$listH->{distribution},$vars,"$indir/XML/distribution/$language",$language,$species,$summary);
	&Intestinal_Pathogen($Data,$listH->{pathogens},$vars,"$indir/XML/pathogens/$language",$language,$species,$summary);
	&get_Intestinal_Flora_pic($vars);

	## 菌群概况(多样性、有益菌、有害菌)+菌群综合指标(有益有害菌)
	# $listH->{"overview"},$listH->{"metasummary"}
	&Micro_Overvwiew($Data,$listH->{overview},$listH->{metasummary},$vars,"$indir/XML/overview/$language","$indir/XML/metasummary/$language",$language,$species,$summary);

	## 菌群代谢(Vitamin 及重要有机小分子代谢)
	# $listH->{"metabolism"}
	# 所有维生素
	my ($Vitamin)=&load_Vitamin("$confdir/tjsjc_Vitamin.conf");
	my ($regulateLow,$badDiseaseTotal)=&Micro_Metabolism($Data,$listH->{metabolism},$vars,"$indir/XML/metabolism/$language",$language,$Vitamin);

	## 膳食方案+肠道调节方案+运动方案
	# $listH->{"advise"};
	&Program($vars,"$indir/XML/advise",$language,$samInfo,$listH->{advise},$listQ,$listP,$listN,$regulateLow,$badDiseaseTotal);

	## 多次检测结果比较
	# $samInfo->{batch}; $samInfo->{'历次检测编号'}
	&Compare_history($vars,$samInfo,$historydir) if ($samInfo->{batch} > 1);

	## 报告日期
	$vars->{reportdate}=strftime("%Y-%m-%d",localtime());
	
	## output xml
	&output_xml($vars,$barcode,$dirout);
	
	return;
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
			$xml->{gaikuang}{curpic} =~ /sun_(.*)\.pdf/;
			$vars->{$batch_str}{gaikuang}{pic}=$1."circle.pdf";
			$xml->{changdaojunqun}{curpic} =~ /sun_(.*)\.pdf/;
			$vars->{$batch_str}{changdaojunqun}{pic}=$1."circle.pdf";
			$xml->{yingyanggongneng}{curpic} =~ /sun_(.*)\.pdf/;
			$vars->{$batch_str}{yingyanggongneng}{pic}=$1."circle.pdf";
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

	$xml->{gaikuang}{curpic} =~ /sun_(.*)\.pdf/;
	$vars->{$batch_str}{gaikuang}{pic}=$1."circle.pdf";
	$xml->{changdaojunqun}{curpic} =~ /sun_(.*)\.pdf/;
	$vars->{$batch_str}{changdaojunqun}{pic}=$1."circle.pdf";
	$xml->{yingyanggongneng}{curpic} =~ /sun_(.*)\.pdf/;
	$vars->{$batch_str}{yingyanggongneng}{pic}=$1."circle.pdf";

	$vars->{$batch_str}{gaikuang}{desc}=\%{$xml->{gaikuang}{desc}};
	&compare_gaikuang($vars, "duoyangxing", $vars->{gaikuang}{desc}{duoyangxing}{per}, $xml->{gaikuang}{desc}{duoyangxing}{per});
	&compare_gaikuang($vars, "youyijun", $vars->{gaikuang}{desc}{youyijun}{per}, $xml->{gaikuang}{desc}{youyijun}{per});
	&compare_gaikuang($vars, "youhaijun", $vars->{gaikuang}{desc}{youhaijun}{per}, $xml->{gaikuang}{desc}{youhaijun}{per});

	$vars->{$batch_str}{changdaojunqun}{fenbu}{goodN}=$xml->{changdaojunqun}{fenbu}{goodN};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{badN}=$xml->{changdaojunqun}{fenbu}{badN};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{good}=$xml->{changdaojunqun}{fenbu}{good};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{bad}=$xml->{changdaojunqun}{fenbu}{bad};

	$vars->{$batch_str}{changdaojunqun}{pathogen}{goodN}=$xml->{changdaojunqun}{pathogen}{goodN};
	$vars->{$batch_str}{changdaojunqun}{pathogen}{badN}=$xml->{changdaojunqun}{pathogen}{badN};
	$vars->{$batch_str}{changdaojunqun}{pathogen}{good}=$xml->{changdaojunqun}{pathogen}{good};
	$vars->{$batch_str}{changdaojunqun}{pathogen}{bad}=$xml->{changdaojunqun}{pathogen}{bad};

	$vars->{$batch_str}{yingyanggongneng}{desc}{goodN}=$xml->{yingyanggongneng}{desc}{goodN};
	$vars->{$batch_str}{yingyanggongneng}{desc}{badN}=$xml->{yingyanggongneng}{desc}{badN};
	$vars->{$batch_str}{yingyanggongneng}{desc}{good}=$xml->{yingyanggongneng}{desc}{good};
	$vars->{$batch_str}{yingyanggongneng}{desc}{bad}=$xml->{yingyanggongneng}{desc}{bad};

}

#&compare_gaikuang($vars, "duoyangxing", $vars->{gaikuang}{desc}{duoyangxing}{per}, $xml->{gaikuang}{desc}{duoyangxing}{per});
sub compare_gaikuang {
	my ($vars, $key, $valnow, $vallast) = @_;

	if (($valnow - $vallast) >= 0.05) {
		$vars->{comparelast}{$key} = "升高";
	}
	elsif (($valnow - $vallast) <= -0.05) {
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
	open my $out, '>:encoding(utf-8)', $outfile || die $!;

	print $out "检测者信息\n检测编号\t姓名\n$vars->{barcode}\t$vars->{name}\n";

#	my $head_str = "";
	my $term_str = "";
	my $result_str = "";
#	my $value_str = "\t";
#	my $range_str = "\t";
#	my $per_str = "\t";

	my $level;
	print $out "菌群概况\n肠道菌群多样性\t有益菌\t有害菌\n";
	$level = &transfer_item_level($vars->{gaikuang}{desc}{duoyangxing}{level});
	$result_str ="$level($vars->{gaikuang}{desc}{duoyangxing}{risk})\t";

	$level = &transfer_item_level($vars->{gaikuang}{desc}{youyijun}{level});
	$result_str .="$level($vars->{gaikuang}{desc}{youyijun}{risk})\t";

	$level = &transfer_item_level($vars->{gaikuang}{desc}{youhaijun}{level});
	$result_str .="$level($vars->{gaikuang}{desc}{youhaijun}{risk})\t";
#	$value_str .= "\t$vars->{gaikuang}{desc}{youhaijun}{value}";
#	$range_str .= "\t$vars->{gaikuang}{desc}{youhaijun}{range}";
#	$per_str .= "\t$vars->{gaikuang}{desc}{youhaijun}{per}";
	$result_str =~s/\s+$//;
	print $out "$result_str\n";

#	$head_str .= "\t菌群分布";
	print $out "菌群分布\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{changdaojunqun}{fenbu}{item}}) {
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

#	$head_str .= "菌群代谢";
	print $out "有机小分子物质代谢\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{yingyanggongneng}{desc}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	my %pathogens_result=();
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		$pathogens_result{$list->{$meta}{cn}} = "未检出(正常)";
	}
	foreach my $item (@{$vars{changdaojunqun}{pathogen}{item}}) {
		$level = &transfer_pathogen_level($item->{level});
		$pathogens_result{$item->{name}} = "$level($item->{risk})";
	}
#	$head_str .= "致病菌含量";
	print $out "致病菌含量\n";
	$term_str = "";
	$result_str = "";
	foreach my $patho_name (sort keys %pathogens_result) {
		$term_str .= "$patho_name\t";
		$result_str .="$pathogens_result{$patho_name}\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	close $out;

	return;
}

#&Program($vars,"$indir/XML",$language,$samInfo,$listH->{diseaseList},$regulateLow,$badDiseaseTotal);
sub Program {
	my ($vars,$dir,$language,$samInfo,$list,$listQ,$listP,$listN,$regulateLow,$badDiseaseTotal)=@_;
	my ($foodKey,$regulateKey,$sportKey,$age)=('-','-','-','-');

	## get food key
	($foodKey)=&get_program_food($samInfo,$vars,$listQ,$listN);
	## get regulate key
	($regulateKey)=&get_program_regulate($samInfo,$vars,$listQ,$listN);
	## get sport key
	($sportKey)=&get_program_food($samInfo,$vars,$listQ,$listN);
	## if no food disease, get age
	($age)=&get_program_age($samInfo->{age});

	## process food and sport
	&process_food_sport($vars,$dir,$foodKey,$sportKey,$age,$language);
	## process regulate
	&process_regulate($vars,$dir,$regulateKey,$regulateLow,$age,$language);
	
	return;
}

#&process_regulate($vars,$dir,$regulateKey,$regulateLow,$language);
sub process_regulate {
	my ($vars,$dir,$regulateKey,$regulate,$age,$language)=@_;
	
	## 补充益生元和益生菌
	my $vd=XMLin("$dir/$language/tjsjctiaojiefangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{regulate}{name} = "营养调节方案";
	#$vars->{fanan}{food}{pic} = "tiaojie.pdf";
	if ($regulateKey=~/condition/) {
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => "调节肠道菌群、补充营养", desc => $vd->{$language}->{microsuggestion}->{$regulateKey}->{$age}->{action}};
	}
	elsif ($regulateKey=~/disease/) {
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => "调节肠道菌群、补充营养", desc => $vd->{$language}->{dieasesuggestion}->{$regulateKey}->{$age}->{action}};
	}
	elsif ($regulateKey=~/age/) {
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => "调节肠道菌群、补充营养", desc => $vd->{$language}->{agesuggestion}->{desc}->{$regulateKey}->{action}};
	}

	## 其他
	push @{$vars->{fangan}{regulate}{desc}{item}},@{$regulate} if (@{$regulate});
	
	return;
}

#&process_food_sport($vars,$dir,$diseaseF,$ageF,$badDiseaseF,$language);
sub process_food_sport {
	my ($vars,$dir,$foodKey,$sportKey,$age,$language)=@_;
	
	## food
	my $vd_food=XMLin("$dir/$language/tjsjcshanshifangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{food}{name} = "膳食方案";
	#$vars->{fangan}{food}{pic} = "shanshi.pdf";
	$vars->{fangan}{food}{desc} = $vd_food->{$language}{dieasesuggestion}{$foodKey}{$age}{action} if ($foodKey=~/disease/);
	$vars->{fangan}{food}{desc} = $vd_food->{$language}{microsuggestion}{$foodKey}{$age}{action} if ($foodKey=~/condition/);
	$vars->{fangan}{food}{desc} = $vd_food->{$language}{agesuggestion}{desc}{$foodKey}{action} if ($foodKey=~/age/);
	## sport
	my $vd_sport=XMLin("$dir/$language/tjsjcyundongfangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{sport}{name} = "运动方案";
	#$vars->{fangan}{sport}{pic} = "yundong.pdf";
	$vars->{fangan}{sport}{desc} = $vd_sport->{$language}{dieasesuggestion}{$sportKey}{$age}{action} if ($sportKey=~/disease/);
	$vars->{fangan}{sport}{desc} = $vd_sport->{$language}{microsuggestion}{$sportKey}{$age}{action} if ($sportKey=~/condition/);
	$vars->{fangan}{sport}{desc} = $vd_sport->{$language}{agesuggestion}{desc}{$sportKey}{action} if ($sportKey=~/age/);
	
	return;
}

#($diseaseF,$diseaseR)=&get_program_regulate($samInfo,$vars->{changdaojunqun}{diseaseHigh},$listQ);
sub get_program_regulate {
	my ($samInfo,$vars,$listQ,$listN)=@_;
	my ($disease,$condition,$age)=('-','-','-');

	### harmful items (多样性、有益菌、有害菌)
	if ($vars->{gaikuang}{desc}{duoyangxing}{level} eq "l1") {
		$condition="condition1";
	}
	elsif ($vars->{gaikuang}{desc}{youyijun}{level} eq "l1") {
		$condition="condition2";
	}
	elsif ($vars->{gaikuang}{desc}{youhaijun}{level} eq "l5") {
		$condition="condition3";
	}
	elsif ($vars->{changdaojunqun}{pathogen}{badN} >= 2) {
		$condition="condition4";
	}
	elsif ($vars->{gaikuang}{desc}{duoyangxing}{level} eq "l2") {
		$condition="condition5";
	}
	elsif ($vars->{gaikuang}{desc}{youyijun}{level} eq "l2") {
		$condition="condition6";
	}
	elsif ($vars->{gaikuang}{desc}{youhaijun}{level} eq "l4") {
		$condition="condition7";
	}
	return ($condition) if ($condition ne '-');

	#### disease
	## sample information
	if ($condition eq '-') {
		foreach my $k (sort {$listQ->{$a} <=> $listQ->{$b}} keys %{$listQ}) {
			if ($samInfo->{disease}{$k}{risk} eq 'Y') {
				$disease="$listN->{$k}";
				return ($disease);
				last;
			}
		}
	}

	### harmful items (多样性、有益菌、有害菌)
	if ($disease eq '-' && $condition eq '-') {
		$age=&get_program_age($samInfo->{age});
		return ($age);
	}

}

#($diseaseF,$diseaseR)=&get_program_disease($samInfo,$vars->{changdaojunqun}{diseaseHigh},$listQ,$listP);
sub get_program_food {
	my ($samInfo,$vars,$listQ,$listN)=@_;
	my ($disease,$condition,$age)=('-','-','-');
	
	#### disease
	## sample information
	foreach my $k (sort {$listQ->{$a} <=> $listQ->{$b}} keys %{$listQ}) {
		if ($samInfo->{disease}{$k}{risk} eq 'Y') {
			$disease="$listN->{$k}";
			return ($disease);
			last;
		}
	}

	### harmful items (多样性、有益菌、有害菌)
	if ($disease eq '-') {
		if ($vars->{gaikuang}{desc}{duoyangxing}{level} eq "l1") {
			$condition="condition1";
		}
		elsif ($vars->{gaikuang}{desc}{youyijun}{level} eq "l1") {
			$condition="condition2";
		}
		elsif ($vars->{gaikuang}{desc}{youhaijun}{level} eq "l5") {
			$condition="condition3";
		}
		elsif ($vars->{changdaojunqun}{pathogen}{badN} >= 2) {
			$condition="condition4";
		}
		elsif ($vars->{gaikuang}{desc}{duoyangxing}{level} eq "l2") {
			$condition="condition5";
		}
		elsif ($vars->{gaikuang}{desc}{youyijun}{level} eq "l2") {
			$condition="condition6";
		}
		elsif ($vars->{gaikuang}{desc}{youhaijun}{level} eq "l4") {
			$condition="condition7";
		}
		return ($condition) if ($condition ne '-');
	}

	### harmful items (多样性、有益菌、有害菌)
	if ($disease eq '-' && $condition eq '-') {
		$age=&get_program_age($samInfo->{age});
		return ($age);
	}

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

#&get_Intestinal_Flora_pic($vars);
sub get_Intestinal_Flora_pic {
	my ($vars)=@_;
	
	##
	$vars->{changdaojunqun}{curpic}="-";
	$vars->{changdaojunqun}{pic}="-";
	if ($vars->{changdaojunqun}{fenbu}{badN} >= 4 || $vars->{changdaojunqun}{pathogen}{badN} >= 1) {
		$vars->{changdaojunqun}{curpic}="sun_red.pdf";
		$vars->{changdaojunqun}{pic}="fenbu_red.pdf";
	}
	elsif ($vars->{changdaojunqun}{fenbu}{badN} >= 1 || $vars->{changdaojunqun}{pathogen}{num} >= 4) {
		$vars->{changdaojunqun}{curpic}="sun_orange.pdf";
		$vars->{changdaojunqun}{pic}="fenbu_orange.pdf";
	}
	else {
		$vars->{changdaojunqun}{curpic}="sun_green.pdf";
		$vars->{changdaojunqun}{pic}="fenbu_green.pdf";
	}
	
	return;
}

#my ($regulateLow,$badDiseaseTotal)=&Digestion_and_Absorption($Data,$listH->{meta}{'消化和吸收'},$vars,"$indir/XML/meta/$language",$language,$Vitamin);
sub Micro_Metabolism {
	my ($Data,$list,$vars,$dir,$language,$vitamin)=@_;
	my @regulateLow=();
	my %badDiseaseTotal=();
	my %TotalbadDiseaseTotal=();
	
	##
	my (@TotalgoodMeta,@TotalbadMeta);                                              ## harmful item of metabolism
	my (@TotalbadDisease,@TotalbadActRule,@TotalbadActDesc);                        ## harmful metabolism related diseases
	$vars->{yingyanggongneng}{name}="肠道菌群营养功能分析";
	$vars->{yingyanggongneng}{ename}="INSUFFICIENCY";

	my @vitaminBad=();
	## 菌群代谢
	$vars->{yingyanggongneng}{desc}{FumeiQBad} = 0;
	my $vitaminBadN=&item_Micro_Metabolism('desc',$Data,$list,$vars,$dir,$language,$vitamin,\@regulateLow,\%badDiseaseTotal,\@TotalbadDisease,\@TotalbadActRule,\@TotalbadActDesc,\@TotalgoodMeta,\@TotalbadMeta,\@vitaminBad);
	&process_vitaminbad($vars,$vitaminBadN,\@vitaminBad);

	## protein、fat、carbon + meta
	($vars->{yingyanggongneng}{explain}{TotalgoodMeta},$vars->{yingyanggongneng}{explain}{TotalbadMeta})=&process_good_bad(\@TotalgoodMeta,\@TotalbadMeta);
	($vars->{yingyanggongneng}{explain}{TotalbadDisease},$vars->{yingyanggongneng}{explain}{TotalbadActRule},$vars->{yingyanggongneng}{explain}{badActDesc})=&process_badDisease_badFood(\@TotalbadDisease,\@TotalbadActRule,\@TotalbadActDesc,\%TotalbadDiseaseTotal);
	
	## get green|yellow|red.pdf
	&get_Micro_Metabolism_pic($vars);
	
	return (\@regulateLow,\%badDiseaseTotal);
}

sub process_vitaminbad {
	my ($vars,$vitaminBadN,$vitaminBad)=@_;

	my $str = "维生素";
	if ($vitaminBadN == 0) {
		$vars->{yingyanggongneng}{desc}{VitaminBadN} = 0;
		$vars->{yingyanggongneng}{desc}{VitaminBad} = "-";
	}
	elsif ($vitaminBadN <= 2) {
		foreach my $vitaminid (@{$vitaminBad}) {
			$vitaminid =~/([a-zA-Z\d]+)/;
			$str .= "$1"."、";
		}
		$str =~s/、$//;
		$vars->{yingyanggongneng}{desc}{VitaminBadN} = $vitaminBadN;
		$vars->{yingyanggongneng}{desc}{VitaminBad} = $str;
	}
	else {
		for (0..1) {
			my $vitaminid = ${$vitaminBad}[$_];
			$vitaminid =~/([a-zA-Z\d]+)/;
			$str .= "$1"."、";
		}
		$str =~s/、$//;
		$vars->{yingyanggongneng}{desc}{VitaminBadN} = $vitaminBadN;
		$vars->{yingyanggongneng}{desc}{VitaminBad} = "$str"."等";
	}
}

#&item_Digestion_and_Absorption('1','zongping',$Data,$list,$vars,$dir,$language,\@Totalgood,\@Totalbad,$vitamin,\@regulateLow,\%badDiseaseTotal,\@TotalbadDisease,\@TotalbadFoodRule,\@TotalbadFoodMore,\@TotalbadFoodLess,\@TotalgoodMeta,\@TotalbadMeta);
sub item_Micro_Metabolism {
	my ($prefix,$Data,$list,$vars,$dir,$language,$vitamin,$regulate,$badDiseaseTotal,$TotalbadDisease,$TotalbadActRule,$TotalbadActDesc,$TotalgoodMeta,$TotalbadMeta,$vitaminBad)=@_;
	
	##
	my $vitaminBadN = 0;
	my (@good,@bad,@gooditem,@baditem);
	my (@goodDesc,@badDesc);
	my (@goodDisease,@badDisease,@goodEffect,@badEffect);
	my (@badActRule,@badActDesc);
	##
	$vars->{yingyanggongneng}{$prefix}{abnormalN}=0;
	$vars->{yingyanggongneng}{$prefix}{goodN}=0;
	$vars->{yingyanggongneng}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}) {
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");

		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		push @{$vars->{yingyanggongneng}{$prefix}{item}},$item;

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{yingyanggongneng}{$prefix}{abnormalN}+=$add;
		$vars->{yingyanggongneng}{$prefix}{goodN}+=$addG;
		$vars->{yingyanggongneng}{$prefix}{badN}+=$addN;

		## get abnormal,badDisease,badFoodRule,badFoodMore,badFoodLess
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			push @baditem,$item;
			### Get risk disease info
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
			push @badActDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
			## get badDesc
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");

		}

		## get goodDesc
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			push @gooditem,$item;
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
		}

		## get regulate
		if ($meta eq "tjsjcfumeiQ") {
			if ($item->{level} eq 'l1' || $item->{level} eq 'l2') {
				$vars->{yingyanggongneng}{$prefix}{FumeiQBad} = 1;
#				push @{$regulate},{name => $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule}, desc => $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc}};
			}
		}
		if (exists $vitamin->{$meta}) {
			if ($item->{level} eq 'l1' || $item->{level} eq 'l2') {
				$vitaminBadN++;
				push @{$vitaminBad}, $vitamin->{$meta};
#				push @{$regulate},{name => $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule}, desc => $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc}};
			}
		}
	}
	#$vars->{yingyanggongneng}{$prefix}{abnormal}{num}=$vars->{yingyanggongneng}{$prefix}{badN};

	##
	push @{$vars->{yingyanggongneng}{$prefix}{abnormalitem}},@baditem,@gooditem;
	## good,bad
	($vars->{yingyanggongneng}{$prefix}{good},$vars->{yingyanggongneng}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	push @{$TotalgoodMeta},@good;
	push @{$TotalbadMeta},@bad;
	push @{$TotalbadDisease},@badDisease;
	push @{$TotalbadActRule},@badActRule;
	push @{$TotalbadActDesc},@badActDesc;

	## badDesc,goodDesc
	($vars->{yingyanggongneng}{$prefix}{goodDesc},$vars->{yingyanggongneng}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);

	## goodDisease,badDisease,goodEffect,badEffect
	($vars->{yingyanggongneng}{$prefix}{goodDisease},$vars->{yingyanggongneng}{$prefix}{badDisease},$vars->{yingyanggongneng}{$prefix}{goodEffect},$vars->{yingyanggongneng}{$prefix}{badEffect})=&process_goodbad_DiseaseEffect(\@goodDisease,\@badDisease,\@goodEffect,\@badEffect);
	## badDisease,badActRule,badActDesc
	## meta 4 part's risk related disease and foodrule
	($vars->{yingyanggongneng}{$prefix}{badDisease},$vars->{yingyanggongneng}{$prefix}{badActRule},$vars->{yingyanggongneng}{$prefix}{badActDesc})=&process_badDisease_badFood(\@badDisease,\@badActRule,\@badActDesc,$badDiseaseTotal);

	return($vitaminBadN);
}

#&get_IMBALANCE_pic($vars);
sub get_Micro_Metabolism_pic {
	my ($vars)=@_;
	
	##
	$vars->{yingyanggongneng}{curpic}="-";
	$vars->{yingyanggongneng}{pic}="-";
	if ($vars->{yingyanggongneng}{desc}{badN} == 0) {
		$vars->{yingyanggongneng}{curpic}="sun_green.pdf";
		$vars->{yingyanggongneng}{pic}="gongneng_green.pdf";
	}
	elsif ($vars->{yingyanggongneng}{desc}{badN} <= 3) {
		$vars->{yingyanggongneng}{curpic}="sun_orange.pdf";
		$vars->{yingyanggongneng}{pic}="gongneng_orange.pdf";
	}
	else {
		$vars->{yingyanggongneng}{curpic}="sun_red.pdf";
		$vars->{yingyanggongneng}{pic}="gongneng_red.pdf";
	}
	
	return;
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

#($vars->{xiaohuahexishou}{desc}{$prefix}{badDisease},$vars->{xiaohuahexishou}{desc}{$prefix}{badActRule},$vars->{xiaohuahexishou}{desc}{$prefix}{badActDesc})=&process_badDisease_badFood(\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,$badDiseaseTotal);
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

#($vars->{xiaohuahexishou}{desc}{$prefix}{badDisease},$vars->{xiaohuahexishou}{desc}{$prefix}{badActRule},$vars->{xiaohuahexishou}{desc}{$prefix}{badActDesc})=&process_badDisease_badFood(\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,$badDiseaseTotal);
sub process_badDisease_badFood {
	my ($Disease,$ActRule,$ActDesc,$badDiseaseTotal)=@_;
	my ($badDisease,$badActRule,$badActDesc)=('-','-','-');

	##
	&split_uniq_hash($Disease,$badDiseaseTotal);
	my ($dis)=&split_uniq_array($Disease);
	my ($rule)=&split_uniq_array($ActRule);
	my ($desc)=&split_uniq_array($ActDesc);
	##
	if (@{$dis}) {
		if (scalar @{$dis} <= 3) {
			$badDisease=join "、",@{$dis};
		}
		else {
			$badDisease=join("、",$dis->[0],$dis->[1],$dis->[2]);
			$badDisease.="等";
		}
		$badActRule=join "、",@{$rule};
		$badActDesc=join "、",@{$desc};
	}

	$badDisease='-' if ($badDisease eq "");
	$badActRule='-' if ($badActRule eq "");
	$badActDesc='-' if ($badActDesc eq "");

	return ($badDisease,$badActRule,$badActDesc);
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
	my ($Data,$list1,$list2,$vars,$dir1,$dir2,$language,$species,$summary)=@_;
	
	##
	$vars->{gaikuang}{name}="肠道菌群概况";
	$vars->{gaikuang}{ename}="OVERVIEW";
	$vars->{gaikuang}{desc}{total}=$summary->{Total_OTU};
	$vars->{gaikuang}{desc}{classified}=$summary->{Classified_OTU};
	$vars->{gaikuang}{desc}{unclassified}=$summary->{Unclassified_OTU};
	$vars->{gaikuang}{desc}{Firmicutes_Bacteroidetes_ratio}=$summary->{Firmicutes_Bacteroidetes_ratio};

	##
	my ($benifit_risk,$harm_risk)=&item_Micro_Overvwiew('desc',$Data,$list1,$vars,$dir1,$language,$species);

	## 
	my $summary_lev = &process_risk_protect($benifit_risk,$harm_risk);
	&item_Micro_Summary('desc',$summary_lev,$list2,$vars,$dir2,$language,$species);

	## get pic
	&get_Micro_Overview_pic($vars);

	return;
}

#&get_Micro_Overview_pic($vars);
sub get_Micro_Overview_pic {
	my ($vars)=@_;
	
	##
	$vars->{gaikuang}{curpic}="-";
	$vars->{gaikuang}{pic}="-";
	if ($vars->{gaikuang}{desc}{badN} == 0) {
		$vars->{gaikuang}{curpic}="sun_green.pdf";
		$vars->{gaikuang}{pic}="gaikuang_green.pdf";
	}
	elsif ($vars->{gaikuang}{desc}{badN} <= 1) {
		$vars->{gaikuang}{curpic}="sun_orange.pdf";
		$vars->{gaikuang}{pic}="gaikuang_orange.pdf";
	}
	else {
		$vars->{gaikuang}{curpic}="sun_red.pdf";
		$vars->{gaikuang}{pic}="gaikuang_red.pdf";
	}
	
	return;
}

#my ($level)=&process_risk_protect($risk,$protect);
sub process_risk_protect {
	my ($protect,$risk)=@_;
	my $level='-';
	
	##
	$protect=&transfer_diseasefactor_level($protect);
	$risk=&transfer_diseasefactor_level($risk);
	##
	$level=&get_summary_level($protect,$risk);
	
	return ($level);
}

#$risk=&transfer_diseasefactor_level($risk);
sub transfer_diseasefactor_level {
	my ($risk)=@_;
	
	##
	if ($risk eq 'l1' || $risk eq 'l2') {
		$risk="低";
	}
	elsif ($risk eq 'l3') {
		$risk="正常";
	}
	else {
		$risk="高";
	}
	
	return ($risk);
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

#&item_Micro_Summary('fenbu',$Data,$list,$vars,$dir,$language,$species,$benifit_risk,$harm_risk);
sub item_Micro_Summary {
	my ($prefix,$lev,$list,$vars,$dir,$language,$species)=@_;
	
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get summary description
		$vars->{gaikuang}{$prefix}{summary}=$vd->{$language}{suggestion}{suggest}{$lev}{part1};
	}
	
	return;
}

#&item_Micro_Overvwiew('fenbu',$Data,$list,$vars,$dir,$language,$species,$benifit_risk,$harm_risk);
sub item_Micro_Overvwiew {
	my ($prefix,$Data,$list,$vars,$dir,$language,$species)=@_;
	
	##
	my (@good,@bad);
	my ($benifit_risk,$harm_risk);
	##
	$vars->{gaikuang}{$prefix}{abnormalN}=0;
	$vars->{gaikuang}{$prefix}{goodN}=0;
	$vars->{gaikuang}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		$vars->{gaikuang}{$prefix}{duoyangxing}=$item if ($meta =~ /junqunduoyangxing/);
		if ($meta =~ /youyijun/) {
			$vars->{gaikuang}{$prefix}{youyijun}=$item;
			$benifit_risk=$vars->{gaikuang}{$prefix}{youyijun}{level};
		}
		if ($meta =~ /youhaijun/) {
			$vars->{gaikuang}{$prefix}{youhaijun}=$item;
			$harm_risk=$vars->{gaikuang}{$prefix}{youhaijun}{level};
		}

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{gaikuang}{$prefix}{abnormalN}+=$add;
		$vars->{gaikuang}{$prefix}{goodN}+=$addG;
		$vars->{gaikuang}{$prefix}{badN}+=$addN;
	}

	## good,bad
	($vars->{gaikuang}{$prefix}{good},$vars->{gaikuang}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	
	return ($benifit_risk,$harm_risk);
}

#&Intestinal_Flora($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Flora {
	my ($Data,$list,$vars,$dir,$language,$species,$summary)=@_;
	
	##
	$vars->{changdaojunqun}{name}="肠道菌群";
	$vars->{changdaojunqun}{ename}="MICROBIOTA";
	$vars->{changdaojunqun}{explain}{majorgenus}=decode("UTF-8",$summary->{predominant_genus});
	$vars->{changdaojunqun}{explain}{majorgenus}.="菌属";
	$vars->{changdaojunqun}{explain}{majorgenus}=decode("UTF-8",$species->{en2cn}{$summary->{predominant_genus}}) if (exists $species->{en2cn}{$summary->{predominant_genus}});
	
	## 分布
	&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language,$species);
	
	return;
}

#&Intestinal_Pathogen($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Pathogen {
	my ($Data,$list,$vars,$dir,$language,$species,$summary)=@_;
	
	##
	$vars->{changdaojunqun}{pathogen}{name}="致病菌";
	$vars->{changdaojunqun}{pathogen}{ename}="PATHOGEN";
	### num,opportunistic,foodborne
	# $vars->{changdaojunqun}{pathogen}{opportunistic}=$summary->{oppor_pathogen};
	# $vars->{changdaojunqun}{pathogen}{foodborne}=$summary->{food_pathogen};
	# $vars->{changdaojunqun}{pathogen}{num}=$summary->{oppor_pathogen}+$summary->{food_pathogen};
	$vars->{changdaojunqun}{pathogen}{num}=0;
	
	## 分布
	# &item_Intestinal_Flora('pathogen',$Data,$list,$vars,$dir,$language,$species);
	&item_Intestinal_Pathogen('pathogen',$Data,$list,$vars,$dir,$language,$species);

	return;
}

#&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language,$species);
sub item_Intestinal_Flora {
	my ($prefix,$Data,$list,$vars,$dir,$language,$species)=@_;
	
	##
	my (@good,@bad,@gooditem,@baditem);
	my (@goodDesc,@badDesc,@goodDisease,@badDisease,@goodEffect,@badEffect);
	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## Filt if the bacteria genus percent is zero
		# next if ($Data->{$meta}{sample} == 0);

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
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
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			
		}
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			push @gooditem,$item;
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
		}
	}

	##
	push @{$vars->{changdaojunqun}{$prefix}{abnormalitem}},@baditem,@gooditem;
	## good,bad
	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	## goodDisease,badDisease,goodEffect,badEffect
	($vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{badEffect})=&process_goodbad_DiseaseEffect(\@goodDisease,\@badDisease,\@goodEffect,\@badEffect);
	## badDesc,goodDesc
	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);
	
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

	my $pathogen_max = 7;
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

	if (@norm_pathogen && $pathogen_number <= $pathogen_max) {
		foreach my $item (@norm_pathogen) {
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

	&undected_pathogen($prefix, $vars) if ($vars->{changdaojunqun}{pathogen}{num} == 0);

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
sub undected_pathogen {
	my ($prefix, $vars) = @_;

	my %undected_item1 = (
		"name" => "脆弱拟杆菌",
		"desc" => "可能导致菌血症、腹内感染、腹膜炎",
		"risk" => "正常",
		"level" => "l0",
	);
	my %undected_item2 = (
		"name" => "艰难梭菌",
		"desc" => "可引起艰难梭菌感染，导致腹痛腹泻等",
		"risk" => "正常",
		"level" => "l0",
	);
	my %undected_item3 = (
		"name" => "肠道沙门氏菌",
		"desc" => "可能引发肠炎，导致腹泻、发烧或腹部痉挛",
		"risk" => "正常",
		"level" => "l0",
	);
	my %undected_item4 = (
		"name" => "空肠弯曲杆菌",
		"desc" => "可导致肠胃炎等",
		"risk" => "正常",
		"level" => "l0",
	);
	push @{$vars->{changdaojunqun}{$prefix}{item}}, (\%undected_item1, \%undected_item2, \%undected_item3, \%undected_item4);
}

#my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
sub get_item {
	my ($vd,$Data,$language,$key)=@_;
	my %item=();
	
	##
	$item{value}=sprintf "%.2f",$Data->{sample}*100;
	if ($key eq "tjsjcjunqunduoyangxing") {
		$item{value}=sprintf "%.2f",$Data->{sample};
		$item{mean}=sprintf "%.2f",$Data->{ref}{mean}*100;
		$item{per100}=sprintf "%.2f",$Data->{per}*100;

		$item{per100} = "1.00" if ($Data->{per}*100 < 1);
		$item{per100} = "99.99" if ($item{per100} eq "100.00");
		## caculate dotted line position (10.2*div_value/div_pic_max);
		$item{valuelocation}=10.2*$item{value}/6;
	}
	if ($key eq "tjsjcyouyijun" || $key eq "tjsjcyouhaijun") {
		$item{per100}=sprintf "%.2f",$Data->{per}*100;

		$item{per100} = "1.00" if ($Data->{per}*100 < 1);
		$item{value} = "0.01" if ($item{value} eq "0.00");
		$item{per100} = "99.99" if ($item{per100} eq "100.00");

		$item{mean}=sprintf "%.2f",$Data->{ref}{mean}*100;
		$item{cnt1}=int($item{per100}/10);

		my $remainder=sprintf "%.2f", (($item{per100}*100)%1000)/100;
		if ($remainder == 0) {
			$item{cnt2}= (10 - $item{cnt1});
			$item{cnt3}= 0;
		}
		else {
			$item{cnt2}= (10 - $item{cnt1} - 1);
			$item{cnt3}=int($remainder+0.5);
			$item{cnt3}=9 if ($item{cnt3} == 10);
		}
	}
	$item{per}=sprintf "%.4f",$Data->{per};
	if ($Data->{per} eq "-") {
		print "$key\n\n";
	}
	$item{level}=$Data->{level};
	$item{range}=$Data->{range};
	#$item{pic}=$key."_dis.png";
	#$item{coordinate}=$Data->{per}*4.98-5.66;
	$item{name}=$vd->{$language}->{title};
	$item{desc}=$vd->{$language}->{summary}->{desc};
	if ($Data->{level} ne '-') {
		$item{risk}=$vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{effect};
		if ($vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{mechanism} ne "") {
			$item{riskdesc}=$vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{mechanism};
			$item{riskdescfull}=$vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{mechanismfull} if ($key eq "tjsjcjunqunduoyangxing");
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
			$good.="等";
		}
	}
	##
	if (@{$badA}) {
		if (scalar @{$badA} <= 5) {
			$bad=join "、",@{$badA};
		}
		else {
			$bad=join("、",$badA->[0],$badA->[1],$badA->[2],$badA->[3],$badA->[4]);
			$bad.="等";
		}
	}
	
	return ($good,$bad);
}

#&check_file($datadir,$indir,$barcode);
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
	
	## check config file
	#
	$file="$confdir/genus_latin2zh";
	$check+=&checkF($file);
	#
	$file="$confdir/tjsjc_genus_fenbu.conf";
	$check+=&checkF($file);
	#
	$file="$confdir/tjsjc_Vitamin.conf";
	$check+=&checkF($file);
	
	## check input data file
	#
	$file="$datadir/$barcode".".extract.xls";
	$check+=&checkF($file);
	#
	$file="$datadir/$barcode".".summary.xls";
	$check+=&checkF($file);
	
	## check sample information file
	#
	$file="$datadir/$barcode".".csv";
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

#my ($Vitamin)=&load_Vitamin("$indir/reference_population/config/Vitamin.conf");
sub load_Vitamin {
	my ($file)=@_;
	my %hash;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^\#/ || /^$/);
		my ($oriid,$en,$cn)=(split /\t/,$_)[0,1,2];
		$hash{$oriid}=$cn;
	}
	close IN;

	return (\%hash);
}

#my ($species)=&load_species("$indir/reference_population/config/genus_fenbu.conf","$indir/reference_population/config/genus_latin2zh");
sub load_species {
	my ($file, $genus2cn)=@_;
	my %species;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^\#/ || /^$/);
		my ($oriid,$en,$cn)=(split /\s+/,$_)[0,1,2];
		$species{meta}{$oriid}=1;
		$species{en2cn}{$en}=$cn;
	}
	close IN;

	open IN,$genus2cn or die $!;
	while (<IN>) {
		chomp;
		next if (/^\#/ || /^$/);
		my ($en,$cn)=(split /\s+/,$_)[0,1];
		$species{en2cn}{$en}=$cn;
	}
	close IN;

	return (\%species);
}

#my ($summary)=&load_summary("$datadir/$barcode.summary.xls");
sub load_summary {
	my ($file)=@_;
	my %summary;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^\#/ || /^$/);
		my ($id,$value)=(split /\s+/,$_)[0,1];
		$summary{$id}=$value;
	}
	close IN;
	
	return (\%summary);
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
				$vars->{$head[$i]}=$unit[$i] if ($head[$i] eq "barcode" || $head[$i] eq "sex" || $head[$i] eq "batch" || $head[$i] eq "receivedate" || $head[$i] eq "testdate");
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
	my ($indir,$samInfo)=@_;
	my %listH=();
	my %term2oriid=();
	
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
		&get_list($file,\%listH,\%term2oriid,"distribution");
	}
	## metabolism.list.order$sex
	$file="$indir/metabolism.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"metabolism");
	}
	## metasummary.list.order$sex
	$file="$indir/metasummary.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"metasummary");
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
	
	return (\%listH,\%term2oriid);
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
#				$listH->{$type}{$section1}{$j}{section}{$oriid}{order}=$i;
#				$listH->{$type}{$section1}{$j}{section}{$oriid}{en}=$en if (defined $en);
#				$listH->{$type}{$section1}{$j}{section}{$oriid}{cn}=$cn;
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

#my ($Data)=&calculate_meta_level_per($datadir,$indir,$barcode,$level_set);
sub calculate_meta_level_per {
	my ($datadir,$indir,$referdir,$barcode,$term2oriid,$level_set)=@_;
	my $Data=();
	
	## load sample and reference data
	($Data)=&load_sample_reference_data($datadir,$indir,$referdir,$barcode,$term2oriid);

	## get level,per,range
	&calculate_level_per($Data,$level_set);
	
	return ($Data);
}

#&calculate_level_per($Data,$level_set)
sub calculate_level_per {
	my ($Data,$level_set)=@_;
	
	##
	foreach my $meta (keys %{$Data}) {

		# get level
		($Data->{$meta}{level},$Data->{$meta}{range})=&get_level($Data->{$meta},$meta,$level_set);

		## get per
		$Data->{$meta}{per}=&get_per($Data->{$meta}, $meta);
	}
	
	return;
}

#$Data->{$meta}{per}=&get_per($Data->{$meta}, $key);
sub get_per {
	my ($meta, $key)=@_;
	my $per='-';

	if ((defined $meta->{ref}{mean} && $meta->{ref}{mean} != 0) && (defined $meta->{ref}{sd} && $meta->{ref}{sd} != 0) && exists $meta->{sample}) {
		## use R pnorm
		$R->set('v', $meta->{sample});
		$R->set('y', $meta->{ref}{mean});
		$R->set('z', $meta->{ref}{sd});
		$R->run(q`prob <- pnorm(v, y, z)`);
		$per = $R->get('prob');
	}

	## use sample data
	elsif (exists $meta->{sample} && exists $meta->{ref}) {
		my $N=scalar @{$meta->{ref}{values}};
		my $n=0;
		foreach my $v (@{$meta->{ref}{values}}) {
			if ($meta->{sample} < $v) {
				last;
			}
			$n++;
		}
		$per=$n/$N if ($N != 0);
	}

	return ($per);
}

#($Data->{$meta}{level},$Data->{$meta}{range})=&get_level($Data->{$meta},$meta,$level_set);
sub get_level {
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

#($Data)=&load_sample_reference_data($datadir,$indir,$barcode,$term2oriid);
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
                          default "$Bin/_template.db2xml";
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
