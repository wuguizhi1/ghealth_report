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
use Math::SymbolicX::Statistics::Distributions ':all';
use Statistics::Distributions qw(uprob);
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my $version = "1.0.0";
my $R = Statistics::R->new();

####
my ($datadir,$indir,$referdir,$dirout,$language,$barcode,$level_set,$historydir,$confdir,$help);
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

my %vars=();
## check file and dir
&check_file($datadir,$indir,$referdir,$confdir,$barcode);

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
);

## load sample information file
my ($samInfo)=&load_sample_info($barcode,$datadir,\%vars);

## load $type.list.order$sex
my ($listH,$term2oriid)=&load_list($indir,$samInfo);

## process relative abundance
&process_relative_abundance($datadir,$indir,$referdir,$confdir,\%vars);

## calculate meta level and per by sample data and ref data
my ($Data)=&calculate_meta_level_per($datadir,$indir,$referdir,$term2oriid,$barcode,$level_set);

## main process
&main_process($Data,$listH,\%{$diseaseOrder{question}},\%{$diseaseOrder{predict}},\%vars,$datadir,$indir,$language,$barcode,$samInfo,$historydir,$confdir);

## sample result csv for erp import
&sample_result_csv(\%vars,$listH->{pathogens}{'致病菌'},$barcode,$dirout);

############################# sub ###############################################

#&main_process($Data,$listH,\%vars,$datadir,$indir,$language,$barcode,$samInfo);
sub main_process {
	my ($Data,$listH,$listQ,$listP,$vars,$datadir,$indir,$language,$barcode,$samInfo,$historydir,$confdir)=@_;
	
	## 消化和吸收
	my ($Vitamin)=&load_Vitamin("$confdir/zx_Vitamin.conf");
	my ($regulateLow,$goodDiseaseTotal,$badDiseaseTotal)=&Absorption_and_Metabolism($Data,$listH,$vars,"$indir/XML",$language,$Vitamin);

	## 炎症和免疫
	&Inflammation_and_Immune($Data,$listH->{giimmune}{'炎症和免疫'},$vars,"$indir/XML/giimmune/$language",$language,$regulateLow);

	## 肠道菌群
	# 所有菌属
	my ($species)=&load_species("$confdir/zx_genus_fenbu.conf", "$confdir/genus_latin2zh");
	# summary信息
	my ($summary)=&load_summary("$datadir/$barcode.summary.xls");
	# 分布+致病菌
	my ($goodDiseaseFlora,$badDiseaseFlora)=&Intestinal_Flora($Data,$listH,$vars,"$indir/XML",$language,$species,$summary);
	# 相关疾病风险
	&Intestinal_Flora_disease($Data,$vars,$listH,"$indir/XML",$language);

	## 粪便状态
	&Stool_Status($listH->{stoolList}{"粪便状态"},$vars,"$indir/XML/stool/$language",$language,$samInfo);

	## 膳食方案+肠道调节方案+运动方案
	&Program($vars,"$indir/XML",$language,$samInfo,$listH->{diseaseList},$listQ,$listP,$regulateLow,$badDiseaseTotal);

	## 多次检测结果比较
	# $samInfo->{batch}; $samInfo->{'历次检测编号'}
	&Compare_history($vars,$samInfo,$historydir) if ($samInfo->{batch} > 1);

	## 报告日期
	$vars->{reportdata}=strftime("%Y-%m-%d",localtime());
	
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
			$xml->{xiaohuahexishou}{curpic} =~ /sun-(.*)\.pdf/;
			$vars->{$batch_str}{xiaohuahexishou}{absorption}{pic}=$1."circle.pdf";
			$vars->{$batch_str}{xiaohuahexishou}{absorption}{piccolor}=$1;
			$xml->{xiaohuahexishou}{metabolism}{curpic} =~ /sun-(.*)\.pdf/;
			$vars->{$batch_str}{xiaohuahexishou}{metabolism}{pic}=$1."circle.pdf";
			$vars->{$batch_str}{xiaohuahexishou}{metabolism}{piccolor}=$1;
			$xml->{yanzhenghemianyi}{curpic} =~ /sun-(.*)\.pdf/;
			$vars->{$batch_str}{yanzhenghemianyi}{pic}=$1."circle.pdf";
			$vars->{$batch_str}{yanzhenghemianyi}{piccolor}=$1;
			$xml->{changdaojunqun}{pathogens}{curpic} =~ /sun-(.*)\.pdf/;
			$vars->{$batch_str}{changdaojunqun}{pathogens}{pic}=$1."circle.pdf";
			$vars->{$batch_str}{changdaojunqun}{pathogens}{piccolor}=$1;
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

	$xml->{xiaohuahexishou}{curpic} =~ /sun-(.*)\.pdf/;
	$vars->{$batch_str}{xiaohuahexishou}{absorption}{pic}=$1."circle.pdf";
	$vars->{$batch_str}{xiaohuahexishou}{absorption}{piccolor}=$1;
	$xml->{xiaohuahexishou}{metabolism}{curpic} =~ /sun-(.*)\.pdf/;
	$vars->{$batch_str}{xiaohuahexishou}{metabolism}{pic}=$1."circle.pdf";
	$vars->{$batch_str}{xiaohuahexishou}{metabolism}{piccolor}=$1;
	$xml->{yanzhenghemianyi}{curpic} =~ /sun-(.*)\.pdf/;
	$vars->{$batch_str}{yanzhenghemianyi}{pic}=$1."circle.pdf";
	$vars->{$batch_str}{yanzhenghemianyi}{piccolor}=$1;
	$xml->{changdaojunqun}{pathogens}{curpic} =~ /sun-(.*)\.pdf/;
	$vars->{$batch_str}{changdaojunqun}{pathogens}{pic}=$1."circle.pdf";
	$vars->{$batch_str}{changdaojunqun}{pathogens}{piccolor}=$1;

	#$vars->{$batch_str}{gaikuang}{desc}=\%{$xml->{gaikuang}{desc}};
	&compare_gaikuang($vars, "duoyangxing", $vars->{changdaojunqun}{gaikuang}{diversity}, $xml->{changdaojunqun}{gaikuang}{diversity});
	#&compare_gaikuang($vars, "youyijun", $vars->{gaikuang}{desc}{youyijun}{per}, $xml->{gaikuang}{desc}{youyijun}{per});
	#&compare_gaikuang($vars, "youhaijun", $vars->{gaikuang}{desc}{youhaijun}{per}, $xml->{gaikuang}{desc}{youhaijun}{per});

	$vars->{$batch_str}{changdaojunqun}{gaikuang}{goodN}=$xml->{changdaojunqun}{gaikuang}{goodN};
	$vars->{$batch_str}{changdaojunqun}{gaikuang}{badN}=$xml->{changdaojunqun}{gaikuang}{badN};
	$vars->{$batch_str}{changdaojunqun}{gaikuang}{good}=$xml->{changdaojunqun}{gaikuang}{good};
	$vars->{$batch_str}{changdaojunqun}{gaikuang}{bad}=$xml->{changdaojunqun}{gaikuang}{bad};

	$vars->{$batch_str}{changdaojunqun}{fenbu}{goodN}=$xml->{changdaojunqun}{fenbu}{goodN};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{badN}=$xml->{changdaojunqun}{fenbu}{badN};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{good}=$xml->{changdaojunqun}{fenbu}{good};
	$vars->{$batch_str}{changdaojunqun}{fenbu}{bad}=$xml->{changdaojunqun}{fenbu}{bad};

	$vars->{$batch_str}{changdaojunqun}{pathogens}{num}=$xml->{changdaojunqun}{pathogens}{num};
	$vars->{$batch_str}{changdaojunqun}{pathogens}{badN}=$xml->{changdaojunqun}{pathogens}{badN};
	$vars->{$batch_str}{changdaojunqun}{pathogens}{bad}=$xml->{changdaojunqun}{pathogens}{bad};

	$vars->{$batch_str}{xiaohuahexishou}{explain}=\%{$xml->{xiaohuahexishou}{explain}};
	$vars->{$batch_str}{xiaohuahexishou}{metabolism}{explain}=\%{$xml->{xiaohuahexishou}{metabolism}{explain}};

	$vars->{$batch_str}{yanzhenghemianyi}{giimmune}{explain}=\%{$xml->{yanzhenghemianyi}{giimmune}{explain}};

}

#&compare_gaikuang($vars, "duoyangxing", $vars->{gaikuang}{desc}{duoyangxing}{per}, $xml->{gaikuang}{desc}{duoyangxing}{per});
sub compare_gaikuang {
	my ($vars, $key, $valnow, $vallast) = @_;

	if (($valnow - $vallast) >= 0.05) {
		$vars->{comparelast}{changdaojunqun}{gaikuang}{$key} = "升高";
	}
	elsif (($valnow - $vallast) <= -0.05) {
		$vars->{comparelast}{changdaojunqun}{gaikuang}{$key} = "降低";
	}
	else {
		$vars->{comparelast}{changdaojunqun}{gaikuang}{$key} = "保持不变";
	}
}

#&sample_result_csv(\%vars,$listH,$barcode,$dirout);
sub sample_result_csv {
	my ($vars,$list,$barcode,$dirout) = @_;


	my $outfile="$dirout/$barcode".".result.csv";
	open my $out, '>:encoding(utf-8)', $outfile || die $!;

	print $out "检测者信息\n检测编号\t姓名\n$vars->{barcode}\t$vars->{name}\n";

	my $term_str = "";
	my $result_str = "";

	my $level;
	print $out "菌群概况\n";
	foreach my $item (@{$vars{changdaojunqun}{gaikuang}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	print $out "菌群分布\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{changdaojunqun}{fenbu}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	print $out "有机小分子物质代谢\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{xiaohuahexishou}{metabolism}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	my %pathogens_result=();
	foreach my $num (sort keys %{$list}) {
		foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}){
			$pathogens_result{$list->{$num}{section}{$meta}{cn}} = "未检出(正常)";
		}
	}
	foreach my $item (@{$vars{changdaojunqun}{pathogens}{items}{item}}) {
		$level = &transfer_pathogen_level($item->{level});
		$pathogens_result{$item->{name}} = "$level($item->{risk})";
	}
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

	print $out "三大营养物质代谢\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{xiaohuahexishou}{zongping}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	foreach my $item (@{$vars{xiaohuahexishou}{protein}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	foreach my $item (@{$vars{xiaohuahexishou}{fat}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	foreach my $item (@{$vars{xiaohuahexishou}{carbon}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	print $out "炎症和免疫\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{yanzhenghemianyi}{giimmune}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	print $out "相关疾病风险\n";
	$term_str = "";
	$result_str = "";
	foreach my $disease (@{$vars{changdaojunqun}{disease}}) {
		foreach my $item (@{$disease->{items}{item}}) {
			$term_str .= "$disease->{name}"."$item->{name}\t";
			$level = &transfer_item_level($item->{level});
			$result_str .="$level($item->{risk})\t";
		}
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	print $out "粪便状态\n";
	$term_str = "";
	$result_str = "";
	foreach my $item (@{$vars{fenbianzhuangtai}{desc}{items}{item}}) {
		$term_str .= "$item->{name}\t";
		$result_str .="$item->{sample}{desc}\t";
	}
	$term_str =~s/\s+$//;
	$result_str =~s/\s+$//;
	print $out "$term_str\n$result_str\n";

	close $out;

	return;
}


#&Program($vars,"$indir/XML",$language,$samInfo,$listH->{diseaseList},$regulateLow,$badDiseaseTotal);
sub Program {
	my ($vars,$dir,$language,$samInfo,$list,$listQ,$listP,$regulateLow,$badDiseaseTotal)=@_;
	my ($badDiseaseF,$badDiseaseR,$diseaseF,$diseaseR,$ageF,$ageR,$regulateC)=('-','-','-','-','-','-','-');

	## get disease
	($diseaseF,$diseaseR)=&get_program_disease($samInfo->{disease},$listQ,$listP);
	## get risk bad disease
	($badDiseaseF,$badDiseaseR)=&get_bad_disease_meta($vars->{changdaojunqun}{diseaseHigh},$badDiseaseTotal,$listP);
	## get regulate condition key
	($regulateC)=&get_regulate_condition($vars);
	## if no food disease, get age
#	if ($diseaseF eq '-' && $badDiseaseF eq '-') {
#		($ageF)=&get_program_age($samInfo->{age});
#	}
#	elsif ($diseaseF eq '-') {
#		$badDiseaseF=$list->{$badDiseaseF};
#	}
#	else {
#		$diseaseF=$list->{$diseaseF};
#	}
#	## if no regulate disease, get age
#	if ($diseaseR eq '-' && $badDiseaseR eq '-') {
#		($ageR)=&get_program_age($samInfo->{age});
#	}
#	elsif ($diseaseR eq '-') {
#		$badDiseaseR=$list->{$badDiseaseR};
#	}
#	else {
#		$diseaseR=$list->{$diseaseR};
#	}

	## give all food and regulate keys to adjust orders of regulation schemes
	$ageF=&get_program_age($samInfo->{age});
	if ($badDiseaseF ne '-') {
		$badDiseaseF=$list->{$badDiseaseF};
	}
	if ($diseaseF ne '-') {
		$diseaseF=$list->{$diseaseF};
	}
	## if no regulate disease, get age
	$ageR=&get_program_age($samInfo->{age});
	if ($badDiseaseR ne '-') {
		$badDiseaseR=$list->{$badDiseaseR};
	}
	if ($diseaseR ne '-') {
		$diseaseR=$list->{$diseaseR};
	}

	## process food and sport
	&process_food_sport($vars,$dir,$diseaseF,$ageF,$badDiseaseF,$regulateC,$language);
	## process regulate
	&process_regulate($vars,$dir,$regulateLow,$diseaseR,$ageR,$badDiseaseR,$regulateC,$language);
	
	return;
}

#&get_regulate_condition($vars);
sub get_regulate_condition {
	my ($vars)=@_;
	my $regulateC='-';

	if ($vars->{changdaojunqun}{gaikuang}{divlevel} eq 'l1') {
		$regulateC="condition1";
	}
	elsif ($vars->{changdaojunqun}{gaikuang}{benifilevel} eq 'l1') {
		$regulateC="condition2";
	}
	elsif ($vars->{changdaojunqun}{gaikuang}{harmlevel} eq 'l5') {
		$regulateC="condition3";
	}
	elsif ($vars->{changdaojunqun}{pathogens}{badN} >= 2) {
		$regulateC="condition4";
	}
	elsif ($vars->{changdaojunqun}{gaikuang}{divlevel} eq 'l2') {
		$regulateC="condition5";
	}
	elsif ($vars->{changdaojunqun}{gaikuang}{benifilevel} eq 'l2') {
		$regulateC="condition6";
	}
	elsif ($vars->{changdaojunqun}{gaikuang}{harmlevel} eq 'l4') {
		$regulateC="condition7";
	}
	return($regulateC);
}

#&process_regulate($vars,$die,$regulateLow,$diseaseR,$ageR,$badDiseaseR,$language);
sub process_regulate {
	my ($vars,$dir,$regulate,$disease,$age,$badDisease,$regulateC,$language)=@_;

	#print "$disease\n$age\n$badDisease\n$regulateC\n";die;
	## 补充益生元和益生菌
	#if ($regulateC ne '-' && $disease ne 'zxchangyan') {
	if ($regulateC ne '-') {
		my $vd=XMLin("$dir/dysbiosis/$language/zxjunqunshihengtiaojie.xml",NoAttr=>1,SuppressEmpty => "");
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{metaregulate}{$regulateC}{$age}{title}, desc => $vd->{$language}{metaregulate}{$regulateC}{$age}{desc}, spotnum => 0};
		#if ($disease eq '-' && $badDisease eq '-') {
		#	push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{metaregulate}{$regulateC}{$age}{nutrititle}, desc => $vd->{$language}{metaregulate}{$regulateC}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{metaregulate}{$regulateC}{$age}{nutridesc} ne "");
		#}
		#elsif ($disease ne '-') {
		#	my $vd=XMLin("$dir/disease/$language/$disease.xml",NoAttr=>1,SuppressEmpty => "");
		#	push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{regulate}{disease}{$age}{nutrititle}, desc => $vd->{$language}{regulate}{disease}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{regulate}{disease}{$age}{nutridesc} ne "");
		#}
		#elsif ($badDisease ne '-') {
		#	my $vd=XMLin("$dir/disease/$language/$badDisease.xml",NoAttr=>1,SuppressEmpty => "");
		#	push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{regulate}{highrisk}{$age}{nutrititle}, desc => $vd->{$language}{regulate}{highrisk}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{regulate}{highrisk}{$age}{nutridesc} ne "");
		#}
		#elsif ($age ne '-') {
		#	my $vd=XMLin("$dir/age/$language/zxnianling.xml",NoAttr=>1,SuppressEmpty => "");
		#	push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{metaregulate}{desc}{$age}{nutrititle}, desc => $vd->{$language}{metaregulate}{desc}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{metaregulate}{desc}{$age}{nutridesc} ne "");
		#}
	}
	elsif ($disease ne '-') {
		my $vd=XMLin("$dir/disease/$language/$disease.xml",NoAttr=>1,SuppressEmpty => "");
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{regulate}{disease}{$age}{title}, desc => $vd->{$language}{regulate}{disease}{$age}{desc}, spotnum => 0};
		#push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{regulate}{disease}{$age}{nutrititle}, desc => $vd->{$language}{regulate}{disease}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{regulate}{disease}{$age}{nutridesc} ne "");
	}
	elsif ($badDisease ne '-') {
		my $vd=XMLin("$dir/disease/$language/$badDisease.xml",NoAttr=>1,SuppressEmpty => "");
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{regulate}{highrisk}{$age}{title}, desc => $vd->{$language}{regulate}{highrisk}{$age}{desc}, spotnum => 0};
		#push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{regulate}{highrisk}{$age}{nutrititle}, desc => $vd->{$language}{regulate}{highrisk}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{regulate}{highrisk}{$age}{nutridesc} ne "");
	}
	elsif ($age ne '-') {
		my $vd=XMLin("$dir/age/$language/zxnianling.xml",NoAttr=>1,SuppressEmpty => "");
#		my @lines=split /\n/,$vd->{$language}->{metaregulate}->{$age}->{desc};
#		#
#		my $spotnum=0;
#		my ($desc,@spot);
#		foreach my $line (@lines) {
#			next if ($line eq "");
#			if ($line =~ /益生元：/ || $line =~ /益生菌：/) {
#				push @spot,$line;
#				$spotnum++;
#			}
#			else {
#				$desc.=$line;
#			}
#		}
#		push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}->{metaregulate}->{$age}->{title}, desc => $desc, spotnum => $spotnum, spot => \@spot};
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{metaregulate}{desc}{$age}{title}, desc => $vd->{$language}{metaregulate}{desc}{$age}{desc}, spotnum => 0};
		#push @{$vars->{fangan}{regulate}{desc}{item}},{name => $vd->{$language}{metaregulate}{desc}{$age}{nutrititle}, desc => $vd->{$language}{metaregulate}{desc}{$age}{nutridesc}, spotnum => 0} if ($vd->{$language}{metaregulate}{desc}{$age}{nutridesc} ne "");
	}

	## 其他
	if ($vars->{xiaohuahexishou}->{metabolism}->{explain}->{fumeiQbad} eq "Y") {
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => "补充辅酶Q", desc => "可以帮助您增强抗氧化能力。",spotnum => 0};
	}
	if ($vars->{xiaohuahexishou}->{metabolism}->{explain}->{vitaminbadN} > 2) {
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => "补充复合维生素", desc => "均衡多种维生素营养，促进碳水化合物、蛋白质、脂类的正常代谢。", spotnum => 0};
	}
	elsif ($vars->{xiaohuahexishou}->{metabolism}->{explain}->{vitaminbadN} > 0) {
		push @{$vars->{fangan}{regulate}{desc}{item}},{name => "补充$vars->{xiaohuahexishou}->{metabolism}->{explain}->{vitaminBad}", desc => "均衡多种维生素营养，促进碳水化合物、蛋白质、脂类的正常代谢。", spotnum => 0};
	}

	return;
}

#&process_food_sport($vars,$dir,$diseaseF,$ageF,$badDiseaseF,$language);
sub process_food_sport {
	my ($vars,$dir,$disease,$age,$badDisease,$regulateC,$language)=@_;

	##
	if ($disease ne '-') {
		my $vd=XMLin("$dir/disease/$language/$disease.xml",NoAttr=>1,SuppressEmpty => "");
		&get_food($vd->{$language}{food}{disease}{$age},$vars);
		&get_sport($vd->{$language}{sport}{disease}{$age},$vars);
	}
	elsif ($regulateC ne '-') {
		my $vd=XMLin("$dir/dysbiosis/$language/zxjunqunshihengtiaojie.xml",NoAttr=>1,SuppressEmpty => "");
		&get_food($vd->{$language}{metafood}{$regulateC}{$age},$vars);
		&get_sport($vd->{$language}{metasport}{$regulateC}{$age},$vars);
	}
	elsif ($badDisease ne '-') {
		my $vd=XMLin("$dir/disease/$language/$badDisease.xml",NoAttr=>1,SuppressEmpty => "");
		&get_food($vd->{$language}{food}{highrisk}{$age},$vars);
		&get_sport($vd->{$language}{sport}{highrisk}{$age},$vars);
	}
	elsif ($age ne '-') {
		my $vd=XMLin("$dir/age/$language/zxnianling.xml",NoAttr=>1,SuppressEmpty => "");
		&get_food($vd->{$language}{metafood}{desc}{$age},$vars);
		&get_sport($vd->{$language}{metasport}{desc}{$age},$vars);
	}
	
	return;
}

#&Intestinal_Flora_disease($Data,$vars,$listH,"$indir/XML",$language);
sub Intestinal_Flora_disease {
	my ($Data,$vars,$list,$dir,$language)=@_;
	
	##
	my (@diseaseName,@diseaseHigh,@diseaseLow,@diseaseNormal);
	my (@actrule,@actmore,@actless,@foodrule,@foodmore,@foodless);
	$vars->{changdaojunqun}{diseaseNum}=0;
	$vars->{changdaojunqun}{diseaseName}="-";
	$vars->{changdaojunqun}{diseaseHighN}=0;
	$vars->{changdaojunqun}{diseaseHigh}="-";
	$vars->{changdaojunqun}{diseaseHighRule}="-";
	$vars->{changdaojunqun}{diseaseHighMore}="-";
	$vars->{changdaojunqun}{diseaseHighLess}="-";
	$vars->{changdaojunqun}{diseaseLowN}=0;
	$vars->{changdaojunqun}{diseaseLow}="-";
	$vars->{changdaojunqun}{diseaseNormalN}=0;
	$vars->{changdaojunqun}{diseaseNormal}="-";
	##
	foreach my $num (sort {$a <=> $b} keys %{$list->{diseasefactor}{'相关疾病风险'}}) {
		my %disease=();
		##
		my ($risk,$protect);
		$vars->{changdaojunqun}{diseaseNum}++;
		$disease{name}=$list->{diseasefactor}{'相关疾病风险'}{$num}{name};
		push @diseaseName,$disease{name};
		my $disXML=XMLin("$dir/disease/$language/$list->{diseaseList}{$disease{name}}.xml",NoAttr=>1,SuppressEmpty => "");
		##
		foreach my $meta (sort {$list->{diseasefactor}{'相关疾病风险'}{$num}{section}{$a}{order} <=> $list->{diseasefactor}{'相关疾病风险'}{$num}{section}{$b}{order}} keys %{$list->{diseasefactor}{'相关疾病风险'}{$num}{section}}) {
			next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});
			
			### get $vd
			my $vd=XMLin("$dir/diseasefactor/$language/$meta.xml",NoAttr=>1,SuppressEmpty => "");
			
			## get one item
			my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
			&rename_disease_factor_name($item);
			push @{$disease{items}{item}},$item;
			
			## get $risk|$protect
			if ($item->{name} eq "危险因子") {
				$risk=$item->{level};
			}
			elsif ($item->{name} eq "保护因子") {
				$protect=$item->{level};
			}
		}
		## get disease risk
		my ($level)=&process_risk_protect($risk,$protect);
		$disease{summary}{line1}=$disXML->{$language}->{suggestion}->{suggest}->{$level}->{part1};
		$disease{summary}{line2}=$disXML->{$language}->{suggestion}->{suggest}->{$level}->{part2};
		$disease{summary}{line3}=$disXML->{$language}->{suggestion}->{suggest}->{$level}->{part3};
		&get_disease_high_middle_low($vars,$disease{name},$disXML->{$language},$level,\@diseaseHigh,\@diseaseLow,\@diseaseNormal,\@actrule,\@actmore,\@actless,\@foodrule,\@foodmore,\@foodless);
	
		##
		push @{$vars->{changdaojunqun}{disease}},\%disease;
	}
	
	##
	$vars->{changdaojunqun}{diseaseName}=join "、", @diseaseName;
	$vars->{changdaojunqun}{diseaseNormal}=join "、", @diseaseNormal if (@diseaseNormal);
	$vars->{changdaojunqun}{diseaseHigh}=join "、", @diseaseHigh if (@diseaseHigh);
	$vars->{changdaojunqun}{diseaseLow}=join "、", @diseaseLow if (@diseaseLow);
	$vars->{changdaojunqun}{diseaseHighActRule}=&process_HighRisk_disease_food(\@actrule);
	$vars->{changdaojunqun}{diseaseHighActMore}=&process_HighRisk_disease_food(\@actmore);
	$vars->{changdaojunqun}{diseaseHighActLess}=&process_HighRisk_disease_food(\@actless);
	$vars->{changdaojunqun}{diseaseHighFoodRule}=&process_HighRisk_disease_food(\@foodrule);
	$vars->{changdaojunqun}{diseaseHighFoodMore}=&process_HighRisk_disease_food(\@foodmore);
	$vars->{changdaojunqun}{diseaseHighFoodLess}=&process_HighRisk_disease_food(\@foodless);
	
	return;
}


#&Intestinal_Flora($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Flora {
	my ($Data,$list,$vars,$dir,$language,$species,$summary)=@_;
	my %goodDiseaseFlora=();
	my %badDiseaseFlora=();
	##
	my (@gooditem,@baditem);											    ## gaikuang|fenbu|pathogens good and bad effect items
	my (@goodDesc,@badDesc);											    ## gaikuang|fenbu|pathogens good and bad items mechanism
	my (@goodDisease,@badDisease);										    ## gaikuang|fenbu|pathogens good and bad items disease related
	my (@goodEffect,@badEffect,@badActRule,@badActDesc);				    ## gaikuang|fenbu|pathogens good and bad items disease related

	$vars->{changdaojunqun}{name}="肠道菌群";
	$vars->{changdaojunqun}{ename}="MICROBIOTA";
	$vars->{changdaojunqun}{abnormal}{num}=0;
	$vars->{changdaojunqun}{explain}{majorgenus}=decode("UTF-8",$summary->{predominant_genus});
	$vars->{changdaojunqun}{explain}{majorgenus}.="菌属";
	$vars->{changdaojunqun}{explain}{majorgenus}=decode("UTF-8",$species->{en2cn}{$summary->{predominant_genus}}) if (exists $species->{en2cn}{$summary->{predominant_genus}});

	## 概况
	&item_Intestinal_Overview('1','gaikuang',$Data,$list->{overview}{"菌群概况"},$vars,"$dir/overview/$language",$language,$species,$summary,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@goodEffect,\@badEffect,\@badActRule,\@badActDesc);

	## 分布
	&item_Intestinal_Flora('1','fenbu',$Data,$list->{distribution}{"菌群分布"},$vars,"$dir/distribution/$language",$language,$species,$summary,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@goodEffect,\@badEffect,\@badActRule,\@badActDesc);

	## 致病菌
	&item_Intestinal_Pathogen('pathogens',$Data,$list->{pathogens}{"致病菌"},$vars,"$dir/pathogens/$language",$language,$species,$summary,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@goodEffect,\@badEffect,\@badActRule,\@badActDesc);

	## good,bad
	($vars->{changdaojunqun}{explain}{good},$vars->{changdaojunqun}{explain}{bad})=&process_good_bad(\@gooditem,\@baditem);
	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActDesc,badDisease,goodDisease
	## get flora overview risk related effect/disease and actrule
	($vars->{changdaojunqun}{explain}{goodDesc},$vars->{changdaojunqun}{explain}{goodEffect},$vars->{changdaojunqun}{explain}{goodDisease},$vars->{changdaojunqun}{explain}{badDesc},$vars->{changdaojunqun}{explain}{badEffect},$vars->{changdaojunqun}{explain}{badDisease},$vars->{changdaojunqun}{explain}{badActRule},$vars->{changdaojunqun}{explain}{badActDesc})=&process_Flora_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActDesc,\%goodDiseaseFlora,\%badDiseaseFlora);

	return(\%goodDiseaseFlora,\%badDiseaseFlora);
}

#&item_Intestinal_Pathogen('pathogens',$Data,$list,$vars,$dir,$language,$species,$summary,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@goodEffect,\@badEffect,\@badActRule,\@badActDesc);
sub item_Intestinal_Pathogen {
	my ($prefix,$Data,$list,$vars,$dir,$language,$species,$summary,$Floragood,$Florabad,$FloragoodDesc,$FlorabadDesc,$FloragoodDisease,$FlorabadDisease,$FloragoodEffect,$FlorabadEffect,$FlorabadActRule,$FlorabadActDesc)=@_;
	my %goodDiseasePatho=();
	my %badDiseasePatho=();

	##
	my (@good,@bad);
	my (@badDesc,@badEffect,@badDisease,@badActRule,@badActDesc);
	my (@goodDesc,@goodEffect,@goodDisease);

	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;
	$vars->{changdaojunqun}{$prefix}{num}=0;

	##
	my @patho_names=();
	my @norm_pathogen=();
	my @bad_pathogen=();
	my %pathogen2oriid=();
	foreach my $num (sort keys %{$list}) {
		foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}){
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
	}

	my $pathogen_max = 7;
	my $pathogen_number = 1;
	if (@bad_pathogen) {
		foreach my $item (@bad_pathogen) {
			if ($pathogen_number <= $pathogen_max) {
				push @patho_names, $item->{name};
				$vars->{changdaojunqun}{$prefix}{num} += 1;
				push @{$vars->{changdaojunqun}{$prefix}{items}{item}},$item;
				my $vd=XMLin("$dir/$pathogen2oriid{$item->{name}}.xml",NoAttr=>1,SuppressEmpty => "");
				## get abnormalN,badN,good,bad
				my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
				$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
				$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
				$vars->{changdaojunqun}{$prefix}{badN}+=$addN;

				###
				my %new_item = %{$item};
				$new_item{index} = $vars->{changdaojunqun}{$prefix}{num};
				push @{$vars->{changdaojunqun}{$prefix}{baditems}{items}{item}},\%new_item;
				## Get risk disease & badDesc info
				push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
				push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
				push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
				push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
				push @badActDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
				$pathogen_number ++;
			}
		}
	}

	if (@norm_pathogen && $pathogen_number <= $pathogen_max) {
		foreach my $item (@norm_pathogen) {
			if ($pathogen_number <= $pathogen_max) {
				push @patho_names, $item->{name};
				$vars->{changdaojunqun}{$prefix}{num} += 1;
				push @{$vars->{changdaojunqun}{$prefix}{items}{item}},$item;
				my $vd=XMLin("$dir/$pathogen2oriid{$item->{name}}.xml",NoAttr=>1,SuppressEmpty => "");
				## get abnormalN,badN,good,bad
				my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
				$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
				$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
				$vars->{changdaojunqun}{$prefix}{badN}+=$addN;

				###
				my %new_item = %{$item};
				$new_item{index} = $vars->{changdaojunqun}{$prefix}{num};
				push @{$vars->{changdaojunqun}{$prefix}{baditems}{items}{item}},\%new_item;
				## Get risk disease & badDesc info
				push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
				push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
				push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
				push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
				push @badActDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
				$pathogen_number ++;
			}
		}
	}

	&undected_pathogen($prefix, $vars) if ($vars->{changdaojunqun}{$prefix}{num} == 0);

	($vars->{changdaojunqun}{$prefix}{nameDesc})=&process_pathogen_names(\@patho_names);

	## good,bad
	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	push @{$Floragood},@good;
	push @{$Florabad},@bad;
	push @{$FloragoodDesc},@goodDesc;
	push @{$FlorabadDesc},@badDesc;
	push @{$FloragoodEffect},@goodEffect;
	push @{$FloragoodDisease},@goodDisease;
	push @{$FlorabadEffect},@badEffect;
	push @{$FlorabadDisease},@badDisease;
	push @{$FlorabadActRule},@badActRule;
	push @{$FlorabadActDesc},@badActDesc;

	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActDesc,badDisease,goodDisease
	## get flora overview risk related effect/disease and actrule
	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDesc},$vars->{changdaojunqun}{$prefix}{badEffect},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{badActRule},$vars->{changdaojunqun}{$prefix}{badActDesc})=&process_Flora_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActDesc,\%goodDiseasePatho,\%badDiseasePatho);

	$vars->{changdaojunqun}{$prefix}{name}="致病菌";
	$vars->{changdaojunqun}{$prefix}{ename}="INFECTION";
	&get_INFECTION_pic($vars);

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
		"pic" => "zxcuiruoniganjun_dis.png",
		"per" => "0",
		"coordinate" => "-5.66",
	);
	my %undected_item2 = (
		"name" => "艰难梭菌",
		"desc" => "可引起艰难梭菌感染，导致腹痛腹泻等",
		"risk" => "正常",
		"level" => "l0",
		"pic" => "zxjiannansuojun_dis.png",
		"per" => "0",
		"coordinate" => "-5.66",
	);
	my %undected_item3 = (
		"name" => "肠道沙门氏菌",
		"desc" => "可能引发肠炎，导致腹泻、发烧或腹部痉挛",
		"risk" => "正常",
		"level" => "l0",
		"pic" => "zxchangdaoshamenshijun_dis.png",
		"per" => "0",
		"coordinate" => "-5.66",
	);
	my %undected_item4 = (
		"name" => "空肠弯曲杆菌",
		"desc" => "可导致肠胃炎等",
		"risk" => "正常",
		"level" => "l0",
		"pic" => "zxkongchangwanquganjun_dis.png",
		"per" => "0",
		"coordinate" => "-5.66",
	);
	push @{$vars->{changdaojunqun}{$prefix}{items}{item}}, (\%undected_item1, \%undected_item2, \%undected_item3, \%undected_item4);
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

#&get_INFECTION_pic($vars);
sub get_INFECTION_pic {
	my ($vars)=@_;
	
	##
	$vars->{changdaojunqun}{pathogens}{curpic}="-";
	$vars->{changdaojunqun}{pathogens}{pic}="-";
	if ($vars->{changdaojunqun}{pathogens}{num} == 0) {
		$vars->{changdaojunqun}{pathogens}{curpic}="sun-green.pdf";
		$vars->{changdaojunqun}{pathogens}{pic}="INFECTION-green.pdf";
	}
	elsif ($vars->{changdaojunqun}{pathogens}{num} <= 3 && $vars->{changdaojunqun}{pathogens}{badN} == 0) {
		$vars->{changdaojunqun}{pathogens}{curpic}="sun-yellow.pdf";
		$vars->{changdaojunqun}{pathogens}{pic}="INFECTION-yellow.pdf";
	}
	else {
		$vars->{changdaojunqun}{pathogens}{curpic}="sun-red.pdf";
		$vars->{changdaojunqun}{pathogens}{pic}="INFECTION-red.pdf";
	}
	
	return;
}

#&item_Intestinal_Flora('1','fenbu',$Data,$list,$vars,$dir,$language,$species,$summary,\@gooditem,\@baditem,\@goodDisease,\@badDisease,\@goodEffect,\@badEffect,\@badActRule,\@badActDesc);
sub item_Intestinal_Flora {
	my ($num,$prefix,$Data,$list,$vars,$dir,$language,$species,$summary,$Floragood,$Florabad,$FloragoodDesc,$FlorabadDesc,$FloragoodDisease,$FlorabadDisease,$FloragoodEffect,$FlorabadEffect,$FlorabadActRule,$FlorabadActDesc)=@_;
	my %goodDiseaseFenbu=();
	my %badDiseaseFenbu=();

	##
	my (@good,@bad);
	my (@badDesc,@badEffect,@badDisease,@badActRule,@badActDesc);
	my (@goodDesc,@goodEffect,@goodDisease);

	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;

	##
	foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		push @{$vars->{changdaojunqun}{$prefix}{items}{item}},$item;

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
		$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
		$vars->{changdaojunqun}{$prefix}{badN}+=$addN;

		## get abnormal
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{changdaojunqun}{$prefix}{badN};
			push @{$vars->{changdaojunqun}{$prefix}{baditems}{items}{item}},\%new_item;
			### Get risk disease info
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
			push @badActDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
		}

		## get gooditems,goodEffect,goodDisease
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{changdaojunqun}{$prefix}{goodN};
			push @{$vars->{changdaojunqun}{$prefix}{gooditems}{items}{item}},\%new_item;
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
		}

	}

	## good,bad
	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	push @{$Floragood},@good;
	push @{$Florabad},@bad;
	push @{$FloragoodDesc},@goodDesc;
	push @{$FlorabadDesc},@badDesc;
	push @{$FloragoodEffect},@goodEffect;
	push @{$FloragoodDisease},@goodDisease;
	push @{$FlorabadEffect},@badEffect;
	push @{$FlorabadDisease},@badDisease;
	push @{$FlorabadActRule},@badActRule;
	push @{$FlorabadActDesc},@badActDesc;

	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActDesc,badDisease,goodDisease
	## get flora overview risk related effect/disease and actrule
	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDesc},$vars->{changdaojunqun}{$prefix}{badEffect},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{badActRule},$vars->{changdaojunqun}{$prefix}{badActDesc})=&process_Flora_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActDesc,\%goodDiseaseFenbu,\%badDiseaseFenbu);

	return;
}

#&item_Intestinal_Overview('1','gaikuang',$Data,$list,$vars,$dir,$language,$species,$summary,\@gooditem,\@baditem,\@goodDisease,\@badDisease,\@goodEffect,\@badEffect,\@badActRule,\@badActDesc);
sub item_Intestinal_Overview {
	my ($num,$prefix,$Data,$list,$vars,$dir,$language,$species,$summary,$Floragood,$Florabad,$FloragoodDesc,$FlorabadDesc,$FloragoodDisease,$FlorabadDisease,$FloragoodEffect,$FlorabadEffect,$FlorabadActRule,$FlorabadActDesc)=@_;
	my %goodDiseaseOverview=();
	my %badDiseaseOverview=();

	##
	my (@good,@bad);
	my (@badDesc,@badEffect,@badDisease,@badActRule,@badActDesc);
	my (@goodDesc,@goodEffect,@goodDisease);

	$vars->{changdaojunqun}{$prefix}{diversity}=0;
	$vars->{changdaojunqun}{$prefix}{divlevel}="-";
	$vars->{changdaojunqun}{$prefix}{divvalue}="-";
	$vars->{changdaojunqun}{$prefix}{benifitper}=0;
	$vars->{changdaojunqun}{$prefix}{benifilevel}="-";
	$vars->{changdaojunqun}{$prefix}{benifivalue}="-";
	$vars->{changdaojunqun}{$prefix}{harmper}=0;
	$vars->{changdaojunqun}{$prefix}{harmlevel}="-";
	$vars->{changdaojunqun}{$prefix}{harmvalue}="-";

	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;

	##
	foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		push @{$vars->{changdaojunqun}{$prefix}{items}{item}},$item;
		if ($meta eq "zxjunqunduoyangxing") {
			$vars->{changdaojunqun}{$prefix}{divlevel}=$item->{level};
			$vars->{changdaojunqun}{$prefix}{diversity}=sprintf "%.2f",$item->{per}*100;
			$vars->{changdaojunqun}{$prefix}{divvalue}=$item->{value};
		}
		if ($meta eq "zxyouyijun") {
			$vars->{changdaojunqun}{$prefix}{benifitper}=sprintf "%.2f",$item->{per}*100;
			$vars->{changdaojunqun}{$prefix}{benifilevel}=$item->{level};
			$vars->{changdaojunqun}{$prefix}{benifivalue}=$item->{value};
		}
		if ($meta eq "zxyouhaijun") {
			$vars->{changdaojunqun}{$prefix}{harmper}=sprintf "%.2f",$item->{per}*100;
			$vars->{changdaojunqun}{$prefix}{harmlevel}=$item->{level};
			$vars->{changdaojunqun}{$prefix}{harmvalue}=$item->{value};
		}

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{changdaojunqun}{$prefix}{abnormalN}+=$add;
		$vars->{changdaojunqun}{$prefix}{goodN}+=$addG;
		$vars->{changdaojunqun}{$prefix}{badN}+=$addN;

		## get abnormal
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{changdaojunqun}{$prefix}{badN};
			push @{$vars->{changdaojunqun}{$prefix}{baditems}{items}{item}},\%new_item;
			### Get risk disease info
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
			push @badActDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
		}

		## get gooditems,goodEffect,goodDisease
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{changdaojunqun}{$prefix}{goodN};
			push @{$vars->{changdaojunqun}{$prefix}{gooditems}{items}{item}},\%new_item;
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
		}

	}

	## good,bad
	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	push @{$Floragood},@good;
	push @{$Florabad},@bad;
	push @{$FloragoodDesc},@goodDesc;
	push @{$FlorabadDesc},@badDesc;
	push @{$FloragoodEffect},@goodEffect;
	push @{$FloragoodDisease},@goodDisease;
	push @{$FlorabadEffect},@badEffect;
	push @{$FlorabadDisease},@badDisease;
	push @{$FlorabadActRule},@badActRule;
	push @{$FlorabadActDesc},@badActDesc;

	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActDesc,badDisease,goodDisease
	## get flora overview risk related effect/disease and actrule
	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDesc},$vars->{changdaojunqun}{$prefix}{badEffect},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{badActRule},$vars->{changdaojunqun}{$prefix}{badActDesc})=&process_Flora_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActDesc,\%goodDiseaseOverview,\%badDiseaseOverview);

	$vars->{changdaojunqun}{$prefix}{total}=$summary->{Total_OTU};
	$vars->{changdaojunqun}{$prefix}{classified}=$summary->{Classified_OTU};
	$vars->{changdaojunqun}{$prefix}{unclassified}=$summary->{Unclassified_OTU};
	$vars->{changdaojunqun}{$prefix}{Firmicutes_Bacteroidetes_ratio}=$summary->{Firmicutes_Bacteroidetes_ratio};
	### get diversity pic
	&get_diversity_pic($vars);

	return;
}

#($vars->{changdaojunqun}{desc}{$prefix}{goodEffect},$vars->{changdaojunqun}{desc}{$prefix}{goodDisease},$vars->{changdaojunqun}{desc}{$prefix}{badEffect},$vars->{changdaojunqun}{desc}{$prefix}{badDisease},$vars->{changdaojunqun}{desc}{$prefix}{badActRule},$vars->{changdaojunqun}{desc}{$prefix}{badActDesc})=&process_Flora_DiseaseEffect(\@goodEffect,\@goodDisease,\@badEffect,\@badDisease,\@badActRule,\@badActDesc);
sub process_Flora_DiseaseEffect {
	my ($goodDesc,$goodEffect,$goodDisease,$badDesc,$badEffect,$badDisease,$badActRule,$badActDesc,$goodDiseaseFlora,$badDiseaseFlora)=@_;
	my ($goodDe,$goodE,$goodD,$badDe,$badE,$badD,$badActR,$badActD)=('-','-','-','-','-','-','-','-','-');

	##
	&split_uniq_hash($goodDisease,$goodDiseaseFlora);
	&split_uniq_hash($badDisease,$badDiseaseFlora);

	my ($gde)=&split_uniq_array($goodDesc);
	my ($ge)=&split_uniq_array($goodEffect);
	my ($gd)=&split_uniq_array($goodDisease);
	my ($bde)=&split_uniq_array($badDesc);
	my ($be)=&split_uniq_array($badEffect);
	my ($bd)=&split_uniq_array($badDisease);

	my ($actrule)=&split_uniq_array($badActRule);
	my ($actdesc)=&split_uniq_array($badActDesc);

	##
	if (@{$gde}) {
		if (scalar @{$gde} <= 3) {
			$goodDe=join "、",@{$gde};
		}
		else {
			$goodDe=join("、",$gde->[0],$gde->[1],$gde->[2]);
			$goodDe.="等";
		}
	}
	if (@{$ge}) {
		if (scalar @{$ge} <= 3) {
			$goodE=join "、",@{$ge};
		}
		else {
			$goodE=join("、",$ge->[0],$ge->[1],$ge->[2]);
			$goodE.="等";
		}
	}
	if (@{$gd}) {
		if (scalar @{$gd} <= 3) {
			$goodD=join "、",@{$gd};
		}
		else {
			$goodD=join("、",$gd->[0],$gd->[1],$gd->[2]);
			$goodD.="等";
		}
	}

	if (@{$bde}) {
		if (scalar @{$bde} <= 3) {
			$badDe=join "、",@{$bde};
		}
		else {
			$badDe=join("、",$bde->[0],$bde->[1],$bde->[2]);
			$badDe.="等";
		}
	}
	if (@{$be}) {
		if (scalar @{$be} <= 3) {
			$badE=join "、",@{$be};
		}
		else {
			$badE=join("、",$be->[0],$be->[1],$be->[2]);
			$badE.="等";
		}
	}
	if (@{$bd}) {
		if (scalar @{$bd} <= 3) {
			$badD=join "、",@{$bd};
		}
		else {
			$badD=join("、",$bd->[0],$bd->[1],$bd->[2]);
			$badD.="等";
		}
	}

	$badActR=join "、",@{$actrule};
	$badActD=join "、",@{$actdesc};

	$goodDe='-' if ($goodDe eq "");
	$goodE='-' if ($goodE eq "");
	$goodD='-' if ($goodD eq "");
	$badE='-' if ($badE eq "");

	$badDe='-' if ($badDe eq "");
	$badD='-' if ($badD eq "");
	$badActR='-' if ($badActR eq "");
	$badActD='-' if ($badActD eq "");

	return ($goodDe,$goodE,$goodD,$badDe,$badE,$badD,$badActR,$badActD);
}

#&get_diversity_pic($vars);
sub get_diversity_pic {
	my ($vars)=@_;
	
	##
	$vars->{changdaojunqun}{gaikuang}{divpic}="-";
	$vars->{changdaojunqun}{pic}="-";
	if ($vars->{changdaojunqun}{gaikuang}{divlevel} eq 'l1' || $vars->{changdaojunqun}{gaikuang}{divlevel} eq 'l2') {
		$vars->{changdaojunqun}{gaikuang}{divpic}="DIVERSITY-low.pdf";
		$vars->{changdaojunqun}{pic}="MICROBIOTA-low.pdf";
	}
	elsif ($vars->{changdaojunqun}{gaikuang}{divlevel} eq 'l3') {
		$vars->{changdaojunqun}{gaikuang}{divpic}="DIVERSITY-middle.pdf";
		$vars->{changdaojunqun}{pic}="MICROBIOTA-middle.pdf";
	}
	else {
		$vars->{changdaojunqun}{gaikuang}{divpic}="DIVERSITY-high.pdf";
		$vars->{changdaojunqun}{pic}="MICROBIOTA-high.pdf";
	}
	
	return;
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

#&Inflammation_and_Immune($Data,$listH->{meta}{'炎症和免疫'},$vars,"$indir/XML/meta/$language",$language,$regulate);
sub Inflammation_and_Immune {
	my ($Data,$list,$vars,$dir,$language,$regulate)=@_;

	##
	$vars->{yanzhenghemianyi}{name}="炎症和免疫";
	$vars->{yanzhenghemianyi}{ename}="INFLAMMATION";
	$vars->{yanzhenghemianyi}{abnormal}{num}=0;
	
	## get item
	&item_Inflammation_and_Immune('1','giimmune',$Data,$list,$vars,$dir,$language,$regulate);
	
	## get green|yellow|red.pdf
	&get_Inflammation_and_Immune_pic($vars);
	
	return;
}

#&item_Inflammation_and_Immune('1','items',$Data,$list,$vars,$dir,$language,$regulate);
sub item_Inflammation_and_Immune {
	my ($num,$prefix,$Data,$list,$vars,$dir,$language,$regulate)=@_;
	my %badDiseaseTotal=();
	my %goodDiseaseTotal=();

	##
	my (@good,@bad,@goodDesc,@badDesc);
	my (@badEffect,@badDisease,@badFoodRule,@badFoodMore,@badFoodLess,@badActRule,@badActMore,@badActLess);
	my (@goodEffect,@goodDisease);
	##
	$vars->{yanzhenghemianyi}{$prefix}{explain}{totalN}=0;
	$vars->{yanzhenghemianyi}{$prefix}{explain}{abnormalN}=0;
	$vars->{yanzhenghemianyi}{$prefix}{explain}{goodN}=0;
	$vars->{yanzhenghemianyi}{$prefix}{explain}{badN}=0;
	##
	foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}) {
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});
		$vars->{yanzhenghemianyi}{$prefix}{explain}{totalN}++;
		
		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		push @{$vars->{yanzhenghemianyi}{$prefix}{items}{item}},$item;
		
		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{yanzhenghemianyi}{$prefix}{explain}{abnormalN}+=$add;
		$vars->{yanzhenghemianyi}{$prefix}{explain}{goodN}+=$addG;
		$vars->{yanzhenghemianyi}{$prefix}{explain}{badN}+=$addN;
		
		## get abnormal
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			##
			my %new_item = %{$item};
			$new_item{index} = $vars->{yanzhenghemianyi}{$prefix}{explain}{badN};
			push @{$vars->{yanzhenghemianyi}{$prefix}{baditems}{items}{item}},\%new_item;
			###
			### Get risk disease info
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badFoodRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} ne "");
			push @badFoodMore,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} ne "");
			push @badFoodLess,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} ne "");
			push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
			push @badActMore,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} ne "");
			push @badActLess,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} ne "");
		}
		## get abnormal
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			##
			my %new_item = %{$item};
			$new_item{index} = $vars->{yanzhenghemianyi}{$prefix}{explain}{goodN};
			push @{$vars->{yanzhenghemianyi}{$prefix}{gooditems}{items}{item}},\%new_item;
			###  Get benifit disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
		}

		## get regulate
		#if ($meta eq "zxmianyizhishu") {
		#	if ($item->{level} eq 'l1' || $item->{level} eq 'l2') {
		#		push @{$regulate},{name => $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddesc}, desc => $vd->{$language}{suggestion}{healthy}{$item->{level}}{aidnutrition}, spotnum => 0};
		#	}
		#}
	}
	$vars->{yanzhenghemianyi}{abnormal}{num}=$vars->{yanzhenghemianyi}{$prefix}{explain}{badN};

	## good,bad
	($vars->{yanzhenghemianyi}{$prefix}{explain}{good},$vars->{yanzhenghemianyi}{$prefix}{explain}{bad})=&process_good_bad(\@good,\@bad);
	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalbadFoodRule,TotalbadFoodMore,TotalbadFoodLess,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActMore,TotalbadActLess,badDiseaseTotal,goodDiseaseTotal
	## get protein fat carbon risk related effect/disease and actrule/foodrule
	($vars->{yanzhenghemianyi}{$prefix}{explain}{goodDesc},$vars->{yanzhenghemianyi}{$prefix}{explain}{goodEffect},$vars->{yanzhenghemianyi}{$prefix}{explain}{goodDisease},$vars->{yanzhenghemianyi}{$prefix}{explain}{badDesc},$vars->{yanzhenghemianyi}{$prefix}{explain}{badEffect},$vars->{yanzhenghemianyi}{$prefix}{explain}{badDisease},$vars->{yanzhenghemianyi}{$prefix}{explain}{badActRule},$vars->{yanzhenghemianyi}{$prefix}{explain}{badActMore},$vars->{yanzhenghemianyi}{$prefix}{explain}{badActLess},$vars->{yanzhenghemianyi}{$prefix}{explain}{badFoodRule},$vars->{yanzhenghemianyi}{$prefix}{explain}{badFoodMore},$vars->{yanzhenghemianyi}{$prefix}{explain}{badFoodLess})=&process_goodbad_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActMore,\@badActLess,\@badFoodRule,\@badFoodMore,\@badFoodLess,\%goodDiseaseTotal,\%badDiseaseTotal);

	return;
}

#&get_Inflammation_and_Immune_pic($vars);
sub get_Inflammation_and_Immune_pic {
	my ($vars)=@_;
	
	##
	$vars->{yanzhenghemianyi}{curpic}="-";
	$vars->{yanzhenghemianyi}{pic}="-";
	if ($vars->{yanzhenghemianyi}{giimmune}{explain}{badN} == 0) {
		$vars->{yanzhenghemianyi}{curpic}="sun-green.pdf";
		$vars->{yanzhenghemianyi}{pic}="INFLAMMATION-green.pdf";
	}
	elsif ($vars->{yanzhenghemianyi}{giimmune}{explain}{badN} <= 2) {
		$vars->{yanzhenghemianyi}{curpic}="sun-yellow.pdf";
		$vars->{yanzhenghemianyi}{pic}="INFLAMMATION-yellow.pdf";
	}
	else {
		$vars->{yanzhenghemianyi}{curpic}="sun-red.pdf";
		$vars->{yanzhenghemianyi}{pic}="INFLAMMATION-red.pdf";
	}
	
	return;
}


#&Stool_Status($listH->{stoolList}{"粪便状态"},$vars,"$indir/XML/stool/$language",$language,$samInfo);
sub Stool_Status {
	my ($meta,$vars,$dir,$language,$samInfo)=@_;
	
	##
	my ($texture,$color)=($samInfo->{"粪便质地"},$samInfo->{"粪便颜色"});
	my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fenbianzhuangtai}{name}="粪便状态";
	$vars->{fenbianzhuangtai}{picdesc}='-';
	$vars->{fenbianzhuangtai}{pic}='-';
	$vars->{fenbianzhuangtai}{colorpic}='-';
	$vars->{fenbianzhuangtai}{desc}{itemnum}=0;
	##
	if (defined $texture && $texture ne "") {
		$vars->{fenbianzhuangtai}{picdesc}=$texture;
		$vars->{fenbianzhuangtai}{pic}=&item_Stool_Status_texture($vd->{$language}{metatexture},$texture,$vars);
	}
	if (defined $color && $color ne "") {
		$vars->{fenbianzhuangtai}{colorpic}=&item_Stool_Status_color($vd->{$language}{metacolor},$color,$vars);
	}
	
	return;
}

#my ($regulateLow,$badDiseaseTotal)=&Absorption_and_Metabolism($Data,$listH->{absorption}{'消化和吸收'},$vars,"$indir/XML/meta/$language",$language);
sub Absorption_and_Metabolism {
	my ($Data,$list,$vars,$dir,$language,$vitamin)=@_;
	my @regulateLow=();
	my %badDiseaseTotal=();
	my %goodDiseaseTotal=();
	my %AbsorbbadDisease=();
	my %AbsorbgoodDisease=();
	my %TotalbadDisease=();
	my %TotalgoodDisease=();

	##
	my (@gooditem,@baditem,@goodDesc,@badDesc);							    ## protein|fat|carbon good and bad effect items
	my (@goodDisease,@badDisease,@badFoodRule,@badFoodMore,@badFoodLess);   ## protein|fat|carbon good and bad items disease related
	my (@goodEffect,@badEffect,@badActRule,@badActMore,@badActLess);        ## protein|fat|carbon good and bad items disease related

	$vars->{xiaohuahexishou}{name}="消化和吸收";
	$vars->{xiaohuahexishou}{ename}="INSUFFICIENCY";
	$vars->{xiaohuahexishou}{explain}{TotalabnormalN}=0;				    ## protein|fat|carbon
	$vars->{xiaohuahexishou}{explain}{TotalbadN}=0;						    ## protein|fat|carbon
	$vars->{xiaohuahexishou}{explain}{TotalgoodN}=0;
	$vars->{xiaohuahexishou}{abnormal}{num}=0;

	$vars->{xiaohuahexishou}{explain}{AbsorbabnormalN}=0;
	$vars->{xiaohuahexishou}{explain}{AbsorbbadN}=0;
	$vars->{xiaohuahexishou}{explain}{AbsorbgoodN}=0;
	$vars->{xiaohuahexishou}{abnormal}{num}=0;
	$vars->{xiaohuahexishou}{baditems}{num}=0;
	$vars->{xiaohuahexishou}{gooditems}{num}=0;

	## 总评
	&item_Digestion_and_Absorption('1','zongping',$Data,$list->{absorption}{'消化和吸收'},$vars,"$dir/absorption/$language",$language,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%goodDiseaseTotal,\%badDiseaseTotal);

	## 蛋白质的消化和吸收
	&item_Digestion_and_Absorption('2','protein',$Data,$list->{absorption}{'消化和吸收'},$vars,"$dir/absorption/$language",$language,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%goodDiseaseTotal,\%badDiseaseTotal);

	## 脂肪的消化和吸收
	&item_Digestion_and_Absorption('3','fat',$Data,$list->{absorption}{'消化和吸收'},$vars,"$dir/absorption/$language",$language,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%goodDiseaseTotal,\%badDiseaseTotal);

	## 碳水化合物的消化和吸收
	&item_Digestion_and_Absorption('4','carbon',$Data,$list->{absorption}{'消化和吸收'},$vars,"$dir/absorption/$language",$language,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%goodDiseaseTotal,\%badDiseaseTotal);

	## protein、fat、carbon共15小项
	($vars->{xiaohuahexishou}{explain}{Absorbgood},$vars->{xiaohuahexishou}{explain}{Absorbbad})=&process_good_bad(\@gooditem,\@baditem);
	($vars->{xiaohuahexishou}{explain}{AbsorbgoodDesc},$vars->{xiaohuahexishou}{explain}{AbsorbgoodEffect},$vars->{xiaohuahexishou}{explain}{AbsorbgoodDisease},$vars->{xiaohuahexishou}{explain}{AbsorbbadDesc},$vars->{xiaohuahexishou}{explain}{AbsorbbadEffect},$vars->{xiaohuahexishou}{explain}{AbsorbbadDisease},$vars->{xiaohuahexishou}{explain}{AbsorbbadActRule},$vars->{xiaohuahexishou}{explain}{AbsorbbadActMore},$vars->{xiaohuahexishou}{explain}{AbsorbbadActLess},$vars->{xiaohuahexishou}{explain}{AbsorbbadFoodRule},$vars->{xiaohuahexishou}{explain}{AbsorbbadFoodMore},$vars->{xiaohuahexishou}{explain}{AbsorbbadFoodLess})=&process_goodbad_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActMore,\@badActLess,\@badFoodRule,\@badFoodMore,\@badFoodLess,\%AbsorbgoodDisease,\%AbsorbbadDisease);

	## 菌群代谢
	&item_Micro_Metabolism('1','metabolism',$Data,$list->{metabolism}{'菌群代谢'},$vars,"$dir/metabolism/$language",$language,$vitamin,\@regulateLow,\@gooditem,\@baditem,\@goodDesc,\@badDesc,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%goodDiseaseTotal,\%badDiseaseTotal);

	## protein、fat、carbon + metabolism
	($vars->{xiaohuahexishou}{explain}{Totalgood},$vars->{xiaohuahexishou}{explain}{Totalbad})=&process_good_bad(\@gooditem,\@baditem);
	($vars->{xiaohuahexishou}{explain}{TotalgoodDesc},$vars->{xiaohuahexishou}{explain}{TotalgoodEffect},$vars->{xiaohuahexishou}{explain}{TotalgoodDisease},$vars->{xiaohuahexishou}{explain}{TotalbadDesc},$vars->{xiaohuahexishou}{explain}{TotalbadEffect},$vars->{xiaohuahexishou}{explain}{TotalbadDisease},$vars->{xiaohuahexishou}{explain}{TotalbadActRule},$vars->{xiaohuahexishou}{explain}{TotalbadActMore},$vars->{xiaohuahexishou}{explain}{TotalbadActLess},$vars->{xiaohuahexishou}{explain}{TotalbadFoodRule},$vars->{xiaohuahexishou}{explain}{TotalbadFoodMore},$vars->{xiaohuahexishou}{explain}{TotalbadFoodLess})=&process_goodbad_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActMore,\@badActLess,\@badFoodRule,\@badFoodMore,\@badFoodLess,\%TotalgoodDisease,\%TotalbadDisease);

	## get green|yellow|red.pdf
	&get_Absorption_pic($vars);
	
	return (\@regulateLow,\%goodDiseaseTotal,\%badDiseaseTotal);
}

#&item_Digestion_and_Absorption('1','zongping',$Data,$list,$vars,$dir,$language,\@gooditem,\@baditem,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%badDiseaseTotal);
sub item_Digestion_and_Absorption {
	my ($num,$prefix,$Data,$list,$vars,$dir,$language,$Totalgood,$Totalbad,$TotalgoodDesc,$TotalbadDesc,$TotalgoodDisease,$TotalbadDisease,$TotalbadFoodRule,$TotalbadFoodMore,$TotalbadFoodLess,$TotalgoodEffect,$TotalbadEffect,$TotalbadActRule,$TotalbadActMore,$TotalbadActLess,$badDiseaseTotal,$goodDiseaseTotal)=@_;
	
	##
	my (@good,@bad,@goodDesc,@badDesc);
	my (@badEffect,@badDisease,@badFoodRule,@badFoodMore,@badFoodLess,@badActRule,@badActMore,@badActLess);
	my (@goodEffect,@goodDisease);
	##
	$vars->{xiaohuahexishou}{$prefix}{explain}{abnormalN}=0;
	$vars->{xiaohuahexishou}{$prefix}{explain}{goodN}=0;
	$vars->{xiaohuahexishou}{$prefix}{explain}{badN}=0;
	##
	foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}) {
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");

		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		$item->{key}=&get_item_key($item) if ($list->{$num}{name} eq "消化和吸收总评");
		push @{$vars->{xiaohuahexishou}{$prefix}{items}{item}},$item;

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{xiaohuahexishou}{$prefix}{explain}{abnormalN}+=$add;
		$vars->{xiaohuahexishou}{$prefix}{explain}{goodN}+=$addG;
		$vars->{xiaohuahexishou}{$prefix}{explain}{badN}+=$addN;
		if ($prefix eq "protein" || $prefix eq "fat" || $prefix eq "carbon") {
			$vars->{xiaohuahexishou}{explain}{TotalabnormalN}+=$add;
			$vars->{xiaohuahexishou}{explain}{TotalbadN}+=$addN;
			$vars->{xiaohuahexishou}{explain}{TotalgoodN}+=$addG;
			$vars->{xiaohuahexishou}{explain}{AbsorbabnormalN}+=$add;
			$vars->{xiaohuahexishou}{explain}{AbsorbbadN}+=$addN;
			$vars->{xiaohuahexishou}{explain}{AbsorbgoodN}+=$addG;
			$vars->{xiaohuahexishou}{abnormal}{num}+=$addN;
			$vars->{xiaohuahexishou}{baditems}{num}+=$addN;
			$vars->{xiaohuahexishou}{gooditems}{num}+=$addG;
		}

		## get baditems,badEffect,badDisease,badFoodRule,badFoodMore,badFoodLess
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			###
			#push @{$vars->{xiaohuahexishou}{$prefix}{baditems}{items}{item}},{index => $vars->{xiaohuahexishou}{desc}{$prefix}{badN}, name => $item->{name}, risk => $item->{level}};
			my %new_item = %{$item};
			$new_item{index} = $vars->{xiaohuahexishou}{$prefix}{explain}{badN};
			push @{$vars->{xiaohuahexishou}{$prefix}{baditems}{items}{item}},\%new_item;
			$new_item{index} = $vars->{xiaohuahexishou}{explain}{AbsorbbadN};
			push @{$vars->{xiaohuahexishou}{baditems}{items}{item}},\%new_item if ($prefix eq "protein" || $prefix eq "fat" || $prefix eq "carbon");
			### Get risk disease info
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badFoodRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} ne "");
			push @badFoodMore,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} ne "");
			push @badFoodLess,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} ne "");
			push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
			push @badActMore,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} ne "");
			push @badActLess,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} ne "");
		}
		## get gooditems,goodEffect,goodDisease
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{xiaohuahexishou}{$prefix}{explain}{goodN};
			push @{$vars->{xiaohuahexishou}{$prefix}{gooditems}{items}{item}},\%new_item;
			$new_item{index} = $vars->{xiaohuahexishou}{explain}{AbsorbbadN};
			push @{$vars->{xiaohuahexishou}{gooditems}{items}{item}},\%new_item if ($prefix eq "protein" || $prefix eq "fat" || $prefix eq "carbon");
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
		}
	}
	$vars->{xiaohuahexishou}{$prefix}{abnormal}{num}=$vars->{xiaohuahexishou}{$prefix}{explain}{badN};
	
	## good,bad
	($vars->{xiaohuahexishou}{$prefix}{explain}{good},$vars->{xiaohuahexishou}{$prefix}{explain}{bad})=&process_good_bad(\@good,\@bad);
	if ($prefix eq "protein" || $prefix eq "fat" || $prefix eq "carbon") {
		push @{$Totalgood},@good;
		push @{$Totalbad},@bad;
		push @{$TotalgoodDesc},@goodDesc;
		push @{$TotalbadDesc},@badDesc;
		push @{$TotalgoodEffect},@goodEffect;
		push @{$TotalgoodDisease},@goodDisease;
		push @{$TotalbadEffect},@badEffect;
		push @{$TotalbadDisease},@badDisease;
		push @{$TotalbadActRule},@badActRule;
		push @{$TotalbadActMore},@badActMore;
		push @{$TotalbadActLess},@badActLess;
		push @{$TotalbadFoodRule},@badFoodRule;
		push @{$TotalbadFoodMore},@badFoodMore;
		push @{$TotalbadFoodLess},@badFoodLess;
	}
	
	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalbadFoodRule,TotalbadFoodMore,TotalbadFoodLess,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActMore,TotalbadActLess,badDiseaseTotal,goodDiseaseTotal
	## get protein fat carbon risk related effect/disease and actrule/foodrule
	($vars->{xiaohuahexishou}{$prefix}{explain}{goodDesc},$vars->{xiaohuahexishou}{$prefix}{explain}{goodEffect},$vars->{xiaohuahexishou}{$prefix}{explain}{goodDisease},$vars->{xiaohuahexishou}{$prefix}{explain}{badDesc},$vars->{xiaohuahexishou}{$prefix}{explain}{badEffect},$vars->{xiaohuahexishou}{$prefix}{explain}{badDisease},$vars->{xiaohuahexishou}{$prefix}{explain}{badActRule},$vars->{xiaohuahexishou}{$prefix}{explain}{badActMore},$vars->{xiaohuahexishou}{$prefix}{explain}{badActLess},$vars->{xiaohuahexishou}{$prefix}{explain}{badFoodRule},$vars->{xiaohuahexishou}{$prefix}{explain}{badFoodMore},$vars->{xiaohuahexishou}{$prefix}{explain}{badFoodLess})=&process_goodbad_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActMore,\@badActLess,\@badFoodRule,\@badFoodMore,\@badFoodLess,$goodDiseaseTotal,$badDiseaseTotal);

	return;
}

#&item_Micro_Metabolism('1','metabolism',$Data,$list,$vars,$dir,$language,$vitamin,\@regulate,\@gooditem,\@baditem,\@goodDisease,\@badDisease,\@badFoodRule,\@badFoodMore,\@badFoodLess,\@goodEffect,\@badEffect,\@badActRule,\@badActMore,\@badActLess,\%badDiseaseTotal);
sub item_Micro_Metabolism {
	my ($num,$prefix,$Data,$list,$vars,$dir,$language,$vitamin,$regulate,$Totalgood,$Totalbad,$TotalgoodDesc,$TotalbadDesc,$TotalgoodDisease,$TotalbadDisease,$TotalbadFoodRule,$TotalbadFoodMore,$TotalbadFoodLess,$TotalgoodEffect,$TotalbadEffect,$TotalbadActRule,$TotalbadActMore,$TotalbadActLess,$badDiseaseTotal,$goodDiseaseTotal)=@_;
	
	##
	#my $vitaminN=0;
	my @vitaminBad=();
	my (@good,@bad,@goodDesc,@badDesc);
	my (@badEffect,@badDisease,@badFoodRule,@badFoodMore,@badFoodLess,@badActRule,@badActMore,@badActLess);
	my (@goodEffect,@goodDisease);
	##
	$vars->{xiaohuahexishou}{$prefix}{explain}{abnormalN}=0;
	$vars->{xiaohuahexishou}{$prefix}{explain}{goodN}=0;
	$vars->{xiaohuahexishou}{$prefix}{explain}{badN}=0;
	$vars->{xiaohuahexishou}{$prefix}{explain}{vitaminbadN}=0;
	$vars->{xiaohuahexishou}{$prefix}{explain}{fumeiQbad} = "N";

	##
	foreach my $meta (sort {$list->{$num}{section}{$a}{order} <=> $list->{$num}{section}{$b}{order}} keys %{$list->{$num}{section}}) {
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");

		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
		$item->{key}=&get_item_key($item) if ($list->{$num}{name} eq "消化和吸收总评");
		push @{$vars->{xiaohuahexishou}{$prefix}{items}{item}},$item;

		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
		$vars->{xiaohuahexishou}{$prefix}{explain}{abnormalN}+=$add;
		$vars->{xiaohuahexishou}{$prefix}{explain}{goodN}+=$addG;
		$vars->{xiaohuahexishou}{$prefix}{explain}{badN}+=$addN;
		$vars->{xiaohuahexishou}{explain}{TotalabnormalN}+=$add;
		$vars->{xiaohuahexishou}{explain}{TotalbadN}+=$addN;
		$vars->{xiaohuahexishou}{explain}{TotalgoodN}+=$addG;

		## get baditems,badEffect,badDisease,badFoodRule,badFoodMore,badFoodLess
		if (defined $item->{risk} && $item->{risk} eq "有害") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{xiaohuahexishou}{$prefix}{explain}{badN};
			push @{$vars->{xiaohuahexishou}{$prefix}{baditems}{items}{item}},\%new_item;
			### Get risk disease info
			push @badDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @badEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @badDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @badFoodRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietdesc} ne "");
			push @badFoodMore,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietmore} ne "");
			push @badFoodLess,$vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{dietless} ne "");
			push @badActRule,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionrule} ne "");
			push @badActMore,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionmore} ne "");
			push @badActLess,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actionless} ne "");
		}
		## get gooditems,goodEffect,goodDisease
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			###
			my %new_item = %{$item};
			$new_item{index} = $vars->{xiaohuahexishou}{$prefix}{explain}{goodN};
			push @{$vars->{xiaohuahexishou}{$prefix}{gooditems}{items}{item}},\%new_item;
			### Get risk disease info
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
		}

		## get regulate
		if ($meta eq "zxfumeiQ") {
			if ($item->{level} eq 'l1' || $item->{level} eq 'l2') {
				$vars->{xiaohuahexishou}{$prefix}{explain}{fumeiQbad} = "Y";
				#push @{$regulate},{name => $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddesc}, desc => $vd->{$language}{suggestion}{healthy}{$item->{level}}{aidnutrition}, spotnum => 0};
			}
		}
		if (exists $vitamin->{$meta}) {
			if ($item->{level} eq 'l1' || $item->{level} eq 'l2') {
				$vars->{xiaohuahexishou}{$prefix}{explain}{vitaminbadN} ++;
				push @vitaminBad,$item->{name};
				#$vitaminN++;
				#push @{$regulate},{name => $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddesc}, desc => $vd->{$language}{suggestion}{healthy}{$item->{level}}{aidnutrition}, spotnum => 0};
			}
		}
	}
	&process_vitaminbad($vars->{xiaohuahexishou}{$prefix}{explain},$vars->{xiaohuahexishou}{$prefix}{explain}{vitaminbadN},\@vitaminBad) if ($prefix eq "metabolism");
	$vars->{xiaohuahexishou}{$prefix}{explain}{abnormal}{num}=$vars->{xiaohuahexishou}{$prefix}{explain}{badN};
	
	## good,bad
	($vars->{xiaohuahexishou}{$prefix}{explain}{good},$vars->{xiaohuahexishou}{$prefix}{explain}{bad})=&process_good_bad(\@good,\@bad);
	push @{$Totalgood},@good;
	push @{$Totalbad},@bad;
	push @{$TotalgoodDesc},@goodDesc;
	push @{$TotalbadDesc},@badDesc;
	push @{$TotalgoodEffect},@goodEffect;
	push @{$TotalgoodDisease},@goodDisease;
	push @{$TotalbadEffect},@badEffect;
	push @{$TotalbadDisease},@badDisease;
	push @{$TotalbadActRule},@badActRule;
	push @{$TotalbadActMore},@badActMore;
	push @{$TotalbadActLess},@badActLess;
	push @{$TotalbadFoodRule},@badFoodRule;
	push @{$TotalbadFoodMore},@badFoodMore;
	push @{$TotalbadFoodLess},@badFoodLess;

	## Totalgood,Totalbad,TotalgoodDisease,TotalbadDisease,TotalbadFoodRule,TotalbadFoodMore,TotalbadFoodLess,TotalgoodEffect,TotalbadEffect,TotalbadActRule,TotalbadActMore,TotalbadActLess,badDiseaseTotal,goodDiseaseTotal
	## get protein fat carbon risk related effect/disease and actrule/foodrule
	($vars->{xiaohuahexishou}{$prefix}{explain}{goodDesc},$vars->{xiaohuahexishou}{$prefix}{explain}{goodEffect},$vars->{xiaohuahexishou}{$prefix}{explain}{goodDisease},$vars->{xiaohuahexishou}{$prefix}{explain}{badDesc},$vars->{xiaohuahexishou}{$prefix}{explain}{badEffect},$vars->{xiaohuahexishou}{$prefix}{explain}{badDisease},$vars->{xiaohuahexishou}{$prefix}{explain}{badActRule},$vars->{xiaohuahexishou}{$prefix}{explain}{badActMore},$vars->{xiaohuahexishou}{$prefix}{explain}{badActLess},$vars->{xiaohuahexishou}{$prefix}{explain}{badFoodRule},$vars->{xiaohuahexishou}{$prefix}{explain}{badFoodMore},$vars->{xiaohuahexishou}{$prefix}{explain}{badFoodLess})=&process_goodbad_DiseaseEffect(\@goodDesc,\@goodEffect,\@goodDisease,\@badDesc,\@badEffect,\@badDisease,\@badActRule,\@badActMore,\@badActLess,\@badFoodRule,\@badFoodMore,\@badFoodLess,$goodDiseaseTotal,$badDiseaseTotal);

	$vars->{xiaohuahexishou}{$prefix}{name}="代谢平衡";
	$vars->{xiaohuahexishou}{$prefix}{ename}="IMBALANCE";
	&get_IMBALANCE_pic($vars);

	return;
}

#&process_vitaminbad($vars,$vitaminbadN,$vitaminBad);
sub process_vitaminbad {
	my ($vars,$vitaminbadN,$vitaminBad)=@_;

	my $str = "维生素";
	if ($vitaminbadN == 0) {
		$vars->{vitaminBad} = "-";
	}
	elsif ($vitaminbadN <= 2) {
		foreach my $vitaminid (@{$vitaminBad}) {
			$vitaminid =~/([a-zA-Z\d]+)/;
			$str .= "$1"."、";
		}
		$str =~s/、$//;
		$vars->{vitaminBad} = $str;
	}
	else {
		for (0..1) {
			my $vitaminid = ${$vitaminBad}[$_];
			$vitaminid =~/([a-zA-Z\d]+)/;
			$str .= "$1"."、";
		}
		$str =~s/、$//;
		$vars->{vitaminBad} = "$str"."等";
	}
}

#&get_Absorption_pic($vars);
sub get_Absorption_pic {
	my ($vars)=@_;
	
	##
	$vars->{xiaohuahexishou}{curpic}="-";
	$vars->{xiaohuahexishou}{pic}="-";
	if ($vars->{xiaohuahexishou}{explain}{AbsorbbadN} == 0) {
		$vars->{xiaohuahexishou}{curpic}="sun-green.pdf";
		$vars->{xiaohuahexishou}{pic}="INSUFFICIENCY-green.pdf";
	}
	elsif ($vars->{xiaohuahexishou}{explain}{AbsorbbadN} <= 3) {
		$vars->{xiaohuahexishou}{curpic}="sun-yellow.pdf";
		$vars->{xiaohuahexishou}{pic}="INSUFFICIENCY-yellow.pdf";
	}
	else {
		$vars->{xiaohuahexishou}{curpic}="sun-red.pdf";
		$vars->{xiaohuahexishou}{pic}="INSUFFICIENCY-red.pdf";
	}
	
	return;
}

#&get_IMBALANCE_pic($vars);
sub get_IMBALANCE_pic {
	my ($vars)=@_;
	
	##
	$vars->{xiaohuahexishou}{metabolism}{curpic}="-";
	$vars->{xiaohuahexishou}{metabolism}{pic}="-";
	if ($vars->{xiaohuahexishou}{metabolism}{explain}{badN} == 0) {
		$vars->{xiaohuahexishou}{metabolism}{curpic}="sun-green.pdf";
		$vars->{xiaohuahexishou}{metabolism}{pic}="IMBALANCE-green.pdf";
	}
	elsif ($vars->{xiaohuahexishou}{metabolism}{explain}{badN} <= 3) {
		$vars->{xiaohuahexishou}{metabolism}{curpic}="sun-yellow.pdf";
		$vars->{xiaohuahexishou}{metabolism}{pic}="IMBALANCE-yellow.pdf";
	}
	else {
		$vars->{xiaohuahexishou}{metabolism}{curpic}="sun-red.pdf";
		$vars->{xiaohuahexishou}{metabolism}{pic}="IMBALANCE-red.pdf";
	}
	
	return;
}

#$vars->{fenbianzhuangtai}{pic}=&item_Stool_Status_texture($vd->{$language}{metatexture},$texture,$vars);
sub item_Stool_Status_texture {
	my ($hash,$texture,$vars)=@_;
	
	##
	$vars->{fenbianzhuangtai}{desc}{itemnum}++;
	my %item;
	## control
	$item{name}="质地";
	$item{control}{alias}="第四型";
	$item{control}{desc}="理想的便型";
	$item{control}{pic}="texture_normal.pdf";
	## sample
	$item{sample}{alias}="第四型";
	$item{sample}{desc}="理想的便型";
	$item{sample}{pic}="texture_part4.pdf";
	foreach my $key (keys %{$hash->{desc}}) {
		next if ($key eq "name");
		if ($hash->{desc}{$key}{texture} eq $texture) {
			$item{sample}{alias}=$texture;
			$item{sample}{desc}=$hash->{desc}{$key}{risk};
			$item{sample}{pic}="texture_".$key.".pdf";
		}
	}
	##
	push @{$vars->{fenbianzhuangtai}{desc}{items}{item}},\%item;
	
	return ($item{sample}{pic});
}

#$vars->{fenbianzhuangtai}{colorpic}=&item_Stool_Status_color($vd->{$language}{metacolor},$color,$vars);
sub item_Stool_Status_color {
	my ($hash,$color,$vars)=@_;
	
	##
	$vars->{fenbianzhuangtai}{desc}{itemnum}++;
	my %item;
	## control
	$item{name}="颜色";
	$item{control}{alias}="棕色";
	$item{control}{desc}="正常";
	$item{control}{pic}="color_part1.pdf";
	## sample
	$item{sample}{alias}="棕色";
	$item{sample}{desc}="正常";
	$item{sample}{pic}="color_part1.pdf";
	foreach my $key (keys %{$hash->{desc}}) {
		next if ($key eq "name");
		if ($hash->{desc}{$key}{color} eq $color) {
			$item{sample}{alias}=$color;
			$item{sample}{desc}=$hash->{desc}{$key}{risk};
			$item{sample}{pic}="color_".$key.".pdf";
		}
	}
	##
	push @{$vars->{fenbianzhuangtai}{desc}{items}{item}},\%item;
	
	return ($item{sample}{pic});
}

#($vars->{xiaohuahexishou}{desc}{$prefix}{goodEffect},$vars->{xiaohuahexishou}{desc}{$prefix}{goodDisease},$vars->{xiaohuahexishou}{desc}{$prefix}{badEffect},$vars->{xiaohuahexishou}{desc}{$prefix}{badDisease},$vars->{xiaohuahexishou}{desc}{$prefix}{badActRule},$vars->{xiaohuahexishou}{desc}{$prefix}{badActMore},$vars->{xiaohuahexishou}{desc}{$prefix}{badActLess},$vars->{xiaohuahexishou}{desc}{$prefix}{badFoodRule},$vars->{xiaohuahexishou}{desc}{$prefix}{badFoodMore},$vars->{xiaohuahexishou}{desc}{$prefix}{badFoodLess})=&process_goodbad_DiseaseEffect(\@goodEffect,\@goodDisease,\@badEffect,\@badDisease,\@badActRule,\@badActMore,\@badActLess,\@badFoodRule,\@badFoodMore,\@badFoodLess,$goodDiseaseTotal,$badDiseaseTotale);
sub process_goodbad_DiseaseEffect {
	my ($goodDesc,$goodEffect,$goodDisease,$badDesc,$badEffect,$badDisease,$badActRule,$badActMore,$badActLess,$badFoodRule,$badFoodMore,$badFoodLess,$goodDiseaseTotal,$badDiseaseTotal)=@_;
	my ($goodDe,$goodE,$goodD,$badDe,$badE,$badD,$badActR,$badActM,$badActL,$badFoodR,$badFoodM,$badFoodL)=('-','-','-','-','-','-','-','-','-','-','-','-');

	##
	&split_uniq_hash($goodDisease,$goodDiseaseTotal);
	&split_uniq_hash($badDisease,$badDiseaseTotal);

	my ($gde)=&split_uniq_array($goodDesc);
	my ($ge)=&split_uniq_array($goodEffect);
	my ($gd)=&split_uniq_array($goodDisease);
	my ($bde)=&split_uniq_array($badDesc);
	my ($be)=&split_uniq_array($badEffect);
	my ($bd)=&split_uniq_array($badDisease);

	my ($actrule)=&split_uniq_array($badActRule);
	my ($actmore)=&split_uniq_array($badActMore);
	my ($actless)=&split_uniq_array($badActLess);

	my ($foodrule)=&split_uniq_array($badFoodRule);
	my ($foodmore)=&split_uniq_array($badFoodMore);
	my ($foodless)=&split_uniq_array($badFoodLess);

	##
	if (@{$gde}) {
		if (scalar @{$gde} <= 3) {
			$goodDe=join "、",@{$gde};
		}
		else {
			$goodDe=join("、",$gde->[0],$gde->[1],$gde->[2]);
			$goodDe.="等";
		}
	}
	if (@{$ge}) {
		if (scalar @{$ge} <= 3) {
			$goodE=join "、",@{$ge};
		}
		else {
			$goodE=join("、",$ge->[0],$ge->[1],$ge->[2]);
			$goodE.="等";
		}
	}
	if (@{$gd}) {
		if (scalar @{$gd} <= 3) {
			$goodD=join "、",@{$gd};
		}
		else {
			$goodD=join("、",$gd->[0],$gd->[1],$gd->[2]);
			$goodD.="等";
		}
	}

		if (@{$bde}) {
		if (scalar @{$bde} <= 3) {
			$badDe=join "、",@{$bde};
		}
		else {
			$badDe=join("、",$bde->[0],$bde->[1],$bde->[2]);
			$badDe.="等";
		}
	}
	if (@{$be}) {
		if (scalar @{$be} <= 3) {
			$badE=join "、",@{$be};
		}
		else {
			$badE=join("、",$be->[0],$be->[1],$be->[2]);
			$badE.="等";
		}
	}
	if (@{$bd}) {
		if (scalar @{$bd} <= 3) {
			$badD=join "、",@{$bd};
		}
		else {
			$badD=join("、",$bd->[0],$bd->[1],$bd->[2]);
			$badD.="等";
		}
	}

	$badActR=join "、",@{$actrule};
	$badActM=join "、",@{$actmore};
	$badActL=join "、",@{$actless};

	$badFoodR=join "、",@{$foodrule};
	$badFoodM=join "、",@{$foodmore};
	$badFoodL=join "、",@{$foodless};

	$goodDe='-' if ($goodDe eq "");
	$goodE='-' if ($goodE eq "");
	$goodD='-' if ($goodD eq "");

	$badDe='-' if ($badDe eq "");
	$badE='-' if ($badE eq "");
	$badD='-' if ($badD eq "");
	$badActR='-' if ($badActR eq "");
	$badActM='-' if ($badActM eq "");
	$badActL='-' if ($badActL eq "");
	$badFoodR='-' if ($badFoodR eq "");
	$badFoodM='-' if ($badFoodM eq "");
	$badFoodL='-' if ($badFoodL eq "");

	return ($goodDe,$goodE,$goodD,$badDe,$badE,$badD,$badActR,$badActM,$badActL,$badFoodR,$badFoodM,$badFoodL);
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


#($vars->{xiaohuahexishou}{desc}{$prefix}{good},$vars->{xiaohuahexishou}{desc}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
sub process_good_bad {
	my ($goodA,$badA)=@_;
	my ($good,$bad)=('-','-');
	
	##
	if (@{$goodA}) {
		if (scalar @{$goodA} <= 3) {
			$good=join "、",@{$goodA};
		}
		else {
			$good=join("、",$goodA->[0],$goodA->[1],$goodA->[2]);
			$good.="等";
		}
	}
	##
	if (@{$badA}) {
		if (scalar @{$badA} <= 3) {
			$bad=join "、",@{$badA};
		}
		else {
			$bad=join("、",$badA->[0],$badA->[1],$badA->[2]);
			$bad.="等";
		}
	}
	
	return ($good,$bad);
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

#$item->{key}=&get_item_key($item);
sub get_item_key {
	my ($item)=@_;
	my $key="-";
	
	##
	if ($item->{name} =~ /蛋白/) {
		$key="蛋白质";
	}
	elsif ($item->{name} =~ /脂肪/) {
		$key="脂肪";
	}
	elsif ($item->{name} =~ /碳水化合物/) {
		$key="碳水化合物";
	}
	else {
	}
	
	return ($key);
}

#$vars->{changdaojunqun}{diseaseHighRule}=&process_HighRisk_disease_food(\@rule);
sub process_HighRisk_disease_food {
	my ($array)=@_;
	my $desc="-";
	
	##
	if (@{$array}) {
		my ($uniq)=&split_uniq_array($array);
		$desc=join "、", @{$uniq};
	}
	
	return ($desc);
}

#my ($item)=&get_item($vd,$Data->{$meta},$language,$meta);
sub get_item {
	my ($vd,$Data,$language,$key)=@_;
	my %item=();
	
	##
	$item{value}=sprintf "%.2f",$Data->{sample}*100;
	$item{value}=sprintf "%.2f",$Data->{sample} if ($key eq "zxjunqunduoyangxing");
	$item{per}=sprintf "%.4f",$Data->{per};
	if ($key eq "zxjunqunduoyangxing" || $key eq "zxyouyijun") {
		$item{per}=0.9999 if ($item{per} == 1);
		$item{per}=0.0001 if ($item{per} == 0);
		$item{value}=0.01 if ($item{value} == 0);
	}
	$item{level}=$Data->{level};
	$item{range}=$Data->{range};
	$item{pic}=$key."_dis.png";
	$item{coordinate}=$Data->{per}*4.98-5.66;
	$item{name}=$vd->{$language}->{title};
	$item{desc}=$vd->{$language}->{summary}->{desc};
	if ($Data->{level} ne '-') {
		$item{risk}=$vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{effect};
		if ($vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{mechanism} ne "") {
			$item{riskdesc}=$vd->{$language}->{suggestion}->{healthy}->{$Data->{level}}->{mechanism};
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

#&process_relative_abundance($datadir,$indir,\%vars);
sub process_relative_abundance {
	my ($datadir,$indir,$referdir,$confdir,$vars)=@_;
	
	##
	my $samFile="$datadir/$vars->{barcode}".".phylum.xls";
	my $refFile="$referdir/refer.phylum.stat.xls";
	my $config="$confdir/phylum_latin2zh";
	## read sample
	$/="\n";
	my %hash;
	open IN,$samFile or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$value)=(split /\t/,$_)[0,1];
		$hash{$id}{sam}=$value;
		$hash{$id}{ref}=0;
	}
	close IN;
	## read ref
	open IN,$refFile or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$value)=(split /\t/,$_)[0,3];
		$hash{$id}{ref}=$value;
		$hash{$id}{sam}=0 if (!exists $hash{$id}{sam});
	}
	close IN;
	## read config
	open IN,$config or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$name)=(split /\t/,$_)[0,1];
		$hash{$id}{name}=decode("UTF-8",$name) if (exists $hash{$id});
	}
	close IN;
	
	## get $vars->{abundance}
	&sort_abundance(\%hash,$vars);
	
	return;
}

#&sort_abundance(\%hash,$vars);
sub sort_abundance {
	my ($hash,$vars)=@_;
	
	##
	my $cnt=0;
	my ($start1,$start2)=(0,0);
	my @sample_phylum=();
	my @refer_phylum=();
	#foreach my $key (sort {$hash->{$b}{ref} <=> $hash->{$a}{ref}} keys %{$hash}) {
	foreach my $key (sort {$hash->{$b}{sam} <=> $hash->{$a}{sam}} keys %{$hash}) {
		next if ($cnt > 6);
		##
		my %item=();
		$item{ename}=$key;
		$item{name}=$key;
		$item{name}=$hash->{$key}{name} if (exists $hash->{$key}{name});
		$item{index}=$cnt;
		#
		#$item{per}{sample}=0;
		#$item{per}{sample}=$hash->{$key}{sam} if (defined $hash->{$key}{sam});
		$item{per}{sample}=$hash->{$key}{sam};
		push @sample_phylum,$hash->{$key}{sam};
		$item{per}{control}=$hash->{$key}{ref};
		push @refer_phylum,$hash->{$key}{ref};
		#
		$item{start}{sample}=$start1;
		$item{start}{control}=$start2;
		#
		$item{end}{sample}=$item{start}{sample}+$item{per}{sample};
		$item{end}{control}=$item{start}{control}+$item{per}{control};
		#
		$start1=$item{end}{sample};
		$start2=$item{end}{control};
		#
		$cnt++;
		push @{$vars->{abundance}->{item}},\%item;
	}
	$vars->{abundance}->{cnt}=$cnt;

	## ger pe
	$R->set('y', \@sample_phylum);
	$R->set('z', \@refer_phylum);
	$R->run(q`pearson_cor <- cor(y, z)`);
	$vars->{abundance}->{pearson_cor} = sprintf "%.2f",$R->get('pearson_cor');

	return;
}

#&get_disease_high_middle_low($vars,$disease{name},$disXML->{$language},$level,\@diseaseHigh,\@diseaseLow,\@diseaseNormal,\@actrule,\@actmore,\@actless,\@foodrule,\@foodmore,\@foodless);
sub get_disease_high_middle_low {
	my ($vars,$name,$hash,$level,$diseaseHigh,$diseaseLow,$diseaseNormal,$actrule,$actmore,$actless,$foodrule,$foodmore,$foodless)=@_;
	
	##
	if ($hash->{suggestion}->{suggest}->{$level}->{compare} eq "一致") {
		push @{$diseaseNormal},$name;
		$vars->{changdaojunqun}{diseaseNormalN}++;
	}
	elsif ($hash->{suggestion}->{suggest}->{$level}->{compare} eq "偏高") {
		push @{$diseaseHigh},$name;
		$vars->{changdaojunqun}{diseaseHighN}++;
		#
		push @{$actrule},$hash->{healthcarerule}->{description} if ($hash->{healthcarerule}->{description} ne "");
		push @{$actmore},$hash->{healthcarerule}->{recommend} if ($hash->{healthcarerule}->{recommend} ne "");
		push @{$actless},$hash->{healthcarerule}->{taboo} if ($hash->{healthcarerule}->{taboo} ne "");
		push @{$foodrule},$hash->{nutritionrule}->{description} if ($hash->{nutritionrule}->{description} ne "");
		push @{$foodmore},$hash->{nutritionrule}->{more} if ($hash->{nutritionrule}->{more} ne "");
		push @{$foodless},$hash->{nutritionrule}->{less} if ($hash->{nutritionrule}->{less} ne "");
	}
	elsif ($hash->{suggestion}->{suggest}->{$level}->{compare} eq "偏低") {
		push @{$diseaseLow},$name;
		$vars->{changdaojunqun}{diseaseLowN}++;
	}
	
	return;
}

#my ($level)=&process_risk_protect($risk,$protect);
sub process_risk_protect {
	my ($risk,$protect)=@_;
	my $level='-';
	
	##
	$risk=&transfer_level($risk);
	$protect=&transfer_level($protect);
	##
	$level=&get_disease_level($risk,$protect);
	
	return ($level);
}

#$level=&get_disease_level($risk,$protect);
sub get_disease_level {
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

#$risk=&transfer_level($risk);
sub transfer_level {
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
	#
	$file="$referdir/refer.phylum.stat.xls";
	$check+=&checkF($file);
	
	## check config file
	#
	$file="$confdir/phylum_latin2zh";
	$check+=&checkF($file);
	#	
	$file="$confdir/genus_latin2zh";
	$check+=&checkF($file);
	#
	$file="$confdir/zx_genus_fenbu.conf";
	$check+=&checkF($file);
	#
	$file="$confdir/zx_Vitamin.conf";
	$check+=&checkF($file);
	
	## check input data file
	#
	$file="$datadir/$barcode".".extract.xls";
	$check+=&checkF($file);
	#
	$file="$datadir/$barcode".".phylum.xls";
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

#&get_sport($vd->{$language}->{sport},$vars);
sub get_sport {
	my ($hash,$vars)=@_;
	
	##
	$vars->{fangan}{sport}{name}=$hash->{title};
	$vars->{fangan}{sport}{desc}=$hash->{desc};
	$vars->{fangan}{sport}{pic}="yundong.pdf";
	
	return;
}

#&get_food($vd->{$language}->{food},$vars);
sub get_food {
	my ($hash,$vars)=@_;
	
	## 主食
	push @{$vars->{fangan}{food}{items}{item}},{name => "主食", pic => "zhushi.pdf", desc => $hash->{staple}};
	## 肉类
	push @{$vars->{fangan}{food}{items}{item}},{name => "肉类", pic => "roulei.pdf", desc => $hash->{meat}};
	## 鱼虾类
	push @{$vars->{fangan}{food}{items}{item}},{name => "鱼虾类", pic => "yuxialei.pdf", desc => $hash->{fish}};
	## 蛋奶、大豆类
	push @{$vars->{fangan}{food}{items}{item}},{name => "蛋奶、大豆类", pic => "danlei.pdf", desc => $hash->{egg}};
	## 蔬菜水果
	push @{$vars->{fangan}{food}{items}{item}},{name => "蔬菜水果", pic => "shucaishuiguo.pdf", desc => $hash->{fruit}};
	## 食用油、盐
	push @{$vars->{fangan}{food}{items}{item}},{name => "食用油、盐", pic => "you.pdf", desc => $hash->{oil}};
	## 其他
	push @{$vars->{fangan}{food}{items}{item}},{name => "其他", pic => "other.pdf", desc => $hash->{other}};
	## 备注
	$vars->{fangan}{food}{remark}=$hash->{remark};
	
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

#($diseaseF,$diseaseR)=&get_program_disease($samInfo,$listQ,$listP);
sub get_program_disease {
	my ($diseaseH,$listQ,$listP)=@_;
	my ($diseaseF,$diseaseR)=('-','-');
	
	#### food
	## sample information
	foreach my $k (sort {$listQ->{$a} <=> $listQ->{$b}} keys %{$listQ}) {
		if ($diseaseH->{$k}{risk} eq 'Y') {
			$diseaseF=$k;
			$diseaseR=$k;
			last;
		}
	}

	return ($diseaseF,$diseaseR);
}

#($badDiseaseF,$badDiseaseR)=&get_bad_disease_meta($vars->{changdaojunqun}{diseaseHigh},$badDiseaseTotal,$samInfo->{disease});
sub get_bad_disease_meta {
	my ($diseaseHigh,$hash,$listP)=@_;
	my ($diseaseF,$diseaseR)=('-','-');

	### high risk disease
	my %diseaseHigh_pre;
	if ($diseaseHigh ne '-') {
		my @d=split /\、/,$diseaseHigh;
		foreach (@d) {
			$diseaseHigh_pre{$_} = 1;
		}
	}
	foreach my $k (sort {$listP->{$a} <=> $listP->{$b}} keys %{$listP}) {
		if (exists $diseaseHigh_pre{$k}) {
			$diseaseF=$k;
			$diseaseR=$k;
			last;
		}
	}

	# disease risk of items
	#if ($diseaseF eq '-') {
	#	foreach my $k (sort {$listP->{$a} <=> $listP->{$b}} keys %{$listP}) {
	#		if (exists $hash->{$k}) {
	#			$diseaseF=$k;
	#			$diseaseR=$k;
	#			last;
	#		}
	#	}
	#}
	
	return ($diseaseF,$diseaseR);
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
		my ($oriid,$en,$cn)=(split /\s+/,$_)[0,1,2];
		$hash{$oriid}=1;
	}
	close IN;

	return (\%hash);
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
				$head[$i]="炎症性肠病" if ($head[$i] =~ /肠炎/ || $head[$i] =~ /炎症性肠病/);
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
	## overview.list.order$sex
	$file="$indir/overview.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"overview");
	}
	## distribution.list.order$sex
	$file="$indir/distribution.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"distribution");
	}
	## pathogens.list.order$sex
	$file="$indir/pathogens.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"pathogens");
	}
	## absorption.list.order$sex
	$file="$indir/absorption.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"absorption");
	}
	## metabolism.list.order$sex
	$file="$indir/metabolism.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"metabolism");
	}
	## giimmune.list.order$sex
	$file="$indir/giimmune.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"giimmune");
	}
	## diseasefactor.list.order$sex
	$file="$indir/diseasefactor.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"diseasefactor");
	}
	## stool.list.order$sex
	$file="$indir/stool.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"stool");
	}

	## disease.list.order$sex
	$file="$indir/disease.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"disease");
	}
	## age.list.order$sex
	$file="$indir/age.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"age");
	}
	## dysbiosis.list.order$sex
	$file="$indir/dysbiosis.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"dysbiosis");
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
			$listH->{$type}{$section1}{$j}{name}=$section2;
			for (my $i=1;$i<@lines;$i++) {
				my @aa=split /\t/,$lines[$i];
				my ($oriid,$en,$cn);
				$oriid=$aa[0];
				$cn=$aa[-1];
				$cn=decode('UTF-8',$cn);
				if (scalar @aa == 3) {
					$en=$aa[1];
				}
				$listH->{$type}{$section1}{$j}{section}{$oriid}{order}=$i;
				$listH->{$type}{$section1}{$j}{section}{$oriid}{en}=$en if (defined $en);
				$listH->{$type}{$section1}{$j}{section}{$oriid}{cn}=$cn;
				$listH->{diseaseList}{$cn}=$oriid if ($type eq "disease");
				$listH->{stoolList}{$cn}=$oriid if ($type eq "stool");
				if (defined $en) {
					$term2oriid->{$en}{$oriid}=1;
				}
			}
		}
	}
	$/="\n";
	close IN;
	
	return;
}

#my ($Data)=&calculate_meta_level_per($datadir,$indir,$barcode);
sub calculate_meta_level_per {
	my ($datadir,$indir,$referdir,$term2oriid,$barcode,$level_set)=@_;
	my $Data=();
	
	## load sample and reference data
	($Data)=&load_sample_reference_data($datadir,$indir,$referdir,$barcode,$term2oriid);
	
	## get level,per,range
	&calculate_level_per($Data,$level_set);
	
	return ($Data);
}

#&calculate_level_per($Data)
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

#$Data->{$meta}{per}=&get_per($Data->{$meta});
sub get_per {
	my ($meta, $key)=@_;
	my $per='-';
	
	##
	if ((defined $meta->{ref}{mean} && $meta->{ref}{mean} != 0) && (defined $meta->{ref}{sd} && $meta->{ref}{sd} != 0) && exists $meta->{sample}) {
		## use R pnorm
		$R->set('v', $meta->{sample});
		$R->set('y', $meta->{ref}{mean});
		$R->set('z', $meta->{ref}{sd});
		$R->run(q`prob <- pnorm(v, y, z)`);
		$per = $R->get('prob');
	}

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

#($Data->{$meta}{level},$Data->{$meta}{range})=&get_level($Data->{$meta});
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
		$range=(sprintf "%.2f",$mean-2*$sd).'~'.(sprintf "%.2f",$mean+2*$sd) if ($k eq "zxjunqunduoyangxing");
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

#($Data)=&load_sample_reference_data($datadir,$indir,$barcode);
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

#&get_sample_values($file,\%Data);
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
		foreach my $key (keys %{$term2oriid->{$id}}) {
			$id = $key;
			$Data->{$id}{sample}=$value;
		}
	}
	close IN;

	return;
}

#&get_ref_mean_sd($file,\%Data);
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
		foreach my $key (keys %{$term2oriid->{$unit[0]}}) {
			$unit[0] = $key;
			$Data->{$unit[0]}{ref}{mean}=$unit[3];
			$Data->{$unit[0]}{ref}{sd}=$unit[4];
		}
	}
	close IN;
	
	return;
}

#&get_ref_values($file,\%Data);
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
		foreach my $key (keys %{$term2oriid->{$unit[0]}}) {
			$unit[0] = $key;
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
	}
	close IN;
	
	return;
}

#&rename_disease_factor_name($item);
sub rename_disease_factor_name {
	my ($item)=@_;
	
	##
	if ($item->{name} =~ /保护因子/) {
		$item->{name}="保护因子";
	}
	elsif ($item->{name} =~ /危险因子/) {
		$item->{name}="危险因子";
	}
	
	return;
}

#my ($diseaseCNtitle2meta)=&get_disease_CN_title_meta($list->{disease}{'相关疾病风险'});
sub get_disease_CN_title_meta {
	my ($list)=@_;
	my %title2meta=();
	
	##
	foreach my $num (keys %{$list}) {
		foreach my $meta (keys %{$list->{$num}{section}}) {
			$title2meta{$list->{$num}{section}{$meta}{cn}}=$meta;
		}
	}
	
	return (\%title2meta);
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
{		#获取指定目录或文件的决定路径
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
      -confdir    <dir>   report conf dir                           optional
                          default "$Bin/../conf/current";
      -historydir <dir>   history test xml result for report dir    optional
                          default "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
      -out        <dir>   output file dir                           optional
                          default "./";
      -language   <str>   result language used for report           optional
                          default "CN";
      -l          <str>   level set to check report correctness (l1~l5)          optional
                          default "N";

      -h          Help

USAGE
	print $usage;
	exit;
}
