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
use JSON::XS;
use MongoDB;

use Statistics::R;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my $version = "1.0.0";
my $R = Statistics::R->new();

###
my ($datadir,$indir,$referdir,$dirout,$language,$barcode,$historydir,$confdir,$level_set,$host,$port,$database,$collection,$help);
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

	"host=s"	=> \$host,
	"port=s"	=> \$port,
	"database=s"	=> \$database,
	"collection=s"	=> \$collection,
	"h"	=> \$help,
) or &USAGE;
&USAGE if ((!defined $datadir || !defined $barcode) || defined $help) ;

########## AbsolutePath and default ##############
$dirout||="./";
$indir||="$Bin/db2xml/current";
$referdir||="$Bin/../../reference_data/qpcr_current";
$confdir||="$Bin/../conf/current";
$historydir||="/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
mkdir $dirout if (! -d $dirout);
$dirout=&AbsolutePath("dir",$dirout);
$datadir=&AbsolutePath("dir",$datadir);
$referdir=&AbsolutePath("dir",$referdir);
$dirout=&AbsolutePath("dir",$dirout);
$historydir=&AbsolutePath("dir",$historydir);
$confdir=&AbsolutePath("dir",$confdir);

$language||="CN";

$level_set||="N";

$host = $host || "10.0.0.8";
$port = $port || "27021";
$database = $database || "result";
$collection = $collection || "tpqpcr";

################# main process ###################
my %vars=();
## check file and dir
&check_file($datadir,$indir,$referdir,$confdir,$barcode);

## load sample information file
my ($samInfo)=&load_sample_info($barcode,$datadir,\%vars);

## load $type.list.order$sex
my ($listH,$term2oriid,$overviewcomposition,$compositionweigth)=&load_list($indir,$samInfo,$confdir);

## calculate meta level and per by sample data and ref data
my ($Data)=&calculate_meta_level_per($datadir,$indir,$referdir,$barcode,$term2oriid,$compositionweigth,$level_set);
&calculate_overview_level_per($Data,$overviewcomposition,$compositionweigth,$level_set);

## main process
#&main_process($Data,$listH,\%{$diseaseOrder{question}},\%{$diseaseOrder{predict}},\%vars,$datadir,$indir,$language,$barcode,$samInfo);
&main_process($Data,$listH,\%vars,$datadir,$indir,$language,$barcode,$samInfo,$historydir);

## sample result csv for erp import
&sample_result_csv(\%vars,$listH->{pathogens},$barcode,$dirout);
my %final_report;
&get_final_report(\%vars,\%final_report,$barcode,$dirout);

## insert final report to mongo db
#my $mongo = MongoDB::MongoClient->new('host' => $host.':'.$port);
#my $db = $mongo->get_database($database);
#my $cl = $db->get_collection($collection);
#my $barcode_find = $cl->find({_id => $barcode});
#
##delete $final_report{page3};
#my $flag = 0;
#while (my $barcode_report = $barcode_find->next()) {
#	$flag = 1;
#	if (!exists $barcode_report->{Report}{page3} && exists $final_report{page3}) {
#		$cl->update({_id => $barcode}, {'$set' => {'Report.page3' => \%{$final_report{page3}} }}, {'upsert' => 1});
#	}
#}
#
#if ($flag == 0) {
#	$cl->update({_id => $barcode}, {'$set' => {'Report' => \%final_report }},{'upsert' => 1});
#}


############################# sub ###############################################
sub get_final_report {
	my ($vars, $final_report,$barcode,$dirout) = @_;

	$final_report->{barcode} = $vars->{barcode};
	$final_report->{name} = $vars->{name};
	$final_report->{age} = $vars->{age};
	$final_report->{receivedate} = $vars->{receivedate};
	$final_report->{sex} = $vars->{sex};
	$final_report->{source} = $vars->{source};

	$final_report->{score} = $vars->{gaikuang}{desc}{youyijun}{per100};
	$final_report->{page1}{part1} = "$vars->{gaikuang}{desc}{youyijun}{actionrule}";
	$final_report->{page1}{part1} .= "益生菌检测总评分越高，说明您的肠道越健康；反之，益生菌检测总评分越低，则说明您的肠道微生态环境越值得关注。";

	foreach my $item (@{$vars{changdaojunqun}{fenbu}{item}}) {
		my %new_item;
		$new_item{name} = $item->{name};
		$item->{name} =~ /(.*)\（([^\）]+)\）/;
		$new_item{cnname} = $1;
		$new_item{ename} = $2;
		$new_item{value} = $item->{value};
		$new_item{threshold} = $item->{threshold};
		$new_item{level} = $item->{level};
		$new_item{risk} = $item->{risk};

		push @{$final_report->{fenbu}{item}}, \%new_item;
	}

	if ($vars->{changdaojunqun}{fenbu}{badN} == 0) {
		$final_report->{page2}{part1} = "恭喜您！您本次检测的4种主要益生菌的含量均处于正常范围内或优于正常范围，其中乳杆菌和双歧杆菌是调节肠道微生态平衡、抵御有害菌、促进营养物质吸收的益生菌，";
		$final_report->{page2}{part1} .= "对维持免疫调节力等都有重要作用；阿克曼氏菌是能够帮助控制体重的“瘦菌”，能促进脂肪代谢，消耗多余能量，维持代谢平衡；柔嫩梭菌是健康人肠道中含量比较高的细菌，能产生大量对人体非常有益的丁酸，";
		$final_report->{page2}{part1} .= "还能够调节人体免疫系统、抑制炎症、调节肠道激素分泌以及人体的代谢平衡。";
	}
	elsif ($vars->{changdaojunqun}{fenbu}{badN} < 4) {
		$final_report->{page2}{part1} = "您本次检测的4种主要益生菌中"."$vars->{changdaojunqun}{fenbu}{bad}"."含量相比参考人群偏低。"."$vars->{changdaojunqun}{fenbu}{badSource}";
	}
	else {
		$final_report->{page2}{part1} = "您本次检测的4种主要益生菌含量相比参考人群均偏低。"."$vars->{changdaojunqun}{fenbu}{badSource}";
	}

	if ($vars->{changdaojunqun}{fenbu}{badEffect} ne '-') {
		$final_report->{page2}{part1} .= "这些菌属含量偏低，有可能导致"."$vars->{changdaojunqun}{fenbu}{badEffect}"."等健康问题";
	}
	if ($vars->{changdaojunqun}{fenbu}{badDisease} ne '-') {
		$final_report->{page2}{part1} .= "，甚至增加"."$vars->{changdaojunqun}{fenbu}{badDisease}"."的风险";
	}
	$final_report->{page2}{part1} .= "。";

	$final_report->{page3}{part1} = "综合您本次肠道益生菌检测，您本次益生菌检测总评分为"."$final_report->{score}"."分，";
	if ($vars->{gaikuang}{desc}{youyijun}{risk} eq "有害") {
		$final_report->{page3}{part1} .= "不利于肠道及人体健康，";
		if ($vars->{changdaojunqun}{fenbu}{badN} == 4) {
			$final_report->{page3}{part1} .= "4种益生菌属的检测指数均偏低，需注意。";
		}
		elsif ($vars->{changdaojunqun}{fenbu}{badN} > 0) {
			$final_report->{page3}{part1} .= "其中，"."$vars->{changdaojunqun}{fenbu}{bad}"."检测指数偏低，需注意。";
		}
	}
	elsif ($vars->{gaikuang}{desc}{youyijun}{risk} eq "正常") {
		$final_report->{page3}{part1} .= "处于正常范围，";
		if ($vars->{changdaojunqun}{fenbu}{badN} == 4) {
			$final_report->{page3}{part1} .= "但4种益生菌属的检测指数均偏低，需注意。";
		}
		elsif ($vars->{changdaojunqun}{fenbu}{badN} > 0) {
			$final_report->{page3}{part1} .= "但其中，"."$vars->{changdaojunqun}{fenbu}{bad}"."检测指数偏低，需注意。";
		}
	}
	else {
		$final_report->{page3}{part1} .= "有利于肠道及人体健康，";
		if ($vars->{changdaojunqun}{fenbu}{badN} == 4) {
			$final_report->{page3}{part1} .= "但4种益生菌属的检测指数均偏低，需注意。";
		}
		elsif ($vars->{changdaojunqun}{fenbu}{badN} > 0) {
			$final_report->{page3}{part1} .= "但其中，"."$vars->{changdaojunqun}{fenbu}{bad}"."检测指数偏低，可适当改善。";
		}
	}
	$final_report->{page3}{part1} .= "建议您：";

	$final_report->{page3}{fangan}{food} = $vars->{fangan}{food};
	$final_report->{page3}{fangan}{prebiotics} = $vars->{fangan}{prebiotics};
	$final_report->{page3}{fangan}{probiotics} = $vars->{fangan}{probiotics};
	#$final_report->{page3}{fangan}{others} = $vars->{fangan}{others};

#	if ($vars->{age} <= 18) {
#		$final_report->{page3}{fangan}{regulate}{head} = "肠道康，更健康！补充适量的益生菌、益生元，助力健康成长。";
#	}
#	else {
#		$final_report->{page3}{fangan}{regulate}{head} = "肠命百岁方可长命百岁！补充适量的益生菌、益生元，助力肠道健康。";
#	}

}

#&main_process($Data,$listH,\%vars,$datadir,$indir,$language,$barcode,$samInfo);
sub main_process {
	my ($Data,$listH,$vars,$datadir,$indir,$language,$barcode,$samInfo,$historydir)=@_;

	## 菌群分布(有益菌属分布)
	# $listH->{"distribution"}; $listH->{"pathogens"}
	&Intestinal_Flora($Data,$listH->{distribution},$vars,"$indir/XML/distribution/$language",$language);

	## 菌群综合指标(有益菌)
	# $listH->{"overview"}
	&Micro_Overvwiew($Data,$listH->{overview},$vars,"$indir/XML/overview/$language",$language);

	## 膳食方案+肠道调节方案+运动方案
	# $listH->{"advise"};
	&Program($vars,"$indir/XML/advise",$language,$samInfo,$listH->{advise});

	## 多次检测结果比较
	# $samInfo->{batch}; $samInfo->{'历次检测编号'}
	# &Compare_history($vars,$samInfo,$historydir) if ($samInfo->{batch} > 1);

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
	open my $out, '>:encoding(utf-8)', $outfile || die $!;

	my $term_str = "检测编号\t姓名\t益生菌\t";
	my $result_str = "$vars->{barcode}\t$vars->{name}\t";

	my $level;
	$level = &transfer_item_level($vars->{gaikuang}{desc}{youyijun}{level});
	$result_str .="$level($vars->{gaikuang}{desc}{youyijun}{risk})\t";

	foreach my $item (@{$vars{changdaojunqun}{fenbu}{page1}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
	}
	foreach my $item (@{$vars{changdaojunqun}{fenbu}{page2}}) {
		$term_str .= "$item->{name}\t";
		$level = &transfer_item_level($item->{level});
		$result_str .="$level($item->{risk})\t";
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
	# &process_regulate($vars,$dir,$ageKey,$samInfo->{sex},$language);
	## process probiotics
	&process_probiotics($vars,$dir,$ageKey,$samInfo->{sex},$language);
	## process prebiotics
	&process_prebiotics($vars,$dir,$ageKey,$samInfo->{sex},$language);
	
	return;
}

#&process_regulate($vars,$dir,$regulateKey,$language);
sub process_regulate {
	my ($vars,$dir,$ageKey,$sex,$language)=@_;
	
	## 补充益生元和益生菌
	my $vd=XMLin("$dir/$language/qpcrjiankangchangdao.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{regulate}{name} = $vd->{$language}{title};
	my @actdescs = split/\n/,$vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
	my $num = 0;
	foreach my $actdesc (@actdescs) {
		$actdesc =~ /(^l\d+)：(.*)/;
		if ($vars->{gaikuang}{desc}{youyijun}{level} eq $1) {
			$num ++;
			$vars->{fangan}{regulate}{desc}{item}{desc} = $2;
			$vars->{fangan}{regulate}{desc}{item}{num} = $num;
		}
	}
	#$vars->{fangan}{regulate}{desc} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};

	return;
}

#&process_prebiotics($vars,$dir,$regulateKey,$language);
sub process_prebiotics {
	my ($vars,$dir,$ageKey,$sex,$language)=@_;
	
	## 补充益生元
	my $vd=XMLin("$dir/$language/tjsjcyishengyuanfangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{prebiotics}{name} = $vd->{$language}{title};
	$vars->{fangan}{prebiotics}{suggest1} = '-';
	$vars->{fangan}{prebiotics}{suggest2} = '-';
	my %prebiotics_sup;
	if ($vars->{changdaojunqun}{fenbu}{prebiotics_sup} eq '-') {
		$prebiotics_sup{num} = 1;
		$prebiotics_sup{desc} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
		$vars->{fangan}{prebiotics}{suggest1} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
	}
	else {
		$prebiotics_sup{num} = 1;
		$prebiotics_sup{desc} = "建议您选择补充"."$vars->{changdaojunqun}{fenbu}{prebiotics_sup}"."等益生元，能帮助"."$vars->{changdaojunqun}{fenbu}{biotics_bad}"."在肠道内增殖，增加有益物质的产生，有利于增强肠道屏障功能，促进肠道健康。";
		$vars->{fangan}{prebiotics}{suggest1} = "建议您选择补充"."$vars->{changdaojunqun}{fenbu}{prebiotics_sup}"."等益生元，能帮助"."$vars->{changdaojunqun}{fenbu}{biotics_bad}"."在肠道内增殖，增加有益物质的产生，有利于增强肠道屏障功能，促进肠道健康。";
	}
	push @{$vars->{fangan}{prebiotics}{desc}{item}}, \%prebiotics_sup;

	if ($vd->{$language}{agesuggestion}{$sex}{$ageKey}{note} ne '') {
		my %prebiotics_note;
		$prebiotics_note{num} = 2;
		$prebiotics_note{desc} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{note};
		$vars->{fangan}{prebiotics}{suggest2} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{note};

		push @{$vars->{fangan}{prebiotics}{desc}{item}}, \%prebiotics_note;
	}
#	my @actdescs = split/\n/,$vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	my $num = 0;
#	foreach my $actdesc (@actdescs) {
#		$actdesc =~ /(^l\d+)：(.*)/;
#		if ($vars->{gaikuang}{desc}{youyijun}{level} eq $1) {
#			$num ++;
#			$vars->{fangan}{prebiotics}{desc}{item}{desc} = $2;
#			$vars->{fangan}{prebiotics}{desc}{item}{num} = $num;
#		}
#	}

	return;
}

#&process_probiotics($vars,$dir,$regulateKey,$language);
sub process_probiotics {
	my ($vars,$dir,$ageKey,$sex,$language)=@_;
	
	## 补充益生菌
	my $vd=XMLin("$dir/$language/tjsjcyishengjunfangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{probiotics}{name} = $vd->{$language}{title};
	$vars->{fangan}{probiotics}{suggest1} = '-';
	$vars->{fangan}{probiotics}{suggest2} = '-';
	my %probiotics_sup;
	if ($vars->{changdaojunqun}{fenbu}{probiotics_sup} eq '-') {
		$probiotics_sup{num} = 1;
		$probiotics_sup{desc} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
		$vars->{fangan}{probiotics}{suggest1} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
	}
	else {
		$probiotics_sup{num} = 1;
		$probiotics_sup{desc} = "建议您选择补充"."$vars->{changdaojunqun}{fenbu}{probiotics_sup}"."等益生菌产品，以增加您肠道内益生菌的含量，并能抑制有害菌的异常增殖，调节肠道菌群平衡。";
		$vars->{fangan}{probiotics}{suggest1} = "建议您选择补充"."$vars->{changdaojunqun}{fenbu}{probiotics_sup}"."等益生菌产品，以增加您肠道内益生菌的含量，并能抑制有害菌的异常增殖，调节肠道菌群平衡。";
	}
	push @{$vars->{fangan}{probiotics}{desc}{item}}, \%probiotics_sup;

	if ($vd->{$language}{agesuggestion}{$sex}{$ageKey}{note} ne '') {
		my %probiotics_note;
		$probiotics_note{num} = 2;
		$probiotics_note{desc} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{note};
		$vars->{fangan}{probiotics}{suggest2} = $vd->{$language}{agesuggestion}{$sex}{$ageKey}{note};

		push @{$vars->{fangan}{probiotics}{desc}{item}}, \%probiotics_note;
	}
#	my @actdescs = split/\n/,$vd->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	my $num = 0;
#	foreach my $actdesc (@actdescs) {
#		$actdesc =~ /(^l\d+)：(.*)/;
#		if ($vars->{gaikuang}{desc}{youyijun}{level} eq $1) {
#			$num ++;
#			$vars->{fangan}{probiotics}{desc}{item}{desc} = $2;
#			$vars->{fangan}{probiotics}{desc}{item}{num} = $num;
#		}
#	}

	return;
}

#&process_food_sport($vars,$dir,$diseaseF,$ageF,$badDiseaseF,$language);
sub process_food {
	my ($vars,$dir,$ageKey,$sex,$language)=@_;
	
	## food
	my $vd_food=XMLin("$dir/$language/tjsjcshanshifangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vars->{fangan}{food}{name} = $vd_food->{$language}{title};
	$vars->{fangan}{food}{suggest1} = '-';
	$vars->{fangan}{food}{suggest2} = '-';
	my %food_sup;
	if ($vars->{changdaojunqun}{fenbu}{food_sup} eq '-') {
		$food_sup{num} = 1;
		$food_sup{desc} = $vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{action};
		$vars->{fangan}{food}{suggest1} = $vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{action};
	}
	else {
		$food_sup{num} = 1;
		$food_sup{desc} = "建议您在日常饮食中选择补充"."$vars->{changdaojunqun}{fenbu}{food_sup}"."等食物或者相关的加工提取食品，为肠道菌群提供更多营养，促进有益共生菌的生长，增加肠道菌群的丰富度及多样性。";
		$vars->{fangan}{food}{suggest1} = "建议您在日常饮食中选择补充"."$vars->{changdaojunqun}{fenbu}{food_sup}"."等食物或者相关的加工提取食品，为肠道菌群提供更多营养，促进有益共生菌的生长，增加肠道菌群的丰富度及多样性。";
	}
	push @{$vars->{fangan}{food}{desc}{item}}, \%food_sup;

	if ($vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{note} ne '') {
		my %food_note;
		$food_note{num} = 2;
		$food_note{desc} = $vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{note};
		$vars->{fangan}{food}{suggest2} = $vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{note};

		push @{$vars->{fangan}{food}{desc}{item}}, \%food_note;
	}

#	#$vars->{fangan}{food}{desc} = $vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	my @food_notes = split/\n/,$vd_food->{$language}{agesuggestion}{$sex}{$ageKey}{action};
#	#$vars->{fangan}{food}{head} = shift @food_notes;
#	my $food_note_index = 0;
#	foreach my $note (@food_notes) {
#		$note=~s/^\s+|\s+$//;
#		next if ($note=~/^$/);
#		$food_note_index ++;
#		#my ($index_name,$index_desc) = split/：/,$note;
#		#$vars->{fangan}{others}{"desc$note_index"} = $note;
#		#push @{$vars->{fangan}{food}{desc}{item}}, {name => $index_name, desc => $index_desc, num => $food_note_index};
#		push @{$vars->{fangan}{food}{desc}{item}}, {desc => $note, num => $food_note_index};
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
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta,$vars->{age});
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

#&Intestinal_Flora($Data,$listH->{meta}{'肠道菌群'},$vars,"$indir/XML/meta/$language",$language,$species,$summary);
sub Intestinal_Flora {
	my ($Data,$list,$vars,$dir,$language)=@_;
	
	##
	$vars->{changdaojunqun}{name}="肠道菌群";
	$vars->{changdaojunqun}{ename}="MICROBIOTA";

	## 分布
	&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language);
	
	return;
}

#&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language);
sub item_Intestinal_Flora {
	my ($prefix,$Data,$list,$vars,$dir,$language)=@_;
	
	##
	my (@good,@bad,@gooditem,@baditem);
	my (@goodDesc,@badDesc,@goodDisease,@goodSource,@badDisease,@goodEffect,@badEffect,@badSource);
	my (@biotics_good,@biotics_bad,@probiotics_sup,@prebiotics_sup);
	my (@food_sup_all,@food_sup_bad1,@food_sup_bad2,@food_sup_bad3);
	##
	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;
	#my $pageindex = 0;
	#my $itemcount = 1;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});

		## Filt if the bacteria genus percent is zero
		# next if ($Data->{$meta}{sample} == 0);

		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		#$pageindex ++ if ($itemcount % 2 == 1);
		#$itemcount ++;
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta,$vars->{age});
		if ($item->{name} =~ /(.*)（(.*)）/) {
			$item->{cnname} = $1;
			$item->{enname} = $2;
		}
		push @{$vars->{changdaojunqun}{$prefix}{item}},$item;
		#push @{$vars->{changdaojunqun}{$prefix}{"page$pageindex"}},$item;

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
			$item->{name}=~/(.*)\（.*/;
			unshift @biotics_bad,$1;
			&process_PreProbiotics_and_food($vars->{age},\@prebiotics_sup,\@probiotics_sup,\@food_sup_all,\@food_sup_bad1,\@food_sup_bad2,\@food_sup_bad3,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc}) if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc} ne "");
		}
		if (defined $item->{risk} && $item->{risk} eq "有益") {
			push @gooditem,$item;
			### Get risk disease info
			push @goodSource,$vd->{$language}{summary}{source} if (defined $vd->{$language}{summary}{source} && $vd->{$language}{summary}{source} ne "");
			push @goodDesc,$vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{mechanism} ne "");
			push @goodDisease,$vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{aiddisease} ne "");
			push @goodEffect,$vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} if (defined $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} && $vd->{$language}{suggestion}{healthy}{$item->{level}}{effectdesc} ne "");
			$item->{name}=~/(.*)\（.*/;
			unshift @biotics_good,$1;
		}
	}

	my $prebiotics_sup_uniq=&uniq_array(\@prebiotics_sup);
	my $probiotics_sup_uniq=&uniq_array(\@probiotics_sup);
	#my $food_sup_uniq=&uniq_array(\@food_sup_all) if ($vars->{changdaojunqun}{$prefix}{badN} <= 1);
	my @food_sup;
	@food_sup = (@food_sup_bad3, @food_sup_all) if ($vars->{changdaojunqun}{$prefix}{badN} <= 2);
	@food_sup = (@food_sup_bad2, @food_sup_all) if ($vars->{changdaojunqun}{$prefix}{badN} == 3);
	@food_sup = (@food_sup_bad1, @food_sup_all) if ($vars->{changdaojunqun}{$prefix}{badN} == 4);
	my $food_sup_uniq = &uniq_array(\@food_sup);
	#$food_sup_uniq=&uniq_array(\@food_sup_bad1) if ($vars->{changdaojunqun}{$prefix}{badN} >= 3);
	#$food_sup_uniq=&uniq_array(\@food_sup_bad2) if ($vars->{changdaojunqun}{$prefix}{badN} == 2);
	#$food_sup_uniq=&uniq_array(\@food_sup_bad3) if ($vars->{changdaojunqun}{$prefix}{badN} <= 1);

	##
	push @{$vars->{changdaojunqun}{$prefix}{abnormalitem}},@baditem,@gooditem;
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

#&process_Pre_and_Probiotics(\@prebiotics_sup,\@probiotics_sup,$vd->{$language}{suggestion}{$item->{level}}{actiondesc})
sub process_PreProbiotics_and_food {
	my ($age,$probiotics_sup,$prebiotics_sup,$food_sup_all,$food_sup_bad1,$food_sup_bad2,$food_sup_bad3,$Actdesc)=@_;

	my @Actarray = split/\n/,$Actdesc;
	foreach my $Act (@Actarray) {
		$Act =~s/^\s+|\s+$//;
		next if ($Act eq "");
		if ($Act =~ /(\d+)\\textasciitilde (\d+)）：/) {
			if ($age >= $1 && $age <= $2) {
				if ($Act =~/^益生元/) {
					$Act=~/^[^：]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$probiotics_sup},@sups;
				}
				if ($Act =~/^益生菌/) {
					$Act=~/^[^：]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$prebiotics_sup},@sups;
				}
				if ($Act =~/^食物/) {
					$Act=~/^[^：]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$food_sup_all},@sups;
					my $bad1 = &randomElem ('1', \@sups);
					push @{$food_sup_bad1}, @{$bad1};
					my $bad2 = &randomElem ('2', \@sups);
					push @{$food_sup_bad2}, @{$bad2};
					my $bad3 = &randomElem ('3', \@sups);
					push @{$food_sup_bad3}, @{$bad3};
				}
			}
		}
		elsif ($Act =~ /(\d+)\+）：/) {
			if ($age >= $1) {
				if ($Act =~/^益生元/) {
					$Act=~/^[^：]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$probiotics_sup},@sups;
				}
				if ($Act =~/^益生菌/) {
					$Act=~/^[^：]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$prebiotics_sup},@sups;
				}
				if ($Act =~/^食物/) {
					$Act=~/^[^：]+：(.*)/;
					my @sups=split/；/,$1;
					unshift @{$food_sup_all},@sups;
					my $bad1 = &randomElem ('1', \@sups);
					push @{$food_sup_bad1}, @{$bad1};
					my $bad2 = &randomElem ('2', \@sups);
					push @{$food_sup_bad2}, @{$bad2};
					my $bad3 = &randomElem ('3', \@sups);
					push @{$food_sup_bad3}, @{$bad3};
				}
			}
		}
		else {
			if ($Act =~/^益生元/) {
				$Act=~/^[^：]+：(.*)/;
				my @sups=split/；/,$1;
				unshift @{$probiotics_sup},@sups;
			}
			if ($Act =~/^益生菌/) {
				$Act=~/^[^：]+：(.*)/;
				my @sups=split/；/,$1;
				unshift @{$prebiotics_sup},@sups;
			}
			if ($Act =~/^食物/) {
				$Act=~/^[^：]+：(.*)/;
				my @sups=split/；/,$1;
				unshift @{$food_sup_all},@sups;
				my $bad1 = &randomElem ('1', \@sups);
				push @{$food_sup_bad1}, @{$bad1};
				my $bad2 = &randomElem ('2', \@sups);
				push @{$food_sup_bad2}, @{$bad2};
				my $bad3 = &randomElem ('3', \@sups);
				push @{$food_sup_bad3}, @{$bad3};
			}
		}
	}
}

# &randomElem($want, $array);
sub randomElem {
	my ($want, $array) = @_ ;
	my (%seen, @ret);
	while ( @ret != $want ) {
	my $num = abs(int(rand($#{$array})));
		if ( ! $seen{$num} ) { 
			++$seen{$num};
			push @ret, $$array[$num];
		}
	}
	return (\@ret);
}

#my ($item)=&get_item($vd,$Data->{$meta},$language,$meta,$age);
sub get_item {
	my ($vd,$Data,$language,$key,$age)=@_;
	my %item=();
	
	##
	#$item{value}=sprintf "%.2f",$Data->{sample}*1000;
	#$item{Val30}=sprintf "%.3f",$Data->{Val30}*100 if (defined $Data->{Val30});
	#$item{Val70}=sprintf "%.3f",$Data->{Val70}*100 if (defined $Data->{Val70});
	$item{value}=(int($Data->{sample}*100000))/100;
	$item{Val30}=$Data->{Val30} if (defined $Data->{Val30});
	$item{Val70}=$Data->{Val70} if (defined $Data->{Val70});
	$item{level}=$Data->{level};

	#$item{per100}=sprintf "%.2f",$Data->{per}*100;
	$item{per100}= int($Data->{per}*100);
	$item{per100} = "1" if ($Data->{per}*100 < 1);
	$item{per100} = "99" if ($item{per100} >= 100);
	if ($key eq "tjsjcjunqunduoyangxing") {
		$item{value}=sprintf "%.2f",$Data->{sample};
		$item{mean}=sprintf "%.2f",$Data->{ref}{mean}*100;

		## caculate dotted line position (10.2*div_value/div_pic_max);
		$item{valuelocation}=10.2*$item{value}/6;
	}
	if ($key eq "tjsjcyouyijun" || $key eq "tjsjcyouhaijun") {

		$item{value} = "0.01" if ($item{value} eq "0.00" && $key eq "tjsjcyouyijun");
		$item{dropcoordinate} = &get_drop_coordinate ($item{per100});
		$item{strawpic} = "straw_"."$item{level}.pdf";
		$item{droppic} = "drop_"."$item{level}.pdf";
		#$item{mean}=sprintf "%.2f",$Data->{ref}{mean}*100;
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

	$item{range}=$Data->{range} if (exists $Data->{range});
	$item{name}=$vd->{$language}{title};
	$item{desc}=$vd->{$language}{summary}{desc};
	if ($vd->{$language}{summary}{descfull} ne "") {
		my @descs = split/\n+/,$vd->{$language}{summary}{descfull};
		my $final_desc = pop @descs;
		if (@descs >= 1) {
			foreach my $desc_line (@descs) {
				$item{descfull} .= $desc_line."\\\\\n";
			}
		}
		$item{descfull} .= $final_desc;
	}

	$item{source}=$vd->{$language}{summary}{source};

	# 转换有益为正常
	if ($vd->{$language}->{suggestion}->{healthy}->{l1}->{effect} eq "有害" && defined $item{Val30}) {
		$item{Val30} = 0.01 if ($item{Val30} == 0);
		if ($item{value} < $item{Val30}) {
			$item{level} = "l2";
		}
		$item{threshold} = "≥ $item{Val30}" ;
		$item{level} = "l3" if ($item{level} eq "l4" or $item{level} eq "l5");
	}
	#if ($vd->{$language}->{suggestion}->{healthy}->{l5}->{effect} eq "有害") {
	#	$item{threshold} = "≤ $item{Val70}";
	#}

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
		if ($vd->{$language}{suggestion}{healthy}{$Data->{level}}{note} ne "") {
			my @notes = split/\n+/,$vd->{$language}{suggestion}{healthy}{$Data->{level}}{note};
			foreach my $note (@notes) {
				if ($note =~ /^（(\d+)\\textasciitilde (\d+)）(.*)/) {
					if ($age >= $1 && $age <= $2) {
						$item{note} .= $3." \\\\\n";
					}
				}
				elsif ($note =~ /^（(\d+)\+）(.*)/) {
					if ($age >= $1) {
						$item{note} .= $2." \\\\\n";
					}
				}
				else {
					$item{note} .= $note." \\\\\n";
				}
			}
		}
	}
	else {
		$item{risk}='-';
		$item{riskdesc}='-';
	}
	
	return (\%item);
}

#$item{dropcoordinate} = &get_drop_coordinate ($item{per100});
sub get_drop_coordinate {
	my ($per) = @_;
	my $coordinate = '-';

	if ($per < 10) {
		$coordinate = $per/9*(1.7-0.4)+0.4;
	}
	elsif ($per < 30) {
		$coordinate = ($per-10)/19*(4.46-1.78)+1.78;
	}
	elsif ($per < 70) {
		$coordinate = ($per-30)/39*(10.06-4.54)+4.54;
	}
	elsif ($per < 90) {
		$coordinate = ($per-70)/19*(12.85-10.12)+10.12;
	}
	else {
		$coordinate = ($per-90)/9*(14.25-12.93)+12.93;
	}

	return($coordinate);
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
	$file="$referdir/qpcr_refer.extract.xls";
	$check+=&checkF($file);
	#
	$file="$referdir/qpcr_refer.extract.stat.xls";
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
	$file="$confdir/qpcr_overview.conf";
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
			elsif ($head[$i] eq "examcode") {
				if ($unit[$i] eq "") {
					$unit[$i]=decode("UTF-8","-");
				}
				else {
					$unit[$i]=decode("UTF-8",$unit[$i]);
				}
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
	my %overviewcomposition=();
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
		&get_list($file,\%listH,\%term2oriid,"distribution");
	}
	## overview.list.order$sex
	$file="$indir/overview.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"overview");
	}
	
	## overviews term composition of distribution oriid and weigth
	$file="$confdir/qpcr_overview.conf";
	if (-f $file) {
		&get_overview_composition($file, \%overviewcomposition, \%compositionweigth);
	}

	return (\%listH,\%term2oriid,\%overviewcomposition,\%compositionweigth);
}

#&get_overview_composition($file,\%overviewcomposition);
sub get_overview_composition {
	my ($file,$overviewcomposition,$compositionweigth)=@_;
	
	##
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($overview_oriid, $distribution_oriid, $weigth) = (split/\t/,$_)[0,1,2];
		$overviewcomposition->{$overview_oriid}{$distribution_oriid} = 1;
		$compositionweigth->{$distribution_oriid} = $weigth;
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

#&calculate_overview_level_per($Data,$overviewcomposition,$level_set);
sub calculate_overview_level_per {
	my ($Data, $overviewcomposition, $compositionweigth, $level_set) = @_;

	foreach my $overview_oriid (keys %{$overviewcomposition}) {
		my $item_score = 0;
		my $item_genus_num = 0;
		my $item_val = 0;
		foreach my $oriid (sort keys %{$overviewcomposition->{$overview_oriid}}) {
			if (exists $Data->{$oriid}) {
				$item_score += $Data->{$oriid}{per}*$compositionweigth->{$oriid};
				$item_val += $Data->{$oriid}{sample};
				$item_genus_num ++;
			}
			else {
				print STDERR "WARN: Unexists of $oriid of overview index $overview_oriid...\n\n";
			}
		}
		$Data->{$overview_oriid}{per} = $item_score/$item_genus_num;
		$Data->{$overview_oriid}{sample} = $item_val;

		$Data->{$overview_oriid}{level} = &get_overview_level($Data->{$overview_oriid}{per}, $level_set);
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
	my ($per,$val30,$val70)=('-',0,0);

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

	## use perl get percentile
	#elsif (exists $meta->{sample} && exists $meta->{ref}) {
	#	my $N=scalar @{$meta->{ref}{values}};
	#	my $n=0;
	#	foreach my $v (@{$meta->{ref}{values}}) {
	#		if ($meta->{sample} < $v) {
	#			last;
	#		}
	#		$n++;
	#	}
	#	$per=$n/$N if ($N != 0);
	#}
	if ($key eq "tjsjcaikemanshijunshu") {
		if ($meta->{sample} == 0) {
			$per = rand(0.25)+0.01;
		}
	}

	$val30 = (int($val30*100000)+1)/100;
	$val70 = (int($val70*100000)+1)/100;

	#my $level=&get_per_level ($per*$weigth, $level_set);
	my $level = "-";
	if ($meta->{sample}*1000 < $val30) {
		$level = "l2";
	}
	if ($meta->{sample}*1000 >= $val30) {
		$level = "l3";
	}
	$level = $level_set if ($level_set ne "N");

	return ($per,$level,$val30,$val70);
}

#&get_overview_level($per,$level)
sub get_overview_level {
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
	$file="$referdir/qpcr_refer.extract.xls";
	&get_ref_values($file,$term2oriid,\%Data);

	## get ref mead && sd
	$file="$referdir/qpcr_refer.extract.stat.xls";
	&get_ref_mean_sd($file,$term2oriid,\%Data);

	return (\%Data);
}

#&get_sample_values($file,$term2oriid,\%Data);
sub get_sample_values {
	my ($file,$term2oriid,$Data)=@_;

	##
	$/="\n";
	#my $benificial_val = 0;
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$value)=(split /\t/,$_)[0,1];
		next if (!exists $term2oriid->{$id});
		$id = $term2oriid->{$id};
		$Data->{$id}{sample}=$value;
		#if ($id ne "tjsjcyouyijun") {
		#	$benificial_val += $value;
		#}
	}
	close IN;

	#$Data->{"tjsjcyouyijun"}{sample} = $benificial_val if (!exists $Data->{"tjsjcyouyijun"});

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
                          default "$Bin/../../reference_data/16s_qpcr_current";
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
