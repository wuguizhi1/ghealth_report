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

my $R = Statistics::R->new();

###
my ($datadir,$db2xml,$level_per,$dirout,$language,$barcode,$historydir,$term2bac,$help);
my ($host,$port,$database,$collection,$updateall);
GetOptions (
	"data=s"	=> \$datadir,
	"db2xml=s"	=> \$db2xml,
	"level_per=s"	=> \$level_per,
	"out=s"		=> \$dirout,
	"language=s"	=> \$language,
	"barcode=s"	=> \$barcode,
	"historydir=s"	=> \$historydir,
	"term2bac=s"	=> \$term2bac,
	"h"		=> \$help,

	"host=s"	=> \$host,
	"port=s"	=> \$port,
	"database=s"	=> \$database,
	"collection=s"	=> \$collection,
) or &USAGE;
&USAGE if((!defined $datadir || !defined $barcode || !$level_per) || defined $help);

#===============AbsolutePath and default
$dirout ||= "./";
$db2xml ||= "$Bin/db2xml/current";
$language ||= "CN";
$historydir ||= "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
$term2bac ||= "$Bin/conf/term2bacteria.list";

$host ||= "10.0.0.204";
$port ||= "27021";
$database ||= "result";
$collection ||= "mmccdjq";

$datadir = &ABSOLUTE_DIR($datadir);
$dirout = &ABSOLUTE_DIR($dirout);
$db2xml = &ABSOLUTE_DIR($db2xml);
$level_per = &ABSOLUTE_DIR($level_per);
$historydir = &ABSOLUTE_DIR($historydir);
$term2bac = &ABSOLUTE_DIR($term2bac);

#==================
my %convertlevel = (
	"L1" => "l1",
	"L2" => "l2",
	"L3" => "l3",
	"L4" => "l4",
	"L5" => "l5",
);

my %final_report;
#===========load sample csv
my ($samInfo) = &load_sample_info($barcode,$datadir,\%final_report);

#===========load sample summary infomation
my ($summary) = &load_summary("$datadir/$barcode.summary.xls");

#===========load $type.list.order$sex
my ($listH,$term2oriid) = &load_list($db2xml,$samInfo);

#===========load term2bacteria
my ($term2bacteria) = &load_term_bacteria($term2bac);

#==========load term level
my ($term2level) = &load_level($level_per);

#=========健康跑分
&total_score_xml($term2level,$listH->{metasummary},"$db2xml/XML/metasummary/",$language,\%final_report);

#==========balance
&balance_xml($summary,$term2level,$term2bacteria,$term2oriid,$listH->{overview},$listH->{distribution},$listH->{pathogens},"$db2xml/XML/overview/","$db2xml/XML/distribution/","$db2xml/XML/pathogens",$language,\%final_report);

#==========Eatperfer
&Eatperfer_xml($term2level,$term2bacteria,$term2oriid,$listH->{metabolism},\%convertlevel,"$db2xml/XML/metabolism/","$db2xml/XML/distribution/",$language,\%final_report);

#==========nutrition
&nutrition_xml($term2level,$term2bacteria,$term2oriid,$listH->{overview},$listH->{metabolism},\%convertlevel,"$db2xml/XML/overview/","$db2xml/XML/metabolism/",$language,\%final_report);

#==========quality
&quality_xml($term2level,$term2bacteria,$term2oriid,$listH->{overview},$listH->{distribution},\%convertlevel,"$db2xml/XML/overview/","$db2xml/XML/distribution/",$language,\%final_report,$samInfo->{age});

#=========emotion
&emotion_xml($term2level,$term2bacteria,$term2oriid,$listH->{overview},$listH->{distribution},$listH->{metabolism},\%convertlevel,"$db2xml/XML/overview/","$db2xml/XML/distribution/","$db2xml/XML/metabolism/",$language,\%final_report);

#=========diseaes
&diseaes_xml($term2level,$term2bacteria,$term2oriid,$listH->{distribution},$listH->{disease},\%convertlevel,"$db2xml/XML/disease/","$db2xml/XML/distribution/",$language,\%final_report);

## insert final report to mongo db
my $mongo = MongoDB::MongoClient->new('host' => $host.':'.$port);
my $db = $mongo->get_database($database);
my $cl = $db->get_collection($collection);
my $barcode_find = $cl->find({_id => $barcode});

$cl->update({_id => $barcode}, {'$set' => {'Report' => \%final_report }},{'upsert' => 1});

#=========print out
&out_json(\%final_report,"$dirout/$barcode.json");

#========sub code
#&out_json(\%final_report,"$dirout/$barcode.json");
sub out_json {
	my ($final_report,$outfile) = @_;

	my $json = encode_json $final_report;
	open OUT,">$outfile";
	print OUT "$json";
	close OUT;
}

sub total_score_xml {
	my ($term2level,$list,$dir,$language,$final_report) = @_;

	$final_report->{total_score}{name} = "健康跑分";
	my $score = sprintf "%.2f",$term2level->{percent}{total_score} * 100;
	$final_report->{total_score}{score} = $score;

	my $xml = XMLin("$dir/$language/$list->{$term2oriid->{'菌群总评分'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($score >= 0 && $score < 10){
		$final_report->{total_score}{summarize} = $xml->{$language}{suggestion}{healthy}{l1}{mechanismfull};
		$final_report->{total_score}{suggest} = $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};	
	}
	if($score >= 10 && $score < 30){
		$final_report->{total_score}{summarize} = $xml->{$language}{suggestion}{healthy}{l2}{mechanismfull};
		$final_report->{total_score}{suggest} = $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	if($score >= 30 && $score < 70){
		$final_report->{total_score}{summarize} = $xml->{$language}{suggestion}{healthy}{l3}{mechanismfull};
		$final_report->{total_score}{suggest} = $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	if($score >= 70 && $score < 90){
		$final_report->{total_score}{summarize} = $xml->{$language}{suggestion}{healthy}{l4}{mechanismfull};
		$final_report->{total_score}{suggest} = $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($score >= 90 && $score <= 100){
		$final_report->{total_score}{summarize} = $xml->{$language}{suggestion}{healthy}{l5}{mechanismfull};
		$final_report->{total_score}{suggest} = $xml->{$language}{suggestion}{healthy}{l5}{actiondesc};
	}
}

sub diseaes_xml {
	my ($term2level,$term2bacteria,$term2oriid,$list1,$list2,$convertlevel,$dir1,$dir2,$language,$final_report) = @_;

	my (%name,%suggest);
	$name{name} = "消化道类疾病风险";
	my ($risk,$level,$xml);

	##炎性肠病
	my %yanxingchangbing;
	&load_part_xml("yanxingchangbing","炎性肠病",$list1,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%yanxingchangbing);
	$risk = $term2level->{level}{yanxingchangbing_risk};
	$xml = XMLin("$dir1/$language/$list2->{$term2oriid->{'炎性肠病'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($risk eq 'L1' || $risk eq 'L2' || $risk eq 'L3'){
		$yanxingchangbing{yanxingchangbing}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	if($risk eq 'L4' || $risk eq 'L5'){
		$yanxingchangbing{yanxingchangbing}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}

	##肠易激综合征
	my %changyiji;
	&load_part_xml("changyiji","肠易激综合征",$list1,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%changyiji);
	$risk = $term2level->{level}{changyiji_risk};
	$xml = XMLin("$dir1/$language/$list2->{$term2oriid->{'肠易激综合征'}}.xml",NoAttr=>1,SuppressEmpty => "");

	if($risk eq 'L1' || $risk eq 'L2' || $risk eq 'L3'){
		$changyiji{changyiji}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	if($risk eq 'L4' || $risk eq 'L5'){
		$changyiji{changyiji}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}

	##结直肠癌
	my %jiezhichangai;
	&load_part_xml("jiezhichangai","结直肠癌",$list1,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%jiezhichangai);
	$risk = $term2level->{level}{jiezhichangai_risk};
	$xml = XMLin("$dir1/$language/$list2->{$term2oriid->{'结直肠癌'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($risk eq 'L1' || $risk eq 'L2' || $risk eq 'L3'){
		$jiezhichangai{jiezhichangai}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	if($risk eq 'L4' || $risk eq 'L5'){
		$jiezhichangai{jiezhichangai}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	
	push @{$final_report->{diseaes}}, \%name, \%yanxingchangbing, \%changyiji, \%jiezhichangai, \%suggest;
}

sub emotion_xml {
	my ($term2level,$term2bacteria,$term2oriid,$list1,$list2,$list3,$convertlevel,$dir1,$dir2,$dir3,$language,$final_report) = @_;

	my (%name,%suggest);
	$name{name} = "情绪管理";
	my ($protect,$risk,$level,$xml);

	my %convert_result = (
		"L1" => "l5",
		"L2" => "l4",
		"L3" => "l3",
		"L4" => "l2",
		"L5" => "l1",
	);

	##焦虑指数
	my %anxietyindex;
	&load_part_xml("anxietyindex","焦虑指数",$list3,$term2bacteria,$term2oriid,$term2level,$dir3,$language,\%anxietyindex);
	my ($level1,$level2);
	foreach my $item(sort keys %{$term2bacteria->{anxietyindex}}){
		foreach my $term(@{$term2bacteria->{anxietyindex}{$item}}){
			my $term =  decode('UTF-8',$term);
			if($term =~/脂多糖/){
				$level1 = $term2level->{level}{$term2oriid->{$term}};
			}
			if($term =~/色氨酸/){
				$level2 = $term2level->{level}{$term2oriid->{$term}};
			}
		}
	}
	$risk = &transfer_emotion_level($level1);
	$protect = &transfer_emotion_level($level2);
	$level = &get_emotion_level($protect,$risk);

	$anxietyindex{anxietyindex}{level} = $level;
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'焦虑指数'}}.xml",NoAttr=>1,SuppressEmpty => "");
	push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{$convertlevel->{$level}}{actiondesc};

	##抗压能力
	my %stresstoler;
	&load_part_xml("stresstoler","抗压能力",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%stresstoler);
	$protect = &transfer_emotion_level($term2level->{level}{stresstoler_protect});
	$risk = &transfer_emotion_level($term2level->{level}{stresstoler_risk});
	$level = &get_emotion_level($protect,$risk);

	$stresstoler{stresstoler}{level} = $level;
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'抗压能力'}}.xml",NoAttr=>1,SuppressEmpty => "");
	push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{$convert_result{$level}}{actiondesc};

	##晴天指数
	my %sunnyindex;
	&load_part_xml("sunnyindex","晴天指数",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%sunnyindex);
	$protect = &transfer_emotion_level($term2level->{level}{sunnyindex_protect});
	$risk = &transfer_emotion_level($term2level->{level}{sunnyindex_risk});
	$level = &get_emotion_level($protect,$risk);

	$sunnyindex{sunnyindex}{level} = $level;
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'晴天指数'}}.xml",NoAttr=>1,SuppressEmpty => "");
	push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{$convert_result{$level}}{actiondesc};

	push @{$final_report->{emotion}}, \%name, \%anxietyindex, \%stresstoler, \%sunnyindex, \%suggest;
	
}

sub transfer_emotion_level {
	my ($in) = @_;
	my $level = "-";

	if($in eq 'L1' || $in eq 'L2'){
		$level = "低";
	}
	if($in eq 'L3'){
		$level = "正常";
	}
	if($in eq 'L4' || $in eq 'L5'){
		$level = "高";
	}

	return ($level);
}
sub get_emotion_level {
	my ($protect,$risk) = @_;
	my $level = "-";

	if($protect eq "正常" && $risk eq "正常" || $protect eq "低" && $risk eq "低" || $protect eq "高" && $risk eq "高"){
		$level = "L3";		
	}
	if($protect eq "正常" && $risk eq "低"){
		$level = "L2";
	}
	if($protect eq "正常" && $risk eq "高"){
		$level = "L4";
	}
	if($protect eq "低" && $risk eq "正常"){
		$level = "L4";
	}
	if($protect eq "低" && $risk eq "高"){
		$level = "L5";
	}
	if($protect eq "高" && $risk eq "低"){
		$level = "L1";
	}
	if($protect eq "高" && $risk eq "正常"){
		$level = "L2";
	}

	return($level);

}

sub quality_xml {
	my ($term2level,$term2bacteria,$term2oriid,$list1,$list2,$convertlevel,$dir1,$dir2,$language,$final_report,$age) = @_;

	my (%name,%suggest);
	my $xml;
	$name{name} = "营养代谢";

	##肠道年龄
	my %gutage;
	$gutage{gutage}{name} = "肠道年龄";
	$gutage{gutage}{age} = $age;
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'肠道年龄'}}.xml",NoAttr=>1,SuppressEmpty => "");

	my $percent = $term2level->{percent}{Bifidobacterium} * 100;
	if($percent >= 0 && $percent < 5){
		$gutage{gutage}{result} = $age + 5;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($percent >= 5 && $percent < 15){
		$gutage{gutage}{result} = $age + 4;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($percent >= 15 && $percent < 25){
		$gutage{gutage}{result} = $age + 3;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($percent >= 25 && $percent < 35){
		$gutage{gutage}{result} = $age + 2;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($percent >= 35 && $percent < 45){
		$gutage{gutage}{result} = $age + 1;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($percent >= 45 && $percent <= 55){
		$gutage{gutage}{result} = $age;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	if($percent > 55 && $percent <= 65){
		$gutage{gutage}{result} = $age - 1;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	if($percent >65 && $percent <= 75){
		$gutage{gutage}{result} = $age - 2;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	if($percent > 75 && $percent <= 85){
		$gutage{gutage}{result} = $age - 3;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	if($percent >85 && $percent <= 95){
		$gutage{gutage}{result} = $age - 4;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	if($percent >95 && $percent <= 100){
		$gutage{gutage}{result} = $age - 5;
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	
	##蹲坑频次
	my %defecatingfreq;
	&load_part_xml("defecatingfreq","蹲坑频次",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%defecatingfreq);

		
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'蹲坑频次'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($term2level->{level}{diarrhea_risk} eq 'L4' || $term2level->{level}{diarrhea_risk} eq 'L5'){
		$defecatingfreq{defecatingfreq}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}

	if($term2level->{level}{diarrhea_risk} eq 'L1' || $term2level->{level}{diarrhea_risk} eq 'L2' || $term2level->{level}{diarrhea_risk} eq 'L3'){
		if($term2level->{level}{constipation_protect} eq 'L3' || $term2level->{level}{constipation_protect} eq 'L4' || $term2level->{level}{constipation_protect} eq 'L5'){
			$defecatingfreq{defecatingfreq}{level} = "L2";
			push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
		}
		if($term2level->{level}{constipation_protect} eq 'L1' || $term2level->{level}{constipation_protect} eq 'L2'){
			$defecatingfreq{defecatingfreq}{level} = "L1";
			push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
		}
	}
	
	##放屁指数
	my %fartpercent;
	&load_part_xml("fartpercent","放屁指数",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%fartpercent);
	
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'放屁指数'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($term2level->{level}{aerogenesis_risk} eq 'L1' || $term2level->{level}{aerogenesis_risk} eq 'L2' || $term2level->{level}{aerogenesis_risk} eq 'L3'){
		$fartpercent{fartpercent}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}
	if($term2level->{level}{aerogenesis_risk} eq 'L4'){
		$fartpercent{fartpercent}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($term2level->{level}{aerogenesis_risk} eq 'L5'){
		$fartpercent{fartpercent}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l5}{actiondesc};
	}

	##肠道月巴指数
	my %obesityindex;
	&load_part_xml("obesityindex","肠道月巴指数",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%obesityindex);

	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'肠道月巴指数'}}.xml",NoAttr=>1,SuppressEmpty => "");

	if($term2level->{level}{obesityindex_risk} eq 'L1' || $term2level->{level}{obesityindex_risk} eq 'L2' || $term2level->{level}{obesityindex_risk} eq 'L3'){
		$obesityindex{obesityindex}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}
	if($term2level->{level}{obesityindex_risk} eq 'L4'){
		$obesityindex{obesityindex}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	if($term2level->{level}{obesityindex_risk} eq 'L5'){
		$obesityindex{obesityindex}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l5}{actiondesc};
	}

	##干吃不胖指数
	my %slimindex;
	&load_part_xml("slimindex","干吃不胖指数",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%slimindex);

	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'干吃不胖指数'}}.xml",NoAttr=>1,SuppressEmpty => "");

	if($term2level->{level}{slimindex_risk} eq 'L3' || $term2level->{level}{slimindex_risk} eq 'L4' || $term2level->{level}{slimindex_risk} eq 'L5'){
		$slimindex{slimindex}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}

	if($term2level->{level}{slimindex_risk} eq 'L2'){
		$slimindex{slimindex}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}
	if($term2level->{level}{slimindex_risk} eq 'L1'){
		$slimindex{slimindex}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}

	push @{$final_report->{quality}}, \%name, \%gutage, \%defecatingfreq, \%fartpercent, \%obesityindex, \%slimindex, \%suggest;
}

sub nutrition_xml {
	my ($term2level,$term2bacteria,$term2oriid,$list1,$list2,$convertlevel,$dir1,$dir2,$language,$final_report) = @_;
	my $xml;

	my (%name,%suggest);
	$name{name} = "营养代谢";
	
	##维生素
	my %vitamin;
	my ($count_vitamin) = &load_part_xml("vitamin","维生素",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%vitamin);
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'维生素'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($count_vitamin->{'维生素'} >= 1 && $count_vitamin->{'维生素'} <= 2){
		$vitamin{vitamin}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}elsif($count_vitamin->{'维生素'} > 2){
		$vitamin{vitamin}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}else{
		$vitamin{vitamin}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	
	
	##氨基酸
	my %aminoacid;
	my ($count_aminoacid) = &load_part_xml("aminoacid","氨基酸",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%aminoacid);
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'氨基酸'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	if($count_aminoacid->{'氨基酸'} >= 1 && $count_aminoacid->{'氨基酸'} <= 2){
		$aminoacid{aminoacid}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}elsif($count_aminoacid->{'氨基酸'} >= 4){
		$aminoacid{aminoacid}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}else{
		$aminoacid{aminoacid}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	
	push @{$final_report->{nutrition}}, \%name, \%vitamin, \%aminoacid, \%suggest;

	return;

}

sub Eatperfer_xml {
	my ($term2level,$term2bacteria,$term2oriid,$list,$convertlevel,$dir1,$dir2,$language,$final_report) = @_;

	my (%name,%suggest);
	$name{name} = "饮食偏好";

	##肉食指数
	my %meatpercent;
	$meatpercent{meatpercent}{name} = "肉食指数";
	my %term_level;
	foreach my $item(sort keys %{$term2bacteria->{meatpercent}}){
		my (@name,@level,$level_info);
		foreach my $term(@{$term2bacteria->{meatpercent}{$item}}){
			my $term =  decode('UTF-8',$term);
			my $level = $term2level->{level}{$term2oriid->{$term}};
			if($level eq "L1" || $level eq "L5"){
				$level_info = "L3";
			}
			if($level eq "L2" || $level eq "L4"){
				$level_info = "L2";
			}
			if($level eq "L3"){
				$level_info = "L1";
			}
			#$term_level{$term} = $level_info;
			my $xml = XMLin("$dir1/$language/$list->{$term2oriid->{$term}}.xml",NoAttr=>1,SuppressEmpty => "");
			
			### suggestion
			push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{$convertlevel->{$level}}{actiondesc};
			
			#my $term =  decode('UTF-8',$term);
			push @name, $term;
			push @level, $level_info;
		}
		
		my ($item_info) = &combine_term(\@name,\@level);
		foreach my $i(sort keys %$item_info){
			push @{$meatpercent{meatpercent}{$item}}, $item_info->{$i};
		}
		 
	}

	if(($term2level->{level}{$term2oriid->{'碳水化合物降解产物'}} eq 'L1' || $term2level->{level}{$term2oriid->{'碳水化合物降解产物'}} eq 'L2' || $term2level->{level}{$term2oriid->{'碳水化合物降解产物'}} eq 'L3') && (($term2level->{level}{$term2oriid->{'蛋白降解产物'}} eq 'L4' || $term2level->{level}{$term2oriid->{'蛋白降解产物'}} eq 'L5') || ($term2level->{level}{$term2oriid->{'脂肪降解产物'}} eq 'L4' || $term2level->{level}{$term2oriid->{'脂肪降解产物'}} eq 'L5'))){
		$meatpercent{meatpercent}{level} = "L3";
	}
	elsif(($term2level->{level}{$term2oriid->{'碳水化合物降解产物'}} eq 'L4' || $term2level->{level}{$term2oriid->{'碳水化合物降解产物'}} eq 'L5') && (($term2level->{level}{$term2oriid->{'蛋白降解产物'}} eq 'L1' || $term2level->{level}{$term2oriid->{'蛋白降解产物'}} eq 'L2' || $term2level->{level}{$term2oriid->{'蛋白降解产物'}} eq 'L3') || ($term2level->{level}{$term2oriid->{'脂肪降解产物'}} eq 'L1' || $term2level->{level}{$term2oriid->{'脂肪降解产物'}} eq 'L2' || $term2level->{level}{$term2oriid->{'脂肪降解产物'}} eq 'L3'))){
		$meatpercent{meatpercent}{level} = "L1";
	}
	else{
		$meatpercent{meatpercent}{level} = "L2";
	}

	push @{$final_report->{Eatperfer}}, \%name, \%meatpercent, \%suggest;

	return;	
}

sub balance_xml {
	my ($summary,$term2level,$term2bacteria,$term2oriid,$list1,$list2,$list3,$dir1,$dir2,$dir3,$language,$final_report) = @_;
	my $xml;

	my (%name,%suggest);
	$name{name} = "菌群平衡";

	##菌群多样性
	my %diversity;
	$diversity{diversity}{name} = "菌群多样性";
	$diversity{diversity}{otunumber} = $summary->{Classified_OTU};
	$diversity{diversity}{score} = sprintf "%.2f",$term2level->{percent}{diversity} * 100;
	$diversity{diversity}{summary} = "本次检测您的肠道菌群有" . $summary->{Classified_OTU} . "种，";
	$diversity{diversity}{summary} .= "多样性检测结果为" . $diversity{diversity}{score} . "分，";
	my $species = $summary->{Classified_OTU};
	my $score = $diversity{diversity}{score};

	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'菌群多样性'}}.xml",NoAttr=>1,SuppressEmpty => "");
	if($diversity{diversity}{score} < 5 ){
		$diversity{diversity}{level} = "L4";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
		$diversity{diversity}{summary} .= "多样性脆弱，需要关注";
	}
	if($diversity{diversity}{score} >= 5 && $diversity{diversity}{score} < 15){
		$diversity{diversity}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
		$diversity{diversity}{summary} .= "多样性单一，需要关注";
	}
	if($diversity{diversity}{score} >= 15 && $diversity{diversity}{score} < 30){
		$diversity{diversity}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
		$diversity{diversity}{summary} .= "多样性较低，需要关注";
	}
	if($diversity{diversity}{score} >= 30){
		$diversity{diversity}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
		$diversity{diversity}{summary} .= "多样性稳定。";
	}

	##有益菌
	my %benifical;
	$benifical{benifical}{refer} = "11.62";
	$benifical{benifical}{percent} = sprintf "%.2f",$term2level->{percent}{$term2oriid->{'有益菌'}} * 100;
	$benifical{benifical}{value} = sprintf "%.2f",$term2level->{value}{$term2oriid->{'有益菌'}} * 100;

	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'有益菌'}}.xml",NoAttr=>1,SuppressEmpty => "");
	if($term2level->{level}{$term2oriid->{'有益菌'}}  eq 'L1'){
		$benifical{benifical}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有益菌'}}  eq 'L2'){
		$benifical{benifical}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有益菌'}}  eq 'L3'){
		$benifical{benifical}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有益菌'}}  eq 'L4'){
		$benifical{benifical}{level} = "L1";		
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有益菌'}}  eq 'L5'){
		$benifical{benifical}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l5}{actiondesc};
	}

	my ($count_benifical) = &load_part_xml("benifical","有益菌",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%benifical);
	
	##有害菌
	my %harmful;
	$harmful{harmful}{refer} = "5.57";
	$harmful{harmful}{percent} = sprintf "%.2f",$term2level->{percent}{$term2oriid->{'有害菌'}} * 100;
	$harmful{harmful}{value} = sprintf "%.2f",$term2level->{value}{$term2oriid->{'有害菌'}} * 100;

	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'有害菌'}}.xml",NoAttr=>1,SuppressEmpty => "");
	if($term2level->{level}{$term2oriid->{'有害菌'}}  eq 'L1'){
		$harmful{harmful}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有害菌'}} eq 'L2'){
		$harmful{harmful}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有害菌'}}  eq 'L3'){
		$harmful{harmful}{level} = "L1";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}

	if($term2level->{level}{$term2oriid->{'有害菌'}}  eq 'L4'){
		$harmful{harmful}{level} = "L2";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l4}{actiondesc};
	}
	
	if($term2level->{level}{$term2oriid->{'有害菌'}} eq 'L5'){
		$harmful{harmful}{level} = "L3";
		push @{$suggest{suggest}}, $xml->{$language}{suggestion}{healthy}{l5}{actiondesc};
	}

	my ($count_harmful) = &load_part_xml("harmful","有害菌",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,\%harmful);

	###致病菌
	my %pathogen;
	&load_pathogen("harmful","items2",$list3,$dir3,$term2level,$term2oriid,\%pathogen,$language,\%harmful);

	$harmful{harmful}{pathogen} = $pathogen{number};
	$harmful{harmful}{bad_pathogen} = $pathogen{bad};

	##菌群紊乱指数
	my %BEpercent;
	$BEpercent{BEpercent}{name} = "菌群紊乱指数";
	$BEpercent{BEpercent}{score} = sprintf "%.2f",$term2level->{percent}{BEpercent};
	$xml = XMLin("$dir1/$language/$list1->{$term2oriid->{'菌群紊乱指数'}}.xml",NoAttr=>1,SuppressEmpty => "");
	
	my $tmp;
	if($term2level->{percent}{BEpercent} < 0.01){
		$BEpercent{BEpercent}{lavel} = "L3";
		$tmp = $xml->{$language}{suggestion}{healthy}{l1}{actiondesc};
	}elsif($term2level->{percent}{BEpercent} >= 0.01 && $term2level->{percent}{BEpercent} < 1){
		$BEpercent{BEpercent}{lavel} = "L2";
		$tmp =  $xml->{$language}{suggestion}{healthy}{l2}{actiondesc};
	}elsif($term2level->{percent}{BEpercent} >= 1){
		$BEpercent{BEpercent}{lavel} = "L1";
		$tmp =  $xml->{$language}{suggestion}{healthy}{l3}{actiondesc};
	}
	$tmp =~s/\{\/\}/\//g;
	push @{$suggest{suggest}},$tmp;

	push @{$final_report->{balance}}, \%name, \%diversity, \%benifical, \%harmful, \%BEpercent, \%suggest;
	
}

sub load_pathogen {
	my ($part,$item,$list,$dir,$term2level,$term2oriid,$pathogen,$language,$final_report) = @_;
	my (@name,@level);	

	foreach my $meta(sort keys %$list){
		next if(!exists $term2level->{level}{$meta});
		next if($term2level->{value}{$meta} == 0);
		$pathogen->{number}++;#print "$term2oriid->{$meta}\n";
		#my $xml = XMLin("$dir/$language/$list->{$term2oriid->{$meta}}.xml",NoAttr=>1,SuppressEmpty => "");print Dumper $xml,"\n";
		if($term2level->{level}{$meta} eq 'L3' || $term2level->{level}{$meta} eq 'L2' || $term2level->{level}{$meta} eq 'L1'){
			push @name,$term2oriid->{$meta};
			push @level,"L1";
		}
		if($term2level->{level}{$meta} eq 'L4'){
			push @name,$term2oriid->{$meta}; 
			push @level,"L2";
			$pathogen->{bad}++;
		}
		if($term2level->{level}{$meta} eq 'L5'){
			push @name,$term2oriid->{$meta};
			push @level,"L3";
			$pathogen->{bad}++;
		}	
	}
	
	my ($item_info) = &combine_term(\@name,\@level);

	foreach my $i(sort keys %$item_info){
		push @{$final_report->{$part}{$item}}, $item_info->{$i};
	}

}

sub load_part_xml {
        my ($part,$name,$list,$term2bacteria,$term2oriid,$term2level,$dir,$language,$final_report) = @_;
	my %count;

	##有害菌等级转换
	my %convert1 = (
		"L1" => "L1",
		"L2" => "L1",
		"L3" => "L1",
		"L4" => "L2",
		"L5" => "L3",
	);

	##有益菌等级转换
	my %convert2 = (
		"L1" => "L3",
		"L2" => "L2",
		"L3" => "L1",
		"L4" => "L1",
		"L5" => "L1",
	);

        $final_report->{$part}{name} = $name;
        foreach my $item(sort keys %{$term2bacteria->{$part}}){
                my (@name,@level);
                foreach my $term(@{$term2bacteria->{$part}{$item}}){
			my $term =  decode('UTF-8',$term);
                        my $level = $term2level->{level}{$term2oriid->{$term}};
			
                        my $xml = XMLin("$dir/$language/$list->{$term2oriid->{$term}}.xml",NoAttr=>1,SuppressEmpty => "");
	
			my $effect = $xml->{$language}{suggestion}{healthy}{l1}{effect};
			
			if($effect =~/有害/){
				$level = $convert2{$level};		
			}else{
				$level = $convert1{$level};
			}

#                        my $term =  decode('UTF-8',$term);
                        push @name, $term;
                        push @level, $level;
                }

		$count{$name} = 0;
		foreach my $m(@level){
			if($m eq 'L2' || $m eq 'L3'){
				$count{$name}++;
			}
		}

                my ($item_info) = &combine_term(\@name,\@level);
                foreach my $i(sort keys %$item_info){
                        push @{$final_report->{$part}{$item}}, $item_info->{$i};
                }
        }

        return(\%count);
}

#my (@item_info) = &combine_term(\@name,\@level,$item);
sub combine_term {
	my ($name,$level) = @_;
	my @item_info = ();
	my %info = ();

	for(my $i=0;$i<@$name;$i++){
		$info{$i}{name} = $name->[$i];
		$info{$i}{level} = $level->[$i];
	}

	return (\%info);
}

#my ($term2level) = &load_level($level_per);
sub load_level {
	my ($level_per) = @_;
	my %term2level = ();

	open IN,"$level_per";
	while(<IN>){
		chomp;
		next if(/^$/ || /^#/);
		my @unit = split/\t/,$_;
		next if($unit[2] =~/-/);
		$term2level{value}{$unit[0]} = $unit[1];
		$term2level{level}{$unit[0]} = $unit[2];
		$term2level{percent}{$unit[0]} = $unit[3];
		if($unit[1] ne 'NA' && $unit[1] == "0" && ($unit[2] eq 'L4' || $unit[2] eq 'L5')){
			$term2level{level}{$unit[0]} = "L3";
		}
	}
	close IN;

	return (\%term2level);
}

#my ($listH,$term2oriid)=&load_list($indir,$samInfo);
sub load_list {
	my ($db2xml,$samInfo) = @_;
	my %listH=();
	my %term2oriid=();

	## pcoress sex
	$samInfo->{sex}='A' if (! defined $samInfo->{sex} || ($samInfo->{sex} ne 'F' && $samInfo->{sex} ne 'M'));

	my $file;
	## disease.list.order$sex
	$file="$db2xml/disease.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"disease");
	}
	## distribution.list.order$sex
	$file="$db2xml/distribution.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"distribution");
	}
	## metabolism.list.order$sex
	$file="$db2xml/metabolism.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"metabolism");
	}
	## metasummary.list.order$sex
	$file="$db2xml/metasummary.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"metasummary");
	}
	## overview.list.order$sex
	$file="$db2xml/overview.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"overview");
	}
	## pathogens.list.order$sex
	$file="$db2xml/pathogens.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"pathogens");
	}

	return (\%listH,\%term2oriid);
}

sub get_list {
	my ($file,$listH,$term2oriid,$type) = @_;

	$/="\#\#\#";
	open IN,"$file" or die $!;
	while(<IN>){
		chomp;
		next if (/^$/);
		my ($section1,$units)=split /\n/,$_,2;
		$section1=(split /\s+/,$section1)[0];
		#$section1=decode('UTF-8',$section1);
		my @unit=split /\#\#/,$units;
		for (my $j=1;$j<@unit;$j++){
			my @lines=split /\n/,$unit[$j];
			my $section2;
			if($lines[0] eq ""){
				$section2="NA";
			}else{
				$section2=(split /\s+/,$lines[0])[0];
				#$section2=decode('UTF-8',$section2);
			}
			#$listH->{$type}{$section1}{$j}{name}=$section2;  #print "$section1\t$section2\n";
			for(my $i=1;$i<@lines;$i++){
				my @aa=split /\t/,$lines[$i];
				my ($oriid,$en,$cn);
				$oriid=$aa[0];
				$cn=$aa[-1];
				$cn=decode('UTF-8',$cn);
				if(scalar @aa == 3){
					$en=$aa[1];
				}
				#$listH->{$type}{$oriid}{order}=$i;
				#$listH->{$type}{$oriid}{cn}=$cn;
				if(defined $en){
					#$listH->{$type}{$oriid}{en}=$en;
					#$term2oriid->{$en}=$oriid;
					$term2oriid->{$cn} = $en;
					$term2oriid->{$en} = $cn;
					$listH->{$type}{$en} = $oriid;
				}
			}
		}
	}
	$/="\n";
	close IN;

	return;
	
}

#my ($term2bacteria) = &load_term_bacteria($term2bac);
sub load_term_bacteria {
	my ($term2bac) = @_;
	my %term2bacteria = ();

	open IN,"$term2bac";
	while(<IN>){
		chomp;
		next if(/^#/ || /^$/);
		my @unit = split/\t/,$_;
		push @{$term2bacteria{$unit[0]}{$unit[1]}}, $unit[2];
	}
	close IN;

	return (\%term2bacteria);
}

#my ($samInfo) = &load_sample_info($barcode,$datadir,$final_report);
sub load_sample_info {
	my ($barcode,$datadir,$final_report) = @_;
	my %samInfo = ();

	my $file="$datadir/$barcode".".csv";
	open IN,$file or die $!;
	my $head=<IN>;chomp($head);
	my @head=split /\s+/,$head;
	while(<IN>){
		chomp;
		next if (/^$/);
		my @unit=split /\t/,$_;
		for (my $i=0;$i<@head;$i++){
			if($head[$i] eq 'barcode' || $head[$i] eq 'sex' || $head[$i] eq 'height' || $head[$i] eq 'weight' || $head[$i] eq 'batch' || $head[$i] eq 'receivedate'){
				$samInfo{$head[$i]} = $unit[$i];
				if($head[$i] eq "barcode"){
					$final_report->{$head[$i]} = $unit[$i];
					#$final_report->{"_id"} = $unit[$i];
				}
			}elsif($head[$i] eq 'name'){
				$unit[$i] = decode("UTF-8",$unit[$i]);
				$samInfo{$head[$i]} = $unit[$i];
				$final_report->{$head[$i]} = $unit[$i];
			}elsif($head[$i] eq 'age'){
				$samInfo{$head[$i]} = $unit[$i];
				#$final_report->{$head[$i]} = $unit[$i];
				$unit[$i] = decode("UTF-8",$unit[$i]);
				if($unit[$i] =~/未知/){
					$samInfo{$head[$i]} = 30;
				#	$final_report->{$head[$i]} = 30;
				}
			}elsif($head[$i] eq 'source'){
				if($unit[$i] eq ""){
					$unit[$i] = decode("UTF-8","-");
				}else{
					$unit[$i] = decode("UTF-8",$unit[$i]);
				}
				$samInfo{$head[$i]} = $unit[$i];
				$final_report->{$head[$i]} = $unit[$i];
			}elsif($head[$i] eq 'version'){
				$unit[$i] = decode("UTF-8",$unit[$i]);
				$samInfo{$head[$i]} = $unit[$i];
				#$final_report->{$head[$i]} = $unit[$i];
				$final_report->{$head[$i]} = "肠道全菌谱检测";
			}
		}
	}
	close IN;

	return (\%samInfo);
}

sub load_summary {
	my ($file) = @_;
	my %summary;

	$/="\n";
	open IN,$file or die $!;
	while (<IN>){
		chomp;
		next if (/^\#/ || /^$/);
		my ($id,$value)=(split /\s+/,$_)[0,1];
		$summary{$id}=$value;
	}
	close IN;

	return (\%summary);
}

sub ABSOLUTE_DIR {
	my $cur_dir = `pwd`;chomp($cur_dir);
	my ($in) = @_;
	my $return = "";

	if(-f $in){
		my $dir = dirname($in);
		my $file = basename($in);
		chdir $dir;$dir = `pwd`;chomp($dir);
		$return = "$dir/$file";
	}
	elsif(-d $in){
		chdir $in;$return = `pwd`;chomp($return);
	}
	else{
		warn "Warning just for file and dir, $in is no exist\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub USAGE {
	my $usage=<<"USAGE";

	Function: result to json and level file
	Usage:
		Options:
		-data		<dir>	sample's report result dir for input	forced

		-db2xml	<dir>	product set database list and xml dir	optional
					default "$Bin/db2xml/current"

		-level_per	<file>	sample's level result			forced

		-out		<dir>	output dir				optional
					default "./"	

		-language	<str>	result language used for report		optional
					default "CN"

		-barcode	<str>	sample's code				forced

		-historydir	<dir>	history test xml result for report dir	optional
					default "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml"

		-h		Help

		-host		<num>	host ID
					default "10.0.0.204"
		-port		<num>	port number
					default "27021"
		-database	<str>	database
					default "result"
		-collection	<str>	default "mmccdjq"

USAGE
	print $usage;
	exit;
}
