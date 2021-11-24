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

use Statistics::R;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my $R = Statistics::R->new();

###
my ($datadir,$db2xml,$level_per,$dirout,$language,$barcode,$historydir,$term2bac,$help);
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
) or &USAGE;
&USAGE if((!defined $datadir || !defined $barcode || !$level_per) || defined $help);

#===============AbsolutePath and default
$dirout ||= "./";
$db2xml ||= "$Bin/db2xml/current";
$language ||= "CN";
$historydir ||= "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
$term2bac ||= "$Bin/conf/term2bacteria.list";

$datadir = &ABSOLUTE_DIR($datadir);
$dirout = &ABSOLUTE_DIR($dirout);
$db2xml = &ABSOLUTE_DIR($db2xml);
$level_per = &ABSOLUTE_DIR($level_per);
$historydir = &ABSOLUTE_DIR($historydir);
$term2bac = &ABSOLUTE_DIR($term2bac);

#==================
#my %convertlevel = (
#	"L1" => "l1",
#	"L2" => "l2",
#	"L3" => "l3",
#	"L4" => "l4",
#	"L5" => "l5",
#);

my %diseaseOrder = (
	"question" => {
		"肠易激综合征" => 0,
		"2型糖尿病" => 3,
		"炎症性肠病" => 1,
		"心脑血管疾病" => 5,
		"结直肠癌" => 4,
		"便秘" => 2,
		"肥胖" => 6,
	},
	"predict" => {
		"肠易激综合征" => 0,
		"2型糖尿病" => 3,
		"炎症性肠病" => 1,
		"心脑血管疾病" => 5,
		"结直肠癌" => 4,
		"便秘" => 2,
		"肥胖" => 6,
	},
	"names" => {
		"肠易激综合征" => "disease0",
		"2型糖尿病" => "disease3",
		"炎症性肠病" => "disease1",
		"心脑血管疾病" => "disease5",
		"结直肠癌" => "disease4",
		"便秘" => "disease2",
		"肥胖" => "disease6",
	},
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
#&total_score_xml($term2level,$listH->{metasummary},"$db2xml/XML/metasummary/",$language,\%final_report);

#==========balance
my %pathogen;
&balance_xml($summary,$term2level,$term2bacteria,$term2oriid,$listH->{overview},$listH->{distribution},$listH->{pathogens},"$db2xml/XML/overview/","$db2xml/XML/distribution/","$db2xml/XML/pathogens",\%pathogen,$language,\%final_report);

#==========quality
&quality_xml($term2level,$term2bacteria,$term2oriid,$listH->{overview},$listH->{distribution},"$db2xml/XML/overview/","$db2xml/XML/distribution/",$language,\%final_report,$samInfo->{age});

#==========advise
&advise_xml($term2level,"$db2xml/XML/advise",$language,$samInfo,$listH->{advise},\%{$diseaseOrder{question}},\%{$diseaseOrder{predict}},\%{$diseaseOrder{names}},\%pathogen,\%final_report);

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

sub advise_xml {
	my ($term2level,$dir,$language,$samInfo,$list,$listQ,$listP,$listN,$pathogen,$final_report) = @_;
	my ($foodKey,$regulateKey,$sportKey,$age)=('-','-','-','-');

	## get food key
	($foodKey)=&get_program_food($samInfo,$term2level,$pathogen,$listQ,$listN);
	## get regulate key
	($regulateKey)=&get_program_regulate($samInfo,$term2level,$pathogen,$listQ,$listN);
	## get sport key
	($sportKey)=&get_program_food($samInfo,$term2level,$pathogen,$listQ,$listN);
	## if no food disease, get age
	($age)=&get_program_age($samInfo->{age});

	## process food and sport
	&process_food_sport($dir,$foodKey,$sportKey,$age,$language,$final_report);
	## process regulate
	&process_regulate($dir,$regulateKey,$age,$language,$final_report);	
}

sub process_regulate {
	my ($dir,$regulateKey,$age,$language,$final_report) = @_;
	
	##补充益生元和益生菌
	my $xml = XMLin("$dir/$language/tjsjctiaojiefangan.xml",NoAttr=>1,SuppressEmpty => "");
	if($regulateKey=~/condition/){
		$final_report->{suggest}{regulate} = $xml->{$language}->{microsuggestion}->{$regulateKey}->{$age}->{action};
	}
	elsif($regulateKey=~/disease/){
		$final_report->{suggest}{regulate} = $xml->{$language}->{dieasesuggestion}->{$regulateKey}->{$age}->{action};
	}
	elsif($regulateKey=~/age/){
		$final_report->{suggest}{regulate} = $xml->{$language}->{agesuggestion}->{desc}->{$regulateKey}->{action};
	}
}

sub process_food_sport {
	my ($dir,$foodKey,$sportKey,$age,$language,$final_report) = @_;
	my $xml;

	##food
	$xml=XMLin("$dir/$language/tjsjcshanshifangan.xml",NoAttr=>1,SuppressEmpty => "");
	$final_report->{suggest}{food} = $xml->{$language}{dieasesuggestion}{$foodKey}{$age}{action} if ($foodKey=~/disease/);
	$final_report->{suggest}{food} = $xml->{$language}{microsuggestion}{$foodKey}{$age}{action} if ($foodKey=~/condition/);
	$final_report->{suggest}{food} = $xml->{$language}{agesuggestion}{desc}{$foodKey}{action} if ($foodKey=~/age/);

	##sport
	$xml=XMLin("$dir/$language/tjsjcyundongfangan.xml",NoAttr=>1,SuppressEmpty => "");
	$final_report->{suggest}{sport} = $xml->{$language}{dieasesuggestion}{$sportKey}{$age}{action} if ($sportKey=~/disease/);
	$final_report->{suggest}{sport} = $xml->{$language}{microsuggestion}{$sportKey}{$age}{action} if ($sportKey=~/condition/);
	$final_report->{suggest}{sport} = $xml->{$language}{agesuggestion}{desc}{$sportKey}{action} if ($sportKey=~/age/);

	return;
}

sub get_program_regulate {
	my ($samInfo,$term2level,$pathogen,$listQ,$listN) = @_;
	my ($disease,$condition,$age)=('-','-','-');

	### harmful items (多样性、有益菌、有害菌)
	if($term2level->{level}{diversity} eq 'L1'){
		$condition="condition1";
	}
	elsif($term2level->{level}{'Benificial bacteria'} eq 'L1'){
		$condition="condition2";
	}
	elsif($term2level->{level}{'Harmful bacteria'} eq 'L5'){
		$condition="condition3";
	}
	elsif($pathogen->{bad} >= 2){
		$condition="condition4";
	}
	elsif($term2level->{level}{diversity} eq 'L2'){
		$condition="condition5";
	}
	elsif($term2level->{level}{'Benificial bacteria'} eq 'L2'){
		$condition="condition6";
	}
	elsif($term2level->{level}{'Harmful bacteria'} eq 'L4'){
		$condition="condition7";
	}
	return ($condition) if ($condition ne '-');

	#### disease
	## sample information
	if ($condition eq '-'){
		foreach my $k (sort {$listQ->{$a} <=> $listQ->{$b}} keys %{$listQ}){
			if (exists $samInfo->{disease}{$k} && $samInfo->{disease}{$k}{risk} eq 'Y'){
				$disease="$listN->{$k}";
				return ($disease);
				last;
			}
		}
	}

	### harmful items (多样性、有益菌、有害菌)
	if ($disease eq '-' && $condition eq '-'){
		$age=&get_program_age($samInfo->{age});
		return ($age);
	}
}

sub get_program_food {
	my ($samInfo,$term2level,$pathogen,$listQ,$listN) = @_;
	my ($disease,$condition,$age)=('-','-','-');

	#### disease
	## sample information
	foreach my $k (sort {$listQ->{$a} <=> $listQ->{$b}} keys %{$listQ}){
		if (exists $samInfo->{disease}{$k} && $samInfo->{disease}{$k}{risk} eq 'Y'){
			$disease="$listN->{$k}";
			return ($disease);
			last;
		}
	}

	### harmful items (多样性、有益菌、有害菌)
	if ($disease eq '-'){
		if($term2level->{level}{diversity} eq 'L1'){
			$condition="condition1";
		}
		elsif($term2level->{level}{'Benificial bacteria'} eq 'L1'){
			$condition="condition2";
		}
		elsif($term2level->{level}{'Harmful bacteria'} eq 'L5'){
			$condition="condition3";
		}
		elsif($pathogen->{bad} >= 2){
			$condition="condition4";
		}
		elsif($term2level->{level}{diversity} eq 'L2'){
			$condition="condition5";
		}
		elsif($term2level->{level}{'Benificial bacteria'} eq 'L2'){
			$condition="condition6";
		}
		elsif($term2level->{level}{'Harmful bacteria'} eq 'L4'){
			$condition="condition7";
		}
		return ($condition) if ($condition ne '-');
	}

	### harmful items (多样性、有益菌、有害菌)
	if ($disease eq '-' && $condition eq '-'){
		$age=&get_program_age($samInfo->{age});
		return ($age);
	}
}

sub get_program_age {
	my ($num)=@_;
	my $age="age5";

	if($num < 1){
		$age="age0";	
	}
	elsif($num < 3){
		$age="age1";
	}
	elsif($num < 7){
		$age="age2";
	}
	elsif($num < 13){
		$age="age3";
	}
	elsif($num < 19){
		$age="age4";
	}
	elsif($num < 45){
		$age="age5";
	}
	elsif($num < 65){
		$age="age6";
	}
	else{
		$age="age7";
	}
	
	return ($age);
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

sub quality_xml {
	my ($term2level,$term2bacteria,$term2oriid,$list1,$list2,$dir1,$dir2,$language,$final_report,$age) = @_;

	##肠道月巴指数
	&load_part_xml("quality","obesityindex","肠道月巴指数",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,$final_report);

	if($term2level->{level}{obesityindex_risk} eq 'L1' || $term2level->{level}{obesityindex_risk} eq 'L2' || $term2level->{level}{obesityindex_risk} eq 'L3'){
		$final_report->{obesityindex}{level} = "L1";
	}
	if($term2level->{level}{obesityindex_risk} eq 'L4'){
		$final_report->{obesityindex}{level} = "L2";
	}
	if($term2level->{level}{obesityindex_risk} eq 'L5'){
		$final_report->{obesityindex}{level} = "L3";
	}

	##干吃不胖指数
	&load_part_xml("quality","slimindex","干吃不胖指数",$list2,$term2bacteria,$term2oriid,$term2level,$dir2,$language,$final_report);

	if($term2level->{level}{slimindex_risk} eq 'L3' || $term2level->{level}{slimindex_risk} eq 'L4' || $term2level->{level}{slimindex_risk} eq 'L5'){
		$final_report->{slimindex}{level} = "L1";
	}

	if($term2level->{level}{slimindex_risk} eq 'L2'){
		$final_report->{slimindex}{level} = "L2";
	}
	if($term2level->{level}{slimindex_risk} eq 'L1'){
		$final_report->{slimindex}{level} = "L3";
	}
}

sub balance_xml {
	my ($summary,$term2level,$term2bacteria,$term2oriid,$list1,$list2,$list3,$dir1,$dir2,$dir3,$pathogen,$language,$final_report) = @_;

	##菌群多样性
	$final_report->{diversity} = sprintf "%.2f",$term2level->{percent}{diversity} * 100;
	
	##有益菌
	$final_report->{benifical}{percent} = sprintf "%.2f",$term2level->{percent}{$term2oriid->{'有益菌'}} * 100;
	$final_report->{benifical}{value} = sprintf "%.2f",$term2level->{value}{$term2oriid->{'有益菌'}} * 100;

	my ($count_benifical) = &load_part_xml("balance","benifical","有益菌",$list2,$term2bacteria,$term2oriid,$term2level,"$dir2",$language,$final_report);

	##有害菌
	$final_report->{harmful}{percent} = sprintf "%.2f",$term2level->{percent}{$term2oriid->{'有害菌'}} * 100;
	$final_report->{harmful}{value} = sprintf "%.2f",$term2level->{value}{$term2oriid->{'有害菌'}} * 100;

	my ($count_harmful) = &load_part_xml("balance","harmful","有害菌",$list2,$term2bacteria,$term2oriid,$term2level,"$dir2",$language,$final_report);

	##致病菌
	&load_pathogen($list3,$dir3,$term2level,$pathogen,$language);

}

sub load_pathogen {
	my ($list,$dir,$term2level,$pathogen,$language) = @_;

	foreach my $meta(sort keys %$list){
		next if(!exists $term2level->{level}{$meta});
		next if($term2level->{value}{$meta} == 0);
		$pathogen->{number}++;
		if($term2level->{level}{$meta} eq 'L4' || $term2level->{level}{$meta} eq 'L5'){
			$pathogen->{bad}++;
		}
	}
}

sub load_part_xml {
        my ($contonts,$part,$name,$list,$term2bacteria,$term2oriid,$term2level,$dir,$language,$final_report) = @_;
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

	##系统对接等级转换
	my %convert2level = (
		"L1" => "l1",
		"L2" => "l2",
		"L3" => "l3",
		"L4" => "l4",
		"L5" => "l5",
	);

        #$final_report->{$contonts}{$part}{name} = $name;
        foreach my $item(sort keys %{$term2bacteria->{$part}}){
                my (@name,@level,@risk);
                foreach my $term(@{$term2bacteria->{$part}{$item}}){
			my $term =  decode('UTF-8',$term);
                        my $level = $term2level->{level}{$term2oriid->{$term}};
			
                        my $xml = XMLin("$dir/$language/$list->{$term2oriid->{$term}}.xml",NoAttr=>1,SuppressEmpty => "");
			my $risk = $xml->{$language}{suggestion}{healthy}{$convert2level{$level}}{mechanism};
			if($risk =~/^$/){
				$risk = "-";
			}
	
			my $effect = $xml->{$language}{suggestion}{healthy}{l1}{effect};
			
			if($effect =~/有害/){
				$level = $convert2{$level};		
			}else{
				$level = $convert1{$level};
			}
			
#                        my $term =  decode('UTF-8',$term);
                        push @name, $term;
                        push @level, $level;
			push @risk, $risk;
                }

		$count{$name} = 0;
		foreach my $m(@level){
			if($m eq 'L2' || $m eq 'L3'){
				$count{$name}++;
			}
		}

                my ($item_info) = &combine_term(\@name,\@level,\@risk);
                foreach my $i(sort keys %$item_info){
                        #push @{$final_report->{$contonts}{$part}{$item}}, $item_info->{$i};
			push @{$final_report->{$part}{$item}}, $item_info->{$i};
                }
        }

        return(\%count);
}

#my (@item_info) = &combine_term(\@name,\@level,$item);
sub combine_term {
	my ($name,$level,$risk) = @_;
	my @item_info = ();
	my %info = ();

	for(my $i=0;$i<@$name;$i++){
		$info{$i}{name} = $name->[$i];
		$info{$i}{value} = $level->[$i];
		$info{$i}{risk} = $risk->[$i];
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
	#$file="$db2xml/disease.list.order".$samInfo->{sex};
	#if (-f $file) {
	#	&get_list($file,\%listH,\%term2oriid,"disease");
	#}
	## distribution.list.order$sex
	$file="$db2xml/distribution.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"distribution");
	}
	## metabolism.list.order$sex
	#$file="$db2xml/metabolism.list.order".$samInfo->{sex};
	#if (-f $file) {
	#	&get_list($file,\%listH,\%term2oriid,"metabolism");
	#}
	## metasummary.list.order$sex
	#$file="$db2xml/metasummary.list.order".$samInfo->{sex};
	#if (-f $file) {
	#	&get_list($file,\%listH,\%term2oriid,"metasummary");
	#}
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
				$final_report->{$head[$i]} = $unit[$i] if ($head[$i] eq "barcode");
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
			}else{
				$head[$i]=decode("UTF-8",$head[$i]);
				$head[$i]="2型糖尿病" if ($head[$i] =~ /糖尿病/);
				$head[$i]="炎症性肠病" if ($head[$i] =~ /肠炎/ || $head[$i] =~ /炎症性肠病/ || $head[$i] =~ /腹泻/);
				$head[$i]="粪便颜色" if ($head[$i] =~ /粪便颜色/ || $head[$i] =~ /颜色/);
				if ($head[$i] eq "粪便质地" || $head[$i] eq "粪便颜色" || $head[$i] =~ /历次检测编号/){
					$samInfo{$head[$i]}=decode("UTF-8",$unit[$i]);
				}else{
					$samInfo{disease}{$head[$i]}{'risk'}=decode("UTF-8",$unit[$i]);
				}
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

USAGE
	print $usage;
	exit;
}
