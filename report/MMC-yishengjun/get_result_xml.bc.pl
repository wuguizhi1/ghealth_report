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
my ($incsv,$infile,$indir,$dirout,$barcode,$language,$help,$host,$port,$database,$collection,$updateadvice,$updateall);
GetOptions (
	"csvfile=s"	=> \$incsv,
	"infile=s"	=> \$infile,
	"in=s"		=> \$indir,
	"out=s"		=> \$dirout,
	"barcode=s"	=> \$barcode,
	"language=s"	=> \$language,

	"h"	=> \$help,
	"host=s"        => \$host,
	"port=s"        => \$port,
	"database=s"    => \$database,
	"collection=s"  => \$collection,
	"updateadvice"  => \$updateadvice,
	"updateall"     => \$updateall,
) or &USAGE;
&USAGE if ((!defined $incsv || !defined $infile) || defined $help) ;

########## AbsolutePath and default ##############
$dirout||="./";
$indir||="$Bin/db2xml/current/";
mkdir $dirout if (! -d $dirout);
$incsv=&AbsolutePath("file",$incsv);
$infile=&AbsolutePath("file",$infile);
$dirout=&AbsolutePath("dir",$dirout);
$language||="CN";
$host = $host || "10.0.0.204";
$port = $port || "27021";
$database = $database || "result";
$collection = $collection || "mmccdjq";

################# main process ###################
my %vars=();
my %json=();
my %XML=();
#$json{_id} = $barcode;
$json{version} = "秒秒测有益菌检测";
## check file and dir
&check_file($incsv,$infile);

## load sample information file
my ($samInfo)=&load_sample_info($incsv,\%vars,\%json,\%XML);

## load $type.list.order$sex
my ($listH,$term2oriid)=&load_list($indir,$samInfo);
##load infile
my($Data)=&load_infile($infile,$term2oriid);
## main process
&main_process($Data,$listH,\%vars,$indir,$language,$barcode,$samInfo,\%json,\%XML);
## sample result csv for erp import
&sample_result_csv(\%vars,\%XML,$barcode,$dirout);
#my %final_report;
#&get_final_report(\%vars,\%final_report,$barcode,$dirout);

## insert final report to mongo db
my $mongo = MongoDB::MongoClient->new('host' => $host.':'.$port);
my $db = $mongo->get_database($database);
my $cl = $db->get_collection($collection);
my $barcode_find = $cl->find({_id => $barcode});
my $flag = 0;
while (my $barcode_report = $barcode_find->next()) {
	$flag = 1;
	if (!exists $barcode_report->{Report}{page3} && exists $json{page3}) {
		$cl->update_one({_id => $barcode}, {'$set' => {'Report.page3' => \%{$json{page3}} }}, {'upsert' => 1});
	}
}
if ($flag == 0 || defined $updateall) {
	$cl->update_one({_id => $barcode}, {'$set' => {'Report' => \%json }},{'upsert' => 1});
}

my $JSON = encode_json \%json;
open JSoN,">$dirout/$barcode.json";
print JSoN"$JSON";
close JSoN;
############################# sub ###############################################
sub load_infile{
	my($infile,$term2oriid) = @_;
	my %Data=();
	open IN,$infile or die $!;
	chomp(my $head = <IN>);
	my @line = split/\t/,$head;
	while(<IN>){
		chomp;
		my @unit = split /\t/,$_;
		if(exists $term2oriid->{$unit[0]}){
			$unit[0] = $term2oriid->{$unit[0]};
		}
		for(my $i=1;$i<@unit;$i++){
			$Data->{$unit[0]}{$line[$i]} = $unit[$i];
		}
	}
	close IN;
	return ($Data);
}


#&main_process($Data,$listH,\%vars,$infile,$indir,$language,$barcode,$samInfo,$json);
sub main_process {
	my ($Data,$listH,$vars,$indir,$language,$barcode,$samInfo,$json,$XML)=@_;

	## 菌群分布(有益菌属分布)
	# $listH->{"distribution"}; $listH->{"pathogens"}
	&Intestinal_Flora($Data,$listH->{distribution},$vars,"$indir/XML/distribution/$language",$language,$json,$XML);

	## 菌群综合指标(有益菌)
	# $listH->{"overview"}
	&Micro_Overvwiew($Data,$listH->{overview},$vars,"$indir/XML/overview/$language",$language,$json,$XML);

	## 膳食方案+肠道调节方案+运动方案
	# $listH->{"advise"};
	&Program($vars,"$indir/XML/advise",$language,$samInfo,$listH->{advise},$json,$XML);

	## 多次检测结果比较
	# $samInfo->{batch}; $samInfo->{'历次检测编号'}
	# &Compare_history($vars,$samInfo,$historydir) if ($samInfo->{batch} > 1);

	## 报告日期
	$XML->{reportdate}=strftime("%Y-%m-%d",localtime());
	
	## output xml
	&output_xml($XML,$barcode,$dirout);
	
	return;
}


#&Program($vars,"$indir/XML",$language,$samInfo,$listH->{diseaseList});
sub Program {
	my ($vars,$dir,$language,$samInfo,$list,$json,$XML)=@_;
	#my ($foodKey,$regulationKey,$sportKey)=('-','-','-');
	my $ageKey=&get_program_age($samInfo->{age});
	## process food and others
	&process_food($vars,$dir,$ageKey,$language,$json,$XML);
	## process regulation
	&process_regulation($vars,$dir,$ageKey,$language,$json,$XML);
	##process exercise
	#&process_exercise($vars,$dir,$ageKey,$language,$json,$XML);
	##summarize
	&process_summarize($vars,$XML,$json);
	## process probiotics
#	&process_probiotics($vars,$dir,$ageKey,$samInfo->{sex},$language);
	## process prebiotics
#	&process_prebiotics($vars,$dir,$ageKey,$samInfo->{sex},$language);
	
	return;
}

#&process_summarize($vars,$json);
sub process_summarize{
	my($vars,$XML,$json) = @_;
	$XML->{fangan}{summary} = "综合您本次检测结果，肠道有益菌总得分为$XML->{gaikuang}{desc}{youyijun}{score}分，";
	if($XML->{gaikuang}{desc}{youyijun}{risk} eq "有害"){
		$XML->{fangan}{summary} .= "不利于肠道及人体健康，";
		if ($vars->{changdaojunqun}{fenbu}{badN} == 4){
			$XML->{fangan}{summary} .= "4种有益菌属的检测指数均偏低，需注意。";
		}
		if ($vars->{changdaojunqun}{fenbu}{badN} > 0 && $vars->{changdaojunqun}{fenbu}{badN} < 4){
			$XML->{fangan}{summary} .= "其中，$vars->{changdaojunqun}{fenbu}{biotics_bad}检测指数偏低，需注意。";
		}
	}
	elsif($XML->{gaikuang}{desc}{youyijun}{risk} eq "正常"){
		$XML->{fangan}{summary} .= "处于正常范围，";
		if ($vars->{changdaojunqun}{fenbu}{badN} == 4){
			$XML->{fangan}{summary} .= "4种有益菌属的检测指数均偏低，需注意。";
		}
		if ($vars->{changdaojunqun}{fenbu}{badN} > 0){
			$XML->{fangan}{summary} .= "但其中$vars->{changdaojunqun}{fenbu}{biotics_bad}检测指数偏低，可重点调理。";
		}
	}
	else{
		$XML->{fangan}{summary} .= "有利于肠道及人体健康，";
		if ($vars->{changdaojunqun}{fenbu}{badN} == 4){
			$XML->{fangan}{summary} .= "4种有益菌属的检测指数均偏低，需注意。";
		}
		if ($vars->{changdaojunqun}{fenbu}{badN} > 0){
			$XML->{fangan}{summary} .= "但其中$vars->{changdaojunqun}{fenbu}{biotics_bad}检测指数偏低，可适当改善。";
		}
	}
	$json->{page3}{summarize} = $XML->{fangan}{summary};

	return;
}

#&process_regulation($vars,$dir,$regulationKey,$language);
sub process_regulation{
	my($vars,$dir,$ageKey,$language,$json,$XML) = @_;
	##process prebiotics
	&process_prebiotics($vars,$dir,$ageKey,$language,$json,$XML);
	##process probiotics
	&process_probiotics($vars,$dir,$ageKey,$language,$json,$XML);
	
	return;
}

sub process_probiotics{
	my ($vars,$dir,$ageKey,$language,$json,$XML)=@_;
	## 补充益生菌	
	my $vd=XMLin("$dir/$language/tjsjcyishengjunfangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vd->{$language}{title} =~ s/方案//;
	$XML->{fangan}{probiotics}{name} = $vd->{$language}{title};
	if($vars->{changdaojunqun}{fenbu}{probiotics_sup} eq '-'){
		push @{$XML->{fangan}{probiotics}{suggest}},$vd->{$language}{agesuggestion}{desc}{$ageKey}{action};
		push @{$json->{page3}{probiotics}} ,$vd->{$language}{agesuggestion}{desc}{$ageKey}{action};
	}
	else{
		push @{$XML->{fangan}{probiotics}{suggest}},"建议您选择补充"."$vars->{changdaojunqun}{fenbu}{probiotics_sup}"."等益生菌产品，以增加您肠道内益生菌的含量，并能抑制有害菌的异常增殖，调节肠道菌群平衡。";
		push @{$json->{page3}{probiotics}},"建议您选择补充"."$vars->{changdaojunqun}{fenbu}{probiotics_sup}"."等益生菌产品，以增加您肠道内益生菌的含量，并能抑制有害菌的异常增殖，调节肠道菌群平衡。";
	}

	return;
}

sub process_prebiotics{
	my ($vars,$dir,$ageKey,$language,$json,$XML)=@_;
	## 补充益生元
	my $vd=XMLin("$dir/$language/tjsjcyishengyuanfangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vd->{$language}{title} =~ s/方案//;
	$XML->{fangan}{prebiotics}{name} = $vd->{$language}{title};
	if($vars->{changdaojunqun}{fenbu}{prebiotics_sup} eq '-'){
		push @{$XML->{fangan}{prebiotics}{suggest}} ,$vd->{$language}{agesuggestion}{desc}{$ageKey}{action};
		push @{$json->{page3}{prebiotics}},$vd->{$language}{agesuggestion}{desc}{$ageKey}{action};
	}
	else{
		push @{$XML->{fangan}{prebiotics}{suggest}},"建议您选择补充"."$vars->{changdaojunqun}{fenbu}{prebiotics_sup}"."等益生元，能帮助"."$vars->{changdaojunqun}{fenbu}{biotics_bad}"."在肠道内增殖，增加有益物质的产生，有利于增强肠道屏障功能，促进肠道健康。";
		push @{$json->{page3}{prebiotics}},"建议您选择补充"."$vars->{changdaojunqun}{fenbu}{prebiotics_sup}"."等益生元，能帮助"."$vars->{changdaojunqun}{fenbu}{biotics_bad}"."在肠道内增殖，增加有益物质的产生，有利于增强肠道屏障功能，促进肠道健康。";
	}
	

	return;
}

#&process_food($vars,$dir,$ageKey,$samInfo->{sex},$language);
sub process_food {
	my ($vars,$dir,$ageKey,$language,$json,$XML) = @_;
	## food
	my $vd_food=XMLin("$dir/$language/tjsjcshanshifangan.xml",NoAttr=>1,SuppressEmpty => "");
	$vd_food->{$language}{title} =~ s/方案//;
	$XML->{fangan}{food}{name} = $vd_food->{$language}{title};
	if ($vars->{changdaojunqun}{fenbu}{food_sup} eq '-'){
		push @{$XML->{fangan}{food}{suggest}} , $vd_food->{$language}{agesuggestion}{desc}{$ageKey}{action};
		push @{$json->{page3}{food}},$vd_food->{$language}{agesuggestion}{desc}{$ageKey}{action};
	}
	else {
		push @{$XML->{fangan}{food}{suggest}} , "建议您在日常饮食中选择补充"."$vars->{changdaojunqun}{fenbu}{food_sup}"."等食物或者相关的加工提取食品，为肠道菌群提供更多营养，促进有益共生菌的生长，增加肠道菌群的丰富度及多样性。";
		push @{$json->{page3}{food}},"建议您在日常饮食中选择补充"."$vars->{changdaojunqun}{fenbu}{food_sup}"."等食物或者相关的加工提取食品，为肠道菌群提供更多营养，促进有益共生菌的生长，增加肠道菌群的丰富度及多样性。";

	}
#	if ($vd_food->{$language}{agesuggestion}{desc}{$ageKey}{note} ne '') {
#		push @{$XML->{fangan}{food}{suggest}} , $vd_food->{$language}{agesuggestion}{desc}{$ageKey}{note};
#
#		push @{$json->{page3}{food}},$XML->{fangan}{food}{suggest};
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

#&Micro_Overvwiew($Data,$listH->{overview},$vars,"$indir/XML/overview/$language",$language,$json);
sub Micro_Overvwiew {
	my ($Data,$list1,$vars,$dir1,$language,$json,$XML)=@_;
	
	##
	$XML->{gaikuang}{name}="肠道菌群概况";
	$XML->{gaikuang}{ename}="OVERVIEW";

	##
	&item_Micro_Overvwiew('desc',$Data,$list1,$vars,$dir1,$language,$json,$XML);

	return;
}


#&item_Micro_Overvwiew('desc',$Data,$list,$vars,$dir,$language,$json);
sub item_Micro_Overvwiew {
	my ($prefix,$Data,$list,$vars,$dir,$language,$json,$XML)=@_;
	
	##
	my (@good,@bad);
	##
#	$vars->{gaikuang}{$prefix}{abnormalN}=0;
#	$vars->{gaikuang}{$prefix}{goodN}=0;
#	$vars->{gaikuang}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta});
		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		
		## get one item
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta,$vars->{age},$json,$XML);
		$XML->{gaikuang}{$prefix}{youyijun}=$item;

	}
}

#&Intestinal_Flora($Data,$listH->{distribution},$vars,"$indir/XML/distribution/$language",$language,$json);
sub Intestinal_Flora {
	my ($Data,$list,$vars,$dir,$language,$json,$XML)=@_;
	##
	$XML->{changdaojunqun}{name}="肠道菌群";
	$XML->{changdaojunqun}{ename}="MICROBIOTA";

	## 分布
	&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language,$json,$XML);
	
	return;
}

#&item_Intestinal_Flora('fenbu',$Data,$list,$vars,$dir,$language);
sub item_Intestinal_Flora {
	my ($prefix,$Data,$list,$vars,$dir,$language,$json,$XML)=@_;
	##
	my (@good,@bad,@gooditem,@baditem);
	my (@goodDesc,@badDesc,@goodDisease,@goodSource,@badDisease,@goodEffect,@badEffect,@badSource);
	my (@biotics_good,@biotics_bad,@probiotics_sup,@prebiotics_sup);
	my (@food_sup_all,@food_sup_bad1,@food_sup_bad2,@food_sup_bad3);
	##
#	$vars->{changdaojunqun}{$prefix}{abnormalN}=0;
#	$vars->{changdaojunqun}{$prefix}{goodN}=0;
	$vars->{changdaojunqun}{$prefix}{badN}=0;
	##
	foreach my $meta (sort {$list->{$a}{order} <=> $list->{$b}{order}} keys %{$list}){
		next if (! exists $Data->{$meta} || ! exists $Data->{$meta}{sample});
		#print "$meta\n";
		## get $vd
		my $vd=XMLin("$dir/$meta.xml",NoAttr=>1,SuppressEmpty => "");
		#$itemcount ++;
		my $effect;
		my ($item)=&get_item($vd,$Data->{$meta},$language,$meta,$vars->{age},$json,$XML);
		if ($item->{name} =~ /(.*)（(.*)）/) {
			$item->{cnname} = $1;
			$item->{enname} = $2;
		}
		push @{$XML->{changdaojunqun}{$prefix}{item}},$item;
		## get abnormalN,badN,good,bad
		my ($add,$addN,$addG)=&get_good_bad_array($item,$Data->{$meta},\@good,\@bad);
#		$vars->{gaikuang}{$prefix}{abnormalN}+=$add;
#		$vars->{gaikuang}{$prefix}{goodN}+=$addG;
		$vars->{changdaojunqun}{$prefix}{badN}+=$addN;
		## get goodDesc,badDesc
		if (defined $Data->{$meta}{effect} && $Data->{$meta}{effect} eq "有害") {
			push @baditem,$item;
			if($item->{name}=~/(.*)\（.*/){
				unshift @biotics_bad,$1;
			}
			else{
				unshift @biotics_bad,$item->{name};
			}
			&process_PreProbiotics_and_food($vars->{age},\@prebiotics_sup,\@probiotics_sup,\@food_sup_all,\@food_sup_bad1,\@food_sup_bad2,\@food_sup_bad3,$vd->{$language}{suggestion}{healthy}{$Data->{$meta}{level}}{actiondesc}) if (defined $vd->{$language}{suggestion}{healthy}{$Data->{$meta}{level}}{actiondesc} && $vd->{$language}{suggestion}{healthy}{$Data->{$meta}{level}}{actiondesc} ne "");
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
#	push @{$vars->{changdaojunqun}{$prefix}{abnormalitem}},@baditem,@gooditem;
	## good,bad
#	($vars->{changdaojunqun}{$prefix}{good},$vars->{changdaojunqun}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
	## goodSource,badSource
#	($vars->{changdaojunqun}{$prefix}{goodSource},$vars->{changdaojunqun}{$prefix}{badSource})=&process_good_badSource(\@goodSource,\@badSource);

	$vars->{changdaojunqun}{$prefix}{biotics_bad}=&process_good_bad(\@biotics_bad);
	($vars->{changdaojunqun}{$prefix}{prebiotics_sup},$vars->{changdaojunqun}{$prefix}{probiotics_sup},$vars->{changdaojunqun}{$prefix}{food_sup})=&process_bad_sups($prebiotics_sup_uniq,$probiotics_sup_uniq,$food_sup_uniq);
	## goodDisease,badDisease,goodEffect,badEffect
#	($vars->{changdaojunqun}{$prefix}{goodDisease},$vars->{changdaojunqun}{$prefix}{badDisease},$vars->{changdaojunqun}{$prefix}{goodEffect},$vars->{changdaojunqun}{$prefix}{badEffect})=&process_goodbad_DiseaseEffect(\@goodDisease,\@badDisease,\@goodEffect,\@badEffect);
	## badDesc,goodDesc
#	($vars->{changdaojunqun}{$prefix}{goodDesc},$vars->{changdaojunqun}{$prefix}{badDesc})=&process_good_badDesc(\@goodDesc,\@badDesc);

	return;
}
#&process_PreProbiotics_and_food($vars->{age},\@prebiotics_sup,\@probiotics_sup,\@food_sup_all,\@food_sup_bad1,\@food_sup_bad2,\@food_sup_bad3,$vd->{$language}{suggestion}{healthy}{$item->{level}}{actiondesc}
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
	my ($vd,$Data,$language,$key,$age,$json,$XML)=@_;
	my %item=();
	if($key eq "tjsjcyouyijun"){
		$item{name} = $vd->{$language}{title};
		$item{score} = int($Data->{per}*100);
		$item{score} = "1" if ($Data->{per}*100 < 1);
		$item{score} = "99" if ($Data->{per}*100 >= 100);
		$json->{page1}{score} = $item{score};
		if($Data->{level} ne '-'){
			$item{risk} = $vd->{$language}{suggestion}{healthy}{$Data->{level}}{effect};
			$item{riskdesc}=$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanism};
			push @{$item{riskdescfull}},"综合您本次肠道有益菌检测，您的粑粑健康跑分为$item{score}分。";
			push @{$item{riskdescfull}},$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanismfull};
			
			my @mechanismfull = split /\s+/,$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanismfull},2;
			push @{$json->{page1}{result}},"综合您本次肠道有益菌检测，您的粑粑健康跑分为$item{score}分。";
			push @{$json->{page1}{result}} , $mechanismfull[0];
			
			$json->{page1}{suggest} = $mechanismfull[1];
		}
	}
	else{
		$Data->{effect} = $vd->{$language}{suggestion}{healthy}{$Data->{level}}{effect};
#		print "$Data{Val30}\n";
		$item{name} = $vd->{$language}{title};
		$item{value} = (int($Data->{sample}*100000))/100;
		$item{threshold} = $Data->{Val30};
		$item{desc} = $vd->{$language}{summary}{desc};
#		$item{descfull} = $vd->{$language}{summary}{descfull};
		
		if ($Data->{Val30} == 0){
			$item{threshold} = 0.01;
			if($Data->{sample} < $Data->{Val30}){
				$Data->{level} = "l2";
			}
		}
		my %part=();
		$part{name} = $item{name};
		$part{value} = $item{value};
		$part{threshold} = $item{threshold};
		if($Data->{level} eq "l3"){
			push @{$item{riskdescfull}},"您本次检测结果显示，您肠道内的$part{name}含量正常。";
		#	$part{result}=$item{riskdescfull};
		}
		else{
			if($item{value}==0){
				push @{$item{riskdescfull}},"您本次检测结果显示，您肠道内的$part{name}未检出。";
				push @{$item{riskdescfull}},$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanismfull};
			}
			else{
				push @{$item{riskdescfull}},"您本次检测结果显示，您肠道内的$part{name}含量较低。";
				push @{$item{riskdescfull}},$vd->{$language}{suggestion}{healthy}{$Data->{level}}{mechanismfull};
			}
		}
		$part{result} = $item{riskdescfull};
		$part{function} = $item{desc} ;
		push @{$json->{page2}},\%part;
	}

	return (\%item);
}

#my ($add,$addN,$addG)=&get_good_bad_array($item,\@good,\@bad);
sub get_good_bad_array {
	my ($item,$Data,$good,$bad)=@_;
	my ($add,$addN,$addG)=(0,0,0);
#	for my $key (%{$item}){
#		print "$key\n";}
	##
	if (! defined $Data->{effect}) {
		$add=1;
	}
	else {
		if ($Data->{effect} ne "正常") {
			$add=1;
			if ($Data->{effect} eq "有益") {
				push @{$good},$item->{name};
				$addG=1;
			}
			elsif ($Data->{effect} eq "有害") {
				push @{$bad},$item->{name};
				$addN=1;
			}
		}
	}
	return ($add,$addN,$addG);
}

#($vars->{yingyanggongneng}{desc}{$prefix}{good},$vars->{yingyanggongneng}{desc}{$prefix}{bad})=&process_good_bad(\@good,\@bad);
sub process_good_bad {
	my ($badA)=@_;
	my ($bad)=('-');
	##
#	if (@{$goodA}) {
#		if (scalar @{$goodA} <= 5) {
#			$good=join "、",@{$goodA};
#		}
#		else {
#			$good=join("、",$goodA->[0],$goodA->[1],$goodA->[2],$goodA->[3],$goodA->[4]);
			#$good.="等";
#		}
#	}
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
	
	return ($bad);
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

#&check_file($infile,$indir,$referdir,$confdir,$barcode);
sub check_file {
	my ($incsv,$infile)=@_;
	
	my $file;
	my $check=0;

	## check infile
	#
	$file=$incsv;
	$check+=&checkF($file);
	#
	$file=$infile;
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

#my ($samInfo)=&load_sample_info($infile,\%vars,\%json);
sub load_sample_info {
	my ($incsv,$vars,$json,$XML)=@_;
	my %samInfo=();
	
	##
	my $fatN=-1;
	$/="\n";
	my $file=$incsv;
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
				$json->{$head[$i]}=$unit[$i] if ($head[$i] eq "barcode");
				$XML->{$head[$i]} = $unit[$i] if ($head[$i] eq "barcode");
				if ($head[$i] eq "sex"){
					$XML->{$head[$i]} = $unit[$i];
					$json->{$head[$i]} = $unit[$i];
				}
			}
			elsif ($head[$i] eq "name") {
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
				$json->{$head[$i]}=$unit[$i];
				$XML->{$head[$i]} = $unit[$i];
			}
			elsif ($head[$i] eq "age") {
				$samInfo{$head[$i]}=$unit[$i];
				$vars->{$head[$i]}=$unit[$i];
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$vars->{'agestr'}=$unit[$i];
				$XML->{$head[$i]} = $unit[$i];
				$json->{$head[$i]} = $unit[$i];
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
				$json->{$head[$i]}=$unit[$i];
				$XML->{$head[$i]} = $unit[$i];
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

#my ($listH,$term2oriid)=&load_list($indir,$samInfo);
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
	## overview.list.order$sex
	$file="$indir/overview.list.order".$samInfo->{sex};
	if (-f $file) {
		&get_list($file,\%listH,\%term2oriid,"overview");
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

sub sample_result_csv {
        my ($vars,$XML,$barcode,$dirout) = @_;


        my $outfile="$dirout/$barcode".".result.csv";
        open my $out, '>:encoding(utf-8)', $outfile || die $!;

        my $term_str = "检测编号\t姓名\t益生菌\t";
        my $result_str = "$vars->{barcode}\t$vars->{name}\t";

        my $level;
        $level = &transfer_item_level($Data->{tjsjcyouyijun}{level});
        $result_str .="$level($XML->{gaikuang}{desc}{youyijun}{risk})\t";

        foreach my $item (@{$vars{changdaojunqun}{fenbu}{page1}}) {
		print "$item\tyes\n";
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



#my $xml=&read_xml($xmlfile);
sub read_xml {
	my ($xmlfile)=@_;
	
	##
	my $xml=XMLin($xmlfile,NoAttr=>1,ForceArray => ['item','abnormalitem','disease','spot'], KeyAttr => []);
	
	return ($xml);
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
      -csvfile      <file>   barcode.csv                              forced

      -barcode    <str>   sample's code                             forced

      -infile	  <file>  sample's report result file for input	    forced

      -in         <dir>   product set database list and xml dir     optional
                          default "$Bin/db2xml/current";
      -out        <dir>   output file dir                           optional
                          default "./";
      -language   <str>   result language used for report           optional
                          default "CN";
      -historydir <str>   history test xml result for report dir    optional
                          default "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/Xml";
      -l          <str>   level set to check report correctness (l1~l5)          optional
                          default "N";

      -h          Help

USAGE
	print $usage;
	exit;
}
