#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(cwd);
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Statistics::Basic qw(:all);
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
my ($infile,$referdir,$outfile,$level_cfg,$method_cfg,$overview_cfg,$help);
GetOptions (
	"infile=s"	=> \$infile,
	"referdir=s"	=> \$referdir,
	"outfile=s"	=> \$outfile,
	"level_cfg=s"	=> \$level_cfg,
	"method_cfg=s"	=> \$method_cfg,
	"overview_cfg=s"	=> \$overview_cfg,
	"h"	=> \$help,
) or &USAGE;
&USAGE if ((!defined $infile || !defined $outfile) || defined $help) ;

########## AbsolutePath and default ##############
$referdir ||= "$Bin/../../reference_data/MMC_current/";
$level_cfg ||= "$Bin/conf/level.conf";
$method_cfg ||= "$Bin/conf/term2method.conf";
$overview_cfg ||= "$Bin/conf/overview.conf";
$infile=&AbsolutePath("file",$infile);
$outfile=&AbsolutePath("file",$outfile);
$referdir=&AbsolutePath("dir",$referdir);
$level_cfg=&AbsolutePath("file",$level_cfg);
$method_cfg=&AbsolutePath("file",$method_cfg);
$overview_cfg=&AbsolutePath("file",$overview_cfg);

#$level_set||="N";

################# main process ###################
## check file and dir
&check_file($infile,$referdir,$level_cfg,$method_cfg,$overview_cfg);


## load method and level config file
my ($method_info) = &load_list($method_cfg);
my ($level_info) = &load_list($level_cfg);

## get level boundary arry
my ($level_val) = &get_val($level_info);

## overviews term composition of distribution oriid and weigth
my ($overviewcomposition,$compositionweigth) = &load_overview($overview_cfg);

## calculate meta level and per by sample data and ref data
#my ($Data)=&calculate_meta_level_per($infile,$referdir,$level_set);
my ($Data)=&calculate_meta_level_per($infile,$referdir,$method_info,$level_info,$level_val);

&calculate_overview_level_per($Data,$overviewcomposition,$compositionweigth,$level_info,$level_val);

&calculate_total_level($Data);

&getresult($Data,$level_val,$outfile,$overviewcomposition);

##get result
#&getresult($Data,$outfile,$overviewcomposition);
sub getresult {
	my ($Data,$level_val,$outfile,$overviewcomposition) = @_;

	open OUT,">$outfile";
	print OUT "#term2bac\tvalue\tlevel\tpercent";
	
	foreach my $m(@$level_val){
		print OUT "\tval$m";
	}
	print OUT "\n";

	foreach my $meta (sort keys %{$Data}){
		if($Data->{$meta}{per} ne "-"){
			if($meta eq 'BEpercent' || exists $overviewcomposition->{$meta}){
				print OUT "$meta\tNA\tNA\t$Data->{$meta}{per}\t";
				print OUT "NA\tNA\tNA\tNA\tNA\n";
			}else{
				print OUT "$meta\t$Data->{$meta}{sample}\t$Data->{$meta}{level}\t$Data->{$meta}{per}\t";
				print OUT join("\t",@{$Data->{$meta}{val}}),"\n";
			}
		}
	}
	close OUT;
}

sub calculate_total_level {
	my ($Data) = @_;
	my $BEpercent;
	my $Enterobacter = 0.00001;	

	if($Data->{Enterobacter}{sample} > 0){
		$BEpercent = sprintf "%.2f",($Data->{Bifidobacterium}{sample} / $Data->{Enterobacter}{sample});
	}else{
		$BEpercent = sprintf "%.2f",($Data->{Bifidobacterium}{sample} / $Enterobacter);
	}

	if($BEpercent > 100){
		$BEpercent = 100;
	}

	$Data->{BEpercent}{per} = $BEpercent;

	foreach my $meta (sort keys %{$Data}){
		#print $meta,"\n";
		
	}
}

#($overviewcomposition,$compositionweigth) = &load_overview($overview_cfg);
sub load_overview {
	my ($overview_cfg) = @_;
	my %overviewcomposition = ();
	my %compositionweigth = ();

	open IN,"$overview_cfg";
	while(<IN>){
		chomp;
		next if (/^$/ || /^\#/);
		my @unit = split/\t/,$_;
		$overviewcomposition{$unit[0]}{$unit[1]} = $unit[2];
		$compositionweigth{$unit[1]} = $unit[3];
	}
	close IN;

	return (\%overviewcomposition,\%compositionweigth);
}

#&load_method($method_cfg,\%method_info);
sub load_list {
	my ($cfg) = @_;
	my %info=();

	open IN,"$cfg";
	while(<IN>){
		chomp;
		next if/^#/;
		my @unit = split/\t/,$_;
		$info{$unit[0]} = $unit[1];
	}
	close IN;

	return (\%info);
}

#&check_file($infile,$referdir,$level_cfg,$method_cfg,$overview_cfg);
sub check_file {
	my ($infile,$referdir,$level_cfg,$method_cfg,$overview_cfg) = @_;
	
	my $file;
	my $check=0;

	## check reference_population
	#
	$file="$referdir/refer.extract.xls";
	$check+=&checkF($file);
	#
	#$file="$referdir/refer.extract.stat.xls";
	#$check+=&checkF($file);
	
	## check input data file
	#
	$file=$infile;
	$check+=&checkF($file);
	
	## check level config file
	#
	$file=$level_cfg;
	$check+=&checkF($file);

	## check method config file
	#
	$file=$method_cfg;
	$check+=&checkF($file);

	## check overview config file
	#
	$file=$overview_cfg;
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

#&calculate_overview_level_per($Data,$overviewcomposition,$compositionweigth,$level_info,$level_val);
sub calculate_overview_level_per {
	my ($Data,$overviewcomposition,$compositionweigth,$level_info,$level_val) = @_;

	foreach my $overview_term(keys %$overviewcomposition){
		my $item_per = 0;
		my $i = 0;
		foreach my $term(sort keys %{$overviewcomposition->{$overview_term}}){
			if(exists $Data->{$term}){
				$i++;
				if($overviewcomposition->{$overview_term}{$term} eq '+'){
					$item_per = $item_per + $Data->{$term}{per} * $compositionweigth->{$term}
				}
				if($overviewcomposition->{$overview_term}{$term} eq '-'){
					$item_per = $item_per + (1 - $Data->{$term}{per} * $compositionweigth->{$term});
				}
			}
		}
		$Data->{$overview_term}{per} = $item_per/$i;
	}
}


#my ($Data)=&calculate_meta_level_per($infile,$referdir,\%method_info,\%level_info,$level_val);
sub calculate_meta_level_per {
	my ($infile,$referdir,$method_info,$level_info,$level_val) = @_;
	my $Data=();
	
	## load sample and reference data
	($Data)=&load_sample_reference_data($infile,$referdir,$method_info);
	
	## get level,per,range
	&calculate_level_per($Data,$method_info,$level_info,$level_val);
	
	return ($Data);
}

#&calculate_level_per($Data,$level_set,$level_val)
sub calculate_level_per {
	my ($Data,$method_info,$level_info,$level_val)=@_;
	
	##
	foreach my $meta (keys %{$Data}) {

		# get percent
		($Data->{$meta}{per},$Data->{$meta}{val})=&get_per_new($Data->{$meta}, $meta, $method_info, $level_info, $level_val);
	
		# get level
		#($Data->{$meta}{level},$Data->{$meta}{range})=&get_level($Data->{$meta},$meta,$level_set);
		($Data->{$meta}{level})=&get_level_new($Data->{$meta},$meta,$level_info);
		
		## get per
		#$Data->{$meta}{per}=&get_per($Data->{$meta}, $meta);
	}
	
	return;
}

#($Data->{$meta}{per},$Data->{$meta}{val)=&get_per_new($Data->{$meta}, $meta, $method_info, $level_info, $level_val);
sub get_per_new {
	my ($meta, $key, $method_info, $level_info, $level_val) = @_;
	my ($per,$val) = ('-',0);
	my @values = ();
	my $weight = 1;
	
	next if (!exists $method_info->{$key});

	## use R pnorm
	if($method_info->{$key} eq 'pnorm' && (defined $meta->{ref}{mean} && $meta->{ref}{mean} != 0) && (defined $meta->{ref}{sd} && $meta->{ref}{sd} != 0) && exists $meta->{sample}){
		$R->set('v', $meta->{sample});
		$R->set('y', $meta->{ref}{mean});
		$R->set('z', $meta->{ref}{sd});
		$R->run(q`prob <- pnorm(v, y, z)`);
		$per = $R->get('prob');
		
		foreach my $v(@$level_val){
			my $weight_val = $v * $weight * 0.01;
			$R->set('m', $weight_val);
			$R->run(q`value <- as.numeric(qnorm(m,y,z))`);
			$val = $R->get('value');
			push @values,$val;
		}
	
	}

	## use ecdf calculate percentile
	if($method_info->{$key} eq 'ecdf' && exists $meta->{sample} && exists $meta->{ref}){
		$R->set('v', $meta->{ref}{values});
		$R->set('y', $meta->{sample});
		$R->run(q`percentile <- ecdf(v)(y)`);
		$per = $R->get('percentile');

		foreach my $v(@$level_val){
			my $weight_val = $v * $weight * 0.01;
			$R->set('m', $weight_val);
			$R->run(q`value <- as.numeric(quantile(v,m))`);
			$val = $R->get('value');
			push @values,$val;
		}
	}

	return ($per,\@values);
}

sub get_val {
	my ($level_info) = @_;
	my @val = ();

	foreach my $k(sort keys %$level_info){
		my $value = (split/-/,$level_info->{$k})[1];
		push @val,$value;
	}

	return (\@val);
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

#($Data->{$meta}{level})=&get_level_new($Data->{$meta},$meta,$level_info);
sub get_level_new {
	my ($meta,$key,$level_info) = @_;
	my $level = "-";

	foreach my $k(sort keys %$level_info){
		my ($L1,$L2) = split/-/,$level_info->{$k};	
		if (exists $meta->{per} && $meta->{per} ne '-'){
			my $per = $meta->{per}*100;
			if($per == $L2){
				$level = $k;
			}
			if($per >= $L1 && $per < $L2){
				$level = $k;
				next;
			}
		}
	}
	
	return ($level);	
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
		$range=(sprintf "%.2f",$mean-2*$sd).'~'.(sprintf "%.2f",$mean+2*$sd) if ($k eq "diversity");
#		$range=&replace($range);
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

#($Data)=&load_sample_reference_data($infile,$referdir,$method_info);
sub load_sample_reference_data {
	my ($infile,$referdir,$method_info) = @_;
	my %Data=();

	my $file;	
	## get sample values
	&get_sample_values($infile,\%Data,$method_info);
	
	## get ref samples values && sort
	$file="$referdir/refer.extract.xls";
	&get_ref_values($file,\%Data,$method_info);

	## get ref mead && sd
	#$file="$referdir/refer.extract.stat.xls";
	#&get_ref_mean_sd($file,\%Data,$method_info);

	return (\%Data); 
}

#&get_sample_values($infile,\%Data);
sub get_sample_values {
	my ($file,$Data,$method_info)=@_;

	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$value)=(split /\t/,$_)[0,1];
		next if (!exists $method_info->{$id});
		$Data->{$id}{sample}=$value;
	}
	close IN;
	
	foreach my $term(sort keys %$method_info){
		if(!exists $Data->{$term}){
			$Data->{$term}{sample} = "0";
		}
	}

	return;
}

#&get_ref_mean_sd($file,\%Data);
sub get_ref_mean_sd {
	my ($file,$Data,$method_info)=@_;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @unit=split /\t/,$_;
		next if (! exists $method_info->{$unit[0]});
		$Data->{$unit[0]}{ref}{mean}=$unit[3];
		$Data->{$unit[0]}{ref}{sd}=$unit[4];
	}
	close IN;
	
	return;
}

#&get_ref_values($file,\%Data);
sub get_ref_values {
	my ($file,$Data,$method_info)=@_;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @unit=split /\t/,$_;
		next if (!exists $method_info->{$unit[0]});
		my $zeroN=0;
		for (my $i=1;$i<@unit;$i++) {
			push @{$Data->{$unit[0]}{ref}{values}},$unit[$i];
			$zeroN++ if ($unit[$i] == 0);
		}
		# array sort
		@{$Data->{$unit[0]}{ref}{values}}=sort {$a<=>$b} @{$Data->{$unit[0]}{ref}{values}};
		# calculate mean and sd
		($Data->{$unit[0]}{ref}{mean},$Data->{$unit[0]}{ref}{sd}) = &STAT_(\@{$Data->{$unit[0]}{ref}{values}}) if ($Data->{$unit[0]}{ref}{values}[-1] > 0);
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

#($Data->{$unit[0]}{ref}{mean},$Data->{$unit[0]}{ref}{sd}) = &STAT_(\@{$Data->{$unit[0]}{ref}{values}}) if ($Data->{$unit[0]}{ref}{values}[-1] > 0);
sub STAT_ {
	my ($values) = @_;

	my $mean = mean(@$values);
	my $sd = stddev(@$values);
	$mean = sprintf "%e", $mean;
	$sd = sprintf "%e", $sd;

	return ($mean,$sd);
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

    Program Function: calculate GIhealth result ;
    Version:    $version
    Usage:
      Options:
	-infile		<file>	sample's report result file			forced

	-outfile	<file>	output file					forced

	-referdir	<dir>	reference pop data dir				optional
                          default "$Bin/../../reference_data/MMC_current/";

	-level_cfg	<file>	level config file				optional
                          default "$Bin/conf/level.conf"

	-method_cfg	<file>	calculate percent method			optional
			  default "$Bin/conf/method.conf"
	
	-ovewview_cfg	<file>	
			  defaul "Bin/conf/overview.conf"

	-h		Help

USAGE
	print $usage;
	exit;
}
