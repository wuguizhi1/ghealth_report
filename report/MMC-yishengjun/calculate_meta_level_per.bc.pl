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
my ($infile,$referdir,$confdir,$level_set,$outfile,$help);

GetOptions (
	"infile=s"	=> \$infile,
	"referdir=s"	=> \$referdir,
	"confdir=s"	=> \$confdir,
	"outfile=s"	=>\$outfile,
	"l=s"		=> \$level_set,

	"h"	=> \$help,
) or &USAGE;
&USAGE if (!defined $infile || defined $help) ;

########## AbsolutePath and default ##############
$referdir||="$Bin/../../reference_data/MMC_qpcr_current/";
$confdir||="$Bin/../conf/current/";
$referdir=&AbsolutePath("dir",$referdir);
$infile=&AbsolutePath("file",$infile);
$outfile=&AbsolutePath("file",$outfile);
$confdir=&AbsolutePath("dir",$confdir);


$level_set||="N";

my %enid;
$enid{"Benificial bacteria"}=1;
$enid{"Lactobacillus"}=2;
$enid{"Bifidobacterium"}=3;
$enid{"Akkermansia"}=4;
$enid{"Faecalibacterium"}=5;

################# main process ###################
## check file and dir
&check_file($infile,$referdir,$confdir);

## load $type.list.order$sex
my ($overviewcomposition,$compositionweigth)=&load_conf($confdir);
## calculate meta level and per by sample data and ref data
my ($Data)=&calculate_meta_level_per($infile,$referdir,$compositionweigth,$level_set,\%enid);
&calculate_overview_level_per($Data,$overviewcomposition,$compositionweigth,$level_set);

##get result;
&get_result($Data,$outfile);

############################# sub ###############################################
sub get_result{
	my ($Data,$outfile) = @_;

	open OUT,">$outfile" or die $!;
	print OUT"#bac\tlevel\tVal30\tVal70\tsample\tper\n";
	for my $bac (keys %{$Data}){
		next if($bac=~/Benificial bacteria/);
		print OUT"$bac\t$Data->{$bac}{level}\t$Data->{$bac}{Val30}\t$Data->{$bac}{Val70}\t$Data->{$bac}{sample}\t$Data->{$bac}{per}\n";
	}
	if(exists $Data->{"Benificial bacteria"}){
		my $bac = "Benificial bacteria";
		if(exists $Data->{$bac}{Val30}){
			print OUT"$bac\t$Data->{$bac}{level}\t$Data->{$bac}{Val30}\t$Data->{$bac}{Val70}\t$Data->{$bac}{sample}\t$Data->{$bac}{per}\n";
		}
	
		else{
			print OUT"$bac\t$Data->{$bac}{level}\tNA\tNA\t$Data->{$bac}{sample}\t$Data->{$bac}{per}\n";
		}
	}
	close OUT;
}


#&check_file($datadir,$indir,$referdir,$confdir);
sub check_file {
	my ($datadir,$referdir,$confdir)=@_;
	
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
	$file=$infile;
	$check+=&checkF($file);
	
	
	## check config file
	#
	$file="$confdir/MMC_qpcr_overview.conf";
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


sub load_conf {
	my ($confdir)=@_;
	my %overviewcomposition=();
	my %compositionweigth=();
	
	
	my $file;
	
	## overviews term composition of distribution oriid and weigth
	$file="$confdir/MMC_qpcr_overview.conf";
	if (-f $file) {
		&get_overview_composition($file, \%overviewcomposition, \%compositionweigth);
	}

	return (\%overviewcomposition,\%compositionweigth);
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
#my ($Data)=&calculate_meta_level_per($datadir,$indir,$level_set);
sub calculate_meta_level_per {
	my ($infile,$referdir,$compositionweigth,$level_set,$bacid)=@_;
	my $Data=();
	
	## load sample and reference data
	($Data)=&load_sample_reference_data($infile,$referdir,$bacid);
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

	if ($key eq "Akkermansia" || $key eq "Lactobacillus") {
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

#($Data)=&load_sample_reference_data($datadir,$indir,$referdir,$term2oriid);
sub load_sample_reference_data {
	my ($infile,$referdir,$bacid)=@_;
	my %Data=();

	my $file;
	## get sample values
	$file=$infile;
	&get_sample_values($file,$bacid,\%Data);
	
	## get ref samples values && sort
	$file="$referdir/qpcr_refer.extract.xls";
	&get_ref_values($file,$bacid,\%Data);

	## get ref mead && sd
	$file="$referdir/qpcr_refer.extract.stat.xls";
	&get_ref_mean_sd($file,$bacid,\%Data);

	return (\%Data);
}

#&get_sample_values($file,$term2oriid,\%Data);
sub get_sample_values {
	my ($file,$bacid,$Data)=@_;

	##
	$/="\n";
	#my $benificial_val = 0;
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my ($id,$value)=(split /\t/,$_)[0,1];
		next if (!exists $bacid->{$id});
#		$id = $term2oriid->{$id};
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
	my ($file,$bacid,$Data)=@_;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @unit=split /\t/,$_;
		next if (!exists $bacid->{$unit[0]});
#		$unit[0] = $term2oriid->{$unit[0]};
		$Data->{$unit[0]}{ref}{mean}=$unit[3];
		$Data->{$unit[0]}{ref}{sd}=$unit[4];
	}
	close IN;
	
	return;
}

#&get_ref_values($file,$term2oriid,\%Data);
sub get_ref_values {
	my ($file,$bacid,$Data)=@_;
	
	##
	$/="\n";
	open IN,$file or die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @unit=split /\t/,$_;
		next if (!exists $bacid->{$unit[0]});
#		$unit[0] = $term2oriid->{$unit[0]};
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

### replace;
sub replace {
	my ($js)=@_;
	my $utfcode = '';
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
      -infile        <file>   barcode.extract.xls          forced

      -referdir   <dir>   reference pop data dir                    optional
                          default "$Bin/../../reference_data/MMC_qpcr_current/";
      -out        <file>   output file                           forced
      
      -confdir    <dir>   report conf dir                           optional
                          default "$Bin/../conf/current/";
      -l          <str>   level set to check report correctness (l1~l5)          optional
                          default "N";

      -h          Help

USAGE
	print $usage;
	exit;
}
