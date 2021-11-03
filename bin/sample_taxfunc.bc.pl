#!/usr/bin/perl -w
#########################################################################
# Author: mengfei
# mail: fred_routine@163.com
# Created Time: Tue May 12 09:57:32 CST 2015
#########################################################################

use strict;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use Data::Dumper;

my $path=`pwd`;
chomp($path);

my %opts;
GetOptions(\%opts,"l=s","left=s","right=s","o=s","merge-m=s","join","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -l       :: sample lable (prefix of result files), required;
                  -left    :: read1 of sequence data, required;
                  -right   :: read2 of sequence data, optional;

                  -o    :: output dir default "./";
                  -join    :: unmerged pair reads joined and used in next analysis;
                  -merge-m :: method for merge pair end reads (flash or usearch), default usearch;



USAGE

if(!defined($opts{l}) || !defined($opts{left}) || defined($opts{h}))
{
	print $usage;
	exit;
}

my $lable=$opts{l};

my $left=&ABSOLUTE_DIR($opts{left});
my $right=&ABSOLUTE_DIR($opts{right}) if (defined $opts{right});
my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

my $merge_method=$opts{"merge-m"} || "usearch";
my $log_file="$out_dir/process.log";
system "touch $log_file";

#################
my $cmd;

# sample microbiota result
$cmd = "perl $Bin/../scripts/sample_fastq2otu_process.bc.pl -l $lable -left $left ";
$cmd .= "-right $right " if (defined $right);
$cmd .= "-o $out_dir/microbiota -m $merge_method ";
$cmd .= "-join" if ($opts{"join"});
&LOG ("sample fastq reads to microbiota..", $cmd, $log_file);
system ($cmd);

# picrust predict sample function
$cmd = "perl $Bin/../scripts/sample_fastq2func.pl -l $lable -fq $out_dir/microbiota/$lable.filter.fastq -o $out_dir/function ";
&LOG ("sample filter fastq reads to function..", $cmd, $log_file);
system ($cmd);

# merge sample's result
$cmd = "perl $Bin/../scripts/combine_json.pl -micro $out_dir/microbiota/result/$lable.microbiota.json -meta $out_dir/function/result/$lable.function.json ";
$cmd .= "-prefix $lable -o $out_dir ";
&LOG ("combine sample results to json..", $cmd, $log_file);
system ($cmd);

# extract sample's test term values and compare with reference
$cmd = "perl $Bin/../scripts/extract_sample.pl -json $out_dir/$lable.result.json -prefix $lable -o $out_dir/report/\n";
&LOG ("extract and compare term value to reference pop ..", $cmd, $log_file);
system ($cmd);

# clean and clear up analysis directory 
&clean_clearup_dir_($lable, $out_dir);
# get sample info from erp system

###########################
sub clean_clearup_dir_ {
	my ($lable, $out_dir) = @_;

	# clear up microbiota dir
	system "rm $out_dir/microbiota/$lable.unmerged_*.fastq";
	system "rm $out_dir/microbiota/$lable.merge.stat $out_dir/microbiota/$lable.merge.fastq";
	system "rm $out_dir/microbiota/$lable.assign_taxon.sort_L*";

	# clear up function dir
	system "rm -rf $out_dir/function/$lable.ko_metagenome_contributions.tab";
	system "rm -rf $out_dir/function/uclust_ref_picked_otus/";
}

sub LOG {
	my ($info,$commond,$log) = @_;

	open OUT,">>$log" || die $!;
	print OUT "#$info\n$commond\n\n";
	close OUT;
}

sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub ABSOLUTE_DIR
{
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";

	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir\n";
		exit;
	}
	chdir $cur_dir;
	return $return;
}

sub sub_format_datetime 
{#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

