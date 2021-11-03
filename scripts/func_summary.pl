#!/usr/bin/perl -w
#########################################################################
# Author: mengfei
# mail: fred_routine@163.com
# Created Time: Tue May 12 09:57:32 CST 2015
#########################################################################

use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use Data::Dumper;
use JSON::XS;

my $path=`pwd`;
chomp($path);

my %opts;
GetOptions(\%opts,"id=s","p=s","o=s","func=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -id    :: formated micro results directory, required;
                  -p     :: outfiles prefix, required;

                  -o     :: output dir default "./";
                  -h     help info;

USAGE

if(!defined($opts{id}) || !defined($opts{p}))
{
	print $usage;
	exit;
}

my $indir=&ABSOLUTE_DIR($opts{id});
my $prefix=$opts{p};

my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

#################
# function stratified stat
my %func_key = (
	"C1" => "Carbohydrate Digestion & Absorption",
	"C2" => "Lipids Digestion & Absorption",
	"C3" => "Proteins Digestion & Absorption",
	"C4" => "Others metabolism"
);

my %sample_func_summary;

my $ko_file = "$indir/$prefix.metapredict_gene.xls";
&LOAD_PRO ($ko_file, $prefix, "KO", \%sample_func_summary);

my $L2_file = "$indir/$prefix.pathway_profile.L2.xls";
&LOAD_PRO ($L2_file, $prefix, "L2", \%sample_func_summary);

my $L3_file = "$indir/$prefix.pathway_profile.L3.xls";
&LOAD_PRO ($L3_file, $prefix, "L3", \%sample_func_summary);

# sample's specific function stat
my $sample_func_json = encode_json \%sample_func_summary;
open OUT,">$out_dir/$prefix.function.json" || die $!;
print OUT "$sample_func_json";
close OUT;

###########################
sub LOAD_PRO {
	my ($in,$pre,$level,$info) = @_;

	open IN,"$in" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($name, $count, $per) = (split/\t+/,$_)[0,1,2];
		#$info -> {metabolism}{$level}{$name}{count} = $count;
		$info -> {metabolism}{$level}{$name}{$pre} = $per;
	}
	close IN;
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

