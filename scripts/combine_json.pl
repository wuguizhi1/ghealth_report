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
use JSON::XS;
use Hash::Merge qw(merge);

my $path=`pwd`;
chomp($path);

my %opts;
GetOptions(\%opts,"micro=s","meta=s","prefix=s","o=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -micro    :: microbiota json result, required;
                  -meta     :: metabolism json result, required;
                  -prefix   :: output files prefix, required;

                  -o     :: output dir default "./";

USAGE

if(!defined($opts{micro}) || !defined($opts{prefix}) || !defined($opts{meta}) || defined($opts{h}))
{
	print $usage;
	exit;
}

my $micro_json=&ABSOLUTE_DIR($opts{micro});
my $meta_json=&ABSOLUTE_DIR($opts{meta});

my $prefix=$opts{prefix};
my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

#################
# load microbiota and metabolism json file
open IN,"$micro_json" || die $!;
my $micro_str;
while (<IN>) {
	chomp;
	s/^\s+|\s+$//;
	next if (/^$/);
	$micro_str .= $_;
}
close IN;

my $microbiota = decode_json $micro_str;

open IN,"$meta_json" || die $!;
my $meta_str;
while (<IN>) {
	chomp;
	s/^\s+|\s+$//;
	next if (/^$/);
	$meta_str .= $_;
}
close IN;

my $metabolism = decode_json $meta_str;

# combine sample's micro and meta result
my %sample_result = %{ merge (\%$microbiota, \%$metabolism)};

my $result = encode_json \%sample_result;
open OUT,">$out_dir/$prefix.result.json" || die $!;
print OUT "$result";
close OUT;

#print Dumper %sample_result;

###########################
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

