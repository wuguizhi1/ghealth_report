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

my $path=`pwd`;
chomp($path);

my %opts;
GetOptions(\%opts,"id=s","p=s","o=s","t=s","conf=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -id    :: formated micro results directory, required;
                  -p     :: input files prefix, required;
                  -conf  :: bacteria conf file type, stat sample's bac percent, required;
                          eg. "/data/bioit/biodata/mengf/Project/16srDNA/release/alpha/conf/good.bacteria";

                  -t     :: output files prefix, required;
                  -o     :: output dir default "./";

                  -h     help info;

USAGE

if(!defined($opts{id}) || !defined($opts{p}) || !defined($opts{conf}) || !defined($opts{t}))
{
	print $usage;
	exit;
}

my $indir=&ABSOLUTE_DIR($opts{id});
my $prefix=$opts{p};
my $conf=&ABSOLUTE_DIR($opts{conf});
my $out=$opts{t};

my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

#################
# load benificial and harmful bacteria config
my %bac_conf;
&LOAD_BAC ($conf, \%bac_conf);

# taxon stratified stat representation of sample's microbiota
my %class_info;

my @levels =sort glob("$indir/$prefix*.profile.xls");

foreach my $level_file (@levels) {
	$level_file =~ /$prefix\.(.*)\.profile\.xls/;
	my $l = $1;
	open IN,"$level_file" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($tax, $per) = (split/\t+/,$_)[0,1];
		next if ($tax =~ /Unclassified/);
		if (exists $bac_conf{$l}{$tax}) {
			$class_info{"total percent"} += $per;
			$class_info{composition}{$l}{$tax} = $per;
		}
	}
	close IN;
}

# output stat

printf ("%s\t%.2f\n","total percent(%)", $class_info{"total percent"}*100);

open OUT,">$out_dir/$out.info.xls" || die $!;
foreach my $level (sort keys %{$class_info{composition}}) {
	print OUT "#$level\n";
	foreach my $tax (sort keys %{$class_info{composition}{$level}}) {
		printf OUT ("%s\t%e\n", $tax, $class_info{composition}{$level}{$tax});
	}
	print OUT "\n";
}
close OUT;

###########################
sub LOAD_BAC {
	my ($in,$bac_l) = @_;

	open IN,"$in" || die $!;
	my $i = 1;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($level, $name) = split/\t+/,$_;
		$bac_l -> {$level}{$name} = $i;

		$i++ ;
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

