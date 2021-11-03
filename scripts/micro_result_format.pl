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
GetOptions(\%opts,"t=s","i=s","p=s","o=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -t     :: otu assign taxon list (qiime), required;
                  -i     :: assign taxon sum results directory, required;
                  -p     :: input files and outfiles prefix, required;

                  -o     :: output dir default "./";
                  -h     help info;

USAGE

if(!defined($opts{t}) || !defined($opts{i}) || !defined($opts{p}))
{
	print $usage;
	exit;
}

my $tax=&ABSOLUTE_DIR($opts{t});
my $indir=&ABSOLUTE_DIR($opts{i});
my $prefix=$opts{p};

my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

# stat otu classfified
my %otu_stat;
$otu_stat{Unclassified} = 0;
open IN,"$tax" || die $!;
while (<IN>) {
	chomp;
	s/^\s+|\s+$//;
	next if (/^$/);
	my ($otu, $taxon) = (split/\t+/,$_)[0,1];
	$otu_stat{Total}++ ;
	if ($taxon !~ /k__/) {
		$otu_stat{Unclassified}++ ;
	}
	if ($taxon =~ /k__/) {
		$otu_stat{Classified}++ ;
	}
}
close IN;

# taxon stratified stat representation of sample's microbiota
my %level_key = (
	"L2" => "Phylum",
	"L3" => "Class",
	"L4" => "Order",
	"L5" => "Family",
	"L6" => "Genus",
	"L7" => "Species"
);

my @levels =sort glob("$indir/$prefix.assign_taxon.sort_*txt");
my %level_rep;
my %level_total;
my %unclassified;

foreach my $level_file (@levels) {
	$level_file =~ /$prefix.assign_taxon.sort_(L\d)\.txt/;
	my $l = $1;
	open IN,"$level_file" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($tax, $per) = (split/\t+/,$_)[0,1];
		if ($tax =~ /__$|Other$/) {
			$unclassified{$l} += $per;
			next;
		}
		$tax =~s/\[|\]//g;
		my @taxs = split/;/, $tax;
		if ($l ne "L7") {
			$taxs[-1] =~ /__(.*)$/;
			$level_rep{$l}{$1} = $per;
			$level_total{$l} += $per;
		}
		else {
			$taxs[-2] =~ /g__(.*)/;
			my $name = "$1";
			$taxs[-1] =~ /s__(.*)/;
			$name .= " $1";
			$level_rep{$l}{$name} = $per;
			$level_total{$l} += $per;
		}
	}
	close IN;
}

# get shannon diversity index
my $div_file = "$indir/$prefix.diversity.txt";
#my $shannon_div = `tail -1 $div_file |cut -f 4`; chomp $shannon_div;
my $shannon_div = `tail -1 $div_file |cut -f 13`; chomp $shannon_div;

# output stat
open OUT,">$out_dir/$prefix.stat.xls" || die $!;
print OUT "#OTU stat\n";
print OUT "Total_OTU\t$otu_stat{Total}\nClassified_OTU\t$otu_stat{Classified}\nUnclassified_OTU\t$otu_stat{Unclassified}\n\n";

print OUT "#Shannon diversity\n";
print OUT "diversity\t$shannon_div\n\n";

print OUT "#Taxon represention\n";
print OUT "#taxon\tclassified(%)\tunclassified(%)\n";

foreach my $l (sort keys %level_key) {
	$unclassified{$l} = 0 if (!defined $unclassified{$l});
	printf OUT ("%s\t%.2f\t%.2f\n", $level_key{$l}, $level_total{$l}*100, $unclassified{$l}*100);
	my $filehand = uc ($l);
	open $filehand, ">$out_dir/$prefix.$level_key{$l}.profile.xls" || die $!;
	foreach my $name (sort {$level_rep{$l}{$b}<=>$level_rep{$l}{$a}} keys %{$level_rep{$l}}) {
		my $standard_per = $level_rep{$l}{$name}/$level_total{$l};
		printf $filehand ("%s\t%e\t%e\n", $name, $level_rep{$l}{$name}, $standard_per);
	}
	printf $filehand ("%s\t%e\t%s\n", "Unclassified", $unclassified{$l}, "--");
	close $filehand;
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

