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
GetOptions(\%opts,"id=s","p=s","o=s","beni=s","harm=s","h");

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
# taxon stratified stat representation of sample's microbiota
my %level_key = (
	"Phylum" => "L2",
	"Class" => "L3",
	"Order" => "L4",
	"Family" => "L5",
	"Genus" => "L6",
	"Species" => "L7"
);

my %sample_summary;

$sample_summary{microbiota}{Genus}{Prevotella}{$prefix} = 0;
$sample_summary{microbiota}{Genus}{Bacteroides}{$prefix} = 0.000001;
$sample_summary{microbiota}{Phylum}{Bacteroidetes}{$prefix} = 0.000001;

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
		$sample_summary{microbiota}{$l}{$tax}{$prefix} = $per;
	}
	close IN;
}

# output stat
open OUT,">$out_dir/$prefix.summary.xls" || die $!;
open IN,"$indir/$prefix.stat.xls" || die $!;
while (<IN>) {
	chomp;
	print OUT "$_\n";
	s/^\s+|\s+$//;
	next if (/^$/ || /^\#/);
	my ($key, $v) = (split/\t+/,$_)[0,1];
	if (exists $level_key{$key}) {
		$sample_summary{summary}{representation}{$key}{$prefix} = $v;
		next;
	}
	$sample_summary{summary}{$key}{$prefix} = $v;
}
close IN;

my $max_genus;
$max_genus = (sort {$sample_summary{microbiota}{Genus}{$b}{$prefix}<=>$sample_summary{microbiota}{Genus}{$a}{$prefix}} keys %{$sample_summary{microbiota}{Genus}})[0];
$sample_summary{summary}{predominant_genus}{$prefix} = $max_genus;
$sample_summary{summary}{predominant_genus_per}{$prefix} = $sample_summary{microbiota}{Genus}{$max_genus}{$prefix};
printf OUT ("\n%s\n%s\n%s\t%.2f\n\n", "#microbiota stat", "Predominant Genus:", $max_genus, $sample_summary{microbiota}{Genus}{$max_genus}{$prefix}*100);

my $F_B_ratio = $sample_summary{microbiota}{Phylum}{Firmicutes}{$prefix}/$sample_summary{microbiota}{Phylum}{Bacteroidetes}{$prefix};
$sample_summary{summary}{"Firmicutes_Bacteroidetes_ratio"}{$prefix} = $F_B_ratio;
printf OUT ("%s%.2f\n\n", "Firmicutes/Bacteroidetes Ratio: ", $F_B_ratio);

my $P_B_ratio = $sample_summary{microbiota}{Genus}{Prevotella}{$prefix}/$sample_summary{microbiota}{Genus}{Bacteroides}{$prefix};
$sample_summary{summary}{"Prevotella_Bacteroides_ratio"}{$prefix} = $P_B_ratio;
printf OUT ("%s%.2f\n\n", "Prevotella/Bacteroides Ratio: ", $P_B_ratio);

#my $B_E_ratio = 0;
#$B_E_ratio = $sample_summary{microbiota}{Genus}{Bifidobacterium}{$prefix}/$sample_summary{microbiota}{Family}{Enterobacteriaceae}{$prefix} if (defined $sample_summary{microbiota}{Family}{Enterobacteriaceae}{$prefix} && defined $sample_summary{microbiota}{Genus}{Bifidobacterium}{$prefix});
#$sample_summary{summary}{"Bifidobacterium_Enterobacteriaceae_ratio"}{$prefix} = $B_E_ratio;
#printf OUT ("%s%.2f\n", "Bifidobacterium/Enterobacteriaceae Ratio: ", $B_E_ratio);

close OUT;

my $sample_summary_json = encode_json \%sample_summary;
open OUT, ">$out_dir/$prefix.microbiota.json" || die $!;
print OUT "$sample_summary_json";
close OUT;

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

