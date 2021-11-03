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
GetOptions(\%opts,"t=s","i=s","sum=s","p=s","o=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -t     :: metagenome prediction gene table (picrust), required;
                  -i     :: prediction categorize ko results directory, required;
                  -sum   :: prediction categorize ko sum prefix (picrust), required;
                  -p     :: outfiles prefix, required;

                  -o     :: output dir default "./";

USAGE

if(!defined($opts{t}) || !defined($opts{i}) || !defined($opts{sum}) || !defined($opts{p}))
{
	print $usage;
	exit;
}

my $tab=&ABSOLUTE_DIR($opts{t});
my $indir=&ABSOLUTE_DIR($opts{i});
my $sum=$opts{sum};
my $prefix=$opts{p};

my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

#################
# stat otu classfified
my %gene_orth;
my %orth_desc;
my $total;
open IN,"$tab" || die $!;
while (<IN>) {
	chomp;
	s/^\s+|\s+$//;
	next if (/^$/ || /^\#/);
	my ($k, $num, $desc) = (split/\t+/,$_)[0,1,2];
	$total += $num;
	$gene_orth{$k} = $num;
	$orth_desc{$k} = $desc;
}
close IN;

open OUT, ">$out_dir/$prefix.metapredict_gene.xls" || die $!;
foreach my $k (sort {$gene_orth{$b} <=> $gene_orth{$a}} keys %gene_orth) {
	if ($gene_orth{$k} ne "0.0") {
		my $per = $gene_orth{$k}/$total;
		printf OUT ("%s\t%d\t%e\t%s\n", $k, $gene_orth{$k}, $per, $orth_desc{$k});
	}
	else {
		print OUT "$k\t0\t0\t$orth_desc{$k}\n";
	}
}

close OUT;

# taxon stratified stat representation of sample's microbiota
my @levels =sort glob("$indir/$sum*L*txt");

foreach my $level_file (@levels) {
	$level_file =~ /$sum.(L\d)\.txt/;
	my $l = $1;
	my $outfile = "$out_dir/$prefix.pathway_profile.$l.xls" || die $!;
	&FORMAT($level_file, $outfile);
}

###########################
sub FORMAT {
	my ($in, $out) = @_;

	my $total_num;
	my %level_info;
	my %level_desc;

	open IN,"$in" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($k, $num, $desc) = (split/\t+/,$_)[0,1,2];
		$total_num += $num;
		$level_info{$k} = $num;
		$level_desc{$k} = $desc;
	}
	close IN;

	open OUT, ">$out" || die $!;
	print OUT "#Path\tCount\tPercent\tDescription\n";
	foreach my $k (sort {$level_info{$b} <=> $level_info{$a}} keys %level_info) {
		my $per = 0;
		if ($level_info{$k} ne "0.0") {
			$per = $level_info{$k}/$total_num;
			printf OUT ("%s\t%d\t%e\t%s\n", $k, $level_info{$k}, $per, $level_desc{$k});
		}
		else {
			print OUT "$k\t0\t0\t$level_desc{$k}\n";
		}
	}

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

