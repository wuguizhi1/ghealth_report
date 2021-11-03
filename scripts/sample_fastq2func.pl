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
GetOptions(\%opts,"l=s","fq=s","o=s","c=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -l       :: sample lable (prefix of result files), required;
                  -fq      :: sample filtered fastq sequence file, required;

                  -o     :: output dir default "./";

              config-file:
                  -c    :: config of pipline parameters;

USAGE

if(!defined($opts{l}) || !defined($opts{fq}))
{
	print $usage;
	exit;
}

my $lable=$opts{l};

my $fq=&ABSOLUTE_DIR($opts{fq});
my $out_dir=$opts{o} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);
my $conf=$opts{c} || "$Bin/../conf/fastq2func.conf";
$conf=&ABSOLUTE_DIR($conf);

my $log_file="$out_dir/process.log";
system "touch $log_file";

my $gg_ref = "$Bin/../database/gg_13_5_otus/rep_set/97_otus.fasta";
my $gg_tax = "$Bin/../database/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt";

#################
my %para;
&LOADPARA($conf, \%para);

#################
my $logfile = "$out_dir/$lable.process.log";
my $cmd;

# fastq format to fasta
$cmd = "perl $Bin/format_fasta.pl $lable $fq $out_dir";
system ($cmd);  

# picking closed reference otus
system "echo pick_otus:enable_rev_strand_match $para{enable_rev_strand_match} >> $out_dir/otu_picking_params.txt";
system "echo pick_otus:similarity $para{similarity} >> $out_dir/otu_picking_params.txt";

my $ref_dir = dirname($gg_ref);
my $tax_dir = dirname($gg_tax);
#$cmd = "pick_closed_reference_otus.py -f -i $out_dir/$lable.filtered.fasta -p $out_dir/otu_picking_params.txt -t $gg_tax ";
#$cmd .= "-r $gg_ref -o $out_dir 2>$logfile";
$cmd = "docker run --name $lable.00 -v $out_dir:/outdir -v $ref_dir:/ref_dir -v $tax_dir:/tax_dir qiime:1.9.1 ";
$cmd .= "pick_closed_reference_otus.py -f -i /outdir/$lable.filtered.fasta -p /outdir/otu_picking_params.txt ";
$cmd .= "-t /tax_dir/97_otu_taxonomy.txt -r /ref_dir/97_otus.fasta -o /outdir";
&LOG ("pick closed reference otus..", $cmd, $log_file);
system ($cmd);

# normalize otu table
$cmd = "docker run --name $lable.01 -v $out_dir:/pircust picrust:1.0.0 normalize_by_copy_number.py -i pircust/otu_table.biom -o pircust/$lable.normalized_otus.biom 2>>$logfile";
&LOG ("derep fulllength of filted fastq reads..", $cmd, $log_file);
system ($cmd);

# Predict Functions For Metagenome
$cmd = "docker run --name $lable.02 -v $out_dir:/pircust picrust:1.0.0 predict_metagenomes.py -i pircust/$lable.normalized_otus.biom -o pircust/$lable.metagenome_predictions.biom -a pircust/$lable.nsti_per_sample.tab ";
$cmd .= "--type_of_prediction $para{type_of_prediction} 2>>$logfile";
&LOG ("predict functions for metagenome of $lable..", $cmd, $log_file);
system ($cmd);

$cmd = "docker run --rm -v $out_dir/:$out_dir/ qiime:1.9.1 biom convert -i $out_dir/$lable.metagenome_predictions.biom -o $out_dir/$lable.metagenome_predictions.txt ";
$cmd .= " --to-tsv --header-key KEGG_Description";
system ($cmd);

# categorize function stat
$cmd = "docker run --name $lable.03 -v $out_dir:/pircust picrust:1.0.0 categorize_by_function.py -f -i pircust/$lable.metagenome_predictions.biom -c KEGG_Pathways -l 3 -o pircust/$lable.metagenome_predictions.L3.txt ";
$cmd .= "&& docker run --name $lable.04 -v $out_dir:/pircust picrust:1.0.0 categorize_by_function.py -f -i pircust/$lable.metagenome_predictions.biom -c KEGG_Pathways -l 2 -o pircust/$lable.metagenome_predictions.L2.txt ";
$cmd .= "&& docker run --name $lable.05 -v $out_dir:/pircust picrust:1.0.0 categorize_by_function.py -f -i pircust/$lable.metagenome_predictions.biom -c KEGG_Pathways -l 1 -o pircust/$lable.metagenome_predictions.L1.txt 2>>$log_file";
&LOG ("categorize kegg function stat..", $cmd, $log_file);
system ($cmd);

# partitions function otus stat
$cmd = "docker run --name $lable.06 -v $out_dir:/pircust picrust:1.0.0 metagenome_contributions.py -i pircust/$lable.normalized_otus.biom -t ko -o pircust/$lable.ko_metagenome_contributions.tab ";
$cmd .= " 2>>$logfile";
&LOG ("partitions function otus stat of $lable..", $cmd, $log_file);
system ($cmd);

# rm docker container
$cmd = "docker ps -a | grep $lable |awk '{print \$1}' | xargs --no-run-if-empty docker rm ";
system ($cmd);

# stat format result
$cmd = "perl $Bin/picrust_result_format.pl -t $out_dir/$lable.metagenome_predictions.txt -i $out_dir -sum $lable.metagenome_predictions ";
$cmd .= "-p $lable -o $out_dir/result ";
$cmd .= "&& perl $Bin/func_summary.pl -id $out_dir/result -p $lable -o $out_dir/result";
&LOG ("format output result..", $cmd, $log_file);
system ($cmd);

###########################
sub LOG {
	my ($info,$commond,$log) = @_;

	open OUT,">>$log" || die $!;
	print OUT "#$info\n$commond\n\n";
	close OUT;
}

sub LOADPARA {
	my ($cfg, $para) = @_;
	open IN,"$cfg" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($p, $v) =(split /\s+/,$_)[0,1];
		$para -> {$p} = $v;
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

