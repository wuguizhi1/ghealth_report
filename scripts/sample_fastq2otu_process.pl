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
GetOptions(\%opts,"l=s","left=s","right=s","o=s","c=s","join","m=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -l       :: sample lable (prefix of result files), required;
                  -left    :: read1 of sequence data, required;
                  -right   :: read2 of sequence data, optional;

                  -o    :: output dir default "./";
                  -join :: unmerged pair reads joined and used in next analysis;
                  -m    :: method for merge pair end reads (flash or usearch), default usearch;

              config-file:
                  -c    :: config of pipline parameters;

USAGE

if(!defined($opts{l}) || !defined($opts{left}))
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
my $conf=$opts{c} || "$Bin/../conf/fastq2otu.conf";
$conf=&ABSOLUTE_DIR($conf);
my $log_file="$out_dir/process.log";
system "touch $log_file";

my $merge_method = $opts{m} || "usearch";

my $chimera_udb = "$Bin/../database/rdp_gold.udb";
my $flash = "$Bin/../utils/flash";

#################
my %para;
&LOADPARA($conf, \%para);

#################
my $logfile = "$out_dir/$lable.process.log";
my $cmd;

## merge pair fastq
if (defined $right) {
	if ($merge_method eq "flash") {
		$cmd .= "$flash";
		$cmd .= " -m $para{flash_m} -M $para{flash_M} -x $para{flash_x}";
		$cmd .= " -d $out_dir -o $lable $left $right && mv $out_dir/$lable.extendedFrags.fastq $out_dir/$lable.merge.fastq 2>$out_dir/$lable.flash.log";
	}
	else {
		$cmd = "usearch8.1 -fastq_mergepairs $left -reverse $right -fastqout $out_dir/$lable.merge.fastq -tabbedout $out_dir/$lable.merge.stat -report $out_dir/$lable.merge.report ";
		$cmd .= "-fastqout_notmerged_fwd $out_dir/$lable.unmerged_1.fastq -fastqout_notmerged_rev $out_dir/$lable.unmerged_2.fastq ";
		$cmd .= "-fastq_maxdiffpct $para{fastq_maxdiffpct} -sample $lable -minhsp $para{minhsp} -fastq_minovlen $para{fastq_minovlen} 2>$logfile";
	}
	&LOG ("merge pair fastq reads..", $cmd, $log_file);
	system ($cmd);
}

# unmerged fastq join together
if (defined $opts{"join"}) {
	$cmd = "usearch8.1 -fastq_join $out_dir/$lable.unmerged_1.fastq -reverse $out_dir/$lable.unmerged_2.fastq -fastqout $out_dir/$lable.unmerged.join.raw.fastq ";
	$cmd .= "-join_padgap \"\" ";
	&LOG ("merge pair fastq reads..", $cmd, $log_file);
	system ($cmd);

	&relable_join_fq_ ("$out_dir/$lable.unmerged.join.raw.fastq", "$out_dir/$lable.unmerged.join.fastq", $lable);
	system ("rm $out_dir/$lable.unmerged.join.raw.fastq");

	system ("cat $out_dir/$lable.merge.fastq $out_dir/$lable.unmerged.join.fastq >$out_dir/$lable.combine.fastq");
}

# filt fastq
$cmd = "usearch8.1 -fastq_filter $left -fastqout $out_dir/$lable.filter.fastq " if (!defined $right);
$cmd = "usearch8.1 -fastq_filter $out_dir/$lable.merge.fastq -fastqout $out_dir/$lable.filter.fastq " if (defined $right && !defined $opts{"join"});
$cmd = "usearch8.1 -fastq_filter $out_dir/$lable.combine.fastq -fastqout $out_dir/$lable.filter.fastq " if (defined $opts{"join"});
$cmd .= "-fastq_truncqual $para{fastq_truncqual} -fastq_minlen $para{fastq_minlen} -fastq_maxns $para{fastq_maxns} -fastq_maxee_rate $para{fastq_maxee_rate} 2>>$logfile";
&LOG ("filter merged fastq reads..", $cmd, $log_file);
system ($cmd);

# derep fulllength
$cmd = "usearch8.1 -derep_fulllength $out_dir/$lable.filter.fastq -fastaout $out_dir/$lable.derep.fasta ";
$cmd .= "-sizeout 2>>$logfile";
&LOG ("derep fulllength of filted fastq reads..", $cmd, $log_file);
system ($cmd);

# cluster_otus
$cmd = "usearch8.1 -cluster_otus $out_dir/$lable.derep.fasta -otus $out_dir/$lable.otus.fasta -uparseout $out_dir/$lable.otus.up -relabel OTU ";
$cmd .= "-minsize $para{minsize} -otu_radius_pct $para{otu_radius_pct} -sizeout 2>>$logfile";
&LOG ("cluster otus of $lable..", $cmd, $log_file);
system ($cmd);

# uchime ref
$cmd = "usearch8.1 -uchime_ref $out_dir/$lable.otus.fasta -db $chimera_udb -uchimeout $out_dir/$lable.otus.uchime -chimeras $out_dir/$lable.otus.chimera.fasta -nonchimeras $out_dir/$lable.otus.nonchimera.fasta ";
$cmd .= "-strand plus 2>>$logfile";
&LOG ("otus exclude chimeras..", $cmd, $log_file);
system ($cmd);

# rarefaction and alpha diversity
#$cmd = "usearch8.1 -fasta_diversity $out_dir/$lable.otus.nonchimera.fasta -output $out_dir/$lable.diversity.txt ";
#$cmd .= "-iters 100 2>>$logfile";
#&LOG ("caculated diversity of $lable..", $cmd, $log_file);
#system ($cmd);

system "sed 's/;size.*//g' $out_dir/$lable.otus.nonchimera.fasta > $out_dir/$lable.otus.nonchimera.fasta.bac";
system "mv $out_dir/$lable.otus.nonchimera.fasta.bac $out_dir/$lable.otus.nonchimera.fasta";

# global search
$cmd = "usearch8.1 -usearch_global $out_dir/$lable.filter.fastq -db $out_dir/$lable.otus.nonchimera.fasta -qout $out_dir/$lable.qout -otutabout $out_dir/$lable.otus.tab.xls ";
$cmd .= "-strand plus -id $para{globalid} 2>>$logfile";
&LOG ("map filted fastq reads to otus..", $cmd, $log_file);
system ($cmd);

# rarefaction and alpha diversity
$cmd = "usearch10.0 -alpha_div $out_dir/$lable.otus.tab.xls -output $out_dir/$lable.diversity.txt ";
$cmd .= "2>>$logfile";
&LOG ("caculated diversity of $lable..", $cmd, $log_file);
system ($cmd);

# qiime assign taxon
$cmd = "docker run --rm -v $out_dir/:$out_dir/ qiime:1.9.1 assign_taxonomy.py -i $out_dir/$lable.otus.nonchimera.fasta -o $out_dir/assigned_taxonomy ";
$cmd .= "-m $para{assign_m}";
&LOG ("assign otus taxon by qiime..", $cmd, $log_file);
system ($cmd);

if ($para{"assign_m"} eq "rdp") {
	system "sed s'/;/; /g' $out_dir/assigned_taxonomy/$lable.otus.nonchimera_tax_assignments.txt >$out_dir/assigned_taxonomy/$lable.otus.nonchimera_tax_assignments.txt.bac";
	system "mv $out_dir/assigned_taxonomy/$lable.otus.nonchimera_tax_assignments.txt.bac $out_dir/assigned_taxonomy/$lable.otus.nonchimera_tax_assignments.txt";
}

$cmd = "perl $Bin/format_qiime_otu_table.pl $out_dir/$lable.otus.tab.xls $out_dir/assigned_taxonomy/$lable.otus.nonchimera_tax_assignments.txt $out_dir/$lable.assign_taxon";
&LOG ("format qiime otu table..", $cmd, $log_file);
system ($cmd);

# summurize taxon and format result
$cmd = "docker run --rm -v $out_dir/:$out_dir/ qiime:1.9.1 sort_otu_table.py -i $out_dir/$lable.assign_taxon.biom -o $out_dir/$lable.assign_taxon.sort.biom ";
$cmd .= " && docker run --rm -v $out_dir/:$out_dir/ qiime:1.9.1 summarize_taxa.py -i $out_dir/$lable.assign_taxon.sort.biom -o $out_dir -L 2,3,4,5,6,7";
print STDOUT "OUTOUT:$cmd\n";
&LOG ("summarize and stats taxon by qiime..", $cmd, $log_file);
system ($cmd);

# stat format result
$cmd = "perl $Bin/micro_result_format.pl -t $out_dir/assigned_taxonomy/$lable.otus.nonchimera_tax_assignments.txt -i $out_dir ";
$cmd .= "-p $lable -o $out_dir/result ";
$cmd .= "&& perl $Bin/micro_summary.pl -id $out_dir/result -p $lable -o $out_dir/result";
&LOG ("format output result..", $cmd, $log_file);
system ($cmd);

###########################
sub relable_join_fq_ {
	my ($in,$out,$l) = @_;

	open IN,"$in" || die $!;
	open OUT,">$out" || die $!;
	while (<IN>) {
		chomp;
		print OUT "$_;sample=$l;\n";
		my $seq = <IN>; chomp $seq;
		my $qual_id = <IN>; chomp $qual_id;
		my $qual = <IN>; chomp $qual;
		print OUT "$seq\n$qual_id\n$qual\n";
	}

	close IN;
	close OUT;
}

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

