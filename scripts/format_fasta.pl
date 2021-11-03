#!/usr/bin/perl -w
use strict;

unless (@ARGV == 3) {
	print "\n\tFunction: format filtered fastq for qiime closed otu picking fasta\n\n";
	print "\t  Usage: <sample_lable> <fastq_in> <out_dir>\n\n";
	exit;
}

my ($sample,$in,$out) = @ARGV;
my %otu;

open IN,"$in" || die $!;
open OUT,">$out/$sample.filtered.fasta" || die $!;
my $i=1;
while (<IN>) {
	chomp;
	s/\s+$//;
	next if (/^$/);
	my $seq = <IN>;chomp $seq;
	$seq =~s/\s+//g;
	<IN>;
	<IN>;
	printf OUT ("%s%s%d%s%s", ">$sample", "_", $i, "\n", "$seq\n");
	$i++ ;
}
close IN;
close OUT;

