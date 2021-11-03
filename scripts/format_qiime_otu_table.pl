#!/usr/bin/perl -w
use strict;
use File::Basename qw(dirname);

unless (@ARGV == 3) {
	print "\n\tFunction: format otu table for qiime diversity analysis\n\n";
	print "\t  Usage: <uearch_otu_table> <uclust_assign_taxonomy_table> <format_out>\n";
	print "\t     usearch_otu_table:\n";
	print "\t        OTUId     Sample1    Sample2    ...\n";
	print "\t     uclust_assign_taxonomy_table:\n";
	print "\t        OTU_69  k__Bacteria; p__Bacteroidetes; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__ 1.00    3\n\n";
	exit;
}

my ($otu_tab,$tax_tab,$out) = @ARGV;
my %otu;

my $make_otu_biom = "make_otu_table.py";

open IN,"$otu_tab" || die $!;
my $head = <IN>; chomp $head;
$head =~ s/\s+$//;
my ($id,$sams) = split/\t+/,$head,2;
while (<IN>) {
	chomp;
	s/\s+$//;
	next if (/^$/);
	my ($otu_id,$count) = split /\t+/,$_,2;
	$otu{$otu_id} = $count;
}
close IN;

open IN,"$tax_tab" || die $!;
open OUT,">$out.xls" || die $!;
print OUT "#Qiime v1.9.1 OTU table\n";
print OUT "#OTU ID\t$sams\ttaxonomy\n";
while (<IN>) {
	chomp;
	s/\s+$//;
	next if (/^$/);
	my ($otu_id,$tax) = split /\t+/,$_;
	next if (!exists $otu{$otu_id});
	print OUT "$otu_id\t$otu{$otu_id}\t$tax\n";
}
close IN;
close OUT;

#+---------------
# convert tab to biom
#+---------------

my $odir = dirname($otu_tab);
system "docker run  --rm -v $odir:$odir qiime:1.9.1 biom convert -i $out.xls -o $out.biom --table-type \"OTU table\" --to-json --process-obs-metadata taxonomy\n";
#print "docker run -it --rm -v $odir:$odir qiime:1.9.1 biom convert -i $out.xls -o $out.biom --table-type \"OTU table\" --to-json --process-obs-metadata taxonomy";
system "docker run  --rm -v $odir:$odir qiime:1.9.1 biom summarize-table -i $out.biom -o $out.biom.summarize";
#print "docker run -it --rm -v $odir:$odir qiime:1.9.1 biom summarize-table -i $out.biom -o $out.biom.summarize\n";
