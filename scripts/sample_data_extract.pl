#!/usr/bin/perl -w
#use strict;
use warnings;
use PerlIO::gzip;
use Data::Dumper;

unless (@ARGV == 7) {
	print "\n\tFunction: extract data from multisamples to one index of 16sRNA rawdata\n\n";
	print "\t\tUsage: <Sam_ids> <barcode.cfg> <barcodes info> <fq1> <fq2> <out_dir> <barcode_log>\n\n";
	exit;
}

my %cutlen = (
	"11" => "3",
	"12" => "0",
	"21" => "4",
	"22" => "2",
	"31" => "5",
	"32" => "3",
	"41" => "6",
	"42" => "6",
);

my ($sam_str,$barcode_cfg,$barcode_info,$fq1,$fq2,$od,$stat) = @ARGV;

my %barcodepairs;
my %splitstat;
&load_barcode_pair_ ($barcode_cfg, \%barcodepairs, \%splitstat);

my %saminfo;
my @samples = split/;/,$sam_str;
&load_sample_info_ ($sam_str, $barcode_info, \%saminfo);

&get_sample_data_ (\@samples, \%saminfo, \%barcodepairs, \%cutlen, $fq1, $fq2, $od, \%splitstat);

&out_log_ ($sam_str,\%splitstat,$stat);

##################
# subs
sub load_barcode_pair_ {
	my ($file,$pair,$splitstat) = @_;

	open IN,"$file" || die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		$_ =~ /F(\d+)R(\d+)\s+([ATCG]+);([ATCG]+)/;
		$pair->{$3}{$4} = $1;
		$splitstat->{"F$1"."R$1"} = 0;
	}
	close IN;
}

sub load_sample_info_ {
	my ($sam_str, $barcode_str, $sample_info) = @_;

	my @sams = split/;/,$sam_str;
	my @barcodes = split/;/,$barcode_str;

	for (0..$#sams) {
		my $num = $_;
		my @barcodeinfo = split/,/,$barcodes[$num];
		foreach my $barcode_pair (@barcodeinfo) {
			$sample_info->{$barcode_pair} = $sams[$num];
		}
	}
}

sub get_sample_data_ {
	my ($sams, $saminfo, $barcodepairs, $cut_len, $fq1, $fq2, $od, $splitstat) = @_;

	foreach my $sam (@{$sams}) {

		my $r1 = "$sam"."1";
		my $r2 = "$sam"."2";
		open $r1,">$od/$sam"."_R1.fastq" || die $!;
		open $r2,">$od/$sam"."_R2.fastq" || die $!;
	}

	open FQ1,"gzip -dc $fq1|" if ($fq1 =~ /\.gz$/);
	open FQ1,"$fq1" if ($fq1 !~ /\.gz$/);
	open FQ2,"gzip -dc $fq2|" if ($fq2 =~ /\.gz$/);
	open FQ2,"<$fq2" if ($fq2 !~ /\.gz$/);

	while (<FQ1>) {
		chomp;
		my $fq1_id = $_;
		my $fq1_seq = <FQ1>; chomp $fq1_seq;
		<FQ1>;
		my $fq1_qual = <FQ1>; chomp $fq1_qual;
		my $F_barcode = substr ($fq1_seq,0,6);

		my $fq2_id = <FQ2>; chomp $fq2_id;
		my $fq2_seq = <FQ2>; chomp $fq2_seq;
		<FQ2>;
		my $fq2_qual = <FQ2>; chomp $fq2_qual;
		my $R_barcode = substr ($fq2_seq,0,6);

		if (exists $barcodepairs->{$F_barcode}{$R_barcode} && exists $saminfo->{$barcodepairs->{$F_barcode}{$R_barcode}}) {
			my $sample = $saminfo->{$barcodepairs->{$F_barcode}{$R_barcode}};
			$splitstat->{"F$barcodepairs->{$F_barcode}{$R_barcode}"."R$barcodepairs->{$F_barcode}{$R_barcode}"} ++;
			my $r1_cut = "$barcodepairs->{$F_barcode}{$R_barcode}"."1";
			my $r2_cut = "$barcodepairs->{$F_barcode}{$R_barcode}"."2";
			my $r1 = "$sample"."1";
			my $r2 = "$sample"."2";
			my $new_fq1 = substr ($fq1_seq,$cut_len->{$r1_cut});
			my $new_fq2 = substr ($fq2_seq,$cut_len->{$r2_cut});
			my $new_qual1 = substr ($fq1_qual,$cut_len->{$r1_cut});
			my $new_qual2 = substr ($fq2_qual,$cut_len->{$r2_cut});
			print $r1 "$fq1_id\n$new_fq1\n+\n$new_qual1\n";
			print $r2 "$fq2_id\n$new_fq2\n+\n$new_qual2\n";
		}
	}

	close FQ1;
	close FQ2;

	foreach my $sam (@{$sams}) {

		my $r1 = "$sam"."1";
		my $r2 = "$sam"."2";
		close $r1;
		close $r2;
	}
}

#&out_log_ ($sam_str,\%splitstat,$stat);
sub out_log_ {
	my ($sam_str, $splitstat, $stat) = @_;


	my $head = 0;
	$head = 1 if (-f $stat);

	open OUT,">>$stat" || die $!;
	my $head_str = "";
	my $val_str = "";
	foreach my $pair (sort keys %{$splitstat}) {
		$head_str .= "$pair\t";
		$val_str .= "$splitstat->{$pair}\t";
	}
	$head_str=~s/\s+$//;
	$val_str=~s/\s+$//;

	print OUT "#Sample\t$head_str\n" if ($head == 0);
	print OUT "$sam_str\t$val_str\n";
	close OUT;
}
