#!/usr/bin/perl -w
#use strict;
use warnings;
use PerlIO::gzip;
use Data::Dumper;

unless (@ARGV == 6) {
	print "\n\tFunction: cut second barcode pair seq of 16sRNA rawdata\n\n";
	print "\t\tUsage: <Sam_ids> <barcode.cfg> <fq1> <fq2> <out_dir> <barcode_log>\n\n";
	exit;
}

my %cutlen = (
	"F1" => "3",
	"R1" => "0",
	"F2" => "4",
	"R2" => "2",
	"F3" => "5",
	"R3" => "3",
	"F4" => "6",
	"R4" => "6",
);

my ($sam_str,$barcode_cfg,$fq1,$fq2,$od,$stat) = @ARGV;

my %barcode;
my %barcodestat;
&load_barcode_pair_ ($barcode_cfg, \%barcode, \%barcodestat);

&get_sample_data_ ($sam_str, \%barcode, \%cutlen, $fq1, $fq2, $od, \%barcodestat);

&out_log_ ($sam_str,\%barcodestat,$stat);

##################
# subs
sub load_barcode_pair_ {
	my ($file,$barcode,$barcodestat) = @_;

	open IN,"$file" || die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		$_ =~ /(F\d+)(R\d+)\s+([ATCG]+);([ATCG]+)/;
		$barcode->{F}{$3} = $1;
		$barcode->{R}{$4} = $2;
		$barcodestat->{$1} = 0;
		$barcodestat->{$2} = 0;
	}
	close IN;
}

sub get_sample_data_ {
	my ($sam_str, $barcode, $cut_len, $fq1, $fq2, $od, $barcodestat) = @_;

	open R1,">$od/$sam_str"."_R1.fastq" || die $!;
	open R2,">$od/$sam_str"."_R2.fastq" || die $!;

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

		if (exists $barcode->{F}{$F_barcode} && exists $barcode->{R}{$R_barcode}) {
			$barcodestat->{$barcode->{F}{$F_barcode}} ++;
			$barcodestat->{$barcode->{R}{$R_barcode}} ++;
			my $new_fq1 = substr ($fq1_seq,$cut_len->{$barcode->{F}{$F_barcode}});
			my $new_fq2 = substr ($fq2_seq,$cut_len->{$barcode->{R}{$R_barcode}});
			my $new_qual1 = substr ($fq1_qual,$cut_len->{$barcode->{F}{$F_barcode}});
			my $new_qual2 = substr ($fq2_qual,$cut_len->{$barcode->{R}{$R_barcode}});
			print R1 "$fq1_id\n$new_fq1\n+\n$new_qual1\n";
			print R2 "$fq2_id\n$new_fq2\n+\n$new_qual2\n";
		}
	}

	close FQ1;
	close FQ2;
	close R1;
	close R2;

}

#&out_log_ ($sam_str,\%barcodestat,$stat);
sub out_log_ {
	my ($sam_str, $barcodestat, $stat) = @_;


	my $head = 0;
	$head = 1 if (-f $stat);

	open OUT,">>$stat" || die $!;
	my $head_str = "";
	my $val_str = "";
	foreach my $barcode_str (sort keys %{$barcodestat}) {
		$head_str .= "$barcode_str\t";
		$val_str .= "$barcodestat->{$barcode_str}\t";
	}
	$head_str=~s/\s+$//;
	$val_str=~s/\s+$//;

	print OUT "#Sample\t$head_str\n" if ($head == 0);
	print OUT "$sam_str\t$val_str\n";
	close OUT;
}
