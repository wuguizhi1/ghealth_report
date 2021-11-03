#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use Data::Dumper;
use Encode;
use utf8;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my $demotime = &GETtime;

my($samID,$data,$out,$conf,$csvdir,$testdate);
GetOptions(
	"samID=s"	=>\$samID,
	"data=s"	=>\$data,
	"out=s"		=>\$out,
	"conf=s"	=>\$conf,
	"csvdir=s"	=>\$csvdir,
	"testdate=s"	=>\$testdate,
);

my $usage =<<USAGE;
Description:
Usage:  perl $0 [options]
	--samID		samples' list file
	--data		sample fastq data
	--out		output dir
			defaut "./";
	--conf		report conf dir
			default "$Bin/conf/report-version.new.conf"
	--csvdir	csv dir
			default "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/SampleInfo"

	--testdate	ex:$demotime
USAGE
die $usage if (!$samID || !$data);
#================
my %order = (
	"1" => "barcode",
	"2" => "source",
	"3" => "version",
	"4" => "素材分组",
);

$out||= "./";
$conf||= "$Bin/conf/report-version.new.conf";
$csvdir||= "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/SampleInfo";
&MKDIR($out);
$testdate ||= $demotime;
$out=&Absolutepath("dir",$out);
$data=&Absolutepath("dir",$data);
$samID=&Absolutepath("file",$samID);
#&Absolutepath("dir",$csvdir);
#&Absolutepath("dir",$conf);
#print "$data\t$out\n";
#my $script = "/data/bioit/biodata/mengf/pipline/16sRNA/amplicon-pip";

my $work_sh = "$out/work_sh";
&MKDIR($work_sh);
$work_sh=&Absolutepath("dir",$work_sh);

my %sample;
&load_sample_andCheck($samID,\%sample,$csvdir);
#print Dumper \%sample;

my %config;
&load_config($conf,\%config);
#print Dumper \%config;
if(!-f "$work_sh/analysis.OK"){
	my $analysis_sh = "$work_sh/analysis.sh";
	&generate_sample_result($data,\%sample,$analysis_sh,$out);

	my $cmd1 = "parallel -j 8 < $analysis_sh >$work_sh/analysis.log";
#	print "parallel -j 8 < $analysis_sh >$work_sh/analysis.log...";
	print "$cmd1...\n";
	system "$cmd1";
	system "touch $work_sh/analysis.OK";
}
if(!-f "$work_sh/cp_result.OK"){
	print "cp csv files and anlysis results...\n";
	&prepair_result($out,\%sample,$csvdir);
	system "touch $work_sh/cp_result.OK";
}
if(!-f "$work_sh/main.OK"){
	my $main_sh = "$work_sh/main.sh";
	my $readme = "$out/Readme";
	print "generating final report ...\n";
	&generate_report_sh(\%sample,$out,\%config,$main_sh,$readme);
	system "parallel -j 8 < $main_sh >$work_sh/main.log";
	system "touch $work_sh/main.OK";
}

#=================sub
sub load_sample_andCheck{
	my ($in,$sample,$csvdir) = @_;
	open IN,$in or die $!;
	while(<IN>){
		chomp;
		next if $_ =~ /#/;
		$_ =~ s/\s+//g;
		$sample->{$_} = 1;
	}
	for my $samID (keys %{$sample}){
		if(!-f "$csvdir/$samID.csv"){
			print STDERR"$csvdir/$samID.csv does not exist!\n";
			delete($sample{$samID});
		}
	}
	return $sample;
}

sub load_config{
	my ($input,$config) = @_;
	open IN,$input or die $!;
	while(<IN>){
		next if ($_ =~ /^#/);
		chomp;
		my @unit = split /\t/;
		$config->{$unit[0]}{$unit[1]} = "$unit[2]\t$unit[3]";
	}
	close IN;
}

sub generate_sample_result{
	my ($data,$sample,$sh,$out) = @_;
	my $OUT = "$out/analysis";
	&MKDIR ($OUT);
#	my $cmd_sh;
	open OUT,">$sh";
	for my $sam (keys %{$sample}){
		print OUT"perl $Bin/../bin/sample_taxfunc.pl -l $sam -left $data/$sam\_R1.fastq -right $data/$sam\_R2.fastq -o $OUT/$sam\n";
	}

#	system"parallel -j 5 < $sh";
}

sub prepair_result{
	my ($out,$sample,$csvdir) = @_;
	my $OUT = "$out/report";
	&MKDIR($OUT);
	for my $samID (keys %{$sample}){
		my $dir = "$OUT/$samID";
		&MKDIR($dir);
		system "cp $out/analysis/$samID/report/* $dir";
		system "cp $csvdir/$samID.csv $dir";
	}
}

sub generate_report_sh{
	my ($sample,$out,$config,$sh,$readme) = @_;
	open OUT,">$sh" or die $!;
	open OUT1,">$readme" or die $!;
	print OUT1"barcode\tsource\tversion\tsucaifenzu\n";
	for my $samID (keys %$sample){
		open IN,"$out/report/$samID/$samID.csv" or die $!;
		chomp(my $line = <IN>);
		my @head = split/\t/,$line;
		my @info;
		my @readmeinfo;
		while(<IN>){
			chomp;
			my @unit = split/\t/,$_;
			for(my $i=0;$i<@head;$i++){
				$head[$i] = decode("utf-8",$head[$i]);
				
				if($head[$i] eq 'version' || $head[$i] eq '素材分组'){
					push @info,$unit[$i];
				}
				if($head[$i] eq 'barcode' || $head[$i] eq 'source' || $head[$i] eq 'version' || $head[$i] eq '素材分组'){
					push @readmeinfo,$unit[$i];
				}
			}
		}
		print OUT1 join"\t",@readmeinfo;
		print OUT1"\n";
		print "@info\n";
#		my $stat = @info;
#		print "statistics:$stat samples\n";
		if(@info == 2){
#			print "@info";
#			print "$config->{$info[0]}\t$config->{$info[0]}{$info[1]}\n";
			if(exists $config->{$info[0]} and exists $config->{$info[0]}{$info[1]}){
				my $version = (split/\t/,$config->{$info[0]}{$info[1]})[0];
				if ($version =~ /MMC-quanjunpu/){
					my $cmd = "perl $Bin/$version/calculate_meta_level_per.pl -infile $out/report/$samID/$samID.extract.xls -outfile $out/report/$samID/$samID.level.xls";
					$cmd .= "  && perl $Bin/$version/get_result_xml.v1.pl -data $out/report/$samID -level_per $out/report/$samID/$samID.level.xls -out $out/report/$samID/pdf -barcode $samID";
					print OUT"$cmd\n";
				}
				elsif ($version =~ /MMC-yishengjun/){
					my $cmd = "perl $Bin/$version/calculate_meta_level_per.pl -infile $out/report/$samID/$samID.extract.xls -out $out/report/$samID/$samID.level.xls";
					$cmd .= " && perl $Bin/$version/get_result_xml.pl -csvfile $out/report/$samID/$samID.csv -barcode $samID -infile $out/report/$samID/$samID.level.xls -out $out/report/$samID/pdf/";
					print OUT"$cmd\n";
				}
				elsif($version =~ /KNYY/){
					my $cmd = "perl $Bin/$version/get_result_xml.pl -data $out/report/$samID -barcode $samID -out $out/report/$samID/pdf -testdate $testdate";
					$cmd .= " && perl $Bin/$version/meta_make_pdf_v2.pl -xml $out/report/$samID/pdf/$samID.xml -barcode $samID -od $out/report/$samID/pdf";
					print OUT"$cmd\n";
				}
				else{
					my $cmd = "perl $Bin/$version/get_result_xml.pl -data $out/report/$samID -barcode $samID -out $out/report/$samID/pdf";
					$cmd .= " && perl $Bin/$version/meta_make_pdf_v2.pl -xml $out/report/$samID/pdf/$samID.xml -barcode $samID -od $out/report/$samID/pdf >$out/report/$samID/$samID.o 2>$out/report/$samID/$samID.e";
					print OUT"$cmd\n";
				}
			}
			else{
				print STDERR"The $samID.csv have some problem,please check1!\n";
			}
		}
		else{
			print STDERR"The $samID.csv have some problem,please check2!\n";
		}
	}
}

sub MKDIR{
	my($dir) = @_;
	rmdir ($dir) if (-d $dir);
	mkdir ($dir) if (!-d $dir);
}

sub Absolutepath{
	my ($type,$in) = @_;
	my $return;
#	my $pwd = `pwd`;
#	chomp $pwd;
	if ($type eq 'dir'){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir $in;
		$return = `pwd`;
		chomp $return;
		chdir $pwd;
	}
	if ($type eq 'file'){
		my $pwd = `pwd`;
                chomp $pwd;

		my $dirname = dirname($in);
		my $filename = basename($in);
#		my $pwd = `pwd`;
		chdir $dirname;
		$return = `pwd`;
		chomp $return;
		$return .="\/".$filename;
		chdir $pwd;
	}
	return $return;
}

sub GETtime{
#	my($in) = @_;
	my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	$mon = $mon + 1;
	$year = $year + 1900;
	my $return = $year."-".$mon."-".$mday;

	return $return;
}
