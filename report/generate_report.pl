#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(cwd);
use FindBin qw($Bin);
use File::Basename qw(basename dirname);
use POSIX;
use Template;
use Data::Dumper;
use File::Slurp;
use Scalar::Util 'refaddr';
use XML::Simple qw(XMLin XMLout);
use Data::Printer colored => 1;
use utf8;
use Encode;

use Statistics::R;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

my $version = "1.0.0";
my $R = Statistics::R->new();

###
my ($datadir,$barcode_list,$dirout,$csvdir,$confdir,$generatepdf,$language,$help);
GetOptions (
	"data=s"	=> \$datadir,
	"barcode=s"	=> \$barcode_list,
	"od=s"		=> \$dirout,
	"csvdir=s"	=> \$csvdir,
	"confdir=s"	=> \$confdir,
	"generatepdf"	=> \$generatepdf,
	"language=s"	=> \$language,
	"h"	=> \$help,
) or &USAGE;
&USAGE if ((!defined $datadir || !defined $barcode_list) || defined $help) ;

########## AbsolutePath and default ##############
$dirout||="./";
$csvdir||="/data/bioit/biodata/mengf/Project/GIhealth_jichuban/SampleInfo";
$confdir||="$Bin/conf";
&MKDIR($dirout);
$dirout=&AbsolutePath("dir",$dirout);
$datadir=&AbsolutePath("dir",$datadir);
$barcode_list=&AbsolutePath("file",$barcode_list);
$confdir=&AbsolutePath("dir",$confdir);

my $work_sh = "$dirout/work_sh";
&MKDIR($work_sh);

$language||="CN";

my %check_key = (
	"name" => "1",
	"sex" => "1",
	"age" => "1",
	"height" => "1",
	"weight" => "1",
	"batch" => "1",
	"receivedate" => "1",
	"粪便质地" => "1",
	"粪便颜色" => "1",
	"version" => "1",
	"sucaifenzu" => "1"
);

my $get_result = "get_result_xml.pl";
my $xml2pdf = "meta_make_pdf_v2.pl";
my $get_result_qpcr = "get_result_xml.qpcr.pl";

################# main process ###################
my %sample;
&load_sample($barcode_list, \%sample);

## load report version conf
my %report_version;
#&load_conf_($confdir, \%report_version);
&load_conf_new($confdir, \%report_version);

my %incorrect_saminfo;
my %sample2report;
my %sample2source;
foreach my $barcode (sort keys %sample) {
	my %samInfo;
	&load_sample_info($barcode, $csvdir, \%samInfo, \%incorrect_saminfo, \%sample2source);
	&check_sample_info(\%check_key, \%samInfo, \%report_version, \%incorrect_saminfo, \%sample2report) if (exists $samInfo{barcode});
	&generate_readme(\%sample, \%sample2report, \%sample2source, $dirout);
}

&incorrect_info_log(\%incorrect_saminfo, $dirout);

&creat_report_dir(\%sample2report, $datadir, $dirout, $csvdir);

my $report_sh = "$work_sh/generate_report.sh";
&generate_report_sh(\%sample2report, \%report_version, $dirout, $report_sh);

if ($generatepdf) {
	system "parallel -j 10 <$report_sh >$dirout/report.log";
	&check_report_pages_num (\%sample2report, \%report_version, $dirout);
}

################# subs #################
sub generate_readme {
	my ($sample, $sample2report, $sample2source, $dirout) = @_;
	
	open OUT,">$dirout/Readme";
	print OUT "#Barcode\tSource\tVersion\tSucaifenzu\n";
	foreach my $barcode (sort keys %{$sample}){
		print OUT "$barcode\t";
		if(exists $sample2source->{$barcode}){
			my $xx = encode('utf8',$sample2source->{$barcode});
			print OUT "$xx\t";
		}else{
			print OUT "-\t";	
		}
		if(exists $sample2report->{$barcode}){
			my $xx = encode('utf8',$sample2report->{$barcode}{version});
			my $yy = encode('utf8',$sample2report->{$barcode}{sucaifenzu});
			print OUT "$xx\t$yy\n"
		}else{
			print OUT "-\t-\n";
		}
	}
}

sub check_report_pages_num {
	my ($sample2report, $report_version, $dirout) = @_;

	my %incorrect_report = ();
	my %sucess_report = ();
	my $error = 0;
	foreach my $barcode (sort keys %{$sample2report}) {
		my $report = "$dirout/$sample2report->{$barcode}{version}/$barcode/pdf/$barcode.pdf";
		my $foot_num = `pdftk $report dump_data|grep NumberOfPages| awk '{print \$2}'`;
		$foot_num=~s/^NumberOfPages: //;
		$foot_num=~s/\s+$//;
		if ($foot_num != $report_version->{$sample2report->{$barcode}{version}}{$sample2report->{$barcode}{sucaifenzu}}{PageNUM}) {
			$error = 1;
			$incorrect_report{$report_version->{$sample2report->{$barcode}{version}}{$sample2report->{$barcode}{sucaifenzu}}{DIR}}{$barcode} = $foot_num;
		}
		else {
			$sucess_report{$report_version->{$sample2report->{$barcode}{version}}{$sample2report->{$barcode}{sucaifenzu}}{DIR}}{$barcode} = 1;
		}
	}

	open OUT,">$dirout/report.error.list" || die $!;
	if ($error == 1) {
		print STDERR "Some sample's report failed (incorrected page number), please check!!!\n";
		foreach my $taocan (sort keys %incorrect_report) {
			print "Incorrect report page Number of $taocan:\n";
			foreach my $barcode (sort keys %{$incorrect_report{$taocan}}) {
				print "\t$barcode\t$incorrect_report{$taocan}{$barcode}\n";
			}
		}
	}
	close OUT;

}

sub incorrect_info_log {
	my ($incorrect_saminfo, $dirout) = @_;

	my $error = 0;
	open OUT,">$dirout/sampleinfo.error.list" || die $!;
	if (keys %{$incorrect_saminfo}) {
		$error = 1;
		foreach my $barcode (sort keys %{$incorrect_saminfo}) {
			print OUT "$barcode:\n";
			print OUT "$incorrect_saminfo->{$barcode}\n";
		}
	}
	close OUT;

	if ($error == 1) {
		print STDERR "WARNING: Some sample's info csv have errors, please check!!!\n\t$dirout/sampleinfo.error.list\n";
	}
	else {
		print STDOUT "Generate report command complete.\n";
	}
}

sub generate_report_sh {
	my ($sample2report, $report_version, $dirout, $report_sh) = @_;

	my $report_cmd;

	foreach my $barcode (sort keys %{$sample2report}) {
		my $scriptsdir = "$Bin/$report_version->{$sample2report->{$barcode}{version}}{$sample2report->{$barcode}{sucaifenzu}}{DIR}";
		$report_cmd .= "perl $scriptsdir/$get_result -data $dirout/$sample2report->{$barcode}{version}/$barcode -barcode $barcode -out $dirout/$sample2report->{$barcode}{version}/$barcode/pdf";
		$report_cmd .= " && perl $scriptsdir/$xml2pdf -xml $dirout/$sample2report->{$barcode}{version}/$barcode/pdf/$barcode.xml -barcode $barcode -od $dirout/$sample2report->{$barcode}{version}/$barcode/pdf";
		$report_cmd .= " >$dirout/$sample2report->{$barcode}{version}/$barcode/$barcode.o 2>$dirout/$sample2report->{$barcode}{version}/$barcode/$barcode.e\n";
	}

	open OUT,">$report_sh" || die $!;
	print OUT "$report_cmd";
	close OUT;

	return;
}

sub creat_report_dir {
	my ($sample2report, $datadir, $dirout, $csvdir) = @_;

#	my @repver = (sort values %{$sample2report});
#	my $uniqversion = &uniq_array(\@repver);
#	foreach my $ver (@{$uniqversion}) {
#		&MKDIR("$dirout/$ver");
#	}

	my $cp_cmd;
	foreach my $barcode (sort keys %{$sample2report}) {
		if (!-d "$datadir/$barcode/report") {
			print "Unexists $barcode report files dir:$datadir/$barcode/report, please Check!!!\n";
			next;
		}
		$cp_cmd .= "mkdir -p $dirout/$sample2report->{$barcode}{version}/$barcode && ";
		$cp_cmd .= "cp $datadir/$barcode/report/$barcode.*.xls $dirout/$sample2report->{$barcode}{version}/$barcode && ";
		$cp_cmd .= "cp $csvdir/$barcode.csv $dirout/$sample2report->{$barcode}{version}/$barcode\n";
	}

	system ($cp_cmd);
}

#my ($dis)=&uniq_array($Disease);
sub uniq_array {
	my ($array)=@_;
	my @dis=();
	
	##
	my $n=1;
	my %hash;
	foreach my $unit (@{$array}) {
		$hash{$unit}=$n++ if (! exists $hash{$unit});
	}
	##
	@dis=sort {$hash{$a} <=> $hash{$b}} keys %hash;
	
	return (\@dis);
}

sub load_sample {
	my ($barcode_list, $sample)=@_;

	open IN,"$barcode_list" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my $code = (split/\t+/,$_)[0];
		$sample->{$code} = 1;
	}
	close IN;
}



sub check_sample_info {
	my ($check_key, $samInfo, $report_version, $incorrect_saminfo, $sample2report) = @_;

	my $complete = 1;
	my $uncomplete_str = '-';
	my $history_str = '-';
	my $check_result = '-';
	my $version_correct = '-';
	my $sucai_correct = '-';

	foreach my $key (sort keys %{$check_key}) {
		if ($key eq "version") {
			if (!exists $report_version->{$samInfo->{$key}}) {
				$version_correct = $samInfo->{$key};
			}
		}
		if ($key eq "sucaifenzu") {
			if (!exists $report_version->{$samInfo->{version}}{$samInfo->{$key}}) {
				$sucai_correct = $samInfo->{$key};
				print "\n\n$key\n$samInfo->{sucaifenzu}\n$samInfo->{barcode}\n\n";
			}
		}
		elsif ($key eq "粪便质地" || $key eq "粪便颜色"){
			$complete = 1;
		}
		elsif ($key ne "batch") {
			if ($samInfo->{$key} eq "") {
				$uncomplete_str .= "$key\t";
				$complete = 0;
			}
		}
		else {
			if ($samInfo->{$key} > 1 && $samInfo->{$key} < 4) {
				my $history_code = $samInfo->{"历次检测编号"};
				my @codes = split/,/,$history_code;
				if (@codes != $samInfo->{$key} - 1) {
					$history_str .= "$samInfo->{$key}\t$history_code";
				}
			}
			elsif ($samInfo->{$key} >= 4) {
				my $history_code = $samInfo->{"历次检测编号"};
				my @codes = split/,/,$history_code;
				if (@codes != 3) {
					$history_str .= "$samInfo->{$key}\t$history_code";
				}
			}
		}
	}

	if ($complete == 0) {
		$check_result .= "$samInfo->{barcode} sample info csv isn't complete, please check $uncomplete_str\n";
	}
	if ($history_str ne '-') {
		$check_result .= "$samInfo->{barcode} sample test history isn't correct, please check $history_str\n";
	}
	if ($version_correct ne '-') {
		$check_result .= "$samInfo->{barcode} report version $version_correct unexists, please check\n";
	}
	if ($sucai_correct ne '-') {
		$check_result .= "$samInfo->{barcode} report sucaifenzu $sucai_correct unexists, please check\n";
	}

	if ($check_result ne '-') {
		$incorrect_saminfo->{$samInfo->{barcode}} = $check_result;
	}
	else {
		$sample2report->{$samInfo->{barcode}}{version} = $samInfo->{version};
		$sample2report->{$samInfo->{barcode}}{sucaifenzu} = $samInfo->{sucaifenzu};
	}

	return;
}

sub load_conf_ {
	my ($confdir, $report_version) = @_;

	open IN,"$confdir/report-version.conf" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($code, $report_v, $page_num) = split/\t+/,$_;
		my @codes = split/,/,$code;
		foreach (@codes) {
			$report_version->{$_}{DIR} = $report_v;
			$report_version->{$_}{PageNUM} = $page_num;
		}
	}
	close IN;
}

sub load_conf_new {
	my ($confdir, $report_version) = @_;

	open IN,"$confdir/report-version.new.conf" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($code, $sucai, $report_v, $page_num) = split/\t+/,$_;
		my @codes = split/,/,$code;
		foreach (@codes) {
			$report_version->{$_}{$sucai}{DIR} = $report_v;
			$report_version->{$_}{$sucai}{PageNUM} = $page_num;
		}
	}
	close IN;
}

#$check+=&checkF($file);
sub checkF {
	my ($file)=@_;
	
	##
	my $check=0;
	if (! -f $file) {
		print STDOUT "file not exists : $file\n";
		$check=1;
	}
	
	return($check);
}

#my ($samInfo)=&load_sample_info($barcode,$datadir,\%vars);
sub load_sample_info {
	my ($barcode,$csvdir,$samInfo,$incorrect_saminfo,$sample2source)=@_;

	##
	my $fatN=-1;
	$/="\n";
	my $file="$csvdir/$barcode".".csv";
	if (!-f $file) {
		$incorrect_saminfo->{$barcode} = "unexists $file, please check\n";
		return;
	}
	open IN,$file or die $!;
	my $head=<IN>;
	my @head=split /\s+/,$head;
	while (<IN>) {
		chomp;
		next if (/^$/);
		my @unit=split /\t/,$_;
		for (my $i=0;$i<@head;$i++) {
			if ($head[$i] eq "barcode" || $head[$i] eq "sex" || $head[$i] eq "age" || $head[$i] eq "weight" || $head[$i] eq "height" || $head[$i] eq "batch" || $head[$i] eq "receivedate") 
			{
				$samInfo->{$head[$i]}=$unit[$i];
			}
			elsif ($head[$i] eq "name") {
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$samInfo->{$head[$i]}=$unit[$i];
			}
			elsif ($head[$i] eq "source") {
				if ($unit[$i] eq "") {
					$unit[$i]=decode("UTF-8","-");
				}
				else {
					$unit[$i]=decode("UTF-8",$unit[$i]);
				}
				$samInfo->{$head[$i]}=$unit[$i];
			}
			elsif ($head[$i] eq "version") {
				$unit[$i]=decode("UTF-8",$unit[$i]);
				$samInfo->{$head[$i]}=$unit[$i];
			}
			else {
				$head[$i]=decode("UTF-8",$head[$i]);
				if ($head[$i] =~ /素材分组/) {
					$unit[$i]=decode("UTF-8",$unit[$i]);
					$samInfo->{'sucaifenzu'}=$unit[$i];
				}
				$head[$i]="2型糖尿病" if ($head[$i] =~ /糖尿病/);
				$head[$i]="炎症性肠病" if ($head[$i] =~ /肠炎/ || $head[$i] =~ /炎症性肠病/ || $head[$i] =~ /腹泻/);
				$head[$i]="粪便质地" if ($head[$i] =~ /质地/);
				$head[$i]="粪便颜色" if ($head[$i] =~ /颜色/);
				if ($head[$i] eq "粪便质地" || $head[$i] eq "粪便颜色" || $head[$i] =~ /历次检测编号/) {
					$samInfo->{$head[$i]}=decode("UTF-8",$unit[$i]);
				}
				else {
					$samInfo->{disease}{$head[$i]}{'risk'}=decode("UTF-8",$unit[$i]);
				}
			}
			$sample2source->{$samInfo->{barcode}} = $samInfo->{source};
		}
	}
	close IN;
	
	## BMI
	my $BMI="NA";
	$BMI=sprintf "%.2f",$samInfo->{weight}/(($samInfo->{height}/100)**2) if ($samInfo->{height} != 0);
	my $fat="肥胖";
	if ($BMI ne "NA" && $BMI > 24) {
		$samInfo->{disease}{$fat}{risk}='Y';
	}
	else {
		$samInfo->{disease}{$fat}{risk}='N';
	}
	
	return;
}

sub MKDIR { # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}

sub AbsolutePath
{       #获取指定目录或文件的决定路径
        my ($type,$input) = @_;

        my $return;
        $/="\n";

        if ($type eq 'dir')
        {
                my $pwd = `pwd`;
                chomp $pwd;
                chdir($input);
                $return = `pwd`;
                chomp $return;
                chdir($pwd);
        }
        elsif($type eq 'file')
        {
                my $pwd = `pwd`;
                chomp $pwd;

                my $dir=dirname($input);
                my $file=basename($input);
                chdir($dir);
                $return = `pwd`;
                chomp $return;
                $return .="\/".$file;
                chdir($pwd);
        }
        return $return;
}

sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
    $wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub USAGE {#
	my $usage=<<"USAGE";

    Program Function: GIhealth result to xml format;
    Version:    $version
    Contact:    fred routine <fred_routine\@163.com> 
    Program Date:   2016.11.15
    Usage:
      Options:
      -data        <dir>   sample's report result dir for input          forced
      -barcode     <str>   sample's code                                 forced

      -csvdir     <dir>   product set database list and xml dir     optional
                          default "/data/bioit/biodata/mengf/Project/GIhealth_jichuban/SampleInfo";
      -confdir    <dir>   report conf dir                           optional
                          default "$Bin/conf";
      -od         <dir>   output file dir                           optional
                          default "./";
      -language   <str>   result language used for report           optional
                          default "CN";

      -generatepdf        run generate report pdf cmd, default "only generate report command";
      -h          Help

USAGE
	print $usage;
	exit;
}
