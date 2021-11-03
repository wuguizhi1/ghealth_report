#!/usr/bin/perl -w
use strict;
use warnings;
use PerlIO::gzip;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use Data::Dumper;

my $BEGIN_TIME = time();
my $version = "1.0.0";
#######################################################################################
# ------------------------------------------------------------------
# GetOptions, Check format, Output info.
# ------------------------------------------------------------------
my ($fc,$data_dir,$lib_cfg,$barcode_cfg,$para_cfg,$basecall,$od,$bmask1,$bmask2);
GetOptions(
		"help|?" =>\&USAGE,
		"fc=s" =>\$fc,
		"data=s" =>\$data_dir,

		"lib-cfg=s" => \$lib_cfg,
		"barcode=s" => \$barcode_cfg,
		"para=s" => \$para_cfg,
		"od=s" => \$od,
		"bmask1=s" => \$bmask1,
		"bmask2=s" => \$bmask2,
	) or &USAGE;
&USAGE unless ($fc and $lib_cfg) ;

$data_dir = $data_dir || "/data/win";
$data_dir = &ABSOLUTE_DIR($data_dir);
$lib_cfg = &ABSOLUTE_DIR($lib_cfg);

$barcode_cfg = $barcode_cfg || "$Bin/../conf/barcode.cfg";
$barcode_cfg = &ABSOLUTE_DIR($barcode_cfg);
$od = $od || "./";
&MKDIR($od);
$od = &ABSOLUTE_DIR($od);
my $work_sh = "$od/work_sh";
&MKDIR ($work_sh);

$para_cfg = $para_cfg || "$Bin/../conf/para.cfg";
$para_cfg = &ABSOLUTE_DIR($para_cfg);

$bmask1 = $bmask1 || "y*,i7,i7,y*";
$bmask2 = $bmask2 || "y*,i10,y*";
############################
my $bcl2fastq = "$Bin/../utils/bcl2fastq";
my $fastqcount = "$Bin/../utils/fastq_count";

######
my $miseq_com = 0;
$miseq_com = &check_RTACom_ ($data_dir, $fc);

die if ($miseq_com != 1);

#### load cfgs
my %para;
&load_para_cfg_($para_cfg, \%para);

my %index2samples;
my %index2uniqsam;
my %onestepindex2samples;
#my %onestepindex2uniqsam;
#my ($uniqsamplesheet, $onestepsamplesheet) = &load_smaple_info_ ($lib_cfg, \%index2samples, \%index2uniqsam, \%onestepindex2samples, \%onestepindex2uniqsam);
my ($uniqsamplesheet, $onestepsamplesheet) = &load_smaple_info_ ($lib_cfg, \%index2samples, \%index2uniqsam, \%onestepindex2samples);

# bcl2fastq step
if (!-f "$work_sh/bcl2fastq.OK") {
	system "mkdir $od/$fc && ln -s $data_dir/$fc/* $od/$fc";
	system "cp $od/$fc/SampleSheet.csv $od/$fc/SampleSheet.csv.bac && rm $od/$fc/SampleSheet.csv" if (-f "$od/$fc/SampleSheet.csv");

	my $bcl2fastq_sh;
	if ($uniqsamplesheet ne "" && !-f "$work_sh/bcl2fastq1.OK") {
		&replace_samplesheet_ ($od, $fc, $uniqsamplesheet);
		$bcl2fastq_sh = "$work_sh/bcl2fastq1.sh" || die $!;
		&basecall_ ($od, $fc, "$od/rawdata_i7i5", $bcl2fastq_sh, $para{"barcode-mismatches"}, $bmask1);
		system "mv $od/$fc/SampleSheet.csv $od/$fc/SampleSheet.csv.i7i5 && touch $work_sh/bcl2fastq1.OK";

	}

	if ($onestepsamplesheet ne "" && !-f "$work_sh/bcl2fastq2.OK") {
		&replace_samplesheet_ ($od, $fc, $onestepsamplesheet);
		$bcl2fastq_sh = "$work_sh/bcl2fastq2.sh" || die $!;
		&basecall_ ($od, $fc, "$od/rawdata_i10", $bcl2fastq_sh, $para{"barcode-mismatches"}, $bmask2);
		system "mv $od/$fc/SampleSheet.csv $od/$fc/SampleSheet.csv.i10 && touch $work_sh/bcl2fastq2.OK";
	}

	system "touch $work_sh/bcl2fastq.OK";
}

# split sample reads by second barcode pair
if (!-f "$work_sh/split_library_reads.OK") {
	#my $barcode_identify_stat = "$od/$fc.barcode.stat.xls";
	my $barcode_pair_stat = "$od/$fc.barcodepair.stat.xls";
	my $barcode_single_stat = "$od/$fc.barcodesingle.stat.xls";


	if ($uniqsamplesheet ne "-") {
		my $split_lib_sh = "$work_sh/split_library_i7i5.sh";
		&generate_split_sh_ ("$od/rawdata_i7i5", \%index2samples, \%index2uniqsam, "$od/sample_data", $split_lib_sh, $barcode_cfg, $barcode_pair_stat);

		system "cat $split_lib_sh | parallel -j 10";
		system "touch $work_sh/split_library_i7i5.OK";
	}

	if ($onestepsamplesheet ne "-") {
		my $cut_lib_sh = "$work_sh/cut_library_i10.sh";
		#&generate_split_sh_ ("$od/rawdata_i10", \%onestepindex2samples, \%onestepindex2uniqsam, "$od/sample_data", $split_lib_sh, $barcode_cfg, $barcode_identify_stat);
		&generate_cut_sh_ ("$od/rawdata_i10", \%onestepindex2samples, "$od/sample_data", $cut_lib_sh, $barcode_cfg, $barcode_single_stat);

		system "cat $cut_lib_sh | parallel -j 10";
		system "touch $work_sh/cut_library_i10.OK";
	}
	system "touch $work_sh/split_library_reads.OK";
}

# stat samples reads number and quality (Q20, Q30)
if (!-f "$work_sh/sample_qual_stat.OK") {
	my $qual_stat_sh = "$work_sh/sample_qual_stat.sh";
	&generate_stat_sh_ ("$od/sample_data", \%index2samples, \%onestepindex2samples, "$od/stat", $qual_stat_sh);

	system "cat $qual_stat_sh | parallel -j 10";

	my $stat_all = "$od/$fc.allsamples.stat.xls";
	my $stat_flowcell = "$od/$fc.summary.xls";
	&cat_stat_ ("$od/stat", "$od/$fc.barcodepair.stat.xls", "$od/$fc.barcodesingle.stat.xls", \%index2samples, \%onestepindex2samples, \%para, $stat_all, $stat_flowcell);

	system "touch $work_sh/sample_qual_stat.OK";
}

# clean and clear up the directory
&clean_clearup_dir_($fc,$od);

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################
#+------------------------------------
# subs
#+------------------------------------
sub clean_clearup_dir_ {
	my ($fc, $od) = @_;

	# clean link dir of miseq run data
	system "rm -r $od/$fc" if (-d "$od/$fc");
	# clean unfiltered rawdata
	#system "rm -r $od/rawdata_*";
}

sub basecall_ {
	my ($rundir, $fc, $out, $sh, $mismatch, $bmask) = @_;

	&MKDIR ($out);
	my $cmd;
	$cmd = "cd $rundir/$fc && $bcl2fastq -i $rundir/$fc/Data/Intensities/BaseCalls -R $rundir/$fc ";
	#$cmd .= "--interop-dir $out/InterOp -o $out --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --barcode-mismatches $mismatch";
	$cmd .= "--interop-dir $out/InterOp -o $out --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --barcode-mismatches $mismatch";
	if ($bmask ne "") {
		$cmd .= " --use-bases-mask $bmask";
	}
	$cmd .= ">$out/bcl2fastq.log 2>&1 || touch bcl2fastq.err\n";

	open OUT,">$sh" || die $!;
	print OUT "$cmd";
	close OUT;

	system ("sh $sh");
}

sub replace_samplesheet_ {
	my ($rundir, $flowcell, $samplesheet) = @_;

	open IN,"$rundir/$flowcell/SampleSheet.csv.bac" || die $!;
	my $newcsv;
	while (<IN>) {
		chomp;
		if ($_ !~ /^Sample_ID,Sample_Name/) {
			$newcsv .= "$_\n";
		}
		else {
			$newcsv .= "$_\n$samplesheet";
			last;
		}
	}
	close IN;

	open OUT,">$rundir/$flowcell/SampleSheet.csv" || die $!;
	print OUT "$newcsv\n";
	close OUT;

}

sub check_RTACom_ {
	my ($rundir, $flowcell) = @_;

	RUNCOMPLETE:
	{
		my $com_file = "$rundir/$flowcell/RTAComplete.txt";
		if (!-f $com_file) {
				print STDERR "NS500 is runing $flowcell data...\n";
				sleep (3600);
				redo RUNCOMPLETE;
		}
		last RUNCOMPLETE if (-f $com_file);
	}

	return("1");
}

sub load_para_cfg_ {
	my ($cfg, $para) = @_;

	open IN,"$cfg" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($key,$val) = split/\s+/,$_;
		$para->{$key} = $val;
	}
	close IN;
}

sub load_smaple_info_ {
	my ($file,$info,$uniq,$onestepinfo) = @_;

	my %i7;
	my %i5;
	my %i10;
	my $onestep = 0;
	my $twostep = 0;
	open IN,"$file" || die $!;
	while (<IN>) {
		chomp;
		next if (/^$/ || /^\#/);
		my @tmp = split /\s+/,$_;
		if (@tmp == 6) {
			$twostep = 1;
			$tmp[1] =~ /F(\d+)R(\d+)/;
			my $index = "$tmp[2]~$tmp[4]";
			while ($tmp[1]=~/F(\d+)R(\d+)/ig) {
				die "library cfg file contain's unpaired $tmp[1]\n" if ($1 ne $2);
				my $barcode_num = $1;
				if (exists $info->{$index}) {
					foreach my $sam (sort keys %{$info->{$index}}) {
						if (exists $info->{$index}{$sam}{$barcode_num}) {
							print "$sam and $tmp[0] conflicting with index and second barcode info, please check\n\n";
							die;
						}
					}
				}
				$info -> {$index}{$tmp[0]}{$barcode_num} = 1;
			}
			$uniq -> {$index} = $tmp[0];
			$i7{$tmp[2]} = $tmp[3];
			$i5{$tmp[4]} = $tmp[5];
		}
		if (@tmp == 4) {
			$onestep = 1;
			my $index = "$tmp[2]";
			#while ($tmp[1]=~/F(\d+)R(\d+)/ig) {
			#	die "library cfg file contain's unpaired $tmp[1]\n" if ($1 ne $2);
			#	my $barcode_num = $1;
			#	if (exists $onestepinfo->{$index}) {
			#		foreach my $sam (sort keys %{$onestepinfo->{$index}}) {
			#			if (exists $onestepinfo->{$index}{$sam}{$barcode_num}) {
			#				print "$sam and $tmp[0] conflicting with index and second barcode info, please check\n\n";
			#				die;
			#			}
			#		}
			#	}
			#	$onestepinfo -> {$index}{$tmp[0]}{$barcode_num} = 1;
			#}
			if (exists $onestepinfo->{$index}) {
				print "$tmp[0] and $onestepinfo->{$index} conflicting with index, please check\n\n";
				die;
			}
			$onestepinfo -> {$index} = $tmp[0];
			$i10{$index} = $tmp[3];
		}
	}
	close IN;

	my $samplesheet = "";
	if ($twostep == 1) {
		foreach my $index_pair (sort keys %$uniq) {
			my ($index_7, $index_5) = split/~/,$index_pair;
			$samplesheet .= "$uniq->{$index_pair},$uniq->{$index_pair},,,$index_7,$i7{$index_7},$index_5,$i5{$index_5},,\n";
		}
	}

	my $onestepsamplesheet = "";
	if ($onestep == 1) {
		foreach my $index (sort keys %$onestepinfo) {
			$onestepsamplesheet .= "$onestepinfo->{$index},$onestepinfo->{$index},,,$index,$i10{$index},,\n";
		}
	}

	return ($samplesheet, $onestepsamplesheet);
}

sub generate_split_sh_ {
	my ($data, $i2sam, $i2uniq, $out_dir, $out_sh, $cfg, $stat) = @_;

	&MKDIR ($out_dir);
	open OUT,">$out_sh" || die $!;
	foreach my $index (sort keys %$i2sam) {
		my ($sam_str,$barcode_info);
		my $read1 = `find $data/$i2uniq->{$index}_*R1*fastq.gz `; chomp $read1;
		my $read2 = `find $data/$i2uniq->{$index}_*R2*fastq.gz `; chomp $read2;
		foreach my $sam (sort keys %{$i2sam->{$index}}) {
			my $barcode_index = $i2sam->{$index}{$sam};
			$sam_str .= "$sam;";
			my $barcode_str = join",",keys %{$i2sam->{$index}{$sam}};
			$barcode_info .= "$barcode_str;";
		}
		$sam_str =~ s/;$//;
		$barcode_info =~ s/;$//;
		print OUT "perl $Bin/../scripts/sample_data_extract.pl \"$sam_str\" $cfg \"$barcode_info\" $read1 $read2 $out_dir $stat\n";
	}
	close OUT;
}

sub generate_cut_sh_ {
	my ($data, $i2sam, $out_dir, $out_sh, $cfg, $stat) = @_;

	&MKDIR ($out_dir);
	open OUT,">$out_sh" || die $!;
	foreach my $index (sort keys %$i2sam) {
		my $sam_str = $i2sam->{$index};
		my $read1 = `find $data/$i2sam->{$index}_*R1*fastq.gz `; chomp $read1;
		my $read2 = `find $data/$i2sam->{$index}_*R2*fastq.gz `; chomp $read2;
		print OUT "perl $Bin/../scripts/sample_mix_cut.pl \"$sam_str\" $cfg $read1 $read2 $out_dir $stat\n";
	}
	close OUT;
}

sub generate_stat_sh_ {
	my ($data, $i2sam, $onestepi2sam, $out_dir, $out_sh) = @_;

	&MKDIR ($out_dir);
	open OUT,">$out_sh" || die $!;
	foreach my $index (sort keys %$i2sam) {
		foreach my $sam (sort keys %{$i2sam->{$index}}) {
			my $sam_pre = "$sam"."_";
			my $read1 = `find $data/$sam_pre*R1*fastq `; chomp $read1;
			my $read2 = `find $data/$sam_pre*R2*fastq `; chomp $read2;
			print OUT "$fastqcount $read1 $read2 -o $out_dir/$sam.stat.xls\n";
		}
	}

	foreach my $index (sort keys %$onestepi2sam) {
		my $sam = $onestepi2sam->{$index};
		my $sam_pre = "$sam"."_";
		my $read1 = `find $data/$sam_pre*R1*fastq `; chomp $read1;
		my $read2 = `find $data/$sam_pre*R2*fastq `; chomp $read2;
		print OUT "$fastqcount $read1 $read2 -o $out_dir/$sam.stat.xls\n";
	}

	close OUT;
}

sub cat_stat_ {
	my ($data, $barcodepair_stat, $barcodesingle_stat, $i2sam, $onestepi2sam, $para, $out1, $out2) = @_;

	my $highdata_threshold = $para->{"highdata-threshold"}; 
	my $lowdata_threshold = $para->{"lowdata-threshold"};
	open OUT,">$out1" || die $!;
	print OUT "#sample\tReadNumber\tBaseNumber\tQ20(%)\tQ30(%)\n";
	my ($total_reads,$total_base,$total_q20,$total_q30);
	my ($total_r1_base,$total_r2_base,$total_r1_q20,$total_r1_q30,$total_r2_q20,$total_r2_q30);
	my @lowdata;
	my $lowdata_num = 0;
	my @highdata;
	my $highdata_num = 0;
	foreach my $index (sort keys %$i2sam) {
		foreach my $sam (sort keys %{$i2sam->{$index}}) {
			open IN,"$data/$sam.stat.xls" || die $!;
			my ($rnum, $bnum, $q20, $q30);
			while (<IN>) {
				chomp;
				next if (/^$/);
				my @tmp = split/\t+/,$_;
				$rnum = $tmp[1];
				$bnum += $tmp[2];
				if ($rnum == 0) {
					print OUT "$sam\t0\t0\tNA\tNA\n";
					push @lowdata, "$sam\t0\t0\tNA\tNA";
					next;
				}
				$q20 += $tmp[2]*$tmp[6]/100;
				$q30 += $tmp[2]*$tmp[7]/100;
				if ($tmp[0]=~/_R1.fastq/) {
					$total_r1_base += $tmp[2];
					$total_r1_q20 += $tmp[2]*$tmp[6]/100;
					$total_r1_q30 += $tmp[2]*$tmp[7]/100;
				}
				if ($tmp[0]=~/_R2.fastq/) {
					$total_r2_base += $tmp[2];
					$total_r2_q20 += $tmp[2]*$tmp[6]/100;
					$total_r2_q30 += $tmp[2]*$tmp[7]/100;
				}
			}
			next if ($bnum == 0);
			$total_reads += $rnum;
			$total_base += $bnum;
			$total_q20 += $q20;
			$total_q30 += $q30;
			my $q20_per = sprintf "%.2f", $q20/$bnum*100;
			my $q30_per = sprintf "%.2f", $q30/$bnum*100;
			print OUT "$sam\t$rnum\t$bnum\t$q20_per\t$q30_per\n";
			if ($rnum < $lowdata_threshold) {
				push @lowdata,"$sam\t$rnum\t$bnum\t$q20_per\t$q30_per";
				$lowdata_num ++;
			}
			if ($rnum > $highdata_threshold) {
				push @highdata,"$sam\t$rnum\t$bnum\t$q20_per\t$q30_per";
				$highdata_num ++;
			}
		}
	}

	foreach my $index (sort keys %$onestepi2sam) {
		my $sam = $onestepi2sam->{$index};
		open IN,"$data/$sam.stat.xls" || die $!;
		my ($rnum, $bnum, $q20, $q30);
		while (<IN>) {
			chomp;
			next if (/^$/);
			my @tmp = split/\t+/,$_;
			$rnum = $tmp[1];
			$bnum += $tmp[2];
			if ($rnum == 0) {
				print OUT "$sam\t0\t0\tNA\tNA\n";
				push @lowdata, "$sam\t0\t0\tNA\tNA";
				next;
			}
			$q20 += $tmp[2]*$tmp[6]/100;
			$q30 += $tmp[2]*$tmp[7]/100;
			if ($tmp[0]=~/_R1.fastq/) {
				$total_r1_base += $tmp[2];
				$total_r1_q20 += $tmp[2]*$tmp[6]/100;
				$total_r1_q30 += $tmp[2]*$tmp[7]/100;
			}
			if ($tmp[0]=~/_R2.fastq/) {
				$total_r2_base += $tmp[2];
				$total_r2_q20 += $tmp[2]*$tmp[6]/100;
				$total_r2_q30 += $tmp[2]*$tmp[7]/100;
			}
		}
		next if ($bnum == 0);
		$total_reads += $rnum;
		$total_base += $bnum;
		$total_q20 += $q20;
		$total_q30 += $q30;
		my $q20_per = sprintf "%.2f", $q20/$bnum*100;
		my $q30_per = sprintf "%.2f", $q30/$bnum*100;
		print OUT "$sam\t$rnum\t$bnum\t$q20_per\t$q30_per\n";
		if ($rnum < $lowdata_threshold) {
			push @lowdata,"$sam\t$rnum\t$bnum\t$q20_per\t$q30_per";
			$lowdata_num ++;
		}
		if ($rnum > $highdata_threshold) {
			push @highdata,"$sam\t$rnum\t$bnum\t$q20_per\t$q30_per";
			$highdata_num ++;
		}
	}
	close OUT;

	open OUT,">$out2" || die $!;
	my $total_q20per = sprintf "%.2f", $total_q20/$total_base*100;
	my $total_q30per = sprintf "%.2f", $total_q30/$total_base*100;
	my $total_r1_q20per = sprintf "%.2f", $total_r1_q20/$total_r1_base*100;
	my $total_r1_q30per = sprintf "%.2f", $total_r1_q30/$total_r1_base*100;
	my $total_r2_q20per = sprintf "%.2f", $total_r2_q20/$total_r2_base*100;
	my $total_r2_q30per = sprintf "%.2f", $total_r2_q30/$total_r2_base*100;
	print OUT "Total Reads:\t$total_reads\nTotal Base:\t$total_base\n\n";
	print OUT "#QuaStat\tQ20(%)\tQ30(%)\n";
	print OUT "Total\t$total_q20per\t$total_q30per\nRead1\t$total_r1_q20per\t$total_r1_q30per\nRead2\t$total_r2_q20per\t$total_r2_q30per\n\n";

	print OUT "Secondary barcode stat:\n";
	if (-f $barcodepair_stat) {
		my %barcode_pair_stat;
		&barcode_stat_ ($barcodepair_stat, \%barcode_pair_stat);
		my ($pair_str, $pair_count);
		foreach my $pair (sort keys %barcode_pair_stat) {
			$pair_str .= "\t$pair";
			$pair_count .= "\t$barcode_pair_stat{$pair}";
		}
		print OUT "barcode_pair$pair_str\nRead_number$pair_count\n\n";
	}
	if (-f $barcodesingle_stat) {
		my %barcode_single_stat;
		&barcode_stat_ ($barcodesingle_stat, \%barcode_single_stat);
		my ($single_str, $single_count);
		foreach my $single (sort keys %barcode_single_stat) {
			$single_str .= "\t$single";
			$single_count .= "\t$barcode_single_stat{$single}";
		}
		print OUT "barcode_pair$single_str\nRead_number$single_count\n\n";
	}

	print OUT "Lowdata Sample Number(Reads number < $lowdata_threshold):\t$lowdata_num\n";
	if ($lowdata_num > 0) {
		print OUT "LowdataSample\tReadNumber\tBaseNumber\tQ20(%)\tQ30(%)\n";
		foreach (@lowdata) {
			print OUT "$_\n";
		}
	}
	print OUT "\nHighdata sample Number(Reads number > $highdata_threshold):\t$highdata_num\n";
	if ($highdata_num > 0) {
		print OUT "HighdataSample\tReadNumber\tBaseNumber\tQ20(%)\tQ30(%)\n";
		foreach (@highdata) {
			print OUT "$_\n";
		}
	}
	close OUT;
}

#&barcode_stat_ ($file, \%barcode_stat);
sub barcode_stat_ {
	my ($file, $barcode_stat) = @_;

	open IN,"$file" || die $!;
	my $head = <IN>; chomp $head;
	my @pairs = split/\t/,$head;
	while (<IN>) {
		chomp;
		s/^\s+\|\s+$//;
		next if (/^$/ || /^\#/);
		my @info = split/\t/,$_;
		for (1..$#info) {
			$barcode_stat->{$pairs[$_]} += $info[$_];
		}
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

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";

    Program Function: split berry's 16srDNA multiple second barcode linrary's reads;
    Version:    $version
    Contact:    fred routine <fred_routine\@163.com> 
    Program Date:   2014.11.24
    Usage:
      Options:
      -lib-cfg    <file>  smaples index and second barcode info file,   forced
                          eg. "$Bin/../conf/library.cfg"
      -fc         <flowcell>  miseq run flowcell info,                  forced
                          eg. "160912_M03342_0010_000000000-AULPU";

      -data       <dir>   sequencing data's directory,                  optical
                          default "/data/win";
      -barcode    <file>  seconde barcode info file,                    optical
                          default "$Bin/../conf/barcode.cfg"
      -para       <file>  parameters for bcl2fastq config,              optical
                          default "$Bin/../conf/para.cfg"

      -od         <dir>   output dir,                                   optical
                          default "./"

      -bmask1     <str>   output dir,                                   optical
                          default "y*,i7,i7,y*"
      -bmask2     <str>   output dir,                                   optical
                          default "y*,i10,y*"

      -h          Help

USAGE
	print $usage;
	exit;
}
