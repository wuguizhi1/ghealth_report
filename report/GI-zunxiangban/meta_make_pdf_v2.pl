#!/usr/bin/perl
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
binmode(STDOUT, ":utf8");
binmode(STDERR, ":utf8");

my $version = "1.0.0";

my ($xmlfile,$templateDir,$referdir,$dirout,$config,$barcode,$simconfig);
GetOptions (
	"xml=s"		=> \$xmlfile,
	"template=s"	=> \$templateDir,
	"referdir=s"	=> \$referdir,
	"config=s"	=> \$config,
	"simconfig=s"	=> \$simconfig,
	"od=s"		=> \$dirout,
	"barcode=s"	=> \$barcode,
);
&USAGE if (!defined $xmlfile || !defined $barcode || !defined $dirout);

########## AbsolutePath and default ##############
mkdir $dirout if (! -d $dirout);
$dirout=&AbsolutePath("dir",$dirout);
$xmlfile=&AbsolutePath("file",$xmlfile);
$config||="$Bin/config.all";
$config=&AbsolutePath("file",$config);
$simconfig||="$Bin/config.simple";
$simconfig=&AbsolutePath("file",$simconfig);
$templateDir||="$Bin/template/current";
$templateDir=&AbsolutePath("dir",$templateDir);
$referdir||="$Bin/../../reference_data/current";
$referdir=&AbsolutePath("dir",$referdir);

my $tpl = Template->new({ENCODING => 'utf8',ABSOLUTE => 1, EVAL_PERL => 1});


################# main process ###################

## read xml file
my $xml=&read_xml($xmlfile);

## cp
&cp($templateDir,$referdir,$dirout);

my (@context1,@context2,@context3,@context4,@context5,@context6);
#### 致辞
&make_front($xml,"zhici",$templateDir,$dirout);
system "cp $dirout/PDF/zhici.pdf $dirout/$barcode"."_zhici.pdf";
#### 扉页
&make_front($xml,"front",$templateDir,$dirout);

#### 检测结果总结
&make_pdf($xml,"zongjie",$templateDir,$dirout,\@context1);

#### 消化和吸收
&make_pdf($xml,"xiaohuahexishou",$templateDir,$dirout,\@context2);
&make_pdf($xml,"xiaohuahexishouP2",$templateDir,$dirout,\@context2);
&make_pdf($xml,"xiaohuahexishouP3",$templateDir,$dirout,\@context2);

#### 炎症和免疫
&make_pdf($xml,"yanzhenghemianyi",$templateDir,$dirout,\@context3);

#### 肠道菌群
&make_pdf($xml,"changdaojunqun",$templateDir,$dirout,\@context4);

#### 致病菌、相关疾病、粪便状态
&make_pdf($xml,"changdaojunqunP2",$templateDir,$dirout,\@context5);
&make_pdf($xml,"changdaojunqunP3",$templateDir,$dirout,\@context5);
&make_pdf($xml,"fenbianzhuangtai",$templateDir,$dirout,\@context5);

### 健康建议
&make_pdf($xml,"zongshuoming",$templateDir,$dirout,\@context6);
&make_pdf($xml,"shanshifangan",$templateDir,$dirout,\@context6);
&make_pdf($xml,"changdaotiaojieyuyundongfangan",$templateDir,$dirout,\@context6);

#### output context
&output_contexts("$dirout/contexts1",\@context1);
&output_contexts("$dirout/contexts2",\@context2);
&output_contexts("$dirout/contexts3",\@context3);
&output_contexts("$dirout/contexts4",\@context4);
&output_contexts("$dirout/contexts5",\@context5);
&output_contexts("$dirout/contexts6",\@context6);

#### final pdf
#&make_finalPDF($templateDir,$dirout,$config,$barcode);
my $personal_prefix = "$barcode";
&make_finalPDF($templateDir,$dirout,$config,$personal_prefix);
my $sim_prefix = "$barcode"."_sim";
&make_finalPDF($templateDir,$dirout,$simconfig,$sim_prefix);

############################# sub ###############################################

#&make_front($xml,"front",$templateDir,$dirout);
sub make_front {
	my ($xml,$key,$tempDir,$dirout)=@_;

	##
	my $tt="$templateDir/$key".".tt";
	my $ttTex=$key.".tex";
	if (-f $tt) {
		$xml->{'reportdate'}=sub_format_datetime(localtime(time()));
		#$xml->{'reportdate'}="2016-11-10";
		&section_tex($tt,$xml,"$dirout/$ttTex");
		chdir $dirout;
		mkdir "./PDF" if (! -d "./PDF");
		`xelatex $ttTex`;
		`mv $key.tex $key.pdf $key.log ./PDF`;
		`rm $key.out $key.aux`;
		print STDOUT "------------------\n$key.pdf\n";
		chdir $Bin;
	}
	
	return;
}

#&make_finalPDF($templateDir,$dirout,$config,$barcode);
sub make_finalPDF {
	my ($templateDir,$dirout,$config,$barcode)=@_;
	
	##
	my $cmd="perl $Bin/report.pl -conf $config -template $templateDir/pages.tt -outdir $dirout -id $barcode  -tables $dirout -contexts $dirout -mulu $templateDir/mulu.tt ";
	print "$cmd\n";
	`$cmd`;

	return;
}

#my $xml=&read_xml($xmlfile);
sub read_xml {
	my ($xmlfile)=@_;
	
	##
	my $xml=XMLin($xmlfile,NoAttr=>1,ForceArray => ['item','abnormalitem','disease','spot'], KeyAttr => []);
	
	return ($xml);
}

#&cp($templateDir,$referdir,$dirout);
sub cp {
	my ($templateDir,$referdir,$dirout)=@_;
	
	##
	if (-f "$templateDir/format.tex") {
		`cp $templateDir/format.tex $dirout`;
	}
	mkdir "$dirout/cores" if (! -d "$dirout/cores");
	if (-d "$templateDir/cores") {
		`cp $templateDir/cores/* $dirout/cores`;
	}
	`cp $referdir/png/* $dirout/cores`;
	
	return;
}

#&make_pdf($xml,"zongjie",$templateDir,$dirout,\@context1);
sub make_pdf {
	my ($xml,$key,$tempDir,$dirout,$context1)=@_;
	
	##
	my $tt="$templateDir/$key".".tt";
	my $ttTex=$key.".tex";
	if (-f $tt) {
		&section_tex($tt,$xml,"$dirout/$ttTex");
		chdir $dirout;
		mkdir "./PDF" if (! -d "./PDF");
		`xelatex $ttTex`;
	
		## get context config
		&push_contextconfig($ttTex,$context1);

		`mv $key.tex $key.pdf $key.log ./PDF`;
		`rm $key.out $key.aux`;
		print STDOUT "------------------\n$key.pdf\n";
		chdir $Bin;
	}
	
	return;
}

#&push_contextconfig("zongjie.tex",$context1);
sub push_contextconfig {
	my ($tex,$context)=@_;
	
	##
	my %config;
	my $basename=basename($tex);$basename=~s/\.tex/.pdf/;
	$config{key}=$basename;
	$config{'pageNum'}=&getPdfPages($tex);
	push @{$context},\%config;

	return;
}

#$config{'pageNum'}=&getPdfPages($tex);
sub getPdfPages {
	my ($tex)=@_;
	my $num;
	
	##
	$tex =~ s/\.tex$/.log/;
	if (-f $tex) {
		my ($found, $data) = (0, '');
		open IN, "<$tex" || die $!;
		while(<IN>){
			chomp;
			if (/^Output written on/) {
				$data = decode('utf8',$_);
				$found = 1;
			}
			elsif ($found == 1) {
				$data = decode('utf8',$_);
			}
			else {
			}
		}
		close IN;
		
		if ($data =~ /\((\d+) pages\)\.$/) {
			$num = $1;
		}
		elsif ($data =~ /\((\d+) page\)\.$/) {
			$num = $1;	
		}
		else {
		}
	}
	
	return($num);
}

sub makeMulu {
	my ($ttf, $datam, $diro) = @_;

	##
	if (-f $ttf) {
		my $texf = $ttf;
		$texf = basename($texf);
		$texf =~ s/\.tt$/\.tex/;
		my $bname = $texf;
		$bname =~ s/\.tex$//;
		$texf = $diro . '/' . $texf;
		&section_tex($ttf, $datam, $texf);
		chdir $diro;
		mkdir "./PDF" if (! -d "./PDF");
		`xelatex $bname.tex`;
		`mv $bname.pdf $bname.tex $bname.log ./PDF`;
		`rm -rf *.aux *.out`;
		print STDOUT "------------------\n$bname.pdf\n";
		chdir $Bin;
	}
	
	return;
}

#&section_tex($tt,\%vars,$tex);
sub section_tex {
	my ($section_tt,$section_d,$section_tex) = @_;
	$tpl->process($section_tt,$section_d,$section_tex,{binmode => ':encoding(UTF-8)'}) || die $tpl->error;
	return;
}

#&output_contexts("$dirout/context1",\@context1);
sub output_contexts {
	my ($outfile,$config)=@_;
	
	##
	open OUT,">$outfile" or die $!;
	binmode(OUT, ':encoding(utf-8)');
	print OUT "pdf\tcnt\ttitle\tsubclass\n";
	for (my $i=0;$i<@{$config};$i++) {
		print OUT "$config->[$i]->{'key'}\t$config->[$i]->{'pageNum'}\t\t\n";
	}
	close OUT;
	
	return;
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
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
    sprintf("%4d-%02d-%02d", $year+1900, $mon+1, $day);
}

sub USAGE {#
	my $usage=<<"USAGE";

    Program Function: generate GIhealth pdf report;
    Version:    $version
    Contact:    fred routine <fred_routine\@163.com> 
    Program Date:   2016.11.15
    Usage:
      Options:
      -xml         <file>  sample's xml format result                    forced
      -barcode     <str>   sample's code                                 forced
      -od          <dir>   output pdf report dir                         forced

      -template   <dir>   latex tt template dir                     optional
                          default "$Bin/template/current";
      -referdir   <dir>   reference pop data dir                    optional
                          default "$Bin/../../reference_data/current";
      -config     <file>  Full zunxiangban report config            optional
                          default "$Bin/config.all";
      -simconfig  <file>  Simple (5 pages) report config            optional
                          default "$Bin/config.simple";

      -h          Help

USAGE
	print $usage;
	exit;
}
