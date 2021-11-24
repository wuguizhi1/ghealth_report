#!/usr/bin/perl -w
use strict;
use Template;
use Getopt::Long;
use Data::Dumper;
use Encode;
use Cwd 'abs_path';
use FindBin qw($Bin);
use File::Basename qw(basename dirname);

use lib '/data/bioit/biodata/chuanj/lib/';
use Genetalks::Utils;

binmode(STDOUT, ':encoding(utf-8)');
binmode(STDERR, ':encoding(utf-8)');

my ($conf, $template, $outdir, $help, %para, %data, $mulu);
GetOptions(
	"conf:s"	=> \$para{'conf'},
	"template:s"	=> \$para{'template'},
	"barcode:s"	=> \$para{'barcode'},
	"id:s"		=> \$para{'id'},
	"outdir:s"	=> \$para{'outdir'},
	"tables:s"	=> \$para{'tables'},
	"contexts:s"	=> \$para{'contexts'},
	"needStyle:s"	=> \$para{'needStyle'},
	"isTitle:s"	=> \$para{'isTitle'},
	"displayPageNo:s"	=> \$para{'displayPageNo'},
	"displayLogo:s"	=> \$para{'displayLogo'},
	"displayBarcode:s"	=> \$para{'displayBarcode'},
	"mulu:s"	=> \$mulu,
	"check"	=> \$para{'check'},
	"help"	=> \$help,
);

die `pod2text $0` if($help);
die `pod2text $0` unless(defined $para{'conf'} && (defined $para{'check'} || (defined $para{'template'} && defined $para{'outdir'} && defined $para{'tables'} && defined $para{'contexts'})));

&absPath(\%para, 'conf');
$Bin = abs_path($Bin);

&getData($para{'conf'}, \%data, 'f', \%para);
push(@{$data{'pages'}}, $data{'last'});
delete $data{'last'};

#print Dumper \%data;

exit(0) if(defined $para{'check'});

Genetalks::Utils->checkDir($para{'outdir'});
&absPath(\%para, 'template', 'outdir', 'tables', 'contexts');
my $driver = Template->new({ENCODING => 'utf8', ABSOLUTE => 1});
#`cp -r templates/ $para{'outdir'}`;
#`cp -r pics/ $para{'outdir'}`;

if(defined $mulu){#added for meta's mulu only
	my $pgn = &getPageNo(\%data);
#print Dumper $pgn;
	&generMulu($Bin, $para{'outdir'}, $driver, $pgn, $mulu);
}

my $i = 0;
my @prefix;
foreach my $part(@{$data{'pages'}}){
	#$driver->process($para{'template'}, \%data, $para{'out'}, {binmode => ':encoding(UTF-8)'});
	$part->{'barcode'}=$para{'barcode'};
	$driver->process($para{'template'}, $part, $para{'outdir'} . '/out' . $i . '.tex', {binmode => ':encoding(UTF-8)'});
	push(@prefix, $para{'outdir'} . '/out' . $i);
	$i++;
}

my $pdfs;
chdir $para{'outdir'};
foreach my $pre(@prefix){
	`xelatex $pre . '.tex'`;
	$pdfs .= " $pre.pdf ";
}

`pdftk $pdfs cat output $para{'id'}.pdf`;

sub getPageNo{
	my ($d) = @_;
	my %pageNo = ();
	foreach my $pages(@{$d->{'pages'}}){
		foreach my $page(@{$pages->{'pages'}}){
			my $key = $page->{'pdf'};
			$key =~ s/\.pdf$//;
			$pageNo{$key} = $page->{'style'}->{'pageNo'};
		}
	}
	return \%pageNo;
}

sub generMulu{
	my ($curd, $outd, $dri, $d, $tt) = @_;
	chdir $outd;
	$dri->process($tt, $d, "$outd/mulu.tex", {binmode => ':encoding(UTF-8)'});
	`xelatex mulu.tex`;
	chdir $curd;
}


#conf %data type %para
sub getData{
	my ($f, $d, $t, $p, $pgStyle) = @_;
	open my $in, "<$f" || die $!;
	my $head;
	my @header;
	while($head = <$in>){
		chomp($head);
		next if($head =~ /^\s*$/ || $head =~ /^#/);
		@header = split(/\t/, $head);
		last;
	}
	while(my $line = <$in>){
		chomp($line);
		next if($line =~ /^\s*$/ || $line =~ /^#/);
		my @tmp = split(/\t/, $line);
		my %page;
		if(defined $pgStyle){
			foreach my $kx(keys %{$pgStyle->{'style'}}){
				$page{'style'}{$kx} = $pgStyle->{'style'}->{$kx};
			}
			$pgStyle->{'style'}->{'isFirstCounter'} = '';
			$pgStyle->{'style'}->{'isRight'} = '';
		}
		for(my $i = 0; $i <= $#tmp; $i++){
			if($i <= 1){
				$page{$header[$i]} = decode('utf-8', $tmp[$i]);
			}else{
				$page{'style'}{$header[$i]} = decode('utf-8', $tmp[$i]);
			}
		}
#		if($tmp[0] eq 'tables' || $tmp[0] eq 'contexts'){
#			&getData($p->{$tmp[0]}, $d, $tmp[0], $p, \%page);
		if($tmp[0] =~ /^(tables)/ || $tmp[0] =~ /^(contexts)/){
#			&getData("/data/bioit/biodata/chuanj/latex/meta/report/rep/test/$tmp[0]", $d, $tmp[0], $p, \%page);
			&getData("$p->{$1}/$tmp[0]", $d, $tmp[0], $p, \%page);
		}else{
			$page{'type'} = $t;
			#总页数增加
			if(defined $d->{'pageCnt'}){
				$d->{'pageCnt'} += $d->{'lastCnt'};
			}else{
				$d->{'pageCnt'} = 1;;
			}
			$page{'style'}{'pageCnt'} = $d->{'pageCnt'};
			&correctPara($p, \%page);#参数重置配置文件

			if(defined $page{'style'}{'isFirstCounter'} && $page{'style'}{'isFirstCounter'} eq 'Y'){#页码重计
				$d->{'pageNo'} = 1;
				$d->{'needBack'} = 0;
			}else{
				#页码增加
				if(defined $d->{'pageNo'}){
					$d->{'pageNo'} += $d->{'lastCnt'};
				}else{
					$d->{'pageNo'} = 1;
				}
				if(defined $d->{'lastIsPageNo'} && $d->{'lastIsPageNo'} eq 'N' && $page{'style'}{'isPageNo'} eq 'Y'){#前面是不计页码的 页码回退 不计的页数清零
					$d->{'pageNo'} -= $d->{'needBack'};
					$d->{'needBack'} = 0;
				}
			}
			$page{'style'}{'pageNo'} = $d->{'pageNo'};
			$d->{'lastCnt'} = $page{'cnt'};#用来计算当前页的页数起始值
			$d->{'lastIsPageNo'} = $page{'style'}{'isPageNo'};#用来判断页码是否需要回退
			$d->{'needBack'} += $page{'cnt'} if(defined $page{'style'}{'isPageNo'} && $page{'style'}{'isPageNo'} eq 'N');#累计页码回退
			my $isOdd = &checkLeftRight(\%page, $p->{'check'});
			if(defined $d->{'last'}->{'isOdd'}){
				if($d->{'last'}->{'isOdd'} != $isOdd){#前后两部分 左右留白距离不一致时拆分
					push(@{$d->{'pages'}}, $d->{'last'});
					$d->{'last'} = ();
					print STDERR "-----new part latex\n";
				}
			}
			$d->{'last'}->{'isOdd'} = $isOdd;
#print Dumper \%page;
			push(@{$d->{'last'}->{'pages'}}, \%page);
		}
	}
	close $in;
}

sub checkLeftRight{
	my ($p, $check) = @_;
	my $cntOdd = $p->{'style'}->{'pageCnt'} % 2;
	my $curOdd = $p->{'style'}->{'pageNo'} % 2;
	my $judge = 0;

	if($cntOdd == $curOdd){
		$judge = 1;
	}else{
		$judge = 0;
	}
	if(defined $check){
		print STDERR "[$p->{'pdf'}]\t$p->{'style'}->{'pageNo'}:$curOdd $p->{'style'}->{'pageCnt'}:$cntOdd\n";
		if(defined $p->{'style'}{'isRight'} && $p->{'style'}{'isRight'} eq 'Y'){
			if($cntOdd == 0){
				print STDERR "[$p->{'pdf'}]\twill not display on the right side\n";
			}elsif($curOdd == 0){
				print STDERR "[$p->{'pdf'}]\twill display on the right side, But the page number is: [$p->{'style'}->{'pageNo'}]\n";
			}else{
				print STDERR "[$p->{'pdf'}]\tOK\n";
			}
		}elsif(defined $p->{'style'}{'isRight'} && $p->{'style'}{'isRight'} eq 'N'){
			if($cntOdd == 1){
				print STDERR "[$p->{'pdf'}]\twill not display on the left side\n";
			}elsif($curOdd == 1){
				print STDERR "[$p->{'pdf'}]\twill display on the left side, But the page number is: [$p->{'style'}->{'pageNo'}]\n";
			}else{
				print STDERR "[$p->{'pdf'}]\tOK\n";
			}
		}
	}
	return $judge;
}

sub correctPara{
	my ($p, $d) = @_;
	foreach my $k('displayLogo', 'displayBarcode', 'displayPageNo'){
		next if(defined $d->{'style'}->{$k} && $d->{'style'}->{$k} eq 'N');
		if(defined $p->{$k}){
			$d->{'style'}->{$k} = $p->{$k};
		}
	}
#	if($d->{'type'} eq 'tables' || $d->{'type'} eq 'contexts'){
	if($d->{'type'} =~ /tables/ || $d->{'type'} =~ /contexts/){
		foreach my $k('displayLogo', 'displayBarcode', 'displayPageNo', 'needStyle', 'isCore', 'isPageNo', 'isTitle'){
			if(defined $p->{$k}){
				$d->{'style'}->{$k} = $p->{$k};
			}else{
				$d->{'style'}->{$k} = 'Y';
			}
		}
#		$d->{'style'}->{'isTitle'} = 'N';
		if(defined $d->{'style'}->{'type'} && $d->{'style'}->{'type'} eq 'ads'){
			$d->{'style'}->{'needStyle'} = 'N';
			$d->{'style'}->{'isCore'} = 'N';
			$d->{'style'}->{'displayPageNo'} = 'N';
			$d->{'style'}->{'displayLogo'} = 'N';
			$d->{'style'}->{'isTitle'} = 'N';
			$d->{'style'}->{'displayBarcode'} = 'N';
			$d->{'style'}->{'isPageNo'} = 'Y';
		}
		$d->{'style'}->{'style'} = $d->{'type'} . '-' . $d->{'style'}->{'pageCnt'};
	}
}

sub absPath{
	my ($p, @ks) = @_;
	foreach my $k(@ks){
		$p->{$k} = abs_path($p->{$k});
	}
}



=head1 Name
	report.pl
Author
	Jun Chuan
Options
	perl report.pl [options]
	-conf	configuration file
	-template	template file
	-outdir	output path
	-id	id
	-barcode	barcode
	-tables	tables file dir
	-contexts	context files dir
	-mulu	mulu.tt file, for meta only
	-check	check the order in the configuration file, and exit with printing information out to STDERR

	style parameters to overwrite the style:
	-needStyle	[Y|N => Y]
	or groups:
	-isTitle	[Y|N => N]
	-displayBarcode	[Y|N => Y]
	-displayLogo	[Y|N => Y]
	-displayPageNo	[Y|N => Y]
Usage
	For run:
		perl report.pl -conf conf -template pages.tt -id id -barcode barcode-outdir ./output/ -tables tablesDir/ -contexts contextsDir/
	For prepare:
		perl report.pl -conf conf -check
=cut
