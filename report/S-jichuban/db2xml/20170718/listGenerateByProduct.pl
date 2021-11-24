#!/usr/bin/perl -w
use strict;
use MongoDB;
use Encode;
use Getopt::Long;
use Data::Dumper;

binmode(STDOUT, ':encoding(utf8)');

my ($customer, $series, $sets, $lan, $outdir, $dbname, $colname, $db, $col, $data, $rec, $cli);

GetOptions(
	"cus=s"	=> \$customer,
	"ser=s"	=> \$series,
	"set=s"	=> \$sets,
	"lan=s"	=> \$lan,
	"out=s"	=> \$outdir,
);

&usage() unless(defined $customer && defined $series && defined $sets && defined $lan && defined $outdir);

$customer = decode('utf8', $customer);
$series = decode('utf8', $series);
$sets = decode('utf8', $sets);
$lan = (split(/,/, $lan))[0];

$cli = MongoDB::MongoClient->new('host' => '10.0.0.204:27021');

$dbname = 'susceptibility';
$colname = 'products';
$db = $cli->get_database($dbname);
$col = $db->get_collection($colname);
my @region;
my %list;
$data = $col->find({'belongsto' => $customer, $lan . '.name' => $series});
while($rec = $data->next){
	my $last = 0;
	foreach my $key(keys %{$rec->{$lan}->{'sets'}}){
		my $name = $rec->{$lan}->{'sets'}->{$key}->{'name'};
		if($sets eq $name){
			@region = @{$rec->{$lan}->{'sets'}->{$key}->{'region'}};
			foreach my $itm(keys %{$rec->{$lan}->{'sets'}->{$key}->{'list'}}){
				$list{$rec->{$lan}->{'sets'}->{$key}->{'list'}->{$itm}} = $itm;
			}
			$last = 1;
			last;
		}
	}
	last if($last eq 1);
}

$dbname = 'susceptibility';
$colname = 'prodata';
$db = $cli->get_database($dbname);
$col = $db->get_collection($colname);
my @itms = keys %list;
&checkdir($outdir);

foreach my $cat(@region){
#	$data = $col->find({'_id' => {'$in' => \@itms}, 'category' => $cat})->sort({$lan . '.subclass.order' => 1, $lan . '.order' => 1});
#	$data = $col->find({'_id' => {'$in' => \@itms}, 'category' => $cat});#->sort({$lan . '.section1.order' => 1, $lan . '.section2.order' => 1, $lan . '.order' => 1});
#	$data = $col->query({'_id' => {'$in' => \@itms}, 'category' => $cat})->fields({'oriid' => 1, 'sex' => 1, 'CN.title' => 1, 'CN.section1' => 1, 'CN.section2' => 1, 'CN.order' => 1});#->sort({$lan . '.section1.order' => 1, $lan . '.section2.order' => 1, $lan . '.order' => 1});
	$data = $col->query({'_id' => {'$in' => \@itms}, 'category' => $cat})->fields({'oriid' => 1, 'sex' => 1, $lan . '.title' => 1, $lan . '.etitle' => 1, $lan . '.section1' => 1, $lan . '.section2' => 1, $lan . '.order' => 1});#->sort({$lan . '.section1.order' => 1, $lan . '.section2.order' => 1, $lan . '.order' => 1});
	my %hash;
	my %wocao;
	while($rec = $data->next){
		my $sex = 'B';
		$sex = $rec->{'sex'} if(defined $rec->{'sex'});
		my @genders;
		if($sex eq 'B'){
			@genders = ('B', 'M', 'F');
		}else{
			@genders = ('B', $sex);
		}
		if(!defined $rec->{$lan}->{'section2'}){
			$rec->{$lan}->{'section2'}->{'name'} = '';
			$rec->{$lan}->{'section2'}->{'order'} = 1;
		}
		foreach my $gender(@genders){
			#push(@{$wocao{$gender}{$rec->{$lan}->{'section1'}->{'order'}}{'subs'}{$rec->{$lan}->{'section2'}->{'order'}}{'subs'}}, $rec->{'oriid'});#$rec->{$lan}->{'title'});
			$rec->{$lan}->{'section2'}->{'order'} = int($rec->{$lan}->{'section2'}->{'order'});
			$rec->{$lan}->{'section1'}->{'order'} = int($rec->{$lan}->{'section1'}->{'order'});
			$rec->{$lan}->{'order'} = int($rec->{$lan}->{'order'});
			$wocao{$gender}{$rec->{$lan}->{'section1'}->{'order'}}{'name'} = $rec->{$lan}->{'section1'}->{'name'};
			$wocao{$gender}{$rec->{$lan}->{'section1'}->{'order'}}{'subs'}{$rec->{$lan}->{'section2'}->{'order'}}{'name'} = $rec->{$lan}->{'section2'}->{'name'};
			$wocao{$gender}{$rec->{$lan}->{'section1'}->{'order'}}{'subs'}{$rec->{$lan}->{'section2'}->{'order'}}{'subs'}{$rec->{$lan}->{'order'}} = $rec->{'oriid'} . "\t" ;
			if (exists $rec->{$lan}->{'etitle'} && $rec->{$lan}->{'etitle'} ne "") {
				$wocao{$gender}{$rec->{$lan}->{'section1'}->{'order'}}{'subs'}{$rec->{$lan}->{'section2'}->{'order'}}{'subs'}{$rec->{$lan}->{'order'}} .= $rec->{$lan}->{'etitle'} . "\t";
			}
			$wocao{$gender}{$rec->{$lan}->{'section1'}->{'order'}}{'subs'}{$rec->{$lan}->{'section2'}->{'order'}}{'subs'}{$rec->{$lan}->{'order'}} .= $rec->{$lan}->{'title'};

		}
	}
	&printOneCatInArray(\%wocao, $outdir, $cat);
}

sub printOneCat{
	my ($onerec, $dir, $cat) = @_;
	foreach my $sex(keys %{$onerec}){
		#print $sex, '-' x 5, "\n";
		my $tag = $sex;
		$tag = 'A' if($tag eq 'B');
		my $file = $dir . '/' . $cat . '.list.order' . $tag;
		open OUT, ">:encoding(utf8)", $file || die $!;
		foreach my $order(sort{$a <=> $b} keys %{$onerec->{$sex}}){
			print OUT "#$onerec->{$sex}->{$order}->{'name'}\n";
			foreach my $suborder(sort{$a <=> $b} keys %{$onerec->{$sex}->{$order}->{'order'}}){
				print OUT "section_$onerec->{$sex}->{$order}->{'order'}->{$suborder}\n";
			}
		}
		close OUT;
	}
}

sub printOneCatInArray{
	my ($all, $dir, $cat) = @_;
	foreach my $sex(keys %{$all}){
		my $tag = $sex;
		$tag = 'A' if($tag eq 'B');
		my $file = $dir . '/' . $cat . '.list.order' . $tag;
		open OUT, ">:encoding(utf8)", $file || die $!;
		foreach my $order(sort{$a <=> $b} keys %{$all->{$sex}}){
			print OUT "###$all->{$sex}->{$order}->{'name'}\n";
			foreach my $suborder(sort{$a <=> $b} keys %{$all->{$sex}->{$order}->{'subs'}}){
				print OUT "##";
				print OUT "$all->{$sex}->{$order}->{'subs'}->{$suborder}->{'name'}" if(defined $all->{$sex}->{$order}->{'subs'}->{$suborder}->{'name'});
				print OUT "\n";
				foreach my $forder(sort{$a <=> $b} keys %{$all->{$sex}->{$order}->{'subs'}->{$suborder}->{'subs'}}){
					print OUT $all->{$sex}->{$order}->{'subs'}->{$suborder}->{'subs'}->{$forder}, "\n";
				}
			}
		}
		close OUT;
	}
}

sub checkdir{
	my ($dir) = @_;
	`mkdir -p $dir` unless(-d $dir);
}

sub usage{
	print STDERR <<P
perl $0 [options]
-cus customer name
-ser series name
-lan lan (CN|EN...)
-set set name
-out outdir
P
;
exit(1);
}
