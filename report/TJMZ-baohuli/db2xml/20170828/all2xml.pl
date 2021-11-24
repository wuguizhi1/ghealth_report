#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use MongoDB;
use Encode;
use XML::Simple;
use Data::Dumper;
use MIME::Base64;
use JSON;

binmode(STDOUT, ':encoding(utf8)');

my ($host, $port, $productsdbname, $productscolname, $prodatadbname, $prodatacolname, $customer, $series, $sets, $language, $outdir);

GetOptions(
	"host=s"=> \$host,
	"port=i"=> \$port,
	"db1=s"	=> \$productsdbname,
        "cl1=s"	=> \$productscolname,
	"db2=s"	=> \$prodatadbname,
	"cl2=s"	=> \$prodatacolname,
        "lan=s"	=> \$language,
	"cus=s"	=> \$customer,
	"ser=s"	=> \$series,
	"set=s"	=> \$sets,
        "out=s"	=> \$outdir,
);

&usage() unless(defined $host && defined $port && defined $productsdbname && defined $productscolname && defined $prodatadbname && defined $prodatacolname && defined $language && defined $customer && defined $series && defined $sets && defined $outdir);

### products collection
my ($cli,$db,$col,$data,$rec);
$cli = MongoDB::MongoClient->new({'host' => $host . ':' . $port});
$db = $cli->get_database($productsdbname);
$col = $db->get_collection($productscolname);

$customer = decode('utf8', $customer);
$series = decode('utf8', $series);
$sets = decode('utf8', $sets);

my %list;

my ($taoxi,$taocan);
$data = $col->find({'belongsto' => $customer, $language . '.name' => $series});
while($rec = $data->next){
	my $last = 0;
	foreach my $key(keys %{$rec->{$language}->{'sets'}}){
		my $name = $rec->{$language}->{'sets'}->{$key}->{'name'};
		if($sets eq $name){
			foreach my $itm(keys %{$rec->{$language}->{'sets'}->{$key}->{'list'}}){
				$list{$rec->{$language}->{'sets'}->{$key}->{'list'}->{$itm}} = $itm;
				my @unit=split /\_/,$rec->{$language}->{'sets'}->{$key}->{'list'}->{$itm};
				$taoxi=join("_",$unit[0],$unit[1]);
				$taocan=$unit[2];
			}
			$last=1;
			last;
		}
	}
	last if($last eq 1);
}

### pagemode collection
my %pagemodes;
$db = $cli->get_database($prodatadbname);
$col = $db->get_collection("pagemodes");
$data = $col->find({});
while ($rec = $data->next) {
	$pagemodes{$rec->{_id}}{itms}=$rec->{itms};
	$pagemodes{$rec->{_id}}{version}=$rec->{version} if (exists $rec->{version});
}
#print Dumper %pagemodes;
#die;

### prodata collection
my %info;my %pic;
$db = $cli->get_database($prodatadbname);
$col = $db->get_collection($prodatacolname);
my $key=join("_",$taoxi,$taocan);
#$data = $col->find({});
$data = $col->query({'_id' => qr/$key/});
while ($rec = $data->next) {
	my $id=$rec->{'_id'};
	next if (! defined $id || ! exists $list{$id});
	my $pagemode=$rec->{pagemode};
	my $categoey=$rec->{category};
	$info{$categoey}{$id}{id}=$id;
	$info{$categoey}{$id}{categoey}=$categoey;
	$info{$categoey}{$id}{oriid}=$rec->{oriid};
	$info{$categoey}{$id}{orititle}=$rec->{$language}->{title};
	$info{$categoey}{$id}{reportNum}=$rec->{reportNum} if (defined $rec->{reportNum});
	foreach my $itm (@{$pagemodes{$pagemode}{itms}}) {
		foreach my $k (keys %{$itm}) {
			next if ($k eq "sex");
			if ($k eq "pic") {
				#$pic{$categoey}{$id}{base64}=decode_base64($rec->{$language}->{$k}->{base64});
				#$pic{$categoey}{$id}{type}=$rec->{$language}->{$k}->{mimetype};
				#$pic{$categoey}{$id}{type}=~s/image\///;
				#$info{$categoey}{$id}{$language}{pic}="section_".$info{$categoey}{$id}{oriid}.".".$pic{$categoey}{$id}{type};
			}
			else {
				if ($k eq "title"  || $k eq "etitle" || $k eq "category" || $k eq "order") {
					if (defined $rec->{$language}->{$k} && $rec->{$language}->{$k} ne "") {
						$info{$categoey}{$id}{$language}{$k}=$rec->{$language}->{$k};
						$info{$categoey}{$id}{$language}{$k}=&replace($info{$categoey}{$id}{$language}{$k});
					}
					else {
						$info{$categoey}{$id}{$language}{$k}="";
					}
				}
				elsif ($k eq "subclass" || $k eq "section1" || $k eq "section2") {
					if (defined $rec->{$language}->{$k}->{name} && $rec->{$language}->{$k}->{name} ne "") {
						$info{$categoey}{$id}{$language}{$k}{name}=$rec->{$language}->{$k}->{name};
						$info{$categoey}{$id}{$language}{$k}{name}=&replace($info{$categoey}{$id}{$language}{$k}{name});
					}
					else {
						$info{$categoey}{$id}{$language}{$k}{name}="";
					}
					$info{$categoey}{$id}{$language}{$k}{order}="";
					$info{$categoey}{$id}{$language}{$k}{order}=$rec->{$language}->{$k}->{order} if (defined $rec->{$language}->{$k}->{order} && $rec->{$language}->{$k}->{order} ne "");
				}
				elsif (($k eq "suggestion" || $k eq "metafood" || $k eq "metaregulate" || $k eq "metasport" || $k eq "agesuggestion" || $k eq "dieasesuggestion" || $k eq "microsuggestion" || $k eq "food" || $k eq "regulate" || $k eq "sport") && exists $pagemodes{$pagemode}{version}) {
					foreach my $unit (@{$itm->{$k}}) {
						foreach my $u (keys %{$unit}) {
							next if ($u eq "node" || $u eq "type");
							if ($u eq "header") {
								if (defined $rec->{$language}->{$k}->{$u} && $rec->{$language}->{$k}->{$u} ne "") {
									$info{$categoey}{$id}{$language}{$k}{$u}=$rec->{$language}->{$k}->{$u};
									$info{$categoey}{$id}{$language}{$k}{$u}=&replace($info{$categoey}{$id}{$language}{$k}{$u});
								}
								else {
									$info{$categoey}{$id}{$language}{$k}{$u}="";
								}
							}
							else {
								foreach my $level (@{$unit->{$u}}) {
									foreach my $l (keys %{$level}) {
										next if ($l eq "node" || $l eq "type");
										foreach my $part (@{$level->{$l}}) {
											foreach my $p (keys %{$part}) {
												next if ($p eq "node" || $p eq "type");
												if (defined $rec->{$language}->{$k}->{$u}->{$l}->{$p} && $rec->{$language}->{$k}->{$u}->{$l}->{$p} ne "") {
													$info{$categoey}{$id}{$language}{$k}{$u}{$l}{$p}=$rec->{$language}->{$k}->{$u}->{$l}->{$p};
													$info{$categoey}{$id}{$language}{$k}{$u}{$l}{$p}=&replace($info{$categoey}{$id}{$language}{$k}{$u}{$l}{$p});
												}
												else {
													$info{$categoey}{$id}{$language}{$k}{$u}{$l}{$p}="";
												}
											}
										}
									}
								}
							}
						}
					}	
				}
				else {
					foreach my $unit (@{$itm->{$k}}) {
						foreach my $u (keys %{$unit}) {
							next if ($u eq "node" || $u eq "type");
							if (defined $rec->{$language}->{$k}->{$u} && $rec->{$language}->{$k}->{$u} ne "") {
								$info{$categoey}{$id}{$language}{$k}{$u}=$rec->{$language}->{$k}->{$u};
								$info{$categoey}{$id}{$language}{$k}{$u}=&replace($info{$categoey}{$id}{$language}{$k}{$u});
							}
							else {
								$info{$categoey}{$id}{$language}{$k}{$u}="";
							}
						}
					}
				}
			}
		}
	}
}
#print Dumper %info;
#die;

### output xml
my $json = JSON->new;
foreach my $type (keys %info) {
	mkdir $outdir if (! -d $outdir);
	my $dir="$outdir/$type";
	mkdir $dir if (! -d $dir);
	mkdir "$dir/$language" if (! -d "$dir/$language");
#	mkdir "$dir/$language/section" if (! -d "$dir/$language/section");
#	mkdir "$dir/$language/pic" if (! -d "$dir/$language/pic");
	foreach my $id (keys %{$info{$type}}) {
		my $oriid=$info{$type}{$id}{oriid};
		$oriid=~s/\'//g;
		$oriid=~s/\`//g;
		$oriid=~s/ //g;
		my $outfile="$dir/$language/$oriid".".xml";
#		my $outpic="$dir/$language/pic/section_"."$oriid"."."."$pic{$type}{$id}{type}";
		open my $out, '>:encoding(utf-8)', $outfile || die $!;
		XMLout($info{$type}{$id},  OutputFile => $out, NoAttr => 1, SuppressEmpty => "", KeyAttr => []);
		close $out;
#		open OUT,">$outpic" || die $!;
#		print OUT "$pic{$type}{$id}{base64}\n";
#		close OUT;
	}
}

### replace
sub replace {
        my ($js)=@_;
        my $utfcode = '';
        ##
        $js=~s/{/\\{/g;
        $js=~s/}/\\}/g;
        $js=~s/%/{\\%}/g;
        $js=~s/\n/\n\n/g;
        $utfcode = decode('utf-8', 'Ⅰ');
        $js=~s/$utfcode/\\RNum{1}/g;
        $utfcode = decode('utf-8', 'Ⅱ');
        $js=~s/$utfcode/\\RNum{2}/g;
        $utfcode = decode('utf-8', 'Ⅲ');
        $js=~s/$utfcode/\\RNum{3}/g;
        $utfcode = decode('utf-8', 'Ⅳ');
        $js=~s/$utfcode/\\RNum{4}/g;
        $utfcode = decode('utf-8', 'Ⅴ');
        $js=~s/$utfcode/\\RNum{5}/g;
        $utfcode = decode('utf-8', 'Ⅵ');
        $js=~s/$utfcode/\\RNum{6}/g;
        $utfcode = decode('utf-8', 'Ⅶ');
        $js=~s/$utfcode/\\RNum{7}/g;
        $utfcode = decode('utf-8', 'Ⅷ');
        $js=~s/$utfcode/\\RNum{8}/g;
        $utfcode = decode('utf-8', 'Ⅸ');
        $js=~s/$utfcode/\\RNum{9}/g;
        $utfcode = decode('utf-8', 'Ⅹ');
        $js=~s/$utfcode/\\RNum{10}/g;
        $utfcode = decode('utf-8', 'Ⅺ');
        $js=~s/$utfcode/\\RNum{11}/g;
        $utfcode = decode('utf-8', 'Ⅻ');
        $js=~s/$utfcode/\\RNum{12}/g;
        $utfcode = decode('utf-8', 'ⅩⅢ');
        $js=~s/$utfcode/\\RNum{13}/g;
        $utfcode = decode('utf-8', 'α');
        $js=~s/$utfcode/\\textalpha /g;
        $utfcode = decode('utf-8', 'β');
        $js=~s/$utfcode/\\textbeta /g;
        $utfcode = decode('utf-8', 'γ');
        $js=~s/$utfcode/\\textgamma /g;
        $utfcode = decode('utf-8', 'μ');
        $js=~s/$utfcode/\\textmu /g;
        $utfcode = decode('utf-8', 'δ');
        $js=~s/$utfcode/\\textdelta /g;
        $utfcode = decode('utf-8', 'κ');
        $js=~s/$utfcode/\\textkappa /g;
        $utfcode = decode('utf-8', 'ε');
        $js=~s/$utfcode/\\textepsilon /g;
        $utfcode = decode('utf-8', 'ƹ');
        $js=~s/$utfcode/\\textepsilon /g;
        $js=~s/\$/\\\$/g;
        $utfcode = decode('utf-8', '≥');
        $js=~s/$utfcode/\$\\geq\$/g;
        $utfcode = decode('utf-8', '≤');
        $js=~s/$utfcode/\$\\leq\$/g;
        $js=~s/~/\\textasciitilde /g;
        $js=~s/_/\\_/g;
        $js=~s/#/\\#/g;
        $js=~s/\//\{\/\}/g;

        while ($js=~/\^(\D+)/) {
                $js=~s/\^/\\\^{}/g;
        }
        while ($js=~/\^(\d+)/) {
                my $ori=$1;
                $js=~s/\^$ori/\$\^\{$ori\}\$/g;
        }

        $utfcode = decode('utf-8', '∮');
        while ($js=~/$utfcode(\d+)/) {
                my $ori=$1;
                $js=~s/$utfcode$ori/\$\_\{$ori\}\$/g;
        }

        $utfcode = decode('utf-8', '【');
        $js=~s/【/\\vskip.6\\baselineskip\\noindent {\\bfseries\\wuhao 【/g;
        $utfcode = decode('utf-8', '】');
        $js=~s/】/】}\\smallskip/g;
        $utfcode = decode('utf-8', '腘');
        $js=~s/腘/\\mbox{\\scalebox{0.5}[1]{月 }\\kern-.15em\\scalebox{0.75}[1]{国}}/g;

        return ($js);
}
=pod
sub replace {
	my ($js)=@_;
	my $utfcode = '';
	##
	$js=~s/{/\\{/g;
	$js=~s/}/\\}/g;
	$js=~s/%/{\\%}/g;
	$js=~s/\n/\n\n/g;
	$utfcode = decode('utf-8', 'Ⅰ');
	$js=~s/$utfcode/\\RNum{1}/g;
	$utfcode = decode('utf-8', 'Ⅱ');
	$js=~s/$utfcode/\\RNum{2}/g;
	$utfcode = decode('utf-8', 'Ⅲ');
	$js=~s/$utfcode/\\RNum{3}/g;
	$utfcode = decode('utf-8', 'Ⅳ');
	$js=~s/$utfcode/\\RNum{4}/g;
	$utfcode = decode('utf-8', 'Ⅴ');
	$js=~s/$utfcode/\\RNum{5}/g;
	$utfcode = decode('utf-8', 'Ⅵ');
	$js=~s/$utfcode/\\RNum{6}/g;
	$utfcode = decode('utf-8', 'Ⅶ');
	$js=~s/$utfcode/\\RNum{7}/g;
	$utfcode = decode('utf-8', 'Ⅷ');
	$js=~s/$utfcode/\\RNum{8}/g;
	$utfcode = decode('utf-8', 'Ⅸ');
	$js=~s/$utfcode/\\RNum{9}/g;
	$utfcode = decode('utf-8', 'Ⅹ');
	$js=~s/$utfcode/\\RNum{10}/g;
	$utfcode = decode('utf-8', 'Ⅺ');
	$js=~s/$utfcode/\\RNum{11}/g;
	$utfcode = decode('utf-8', 'Ⅻ');
	$js=~s/$utfcode/\\RNum{12}/g;
	$utfcode = decode('utf-8', 'ⅩⅢ');
	$js=~s/$utfcode/\\RNum{13}/g;
	$utfcode = decode('utf-8', 'α');
	$js=~s/$utfcode/\\textalpha /g;
	$utfcode = decode('utf-8', 'β');
	$js=~s/$utfcode/\\textbeta /g;
	$utfcode = decode('utf-8', 'γ');
	$js=~s/$utfcode/\\textgamma /g;
	$utfcode = decode('utf-8', 'μ');
	$js=~s/$utfcode/\\textmu /g;
	$utfcode = decode('utf-8', 'δ');
	$js=~s/$utfcode/\\textdelta /g;
	$utfcode = decode('utf-8', 'κ');
	$js=~s/$utfcode/\\textkappa /g;
	$utfcode = decode('utf-8', 'ε');
	$js=~s/$utfcode/\\textepsilon /g;
	$utfcode = decode('utf-8', 'ƹ');
	$js=~s/$utfcode/\\textepsilon /g;
	$js=~s/\$/\\\$/g;
	$utfcode = decode('utf-8', '≥');
	$js=~s/$utfcode/\$\\geq\$/g;
	$utfcode = decode('utf-8', '≤');
	$js=~s/$utfcode/\$\\leq\$/g;
	$js=~s/~/\\textasciitilde /g;
	$js=~s/_/\\_/g;
	$js=~s/#/\\#/g;
	$js=~s/\//\{\/\}/g;
	
	while ($js=~/\^(\D+)/) {
		$js=~s/\^/\\\^{}/g;
	}
	while ($js=~/\^(\d+)/) {
		my $ori=$1;
		$js=~s/\^$ori/\$\^\{$ori\}\$/g;
	}

	$utfcode = decode('utf-8', '∮');
	while ($js=~/$utfcode(\d+)/) {
		my $ori=$1;
		$js=~s/$utfcode$ori/\$\_\{$ori\}\$/g;
	}
	
	$utfcode = decode('utf-8', '【');
	$js=~s/【/\\vskip.6\\baselineskip\\noindent {\\bfseries\\wuhao 【/g;
	$utfcode = decode('utf-8', '】');
	$js=~s/】/】}\\smallskip/g;
	$utfcode = decode('utf-8', '腘');
	$js=~s/腘/\\mbox{\\scalebox{0.5}[1]{月 }\\kern-.15em\\scalebox{0.75}[1]{国}}/g;

	return ($js);
}
=cut

###############################################################################
sub usage{
	print STDERR <<P
perl $0 [options]
-host	host
-port	port
-db1	products database name
-cl1	products collection name
-db2	prodata database name
-cl2	prodata collection name
-lan	which language to export
-cus	customer name
-ser	series name
-set	set name
-out	outdir
P
}
