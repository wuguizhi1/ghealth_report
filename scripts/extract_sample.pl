#!/usr/bin/perl -w
#########################################################################
# Author: mengfei
# mail: fred_routine@163.com
# Created Time: Tue May 12 09:57:32 CST 2015
#########################################################################

use strict;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Getopt::Long;
use Data::Dumper;
use JSON::XS;
use XML::Simple qw(XMLin XMLout);

my $path=`pwd`;
chomp($path);

my %opts;
GetOptions(\%opts,"json=s","prefix=s","o=s","d=s","h");

my $usage = <<"USAGE";
        Program : $0
        Contact : mengfei   fred_routine\@163.com
        Usage   : $0 [options]
        Option  :
                  -json     :: results json file, required;
                  -prefix   :: output files prefix, required;

                  -o    :: output dir default "./";

              conf:
                  conf's: config files dir,
                          conf's file directory, default "$Bin/../reference_data";
                   -d   :: directory of config files;

USAGE

if(!defined($opts{"json"}) || !defined($opts{"prefix"}) || defined($opts{"h"}))
{
	print $usage;
	exit;
}

my $result_json=&ABSOLUTE_DIR($opts{"json"});
my $prefix=$opts{"prefix"};
my $out_dir=$opts{"o"} || "./";
&MKDIR($out_dir);
$out_dir=&ABSOLUTE_DIR($out_dir);

my $conf_d=$opts{"d"} || "$Bin/../reference_data";
$conf_d=&ABSOLUTE_DIR($conf_d);

my $term_conf="$conf_d/dbconf/term.conf";
my $term2bac_conf="$conf_d/dbconf/term2bac.conf";
my $pathogen_conf="$conf_d/dbconf/pathogen.conf";
my $factor2bac_conf="$conf_d/dbconf/disease_factor2bac.conf";
my $diseaserisk_conf="$conf_d/dbconf/disease_risk2bac.conf";
my $mmc_conf="$conf_d/dbconf/mmc_factor2bac.conf";

#############################
# load conf files for result extraction
my %term;
&LOAD_TERM ($term_conf, \%term);

my %term2bac;
&LOAD_BAC ($term2bac_conf, \%term2bac);

my %pathogen;
&LOAD_BAC ($pathogen_conf, \%pathogen);

my %disease;
&LOAD_F2B ($factor2bac_conf, \%disease);
&LOAD_F2B ($mmc_conf, \%disease);

my %disease_risk;
&LOAD_DR ($diseaserisk_conf, \%disease_risk);

# load microbiota and function result json
open IN,"$result_json" || die $!;
my $result_str;
while (<IN>) {
	chomp;
	s/^\s+|\s+$//;
	next if (/^$/);
	$result_str .= $_;
}
close IN;
my $refer = decode_json $result_str;
my @samples = sort keys %{$refer->{"summary"}{"diversity"}};

#### get microbiota, pathway, ko results
my %microbiota_result; my %pathway_result; my %ko_result;
&REFER_SPLIT ($refer, \%term, \%pathogen, \%microbiota_result, \%pathway_result, \%ko_result);

#### extract bacteria or function to terms
my %term2bac_result;
&EXTRACT_T2B (\%term2bac, $refer, \%term2bac_result);

my %diseasefactor_result;
&EXTRACT_F2B (\%disease, $refer, \%diseasefactor_result);

my %diseaserisk_result;
&EXTRACT_DR (\%disease_risk, $refer, \%diseaserisk_result);

&output_sam_data (\@samples, \%microbiota_result, \%pathway_result, \%ko_result, \%term2bac_result, \%diseasefactor_result, \%diseaserisk_result, $out_dir, $prefix, "extract");

#### output phylum & summary result
&output_phylum (\@samples, $refer, $out_dir, $prefix);
&output_summary (\@samples, $refer, $out_dir, $prefix);

###########################
sub LOAD_DR {
	my ($in,$comp_l) = @_;

	open IN,"$in" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($term, $category, $level, $name) = split/\t+/,$_;
		$comp_l -> {$term}{$category}{$level}{$name} = 1;

	}
	close IN;
}

#&output_xml(\%vars,$dirout,$section);
sub output_xml {
	my ($vars,$dirout,$section)=@_;

	my $outfile="$dirout/$section".".xml";
	open my $out, '>:encoding(utf-8)', $outfile || die $!;
	XMLout($vars,  OutputFile => $out,NoAttr => 1,SuppressEmpty => "",KeyAttr => []);#, NoAttr => 1,SuppressEmpty => "");

	return;
}

sub EXTRACT_T2B {
	my ($in,$ref,$result1) = @_;

	foreach my $t (sort keys %$in) {
		$result1 -> {$t} = 0;
		foreach my $l (sort keys %{$$in{$t}}) {
			foreach my $n (sort keys %{$$in{$t}{$l}}) {
				if (exists $$ref{microbiota}{$l}{$n}) {
					foreach my $s (sort keys %{$$ref{"microbiota"}{$l}{$n}}) {
						$result1 -> {$t} += $$ref{"microbiota"}{$l}{$n}{$s};
					}
				}
			}
		}
	}
}

sub EXTRACT_F2B {
	my ($in,$ref,$result1) = @_;

	foreach my $t (sort keys %$in) {
		foreach my $e (sort keys %{$$in{$t}}) {
			my $key = "$t"."_$e";
			$result1 -> {$key} = 0;
			foreach my $l (sort keys %{$$in{$t}{$e}}) {
				foreach my $n (sort keys %{$$in{$t}{$e}{$l}}) {
					if (exists $$ref{microbiota}{$l}{$n}) {
						foreach my $s (sort keys %{$$ref{"microbiota"}{$l}{$n}}) {
							$result1 -> {$key} += $$ref{"microbiota"}{$l}{$n}{$s};
						}
					}
				}
			}
		}
	}
}

sub EXTRACT_DR {
	my ($in,$ref,$result1) = @_;

	foreach my $t (sort keys %$in) {
		$result1 -> {$t} = 0;
		foreach my $c (sort keys %{$$in{$t}}) {
			foreach my $l (sort keys %{$$in{$t}{$c}}) {
				foreach my $n (sort keys %{$$in{$t}{$c}{$l}}) {
					if (exists $$ref{$c}{$l}{$n}) {
						foreach my $s (sort keys %{$$ref{$c}{$l}{$n}}) {
							$result1 -> {$t} += $$ref{$c}{$l}{$n}{$s};
						}
					}
				}
			}
		}
	}
}

sub LOAD_TERM {
	my ($in,$term_l) = @_;

	open IN,"$in" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my @tmp = split/\t+/,$_;
		if ($#tmp == 1) {
			$term_l -> {$tmp[0]}{$tmp[1]} = 1;
		}
		if ($#tmp == 2) {
			$term_l -> {$tmp[0]}{$tmp[1]}{$tmp[2]} = 1;
		}

	}
	close IN;
}

sub LOAD_BAC {
	my ($in,$bac_l) = @_;

	open IN,"$in" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($term, $level, $name) = split/\t+/,$_;
		$bac_l -> {$term}{$level}{$name} = 1;

	}
	close IN;
}

sub LOAD_F2B {
	my ($in,$bac_l) = @_;

	open IN,"$in" || die $!;
	while (<IN>) {
		chomp;
		s/^\s+|\s+$//;
		next if (/^$/ || /^\#/);
		my ($term, $effect, $level, $name) = split/\t+/,$_;
		$bac_l -> {$term}{$effect}{$level}{$name} = 1;

	}
	close IN;
}

sub output_pop_data {
	my ($sam, $data, $od, $prefix, $type) = @_;

	my $sam_str = join "\t",@$sam;
	open OUT,">$od/$prefix.$type.xls" || die $!;
	print OUT "#$type\t$sam_str\n";
	foreach my $name (sort keys %{$data}) {
		print OUT "$name";
		my $data_str;
		foreach (@$sam) {
			$data_str .= "\t$data->{$name}{$_}";
		}
		print OUT "$data_str\n";
	}
	close OUT;

}

sub output_sam_data {
	my ($sam, $microbiota, $pathway, $ko, $t2b, $d2f, $dr, $od, $prefix, $type) = @_;

	my $sam_str = join "\t",@$sam;
	open OUT,">$od/$prefix.$type.xls" || die $!;
	print OUT "#Term\t$sam_str\n";
	foreach my $name (sort keys %{$microbiota}) {
		print OUT "$name";
		my $data_str;
		foreach (@$sam) {
			$data_str .= "\t$microbiota->{$name}{$_}";
		}
		print OUT "$data_str\n";
	}
	foreach my $name (sort keys %{$pathway}) {
		print OUT "$name";
		my $data_str;
		foreach (@$sam) {
			$data_str .= "\t$pathway->{$name}{$_}";
		}
		print OUT "$data_str\n";
	}
	foreach my $name (sort keys %{$ko}) {
		print OUT "$name";
		my $data_str;
		foreach (@$sam) {
			$data_str .= "\t$ko->{$name}{$_}";
		}
		print OUT "$data_str\n";
	}

	foreach my $name (sort keys %{$t2b}) {
		print OUT "$name\t$t2b->{$name}\n";
	}
	foreach my $name (sort keys %{$d2f}) {
		print OUT "$name\t$d2f->{$name}\n";
	}

	foreach my $name (sort keys %{$dr}) {
		print OUT "$name\t$dr->{$name}\n";
	}

	close OUT;

}

sub output_sam_factor {
	my ($sam, $data, $od, $prefix, $type) = @_;

	my $sam_str = join "\t",@$sam;
	open OUT,">$od/$prefix.$type.xls" || die $!;
	print OUT "#$type\t$sam_str\n";
	foreach my $name (sort keys %{$data}) {
		print OUT "$name\t$data->{$name}\n";
	}
	close OUT;

}

sub REFER_SPLIT {
	my ($ref, $term, $pathogen, $microbiota_result, $pathway_result, $ko_result) = @_;

	my @samples = keys %{$ref->{"summary"}{"diversity"}};
	foreach my $key (sort keys %{$ref->{"summary"}}) {
		if (exists $term->{"summary"}{$key}) {
			$microbiota_result->{$key} = $ref->{"summary"}{$key};
		}
	}
	foreach my $key (sort keys %{$term->{"summary"}}) {
		if (!exists $microbiota_result->{$key}) {
			print STDERR "There no summary $key data in reference population, Set to 0\n";
			foreach my $sam (@samples) {
				$microbiota_result->{$key} = 0;
			}
		}
	}

	foreach my $level (sort keys %{$ref->{"microbiota"}}) {
		foreach my $latin (sort keys %{$ref->{"microbiota"}{$level}}) {
			$microbiota_result->{$latin} = $ref->{"microbiota"}{$level}{$latin};
		}
	}
	foreach my $level (sort keys %{$term->{"microbiota"}}) {
		#print STDERR "$level list of no data in reference population, Set to 0\n";
		foreach my $latin (sort keys %{$term->{"microbiota"}{$level}}) {
			if (!exists $microbiota_result->{$latin}) {
				#print STDERR "\t$latin\n";
				foreach my $sam (@samples) {
					$microbiota_result->{$latin}{$sam} = 0;
				}
			}
		}
	}
	foreach my $class (sort keys %{$pathogen}) {
		#print STDERR "Pathogen list of no data in reference population, Set to 0\n";
		foreach my $level (sort keys %{$pathogen->{$class}}) {
			foreach my $latin (sort keys %{$pathogen->{$class}{$level}}) {
				if (!exists $microbiota_result->{$latin}) {
					foreach my $sam (@samples) {
						$microbiota_result->{$latin}{$sam} = 0;
					}
				}
			}
		}
	}


	foreach my $level (sort keys %{$ref->{"metabolism"}}) {
		foreach my $entry (sort keys %{$ref->{"metabolism"}{$level}}) {
			$pathway_result->{$entry} = $ref->{"metabolism"}{$level}{$entry} if ($level =~ /^L/);
			$ko_result->{$entry} = $ref->{"metabolism"}{$level}{$entry} if ($level =~ /^KO/);
		}
	}
	foreach my $level (sort keys %{$term->{"metabolism"}}) {
		foreach my $entry (sort keys %{$term->{"metabolism"}{$level}}) {
			if (!exists $pathway_result->{$entry} && !exists $ko_result->{$entry}) {
				print STDERR "There no $level $entry data in reference population, Set to 0\n";
				foreach my $sam (@samples) {
					$pathway_result->{$entry}{$sam} = 0 if ($level =~ /^L/);
					$ko_result->{$entry}{$sam} = 0 if ($level =~ /^KO/);
				}
			}
		}
	}
}

sub output_phylum {
	my ($sam, $refer, $out_dir, $prefix) = @_;

	my $sam_str = join "\t",@$sam;
	open OUT,">$out_dir/$prefix.phylum.xls" || die $!;
	print OUT "#phylum\t$sam_str\n";
	foreach my $phy (sort keys %{$refer->{"microbiota"}{"Phylum"}}) {
		my $str = $phy;
		foreach (@$sam) {
			$str .= "\t$refer->{microbiota}{Phylum}{$phy}{$_}";
		}
		print OUT "$str\n";
	}
	close OUT;
}

sub output_summary {
	my ($sam, $refer, $out_dir, $prefix) = @_;

	my $sam_str = join "\t",@$sam;
	open OUT,">$out_dir/$prefix.summary.xls" || die $!;
	print OUT "#summary\t$sam_str\n";
	foreach my $sum (sort keys %{$refer->{"summary"}}) {
		next if ($sum =~ /representation/);
		my $str = $sum;
		foreach (@$sam) {
			$str .= "\t$refer->{summary}{$sum}{$_}";
		}
		print OUT "$str\n";
	}
	close OUT;
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

sub sub_format_datetime 
{#Time calculation subroutine
	my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
	sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

