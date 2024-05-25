#!/usr/bin/perl
use strict;
use warnings;

my $usage = <<"USAGE";
	Usage :perl $0 kegg.need.list kegg.category.xls output
		 kegg.need.list          a list of KO id what you need,one column,required
		 kegg.category.xls       kegg.category.xls from kegg annnotation,required
		 output			 function2gene.txt, required
USAGE

die "$usage\n" unless(@ARGV == 3);
my ($lst,$category,$output) = @ARGV;
my %list;
my %h;
my $k1;
my $k2;

open INL,"$lst" or die "$!\n";
while(<INL>){
	chomp;my @a=split/\t/;
	$list{$a[0]}=1;
}
close INL;
open INK,"$category" or die "$!\n";
while(<INK>){
	chomp;my @a=split/\t/;
	if(exists $list{$a[1]}){
		$h{$a[1]}{$a[0]}=1;
	}
}
close INK;
open OUT,"> $output" or die "$!\n";
foreach $k1(sort keys %h){
	print OUT "$k1";
	foreach $k2(sort keys %{$h{$k1}}){
		print OUT "\t$k2";
	}
	print OUT "\n";

}
close OUT;
