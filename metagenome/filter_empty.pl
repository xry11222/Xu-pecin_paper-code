#!/usr/bin/perl -w
use strict;
use warnings;

die "usage: perl $0 input.table whether-metadata[T|F] output\n" unless(@ARGV == 3);
my ($input,$flag,$output) = @ARGV;

my %abun;
open IN,$input or die "$!\n";
while(<IN>){
	chomp;
	next if(/^#/ or $. == 1);
	my @temp = split /\t/;
	pop @temp if($flag eq "T");
	for(my $i = 1;$i < @temp;$i++){
		$abun{$i} += $temp[$i];
	}
}
close IN;

open IN2,$input or die "$!\n";
open OUT,"> $output" or die "$!\n";
while(<IN2>){
	chomp;
	my @temp = split /\t/;
	my $metadata = pop @temp if($flag eq "T");
	my $num = 0;
	my $mess = $temp[0];
	for(my $i = 1;$i < @temp;$i++){
		die "error: $_\n" unless(exists $abun{$i});
		next if($abun{$i} == 0);
		$num += $temp[$i] unless($. == 1);
		$mess .= "\t$temp[$i]";
	}
	if($flag eq "T"){
		print OUT "$mess\t$metadata\n" if($num > 0 or $. == 1);
	}else{
		print OUT "$mess\n" if($num > 0 or $. == 1);
	}
}
close IN2;
close OUT;
