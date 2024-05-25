#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw/$Bin/;

my $opt;
my ($geneProfile,$list,$taxFile,$taxFlag,$outDir);
while($opt = shift){
	if($opt eq "-p"){
		$geneProfile = shift;
	}elsif($opt eq "-l"){
		$list = shift;
	}elsif($opt eq "-f"){
		$taxFlag = shift;
	}elsif($opt eq "-t"){
		$taxFile = shift;
	}elsif($opt eq "-o"){
		$outDir = shift;
	}elsif($opt eq "-h"){
		&usage;
		exit;
	}else{
		print STDERR "unknown parameter: $opt\n";
		&usage;
		exit;
	}
}

unless($geneProfile and $list and $taxFile and $outDir){
	&usage;
	exit;
}

$taxFlag = "P" unless($taxFlag);
$outDir =~ s/\/$//;
`mkdir $outDir` unless(-e $outDir);

my (@order,%abun,%tax,%num2sam);

print "read gene profile ...\n";
open INP,$geneProfile or die "$!\n";
while(<INP>){
	chomp;
	my @temp = split /\t/;
	for(my $i = 1;$i < @temp;$i++){
		if($. == 1){
			push @order,$temp[$i];
			$num2sam{$i} = $temp[$i];
		}else{
			$abun{$temp[0]}{$num2sam{$i}} = $temp[$i];
		}
	}
}
close INP;

print "read taxonomy ...\n";
open INT,$taxFile or die "$!\n";
while(<INT>){
	chomp;
	my @temp = split /\t/;
	next if(/^#/);
	if($taxFlag eq "P"){
		$temp[0] = $1 if($temp[0] =~ /(\S+)\_/);
	}
	$tax{$temp[0]} = $temp[1];
}
close INT;

my (%mess,@function,%used_tax);
open INL,$list or die "$!\n";
while(<INL>){
	chomp;
	next if(/^#/);
	my @temp = split /\t/;
	push @function,$temp[0];
	for(my $i = 1;$i < @temp;$i++){
		if($taxFlag =~ /P/i){
			$temp[$i] = $1 if($temp[$i] =~ /(\S+)\_/);
		}
		foreach my $sam(@order){
			die "error: $temp[$i] >> $sam\n" unless(exists $abun{$temp[$i]}{$sam});
			if(exists $tax{$temp[$i]}){
				$mess{$sam}{$temp[0]}{$tax{$temp[$i]}} += $abun{$temp[$i]}{$sam};
				$used_tax{$tax{$temp[$i]}}++;
			}else{
				$mess{$sam}{$temp[0]}{"Unclassified"} += $abun{$temp[$i]}{$sam};
				$used_tax{"Unclassified"}++;
			}
		}
	}
}
close INL;

`mkdir $outDir/sample` unless(-e "$outDir/sample");
`mkdir $outDir/function` unless(-e "$outDir/function");
`mkdir $outDir/taxonomy` unless(-e "$outDir/taxonomy");
foreach my $s(@order){
	open OUS,"> $outDir/sample/$s.all.xls" or die "$!\n";
	print OUS "#Category\t";
	print OUS join "\t",@function;
	print OUS "\n";
	foreach my $t(sort keys %used_tax){
		print OUS $t;
		my $temp_abun = 0;
		my $temp_mess;
		foreach my $f(@function){
			$mess{$s}{$f}{$t} = 0 unless(exists $mess{$s}{$f}{$t});
			$temp_abun += $mess{$s}{$f}{$t};
			$temp_mess .= "\t$mess{$s}{$f}{$t}";
		}
#		print OUS "$temp_mess\n" if($temp_abun > 0);
		print OUS "$temp_mess\n";
	}
	close OUS;
	`perl $Bin/filter_empty.pl $outDir/sample/$s.all.xls F $outDir/sample/$s.xls`;
	`rm $outDir/sample/$s.all.xls`;
}

my $temp_num = 0;
open OUR,"> $outDir/function/num2function.txt" or die "$!\n";
foreach my $f(@function){
	$temp_num++;
	print OUR "$temp_num\t$f\n";
	open OUF,"> $outDir/function/$temp_num.all.xls" or die "$!\n";
	print OUF "#Category\t";
	print OUF join "\t",@order;
	print OUF "\n";
	foreach my $t(sort keys %used_tax){
		print OUF $t;
		my $temp_abun = 0;
		my $temp_mess;
		foreach my $s(@order){
			$mess{$s}{$f}{$t} = 0 unless(exists $mess{$s}{$f}{$t});
			$temp_abun += $mess{$s}{$f}{$t};
			$temp_mess .= "\t$mess{$s}{$f}{$t}";
		}
#		print OUF "$temp_mess\n" if($temp_abun > 0);
		print OUF "$temp_mess\n";
	}
	close OUF;
	`perl $Bin/filter_empty.pl $outDir/function/$temp_num.all.xls F $outDir/function/$temp_num.xls`;
	`rm $outDir/function/$temp_num.all.xls`;
}
close OUR;

$temp_num = 0;
open OUTR,"> $outDir/taxonomy/num2tax.txt" or die "$!\n";
foreach my $t(sort keys %used_tax){
	$temp_num++;
	print OUTR "$temp_num\t$t\n";
	open OUTT,"> $outDir/taxonomy/$temp_num.all.xls" or die "$!\n";
	print OUTT "#Category\t";
	print OUTT join "\t",@order;
	print OUTT "\n";
	foreach my $f(@function){
		print OUTT $f;
		my $temp_abun = 0;
		my $temp_mess;
		foreach my $s(@order){
			$mess{$s}{$f}{$t} = 0 unless(exists $mess{$s}{$f}{$t});
			$temp_abun += $mess{$s}{$f}{$t};
			$temp_mess .= "\t$mess{$s}{$f}{$t}";
		}
#		print OUTT "$temp_mess\n" if($temp_abun > 0);
		print OUTT "$temp_mess\n";
	}
	close OUTT;
	`perl $Bin/filter_empty.pl $outDir/taxonomy/$temp_num.all.xls F $outDir/taxonomy/$temp_num.xls`;
	`rm $outDir/taxonomy/$temp_num.all.xls`;
}
close OUTR;

sub usage{
	print <<EOD
	usage: perl $0 -p geneProfile.xls -l function2gene.list -f taxSeq.flag -t tax.xls -o outDir
		-p	gene profile, required
		-l	function2gene list, required
		-f	taxononmy sequence type[N:nuc|P:pro], default P
		-t	gene taxonomy file, required
		-o	out direction, required
EOD
}
