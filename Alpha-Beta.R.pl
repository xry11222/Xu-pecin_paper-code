#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
my %opts;
my $VERSION = "v2.1-20140214"; 
my $UpdatedVersion = "V2.20170630";
GetOptions( \%opts,"i=s","o=s");
my $usage = <<"USAGE";
              Usage :perl $0 [options]		
			-i	 * input distance matrix 
			-o * output prefix
	   Example:$0 

USAGE

die $usage if(!($opts{i}&&$opts{o}));

$opts{i}=~/([^\/]+)$/;
my $bn=$1;
		
open CMD,">$opts{o}.cmd.r";
print CMD "

library(vegan)
dr <- read.table(\"$opts{i}\",header=TRUE, row.names = 1)
dr <- dr[which(rowSums(dr) > 0),]
dr <- t(dr)
beta <- vegdist(dr, method= 'bray')
beta<-as.matrix(beta)
shannon_index <- diversity(dr, index = 'shannon',MARGIN = 1, base = exp(1))
richness <- rowSums(dr > 0)
res<-cbind(shannon_index,richness) 
write.csv(res,\'$opts{o}_alpha.csv\', quote = FALSE)
write.csv(beta,\'$opts{o}_beta.csv\', quote = FALSE)

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
    data.frame(
	    row = rownames(cormat)[row(cormat)[ut]],
		    column = rownames(cormat)[col(cormat)[ut]],
			    dist  =(cormat)[ut]
				  )
				  }
MM<-flattenCorrMatrix(beta)
write.table(MM,\'$opts{o}_beta.1v1.csv\', sep=\'\,\')
";

`R --restore --no-save < $opts{o}.cmd.r`;
