# Xu-pecin_paper-code
The scripts of paper "Pectin supplementation accelerates post-antibiotic gut microbiome reconstitution orchestrated with reduced gut redox potential"

Metagenome: (use A-3-JS as an example)
Filtering of reads；use trimmomatic 
java -jar -Xms20G -Xmx20G trimmomatic-0.33.jar PE -threads 8 -phred33 rawData/A-3-JS_R1.fq.gz rawData/A-3-JS_R2.fq.gz A-3-JS.clip.1.fq.gz A-3-JS.single.R1.fastq.gz A-3-JS.clip.2.fq.gz A-3-JS.single.R2.fastq.gz ILLUMINACLIP:/mnt/sdb/bin/Trimmomatic-0.33/adapters/merge.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75

remove host contamination: use bwa
bwa mem -t 40 -M 00.ref/host.fa A-3-JS.clip.1.fq.gz A-3-JS.clip.2.fq.gz  |awk '$3!="*"' > A-3-JS.host.sam
perl /mnt/sdb/bin/remove-host.pl A-3-JS.host.sam A-3-JS