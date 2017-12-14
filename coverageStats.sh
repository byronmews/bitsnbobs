#!/bin/bash

#
# Simple coverge statistics from a reference based mapping.
# Input a sorted bam file.
#

# Convert sorted bam output to bed format.
bedtools bamtobed -i $1 > $1.bed

# Index file and parse first two columns to generate a bedtools genome file. Ignore last line.
samtools idxstats $1 | grep -v '*' | awk '{print $1 "\t" $2}' > $1.genome

## This keeps multimapped reads within the metrics, not an issue so far
# Show base coverage by genomic position
echo "Running genome coverage for $1."
bedtools genomecov -d -i $1.bed -g $1.genome > $1.bed_coverage

# Show base coverage by genomic position removing multimapped reads, if the multiple mapped reads need to be removed
#samtools depth -Q 1 $1 > $1.bed_coverage

# Calculate number of bases, genome length and mean fold coverage/depth
awk '{sum+=$3} END {print "bases= " sum "\n" "genome_length= "NR "\n" "foldX= "sum/NR}' $1.bed_coverage > $1.coverage_1

# Calculate number of bases with 1+ base coverage, show as % coverage
awk '$3>=1 {count++} END {print "percent_coverage= "((count/NR)*100)}' $1.bed_coverage > $1.coverage_2

# Tidy up the coverage files, deleting any previous runs
> $1.genome_coverage.txt
cat $1.coverage_1 $1.coverage_2 >> $1.genome_coverage.txt
rm  $1.coverage_*


echo "Completed mapping for $1"
echo "Output file: $1.genome_coverage.txt"
