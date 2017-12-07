#!/bin/bash

#
# Parsing of the output from the groovy script [mapping_process_single_sam.bpipe]
#
# Takes all *.bam_report.txt and *.bam.genome_coverage.txt files with the current
# directory and creates a single tabulated report table.
#
# Currently plugged into pipelines after writing mapping reports for multiple sam files
#


# Table header labels
echo -e "sample \t reads_mapped \t reads_mapped_% \t bases_mapped \t ref_base_length \t fold_coverage \t genome_coverage_%"

# Loop over any report files in the directory and parse into one tab delim table
for file in *.bam.report.txt; do
	
	# Parse text from *.bam.report.txt file
	# Base string to use as output
	sample_name=`echo $file | sed 's/.bam.report.txt//'`
	
	cov_report=`echo $file | sed 's/.bam.report.txt/.bam.genome_coverage.txt/'`
	mapped=`grep "Mapped" $file | awk '{print $3}'`
	mapped_pct=`grep "Mapped" $file | awk '{print $4}' | sed 's/(\|)\|%//g'`
	
	# Parse text from *.bam.genome_coverage.txt file (base string sync'd with *.bam.report.txt file)
	bases=`grep -w "bases" $cov_report | awk '{print $2}'`
	genome=`grep -w "genome_length" $cov_report | awk '{print $2}'`
	fold_cov=`grep -w "foldX" $cov_report | awk '{print $2}'`
	base_cov=`grep -w "percent_coverage" $cov_report | awk '{print $2}'`

	# Results table output
	echo -e "$sample_name \t $mapped \t $mapped_pct \t $bases \t $genome \t $fold_cov \t $base_cov"
done


