#!/bin/bash

##
## Parser for qc trim log files from 16S qc run using merged fastq files
## Log file from qc_trim used as input (nohup concatenated file), all screening logs collected from folder.
##

## Headers for matrix. Hard set screen databases.
echo -e "sample \t trimmed_input_pairs \t trimmed_pairs_remaining \t trimmed_pairs_remaining_% \t trimmed_forward_only \t trimmed_forward_only_% \t trimmed_reverse_only \t trimmed_reverse_only_% \t trimmed_discarded \t trimmed_discarded_% \t screen_input_pairs \t screen_univec_pairs \t screen_univec_pairs_% \t pe_reads_passed_qc"

# Cycle over single run nohup output file. Remove all () for parsing.
grep "^Starting" $1 | sed 's/Starting //' > temp1_sample_id
grep "Input Read Pairs:" $1 | sed 's/(\|)\|%//g' | awk '{print $4,$7,$8,$12,$13,$17,$18,$20,$21}' > temp2_trimmed_metrics

## Set at screening two databases. Hard set databases currently
# $2: input pairs
# $3: unmapped reads against library screen genome
for file in *_screen.txt; do
	awk '{if(NR==3) print $2,($2-$3),((($2-$3)/($2))*100)}' $file >> temp3_screen_1
done

# Count all sequences in final trimmed and screen PE fastq file. Current version based on no PE joining naming.
for file in *R1_001_PE_no_hits_file.1.fastq; do
	seqs=`wc -l $file | awk '{print $1/4}'`; echo $seqs >> temp4_seqs;
done

# Join all temp files to single flat file. All files syncronised by sample name.
paste temp[1234]*

## Tidy up temp files
rm temp[1234]_*


