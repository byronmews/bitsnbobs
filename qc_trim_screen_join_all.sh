#!/bin/bash

# qc_trim_screen_join_all.sh
#
# Fastq QC, with stages consisting of: FastQC > Trimmomatic_1 > FastQC > Fastq Screen > FLASH  > Trimmomatic_2 > FastQC
# Uses multiple configuration files for host screen depending on specified sample type.
# Script geared for longer 500/600 cycle MiSeq PE runs for amplicon 16S sequencing. Joins overlapping PE reads.
# Various runtime messsages to stout, can be saved for log file parsing.
#
# GRU01 server config. * only human and adapter paths setup.
#
# To run on 16S samples use the below command:
#
# 	bash qc_trim_screen_join_all.sh adapter
# 
# Author: Graham Rose
#




## Loop over all PE libraries in directory
for file_r1 in *R1_001.fastq.gz; do
	file_r2=`echo $file_r1 | sed 's/_R1_/_R2_/'`

	## Input screen config file name, if missing then exit
	if [ $# -ne 1 ] ; then
		echo 'Input: [fastq_screen.conf]'
		exit 0
	fi


	base_1=`echo $file_r1 | sed 's/_R1_.*.fastq.gz//'`
	base_2=`echo $file_r2 | sed 's/_R2_.*.fastq.gz//'`
	length=$2; # useful for amplicon projects eg. Illumina 16S V3/V4 expected length ~460bp


	# Catch fastqscreen genome(s) dbs to use
	if [ $1 == "nibsc" ]
	then
		screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_nibsc.conf
	elif [ $1 == "unknown" ]
	then
		screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_unknowns.conf
	elif [ $1 == "pig" ]
	then
		screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_pig.conf
	elif [ $1 == "human" ]
	then
		screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_human.conf
	elif [ $1 == "ebola" ]
	then
		screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_ebola.conf
	elif [ $1 == "adapter" ]
	then
       		screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_adapter_only.conf
	else
		printf "\nEnter library fastq_screen conf name to screen fastq: nibsc, unknown, pig, human, ebola, adapter \n\n"
		exit 0;
	fi


	mkdir fastqc_pre_qc
	mkdir fastqc_post_qc_1
	mkdir fastqc_post_qc_2

	# Start time
	printf "Starting $base_1\n"
	date

	echo "FastQC pre trimmed fastqs"
	fastqc -k 8 --nogroup $file_r1 -t 12 -o fastqc_pre_qc
	fastqc -k 8 --nogroup $file_r2 -t 12 -o fastqc_pre_qc

	echo "Trimming fastqs - soft trim pre joining ends"
	# Contaminant list originates from fastq library, with any solid  platform primers removed from list. Prob redundant for qualilty aware joiners. Hexamer bias from 5' end kept in this trimming round. No trimming within the sequence either.
	java -jar /usr/local/src/trimmomatic-0.32.jar PE -phred33 $file_r1 $file_r2 $base_1"_R1_001_PE.fastq" $base_1"_R1_001_SE.fastq" $base_2"_R2_001_PE.fastq" $base_2"_R2_001_SE.fastq" ILLUMINACLIP:/srv/data0/dbs/trimmomatic_db/contaminant_list.fasta:2:30:10 LEADING:10 TRAILING:10 MINLEN:40

	echo "FastQC of trimmed fastqs"
	# Fastqc post qc trim
	fastqc -k 8 --nogroup $base_1"_R1_001_PE.fastq" -t 12 -o fastqc_post_qc_1
	fastqc -k 8 --nogroup $base_2"_R2_001_PE.fastq" -t 12 -o fastqc_post_qc_1

	echo "Screening qc passed fastqs"
	# Screen fastq against set libraries, recover those that do not map to these
	fastq_screen --conf $screen_conf --aligner bowtie2 --threads 8 --nohits --paired $base_1"_R1_001_PE.fastq" $base_2"_R2_001_PE.fastq"

	echo "Joining soft trimmed pairs - 100bp-200bp overlap expected"
	flash --max-overlap 200 --threads 12 --output-prefix $base_1"_no_hits" $base_1"_R1_001_PE_no_hits_file.1.fastq" $base_2"_R2_001_PE_no_hits_file.2.fastq"

	echo "Trimming fastqs - hard trim of joined pairs"
	# Contaminant list originates from fastq library, with any solid  platform primers removed from list. Prob overkill for the projects though. Hexamer bias from 5' end kept in this trimming round. Trimming within sequence.
	java -jar /usr/local/src/trimmomatic-0.32.jar  SE -phred33 $base_1"_no_hits.extendedFrags.fastq" $base_1"_no_hits.extendedFrags.trimmed.fastq" ILLUMINACLIP:/srv/data0/dbs/trimmomatic_db/contaminant_list.fasta:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

	echo "FastQC of joined and trimmed fastqs"
	# Fastqc post qc trim
	fastqc -k 8 --nogroup $base_1"_no_hits.extendedFrags.trimmed.fastq" -t 12 -o fastqc_post_qc_2

	echo "Done"

	# End time
	date

done


