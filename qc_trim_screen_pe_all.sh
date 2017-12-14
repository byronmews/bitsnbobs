#!/bin/bash

#
# Fastq QC stages consisting of: FastQC > Trimmomatic > FastQC > Fastq Screen
#
# Uses multiple configuration files for host screen depending on specified sample type
# This all version takes all PE files within the folder and runs the analysis. Only needs screen config file as argument.

# Only human and adapter setup for GRU01 server


## Loop over all PE libraries in folder
for file_r1 in *R1_001.fastq.gz; do

	file_r2=`echo $file_r1 | sed 's/_R1_/_R2_/'`

	## Input screen config file name, if missing then exit
	if [ $# -ne 1 ] ; then
    		echo 'Input: [fastq_screen.conf]'
			exit 0
	fi


	base_1=`echo $file_r1 | sed 's/_R1_.*.fastq.gz//'`
	base_2=`echo $file_r2 | sed 's/_R2_.*.fastq.gz//'`

	# Catch eukaryotic screen to use
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
	mkdir fastqc_post_qc

	# Start time
	echo "Starting $base_1"
	date

	echo "FastQC pre trimmed fastqs"
	#fastqc pre qc trim
	fastqc -k 8 --nogroup $file_r1 -t 10 -o fastqc_pre_qc
	fastqc -k 8 --nogroup $file_r2 -t 10 -o fastqc_pre_qc

	echo "Trimming fastqs, no headcrop set (15bp)"
	# Contaminant list originates from fastq library, with any solid platform primers removed. Prob overkill for the projects though. Crop 5' 12bp to remove common hexamer bias issue currently disabled.
	java -jar /usr/local/src/trimmomatic-0.32.jar PE -phred33 $file_r1 $file_r2 $base_1"_R1_001_PE.fastq" $base_1"_R1_001_SE.fastq" $base_2"_R2_001_PE.fastq" $base_2"_R2_001_SE.fastq" ILLUMINACLIP:/srv/data0/dbs/trimmomatic_db/contaminant_list.fasta:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

	echo "FastQC of trimmed fastqs"
	#fastqc post qc trim
	fastqc -k 8 --nogroup $base_1"_R1_001_PE.fastq" -t 10 -o fastqc_post_qc
	fastqc -k 8 --nogroup $base_2"_R2_001_PE.fastq" -t 10 -o fastqc_post_qc

	echo "Screening qc passed fastqs"
	# screen fastq against set libraries, recover those that do not map to these
	fastq_screen --conf $screen_conf --aligner bowtie2 --threads 12 --nohits --paired $base_1"_R1_001_PE.fastq" $base_2"_R2_001_PE.fastq"

	### Parse all qc files to output summary matrix
	# todo
	
	echo "Done"

	# End time
	echo "Ending file $base_1"
	date

done
