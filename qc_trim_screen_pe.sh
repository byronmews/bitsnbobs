#!/bin/bash

#
# Fastq QC stages consisting of: FastQC > Trimmomatic > FastQC > Fastq Screen
#
# Uses multiple configuration files for host screen depending on specifed sample type
#
# Only human setup for GRU01 server paths

## Input fastq
if [ $# -ne 3 ] ; then
    echo 'Input: [PE R1 fastq.gz] [PE R2 fastq.gz] [fastq_screen.conf]'
    exit 0
fi

base_1=`echo $1 | sed 's/_R1_.*.fastq.gz//'`
base_2=`echo $2 | sed 's/_R2_.*.fastq.gz//'`

# Catch screen dbs genomes to use
        if [ $3 == "nibsc" ]
        then
                screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_nibsc.conf
        elif [ $3 == "unknown" ]
        then
                screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_unknowns.conf
        elif [ $3 == "pig" ]
        then
                screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_pig.conf
        elif [ $3 == "human" ]
        then
                screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_human.conf
        elif [ $3 == "ebola" ]
        then
                screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_ebola.conf
        elif [ $3 == "adapter" ]
        then
                screen_conf=/usr/local/etc/fastq_screen_v0.4.4_config/fastq_screen_adapter_only.conf
        else
                printf "\nEnter library fastq_screen conf name to screen fastq: nibsc, unknown, pig, human, ebola, adapter \n\n"
                exit 0;
        fi

mkdir fastqc_pre_qc
mkdir fastqc_post_qc

# Start time
date

echo "FastQC pre trimmed fastqs"
#fastqc pre qc trim
fastqc -k 8 --nogroup $1 -t 6 -o fastqc_pre_qc
fastqc -k 8 --nogroup $2 -t 6 -o fastqc_pre_qc

echo "Trimming fastqs, no headcrop set (15bp)"
# Contaminant list originates from fastq library, with any solid platform primers removed. Prob overkill for the projects though. Crop 5' 12bp to remove common hexamer bias issue currently disabled.
java -jar /usr/local/src/trimmomatic-0.32.jar PE -phred33 $1 $2 $base_1"_R1_001_PE.fastq" $base_1"_R1_001_SE.fastq" $base_2"_R2_001_PE.fastq" $base_2"_R2_001_SE.fastq" ILLUMINACLIP:/srv/data0/dbs/trimmomatic_db/contaminant_list.fasta:2:30:10 LEADING:20 TRAILING:20 MINLEN:40

echo "FastQC of trimmed fastqs"
#fastqc post qc trim
fastqc -k 8 --nogroup $base_1"_R1_001_PE.fastq" -t 12 -o fastqc_post_qc
fastqc -k 8 --nogroup $base_2"_R2_001_PE.fastq" -t 12 -o fastqc_post_qc

echo "Screening qc passed fastqs"
# screen fastq against set libraries, recover those that do not map to these
fastq_screen --conf $screen_conf --aligner bowtie2 --threads 12 --nohits --paired $base_1"_R1_001_PE.fastq" $base_2"_R2_001_PE.fastq"

echo "Done"

# End time
date
