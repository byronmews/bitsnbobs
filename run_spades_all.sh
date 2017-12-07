#!/bin/bash

#################################################################################################
# Pipeline: non_NEC_metagenomic_annotation_analysis.sh 						#
# Author: graham.rose@phe.gov.uk								#
# Project: Simon Kroll collaboration								#
#  NEC metagenomic shotgun datasets AMR analysis. Runs on QC passed non binned reads		#
# 												#
#################################################################################################

# Run as test case already: Q89_111025_S13

# nonNEC samples miseq repeats (N=3)
samples=(
Q87_S20
)


# Cycle all samples, setup dirs, assemble, check, annnotate
for i in "${samples[@]}"
do

	# Setup variables
	fastqR1=$i"_R1_001_PE_no_hits_file.1.fastq"
        fastqR2=$i"_R2_001_PE_no_hits_file.2.fastq"
	base_name=$i
        
	echo "Running $i"
        
#	mkdir $i
#	
#	# Symlink relevant qc passed fastqs to dir
#	ln -s /srv/data0/graham/projects/simon_kroll_lab/AMR_dataset_N11_plusreps/qc/$i"_R1_001_PE_no_hits_file.1.fastq" -t /srv/data0/graham/projects/simon_kroll_lab/AMR_dataset_N11_plusreps/spades/$i
#	ln -s /srv/data0/graham/projects/simon_kroll_lab/AMR_dataset_N11_plusreps/qc/$i"_R2_001_PE_no_hits_file.2.fastq" -t /srv/data0/graham/projects/simon_kroll_lab/AMR_dataset_N11_plusreps/spades/$i
#
#	# Assemble non binned PE reads, using set kmer range
#        echo "Assembling using spades with --meta option"
#	spades_v3.7.1 --meta -1 $i/$fastqR1 -2 $i/$fastqR2 -o $i/$base_name"_spades_assembly"
#	
#	# Run QUAST on assembly
#	echo "Running QUAST on assembly $i"
#	quast.py $i/$base_name"_spades_assembly"/contigs.fasta -o $i/$base_name"_spades_assembly"/quast
#	
	# Run PROKKA on assembly
	echo "Running PROKKA on assembly $i"
	mkdir $base_name/$base_name"_prokka"
        # setup new fasta with headers <20 chars
        echo "Copying and renaming fasta file"
        sed 's/_length.*//' $base_name/$base_name"_spades_assembly/contigs.fasta" > $base_name/$base_name"_spades_assembly_contigs.prokka_ready.fasta"
        prokka --locustag $base_name --force --metagenome --cpus 12 --outdir $base_name/$base_name"_prokka" --prefix $base_name $base_name/$base_name"_spades_assembly_contigs.prokka_ready.fasta"
	



echo "Done"
done

## End
