#!/bin/bash

#
# Run bbsplit (bbmap aligner) to bin a PE fastq into separate files, ignoring reads mapping
# to multiple genomes. Requires previous indexing of the genome fasta files.
#

OUT_1=`echo $1 | sed 's/_no_hits_.*.fastq/_bbsplit_unmapped.fastq/'`
OUT_2=`echo $2 | sed 's/_no_hits_.*.fastq/_bbsplit_unmapped.fastq/'`

~/programs/bbmap/bbsplit.sh build=2 path=~/data/data3/graham/nibsc/nibsc_references/dna_genomes/bbmap in1=$1 in2=$2 minratio=0.56 minhits=1 maxindel=16000 ambig2=split basename=split_%.sam outu1=$OUT_1 outu2=$OUT_2
