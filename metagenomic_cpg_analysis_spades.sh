#!/bin/bash

#################################################################################################
# Pipeline: metagenomic_cpg_analysis.sh 							#
# Author: graham.rose@phe.gov.uk								#
# Project: Simon Kroll collaboration								#
#  NEC metagenomic shotgun datasets CpG analysis. Runs after megan binning, and uses tab delim	#
#  file of species level fastq read annotations.						#
# 												#
#################################################################################################


## Arguments. Defaults set here.
ASSEMBLE_OPT="no_assembly" # default is no

## Help file function
function Help
{
	printf "\n"
	printf "Usage:\tmetagenomics_cpg_analysis.sh [options]\n\n"
	printf "Options:\t-i R1 fastq file. Must be input.\n"
	printf "	\tfastq format: eg. 142_S1_L001_R1_001_PE_no_hits_file.1.fastq\n"
	printf "	\t-a optional split fastqs and run spades assembly (default is skip to motif search)\n"
	printf "	\t-m megan read species annotation file. Must be input if -a set\n"
	printf "	\t-h help\n"
	printf "\n"
}

NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
	Help
	exit 1
fi

## Argument handling
while getopts ":ahm:i:" OPTION; do
	case $OPTION in
	a) ASSEMBLE_OPT="run_assembly"; echo "Assembly -a flag set";;
	h) Help; exit 1;;
	i) INPUT_FASTQ_OPT=$OPTARG;  printf "Input R1 fastq file: $OPTARG\n";;
	m) MEGAN_ANNOTATED_READS_OPT=$OPTARG; printf "Input megan annotation file: $OPTARG\n";; # todo
	\?) Help; exit 1;;
	--) echo "Missing option"; Help; exit 1;;
	:) echo "Missing option"; Help; exit 1;;
	*) echo "Unknown option"; Help; exit 1;; 	
	esac
done

# Setup global variables, most transferred from getopts input. Light error checking.
FASTQ_R1=$INPUT_FASTQ_OPT
FASTQ_R2=`echo $FASTQ_R1 | sed 's/_R1_/_R2_/' | sed 's/_file.1.fastq/_file.2.fastq/'`
MEGAN_ANNOTATED_READS=$MEGAN_ANNOTATED_READS_OPT
FILE_BASE_NAME=`echo $FASTQ_R1 | sed 's/_R1_001_PE.*.fastq//'`


## If -a flag used, run fastq splitter
if [ $ASSEMBLE_OPT == "run_assembly" ]
then

	# Run once only. Command line flag enables this.
	# Separate the megan read to taxon_id file into species (taxon_id) bins named by taxon_id
	printf "Splitting annotated reads within sample: $FILE_BASE_NAME\n"
	awk -F " " '{print > "taxon_id_"$2".txt"}' $MEGAN_ANNOTATED_READS
	printf "Reconstructing PE species fastqs\n"
	## Build new fastq files looping over each read_id bin
	for file in taxon_id_*.txt; do
		base=`echo $file | sed 's/taxon_id_//' | sed 's/.txt//'`
		printf "Extracting species taxon_id: $base\n"
		awk '{print $1}' $file | seqtk subseq $FASTQ_R1 - > "taxon_id_"$base"_R1.fastq"
		awk '{print $1}' $file | seqtk subseq $FASTQ_R2 - > "taxon_id_"$base"_R2.fastq"
	done
fi

## If -a flag used, run assembler
if [ $ASSEMBLE_OPT == "run_assembly" ]
then
	# Remove any previous assembly run directories
	#printf "Removing previous runs\n"
	#rm -r taxon_id_*assembly
	# Assemble binned species PE reads, using set kmer range
	printf "Assembling using spades with kmer lengths 21,33,55 using --sc option for uneven depths\n"
	for file in taxon_id_*.txt; do
		base=`echo $file | sed 's/taxon_id_//' | sed 's/.txt//'`
		printf "Assembling species taxon_id: $base\n"
		spades.py --sc -k 21,33,55 --careful -1 "taxon_id_"$base"_R1.fastq" -2 "taxon_id_"$base"_R2.fastq" -o "taxon_id_"$base"_spades_assembly"
		quast.py "taxon_id_"$base"_spades_assembly/contigs.fasta" -o "taxon_id_"$base"_spades_assembly/quast"
	done
fi

## Motif search for CpG motif hexamer within sequences (using NNCGNN).
# Message if skipping fastq splitter and assebly stages.
if [ $ASSEMBLE_OPT == "no_assembly" ]
then
	printf "\nSkipping fastq split and spades assemblies\n\n"
elif [ $ASSEMBLE_OPT == "run_assembly" ]
then
	printf "Fuzznuc searching assemblies for generic CpG hexamer NNCGNN\n"
	## Motif search for CpG motif hexamer within sequences (NNCGNN)
	for file in taxon_id_*_spades_assembly; do
		fuzznuc -sequence $file"/contigs.fasta" -pattern NNCGNN -stdout -auto -complement > $file".spades_fuzznuc_nncgnn.out"
	done
fi

## Cycle primary NNCGNN file for all possible 6-mer motifs with central CG dinucleotide
# Iterations kept in sync by filename. 
#
# Motifs array to lookup against primary motif output per species assembly. Column 7 lists string matched.
# GTCGTT=Human (Hartnann, 2000)
# GACGTT=Murine (Yi et al, 1998)

cpg_motifs=(AACGAA
AACGAT
AACGAG
AACGAC
AACGTA
AACGTT
AACGTG
AACGTC
AACGGA
AACGGT
AACGGG
AACGGC
AACGCA
AACGCT
AACGCC
ATCGAA
ATCGAT
ATCGAG
ATCGAC
ATCGTA
ATCGTT
ATCGTG
ATCGTC
ATCGGA
ATCGGT
ATCGGG
ATCGGC
ATCGCA
ATCGCT
ATCGCC
AGCGAA
AGCGAT
AGCGAG
AGCGAC
AGCGTA
AGCGTT
AGCGTG
AGCGTC
AGCGGA
AGCGGT
AGCGGG
AGCGGC
AGCGCA
AGCGCT
AGCGCC
ACCGAA
ACCGAT
ACCGAG
ACCGAC
ACCGTA
ACCGTT
ACCGTG
ACCGTC
ACCGGA
ACCGGT
ACCGGG
ACCGGC
ACCGCA
ACCGCT
ACCGCC
TACGAA
TACGAT
TACGAG
TACGAC
TACGTA
TACGTT
TACGTG
TACGTC
TACGGA
TACGGT
TACGGG
TACGGC
TACGCA
TACGCT
TACGCC
TTCGAA
TTCGAT
TTCGAG
TTCGAC
TTCGTA
TTCGTT
TTCGTG
TTCGTC
TTCGGA
TTCGGT
TTCGGG
TTCGGC
TTCGCA
TTCGCT
TTCGCC
TGCGAA
TGCGAT
TGCGAG
TGCGAC
TGCGTA
TGCGTT
TGCGTG
TGCGTC
TGCGGA
TGCGGT
TGCGGG
TGCGGC
TGCGCA
TGCGCT
TGCGCC
TCCGAA
TCCGAT
TCCGAG
TCCGAC
TCCGTA
TCCGTT
TCCGTG
TCCGTC
TCCGGA
TCCGGT
TCCGGG
TCCGGC
TCCGCA
TCCGCT
TCCGCC
GACGAA
GACGAT
GACGAG
GACGAC
GACGTA
GACGTT
GACGTG
GACGTC
GACGGA
GACGGT
GACGGG
GACGGC
GACGCA
GACGCT
GACGCC
GTCGAA
GTCGAT
GTCGAG
GTCGAC
GTCGTA
GTCGTT
GTCGTG
GTCGTC
GTCGGA
GTCGGT
GTCGGG
GTCGGC
GTCGCA
GTCGCT
GTCGCC
GGCGAA
GGCGAT
GGCGAG
GGCGAC
GGCGTA
GGCGTT
GGCGTG
GGCGTC
GGCGGA
GGCGGT
GGCGGG
GGCGGC
GGCGCA
GGCGCT
GGCGCC
GCCGAA
GCCGAT
GCCGAG
GCCGAC
GCCGTA
GCCGTT
GCCGTG
GCCGTC
GCCGGA
GCCGGT
GCCGGG
GCCGGC
GCCGCA
GCCGCT
GCCGCC
CACGAA
CACGAT
CACGAG
CACGAC
CACGTA
CACGTT
CACGTG
CACGTC
CACGGA
CACGGT
CACGGG
CACGGC
CACGCA
CACGCT
CACGCC
CTCGAA
CTCGAT
CTCGAG
CTCGAC
CTCGTA
CTCGTT
CTCGTG
CTCGTC
CTCGGA
CTCGGT
CTCGGG
CTCGGC
CTCGCA
CTCGCT
CTCGCC
CCCGAA
CCCGAT
CCCGAG
CCCGAC
CCCGTA
CCCGTT
CCCGTG
CCCGTC
CCCGGA
CCCGGT
CCCGGG
CCCGGC
CCCGCA
CCCGCT
CCCGCC
)

for file in taxon_id_*_spades_assembly.spades_fuzznuc_nncgnn.out; do
	base=`echo $file | sed 's/_spades_assembly.*//'`
	printf "$base \n" >> "headers_cpg_motifs.temp"
done

for i in "${cpg_motifs[@]}"; do
	printf "Counting motif: $i\n"
	for file in taxon_id_*_spades_assembly.spades_fuzznuc_nncgnn.out; do 
		matches=`grep -wc $i $file`
		#printf "$i=$matches \n" >> $i".single_cpg_motifs.temp" #Debugging removed motif string from hits#
		printf "$matches \n" >> $i".single_cpg_motifs.temp"
	done
done

# Output to file and clean up temps
paste headers_cpg_motifs.temp *.single_cpg_motifs.temp > $FILE_BASE_NAME".contigs_cpg_motifs_all.temp"
rm *single_cpg_motifs.temp
rm headers_cpg_motifs.temp

## Parse all generated files, join tables at end, based on taxon_id sting and output to single tab delim file
#
# Temp files generated:
# cpg_motifs_all.temp
# contigs_gc.temp
# read_number.temp
# contigs_cpg_sites.temp 

## Calc number of reads per fastq file for GC normalisation. Output to temp.
printf "Calculating reads input per fastq bin \n"
for file in taxon_id*_R1.fastq; do
	base=`echo $file | sed 's/_R1.fastq//'`
	fastq_lines=`wc -l $file | cut -f1 -d ' '`
	fasta_denominator=4
	lines=`expr $fastq_lines / $fasta_denominator`
	printf "$base $lines \n" >> $FILE_BASE_NAME".read_number.temp"
done

## Output base nuc freq table for CpG normalisation. AT also calculated for checksum.
# Calculate frequencies as pseudo concatenated sequence.
printf "Calculating GC per assembly \n"
for file in taxon_id_*_spades_assembly/contigs.fasta; do 
	base=`echo $file | sed 's/_spades_assembly.*//'` 
	bases_atgc=`grep -v ">" $file | awk 'BEGIN {a=0; t=0; g=0; c=0}; {a+=gsub("A",""); t+=gsub("T",""); g+=gsub("G",""); c+=gsub("C","")} END {print a"\t"t"\t"g"\t"c}'`
	printf "$base \t $bases_atgc \n" >> $FILE_BASE_NAME".contigs_gc.temp"
done

## Output CpG sites table.
printf "Parsing fuzznuc files\n"
for file in taxon_id_*_spades_assembly.spades_fuzznuc_nncgnn.out; do 
	v1=`echo $file | sed 's/_spades_assembly.*//'`
	v2=`grep "Total_sequences" $file | awk '{print $3}'`
	v3=`grep "Total_length" $file | awk '{print $3}'`
	v4=`grep "Reported_sequences" $file | awk '{print $3}'`
	v5=`grep "Reported_hitcount" $file | awk '{print $3}'`
	printf "$v1\t$v2\t$v3\t$v4\t$v5\n" >> $FILE_BASE_NAME".contigs_cpg_sites.temp"
done

## Pool all temp files using taxon column 1 as key.
# Write header row first.
# Initialise variable
echo -n "" > $FILE_BASE_NAME".master.txt"
# Temp file 1: read pairs per species
printf	"taxon_id\t" >> $FILE_BASE_NAME".master.txt"
printf	"read_pairs\t" >> $FILE_BASE_NAME".master.txt"
# Temp file 2: atgc base frequencies per species of assembled contigs
printf	"a_bases\t" >> $FILE_BASE_NAME".master.txt"
printf	"t_bases\t" >> $FILE_BASE_NAME".master.txt"
printf	"g_bases\t" >> $FILE_BASE_NAME".master.txt"
printf	"c_bases\t" >> $FILE_BASE_NAME".master.txt"
# Temp file 3: cpg generic hexamers per assembly
printf	"total_sequences\t" >> $FILE_BASE_NAME".master.txt"
printf	"total_length\t" >> $FILE_BASE_NAME".master.txt"
printf	"reported_sequences\t" >> $FILE_BASE_NAME".master.txt"
printf	"reported_cg_hitcount\t" >> $FILE_BASE_NAME".master.txt"
# The motif header row
printf  "AACGAA	AACGAC	AACGAG	AACGAT	AACGCA	AACGCC	AACGCT	AACGGA	AACGGC	AACGGG	AACGGT	AACGTA	AACGTC	AACGTG	AACGTT	ACCGAA	ACCGAC	ACCGAG	ACCGAT	ACCGCA	ACCGCC	ACCGCT	ACCGGA	ACCGGC	ACCGGG	ACCGGT	ACCGTA	ACCGTC	ACCGTG	ACCGTT	AGCGAA	AGCGAC	AGCGAG	AGCGAT	AGCGCA	AGCGCC1	AGCGCT	AGCGGA	AGCGGC	AGCGGG	AGCGGT	AGCGTA	AGCGTC	AGCGTG	AGCGTT	ATCGAA	ATCGAC	ATCGAG	ATCGAT	ATCGCA	ATCGCC	ATCGCT	ATCGGA	ATCGGC	ATCGGG	ATCGGT	ATCGTA	ATCGTC	ATCGTG	ATCGTT	CACGAA	CACGAC	CACGAG	CACGAT	CACGCA	CACGCC	CACGCT	CACGGA	CACGGC	CACGGG	CACGGT	CACGTA	CACGTC	CACGTG	CACGTT	CCCGAA	CCCGAC	CCCGAG	CCCGAT	CCCGCA	CCCGCC	CCCGCT	CCCGGA	CCCGGC	CCCGGG	CCCGGT	CCCGTA	CCCGTC	CCCGTG	CCCGTT	CTCGAA	CTCGAC	CTCGAG	CTCGAT	CTCGCA	CTCGCC	CTCGCT	CTCGGA	CTCGGC	CTCGGG	CTCGGT	CTCGTA	CTCGTC	CTCGTG	CTCGTT	GACGAA	GACGAC	GACGAG	GACGAT	GACGCA	GACGCC	GACGCT	GACGGA	GACGGC	GACGGG	GACGGT	GACGTA	GACGTC	GACGTG	GACGTT	GCCGAA	GCCGAC	GCCGAG	GCCGAT	GCCGCA	GCCGCC	GCCGCT	GCCGGA	GCCGGC	GCCGGG	GCCGGT	GCCGTA	GCCGTC	GCCGTG	GCCGTT	GGCGAA	GGCGAC	GGCGAG	GGCGAT	GGCGCA	GGCGCC	GGCGCT1	GGCGGA	GGCGGC	GGCGGG	GGCGGT	GGCGTA	GGCGTC	GGCGTG	GGCGTT	GTCGAA	GTCGAC	GTCGAG	GTCGAT	GTCGCA	GTCGCC	GTCGCT	GTCGGA	GTCGGC	GTCGGG	GTCGGT	GTCGTA	GTCGTC	GTCGTG	GTCGTT	TACGAA	TACGAC	TACGAG	TACGAT	TACGCA	TACGCC	TACGCT	TACGGA	TACGGC	TACGGG	TACGGT	TACGTA	TACGTC	TACGTG	TACGTT	TCCGAA	TCCGAC	TCCGAG	TCCGAT	TCCGCA	TCCGCC	TCCGCT	TCCGGA	TCCGGC	TCCGGG	TCCGGT	TCCGTA	TCCGTC	TCCGTG	TCCGTT	TGCGAA	TGCGAC	TGCGAG	TGCGAT	TGCGCA	TGCGCC	TGCGCT	TGCGGA	TGCGGC	TGCGGG	TGCGGT	TGCGTA	TGCGTC	TGCGTG	TGCGTT	TTCGAA	TTCGAC	TTCGAG	TTCGAT	TTCGCA	TTCGCC	TTCGCT	TTCGGA	TTCGGC	TTCGGG	TTCGGT	TTCGTA	TTCGTC	TTCGTG	TTCGTT" >> $FILE_BASE_NAME".master.txt"


printf	"\n" >> $FILE_BASE_NAME".master.txt"

# Merge tables using coloumn 1 as key, and print as earch species bin per row. Columns increase by motif array size.
printf "Merging tables\n"

perl ~/bin/merge.pl -k -e NA $FILE_BASE_NAME".read_number.temp" $FILE_BASE_NAME".contigs_gc.temp" $FILE_BASE_NAME".contigs_cpg_sites.temp" $FILE_BASE_NAME".contigs_cpg_motifs_all.temp" >> $FILE_BASE_NAME".master.txt"

# Moving to more later version with complete scan of all 6-mer motifs
mv $FILE_BASE_NAME".master.txt" $FILE_BASE_NAME".spades_master_225.txt"

# Main tab delim output file
printf "Output table: $FILE_BASE_NAME".spades_master_225.txt""

## Clean up remaining temp files - end of script.
rm *.temp

printf "\n\nComplete\n\n"

## End
