#!/bin/sh

#  GATK.sh
#  Pipeline for the GATK workflow for calling variants.
#  Uses HaplotypeCaller.
#



##################################################################
#   SOFTWARES
# Provide locations of the softwares to be used. 

PICARD="java -jar /usr/local/src/picard-tools-1.119"
SAMTOOLS="samtools"
GATK="java -jar /usr/local/src/gatk/GenomeAnalysisTK.jar"


#################################################################
#   FILES 
# Aligned BAM file, reference FASTA file and VCF file of known variants. 

MODE=$1
REFERENCE=$2
BAM=$3
READGROUP_SMTAG=$4


if [ $# -ne 4 ] ; then
    echo 'Input: [run mode: single/multi] [reference.fasta] [sample.bam] [sample_id] '
    exit 0
fi



######################################################################
## If MODE flag used, run as single or multimode
if [ $MODE == "single" ]
then

	#####################################################################
	#   Check if the BAM file satisfies the requirements of the GATK

	$SAMTOOLS view -H $BAM > $READGROUP_SMTAG"_alpha.txt"
	head_file=$READGROUP_SMTAG"_alpha.txt"
	grep '^@SQ' $head_file > $READGROUP_SMTAG"_beta.txt"
	file=$READGROUP_SMTAG"_beta.txt"

	echo "Created $file"

	#less $file
	# check if the file is sorted.
	T="sort -c -t':' -nk2 $file"
	if [ "$T" ]; then
    	echo "The file is sorted!"

	else
    	echo "The file is not sorted."
    	echo "Sorting........."
    	$PICARD/SortSam.jar INPUT=$BAM OUTPUT="${BAM%.bam}.sorted.bam" SORT_ORDER=coordinate

	fi


#   check if the file contains RG information.
	if grep -q '^@RG' $head_file; then
		echo "The file contains RG information."
	else
		echo "The file does not contain RG information and the GATK will not work!"
		echo "Attempting fix..."
	
		cp $BAM $BAM".bkup"	

		$PICARD/AddOrReplaceReadGroups.jar \
		INPUT= $BAM \
		OUTPUT=$BAM".temp" RGID=1 RGLB=Library1 RGPL=ILLUMINA RGPU=1 RGSM=$READGROUP_SMTAG RGCN=AA RGDS=AD

		# Renaming BAM variable to new RG BAM
		mv $BAM".temp" $BAM
	fi


	###################################################################
	#   Prepping reference geonome fasta file for GATK
	echo "Prepping reference geonome fasta file for GATK....."
	
	# Create sequence dictionary using Picard Tools.
	# the following command produces a SAM-style header file describing the contents of our fasta file.
	$PICARD/CreateSequenceDictionary.jar \
	REFERENCE=$REFERENCE \
	OUTPUT=$REFERENCE".dict"

	# Looking for reference without fasta suffix
	NEWLIB=`echo $REFERENCE".dict" | sed 's/.fasta//'`
	mv $REFERENCE".dict" $NEWLIB

	echo "created sequence dictionary for the reference genome."
	echo "indexing the reference genome...."

	# Create the fasta index file.
	# The index file describes byte offset in the fasta file for each contig. It is a text file with one record
	# per line for each of the fasta contigs. Each record is of the type -
	# contig, size, location, basePerLine, bytesPerLine
	$SAMTOOLS faidx $REFERENCE

	echo "Reference genome is now ready for GATK."


	###############################################################
	## Summary Statistics

	$PICARD/MeanQualityByCycle.jar \
	INPUT=$BAM \
	CHART_OUTPUT=$READGROUP_SMTAG"_mean_quality_by_cycle.pdf" \
	OUTPUT=$READGROUP_SMTAG"_read_quality_by_cycle.txt" \
	REFERENCE_SEQUENCE=$REFERENCE


	$PICARD/QualityScoreDistribution.jar \
	INPUT=$BAM \
	CHART_OUTPUT=$READGROUP_SMTAG"_mean_quality_overall.pdf" \
	OUTPUT=$READGROUP_SMTAG"_read_quality_overall.txt" \
	REFERENCE_SEQUENCE=$REFERENCE


	$PICARD/CollectWgsMetrics.jar \
	INPUT=$BAM OUTPUT=$READGROUP_SMTAG"_stats_picard.txt" \
	REFERENCE_SEQUENCE=$REFERENCE \
	MINIMUM_MAPPING_QUALITY=20 \
	MINIMUM_BASE_QUALITY=20

	
	#############################################################
	# Mark duplicate reads.

	echo "mark the duplicates in the bam file."

	$PICARD/MarkDuplicates.jar INPUT=$BAM OUTPUT="${BAM%.bam}_dups_marked.bam" \
	METRICS_FILE="${BAM%.bam}_dups_metrics.txt" REMOVE_DUPLICATES=false

	echo "Index the dup-marked bam file,${BAM%.bam}_dups_marked.bam"

	$SAMTOOLS index "${BAM%.bam}_dups_marked.bam"


	BAM_FILE="${BAM%.bam}_dups_marked.bam"

	#############################################################
	## GATK Data Pre-Processing

	# Step 1 - Local realignment around indels.
	# Create a target list of intervals to be realigned.

	echo "Creating a target list of intervals to be realigned...."

	$GATK \
	-T RealignerTargetCreator \
	-R $REFERENCE \
	-I $BAM_FILE \
	-o "${BAM%.bam}_target_intervals.list"

	# do the local realignment.
	echo "local realignment..."

	$GATK \
	-T IndelRealigner \
	-R $REFERENCE \
	-I $BAM_FILE \
	-targetIntervals "${BAM%.bam}_target_intervals.list" \
	-o "${BAM%.bam}_realigned_reads.bam"

	echo "indexing the realigned bam file..."

	# Create a new index file.
	$SAMTOOLS index "${BAM%.bam}_realigned_reads.bam"

	# Step 2 - Base recalibration (fixes them so they better reflect the probability of mismatching the genome).
	# Analyze patterns of covariation in the sequence.
	echo "base recalibration...skipped"

	#$GATK \
	#-T BaseRecalibrator \
	#-R $REFERENCE \
	#-I "${BAM%.bam}_realigned_reads.bam" \
	#-o "${BAM%.bam}_recal_data.table"

	# Apply recalibration to the sequence data.

	echo "recalibrating the sequence data...skipped"

	#$GATK \
	#-T PrintReads \
	#-R $REFERENCE \
	#-I "${BAM%.bam}_realigned_reads.bam" \
	#-BQSR "${BAM%.bam}_recal_data.table" \
	#-o "${BAM%.bam}_recal_reads.bam"




	###########################################################################
	# GATK Variant Calling -  HaplotypeCaller

	minBaseScore=20	#Minimum Phred base score to count a base (20 = 0.01 error, 30=0.001 error, etc)
	stand_call_conf=30
	targets_interval_list="${BAM%.bam}_target_intervals.list"

	echo "calling variants...."

	$GATK \
	-T HaplotypeCaller \
	-R $REFERENCE \
	-I "${BAM%.bam}_realigned_reads.bam" \
	-mbq $minBaseScore \
	-stand_call_conf $stand_call_conf \
	-L $targets_interval_list \
	-nct 15 \
	--emitRefConfidence GVCF \
	-o "${BAM%.bam}_output.raw.snps.indels.g.vcf"

	echo "sample bam GVCF file created...."

	echo "run GenotypeGVCFs across all files next..."

	echo "eg. java -jar /usr/local/src/gatk/GenomeAnalysisTK.jar -T GenotypeGVCFs -R BK_polyomavirus.fasta -V FILE1 -V FILE2 -o out.vcf"

fi


######################################################################
## If MODE flag used, run as single or multimode
if [ $MODE == "multi" ]
then

	for file in *_output.raw.snps.indels.g.vcf
	do
		echo $file >> genotypeGVCFs.list
	done
	
	## GenotypeGVCFs
	$GATK \
	-T GenotypeGVCFs \
	-R $REFERENCE \
	--variant genotypeGVCFs.list \
	-o genotype.pooled.raw.snps.indels.vcf


	#############################################################
	## GenotypeGVCFs filtering
	#

	### SNPs
	## Variant filtration on SNPs. *Note complete dataset too small for truth filtering with VQSR
	java -jar /usr/local/src/gatk/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $REFERENCE \
	-V genotype.pooled.raw.snps.indels.vcf \
	-selectType SNP \
	-o genotype.pooled.raw.snps.vcf

	## Filter SNPs by hard set thresholds
	$GATK \
	-T VariantFiltration \
	-R $REFERENCE \
	--variant genotype.pooled.raw.snps.vcf \
	--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "gatk_snp_filter" \
	-o genotype.pooled.filtered.snps.vcf 
	
	### INDELS
	## Variant filtration on SNPs. *Note complete dataset too small for truth filtering with VQSR
	java -jar /usr/local/src/gatk/GenomeAnalysisTK.jar \
	-T SelectVariants \
	-R $REFERENCE \
	-V genotype.pooled.raw.snps.indels.vcf \
	-selectType INDEL \
	-o genotype.pooled.raw.indels.vcf
	
	## Filter Indels but hard set thresholds
	$GATK \
	-T VariantFiltration \
	-R $REFERENCE \
	--variant genotype.pooled.raw.indels.vcf \
	--filterExpression "QD < 2.0 || FS > 200.0 || MQ < 40.0 || ReadPosRankSum < -20.0" \
	--filterName "gatk_indel_filter" \
	-o genotype.pooled.filtered.indels.vcf


	#############################################################
	## Merge separate filtered snp and indel vcf
	$GATK \
	-T CombineVariants \
	-R BK_polyomavirus.fasta \
	--variant genotype.pooled.filtered.snps.vcf \
	--variant genotype.pooled.filtered.indels.vcf \
	-o genotype.pooled.filtered.snps.indels.vcf \
	--genotypemergeoption UNSORTED
	
fi

#############################################################
## Variants to table
#java -jar /usr/local/src/gatk/GenomeAnalysisTK.jar -T VariantsToTable -R CP010376.fasta -V E-cloacae_all_snps_filtered.vcf -F CHROM -F POS -F REF -F ALT -GF GT -F QUAL -F AC -F HET   -o E-cloacae_all_snps_filtered.csv


