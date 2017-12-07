#¬/bin/bash


#  Pipeline for the samtools workflow for calling variants.
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


if [ $# -ne 3 ] ; then
    echo 'Input: [run mode: single/multi] [reference.fasta] [sample.bam]'
    exit 0
fi

READGROUP_SMTAG=`echo $BAM | sed 's/_L001_R1.*.bam//'`

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
	# Samtools Variant Calling
	
	echo "running mpileup..."

	$SAMTOOLS mpileup -ugf $REFERENCE "${BAM%.bam}_realigned_reads.bam" | bcftools call -vmO v -o "${BAM%.bam}_realigned_reads.vcf"

	echo "tabix..."

	/usr/local/bin/tabix -p vcf "${BAM%.bam}_realigned_reads.vcf"
	
	echo "plotting snps..."

	bcftools stats -F $REFERENCE -s - "${BAM%.bam}_realigned_reads.vcf" > "${BAM%.bam}_realigned_reads.vcf.stats"
	mkdir plots
	plot-vcfstats -p plots/ "${BAM%.bam}_realigned_reads.vcf.stats"

	# Filtering
	bcftools filter -O v -o "${BAM%.bam}_realigned_reads.filtered.vcf" -s LOWQUAL -i'%QUAL>30 && DP>5' "${BAM%.bam}_realigned_reads.vcf"

	
	# Convert to fasta file with gentypes inserted
	cat  | bcftools convert -O b -o $file.bcf
	
	echo "done..."

fi
