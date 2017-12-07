// Mapping PE fastq to de novo assembled contigs.
// Input: [REF] [PE1] [PE2]

@transform("sam")
align={
	doc "Aligning reads to single ref using bwa mem"
	exec "bwa mem -t 7 $input1 $input2 $input3 > $output.sam"
}

@transform("bam")
sam2SortedBam={
	exec """
		java -Xmx2g -jar ~/programs/picard-tools-1.119/SortSam.jar
		VALIDATION_STRINGENCY=LENIENT
		INPUT=$input.sam
		OUTPUT=$output.bam
		SORT_ORDER=coordinate
	"""
	exec "rm $input.sam"
}

indexBam = {
	transform("bam") to ("bam.bai") 
	{
		exec "samtools index $input.bam"
	}
}

// Remove potential reads mapping to multiple positions. Removes 0 quality alignments scores which bwa mem uses to assign reads with multiple primary mappings
multiMapReadsFilter={
	transform("bam") to ("unique.alignment.sam")
		{
			exec "samtools view -q1 -h $input.bam > $output.sam"	
		}
		transform("unique.alignment.sam") to ("unique.alignment.bam")
		{
			exec "samtools view -bS $input.sam > $output.bam"
			exec "rm $input.sam"
		}
}

indexUniqueBam={
	transform("unique.alignment.bam") to ("unique.alignment.bam.bai")
	{
		exec "samtools index $input.bam"
	}
}

mappingReport= {
	transform("bam") to ("bam.report.txt")
	{
		exec "~/programs/bamtools/bin/bamtools stats -in $input.bam > $output.txt"
	}
}

coverageReport = {
	transform("bam") to ("bam.genome_coverage.txt")
	{
		exec "bash ~/graham/pipelines/shell/coverageStats.sh $input.bam"
	}
}

// Final combined report currently run as a separate script: parse_bwa_mapping_single_ref.sh


// Run
Bpipe.run {
 align + sam2SortedBam + indexBam + multiMapReadsFilter + indexUniqueBam + mappingReport + coverageReport
}
