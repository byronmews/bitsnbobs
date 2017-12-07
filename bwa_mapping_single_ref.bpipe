// Mapping PE fastq to a single reference.
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

mappingReport= {
	transform("bam") to ("bam.report.txt")
	{
		exec "~/programs/bamtools/bin/bamtools stats -in $input.bam > $output.txt"
	}
}

coverageReport = {
	transform("bam") to ("bam.genome_coverage.txt")
	{
		exec "bash ~/pipelines/shell/coverageStats.sh $input.bam"
	}
}

// Final combined report currently run as a separate script: parse_bwa_mapping_single_ref.sh


// Run
Bpipe.run {
	indexBam + mappingReport + coverageReport
}
