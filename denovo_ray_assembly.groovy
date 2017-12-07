// Denovo assembly pipeline using RayMeta with downstream annotation of
// contigs with blastn, bwa and kronaTools. Outputs classifications and
// abundances of assembled contigs
// 
// Input PE fastq files of unmapped reads from primary reference based mapping

// Denovo assembly

rayMetaAssemble={
		doc "Assembling PE fastq reads with RayMeta using kmer length 31 (from std K41)"
		output.dir="ray_output"
		produce("Contigs.fasta")
		{
			exec "mpiexec -n 4 Ray -k 31 -p $input1.fastq $input2.fastq -output ray_temp"
			exec "cp -r ray_temp/* $output.dir"
			exec "rm -r ray_temp"
		}
}

cp={
	doc "Set up new dir"
	output.dir = "mapped_contigs_output"
	exec "cp $input.fasta $output.raymeta.fasta"
}

indexFasta={
		doc "Create new mapped contigs output folder. Index assembled contigs reads for bwa mapping"
		output.dir = "mapped_contigs_output"
		exec "bwa index $input.raymeta.fasta"
}


// Mapping of reads back onto assembled contigs

@transform("sam")
align={
	output.dir = "mapped_contigs_output"
	transform("fasta") to ("sam")
	{
		doc "Aligning input unmapped fastq reads to assembled contigs with bwa-mem"
		exec "bwa mem -t 7 -M $input.raymeta.fasta $input1.fastq $input2.fastq > $output.sam"
	}
}


@transform("bam")
sam2SortedBam={
		doc "Picard sam to sorted bam"
		output.dir = "mapped_contigs_output"
		exec """
			java -Xmx2g -jar ~/programs/picard-tools-1.119/SortSam.jar
			VALIDATION_STRINGENCY=LENIENT
			INPUT=$input.sam
			OUTPUT=$output.bam
			SORT_ORDER=coordinate
			"""
		exec "rm $input.sam"
}


indexBam={
	output.dir = "mapped_contigs_output"
		doc "Indexing bam file"
		exec "samtools index $input.bam"
	
		forward input.bam
}


// File output reporting of assembly and mapping stages

mappingReport={
	output.dir = "mapped_contigs_output"
	transform("raymeta.bam") to ("raymeta.mapping_report.txt")
	{
		doc "Bamtools stats report"
		exec "~/programs/bamtools/bin/bamtools stats -in $input > $output"
	}
	forward input.bam
}

contigReport={
	output.dir = "mapped_contigs_output"
	transform("raymeta.bam") to ("raymeta.contigs_report.txt")
	{
		doc "Calculating contig mapping stats using shell script: assemblyRayMetaCoverageStats.sh"
		exec "bash ~/pipelines/shell/assemblyRayMetaCoverageStats.sh $input > $output"
	}
}


// Blastn of assembled contigs against ncbi nt with xml output format (outfmt 5). Use xml as base format
contigBlast={
	output.dir = "blast_output"
	transform("raymeta.fasta") to ("raymeta.blastn.nt.xml.out")
	{
		exec "~/programs/ncbi-blast-2.2.29+/bin/blastn -outfmt 5 -query $input.raymeta.fasta -db nt -num_threads 8 -out $output"
	}
}



// Convert blast xml to tab format (outfmt 6) for krona input
convertContigBlast={
	output.dir = "blast_output"
	transform("raymeta.blastn.nt.xml.out") to ("raymeta.blastn.nt.tab.out")
	{
		exec "python ~/bin/blastxml_to_tabular.py -c std $input -o $output"
	}
}

// Filter BLASTn report. Use parameters: pct_id >= 90, eval >= 1e-8, bit_score >= 50 (as megan5 bit cutoff)
filterContigBlast={
	output.dir = "blast_output"
	transform("raymeta.blastn.nt.tab.out") to ("raymeta.blastn.nt.tab.filtered.out")
	{
		doc ""
		exec "awk '\$3 >= 90 && \$11 <= 1e-8 && \$12 >= 50 {print}' $input > $output"
	}

}

// Parse contigReport for krona magnitudes import. Uses query_id and foldX_read_depth column
readyMagnitudes={
	output.dir = "mapped_contigs_output"
	transform("raymeta.contigs_report.txt") to ("raymeta.contigs.krona_magnitudes.txt")
	{
		doc ""
		exec "awk '{print \$1 \"\t\" \$11}' $input > $output"
	}
}


// ktools visulisation
ktImportBlast={
	output.dir = "blast_output"
	transform("raymeta.blastn.nt.tab.filtered.out") to ("raymeta.blastn.nt.tab.filtered.krona.html")
	{
		doc "Convert blast results in krona plot, using LCA option and magnitudes"
		exec """
		ktImportBLAST $input:$input.cp.raymeta.contigs.krona_magnitudes.txt -o $output
		"""
	}
}

// ktools classify blast as above but output tab file
ktClassifyBlast={
	output.dir = "blast_output"
	transform("raymeta.blastn.nt.tab.filtered.out") to ("raymeta.blastn.nt.filtered.taxonomy.tab")
	{
			exec "ktClassifyBLAST $input -o $output"
	}
}

// Use ktools tab file and bioperl to generate taxon report
taxonSummaryReport={
	output.dir = "taxon_report"
	transform("raymeta.blastn.nt.filtered.taxonomy.tab") to ("raymeta.blastn.nt.taxonomy_report.txt")
	{
		exec "perl ~/bin/taxonId2name.pl $input > $output"
	}
}

// Merge taxonSummaryReport output with contigReport output. Write NA on keys with no records (contigs without blast hits)
taxonMappingSummaryReport={
	output.dir = "taxon_report"
	transform("raymeta.blastn.nt.taxonomy_report.txt") to ("rayameta.blastn.nt.taxon_mapping_report.txt")
	{
		exec "~/bin/merge.pl -k -e NA $input.cp.raymeta.contigs_report.txt $input.cp.raymeta.blastn.nt.taxonomy_report.txt > $output"
	}

}



// Run
Bpipe.run{ 
		rayMetaAssemble + cp + indexFasta + align + sam2SortedBam + indexBam + mappingReport + contigReport + contigBlast + convertContigBlast + filterContigBlast + readyMagnitudes + ktImportBlast + ktClassifyBlast + taxonSummaryReport + taxonMappingSummaryReport
	}
