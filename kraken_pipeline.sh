#!/bin/bash

#
# Kraken kmer mapping of reads
#
# GRU01 version

# Input fastq
if [ $# -ne 2 ] ; then
    echo 'Input paired fastq files: R1 R2'
    exit 0
fi

base=`echo $1 | sed 's/.fastq//'`
krakenDB=/srv/data0/dbs/kraken/standard_210317

# Run kraken with PE reads joined with N linkers.
kraken \
--db $krakenDB \
--threads 12 \
--fastq-input \
--paired \
$1 $2 \
--output $base".kraken.out"

# Report metaphlan format pre filter, ready for krona input
kraken-mpa-report \
--db $krakenDB $base".kraken.out" \
> $base".kraken.mpa_report"

# Convert to ktTools ready format and import taxonomy for krona visual
cut -f2,3 $base".kraken.out" | ktImportTaxonomy - -o $base".kraken.krona.html"

#python /home/graham/programs/nsegata-metaphlan/conversion_scripts/metaphlan2krona.py \
#-p $base".kraken.mpa_report" \
#-k $base".kraken.krona"
