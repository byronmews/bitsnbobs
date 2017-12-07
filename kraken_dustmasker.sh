#!/bin/sh

#
#

# Run Dustmakser on each genome and gene file replacing the genome file  
for file in `find -name '*.fna' -o -name '*.fa'`
do
	# Check if run already
	if [[ $file == *.masked.* ]]
	then
		echo "Files already dusted, overwriting existing ones"
		exit;
	else
	
		echo "Running: $file"
    		dustmasker -in $file -infmt fasta -outfmt fasta | sed -e '/>/!s/a\|t\|g\|c/N/g' > "${file%.fna}.masked.fna"
	fi
done

