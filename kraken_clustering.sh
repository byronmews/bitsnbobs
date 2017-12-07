#!/bin/bash

#
# Feeding output from kraken.report into metaphian utils.
# Uses hclust etc for sample visualisation
#

# Normalise kraken map report counts by library size
awk '{print $1,(($2/[library_size])*1000000)}' [mpa_report]

# Build matrix from kraken mpa report format
~/programs/nsegata-metaphlan/utils/merge_metaphlan_tables.py [normalised_mpa_reports]

# Format matrix
sed 's/ /\t/g' [mpa_report_matrix]

# Plot matrix
~/programs/nsegata-metaphlan/plotting_scripts/metaphlan_hclust_heatmap.py --top 40 --minv 0.1 --tax_lev s -s log --in [mpa_report_matrix] --out out.png
