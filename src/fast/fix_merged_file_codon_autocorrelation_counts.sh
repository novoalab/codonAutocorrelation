#!/bin/bash

# What for: Fix the "merged" datasets that may not have 197 lines (some pairs missing in the case of individual seq analysis).

if [[ $1 == "-h" ]];then
	echo -e "Usage: $0 <merged_codon_file>"
	echo -e "Uses file $HOME/scripts/ADAT_codons_combinations_with_zeroes.txt as reference set"
	exit
fi

# Take input
input_file=$1
reference_file="$HOME/scripts/ADAT_codons_combinations_with_zeroes.txt"

# Compare two files
awk {'print $1"_"$2"_"$3'} $input_file > tmp_input_file.txt
awk {'print $1"_"$2"_"$3'} $reference_file > tmp_reference_file.txt
~/scripts/compare_lists.py tmp_reference_file.txt tmp_input_file.txt  # outputs common.txt diff.txt

# Build merged file with missing 
cat diff.txt | sed 's/ /_/g'| while read line; do cat $reference_file |sed 's/ /_/g'| grep $line >> diff_tmp.txt; done
cat diff_tmp.txt | sed 's/_/ /g' > diff_final.txt

cat $input_file diff_final.txt | sort -g | awk {'print $1,$2,$3,$4,$5,$6'}|  sort -g  > $1.FIXED

# Remove unnecessary files
rm common.txt
rm tmp_input_file.txt
rm tmp_reference_file.txt
rm diff_final.txt
rm diff.txt
rm diff_tmp.txt


