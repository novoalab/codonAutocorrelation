#!/bin/bash

if [[ $1 == "-h" ]];then
	echo -e "Usage: $0 <paired_codons_file> <total_codons_file>"
	echo -e "\tBoth infiles have been generated with the script codon_autocorrelation_multiple_sequences.py"
	echo -e "\tOutput merges both files to have the individual values"
	exit
fi


# Input 
paired_codons_file=$1
total_codons_file=$2
name=`echo $1 | sed 's/\.paired_codons//'`


# Merge files
echo -e "Codon1\tCodon2\tPaired_total\tTotal_codon1\tTotal_codon2" > $name.codon_counts.merged
cat $paired_codons_file | while read line; do
	codon1=`echo $line | awk {'print $1'}`
	total_codon1=`cat $total_codons_file | grep $codon1 | awk {'print $2'}`
	codon2=`echo $line | awk {'print $2'}`
	total_codon2=`cat $total_codons_file | grep $codon2 | awk {'print $2'}`
	echo -e "$line\t$total_codon1\t$total_codon2" >> $name.codon_counts.merged
done
	
