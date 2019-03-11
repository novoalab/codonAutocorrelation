#!/bin/bash -l

###########################################################################
## MAIN SCRIPT TO RUN CODON AUTOCORRELATION
# Eva Maria Novoa, Feb 2019
# Needs as input a fasta file (can be multifasta) containing CDS sequences 
# i.e. each sequence should start with ATG and end with ter codon
###########################################################################


## READ INPUT
input_fasta=$1

## STEPS 
# 1. Build total_codons and paired_codons from fasta
echo "##" 
echo "## STEP 1: Extracting total codons and paired codons..."
echo "##" 
src/fast/codon_autocorrelation_multiple_sequences.py $input_fasta DNA

# 2. Build merged from total_codons and paired_codons
echo "##" 
echo "## STEP 2: Merging total_codons and paired_codons files..."
echo "##" 
src/fast/parse_codon_autocorrelation_output.sh $input_fasta.paired_codons $input_fasta.total_codons

# 3. Fix merged with missing pairs
echo "##" 
echo "## STEP 3: Add missing codon pairs to merged file..."
echo "##" 
src/fast/fix_merged_file_codon_autocorrelation_counts.sh $input_fasta.codon_counts.merged;                                                                     

# 4. Build R code with specific input and get SDEVS and plots
echo "##" 
echo "## STEP 4: Computing codon autocorrelation and RSCPU values..."
cat src/fast/parse_codon_autocorrelation_merged.R | sed "s/INPUTFILE/$input_fasta.codon_counts.merged.FIXED/g" > $input_fasta.codon_counts.merged.FIXED.R;
echo "##" 
R --save < $input_fasta.codon_counts.merged.FIXED.R

## FINISH
echo "##" 
echo "## STEP 5: Cleaning up..."
echo "##" 
echo
mv $input_fasta.* results
echo "--> Done! :)"
echo "--> You will find the results in the 'results' folder. Enjoy!"


