# codonAutocorrelation
From fasta, compute codon autocorrelation (aka codon covariation, codon reuse) 
relative to 


## What is included
- Sets of scripts to analyze codon autocorrelation in fasta sequences (please note, current implementation is limited to computing codon covariation ONLY within synonymous codons)
- Outputs 3 types of measurements: number of standard deviations from expected, effect vector and  relative synonymous codon pair usage (RSCPU)

## Quick start
If you want to go directly from a fasta to relative codon pair usage, you can directly run:
''' (either in the form of the Relative Synonymos Codon Pair Usage (RSCPU) and/or codon pair usage standard deviation from expected 

## How to run the software: step by step

### STEP 1: Download genome CDS fasta sequences
You can get them, for example, from EMBL CDS database: ftp://ftp.ebi.ac.uk/

### STEP 2: Get frequencies of codon pairs 

```
codon_autocorrelation_multiple_sequences.py <FILE.fasta> <mode>
where: <mode> can be DNA or RNA
```
This will generate two files: FILE.total_codons an FILE.paired_codons

### STEP 3: Merge two files from step 2
```
parse_codon_autocorrelation_output.sh <FILE.paired_codons> <FILE.total_codons>

```
This will generate a merged file: <INPUT>.merged

### STEP 4: Fix missing values for codon pairs that did not exist
```
fix_merged_file_codon_autocorrelation_counts_ALLaminoacids.sh <FILE.merged>
```
This will generate a merged file: <INPUT>.FIXED

### STEP 5: COMPUTE CODON STATISTICS OF PAIR USAGE, RELATIVE TO INDIVIDUAL CODON USAGE
```
  parse_codon_autocorrelation_merged.R --save <FILE.FIXED>
```

## Future work
Code can be modified to allow for the calculation of codon autocorrelation among non-synonymous codons (ongoing work)

## Citing this work

If you find this code useful, please cite: Novoa EM, Jungreis I, Jaillon O, Kellis M. Elucidation of codon usage signatures across the domains of life. bioRxiv 2018. https://doi.org/10.1101/421487

## Contact

If you have any doubts/questions/concerns, please contact: enovoa@mit.edu. Thanks!


