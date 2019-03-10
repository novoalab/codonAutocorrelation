# codonAutocorrelation
Compute codon autocorrelation (aka codon covariation, codon reuse, codon pair usage) from fasta sequences

## What's included
- Scripts to compute and analyze codon autocorrelation from fasta sequences 
- Output: 3 types of measurements: 
* (i) SDEVS: number of standard deviations from expected, as defined in Cannarozzi et al, Cell 2010 (https://www.cell.com/abstract/S0092-8674(10)00189-3) 
  --> output file: *.sdevs.txt
  * (ii) EFFECT MATRIX: normalized sdevs, for each amino acid 
  --> output files: *effect_matrix.txt  
  * (iii) RSCPU: relative synonymous codon pair usage, as defined in Novoa et al., bioRxiv 2018. (https://doi.org/10.1101/421487) 
   --> output files *rscpu_way1.txt and *rscpu_way2.txt

## Quick start
If you want to go directly from a fasta to relative codon pair usage, you can directly run:
```
fromFasta2CodonAutocorrelationFAST.sh <FILE.fasta>
```
Example: 
```
fromFasta2CodonAutocorrelationFAST.sh Saccharomyces_cerevisiae.CDS.fasta
```
Results will be placed in 'results' folder

## Comprehensive start
To obtain ALL vs ALL codons, to compute codon autocorrelation within codons that encode for same amino acid, but also across amino acids, a much slower version of this script should be used:
```
fromFasta2CodonAutocorrelation.sh <FILE.fasta>
```
Example: 
```
fromFasta2CodonAutocorrelation.sh Saccharomyces_cerevisiae.CDS.fasta
```
Results will be placed in 'results' folder

## How to run the software: step by step

### STEP 1: Download genome CDS fasta sequences
You can download them, for example, from the EMBL CDS database: ftp://ftp.ebi.ac.uk/

### STEP 2: Get frequencies of codon pairs 

```
codon_autocorrelation_multiple_sequences.py <FILE.fasta> 
```
This will generate two files: <INPUT>.total_codons an <INPUT>.paired_codons

### STEP 3: Merge two files from step 2
```
parse_codon_autocorrelation_output.sh <FILE.paired_codons> <FILE.total_codons>

```
This will generate a merged file: <INPUT>.merged

### STEP 4: Fix missing values for codon pairs that did not exist
```
fix_merged_file_codon_autocorrelation_counts.sh <FILE.merged>
```
This will generate a merged file: <INPUT>.FIXED

### STEP 5: Compute codon pair usage statistics, relative to individual codon usage
```
  parse_codon_autocorrelation_merged.R --save < FILE.FIXED
```

## Future work
- Code can be modified to allow for the calculation of codon autocorrelation among non-synonymous codons (ongoing work)
- Script that computes all VS all codons needs further tuning in terms of speed (ongoing work)

## Citing this work

If you find this code useful, please cite: Novoa EM, Jungreis I, Jaillon O, Kellis M. Elucidation of codon usage signatures across the domains of life. bioRxiv 2018. https://doi.org/10.1101/421487

## Contact

If you have any doubts/questions/concerns, please contact: eva.novoa@crg.eu. Thanks!

