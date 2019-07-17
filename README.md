# codonAutocorrelation
Compute codon autocorrelation (aka codon covariation, codon reuse, codon pair usage) from fasta sequences.

This code was used in the analyses of the manuscript *'Elucidation of codon usage signatures across the domains of life'* (Novoa et al., Mol Biol Evol 2019), available here: https://doi.org/10.1093/molbev/msz124

## What is codon autocorrelation? 
Codon autocorrelation reflects the non-random distribution of consecutive codon occurrences throughout a transcript. Previous studies in yeast have shown that once a particular codon has been used, subsequent occurrences of the same amino acid in the same transcript are not random (Cannarozzi, et al. 2010), a phenomenon termed as ‘codon autocorrelation’ or ‘codon covariation’. Mechanistically, it was argued that tRNA recycling was the driving force causing the observed biased distribution of synonymous codons along a sequence, i.e., codons that would reuse the same tRNA would be favored as a means to increase the speed of translation (Cannarozzi, et al. 2010). A subsequent study re-examined this question, and compared the autocorrelation between codons encoding the same amino acids to those encoding different ones (Hussmann and Press 2014). Intriguingly, this second study found that covariation between codons encoding different amino acids was as strong as covariation between codons encoding the same amino acid, concluding that there was insufficient evidence to claim that tRNA recycling is the force driving codon autocorrelation. Despite the uncertain cause of codon covariation, both studies show that the probability of observing a specific codon is dependent on previous codon occurrences.

## What's included
- Scripts to compute and analyze codon autocorrelation from fasta sequences 
- Output: 3 types of measurements: 
  * *SDEVS*: number of standard deviations from expected, as defined in Cannarozzi et al, Cell 2010 (https://www.cell.com/abstract/S0092-8674(10)00189-3) 
  --> output file: *.sdevs.txt
  * *EFFECT MATRIX*: normalized sdevs, for each amino acid 
  --> output files: *effect_matrix.txt  
  * *RSCPU*: relative synonymous codon pair usage, as defined in Novoa et al., bioRxiv 2018. (https://doi.org/10.1101/421487) 
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
It is important to note that this mode will also consider consecutive codons (i.e. dicodon frequencies), which are excluded in the case of the slower version, to avoid biases that exist naturally as "dicodon".  Results will be placed in 'results' folder.

## Comprehensive start
This version *excludes* consecutive codons (i.e. dicodons) from the analysis. Furthermore this version also computes codon autocorrelation between codons that belong to different amino acids. For these two reasons, this code is much slower than the "quick start" one. In our hands, we have found that excluding consecutive codons does not significantly vary the results (at least in the few species that we have checked this), however, this may not be the case for all species. 
Overall, the fast version (i.e. fromFasta2CodonAutocorrelationFAST.sh) is recommended if only codon autocorrelation across codons that encode for the same amino acid is needed.

Usage: 
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

If you find this code useful, please cite: Novoa EM, Jungreis I, Jaillon O, Kellis M. Elucidation of codon usage signatures across the domains of life. Mol Biol Evol 2019. https://doi.org/10.1093/molbev/msz1247

## Contact

If you have any doubts/questions/concerns, please contact: evamaria.novoa@gmail.com. Thanks!

