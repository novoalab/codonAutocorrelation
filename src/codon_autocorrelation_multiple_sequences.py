#!/usr/bin/python

##################################
## CODON AUTOCORRELATION #########
## Eva Maria Novoa, March 2015  ##
##################################


import sys
import re


## USAGE

if sys.argv[1]=="-h" or len(sys.argv)<3:
	print "\nUsage:",sys.argv[0],"<file.fasta> <option>"
	print "\twhere file.fasta is a fasta file, which can have multiple sequences"
	print "\twhere <option> has to be DNA or RNA\n"
	sys.exit()
	
## SPLIT CODONS FUNCTION

def split_codons(str,num):
	return [str[start:start+num] for start in range(0, len(str), num) ]
	

## SELECT OPTION TO ASSIGN GENETIC CODE

option=sys.argv[2]

if option == "RNA":
	ala=['GCU','GCG','GCC','GCA']
	arg=['CGU','CGG','CGC','CGA','AGA','AGG']
	asn=['AAU','AAC']
	asp=['GAU','GAC']
	cys=['UGU','UGC']
	gln=['CAG','CAA']
	glu=['GAG','GAA']
	gly=['GGU','GGG','GGC','GGA']
	his=['CAU','CAC']
	ile=['AUU','AUC','AUA']
	leu=['CUU','UUG','CUG','CUC','UUA','CUA']
	lys=['AAG','AAA']
	met=['AUG']
	phe=['UUU','UUC']
	pro=['CCU','CCG','CCC','CCA']
	ser=['AGU','UCU','UCG','AGC','UCC','UCA']
	ter=['UAG','UGA','UAA']
	thr=['ACU','ACG','ACC','ACA']
	trp=['UGG']
	tyr=['UAU','UAC']
	val=['GUU','GUG','GUC','GUA']


elif option == "DNA":
	ala=['GCA','GCC','GCG','GCT']
	arg=['CGA','CGC','CGG','CGT','AGA','AGG']
	asn=['AAT','AAC']
	asp=['GAT','GAC']
	cys=['TGT','TGC']
	gln=['CAG','CAA']
	glu=['GAG','GAA']
	gly=['GGT','GGG','GGC','GGA']
	his=['CAT','CAC']
	ile=['ATT','ATC','ATA']
	leu=['CTT','TTG','CTG','CTC','TTA','CTA']
	lys=['AAG','AAA']
	met=['ATG']
	phe=['TTT','TTC']
	pro=['CCT','CCG','CCC','CCA']
	ser=['AGT','TCT','TCG','AGC','TCC','TCA']
	ter=['TAG','TGA','TAA']
	thr=['ACT','ACG','ACC','ACA']
	trp=['TGG']
	tyr=['TAT','TAC']
	val=['GTT','GTG','GTC','GTA']
	all_aa=['GCA','GCC','GCG','GCT','CGA','CGC','CGG','CGT','AGA','AGG','AAT','AAC','GAT','GAC','TGT','TGC','CAG','CAA','GAG','GAA','GGT','GGG','GGC','GGA','CAT','CAC','ATT','ATC','ATA','CTT','TTG','CTG','CTC','TTA','CTA','AAG','AAA','ATG','TTT','TTC','CCT','CCG','CCC','CCA','AGT','TCT','TCG','AGC','TCC','TCA','TAG','TGA','TAA','ACT','ACG','ACC','ACA','TGG','TAT','TAC','GTT','GTG','GTC','GTA']
	
	all_aa_list_of_lists=[ala,arg,asn,asp,cys,gln,glu,gly,his,ile,leu,lys,met,phe,pro,ser,thr,trp,tyr,val]
	#print all_aa_list_of_lists

	
else:
	print "\nUsage:",sys.argv[0],"<file.fasta> <option>"
	print "\twhere file.fasta is a fasta file, which can have multiple sequences"
	print "\twhere <option> has to be DNA or RNA\n"
	sys.exit()





## SCRIPT

# 1. Open fasta file
infile=open(sys.argv[1],'r')

# 2. Get list of codons of the fasta file
first=True
for i in infile:
	line=i.strip().replace("\n","")
	line=line.upper()
	if re.match(">",line):
		if first==True: # so first sequence
			seqCodonList=[]
			first=False
			myseq=""
		else:
			codonList=split_codons(myseq,3)
			seqCodonList.append(codonList)
			myseq=""
	else:
		myseq=myseq+line
codonList=split_codons(myseq,3)
seqCodonList.append(codonList)
			

# 3. Count pairs of codon occurrences for allcodons

# Build a list for codons for each PAIR of amino acids that we will find in the sequence 
part_PairList=[]
PairList=[]
# Output file
name=sys.argv[1]+".paired_codons"
outputfile2=open(name,'w')


# LOOP TO FILL CODON COMBS FOR EACH GROUP OF CODONS WITHIN 2 AA: 
for seq in range(len(seqCodonList)): # each sequence analyzed
	for aa1 in all_aa: # double loop to get two pairs of aa (with its codons, and build merged list
		for c in range(len(all_aa_list_of_lists)):
			for aa2 in all_aa_list_of_lists[c]:
				if aa1 in seqCodonList[seq] and aa2 in seqCodonList[seq]:
					#print aa1,aa2
					set=[]
					set.append(aa1)
					set=set+all_aa_list_of_lists[c]
					for i in range(len(seqCodonList[seq])): # the list of codons in the given sequence	
						if seqCodonList[seq][i] in set:
							part_PairList.append(seqCodonList[seq][i])	
					for k in range(len(part_PairList)-1):		
						if part_PairList[k]==aa1 and part_PairList[k+1]==aa2:
							#print part_PairList
							pair=(part_PairList[k],part_PairList[k+1])	
							PairList.append(pair)
					part_PairList=[]			

		

# Write file of pair co-occurrences	
#print "This is the whole PairList:",PairList
#PairList.sort()	
Seen=[]
for el in PairList:
	if el in Seen:
		pass
	else:
		Seen.append(el)
		text=str(el[0])+"\t"+str(el[1])+"\t"+str(PairList.count(el))+"\n"
		outputfile2.write(text)
			
# 4. Count occurrences for allcodons

# Build a list for codons for each amino acid that we will find in the sequence
alaList=[] 
part_alaList=[]

for seq in range(len(seqCodonList)): # each sequence analyzed
	for i in range(len(seqCodonList[seq])): # the list of codons in the given sequence
		for codon in all_aa: # codons of alanine amino acid
			if seqCodonList[seq][i]==codon:
				part_alaList.append(seqCodonList[seq][i])
	alaList.append(part_alaList)
	part_alaList=[]


# 5. Compute how many times each individual codon occurs

part_count_ala=[]
count_ala=[]
for seq in range(len(alaList)):
	for codon in all_aa:
		counter=alaList[seq].count(codon)
		part_count_ala.append(counter)
	count_ala.append(part_count_ala)
	part_count_ala=[]



# 6. Compute the SUMS of all individual codons for all proteins

# Output file 1 --> total number of individual codons
name=sys.argv[1]+".total_codons"
outputfile1=open(name,'w')

tmp_ala=0
total_ala=[]
for j in range(len(all_aa)):
	for i in range(len(count_ala)):
		tmp_ala=int(count_ala[i][j])+tmp_ala
	total_ala.append(tmp_ala)
	tmp_ala=0
for i in range(len(all_aa)):
	text=all_aa[i]+"\t"+str(total_ala[i])+"\n"
	outputfile1.write(text)

# 4. Check results in lists			
print "\nShort extract of results to validate the run\n"
print "\tCodons in first sequence",alaList[0]
print "\tCodons in last sequence",alaList[-1]
print "\tTotal codons in first sequence:",count_ala[0]
print "\tTotal codons in last sequence:",count_ala[-1],"\n"
print  "\tTotal codons summed from all sequences:",total_ala
print "\tPrinted out output file 1 named as",sys.argv[1]+".total_codons, containing total counts of codons "
print "\tPrinted out output file 2 named as",sys.argv[1]+".paired_codons, containing total pairs of codons\n\n"
