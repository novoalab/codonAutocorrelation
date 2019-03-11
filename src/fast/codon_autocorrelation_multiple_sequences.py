#!/usr/bin/python

##################################
## CODON AUTOCORRELATION #########
## Eva Maria Novoa, August 2013 ##
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
			

# 3. Count pairs of codon occurrences for a same amino acid

# Build a list for codons for each amino acid that we will find in the sequence (only has sense for 3,4,6 box codons)
alaList=[]
argList=[]
glyList=[]
ileList=[]
leuList=[]
proList=[]
serList=[]
thrList=[]
valList=[]
part_alaList=[]
part_argList=[]
part_glyList=[]
part_ileList=[]
part_leuList=[]
part_proList=[]
part_serList=[]
part_thrList=[]
part_valList=[]

for seq in range(len(seqCodonList)): # each sequence analyzed
	if seq==0:
		print seqCodonList[0]
	for i in range(len(seqCodonList[seq])-1): # the list of codons in the given sequence
		for ala_codon in ala: # codons of alanine amino acid
			if seqCodonList[seq][i]==ala_codon:
				part_alaList.append(seqCodonList[seq][i])
		for arg_codon in arg: 
			if seqCodonList[seq][i]==arg_codon:
				part_argList.append(seqCodonList[seq][i])
		for gly_codon in gly: 
			if seqCodonList[seq][i]==gly_codon:
				part_glyList.append(seqCodonList[seq][i])
		for ile_codon in ile: 
			if seqCodonList[seq][i]==ile_codon:
				part_ileList.append(seqCodonList[seq][i])
		for leu_codon in leu: 
			if seqCodonList[seq][i]==leu_codon:
				part_leuList.append(seqCodonList[seq][i])
		for pro_codon in pro: 
			if seqCodonList[seq][i]==pro_codon:
				part_proList.append(seqCodonList[seq][i])
		for ser_codon in ser: 
			if seqCodonList[seq][i]==ser_codon:
				part_serList.append(seqCodonList[seq][i])
		for thr_codon in thr: 
			if seqCodonList[seq][i]==thr_codon:
				part_thrList.append(seqCodonList[seq][i])
		for val_codon in val: 
			if seqCodonList[seq][i]==val_codon:
				part_valList.append(seqCodonList[seq][i])
	alaList.append(part_alaList)
	argList.append(part_argList)
	glyList.append(part_glyList)
	ileList.append(part_ileList)
	leuList.append(part_leuList)
	proList.append(part_proList)
	serList.append(part_serList)
	thrList.append(part_thrList)
	valList.append(part_valList)
	part_alaList=[]
	part_argList=[]
	part_glyList=[]
	part_ileList=[]
	part_leuList=[]
	part_proList=[]
	part_serList=[]
	part_thrList=[]
	part_valList=[]	
	
			
# 4. Check results in lists			
print "\nShort extract of results to validate the run\n"
print "\tAla codons in genetic code:",ala
print "\tAla codons in first sequence",alaList[0]
print "\tAla codons in last sequence",alaList[-1]

#print "arg codons:",arg
#print argList

#print "gly codons:",gly
#print glyList

#print "ile codons:",ile
#print ileList

#print "leu codons:",leu
#print leuList

#print "pro codons:",pro
#print proList

#print "ser codons:",ser
#print serList

#print "thr codons:",thr
#print thrList

#print "val codons:",val
#print valList


# 5. Compute how many times each individual codon occurs

#ala
part_count_ala=[]
count_ala=[]
for seq in range(len(alaList)):
	for codon in ala:
		counter=alaList[seq].count(codon)
		part_count_ala.append(counter)
	count_ala.append(part_count_ala)
	part_count_ala=[]
print "\tTotal Ala codons in first sequence:",count_ala[0]
print "\tTotal Ala codons in last sequence:",count_ala[-1],"\n"

#arg
part_count_arg=[]
count_arg=[]
for seq in range(len(argList)):
	for codon in arg:
		counter=argList[seq].count(codon)
		part_count_arg.append(counter)
	count_arg.append(part_count_arg)
	part_count_arg=[]
#print "Arg:",count_arg


#gly
part_count_gly=[]
count_gly=[]
for seq in range(len(glyList)):
	for codon in gly:
		counter=glyList[seq].count(codon)
		part_count_gly.append(counter)
	count_gly.append(part_count_gly)
	part_count_gly=[]
#print "Gly:",count_gly

#ile
part_count_ile=[]
count_ile=[]
for seq in range(len(ileList)):
	for codon in ile:
		counter=ileList[seq].count(codon)
		part_count_ile.append(counter)
	count_ile.append(part_count_ile)
	part_count_ile=[]
#print "Ile:",count_ile

#leu
part_count_leu=[]
count_leu=[]
for seq in range(len(leuList)):
	for codon in leu:
		counter=leuList[seq].count(codon)
		part_count_leu.append(counter)
	count_leu.append(part_count_leu)
	part_count_leu=[]
#print "Leu:",count_leu
	
#pro
part_count_pro=[]
count_pro=[]
for seq in range(len(proList)):
	for codon in pro:
		counter=proList[seq].count(codon)
		part_count_pro.append(counter)
	count_pro.append(part_count_pro)
	part_count_pro=[]
#print "Pro:",count_pro		
	
#ser
part_count_ser=[]
count_ser=[]
for seq in range(len(serList)):
	for codon in ser:
		counter=serList[seq].count(codon)
		part_count_ser.append(counter)
	count_ser.append(part_count_ser)
	part_count_ser=[]
#print "Ser:",count_ser
	
#thr
part_count_thr=[]
count_thr=[]
for seq in range(len(thrList)):
	for codon in thr:
		counter=thrList[seq].count(codon)
		part_count_thr.append(counter)
	count_thr.append(part_count_thr)
	part_count_thr=[]	
#print "Thr:",count_thr

#val
part_count_val=[]
count_val=[]
for seq in range(len(valList)):
	for codon in val:
		counter=valList[seq].count(codon)
		part_count_val.append(counter)
	count_val.append(part_count_val)
	part_count_val=[]
#print "Val:",count_val

# 6. Compute the SUMS of all individual codons for all proteins

# Output file 1 --> total number of individual codons
name=sys.argv[1]+".total_codons"
outputfile1=open(name,'w')

#ala
tmp_ala=0
total_ala=[]
for j in range(len(ala)):
	for i in range(len(count_ala)):
		tmp_ala=int(count_ala[i][j])+tmp_ala
	total_ala.append(tmp_ala)
	tmp_ala=0
for i in range(len(ala)):
	text="Ala\t"+ala[i]+"\t"+str(total_ala[i])+"\n"
	outputfile1.write(text)


#arg
tmp_arg=0
total_arg=[]
for j in range(len(arg)):
	for i in range(len(count_arg)):
		tmp_arg=int(count_arg[i][j])+tmp_arg
	total_arg.append(tmp_arg)
	tmp_arg=0
for i in range(len(arg)):
	text="Arg\t"+arg[i]+"\t"+str(total_arg[i])+"\n"
	outputfile1.write(text)

#gly
tmp_gly=0
total_gly=[]
for j in range(len(gly)):
	for i in range(len(count_gly)):
		tmp_gly=int(count_gly[i][j])+tmp_gly
	total_gly.append(tmp_gly)
	tmp_gly=0
for i in range(len(gly)):
	text="Gly\t"+gly[i]+"\t"+str(total_gly[i])+"\n"
	outputfile1.write(text)

#ile
tmp_ile=0
total_ile=[]
for j in range(len(ile)):
	for i in range(len(count_ile)):
		tmp_ile=int(count_ile[i][j])+tmp_ile
	total_ile.append(tmp_ile)
	tmp_ile=0
for i in range(len(ile)):
	text="Ile\t"+ile[i]+"\t"+str(total_ile[i])+"\n"
	outputfile1.write(text)

#leu
tmp_leu=0
total_leu=[]
for j in range(len(leu)):
	for i in range(len(count_leu)):
		tmp_leu=int(count_leu[i][j])+tmp_leu
	total_leu.append(tmp_leu)
	tmp_leu=0
for i in range(len(leu)):
	text="Leu\t"+leu[i]+"\t"+str(total_leu[i])+"\n"
	outputfile1.write(text)
#pro
tmp_pro=0
total_pro=[]
for j in range(len(pro)):
	for i in range(len(count_pro)):
		tmp_pro=int(count_pro[i][j])+tmp_pro
	total_pro.append(tmp_pro)
	tmp_pro=0
for i in range(len(pro)):
	text="Pro\t"+pro[i]+"\t"+str(total_pro[i])+"\n"
	outputfile1.write(text)

#ser
tmp_ser=0
total_ser=[]
for j in range(len(ser)):
	for i in range(len(count_ser)):
		tmp_ser=int(count_ser[i][j])+tmp_ser
	total_ser.append(tmp_ser)
	tmp_ser=0
for i in range(len(ser)):
	text="Ser\t"+ser[i]+"\t"+str(total_ser[i])+"\n"
	outputfile1.write(text)
#thr
tmp_thr=0
total_thr=[]
for j in range(len(thr)):
	for i in range(len(count_thr)):
		tmp_thr=int(count_thr[i][j])+tmp_thr
	total_thr.append(tmp_thr)
	tmp_thr=0
for i in range(len(thr)):
	text="Thr\t"+thr[i]+"\t"+str(total_thr[i])+"\n"
	outputfile1.write(text)
#val
tmp_val=0
total_val=[]
for j in range(len(val)):
	for i in range(len(count_val)):
		tmp_val=int(count_val[i][j])+tmp_val
	total_val.append(tmp_val)
	tmp_val=0
for i in range(len(val)):
	text="Val\t"+val[i]+"\t"+str(total_val[i])+"\n"
	outputfile1.write(text)
# 7. Loops around each amino acid to count for pairs 

# Output file 2 --> total number of paired codons
name=sys.argv[1]+".paired_codons"
outputfile2=open(name,'w')

#Ala
alaPairList=[]
for seq in range(len(alaList)):
	for i in range(len(alaList[seq])-1):
		pair=(alaList[seq][i],alaList[seq][i+1])
		alaPairList.append(pair)		
alaPairList.sort()
alaSeen=[]
for el in alaPairList:
	if el in alaSeen:
		pass
	else:
		alaSeen.append(el)
		text="Ala\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(alaPairList.count(el))+"\n"
		outputfile2.write(text)
#Arg	
argPairList=[]
for seq in range(len(argList)):
	for i in range(len(argList[seq])-1):
		pair=(argList[seq][i],argList[seq][i+1])
		argPairList.append(pair)
argPairList.sort()
argSeen=[]
for el in argPairList:
	if el in argSeen:
		pass
	else:
		argSeen.append(el)
		text="Arg\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(argPairList.count(el))+"\n"
		outputfile2.write(text)	
		
#Gly	
glyPairList=[]
for seq in range(len(glyList)):
	for i in range(len(glyList[seq])-1):
		pair=(glyList[seq][i],glyList[seq][i+1])
		glyPairList.append(pair)
glyPairList.sort()
glySeen=[]
for el in glyPairList:
	if el in glySeen:
		pass
	else:
		glySeen.append(el)
		text="Gly\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(glyPairList.count(el))+"\n"
		outputfile2.write(text)

#Ile	
ilePairList=[]
for seq in range(len(ileList)):
	for i in range(len(ileList[seq])-1):
		pair=(ileList[seq][i],ileList[seq][i+1])
		ilePairList.append(pair)
ilePairList.sort()
ileSeen=[]
for el in ilePairList:
	if el in ileSeen:
		pass
	else:
		ileSeen.append(el)
		text="Ile\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(ilePairList.count(el))+"\n"
		outputfile2.write(text)

#Leu	
leuPairList=[]
for seq in range(len(leuList)):
	for i in range(len(leuList[seq])-1):
		pair=(leuList[seq][i],leuList[seq][i+1])
		leuPairList.append(pair)
leuPairList.sort()
leuSeen=[]
for el in leuPairList:
	if el in leuSeen:
		pass
	else:
		leuSeen.append(el)
		text="Leu\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(leuPairList.count(el))+"\n"
		outputfile2.write(text)

#Pro	
proPairList=[]
for seq in range(len(proList)):
	for i in range(len(proList[seq])-1):
		pair=(proList[seq][i],proList[seq][i+1])
		proPairList.append(pair)
proPairList.sort()
proSeen=[]
for el in proPairList:
	if el in proSeen:
		pass
	else:
		proSeen.append(el)
		text="Pro\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(proPairList.count(el))+"\n"
		outputfile2.write(text)
		
#Ser	
serPairList=[]
for seq in range(len(serList)):
	for i in range(len(serList[seq])-1):
		pair=(serList[seq][i],serList[seq][i+1])
		serPairList.append(pair)
serPairList.sort()
serSeen=[]
for el in serPairList:
	if el in serSeen:
		pass
	else:
		serSeen.append(el)
		text="Ser\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(serPairList.count(el))+"\n"
		outputfile2.write(text)	
		
#Thr	
thrPairList=[]
for seq in range(len(thrList)):
	for i in range(len(thrList[seq])-1):
		pair=(thrList[seq][i],thrList[seq][i+1])
		thrPairList.append(pair)
thrPairList.sort()
thrSeen=[]
for el in thrPairList:
	if el in thrSeen:
		pass
	else:
		thrSeen.append(el)
		text="Thr\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(thrPairList.count(el))+"\n"
		outputfile2.write(text)
		
#Val	
valPairList=[]
for seq in range(len(valList)):
	for i in range(len(valList[seq])-1):
		pair=(valList[seq][i],valList[seq][i+1])
		valPairList.append(pair)
valPairList.sort()
valSeen=[]
for el in valPairList:
	if el in valSeen:
		pass
	else:
		valSeen.append(el)
		text="Val\t"+str(el[0])+"\t"+str(el[1])+"\t"+str(valPairList.count(el))+"\n"
		outputfile2.write(text)
		
		
# Printing results
print "Results computed only for the following amino acids: Ala, Arg, Gly, Ile, Leu, Ser, Pro, Thr, Val\n"
print "\tPrinted out output file 1 named as",sys.argv[1]+".total_codons, containing total counts of codons for each amino acid"
print "\tPrinted out output file 2 named as",sys.argv[1]+".paired_codons, containing total pairs of codons for each amino acid (all combinations)\n\n"
