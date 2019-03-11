#################################################
###  SCRIPT TO ANALYZE CODON AUTOCORRELATION  ###
#################################################
# Eva Maria Novoa, August 2013


# STEPS:
# 1. BUILD OBSERVED CODON CO-OCURRENCE MATRICES AND CODON COUNT VECTORS for each AA ###
# 2. METHOD 1: COMPUTES for each pair, STANDARD DEVIATIONS from EXPECTED.
# 3. METHOD 2: COMPUTES LOG LIKELIHOOD RATIOS & MAXIMUM LIKELIHOOD TEST. Outputs a p-value


# INFO:
## Input: merged files from "codon counts" and "codon pairs"
## Analysis only for codons from the AA:  Ala, Arg, Gly, Ile, Leu, Ser, Pro, Thr, Val 
## Thus, expects a total of 197 lines (all possible pair combinations within each amino acid) plus a header line
## Header line expected: AA	Codon1	Codon2	Paired_total	Total_codon1	Total_codon2

#####  1. READ DATA
data<-read.table("Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED",header=T) # file must be "***.merged"
head(data)
ala_data<-data[1:16,]
arg_data<-data[17:52,]
gly_data<-data[53:68,]
ile_data<-data[69:77,]
leu_data<-data[78:113,]
pro_data<-data[114:129,]
ser_data<-data[130:165,]
thr_data<-data[166:181,]
val_data<-data[182:197,]

##### 2. TOTAL SINGLE CODON READS for each codon and Amino acid (single sum)
ala_per_codon<-as.matrix(ala_data[1:4,6]) 
arg_per_codon<-as.matrix(arg_data[1:6,6])
gly_per_codon<-as.matrix(gly_data[1:4,6])
ile_per_codon<-as.matrix(ile_data[1:3,6])
leu_per_codon<-as.matrix(leu_data[1:6,6])
pro_per_codon<-as.matrix(pro_data[1:4,6])
ser_per_codon<-as.matrix(ser_data[1:6,6])
thr_per_codon<-as.matrix(thr_data[1:4,6])
val_per_codon<-as.matrix(val_data[1:4,6])

rownames(ala_per_codon)<-ala_data[1:4,3]
rownames(arg_per_codon)<-arg_data[1:6,3]
rownames(gly_per_codon)<-gly_data[1:4,3]
rownames(ile_per_codon)<-ile_data[1:3,3]
rownames(leu_per_codon)<-leu_data[1:6,3]
rownames(pro_per_codon)<-pro_data[1:4,3]
rownames(ser_per_codon)<-ser_data[1:6,3]
rownames(thr_per_codon)<-thr_data[1:4,3]
rownames(val_per_codon)<-val_data[1:4,3]

##### 3. TOTAL PAIRED CODON READS for each amino acid (build combination by combination and then bind to have a matrix)

#Codon1 (A-ended)
ala_pairs_codon1<-as.matrix(ala_data[1:4,4])
arg_pairs_codon1<-as.matrix(arg_data[1:6,4])
gly_pairs_codon1<-as.matrix(gly_data[1:4,4])
ile_pairs_codon1<-as.matrix(ile_data[1:3,4])
leu_pairs_codon1<-as.matrix(leu_data[1:6,4])
pro_pairs_codon1<-as.matrix(pro_data[1:4,4])
ser_pairs_codon1<-as.matrix(ser_data[1:6,4])
thr_pairs_codon1<-as.matrix(thr_data[1:4,4])
val_pairs_codon1<-as.matrix(val_data[1:4,4])
rownames(ala_pairs_codon1)<-ala_data[1:4,3]
rownames(arg_pairs_codon1)<-arg_data[1:6,3]
rownames(gly_pairs_codon1)<-gly_data[1:4,3]
rownames(ile_pairs_codon1)<-ile_data[1:3,3]
rownames(leu_pairs_codon1)<-leu_data[1:6,3]
rownames(pro_pairs_codon1)<-pro_data[1:4,3]
rownames(ser_pairs_codon1)<-ser_data[1:6,3]
rownames(thr_pairs_codon1)<-thr_data[1:4,3]
rownames(val_pairs_codon1)<-val_data[1:4,3]

#Codon2 (C-ended)
ala_pairs_codon2<-as.matrix(ala_data[5:8,4])
arg_pairs_codon2<-as.matrix(arg_data[7:12,4])
gly_pairs_codon2<-as.matrix(gly_data[5:8,4])
ile_pairs_codon2<-as.matrix(ile_data[4:6,4])
leu_pairs_codon2<-as.matrix(leu_data[7:12,4])
pro_pairs_codon2<-as.matrix(pro_data[5:8,4])
ser_pairs_codon2<-as.matrix(ser_data[7:12,4])
thr_pairs_codon2<-as.matrix(thr_data[5:8,4])
val_pairs_codon2<-as.matrix(val_data[5:8,4])
rownames(ala_pairs_codon2)<-ala_data[1:4,3]
rownames(arg_pairs_codon2)<-arg_data[1:6,3]
rownames(gly_pairs_codon2)<-gly_data[1:4,3]
rownames(ile_pairs_codon2)<-ile_data[1:3,3]
rownames(leu_pairs_codon2)<-leu_data[1:6,3]
rownames(pro_pairs_codon2)<-pro_data[1:4,3]
rownames(ser_pairs_codon2)<-ser_data[1:6,3]
rownames(thr_pairs_codon2)<-thr_data[1:4,3]
rownames(val_pairs_codon2)<-val_data[1:4,3]

#Codon3 (G-ended)
ala_pairs_codon3<-as.matrix(ala_data[9:12,4])
arg_pairs_codon3<-as.matrix(arg_data[13:18,4])
gly_pairs_codon3<-as.matrix(gly_data[9:12,4])
#ile_pairs_codon3<-as.matrix(ile_data[7:9,4]) --> doesn't exist G-ended in Ile
leu_pairs_codon3<-as.matrix(leu_data[13:18,4])
pro_pairs_codon3<-as.matrix(pro_data[9:12,4])
ser_pairs_codon3<-as.matrix(ser_data[13:18,4])
thr_pairs_codon3<-as.matrix(thr_data[9:12,4])
val_pairs_codon3<-as.matrix(val_data[9:12,4])
rownames(ala_pairs_codon3)<-ala_data[1:4,3]
rownames(arg_pairs_codon3)<-arg_data[1:6,3]
rownames(gly_pairs_codon3)<-gly_data[1:4,3]
#rownames(ile_pairs_codon3)<-ile_data[1:3,3]
rownames(leu_pairs_codon3)<-leu_data[1:6,3]
rownames(pro_pairs_codon3)<-pro_data[1:4,3]
rownames(ser_pairs_codon3)<-ser_data[1:6,3]
rownames(thr_pairs_codon3)<-thr_data[1:4,3]
rownames(val_pairs_codon3)<-val_data[1:4,3]

#Codon4 (T-ended)
ala_pairs_codon4<-as.matrix(ala_data[13:16,4])
arg_pairs_codon4<-as.matrix(arg_data[19:24,4])
gly_pairs_codon4<-as.matrix(gly_data[13:16,4])
ile_pairs_codon4<-as.matrix(ile_data[7:9,4])
leu_pairs_codon4<-as.matrix(leu_data[19:24,4])
pro_pairs_codon4<-as.matrix(pro_data[13:16,4])
ser_pairs_codon4<-as.matrix(ser_data[19:24,4])
thr_pairs_codon4<-as.matrix(thr_data[13:16,4])
val_pairs_codon4<-as.matrix(val_data[13:16,4])
rownames(ala_pairs_codon4)<-ala_data[1:4,3]
rownames(arg_pairs_codon4)<-arg_data[1:6,3]
rownames(gly_pairs_codon4)<-gly_data[1:4,3]
rownames(ile_pairs_codon4)<-ile_data[1:3,3]
rownames(leu_pairs_codon4)<-leu_data[1:6,3]
rownames(pro_pairs_codon4)<-pro_data[1:4,3]
rownames(ser_pairs_codon4)<-ser_data[1:6,3]
rownames(thr_pairs_codon4)<-thr_data[1:4,3]
rownames(val_pairs_codon4)<-val_data[1:4,3]

#Codon5 (only for 6box)
arg_pairs_codon5<-as.matrix(arg_data[25:30,4])
leu_pairs_codon5<-as.matrix(leu_data[25:30,4])
ser_pairs_codon5<-as.matrix(ser_data[25:30,4])
rownames(arg_pairs_codon5)<-arg_data[1:6,3]
rownames(leu_pairs_codon5)<-leu_data[1:6,3]
rownames(ser_pairs_codon5)<-ser_data[1:6,3]

#Codon6 (only for 6box)
arg_pairs_codon6<-as.matrix(arg_data[31:36,4])
leu_pairs_codon6<-as.matrix(leu_data[31:36,4])
ser_pairs_codon6<-as.matrix(ser_data[31:36,4])
rownames(arg_pairs_codon6)<-arg_data[1:6,3]
rownames(leu_pairs_codon6)<-leu_data[1:6,3]
rownames(ser_pairs_codon6)<-ser_data[1:6,3]

# Merge all into matrix of codon pairs
ala_pairs_all<-cbind(ala_pairs_codon1,ala_pairs_codon2,ala_pairs_codon3,ala_pairs_codon4)
arg_pairs_all<-cbind(arg_pairs_codon1,arg_pairs_codon2,arg_pairs_codon3,arg_pairs_codon4,arg_pairs_codon5,arg_pairs_codon6)
gly_pairs_all<-cbind(gly_pairs_codon1,gly_pairs_codon2,gly_pairs_codon3,gly_pairs_codon4)
ile_pairs_all<-cbind(ile_pairs_codon1,ile_pairs_codon2,ile_pairs_codon4)
leu_pairs_all<-cbind(leu_pairs_codon1,leu_pairs_codon2,leu_pairs_codon3,leu_pairs_codon4,leu_pairs_codon5,leu_pairs_codon6)
pro_pairs_all<-cbind(pro_pairs_codon1,pro_pairs_codon2,pro_pairs_codon3,pro_pairs_codon4)
ser_pairs_all<-cbind(ser_pairs_codon1,ser_pairs_codon2,ser_pairs_codon3,ser_pairs_codon4,ser_pairs_codon5,ser_pairs_codon6)
thr_pairs_all<-cbind(thr_pairs_codon1,thr_pairs_codon2,thr_pairs_codon3,thr_pairs_codon4)
val_pairs_all<-cbind(val_pairs_codon1,val_pairs_codon2,val_pairs_codon3,val_pairs_codon4)
colnames(ala_pairs_all)<-rownames(ala_pairs_codon4)
colnames(arg_pairs_all)<-rownames(arg_pairs_codon4)
colnames(gly_pairs_all)<-rownames(gly_pairs_codon4)
colnames(ile_pairs_all)<-rownames(ile_pairs_codon4)
colnames(leu_pairs_all)<-rownames(leu_pairs_codon4)
colnames(pro_pairs_all)<-rownames(pro_pairs_codon4)
colnames(ser_pairs_all)<-rownames(ser_pairs_codon4)
colnames(thr_pairs_all)<-rownames(thr_pairs_codon4)
colnames(val_pairs_all)<-rownames(val_pairs_codon4)


##### 4. GET TOTALS PER COLUMN, PER ROW AND PER AMINO ACID

ala_row_sums<-rowSums(ala_pairs_all)
ala_col_sums<-colSums(ala_pairs_all)
ala_total_sum<-sum(colSums(ala_pairs_all))

arg_row_sums<-rowSums(arg_pairs_all)
arg_col_sums<-colSums(arg_pairs_all)
arg_total_sum<-sum(colSums(arg_pairs_all))

gly_row_sums<-rowSums(gly_pairs_all)
gly_col_sums<-colSums(gly_pairs_all)
gly_total_sum<-sum(colSums(gly_pairs_all))

ile_row_sums<-rowSums(ile_pairs_all)
ile_col_sums<-colSums(ile_pairs_all)
ile_total_sum<-sum(colSums(ile_pairs_all))

leu_row_sums<-rowSums(leu_pairs_all)
leu_col_sums<-colSums(leu_pairs_all)
leu_total_sum<-sum(colSums(leu_pairs_all))

pro_row_sums<-rowSums(pro_pairs_all)
pro_col_sums<-colSums(pro_pairs_all)
pro_total_sum<-sum(colSums(pro_pairs_all))

ser_row_sums<-rowSums(ser_pairs_all)
ser_col_sums<-colSums(ser_pairs_all)
ser_total_sum<-sum(colSums(ser_pairs_all))

thr_row_sums<-rowSums(thr_pairs_all)
thr_col_sums<-colSums(thr_pairs_all)
thr_total_sum<-sum(colSums(thr_pairs_all))

val_row_sums<-rowSums(val_pairs_all)
val_col_sums<-colSums(val_pairs_all)
val_total_sum<-sum(colSums(val_pairs_all))


#For a given amino acid, there will be a difference between the total number of pairs and the total number of codons that corresponds to the number of sequences analyzed (because for example, if the protein has 7 Ala codons, that will give 6 Ala codon pairs).
ala_per_codon
colSums(ala_per_codon)
ala_pairs_all
colSums(ala_pairs_all)

ala_num_of_first_codons<-colSums(ala_per_codon)-ala_total_sum
arg_num_of_first_codons<-colSums(arg_per_codon)-arg_total_sum
gly_num_of_first_codons<-colSums(gly_per_codon)-gly_total_sum
ile_num_of_first_codons<-colSums(ile_per_codon)-ile_total_sum
leu_num_of_first_codons<-colSums(leu_per_codon)-leu_total_sum
ser_num_of_first_codons<-colSums(ser_per_codon)-ser_total_sum
pro_num_of_first_codons<-colSums(pro_per_codon)-pro_total_sum
thr_num_of_first_codons<-colSums(thr_per_codon)-thr_total_sum
val_num_of_first_codons<-colSums(val_per_codon)-val_total_sum

### PROBABILITIES

# Codon counts probability vectors
prob_ala_codon<-ala_per_codon/sum(ala_per_codon)
prob_arg_codon<-arg_per_codon/sum(arg_per_codon)
prob_gly_codon<-gly_per_codon/sum(gly_per_codon)
prob_ile_codon<-ile_per_codon/sum(ile_per_codon)
prob_leu_codon<-leu_per_codon/sum(leu_per_codon)
prob_pro_codon<-pro_per_codon/sum(pro_per_codon)
prob_ser_codon<-ser_per_codon/sum(ser_per_codon)
prob_thr_codon<-thr_per_codon/sum(thr_per_codon)
prob_val_codon<-val_per_codon/sum(val_per_codon)

# Codon pairs probability matrices (normalized per row)
prob_ala_pairs<-ala_pairs_all/ala_row_sums 
prob_arg_pairs<-arg_pairs_all/arg_row_sums
prob_gly_pairs<-gly_pairs_all/gly_row_sums
prob_ile_pairs<-ile_pairs_all/ile_row_sums
prob_leu_pairs<-leu_pairs_all/leu_row_sums
prob_pro_pairs<-pro_pairs_all/pro_row_sums
prob_ser_pairs<-ser_pairs_all/ser_row_sums
prob_thr_pairs<-thr_pairs_all/thr_row_sums
prob_val_pairs<-val_pairs_all/val_row_sums



##### 5. ANALYSIS METHOD 1: Z-SCORES AND STANDARD DEVIATIONS FROM EXPECTED

## 5.1 OBTAIN EXPECTED CO-OCCURENCE MATRICES

ala_exp<-matrix(0,ncol=length(ala_row_sums),nrow=length(ala_row_sums)) 
for (i in 1:length(ala_row_sums)) {
	for (j in 1:length(ala_row_sums)) {
		ala_exp[i,j]<-prob_ala_codon[i]*prob_ala_codon[j]* colSums(ala_per_codon)
	}
}
colnames(ala_exp)<-rownames(ala_pairs_all)
rownames(ala_exp)<-rownames(ala_pairs_all)

arg_exp<-matrix(0,ncol=length(arg_row_sums),nrow=length(arg_row_sums)) 
for (i in 1:length(arg_row_sums)) {
	for (j in 1:length(arg_row_sums)) {
		arg_exp[i,j]<-prob_arg_codon[i]*prob_arg_codon[j]* colSums(arg_per_codon)
	}
}
colnames(arg_exp)<-rownames(arg_pairs_all)
rownames(arg_exp)<-rownames(arg_pairs_all)

gly_exp<-matrix(0,ncol=length(gly_row_sums),nrow=length(gly_row_sums)) 
for (i in 1:length(gly_row_sums)) {
	for (j in 1:length(gly_row_sums)) {
		gly_exp[i,j]<-prob_gly_codon[i]*prob_gly_codon[j]* colSums(gly_per_codon)
	}
}
colnames(gly_exp)<-rownames(gly_pairs_all)
rownames(gly_exp)<-rownames(gly_pairs_all)

ile_exp<-matrix(0,ncol=length(ile_row_sums),nrow=length(ile_row_sums)) 
for (i in 1:length(ile_row_sums)) {
	for (j in 1:length(ile_row_sums)) {
		ile_exp[i,j]<-prob_ile_codon[i]*prob_ile_codon[j]* colSums(ile_per_codon)
	}
}
colnames(ile_exp)<-rownames(ile_pairs_all)
rownames(ile_exp)<-rownames(ile_pairs_all)

leu_exp<-matrix(0,ncol=length(leu_row_sums),nrow=length(leu_row_sums)) 
for (i in 1:length(leu_row_sums)) {
	for (j in 1:length(leu_row_sums)) {
		leu_exp[i,j]<-prob_leu_codon[i]*prob_leu_codon[j]* colSums(leu_per_codon)
	}
}
colnames(leu_exp)<-rownames(leu_pairs_all)
rownames(leu_exp)<-rownames(leu_pairs_all)

pro_exp<-matrix(0,ncol=length(pro_row_sums),nrow=length(pro_row_sums)) 
for (i in 1:length(pro_row_sums)) {
	for (j in 1:length(pro_row_sums)) {
		pro_exp[i,j]<-prob_pro_codon[i]*prob_pro_codon[j]* colSums(pro_per_codon)
	}
}
colnames(pro_exp)<-rownames(pro_pairs_all)
rownames(pro_exp)<-rownames(pro_pairs_all)

ser_exp<-matrix(0,ncol=length(ser_row_sums),nrow=length(ser_row_sums)) 
for (i in 1:length(ser_row_sums)) {
	for (j in 1:length(ser_row_sums)) {
		ser_exp[i,j]<-prob_ser_codon[i]*prob_ser_codon[j]* colSums(ser_per_codon)
	}
}
colnames(ser_exp)<-rownames(ser_pairs_all)
rownames(ser_exp)<-rownames(ser_pairs_all)

thr_exp<-matrix(0,ncol=length(thr_row_sums),nrow=length(thr_row_sums)) 
for (i in 1:length(thr_row_sums)) {
	for (j in 1:length(thr_row_sums)) {
		thr_exp[i,j]<-prob_thr_codon[i]*prob_thr_codon[j]* colSums(thr_per_codon)
	}
}
colnames(thr_exp)<-rownames(thr_pairs_all)
rownames(thr_exp)<-rownames(thr_pairs_all)

val_exp<-matrix(0,ncol=length(val_row_sums),nrow=length(val_row_sums)) 
for (i in 1:length(val_row_sums)) {
	for (j in 1:length(val_row_sums)) {
		val_exp[i,j]<-prob_val_codon[i]*prob_val_codon[j]* colSums(val_per_codon)
	}
}
colnames(val_exp)<-rownames(val_pairs_all)
rownames(val_exp)<-rownames(val_pairs_all)


## 5.2 COMPUTE Z-SCORES (=observed-expected/SE)

#Observed: e.g. ala_pairs_all
#Expected: e.g. ala_exp
#Standard dev from expected = (obs-exp)/sqrt(exp*(1-exp/total)

ala_sdevs<-matrix(0,ncol=length(ala_row_sums),nrow=length(ala_row_sums)) 
for (i in 1:length(ala_row_sums)) { # i = row pos
	for (j in 1:length(ala_row_sums)) { # j = col pos
		ala_sdevs[i,j]<-(ala_pairs_all[i,j]-ala_exp[i,j])/(sqrt (ala_exp[i,j]*(1-ala_exp[i,j]/ala_total_sum)))
	}
}
colnames(ala_sdevs)<-rownames(ala_pairs_all)
rownames(ala_sdevs)<-rownames(ala_pairs_all)


arg_sdevs<-matrix(0,ncol=length(arg_row_sums),nrow=length(arg_row_sums)) 
for (i in 1:length(arg_row_sums)) { # i = row pos
	for (j in 1:length(arg_row_sums)) { # j = col pos
		arg_sdevs[i,j]<-(arg_pairs_all[i,j]-arg_exp[i,j])/(sqrt (arg_exp[i,j]*(1-arg_exp[i,j]/arg_total_sum)))
	}
}
colnames(arg_sdevs)<-rownames(arg_pairs_all)
rownames(arg_sdevs)<-rownames(arg_pairs_all)


gly_sdevs<-matrix(0,ncol=length(gly_row_sums),nrow=length(gly_row_sums)) 
for (i in 1:length(gly_row_sums)) { # i = row pos
	for (j in 1:length(gly_row_sums)) { # j = col pos
		gly_sdevs[i,j]<-(gly_pairs_all[i,j]-gly_exp[i,j])/(sqrt (gly_exp[i,j]*(1-gly_exp[i,j]/gly_total_sum)))
	}
}
colnames(gly_sdevs)<-rownames(gly_pairs_all)
rownames(gly_sdevs)<-rownames(gly_pairs_all)


ile_sdevs<-matrix(0,ncol=length(ile_row_sums),nrow=length(ile_row_sums)) 
for (i in 1:length(ile_row_sums)) { # i = row pos
	for (j in 1:length(ile_row_sums)) { # j = col pos
		ile_sdevs[i,j]<-(ile_pairs_all[i,j]-ile_exp[i,j])/(sqrt (ile_exp[i,j]*(1-ile_exp[i,j]/ile_total_sum)))
	}
}
colnames(ile_sdevs)<-rownames(ile_pairs_all)
rownames(ile_sdevs)<-rownames(ile_pairs_all)


leu_sdevs<-matrix(0,ncol=length(leu_row_sums),nrow=length(leu_row_sums)) 
for (i in 1:length(leu_row_sums)) { # i = row pos
	for (j in 1:length(leu_row_sums)) { # j = col pos
		leu_sdevs[i,j]<-(leu_pairs_all[i,j]-leu_exp[i,j])/(sqrt (leu_exp[i,j]*(1-leu_exp[i,j]/leu_total_sum)))
	}
}
colnames(leu_sdevs)<-rownames(leu_pairs_all)
rownames(leu_sdevs)<-rownames(leu_pairs_all)


pro_sdevs<-matrix(0,ncol=length(pro_row_sums),nrow=length(pro_row_sums)) 
for (i in 1:length(pro_row_sums)) { # i = row pos
	for (j in 1:length(pro_row_sums)) { # j = col pos
		pro_sdevs[i,j]<-(pro_pairs_all[i,j]-pro_exp[i,j])/(sqrt (pro_exp[i,j]*(1-pro_exp[i,j]/pro_total_sum)))
	}
}
colnames(pro_sdevs)<-rownames(pro_pairs_all)
rownames(pro_sdevs)<-rownames(pro_pairs_all)


ser_sdevs<-matrix(0,ncol=length(ser_row_sums),nrow=length(ser_row_sums)) 
for (i in 1:length(ser_row_sums)) { # i = row pos
	for (j in 1:length(ser_row_sums)) { # j = col pos
		ser_sdevs[i,j]<-(ser_pairs_all[i,j]-ser_exp[i,j])/(sqrt (ser_exp[i,j]*(1-ser_exp[i,j]/ser_total_sum)))
	}
}
colnames(ser_sdevs)<-rownames(ser_pairs_all)
rownames(ser_sdevs)<-rownames(ser_pairs_all)


thr_sdevs<-matrix(0,ncol=length(thr_row_sums),nrow=length(thr_row_sums)) 
for (i in 1:length(thr_row_sums)) { # i = row pos
	for (j in 1:length(thr_row_sums)) { # j = col pos
		thr_sdevs[i,j]<-(thr_pairs_all[i,j]-thr_exp[i,j])/(sqrt (thr_exp[i,j]*(1-thr_exp[i,j]/thr_total_sum)))
	}
}
colnames(thr_sdevs)<-rownames(thr_pairs_all)
rownames(thr_sdevs)<-rownames(thr_pairs_all)


val_sdevs<-matrix(0,ncol=length(val_row_sums),nrow=length(val_row_sums)) 
for (i in 1:length(val_row_sums)) { # i = row pos
	for (j in 1:length(val_row_sums)) { # j = col pos
		val_sdevs[i,j]<-(val_pairs_all[i,j]-val_exp[i,j])/(sqrt (val_exp[i,j]*(1-val_exp[i,j]/val_total_sum)))
	}
}
colnames(val_sdevs)<-rownames(val_pairs_all)
rownames(val_sdevs)<-rownames(val_pairs_all)

sink("Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED_sdevs.txt",append=FALSE, split=FALSE)
ala_sdevs
arg_sdevs
gly_sdevs
ile_sdevs
leu_sdevs
ser_sdevs
pro_sdevs
thr_sdevs
val_sdevs
sink()


##### 6. ANALYSIS THROUGH LIKELIHOOD RATIO MODELS

## 6.1.  OBTAIN PROBABILITY MATRICES

# Codon pairs probability matrices (normalized per matrix) ; ala_total_sum --> total sum of pairs --> "OBSERVED" MATRIX NORMALIZED TO 1
prob_matrix_ala_pairs<-ala_pairs_all/ala_total_sum 
prob_matrix_arg_pairs<-arg_pairs_all/arg_total_sum
prob_matrix_gly_pairs<-gly_pairs_all/gly_total_sum
prob_matrix_ile_pairs<-ile_pairs_all/ile_total_sum
prob_matrix_leu_pairs<-leu_pairs_all/leu_total_sum
prob_matrix_pro_pairs<-pro_pairs_all/pro_total_sum
prob_matrix_ser_pairs<-ser_pairs_all/ser_total_sum
prob_matrix_thr_pairs<-thr_pairs_all/thr_total_sum
prob_matrix_val_pairs<-val_pairs_all/val_total_sum


# Codon pairs probability matrices (normalized per matrix) --> "EXPECTED" MATRIX NORMALIZED TO 1
prob_matrix_ala_exp<-ala_exp/colSums(ala_per_codon)
prob_matrix_arg_exp<-arg_exp/colSums(arg_per_codon)
prob_matrix_gly_exp<-gly_exp/colSums(gly_per_codon)
prob_matrix_ile_exp<-ile_exp/colSums(ile_per_codon)
prob_matrix_leu_exp<-leu_exp/colSums(leu_per_codon)
prob_matrix_pro_exp<-pro_exp/colSums(pro_per_codon)
prob_matrix_ser_exp<-ser_exp/colSums(ser_per_codon)
prob_matrix_thr_exp<-thr_exp/colSums(thr_per_codon)
prob_matrix_val_exp<-val_exp/colSums(val_per_codon)

## 6.2. LOG LIKELIHOOD RATIOS AND THE P-VALUE of THE LOG LIKELIHOOD STATISTIC

# Step 1 :
# Compute the log likelihood ratios
# log P(null)/P(alt) = log P(null) - log P(alt)
# P(null)= P(codon1)^num(codon1)*P(codon2)^num(codon2)....
# log P(null) = num(codon1)*log P(codon1) + num(codon2)*log P(codon2)....
# log P(alt)= log P(codon pairs) + log P(first codon)
# log P(alt)= num(

# Step 2:
# Compute the lrt between two log likelihoods (Gerald Quon)

lrt <- function(lnull,lalt,df=1) {
	stat=-2*(lnull-lalt); 
	pval=1-pchisq(q=stat,df,ncp=0); 
	pval;
};

# Remember to change df for each case!! 
# 3box --> df(null)=2;df(alt)=6+2 -->df2-1--> 6
# 4box --> df(null)=3;df(alt)=12+3 -->df2-1--> 12
# 6box --> df(null)=5;df(alt)=30+5 -->df2-1--> 30

#Ala
log_Pnull=0  # individual codons probabilities
for (i in 1:length(ala_row_sums)) {
	log_Pnull<-log_Pnull+ala_per_codon[i]*log(prob_ala_codon[i])
}

log_Palt_part1=0  # first codon probabilities --> assumed from general codon abundances
for (i in 1:length(ala_row_sums)) {
	log_Palt_part1<-log_Palt_part1+prob_ala_codon[i]*ala_num_of_first_codons*log(prob_ala_codon[i])
}

log_Palt_part2=0  # pair probability
for (i in 1:length(ala_row_sums)) {
	for (j in 1:length(ala_row_sums)) {
		log_Palt_part2<-log_Palt_part2+ala_pairs_all[i,j]*log(prob_ala_pairs[i,j]) # prob_ala_pairs OR prob_matrix_ala_pairs?
	}
}
log_Palt<-log_Palt_part1+log_Palt_part2
pval_ala<-lrt(log_Pnull,log_Palt,12)

colSums(ala_per_codon) # num of ind. ala codons
ala_total_sum  # num of ala pairs
ala_num_of_first_codons
log_Pnull
log_Palt
log_Palt_part1
log_Palt_part2
pval_ala
prob_ala_codon
ala_num_of_first_codons
ala_per_codon
ala_exp
prob_ala_pairs
ala_pairs_all

# Ala codon pair by codon pair

#a) Using number of codons multiplying the probabilities
ala_pval_matrix<-matrix(0,ncol=length(ala_row_sums),nrow=length(ala_row_sums)) 
colnames(ala_pval_matrix)<-rownames(ala_pairs_all)
rownames(ala_pval_matrix)<-rownames(ala_pairs_all)

for (i in 1:length(ala_row_sums)) {
	for (j in 1:length(ala_row_sums)) {	
		log_Pnull<-ala_per_codon[i]*log(prob_ala_codon[i])+ala_per_codon[j]*log(prob_ala_codon[j])
		log_Palt<-ala_pairs_all[i,j]*log(prob_matrix_ala_pairs[i,j]) 
		ala_pval_matrix[i,j]<-lrt(log_Pnull,log_Palt,0) # df? 
	}
}

log_Pnull
log_Palt
ala_per_codon
prob_ala_codon
ala_pairs_all
prob_matrix_ala_pairs
ala_pval_matrix

#b) Only using the probabilities
for (i in 1:length(ala_row_sums)) {
	for (j in 1:length(ala_row_sums)) {	
		log_Pnull<-log(prob_ala_codon[i])+log(prob_ala_codon[j])
		log_Palt<-log(prob_matrix_ala_pairs[i,j]) 
		ala_pval_matrix[i,j]<-lrt(log_Pnull,log_Palt,0) # df? 
	}
}

log_Pnull
log_Palt
ala_per_codon
prob_ala_codon
ala_pairs_all
prob_matrix_ala_pairs
ala_pval_matrix

#Val
#b) Only using the probabilities
val_pval_matrix<-matrix(0,ncol=length(val_row_sums),nrow=length(val_row_sums)) 
colnames(val_pval_matrix)<-rownames(val_pairs_all)
rownames(val_pval_matrix)<-rownames(val_pairs_all)

for (i in 1:length(val_row_sums)) {
	for (j in 1:length(val_row_sums)) {	
		log_Pnull<-log(prob_val_codon[i])+log(prob_val_codon[j])
		log_Palt<-log(prob_matrix_val_pairs[i,j]) 
		val_pval_matrix[i,j]<-lrt(log_Pnull,log_Palt,1) # df? 
	}
}

log_Pnull
log_Palt
prob_val_codon
prob_matrix_val_pairs
val_pval_matrix


## NOT WORKING YET

##### 7. TO ADJUST FOR EFFECT SIZE (COMPARISON ACROSS GENOMES) BUILD MATRICES OF (P(obs)-P(exp))/P(exp)  -- Irwin recomendation

#prob_matrix_ala_pairs --> observed matrix (norm to sum 1)
#prob_matrix_ala_exp --> expected matrix (norm to sum 1)

ala_effect_matrix<-(prob_matrix_ala_pairs-prob_matrix_ala_exp)/prob_matrix_ala_exp
arg_effect_matrix<-(prob_matrix_arg_pairs-prob_matrix_arg_exp)/prob_matrix_arg_exp
gly_effect_matrix<-(prob_matrix_gly_pairs-prob_matrix_gly_exp)/prob_matrix_gly_exp
ile_effect_matrix<-(prob_matrix_ile_pairs-prob_matrix_ile_exp)/prob_matrix_ile_exp
leu_effect_matrix<-(prob_matrix_leu_pairs-prob_matrix_leu_exp)/prob_matrix_leu_exp
pro_effect_matrix<-(prob_matrix_pro_pairs-prob_matrix_pro_exp)/prob_matrix_pro_exp
ser_effect_matrix<-(prob_matrix_ser_pairs-prob_matrix_ser_exp)/prob_matrix_ser_exp
thr_effect_matrix<-(prob_matrix_thr_pairs-prob_matrix_thr_exp)/prob_matrix_thr_exp
val_effect_matrix<-(prob_matrix_val_pairs-prob_matrix_val_exp)/prob_matrix_val_exp

sink("Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED_effect_matrix.txt",append=FALSE, split=FALSE)
ala_effect_matrix
arg_effect_matrix
gly_effect_matrix
ile_effect_matrix
leu_effect_matrix
ser_effect_matrix
pro_effect_matrix
thr_effect_matrix
val_effect_matrix
sink()

# Convert into appended vectors, to have one vector per species
ala_effect_vector<-as.vector(t(ala_effect_matrix))
arg_effect_vector<-as.vector(t(arg_effect_matrix))
gly_effect_vector<-as.vector(t(gly_effect_matrix))
ile_effect_vector<-as.vector(t(ile_effect_matrix))
leu_effect_vector<-as.vector(t(leu_effect_matrix))
ser_effect_vector<-as.vector(t(ser_effect_matrix))
pro_effect_vector<-as.vector(t(pro_effect_matrix))
thr_effect_vector<-as.vector(t(thr_effect_matrix))
val_effect_vector<-as.vector(t(val_effect_matrix))

Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED<-c(ala_effect_vector, arg_effect_vector, gly_effect_vector, ile_effect_vector, leu_effect_vector, ser_effect_vector, pro_effect_vector, thr_effect_vector, val_effect_vector)  # all_all_vector

write.table(data.frame(Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED), "Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED_effect_vector.txt", row.names=F)

##### 8. BUILD NEW INDEX: Relative synonymous codon pair usage (RSCPU)  (P(obs)-/P(exp)  

## Way A) (P(obs)-P(exp))/P(exp)
ala_rscpu_matrix<-prob_matrix_ala_pairs/prob_matrix_ala_exp
arg_rscpu_matrix<-prob_matrix_arg_pairs/prob_matrix_arg_exp
gly_rscpu_matrix<-prob_matrix_gly_pairs/prob_matrix_gly_exp
ile_rscpu_matrix<-prob_matrix_ile_pairs/prob_matrix_ile_exp
leu_rscpu_matrix<-prob_matrix_leu_pairs/prob_matrix_leu_exp
pro_rscpu_matrix<-prob_matrix_pro_pairs/prob_matrix_pro_exp
ser_rscpu_matrix<-prob_matrix_ser_pairs/prob_matrix_ser_exp
thr_rscpu_matrix<-prob_matrix_thr_pairs/prob_matrix_thr_exp
val_rscpu_matrix<-prob_matrix_val_pairs/prob_matrix_val_exp


# Convert into appended vectors, to have one vector per species
ala_rscpu_vector<-as.vector(t(ala_rscpu_matrix))
arg_rscpu_vector<-as.vector(t(arg_rscpu_matrix))
gly_rscpu_vector<-as.vector(t(gly_rscpu_matrix))
ile_rscpu_vector<-as.vector(t(ile_rscpu_matrix))
leu_rscpu_vector<-as.vector(t(leu_rscpu_matrix))
ser_rscpu_vector<-as.vector(t(ser_rscpu_matrix))
pro_rscpu_vector<-as.vector(t(pro_rscpu_matrix))
thr_rscpu_vector<-as.vector(t(thr_rscpu_matrix))
val_rscpu_vector<-as.vector(t(val_rscpu_matrix))

Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED<-c(ala_rscpu_vector, arg_rscpu_vector, gly_rscpu_vector, ile_rscpu_vector, leu_rscpu_vector, ser_rscpu_vector, pro_rscpu_vector, thr_rscpu_vector, val_rscpu_vector)  # all_all_vector

write.table(data.frame(Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED), "Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED_rscpu_vector.txt", row.names=F)


## Way B) (P(obs)-P(exp))/P(exp)


##Relative synonymous codon usage
#RSCU(codon)= Obs frequency of codon(compared to total for that AA)/ Exp frequency of codon(compared to total for that AA).
#Ex. Ala(AGC)= obs frequency/0.25

##Relative synonumous codon pair usage
#RSCU(codon pair) = Obs frequency of codon(compared to total for that AA)/ Exp frequency of codon(compared to total for that AA
#Ex. Ala(AGC-AGC)= obs frequency/0.0625


ala_rscpu2_matrix<-prob_matrix_ala_pairs/(1/16)
arg_rscpu2_matrix<-prob_matrix_arg_pairs/(1/36)
gly_rscpu2_matrix<-prob_matrix_gly_pairs/(1/16)
ile_rscpu2_matrix<-prob_matrix_ile_pairs/(1/9)
leu_rscpu2_matrix<-prob_matrix_leu_pairs/(1/36)
pro_rscpu2_matrix<-prob_matrix_pro_pairs/(1/16)
ser_rscpu2_matrix<-prob_matrix_ser_pairs/(1/36)
thr_rscpu2_matrix<-prob_matrix_thr_pairs/(1/16)
val_rscpu2_matrix<-prob_matrix_val_pairs/(1/16)

# Convert into appended vectors, to have one vector per species
ala_rscpu2_vector<-as.vector(t(ala_rscpu2_matrix))
arg_rscpu2_vector<-as.vector(t(arg_rscpu2_matrix))
gly_rscpu2_vector<-as.vector(t(gly_rscpu2_matrix))
ile_rscpu2_vector<-as.vector(t(ile_rscpu2_matrix))
leu_rscpu2_vector<-as.vector(t(leu_rscpu2_matrix))
ser_rscpu2_vector<-as.vector(t(ser_rscpu2_matrix))
pro_rscpu2_vector<-as.vector(t(pro_rscpu2_matrix))
thr_rscpu2_vector<-as.vector(t(thr_rscpu2_matrix))
val_rscpu2_vector<-as.vector(t(val_rscpu2_matrix))

Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED<-c(ala_rscpu2_vector, arg_rscpu2_vector, gly_rscpu2_vector, ile_rscpu2_vector, leu_rscpu2_vector, ser_rscpu2_vector, pro_rscpu2_vector, thr_rscpu2_vector, val_rscpu2_vector)  # all_all_vector

write.table(data.frame(Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED), "Saccharomyces_cerevisiae.CDS.fa.codon_counts.merged.FIXED_rscpu_vector_way2.txt", row.names=F)
