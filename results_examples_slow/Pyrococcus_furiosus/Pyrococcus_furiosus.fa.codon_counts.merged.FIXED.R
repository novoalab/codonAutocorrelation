### FROM FIXED merged AUTOCORRELATION COUNTS, get MATRIX ####

	
# Genetic code	
ala=c('GCA','GCC','GCG','GCT')
arg=c('CGA','CGC','CGG','CGT','AGA','AGG')
asn=c('AAT','AAC')
asp=c('GAT','GAC')
cys=c('TGT','TGC')
gln=c('CAG','CAA')
glu=c('GAG','GAA')
gly=c('GGT','GGG','GGC','GGA')
his=c('CAT','CAC')
ile=c('ATT','ATC','ATA')
leu=c('CTT','TTG','CTG','CTC','TTA','CTA')
lys=c('AAG','AAA')
met=c('ATG')
phe=c('TTT','TTC')
pro=c('CCT','CCG','CCC','CCA')
ser=c('AGT','TCT','TCG','AGC','TCC','TCA')
ter=c('TAG','TGA','TAA')
thr=c('ACT','ACG','ACC','ACA')
trp=c('TGG')
tyr=c('TAT','TAC')
val=c('GTT','GTG','GTC','GTA')


all_aa=c('GCA','GCC','GCG','GCT','CGA','CGC','CGG','CGT','AGA','AGG','AAT','AAC','GAT','GAC','TGT','TGC','CAG','CAA','GAG','GAA','GGT','GGG','GGC','GGA','CAT','CAC','ATT','ATC','ATA','CTT','TTG','CTG','CTC','TTA','CTA','AAG','AAA','ATG','TTT','TTC','CCT','CCG','CCC','CCA','AGT','TCT','TCG','AGC','TCC','TCA','TAG','TGA','TAA','ACT','ACG','ACC','ACA','TGG','TAT','TAC','GTT','GTG','GTC','GTA')
all_aa.sorted<-sort(all_aa)

## FUNCTIONS ##

from_codonFIXED_to_matrix_FINAL<-function(FILE) {
	# Read data
	data<-read.table(FILE)
	matrix<-matrix(data[,3],ncol=64,nrow=64, byrow = TRUE)  # so 1st codon in rows; 2nd codon in columns !!!!!
	colnames(matrix)<-all_aa.sorted
	row.names(matrix)<-all_aa.sorted
	# Reoder matrix to be by AA
	matrix.ordered_rows<-matrix[all_aa,,drop=FALSE]
	tmatrix<-t(matrix.ordered_rows)
	tmatrix.ordered_rows<-tmatrix[all_aa,,drop=FALSE]
	matrix.CLEAN<-t(tmatrix.ordered_rows)
	# Remove stop codon rows and columns
	matrix.FINAL<- matrix.CLEAN[-51:-53,-51:-53]
	return(matrix.FINAL)
}

plot_matrix_levelplot<-function(mat,my_title) {
	# Plot
	library(latticeExtra)
	x<-t(mat)
	rgb.palette <- colorRampPalette(c("white","yellow","orange","red","darkred" ), space = "rgb")
	levelplot(x,  
	  	main=my_title,
      	aspect = "fill",
      	#scales = list(x = list(rot = 90)),
      	scales=list(x=list(cex=0.5, rot = 90),y=list(cex=0.5)),
      	colorkey = list(space = "right"),
      	col.regions=rgb.palette(120), 
     	#cuts=100, 
      	#at=seq(0,0.01,0.0001),
      	xlab="",
     	ylab="")
} 

plot_matrix_sdev<-function(mat,my_title) {
	# Plot
	library(latticeExtra)
	mat<-replace(mat, is.nan(mat),0)
	x<-t(mat)
	my_value<-round(max(abs(mat))+0.5)
	my_step=my_value*2/100
	rgb.palette <- colorRampPalette(c("red","black","green"), space = "rgb")
	levelplot(x,  
	  	main=my_title,
      	aspect = "fill",
      	#scales = list(x = list(rot = 90)),
      	scales=list(x=list(cex=0.5, rot = 90),y=list(cex=0.5)),
      	colorkey = list(space = "right"),
      	col.regions=rgb.palette(120), 
     	cuts=100, 
      	at=seq(-my_value,my_value,my_step),
      	xlab="",
     	ylab="")
} 

from_obs_to_exp<-function(obs) {
	exp<-obs
	for (i in 1:dim(obs)[1]){
		for (j in 1:dim(obs)[2]){ 
			exp[i,j]<-(colSums(obs)[j]/sum(obs))*sum(obs)*(rowSums(obs)[i]/sum(obs))
		}
	}
	return(exp)
}

from_obs_and_exp_to_sdev<-function(ala_obs,ala_exp) {
	ala_sdev<-ala_obs
	for (i in 1:dim(ala_obs)[1]){
		for (j in 1:dim(ala_obs)[2]){ 
			ala_sdev[i,j]<-(ala_obs[i,j]-ala_exp[i,j])/sqrt(ala_exp[i,j]*(1-colSums(ala_obs)[j]/sum(ala_obs)*rowSums(ala_obs)[i]/sum(ala_obs)))
		}
	}
	return(ala_sdev)
}

from_obs_to_rscpu_way1<-function(ala_obs) {
	ala_exp_prob<-ala_obs
	my_dim<-dim(ala_obs)[1]*dim(ala_obs)[2]
	for (i in 1:dim(ala_obs)[1]){
		for (j in 1:dim(ala_obs)[2]){ 
			ala_obs_prob<-ala_obs/sum(ala_obs)
			ala_exp_prob[i,j]<-((colSums(ala_obs)[j]/sum(ala_obs))*(rowSums(ala_obs)[i]/sum(ala_obs)))			
		}
	}
	ala_rscpu1<-ala_obs_prob/ala_exp_prob
	return(ala_rscpu1)
}

from_obs_to_rscpu_way2<-function(ala_obs) {
	ala_exp_prob<-ala_obs
	my_dim<-dim(ala_obs)[1]*dim(ala_obs)[2]
	for (i in 1:dim(ala_obs)[1]){
		for (j in 1:dim(ala_obs)[2]){ 
			ala_obs_prob<-ala_obs/sum(ala_obs)
			ala_exp_prob[i,j]<-((colSums(ala_obs)[j]/sum(ala_obs))*(rowSums(ala_obs)[i]/sum(ala_obs)))			
		}
	}
	ala_rscpu2<-ala_obs_prob/(1/my_dim)
	return(ala_rscpu2)
}


########################################
### Load data and get OBS, EXP, SDEV
########################################


matrix_obs<-from_codonFIXED_to_matrix_FINAL("Pyrococcus_furiosus.fa.codon_counts.merged.FIXED")

# Expected
matrix_exp<-from_obs_to_exp(matrix_obs)

# Sdev
matrix_sdev<-from_obs_and_exp_to_sdev(matrix_obs,matrix_exp)


# Plots
pdf(file="Pyrococcus_furiosus.fa.codon_counts.merged.FIXED_plots.pdf", height=7, width=7)
plot_matrix_levelplot(matrix_obs,"OBS")
plot_matrix_levelplot(matrix_exp,"EXP")
plot_matrix_sdev(matrix_sdev,"SDEV")
dev.off()

########################################
### Subdivide By AA
########################################

ala_obs<-matrix_obs[1:4,1:4]
arg_obs<-matrix_obs[5:10,5:10]
gly_obs<-matrix_obs[21:24,21:24]
ile_obs<-matrix_obs[27:29,27:29]
leu_obs<-matrix_obs[30:35,30:35]
pro_obs<-matrix_obs[41:44,41:44]
ser_obs<-matrix_obs[45:50,45:50]
thr_obs<-matrix_obs[51:54,51:54]
val_obs<-matrix_obs[58:61,58:61]

ala_sdev<-from_obs_and_exp_to_sdev(ala_obs,from_obs_to_exp(ala_obs))
arg_sdev<-from_obs_and_exp_to_sdev(arg_obs,from_obs_to_exp(arg_obs))
gly_sdev<-from_obs_and_exp_to_sdev(gly_obs,from_obs_to_exp(gly_obs))
ile_sdev<-from_obs_and_exp_to_sdev(ile_obs,from_obs_to_exp(ile_obs))
leu_sdev<-from_obs_and_exp_to_sdev(leu_obs,from_obs_to_exp(leu_obs))
pro_sdev<-from_obs_and_exp_to_sdev(pro_obs,from_obs_to_exp(pro_obs))
ser_sdev<-from_obs_and_exp_to_sdev(ser_obs,from_obs_to_exp(ser_obs))
thr_sdev<-from_obs_and_exp_to_sdev(thr_obs,from_obs_to_exp(thr_obs))
val_sdev<-from_obs_and_exp_to_sdev(val_obs,from_obs_to_exp(val_obs))

write.table(as.vector(t(arg_sdev)), "Pyrococcus_furiosus.fa.codon_counts.merged.FIXED_arg_sdev.vector.txt", row.names=F)

sink("Pyrococcus_furiosus.fa.codon_counts.merged.FIXED_sdevs.txt",append=FALSE, split=FALSE)
ala_sdev
arg_sdev
gly_sdev
ile_sdev
leu_sdev
ser_sdev
pro_sdev
thr_sdev
val_sdev
sink()


pdf(file="Pyrococcus_furiosus.fa.codon_counts.merged.FIXED_plots_byAA.pdf", height=7, width=7)
plot_matrix_sdev(ala_sdev,"Ala")
plot_matrix_sdev(arg_sdev,"Arg")
plot_matrix_sdev(gly_sdev,"Gly")
plot_matrix_sdev(ile_sdev,"Ile")
plot_matrix_sdev(leu_sdev,"Leu")
plot_matrix_sdev(pro_sdev,"Pro")
plot_matrix_sdev(ser_sdev,"Ser")
plot_matrix_sdev(thr_sdev,"Thr")
plot_matrix_sdev(val_sdev,"Val")
dev.off()


########################################
### RSCPU --> P(obs)/P(exp)
########################################

##Relative synonymous codon usage
#RSCU(codon)= Obs frequency of codon(compared to total for that AA)/ Exp frequency of codon(compared to total for that AA).
#Ex. Ala(AGC)= obs frequency/0.25

##Relative synonumous codon pair usage
#RSCU(codon pair) = Obs frequency of codon(compared to total for that AA)/ Exp frequency of codon(compared to total for that AA
#Ex. Ala(AGC-AGC)= obs frequency/0.0625

## Way 1 --> "real expected" based on the multiplication of the individual two codons
## Way 2 --> "mathematically expected", based on the equal probabilities of codons (Ex. 0.25*0.25 --> 0.0625)
########################################

#ala_rscpu1_matrix<-from_obs_to_rscpu_way1(ala_obs)
#arg_rscpu1_matrix<-from_obs_to_rscpu_way1(arg_obs)
#gly_rscpu1_matrix<-from_obs_to_rscpu_way1(gly_obs)
#ile_rscpu1_matrix<-from_obs_to_rscpu_way1(ile_obs)
#leu_rscpu1_matrix<-from_obs_to_rscpu_way1(leu_obs)
#ser_rscpu1_matrix<-from_obs_to_rscpu_way1(ser_obs)
#pro_rscpu1_matrix<-from_obs_to_rscpu_way1(pro_obs)
#thr_rscpu1_matrix<-from_obs_to_rscpu_way1(thr_obs)
#val_rscpu1_matrix<-from_obs_to_rscpu_way1(val_obs)

#ala_rscpu2_matrix<-from_obs_to_rscpu_way2(ala_obs)
#arg_rscpu2_matrix<-from_obs_to_rscpu_way2(arg_obs)
#gly_rscpu2_matrix<-from_obs_to_rscpu_way2(gly_obs)
#ile_rscpu2_matrix<-from_obs_to_rscpu_way2(ile_obs)
#leu_rscpu2_matrix<-from_obs_to_rscpu_way2(leu_obs)
#ser_rscpu2_matrix<-from_obs_to_rscpu_way2(ser_obs)
#pro_rscpu2_matrix<-from_obs_to_rscpu_way2(pro_obs)
#thr_rscpu2_matrix<-from_obs_to_rscpu_way2(thr_obs)
#val_rscpu2_matrix<-from_obs_to_rscpu_way2(val_obs)

all_rscpu1_matrix<-from_obs_to_rscpu_way1(matrix_obs)
all_rscpu2_matrix<-from_obs_to_rscpu_way2(matrix_obs)

# Convert into appended vectors, to have one vector per species
#ala_rscpu1_vector<-as.vector(t(ala_rscpu1_matrix))
#arg_rscpu1_vector<-as.vector(t(arg_rscpu1_matrix))
#gly_rscpu1_vector<-as.vector(t(gly_rscpu1_matrix))
#ile_rscpu1_vector<-as.vector(t(ile_rscpu1_matrix))
#leu_rscpu1_vector<-as.vector(t(leu_rscpu1_matrix))
#ser_rscpu1_vector<-as.vector(t(ser_rscpu1_matrix))
#pro_rscpu1_vector<-as.vector(t(pro_rscpu1_matrix))
#thr_rscpu1_vector<-as.vector(t(thr_rscpu1_matrix))
#val_rscpu1_vector<-as.vector(t(val_rscpu1_matrix))

#ala_rscpu2_vector<-as.vector(t(ala_rscpu2_matrix))
#arg_rscpu2_vector<-as.vector(t(arg_rscpu2_matrix))
#gly_rscpu2_vector<-as.vector(t(gly_rscpu2_matrix))
#ile_rscpu2_vector<-as.vector(t(ile_rscpu2_matrix))
#leu_rscpu2_vector<-as.vector(t(leu_rscpu2_matrix))
#ser_rscpu2_vector<-as.vector(t(ser_rscpu2_matrix))
#pro_rscpu2_vector<-as.vector(t(pro_rscpu2_matrix))
#thr_rscpu2_vector<-as.vector(t(thr_rscpu2_matrix))
#val_rscpu2_vector<-as.vector(t(val_rscpu2_matrix))

library("gdata")
all_rscpu1_vector<-unmatrix(all_rscpu1_matrix, byrow=FALSE)
all_rscpu2_vector<-unmatrix(all_rscpu2_matrix, byrow=FALSE)


#Pyrococcus_furiosus.fa.codon_counts.merged.FIXED<-c(ala_rscpu1_vector, arg_rscpu1_vector, gly_rscpu1_vector, ile_rscpu1_vector, leu_rscpu1_vector, ser_rscpu1_vector, pro_rscpu1_vector, thr_rscpu1_vector, val_rscpu1_vector)  # all_all_vector
Pyrococcus_furiosus.fa.codon_counts.merged.FIXED<-c(all_rscpu1_vector)
write.table(data.frame(Pyrococcus_furiosus.fa.codon_counts.merged.FIXED), "Pyrococcus_furiosus.fa.codon_counts.merged.FIXED_rscpu_vector_way1.txt", row.names=T,quote=FALSE)

#Pyrococcus_furiosus.fa.codon_counts.merged.FIXED<-c(ala_rscpu2_vector, arg_rscpu2_vector, gly_rscpu2_vector, ile_rscpu2_vector, leu_rscpu2_vector, ser_rscpu2_vector, pro_rscpu2_vector, thr_rscpu2_vector, val_rscpu2_vector)  # all_all_vector
Pyrococcus_furiosus.fa.codon_counts.merged.FIXED<-c(all_rscpu2_vector)
write.table(data.frame(Pyrococcus_furiosus.fa.codon_counts.merged.FIXED), "Pyrococcus_furiosus.fa.codon_counts.merged.FIXED_rscpu_vector_way2.txt", row.names=T, quote=FALSE)

