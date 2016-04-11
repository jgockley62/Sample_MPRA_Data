Rep1_mRNA<-read.table(file="Dropbox/Sample_MPRA_Data/GSM792095_CRE_Multi_100uM_Rep1_mRNA.counts.txt", header=T, sep="\t")
Rep1_pDNA<-read.table(file="Dropbox/Sample_MPRA_Data/GSM792096_CRE_Multi_100uM_Rep1_Plasmid.counts.txt", header=T, sep="\t")

Rep2_mRNA<-read.table(file="Dropbox/Sample_MPRA_Data/GSM792097_CRE_Multi_100uM_Rep2_mRNA.counts.txt", header=T, sep="\t")
Rep2_pDNA<-read.table(file="Dropbox/Sample_MPRA_Data/GSM792098_CRE_Multi_100uM_Rep2_Plasmid.counts.txt", header=T, sep="\t")

row.names(Rep1_mRNA)<-Rep1_mRNA[,1]
  
#MAST<-cbind(Rep1_mRNA[,1:4],Rep1_pDNA[,1:4],Rep2_mRNA[,1:4],Rep2_pDNA[,1:4])
MAST<-cbind(Rep1_mRNA[,1:4],Rep1_pDNA[,1:4],((Rep1_mRNA[,4]/sum(Rep1_mRNA[,4]))/(Rep1_pDNA[,4]/sum(Rep1_pDNA[,4]))),Rep2_mRNA[,1:4],Rep2_pDNA[,1:4],((Rep2_mRNA[,4]/sum(Rep2_mRNA[,4]))/(Rep2_pDNA[,4]/sum(Rep2_pDNA[,4]))))

CureMAST_a<-MAST[MAST[,18] != Inf ,]
CureMAST_b<-CureMAST_a[CureMAST_a[,9] != Inf ,]
CureMAST_c<-CureMAST_b[!is.na(CureMAST_b[,9]),]
CureMAST_d<-CureMAST_c[!is.na(CureMAST_c[,18]),]

#Plot reproduibility
cor(CureMAST_d[,9],CureMAST_d[,18])

setwd("~/Dropbox/Sample_MPRA_Data/Down_Sampling/")
setEPS()
postscript("100Percent.eps")

m <- rbind(c(1, 1), c(2, 3))
print(m)
layout(m)

plot(CureMAST_d[,9], CureMAST_d[,18], pch =16, cex =.1, xlab="Replicate 1", ylab="Replicate 2", las=1, bty="n")
abline(lm(CureMAST_d[,18]~CureMAST_d[,9]), col="red")
mtext(cor(CureMAST_d[,9],CureMAST_d[,18]), side=3)

hist(CureMAST_d[,9], breaks=200, las=1, bty="n", main="Rep 1", xlab="Fold Enrichment")
hist(CureMAST_d[,18], breaks=200, las=1, bty="n",  main="Rep 2", xlab="Fold Enrichment")

dev.off()

setwd("~/Dropbox/Sample_MPRA_Data/")
setEPS()
postscript("Reproduce_WholeData.eps")

plot(CureMAST_d[,9], CureMAST_d[,18], pch =16, cex =.1, xlab="Replicate 1", ylab="Replicate 2", las=1, bty="n")
abline(lm(CureMAST_d[,18]~CureMAST_d[,9]), col="red")

dev.off()

plot.new()

setwd("~/Dropbox/Sample_MPRA_Data/")
setEPS()
postscript("HistFold_WholeData.eps")

par(mfrow=c(1,2))
hist(CureMAST_d[,9], breaks=200, las=1, bty="n", main="Rep 1", xlab="Fold Enrichment")
hist(CureMAST_d[,18], breaks=200, las=1, bty="n",  main="Rep 2", xlab="Fold Enrichment")
dev.off()

##########DownSampling each pool
###Create tables with proportion of tag counts (Column 5) and a whole number value of that proportion (Column6)
#The downsampling script will use Col6 to seed an array to sample from with the correct probability for each tag

Rep1_mRNA_COUNTS<-sum(Rep1_mRNA[,4])
Rep1_mRNA_foo<-cbind(Rep1_mRNA[,1:4],Rep1_mRNA[,4]/sum(Rep1_mRNA[,4]),round(Rep1_mRNA[,4]/sum(Rep1_mRNA[,4])*Rep1_mRNA_COUNTS))

Rep1_pDNA_COUNTS<-sum(Rep1_pDNA[,4])
Rep1_pDNA_foo<-cbind(Rep1_pDNA[,1:4],Rep1_pDNA[,4]/sum(Rep1_pDNA[,4]),round(Rep1_pDNA[,4]/sum(Rep1_pDNA[,4])*Rep1_pDNA_COUNTS))


Rep2_mRNA_COUNTS<-sum(Rep2_mRNA[,4])
Rep2_mRNA_foo<-cbind(Rep2_mRNA[,1:4],Rep2_mRNA[,4]/sum(Rep2_mRNA[,4]),round(Rep2_mRNA[,4]/sum(Rep2_mRNA[,4])*Rep2_mRNA_COUNTS))


Rep2_pDNA_COUNTS<-sum(Rep2_pDNA[,4])
Rep2_pDNA_foo<-cbind(Rep2_pDNA[,1:4],Rep2_pDNA[,4]/sum(Rep2_pDNA[,4]),round(Rep2_pDNA[,4]/sum(Rep2_pDNA[,4])*Rep2_pDNA_COUNTS))

####Write XXXX_foo tables to file for input into perl downsampling script
write.table(Rep1_mRNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep1_mRNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(Rep1_pDNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep1_pDNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)

write.table(Rep2_mRNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep2_mRNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(Rep2_pDNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep2_pDNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)

#########
#Load the down sampled files
Rep1_mRNA_DS<-read.table(file="Rep1_mRNA_DownSampled.txt", header=T, sep="\t")
Rep1_pDNA_DS<-read.table(file="Rep1_pDNA_DownSampled.txt", header=T, sep="\t")

Rep2_mRNA_DS<-read.table(file="Rep2_mRNA_DownSampled.txt", header=T, sep="\t")  
Rep2_pDNA_DS<-read.table(file="Rep2_pDNA_DownSampled.txt", header=T, sep="\t")

#Calc fold increase
Rep1_DS<-data.frame(cbind(as.character(Rep1_mRNA_DS$Name), (Rep1_mRNA_DS[,6]/sum(Rep1_mRNA_DS[,6]))/(Rep1_pDNA_DS[,6]/sum(Rep1_pDNA_DS[,6])), (Rep1_mRNA_DS[,7]/sum(Rep1_mRNA_DS[,7]))/(Rep1_pDNA_DS[,7]/sum(Rep1_pDNA_DS[,7])), (Rep1_mRNA_DS[,8]/sum(Rep1_mRNA_DS[,8]))/(Rep1_pDNA_DS[,8]/sum(Rep1_pDNA_DS[,8])), (Rep1_mRNA_DS[,9]/sum(Rep1_mRNA_DS[,9]))/(Rep1_pDNA_DS[,9]/sum(Rep1_pDNA_DS[,9])), (Rep1_mRNA_DS[,10]/sum(Rep1_mRNA_DS[,10]))/(Rep1_pDNA_DS[,10]/sum(Rep1_pDNA_DS[,10])), (Rep1_mRNA_DS[,11]/sum(Rep1_mRNA_DS[,11]))/(Rep1_pDNA_DS[,11]/sum(Rep1_pDNA_DS[,11]))))
Rep2_DS<-data.frame(cbind(as.character(Rep2_mRNA_DS$Name), (Rep2_mRNA_DS[,6]/sum(Rep2_mRNA_DS[,6]))/(Rep2_pDNA_DS[,6]/sum(Rep2_pDNA_DS[,6])), (Rep2_mRNA_DS[,7]/sum(Rep2_mRNA_DS[,7]))/(Rep2_pDNA_DS[,7]/sum(Rep2_pDNA_DS[,7])), (Rep2_mRNA_DS[,8]/sum(Rep2_mRNA_DS[,8]))/(Rep2_pDNA_DS[,8]/sum(Rep2_pDNA_DS[,8])), (Rep2_mRNA_DS[,9]/sum(Rep2_mRNA_DS[,9]))/(Rep2_pDNA_DS[,9]/sum(Rep2_pDNA_DS[,9])), (Rep2_mRNA_DS[,10]/sum(Rep2_mRNA_DS[,10]))/(Rep2_pDNA_DS[,10]/sum(Rep2_pDNA_DS[,10])), (Rep2_mRNA_DS[,11]/sum(Rep2_mRNA_DS[,11]))/(Rep2_pDNA_DS[,11]/sum(Rep2_pDNA_DS[,11]))))

#Table Specific Down sampling
HUNPercent<-data.frame(cbind(as.character(Rep1_DS[,1]), as.character(Rep1_DS[,2]), as.character(Rep2_DS[,2])))
STFPercent<-data.frame(cbind(as.character(Rep1_DS[,1]), as.character(Rep1_DS[,3]), as.character(Rep2_DS[,3])))
FIFPercent<-data.frame(cbind(as.character(Rep1_DS[,1]), as.character(Rep1_DS[,4]), as.character(Rep2_DS[,4])))
TWFPercent<-data.frame(cbind(as.character(Rep1_DS[,1]), as.character(Rep1_DS[,5]), as.character(Rep2_DS[,5])))
TENPercent<-data.frame(cbind(as.character(Rep1_DS[,1]), as.character(Rep1_DS[,6]), as.character(Rep2_DS[,6])))
FIVPercent<-data.frame(cbind(as.character(Rep1_DS[,1]), as.character(Rep1_DS[,7]), as.character(Rep2_DS[,7])))

#Clean Each DS for ploting
CleanDS <- function(FRAME){
  FRAME<-FRAME[FRAME[,2] != Inf ,]
  FRAME<-FRAME[FRAME[,3] != Inf ,]
  FRAME<-FRAME[!is.na(FRAME[,2]),]
  FRAME<-FRAME[!is.na(FRAME[,3]),]
  return(FRAME)
}

HUNPercent<-CleanDS(HUNPercent)
STFPercent<-CleanDS(STFPercent)
FIFPercent<-CleanDS(FIFPercent)
TWFPercent<-CleanDS(TWFPercent)
TENPercent<-CleanDS(TENPercent)
FIVPercent<-CleanDS(FIVPercent)


setwd("~/Dropbox/Sample_MPRA_Data/Down_Sampling/")
setEPS()
postscript("100Percent.eps")

m <- rbind(c(1, 1), c(2, 3))
print(m)
layout(m)

plot(HUNPercent[,2], HUNPercent[,3], pch =16, cex =.1, xlab="Replicate 1", ylab="Replicate 2", las=1, bty="n")
abline(lm(HUNPercent[,3]~HUNPercent[,2]), col="red")
mtext(cor(HUNPercent[,2],HUNPercent[,3]), side=3)

hist(CureMAST_d[,9], breaks=200, las=1, bty="n", main="Rep 1", xlab="Fold Enrichment")
hist(CureMAST_d[,18], breaks=200, las=1, bty="n",  main="Rep 2", xlab="Fold Enrichment")

dev.off()









