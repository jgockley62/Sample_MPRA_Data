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
Rep1_mRNA_foo<-cbind(Rep1_mRNA[,1:4],Rep1_mRNA[,4]/sum(Rep1_mRNA[,4]),round(Rep1_mRNA[,4]/sum(Rep1_mRNA[,4])*1e10))

Rep1_pDNA_COUNTS<-sum(Rep1_pDNA[,4])
Rep1_pDNA_foo<-cbind(Rep1_pDNA[,1:4],Rep1_pDNA[,4]/sum(Rep1_pDNA[,4]),round(Rep1_pDNA[,4]/sum(Rep1_pDNA[,4])*1e10))


Rep2_mRNA_COUNTS<-sum(Rep2_mRNA[,4])
Rep2_mRNA_foo<-cbind(Rep2_mRNA[,1:4],Rep2_mRNA[,4]/sum(Rep2_mRNA[,4]),round(Rep2_mRNA[,4]/sum(Rep2_mRNA[,4])*1e10))


Rep2_pDNA_COUNTS<-sum(Rep2_pDNA[,4])
Rep2_pDNA_foo<-cbind(Rep2_pDNA[,1:4],Rep2_pDNA[,4]/sum(Rep2_pDNA[,4]),round(Rep2_pDNA[,4]/sum(Rep2_pDNA[,4])*1e10))

####Write XXXX_foo tables to file for input into perl downsampling script
write.table(Rep1_mRNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep1_mRNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(Rep1_pDNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep1_pDNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)

write.table(Rep2_mRNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep2_mRNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)
write.table(Rep2_pDNA_foo, file = "~/Dropbox/Sample_MPRA_Data/Down_Sampling/Rep2_pDNA.txt", sep = "\t", quote=F, row.names = FALSE, col.names = FALSE)


















