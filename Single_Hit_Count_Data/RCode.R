#Single Hit Tag Code

Rep1<-read.table(file="~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/Tags_Replicate_1.txt", header = T, sep = "\t")
Rep2<-read.table(file="~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/Tags_Replicate_2.txt", header = T, sep = "\t")

Rep1_Sig<-cbind(Rep1[,1:4], (Rep1[,3]/sum(Rep1[,3]))/(Rep1[,4]/sum(Rep1[,4])))
colnames(Rep1_Sig)<- c(colnames(Rep1), "Signal")
Rep2_Sig<-cbind(Rep1[,1:4], (Rep2[,3]/sum(Rep2[,3]))/(Rep2[,4]/sum(Rep2[,4])))
colnames(Rep2_Sig)<- c(colnames(Rep2), "Signal")

#Clean Each DS for ploting
CleanDS <- function(FRAME){
  FRAME<-FRAME[as.numeric(as.character(FRAME[,5])) != Inf ,]
  FRAME<-FRAME[as.numeric(as.character(FRAME[,10])) != Inf ,]
  FRAME<-FRAME[as.numeric(as.character(FRAME[,5])) > 0 ,]
  FRAME<-FRAME[as.numeric(as.character(FRAME[,10])) > 0 ,]
  FRAME<-FRAME[!is.na(as.numeric(as.character(FRAME[,5]))),]
  FRAME<-FRAME[!is.na(as.numeric(as.character(FRAME[,10]))),]
  return(FRAME)
}

Total<-cbind(Rep1_Sig[,1:5], Rep2_Sig[,1:5])
Final_Total<-CleanDS(Total)

cor(Final_Total[,5],Final_Total[,10])
#0.7933779
plot(Final_Total[,5],Final_Total[,10], pch=16, cex=.1)

#########
###Look at Tag Variation

Frag_mRNA_Rep1<-read.table(file="~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/SingleHit_mRNA_Rep1_Frag.txt", header=T, sep="\t", row.names=1)
Frag_mRNA_Rep2<-read.table(file="~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/SingleHit_mRNA_Rep2_Frag.txt", header=T, sep="\t", row.names=1)

Frag_pDNA_Rep1<-read.table(file="~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/SingleHit_pDNA_Rep1_Frag.txt", header=T, sep="\t", row.names=1)
Frag_pDNA_Rep2<-read.table(file="~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/SingleHit_pDNA_Rep2_Frag.txt", header=T, sep="\t", row.names=1)


CompileDS <- function(FRAME){
  FRAME<-as.matrix(FRAME)
  FRAME[,1:13]<-as.numeric(FRAME[,1:13])
  FRAME<-cbind(FRAME[,1:13], apply(FRAME, 1, mean), apply(FRAME, 1, sd))
  FRAME<-cbind(FRAME[,1:15], FRAME[,14]/sum(FRAME[,14]))
  return(FRAME)
}

Frag_mRNA_Rep1<-CompileDS(Frag_mRNA_Rep1)
Frag_mRNA_Rep2<-CompileDS(Frag_mRNA_Rep2)
Frag_pDNA_Rep1<-CompileDS(Frag_pDNA_Rep1)
Frag_pDNA_Rep2<-CompileDS(Frag_pDNA_Rep2)

plot(Frag_mRNA_Rep1[,16]/Frag_pDNA_Rep1[,16], Frag_mRNA_Rep2[,16]/Frag_pDNA_Rep2[,16], pch =16, cex=.3)
cor(Frag_mRNA_Rep1[,16]/Frag_pDNA_Rep1[,16], Frag_mRNA_Rep2[,16]/Frag_pDNA_Rep2[,16])


###0.5 SD
foo<-Frag_mRNA_Rep1
for(r in 1:999){
  for(c in 1:13){
    if(foo[r,c] > (foo[r,14]-.5*foo[r,15]) & foo[r,c] < (foo[r,14]+.5*foo[r,15])){
      
    }else{
      foo[r,c]<-NA
    }
    
  }  
  
}
#Loose 7693/12987 tags total (59.2%) THIS IS a bit LOW exp = 62.8%



###ONE SD
foo<-Frag_mRNA_Rep1
for(r in 1:999){
  for(c in 1:13){
    if(foo[r,c] > (foo[r,14]-foo[r,15]) & foo[r,c] < (foo[r,14]+foo[r,15])){
      
    }else{
      foo[r,c]<-NA
    }
    
  }  
  
}
#Loose 3182/12987 tags total (24.5%) THIS IS LOW exp = 32%

###TWO SD
foo<-Frag_mRNA_Rep1
for(r in 1:999){
  for(c in 1:13){
    if(foo[r,c] > (foo[r,14]-2*foo[r,15]) & foo[r,c] < (foo[r,14]+2*foo[r,15])){
      
    }else{
      foo[r,c]<-NA
    }
    
  }  
  
}
#Loose 758/12987 tags total (5.8%) THIS IS a bit high exp = 4.46%

###Three SD
foo<-Frag_mRNA_Rep1
for(r in 1:999){
  for(c in 1:13){
    if(foo[r,c] > (foo[r,14]-3*foo[r,15]) & foo[r,c] < (foo[r,14]+3*foo[r,15])){
      
    }else{
      foo[r,c]<-NA
    }
    
  }  
  
}
#Loose 88/12987 tags total (0.6%) THIS IS high exp = 0.27%


#Look at correlation after filtering 2SD's
##Filter function to set values above or below 2SD's of tag means to zero
FILTER<- function(FRAME){
  for(r in 1:999){
    for(c in 1:13){
      if(FRAME[r,c] > (FRAME[r,14]-2*FRAME[r,15]) & FRAME[r,c] < (FRAME[r,14]+2*FRAME[r,15])){
        
      }else{
        FRAME[r,c]<-NA
      }
      
    }  
    
  }
  for(r in 1:999){
    FRAME[r,14]<-mean(na.omit(FRAME[r,1:13]))
    FRAME[r,15]<-sd(na.omit(FRAME[r,1:13]))
    
  }
  Total<-sum(FRAME[,14])
  for(r in 1:999){
    FRAME[r,16]<-FRAME[r,14]/Total    
  }
  return(FRAME)
}



Frag_mRNA_Rep1_filt<-FILTER(Frag_mRNA_Rep1)
Frag_pDNA_Rep1_filt<-FILTER(Frag_pDNA_Rep1)

Frag_mRNA_Rep2_filt<-FILTER(Frag_mRNA_Rep2)
Frag_pDNA_Rep2_filt<-FILTER(Frag_pDNA_Rep2)

plot(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16], pch =16, cex=.3)
cor(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16])


##Filter function to set values above or below 3SD's of tag means to zero
FILTER<- function(FRAME, Num){
  for(r in 1:999){
    for(c in 1:13){
      if(FRAME[r,c] > (FRAME[r,14]-Num*FRAME[r,15]) & FRAME[r,c] < (FRAME[r,14]+Num*FRAME[r,15])){
        
      }else{
        FRAME[r,c]<-NA
      }
      
    }  
    
  }
  for(r in 1:999){
    FRAME[r,14]<-mean(na.omit(FRAME[r,1:13]))
    FRAME[r,15]<-sd(na.omit(FRAME[r,1:13]))
    
  }
  Total<-sum(FRAME[,14])
  for(r in 1:999){
    FRAME[r,16]<-FRAME[r,14]/Total    
  }
  return(FRAME)
}

setwd("~/Dropbox/Sample_MPRA_Data/Single_Hit_Count_Data/")
setEPS()
postscript("Correlation_by_Fragment.eps")

par(mfrow=c(6,1), mar=c(2, 3, 2, 1))

#Without Filtering:
plot(Frag_mRNA_Rep1[,16]/Frag_pDNA_Rep1[,16], Frag_mRNA_Rep2[,16]/Frag_pDNA_Rep2[,16], pch =16, cex=.3, bty='n', xlab = "Rep1", ylab="Rep2", main="All Data")
tmp1<-as.numeric(Frag_mRNA_Rep2[,16])/as.numeric(Frag_pDNA_Rep2[,16])
tmp2<-as.numeric(Frag_mRNA_Rep1[,16])/as.numeric(Frag_pDNA_Rep1[,16])
abline(lm(tmp1~tmp2), col="red")
mtext(round(cor(Frag_mRNA_Rep1[,16]/Frag_pDNA_Rep1[,16], Frag_mRNA_Rep2[,16]/Frag_pDNA_Rep2[,16]),3), adj=1, line=0)

Frag_mRNA_Rep1_filt<-FILTER(Frag_mRNA_Rep1,1)
Frag_pDNA_Rep1_filt<-FILTER(Frag_pDNA_Rep1,1)
Frag_mRNA_Rep2_filt<-FILTER(Frag_mRNA_Rep2,1)
Frag_pDNA_Rep2_filt<-FILTER(Frag_pDNA_Rep2,1)

plot(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16], pch =16, cex=.3, bty='n', xlab = "Rep1", ylab="Rep2", main="1 SD")
tmp1<-as.numeric(Frag_mRNA_Rep2_filt[,16])/as.numeric(Frag_pDNA_Rep2_filt[,16])
tmp2<-as.numeric(Frag_mRNA_Rep1_filt[,16])/as.numeric(Frag_pDNA_Rep1_filt[,16])
abline(lm(tmp1~tmp2), col="red")
mtext(round(cor(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16]),3), adj=1, line=0)

Frag_mRNA_Rep1_filt<-FILTER(Frag_mRNA_Rep1,2)
Frag_pDNA_Rep1_filt<-FILTER(Frag_pDNA_Rep1,2)
Frag_mRNA_Rep2_filt<-FILTER(Frag_mRNA_Rep2,2)
Frag_pDNA_Rep2_filt<-FILTER(Frag_pDNA_Rep2,2)

plot(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16], pch =16, cex=.3, bty='n', xlab = "Rep1", ylab="Rep2", main="2 SDs")
tmp1<-as.numeric(Frag_mRNA_Rep2_filt[,16])/as.numeric(Frag_pDNA_Rep2_filt[,16])
tmp2<-as.numeric(Frag_mRNA_Rep1_filt[,16])/as.numeric(Frag_pDNA_Rep1_filt[,16])
abline(lm(tmp1~tmp2), col="red")
mtext(round(cor(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16]),3), adj=1, line=0)

Frag_mRNA_Rep1_filt<-FILTER(Frag_mRNA_Rep1,3)
Frag_pDNA_Rep1_filt<-FILTER(Frag_pDNA_Rep1,3)
Frag_mRNA_Rep2_filt<-FILTER(Frag_mRNA_Rep2,3)
Frag_pDNA_Rep2_filt<-FILTER(Frag_pDNA_Rep2,3)

plot(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16], pch =16, cex=.3, bty='n', xlab = "Rep1", ylab="Rep2", main="3 SDs")
tmp1<-as.numeric(Frag_mRNA_Rep2_filt[,16])/as.numeric(Frag_pDNA_Rep2_filt[,16])
tmp2<-as.numeric(Frag_mRNA_Rep1_filt[,16])/as.numeric(Frag_pDNA_Rep1_filt[,16])
abline(lm(tmp1~tmp2), col="red")
mtext(round(cor(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16]),3), adj=1, line=0)

Frag_mRNA_Rep1_filt<-FILTER(Frag_mRNA_Rep1,4)
Frag_pDNA_Rep1_filt<-FILTER(Frag_pDNA_Rep1,4)
Frag_mRNA_Rep2_filt<-FILTER(Frag_mRNA_Rep2,4)
Frag_pDNA_Rep2_filt<-FILTER(Frag_pDNA_Rep2,4)

plot(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16], pch =16, cex=.3, bty='n', xlab = "Rep1", ylab="Rep2", main="4 SDs")
tmp1<-as.numeric(Frag_mRNA_Rep2_filt[,16])/as.numeric(Frag_pDNA_Rep2_filt[,16])
tmp2<-as.numeric(Frag_mRNA_Rep1_filt[,16])/as.numeric(Frag_pDNA_Rep1_filt[,16])
abline(lm(tmp1~tmp2), col="red")
mtext(round(cor(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16]),3), adj=1, line=0)

Frag_mRNA_Rep1_filt<-FILTER(Frag_mRNA_Rep1,5)
Frag_pDNA_Rep1_filt<-FILTER(Frag_pDNA_Rep1,5)
Frag_mRNA_Rep2_filt<-FILTER(Frag_mRNA_Rep2,5)
Frag_pDNA_Rep2_filt<-FILTER(Frag_pDNA_Rep2,5)

plot(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16], pch =16, cex=.3, bty='n', xlab = "Rep1", ylab="Rep2", main="5 SDs")
tmp1<-as.numeric(Frag_mRNA_Rep2_filt[,16])/as.numeric(Frag_pDNA_Rep2_filt[,16])
tmp2<-as.numeric(Frag_mRNA_Rep1_filt[,16])/as.numeric(Frag_pDNA_Rep1_filt[,16])
abline(lm(tmp1~tmp2), col="red")
mtext(round(cor(Frag_mRNA_Rep1_filt[,16]/Frag_pDNA_Rep1_filt[,16], Frag_mRNA_Rep2_filt[,16]/Frag_pDNA_Rep2_filt[,16]),3), adj=1, line=0)

dev.off()















