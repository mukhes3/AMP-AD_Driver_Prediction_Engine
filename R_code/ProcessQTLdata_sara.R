#download.file("http://mostafavilab.stat.ubc.ca/xqtl/eQTLs.txt",destfile = '../DataFiles/eQTLS_sara.txt')
#mQTL <- read.delim('../DataFiles/eQTLS_sara.txt', stringsAsFactors = F)
mQTL <- read.delim('../DataFiles/SWB_Okbay_Full.txt', stringsAsFactors = F)
mQTL$SNPid <- mQTL$MarkerName
mQTL$pValue <- mQTL$Pval
source('convertSnpsToGenes.R')

c <- 1000

Dat <- convertSnpsToGenes(mQTL$SNPid[1:c])
#Dat <- convertSnpsToGenes(mQTL$SNP[1:c])
Dat <- Dat[Dat$distance_to_transcript<3000,]

for (i in 2:(length(mQTL$SNPid)/c)){
  
  temp <- convertSnpsToGenes(mQTL$SNPid[(i-1)*c+1:i*c])
  Dat <- rbind(Dat,temp[temp$distance_to_transcript<10000,])
  
  if (i%%10 == 0){
    print(i)
  }
  
}

RmNAs <- function(Dat){
  
  Dat$distance_to_transcript <- NULL
  Dat$sift_score <- NULL
  Dat$clinical_significance <- NULL
  Dat$polyphen_score <- NULL
  Dat$phenotype_name <- NULL
  Dat$p_value <- NULL
  In <- is.na(Dat$refsnp_id)
  Dat <- Dat[which(!In),]
  
  return(Dat)

}


Dat <- RmNAs(Dat)

#Dat <- Dat_mQTL

#Read list of ENSG genes 
GeneID <- read.csv('ResponseVec_040318.csv', stringsAsFactors = F)
GeneID <- GeneID$GeneID

GenFeatStats <- function(Dat,mQTL){
  
  library(moments)
  GL <- unique(Dat$ensembl_gene_stable_id)
  
  mn <- rep(0,length(GL))
  m <- rep(0,length(GL))
  v <- rep(0,length(GL))
  std <- rep(0,length(GL))
  s <- rep(0,length(GL))
  n <- rep(0,length(GL))
  
  mQTL$pValue[mQTL$pValue==0] = min(mQTL$pValue[mQTL$pValue>0])


  for(i in 1:length(GL)){
    
    In1 <- which(Dat$ensembl_gene_stable_id==GL[i])
    In2 <- which(mQTL$SNPid %in% unique(Dat$refsnp_id[In1]))
    n[i] <- length(In2)
    m[i] <- mean(-log10(mQTL$pValue[In2]))
    mn[i] <- max(-log10(mQTL$pValue[In2]))
    
    if (n[i]>1){
      
      v[i] <- var(-log10(mQTL$pValue[In2]))
      std[i] <- sd(-log10(mQTL$pValue[In2]))
      s[i] <- skewness(-log10(mQTL$pValue[In2]))
    }
  }
  
  l <- list()
  l$GL <- GL
  l$m <- m
  l$v <- v
  l$std <- std
  l$s <- s
  l$n <- n
  l$mn <- mn 
  
  l <- as.data.frame(l)
  
  return(l)
}


Dat2 <- GenFeatStats(Dat,mQTL)


