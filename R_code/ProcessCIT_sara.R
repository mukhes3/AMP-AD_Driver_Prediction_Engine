CIT <- read.delim('../DataFiles/CIT.txt', stringsAsFactors = F)
GL <- unique(CIT$gene)
GL_e <- convertHgncToEnsembl(GL)

CIT$ENSG <- CIT$gene



for ( i in 1:length(GL_e$ensembl_gene_id)){
  In <- which(CIT$gene == GL_e$external_gene_name[i])
  
  if(length(In)>0){
    CIT$ENSG[In] = GL_e$ensembl_gene_id[i]
  }
}


CIT <- read.csv('../DataFiles/CIT_ensg.csv')
GL <- unique(CIT$ENSG)
ENSG <- CIT$ENSG
CIT$X <- NULL
CIT$snp <- NULL
CIT$nProbe <- NULL
CIT$nPeak <- NULL
CIT$probes <- NULL
CIT$gene <- NULL
CIT$ENSG <- NULL
CIT$peaks <- NULL


GenFeatCIT <- function(GL,ENSG,x,sfx){
  
  library(moments)
  
  mn <- rep(0,length(GL))
  n <- rep(0,length(GL))
  m <- rep(0,length(GL))
  v <- rep(0,length(GL))
  std <- rep(0,length(GL))
  s <- rep(0,length(GL))
  x[x==0] = min(x[x>0])
  
  for(i in 1:length(GL)){
    
    In <- which(ENSG==GL[i])
    n[i] <- length(In)
    #print(length(GL))
    #print(x[In])
    m[i] <- mean(-log10(x[In]))
    mn[i] <- min(-log10(x[In]))
    
    
    if(n[i]>1){
      
      v[i] <- var(-log10(x[In]))
      std[i] <- sd(-log10(x[In]))
      s[i] <- skewness(-log10(x[In]))
      
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
  
  colnames(l) <- paste(colnames(l),sfx,sep='')
  
  return(l)
  
}

CN <- colnames(CIT)
temp <- GenFeatCIT(GL,ENSG,CIT[[CN[1]]],CN[1])

for(i in 2:length(CN)){
  
  temp2 <- GenFeatCIT(GL,ENSG,CIT[[CN[i]]],CN[i])
  temp2 <- temp2[,2:dim(temp2)[2]]
  temp <- cbind(temp,temp2)
  
}

write.csv(temp,'CIT_feat.csv')


