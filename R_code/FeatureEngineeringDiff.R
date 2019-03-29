Get.DE.Features <- function(Dat, FSet, validComps = "AD-CONTROL", 
                            validSex = 'ALL', validModel = "Diagnosis", 
                            validTissue){


#  Dat <- read.csv('E:/SageDocs/PredictingDriverGenes/DataFiles/geneExprData.csv', 
#                  stringsAsFactors = F)
#  FSet <- read.csv('E:/SageDocs/PredictingDriverGenes/MiscFiles/amp_ad_rna_seq_feature_set.csv', 
#                   stringsAsFactors = F)
  GeneNames <- FSet$GeneID
  rm(FSet)
  
  
#  validComps <- "AD-CONTROL"
#  validTissue <- "CBE"
#  validSex <- 'ALL'
#  validModel <- "Diagnosis"
  
  AvgExpr <- rep(0,length(GeneNames))
  LogFC <- rep(0,length(GeneNames))
  LogPV <- rep(0,length(GeneNames))
  t <- rep(0,length(GeneNames))
  
  Dat <- Dat[Dat$Comparison == validComps,]
  Dat <- Dat[Dat$Tissue == validTissue,]
  Dat <- Dat[Dat$Sex == validSex,]
  Dat <- Dat[Dat$Model == validModel,]
  
  for (i in 1:length(GeneNames)){
    
    In <- (Dat$ensembl_gene_id == GeneNames[i])
    
    if (sum(In) > 0){
      In = which(In)
      In = In[1]
      AvgExpr[i] <- Dat$AveExpr[In]
      LogFC[i] <- Dat$logFC[In]
      LogPV[i] <- Dat$neg.log10.adj.P.Val[In]
      t[i] <- Dat$t[In]
    }
    
    
  }
  
  Temp <- data.frame(AvgExpr, LogFC, LogPV, t)
  colnames(Temp) <- unlist(lapply(colnames(Temp), function(x) paste(x,validTissue,sep=".")))
  rownames(Temp) <- GeneNames
  
  return(Temp)

}


Bind.All.Dfs <- function(l){
  N <- names(l)
  t <- l[[N[1]]]
  
  for (i in 2:length(l)){
    t <- cbind(t,l[[N[i]]])
  }
  
  return(t)
}

Run.All.BRs <- function(Dat,FSet,validComps = "AD-CONTROL", 
                        validSex = 'ALL', validModel = "Diagnosis"){
  
  temp <- unique(Dat$Tissue)
  print(temp)
  l <- list()
  
  for(i in 1:length(temp)){
    l[[temp[i]]] <- Get.DE.Features(Dat,FSet,validComps, validSex, 
                                  validModel,validTissue = temp[i])
  }
  
  #l2 <- Bind.All.Dfs(l)
  
  return(l)
  
}

