Process.Feature <- function(X){
  
  #Save the gene names and drop the gene name column 
  GeneId <- X$GeneID
  X$GeneID <- NULL
  X$X <- NULL
  
  #Impute missing values with 0 
  X[is.na(X)] <- 0 
  
  #scale columns to between -1 and 1 
  d <- dim(X)[2]
  
  for(i in 1:d){
    X[,i] <- X[,i]/max(abs(X[,i]))
  }
  
  return(X)
}

Get.Gene.Symbs <- function(GeneList,NetDat){
  
  temp <- c()
  
  for (i in 1:length(GeneList)){
    In <- which(NetDat$GeneID==GeneList[i])
    if (length(In)==0){
      temp <- c(temp,str(i))
    }else{
      temp <- c(temp,NetDat$external_gene_name[In[1]])
    }
  }
  return(temp)
}

Get.Pval <- function(GeneList,IGAP_list,mode = 'min'){
  
  pv <- c()
  for(i in 1:length(GeneList)){
    In <- which(IGAP_list$Names==GeneList[i])
    if(length(In)==0){
      pv <- c(pv,1)
    } else{
      if(mode=='min'){
        pv <- c(pv,IGAP_list$Min[In[1]])
      }else{
        pv <- c(pv,IGAP_list$Mean[In[1]])
      }
    }
  }
  
  return(pv)
}