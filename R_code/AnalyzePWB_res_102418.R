Dat <- read.csv('../DataFiles/WingoGWAS_sexPC_summary.csv',stringsAsFactors = F)
Y_pred <- read.csv('./EC2_ConsUnion1_NoDeProbs_LS.csv',stringsAsFactors = F)
Dat$gene <- Dat$GL
Dat$pval_affy <- Dat$mn

source('convertEnsemblToHgnc.R')
#Yconv = convertEnsemblToHgnc(Y_pred$Gene)

ResTT <- function(GL,Yconv,Dat, ENSG_conv = T){
  
  if (ENSG_conv){
    GL2 <- unique(GL)
    In1 <- which(Dat$gene %in% GL2)
    print(length(In1))
    In0 <- which(!(Dat$gene %in% GL2))
    print(length(In0))
    
  } else {
    In <- which(Yconv$ensembl_gene_id %in% GL)
    GL2 <- unique(Yconv$external_gene_name[In])
    
    In1 <- which(Dat$gene %in% GL2)
    print(length(In1))
    In0 <- which(!(Dat$gene %in% GL2))
    print(length(In0))
  }
  
  
  
  res <- t.test(-log10(Dat$pval_affy[In1]),-log10(Dat$pval_affy[In0]))
  
  print(res$p.value)
  
  return(res$p.value)
}


#ResTT(Y_pred$Gene[which(Y_pred$Y0a>0.5)],Yconv,Dat)
ResTT(Y_pred$Gene[which(Y_pred$Y1a>0.5)],Yconv,Dat,ENSG_conv = T)
#ResTT(Y_pred$Gene[which(Y_pred$Y2a>0.5)],Yconv,Dat)
ResTT(Y_pred$Gene[which(Y_pred$Y3a>0.5)],Yconv,Dat,ENSG_conv = T)
ResTT(Y_pred$Gene[which(Y_pred$Y3a>0.5 | Y_pred$Y1a>0.5)],Yconv,Dat,ENSG_conv = T)