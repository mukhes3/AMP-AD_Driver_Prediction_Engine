setwd('Documents/DriverPrediction/MiscFiles/')

#creating file to process feature files 
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

#Read and process different feature files
X0 = Process.Feature(read.csv('./amp_ad_de_feature_set.csv'))
X1 = Process.Feature(read.csv('./amp_ad_agg_feature_set.csv'))
X2 = Process.Feature(read.csv('./amp_ad_deVal_feature_set.csv'))
X3 = Process.Feature(read.csv('./amp_ad_subNet_feature_set.csv'))

#Read original response vector 
Y = read.csv('ResponseVec_040318.csv',stringsAsFactors = F)
GeneId <- Y$GeneID
Y <- Y$Y



#perform principal component analysis on scales features 
PC0 <- prcomp(X0, center = TRUE)
PC1 <- prcomp(X1, center = TRUE)
PC2 <- prcomp(X2, center = TRUE)
PC3 <- prcomp(X3, center = TRUE)


#make scatter plot of principal components 
par(mfrow = c(2,2))
plot(PC0$x[,1],PC0$x[,2], col = Y+1, pch=16, cex = Y + .1
     , main = 'DE')
plot(PC1$x[,1],PC1$x[,2], col = Y+1, pch=16, cex = Y + .1
     , main = 'Mod')
plot(PC2$x[,1],PC2$x[,2], col = Y+1, pch=16, cex = Y + .1,
     , main = 'DE_val')
plot(PC3$x[,1],PC3$x[,2], col = Y+1, pch=16, cex = Y + .1, 
     , main = 'Global_net')

#Get predicted response vector
Y_pred_df <- read.csv('EC2_ConsUnion1_NoDeProbs_0814.csv')
Y_pred <- (Y_pred_df$Y1a + Y_pred_df$Y1b + 
             Y_pred_df$Y3a + Y_pred_df$Y3b)/4.0

#Get neighbors of other genes in coexpression module 
NetDat <- read.csv('AggregateMods_March2018.csv',stringsAsFactors = F)
br <- 'DLPFC'
GeneName <- GeneId[15]
In <- which((NetDat$brainRegion==br) & (NetDat$GeneID==GeneName))
In1 <- which((NetDat$ModuleName == NetDat$ModuleName[In]) &
               (NetDat$GeneID != GeneName))
In2 <- which((GeneId %in% NetDat$GeneID[In1]) & (Y_pred>0.5))

#visualize the neighbors of the coexpression module on graph
IGAP_list <- read.csv('IGAP_gene_summary.csv')
library(igraph)
par(mfrow = c(1,1))
g_obj <- make_star(length(In2)+1, mode = 'undirected')
V_list <- c(GeneName, GeneId[In2])
V_list <- Get.Gene.Symbs(V_list,NetDat)
V_size <- -log10(Get.Pval(V_list,IGAP_list))
V_size[V_size>10] <- 10
V_size <- 5 + (V_size/10)*10
V(g_obj)$label <- V_list
V(g_obj)$size <-  V_size
plot(g_obj, vertex.label.cex = 0.5, main = br)


#visualize a given gene on the PCA plot
Y_pc = (Y_pred>0.5) + 0.0
Y_pc[In] <- 2.0
Y_p = (Y_pred>0.5) + 0.0
Y_p = Y_p*.5
Y_p[In] <- 2.0

par(mfrow = c(2,2))
plot(PC0$x[,1],PC0$x[,2], col = Y_pc+1, pch=16, cex = Y_p + .1
     , main = 'DE')
plot(PC1$x[,1],PC1$x[,2], col = Y_pc+1, pch=16, cex = Y_p + .1
     , main = 'Mod')
plot(PC2$x[,1],PC2$x[,2], col = Y_pc+1, pch=16, cex = Y_p + .1,
     , main = 'DE_val')
plot(PC3$x[,1],PC3$x[,2], col = Y_pc+1, pch=16, cex = Y_p + .1, 
     , main = 'Global_net')

