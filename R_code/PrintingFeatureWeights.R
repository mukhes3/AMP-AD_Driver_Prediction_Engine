library(R.matlab)

#creating file to process feature files 
Get.Feature <- function(X){
  
  #Save the gene names and drop the gene name column 
  GeneId <- X$GeneID
  X$GeneID <- NULL
  X$X <- NULL
  
  #Impute missing values with 0 

  return(colnames(X))
}


library(ggplot2)

Get.WeightDF <- function(X,W, DispNo = 10){
  
  l <- list()
  l$Features <- X
  l$L1 <- as.vector(abs(W$L1))
  l$L2 <- as.vector(abs(W$L2))
  
  l$L1 <- l$L1[2:length(l$L1)]
  l$L2 <- l$L2[2:length(l$L2)]
  
  df <- as.data.frame(l)
  
  df <- df[order(-df$L1),]

  df <- df[1:DispNo,]
  
  df2 <- df
  
  df2$Features <- NULL
    
  counts <- as.matrix(df2)
  #print(counts)
  row.names(counts) <- df$Features
  #print(counts)


  barplot(t(counts), main="Feature weights",
          #names.arg = df$Features,
          #space = c(1,2),
          legend.text = c('L1','L2'),
          col=c("darkblue","red"),
          las=2, beside = T, cex.names = .4)
   
  
#  return(df)
  
}


X0 = Get.Feature(read.csv('./amp_ad_de_feature_set.csv'))
X1 = Get.Feature(read.csv('./amp_ad_agg_feature_set.csv'))
X2 = Get.Feature(read.csv('./amp_ad_deVal_feature_set2.csv'))
X3 = Get.Feature(read.csv('./amp_ad_subNet_feature_set.csv'))


W0 = readMat('ConsUnion_Wts_0904_x0.mat')
W1 = readMat('ConsUnion_Wts_0904_x1.mat')
W2 = readMat('ConsUnion_Wts_0904_x2.mat')
W3 = readMat('ConsUnion_Wts_0904_x3.mat')

dev.off()
#jpeg("ConsUnion_Fwt0.jpg")
Get.WeightDF(X0,W0, DispNo = 25)


#dev.off()
Get.WeightDF(X1,W1, DispNo = 25)
#jpeg("ConsUnion_Fwt1.jpg", width = 350, height = "350")

#dev.off()
Get.WeightDF(X2,W2, DispNo = 25)
#jpeg("ConsUnion_Fwt2.jpg", width = 350, height = "350")

#dev.off()
Get.WeightDF(X3,W3, DispNo = 25)
#jpeg("ConsUnion_Fwt3.jpg", width = 350, height = "350")

