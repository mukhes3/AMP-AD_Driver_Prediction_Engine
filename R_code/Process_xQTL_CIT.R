setwd('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/')

Xm <- read.csv('../DataFiles/mQTL_feat.csv',stringsAsFactors = F)
Xm$X <- NULL
colnames(Xm)[2:dim(Xm)[2]] <- paste(colnames(Xm)[2:dim(Xm)[2]],'_m',sep='')

Xh <- read.csv('../DataFiles/hQTL_feat.csv',stringsAsFactors = F)
Xh$X <- NULL
colnames(Xh)[2:dim(Xh)[2]] <- paste(colnames(Xh)[2:dim(Xh)[2]],'_h',sep='')


Xe <- read.csv('../DataFiles/eQTL_feat.csv',stringsAsFactors = F)
Xe$X <- NULL
colnames(Xe)[2:dim(Xe)[2]] <- paste(colnames(Xe)[2:dim(Xe)[2]],'_e',sep='')

Xc <- read.csv('../DataFiles/CIT_feat.csv',stringsAsFactors = F)
Xc$X <- NULL
Xc$GL <- Xc$GLpCausal
Xc$GLpCausal <- NULL

#merge all dataframes 
Dat <- merge(Xm,Xh,by = 'GL' ,all.x = T,all.y = T)
Dat <- merge(Dat,Xe,by = 'GL' ,all.x = T,all.y = T)
Dat <- merge(Dat,Xc,by = 'GL' ,all.x = T,all.y = T)
Dat[is.na(Dat)] <- 0.0


#Create a larger merged dataset with all other genes as blank
GeneID <- read.csv('ResponseVec_040318.csv', stringsAsFactors = F)
GeneID$Y <- NULL
GeneID$X <- NULL
GeneID$X.1 <- NULL
GeneID$X.2 <- NULL
GeneID$GL <- GeneID$GeneID
GeneID$GeneID <- NULL

Dat2 <- merge(Dat,GeneID,by = 'GL',all.y = T)
Dat2[is.na(Dat2)] = 0.0
Dat2$GeneID <- Dat2$GL
Dat2$GL <- NULL

write.csv(Dat2,'amp_ad_xQTL_CIT_feature_set_DLPFC.csv')

