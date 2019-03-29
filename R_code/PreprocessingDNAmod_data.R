Dat <- read.table('/Users/sumitmukherjee/Documents/DriverPrediction/DataFiles/CIT.txt', 
                  header = T, stringsAsFactors = F)

Dat2 <- read.csv('amp_ad_de_feature_set.csv',stringsAsFactors = F)
source('convertHgncToEnsembl.R')


GeneList <- unique(Dat$gene)
GeneTemp <- convertHgncToEnsembl(GeneList)

l <- list()
l$GeneID <- Dat2$GeneID

for (i in 1:length(Dat2$GeneID)){
  
  
}