synapser::synLogin()

#pull modules from synapse
aggregateModules <- synapser::synTableQuery("select * from syn11182793")
aggMods <-aggregateModules$asDataFrame()

#pull network from synapse
foo <- synapser::synTableQuery('select * from syn8681664 where method=\'bic\' and assay =\'RNAseq\'')
foo <- foo$asDataFrame()

loadBic <- function(synId){
  foo<-synapser::synGet(synId)
  load(foo$path)
  net <- as.matrix(bicNetworks$network)
  net <- net+t(net)
  return(net)
}

ampNetworks <- lapply(foo$id,loadBic)
names(ampNetworks) <- foo$tissueTypeAbrv

rownames(foo) <- names(ampNetworks)


nedges <- sapply(ampNetworks,function(x) sum(x))
foo$nedges <- nedges
foo5 <- foo

makeAggregateNetworks <- function(br,aggMod,aggModMani,networks){
  net <- networks[[br]]
  modDef <- dplyr::filter(aggModMani,ModuleNameFull == aggMod)$GeneID
  net <- net[modDef,modDef]
  net <- igraph::graph_from_adjacency_matrix(net,mode='undirected')
  return(net)
}
mani <- dplyr::group_by(aggMods,brainRegion,ModuleNameFull)
mani2 <- dplyr::summarise(mani,ngenes=length(GeneID))
aggnets <- mapply(makeAggregateNetworks,
               mani2$brainRegion,
               mani2$ModuleNameFull,
               MoreArgs = list(aggModMani=aggMods,
                               networks=ampNetworks),SIMPLIFY = FALSE)
names(aggnets) <- mani2$ModuleNameFull


networkGraphs <- lapply(ampNetworks,
                        igraph::graph_from_adjacency_matrix,
                        mode='undirected')
#rosmap <- igraph::graph_from_adjacency_matrix(bicNetworks$network,
#                                              mode='undirected')
rm(ampNetworks)
gc()

collateGraphStatistics <- function(graph){
  model <- list()
  cat('Computing Degree...\n')
  model$degree <- igraph::degree(graph)
  #model$alpha_centrality <- igraph::alpha_centrality(graph,tol=1e-14)
  cat('Computing Authority Score...\n')
  model$authority_score <- igraph::authority_score(graph)$vector
  
  cat('Computing Closeness...\n')
  model$closeness <- igraph::closeness(graph)
  
  cat('Computing Eccentricity...\n')
  model$eccentricity <- igraph::eccentricity(graph)
  
  cat('Computing Eigenvector Centrality...\n')
  model$eigen_centrality <- igraph::eigen_centrality(graph)$vector
  
  cat('Computing Betweeness Centrality...\n')
  model$centr_betw <- igraph::betweenness(graph)
  
  cat('Computing page rank...\n')
  model$pagerank <- igraph::page.rank(graph)$vector
  
  cat('Computing Transitivity...\n')
  model$transitivity <- igraph::transitivity(graph,
                                             type='undirected')
  
  
  return(model)
}


aggNetworkProperties <- lapply(aggnets,collateGraphStatistics)
aggNetworkProperties2 <- lapply(aggNetworkProperties,
                             function(x){y <- data.frame(x[-which(names(x)=='transitivity')],stringsAsFactors=F);y$GeneID <- names(x$degree);return(y);})


networkProperties <- lapply(networkGraphs,collateGraphStatistics)
str(networkProperties)

networkProperties2 <- lapply(networkProperties,
                             function(x){y <- data.frame(x[-which(names(x)=='transitivity')],stringsAsFactors=F);y$GeneID <- names(x$degree);return(y);})

dummyFun <- function(x,y){
  merge(x,y,by='GeneID')
}

networkProperties3 <- mapply(function(x,y){
  colnames(y)[1:7] <- paste0(x,colnames(y)[1:7]);return(y)},
  names(networkProperties2),
  networkProperties2,SIMPLIFY = F)

aggNetworkProperties3 <- mapply(function(x,y){
  colnames(y)[1:7] <- paste0(x,colnames(y)[1:7]);return(y)},
  names(aggNetworkProperties2),
  aggNetworkProperties2,SIMPLIFY = F)

networkProperties4 <- Reduce(merge,networkProperties3)
networkProperties5 <- dplyr::select(networkProperties4,-GeneID)

fxn34 <- function(x,y){
  return(merge(x,y,all=TRUE))
}

aggNetworkProperties4 <- Reduce(fxn34,aggNetworkProperties3)
aggNetworkProperties5 <- dplyr::select(aggNetworkProperties4,-GeneID)



genesets1 <- synapser::synGet('syn5923958')
load(genesets1$path)
adList <- GeneSets$Alzheimers$`AD:GeneticLoci`
adList <- c(adList,'HLA-DRB5','HLA-DRB1')
adList <- adList[-which(adList=='HLA-DRB5-DRB1')]

adList2 <- list(ad_gwas=adList,
                dummyList=c('VEGF','APOE'))
#ng 2017 paper
adList2$ad_gwas <- c(adList2$ad_gwas,'ABI3','PLCG2')
#alison goate's paper
adList2$ad_gwas <- c(adList2$ad_gwas,'SPI1')
#biorxiv paper
adList2$ad_gwas <- c(adList2$ad_gwas,'ADAM10','ADAMTS4','ACE')

#solly's coloc eqtl analysis
#adList2$ad_gwas <- c(adList2$ad_gwas,'CR1','SLC4A9','C6orf48','HLA-DRB5','PSMB9','HLA-DPB1','CD2AP','SLC39A13','AREL1','RPS6KL1','AL162171.1','RIN3','USP6','ZNF594','RABEP1','NUP88','AC004148.2')
adList2$ad_gwas <- c(adList2$ad_gwas,'RIN3','AREL1','SLC4A9')


#add in loci from Jun et al.
adList2$ad_gwas <- c(adList2$ad_gwas,'HBEGF')

#add in loci from Preuss et al.
adList2$ad_gwas <- c(adList2$ad_gwas,'SNX1')

adList2$ad_gwas <- unique(adList2$ad_gwas)


#fb<-readLines('ng_gene')
#adList2$ad_gwas <- union(adList2$ad_gwas,fb)

convertAd <- utilityFunctions::convertHgncToEnsembl(adList2$ad_gwas)
adGene <- networkProperties4$GeneID%in% convertAd$ensembl_gene_id



#dv <- rowMeans(dplyr::select(networkProperties5,dplyr::ends_with('degree')))
#ev <- rowMeans(dplyr::select(networkProperties5,dplyr::ends_with('authority_score')))
#tcx <- dplyr::select(networkProperties5,dplyr::starts_with('TCX'))
#tcx <- dplyr::select(tcx,-TCXeigen_centrality)
#summary(glm(adGene ~ .,data=tcx))
fun <- function(x,y){
  return(summary(glm(y~x,family='binomial'))$coef[2,4])
}
fun2 <- function(x,y){
  return(summary(glm(y~x,family='binomial'))$coef[2,1])
}
networkProperties6 <- dplyr::select(networkProperties5,-dplyr::ends_with('eigen_centrality'))
networkProperties6 <- scale(networkProperties6)


####build deg set
degResObj <- synapser::synGet("syn10496554")
load(degResObj$path)

foo <- utilityFunctions::list2df(amp.ad.de.geneSets)
View(foo)

#get rid of missing values
foo <- dplyr::filter(foo,value!='')
View(foo)

#Add ENSEMBL Ids back
map <- utilityFunctions::convertHgncToEnsembl(foo$value)
foo <- dplyr::left_join(foo,map,by=c('value'='external_gene_name'))
foo <- foo[!duplicated(foo),]
foo2 <- tidyr::drop_na(foo)
foo2 <- dplyr::select(foo2,key,ensembl_gene_id)
foo2$value <- rep(1,nrow(foo2))
foo3 <- tidyr::spread(foo2,key,value,fill=0)

combinedFeatureSet <- dplyr::left_join(networkProperties4,foo3,by=c('GeneID'='ensembl_gene_id'))
combinedFeatureSet <- dplyr::left_join(combinedFeatureSet,aggNetworkProperties4,by=c('GeneID'))
combinedFeatureSet2 <- dplyr::select(combinedFeatureSet,-dplyr::ends_with('eigen_centrality'))
combinedFeatureSet2 <- dplyr::select(combinedFeatureSet2,-GeneID)
combinedFeatureSet2[is.na(combinedFeatureSet2)] <- 0
adGene2 <- combinedFeatureSet$GeneID %in% convertAd$ensembl_gene_id
wZ <- which(colSums(combinedFeatureSet2)<3 & (1:ncol(combinedFeatureSet2))%in%(43:107))
combinedFeatureSet2 <- combinedFeatureSet2[,-wZ]

store_on_synapse_fs <- combinedFeatureSet2

permLink <- githubr::getPermlink(repository = 'blogsdon/AMP-AD_Driver_Prediction_Engine',
                                 ref = 'branch',
                                 refName = 'dev',
                                 repositoryPath = 'ad_causal_gene_predictor_ver0_0_1.R')
fsObj<-rSynapseUtilities::pushDf2Synapse(df = store_on_synapse_fs,
                                               fileName = 'amp_ad_rna_seq_feature_set.csv',
                                               synapseFolderId = 'syn7525089',
                                               annos = c('analysisType' = 'featureLearning',
                                                         'assay' = 'RNAseq',
                                                         'consortium' = 'AMP-AD',
                                                         'dataSubType' = 'geneExp',
                                                         'dataType' = 'analysis',
                                                         'fileType' = 'csv',
                                                         'normalizationStatus' = 'TRUE',
                                                         'normalizationType' = 'CPM',
                                                         'organism' = 'HomoSapiens',
                                                         'summaryLevel' = 'gene'),
                                               comment = 'first version of network + deg features used for prediction engine',
                                               usedVector = c(foo5$id,'syn10496554'),
                                               executedVector = permLink,
                                               activityName1 = 'feature building',
                                               activityDescription1 = 'build network and deg features')





combinedFeatureSet2 <- scale(combinedFeatureSet2)

####fit models
res_cf <- apply(combinedFeatureSet2,2,fun,adGene2)
res_cf3 <- apply(combinedFeatureSet2,2,fun,adGene2)

test_lasso <- glmnet::glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial')

set.seed(1)
lasso <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial',lambda=test_lasso$lambda[1:40])
set.seed(1)
lassoA <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2[,1:42],family='binomial',lambda=test_lasso$lambda[1:40])
set.seed(1)
lassoB <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2[,43:96],family='binomial',lambda=test_lasso$lambda[1:40])
set.seed(1)
lassoC <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2[,97:348],family='binomial',lambda=test_lasso$lambda[1:40])
no_cores <- parallel::detectCores()
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
set.seed(1)
system.time(lasso <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial',lambda=test_lasso$lambda[1:40]))
set.seed(1)
system.time(lasso <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial',lambda=test_lasso$lambda[1:40],parallel = TRUE))
set.seed(1)
lassoLOOV <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial',lambda=test_lasso$lambda[1:40],nfolds = length(adGene2),parallel = TRUE)


#lasso3 <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial')
#ridge <- glmnet::cv.glmnet(y=adGene2,x=combinedFeatureSet2,family='binomial',alpha=0)
betaScore <- lasso$glmnet.fit$beta[,which(lasso$lambda==lasso$lambda.min)]
score3 <- combinedFeatureSet2[,names(betaScore)]%*%betaScore
names(score3) <- networkProperties4$GeneID
score3Df <- data.frame(gene=names(score3),score5=score3,stringsAsFactors = F)

mapTable2 <- utilityFunctions::convertEnsemblToHgnc(score3Df$gene)
score3Df <- dplyr::left_join(score3Df,mapTable2,by=c('gene'='ensembl_gene_id'))
score3Df <- dplyr::arrange(score3Df,desc(score5))
colnames(score3Df)[2] <- 'adDriverScore'
View(score3Df)


permLink <- githubr::getPermlink(repository = 'blogsdon/AMP-AD_Driver_Prediction_Engine',
                                 ref = 'branch',
                                 refName = 'dev',
                                 repositoryPath = 'ad_causal_gene_predictor_ver0_0_1.R')

edgeListObj<-rSynapseUtilities::pushDf2Synapse(df = score3Df,
                                               fileName = 'adDriverScore.csv',
                                               synapseFolderId = 'syn7525089',
                                               annos = c('analysisType' = 'supervisedLearning',
                                                         'assay' = 'RNAseq',
                                                         'consortium' = 'AMP-AD',
                                                         'dataSubType' = 'geneExp',
                                                         'dataType' = 'analysis',
                                                         'fileType' = 'csv',
                                                         'method' = 'lasso',
                                                         'normalizationStatus' = 'TRUE',
                                                         'normalizationType' = 'CPM',
                                                         'organism' = 'HomoSapiens',
                                                         'summaryLevel' = 'gene'),
                                               comment = 'add Preuss et al snps',
                                               usedVector = c(foo5$id,'syn5923958','syn10496554'),
                                               executedVector = permLink,
                                               activityName1 = 'AD gene ranking',
                                               activityDescription1 = 'build network and deg features, fit a predictive model of known LOAD genes')

synapser::synLogin()
driverScoreObj <- synapser::synGet('syn11688680')
driverScore<-read.csv(driverScoreObj$path,stringsAsFactors = F)


genecardsObj <- synapser::synGet('syn10507702')
genecards <- data.table::fread(genecardsObj$path,data.table=F)


driverScore2 <- dplyr::left_join(driverScore,genecards,by=c('external_gene_name' = 'Gene Symbol'))


adList <- list()
adList$genecards <- genecards$`Gene Symbol`



#load igap loci
igapObj<-synapser::synGet('syn10008575')
igap <- data.table::fread(igapObj$path,data.table=F)
View(igap)
igapAnnotated <- utilityFunctions::convertSnpsToGenes(igap$MarkerName)

####merge with driver score
combinedAnnos<-dplyr::left_join(igapAnnotated,driverScore,by=c('ensembl_gene_stable_id'='gene'))
dups <- duplicated(combinedAnnos[,c(1,2)])
combinedAnnos2 <- combinedAnnos[-which(dups),]

combinedAnnos2 <- dplyr::left_join(combinedAnnos2,igap,by=c('refsnp_id'='MarkerName'))
combinedAnnos3 <- combinedAnnos2[-which(duplicated(combinedAnnos2$refsnp_id)),]
combinedAnnos3 <- dplyr::arrange(combinedAnnos3,desc(adDriverScore))
combinedAnnos4 <- combinedAnnos3[1:535,]
combinedAnnos4 <- dplyr::filter(combinedAnnos4,Pvalue < 0.05/nrow(combinedAnnos4))


unique(combinedAnnos4$external_gene_name)
gap::qqunif(combinedAnnos3$Pvalue,xlim=c(0,4),ylim=c(0,42))
par(new=T)
gap::qqunif(combinedAnnos4$Pvalue,xlim=c(0,4),ylim=c(0,42),col='red')




#iga full

igapObj<-synapser::synGet('syn10008574')
igap <- data.table::fread(igapObj$path,data.table=F)
View(igap)
igapAnnotated <- utilityFunctions::convertSnpsToGenes(igap$MarkerName)

####merge with driver score
combinedAnnos<-dplyr::left_join(igapAnnotated,driverScore,by=c('ensembl_gene_stable_id'='gene'))
dups <- duplicated(combinedAnnos[,c(1,2)])
combinedAnnos2 <- combinedAnnos[-which(dups),]

combinedAnnos2 <- dplyr::left_join(combinedAnnos2,igap,by=c('refsnp_id'='MarkerName'))
combinedAnnos3 <- combinedAnnos2[-which(duplicated(combinedAnnos2$refsnp_id)),]
combinedAnnos3<-dplyr::arrange(combinedAnnos3,desc(adDriverScore))

combinedAnnos4 <- combinedAnnos3[1:535,]
combinedAnnos4 <- dplyr::filter(combinedAnnos4,Pvalue < 0.05/nrow(combinedAnnos4))

unique(combinedAnnos4$external_gene_name)



gap::qqunif(combinedAnnos3$Pvalue,xlim=c(0,4),ylim=c(0,42))
par(new=T)
gap::qqunif(combinedAnnos4$Pvalue,xlim=c(0,4),ylim=c(0,42),col='red')
