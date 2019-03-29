#Creating tool for drug gene enrichment analysis 
dg <- readRDS('drug_target_associations_v2.rds')

#Reading gene names and keeping only those that are predicted to be drivers 
res <- readRDS('AD_UnionModelResults.rds')

#Getting the genes that are predicted drivers 
GenePred <- as.vector(unique(res$GeneSymb[res$Y_union==1]))

#Keeping only drug interactions with specificity above a threshold
th <- 0.1 
dg <- dg[dg$known_selectivity_index>=th,]

#Keeping only interactions with predicted driver genes 
dg <- dg[which(dg$hugo_gene %in% GenePred),]

#Keep only drugs which interact with more than one predicted driver gene 
dg <- dg[ dg$common_name %in%  names(table(dg$common_name))[table(dg$common_name) >1] , ]

#Keep only genes which interact with one than one drug
dg <- dg[ dg$hugo_gene %in%  names(table(dg$hugo_gene))[table(dg$hugo_gene) >1] , ]

#creating matrix of undirected weights 
library(Matrix)
g <- unique(dg$hugo_gene)
c <- unique(dg$common_name)

M <- matrix(rep(0,length(g)*length(c)),nrow = length(g), ncol = length(c))

for (i in 1:length(g)){
  In <- which(c %in% dg$common_name[dg$hugo_gene==g[i]])
  M[i,In] <- 1
}

rownames(M) <- g
colnames(M) <- 1:length(c)

#Performing hierarchical clustering on drug-gene matrix 
library(pheatmap)
pheatmap(M, cluster_cols = T, clustering_distance_rows = 'binary',
         clustering_distance_cols = 'binary')

