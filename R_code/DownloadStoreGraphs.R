synapser::synLogin(email = 'mukhes3@uw.edu', password = 'Nirvana99')
foo <- synapser::synTableQuery("SELECT * FROM syn8681664 WHERE ( ( method = 'bic' ) AND ( study = 'MSBB' OR study = 'MayoRNAseq' OR study = 'ROSMAP' ) )")
foo <- foo$asDataFrame()
ext <- '.edgelist'

library(synapser)
library(igraph)
for (i in 1:7){
  bla <- synGet(foo$id[i])
  bla <- load(bla$path)
  
  tmp <- paste(foo$tissueTypeAbrv[i],ext,sep = '')
  graph = graph.adjacency(bicNetworks$network, mode = "undirected")
  write_graph(graph,tmp,format="edgelist")
  #writeMM(bicNetworks$network, tmp)
  #Genes = colnames(bicNetworks$network)
  #Genes = as.data.frame(Genes)
  #write.csv(Genes, tmp)

}