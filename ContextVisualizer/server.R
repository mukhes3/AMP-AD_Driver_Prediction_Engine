library(shiny)
library(igraph)
library(ggplot2)

#library(datasets)
#setwd('/Users/sumitmukherjee/Documents/
#               DriverPrediction/MiscFiles/')
#source('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/VisualizerHelperFiles.R')
Dat <- readRDS('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/AD_UnionModelResults.rds')

#functions
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
X0 = Process.Feature(read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/amp_ad_de_feature_set.csv'))
X1 = Process.Feature(read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles//amp_ad_agg_feature_set.csv'))
X2_s = read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/amp_ad_deVal_feature_set2.csv')
X2 = Process.Feature(read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/amp_ad_deVal_feature_set2.csv'))
X3 = Process.Feature(read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/amp_ad_subNet_feature_set.csv'))

NetDat <- read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/AggregateMods_March2018.csv',stringsAsFactors = F)

IGAP_list <- read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/IGAP_gene_summary.csv')


#Read original response vector 
Y = read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/ResponseVec_040318.csv',stringsAsFactors = F)
GeneId <- Y$GeneID
Y <- Y$Y

#Get predicted response vector
Y_pred_df <- read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/EC2_ConsUnion1_NoDeProbs_0814.csv')
Y_pred <- (Y_pred_df$Y1a + Y_pred_df$Y1b + 
             Y_pred_df$Y3a + Y_pred_df$Y3b)/4.0

#print('bla1')
#head(TempEnr)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  
  Feat <- reactive(input$FeatureList)
  Gene <- reactive(input$Genes)
  BR <- reactive(input$BR)
  FV <- reactive(input$FeatureViz)
  
  tmp2 <- reactive({
    t <- list()
    t$l <- c(1:5)
    #print('bla3')
    
    t
  })
  
  tmp <- reactive({
    
    
    
    if(is.null(Feat())){
      return(NULL)
    }
    

    if(is.null(Gene())){
      return(NULL)
    }
    
    if(is.null(BR())){
      return(NULL)
    }
    
    
    
    Feat <- Feat()
    Gene <- Gene()
    #GeneName <- Dat$GenesENSG[which(Dat$GeneSymb==Gene)][1]
    BR <- BR()
    
    #Get neighbors of other genes in coexpression module 
    In <- which((NetDat$brainRegion==BR) & (NetDat$external_gene_name==Gene))
    In1 <- which((NetDat$ModuleName == NetDat$ModuleName[In]) &
                   (NetDat$external_gene_name != Gene))
    In2 <- which((GeneId %in% NetDat$GeneID[In1]) & (Y_pred>0.5))
    
    if(Feat == 'DE'){
      PC1 <- Dat$PC_de1
      PC2 <- Dat$PC_de2
      X <- X0
    } else if(Feat == 'Module-Network'){
      PC1 <- Dat$PC_mod1
      PC2 <- Dat$PC_mod2
      X <- X1
    } else if(Feat == 'DE-Values'){
      PC1 <- Dat$PC_deVal1
      PC2 <- Dat$PC_deVal2
      X <- X2
      
    } else {
      PC1 <- Dat$PC_global1
      PC2 <- Dat$PC_global2
      X <- X3
    }
    
    
    #visualize a given gene on the PCA plot
    In_gene <- which(Dat$GeneSymb==Gene)
    Y_pc = (Y_pred>0.5) + 0.0
    Y_pc[In_gene] <- 2.0
    Y_p = (Y_pred>0.5) + 0.0
    Y_p = Y_p*.5
    Y_p[In_gene] <- 2.0
    
    #getting butterfly plot parameters 
    LogFC <- X2_s[[paste('LogFC',BR,sep='.')]]
    LogPV <- X2_s[[paste('LogPV',BR,sep='.')]]
    
    #creating appropriate graph object 
    par(mfrow = c(1,1))
    g_obj <- make_star(length(In2)+1, mode = 'undirected')
    #V_list <- c(GeneName, GeneId[In2])
    V_list <- c(Gene, as.vector(Dat$GeneSymb[In2]))
    #V_list <- Get.Gene.Symbs(V_list,NetDat)
    V_size <- -log10(Get.Pval(V_list,IGAP_list))
    V_size[V_size>20] <- 20
    palf <- colorRampPalette(c(rgb(1,1,1,.2),rgb(0.8,0,0,.7)),
                             alpha=TRUE)
    V(g_obj)$color <- palf(V_size)
    V(g_obj)$label <- V_list
    
    print('bla2')
    
    In_Y <- which(Dat$Y_union > 0)
    S <- log10(Dat$Mod[In_Y]/(1-Dat$Mod[In_Y]))
         + log10(Dat$Global[In_Y]/(1-Dat$Global[In_Y]))  - Dat$Log10Mean[In_Y]
    
    Sx <- sort(S, decreasing = T)
    
    print(Sx)

    t <- list()
    t$PC1 <- PC1
    t$PC2 <- PC2 
    t$Y_p <- Y_p 
    t$Y_pc <- Y_pc 
    t$LogFC <- LogFC
    t$LogPV <- LogPV 
    t$g_obj <- g_obj
    t$BR <- BR
    t$Feat <- Feat
    t$Gene <- Gene
    t$X <- X
    t$In_gene <- In_gene
    t$Score <- log10(Dat$Mod[In_gene]/(1-Dat$Mod[In_gene]))
    + log10(Dat$Global[In_gene]/(1-Dat$Global[In_gene]))- Dat$Log10Mean[In_gene]
    t$Rank <- which(Sx == t$Score)

    print(In_gene %in% In_Y)
    
    t
    

  })
  
  
  output$VizFeature <- renderUI({
    
    if(is.null(tmp())){
      return(NULL)
    }
    
    t <- tmp()
    
    FeatList <- colnames(t$X)
    selectInput("FeatureViz", "Choose Feature", FeatList)
  })
  
  output$PCA_plot <- renderPlot({

    
    if(is.null(tmp())){
      return(NULL)
    }
    
    t <- tmp()
    
    
    plot(t$PC1,t$PC2, col = t$Y_pc+1, pch=16, cex = t$Y_p + .1, main = t$Feat)
    
    
  })
  
  output$Butterfly_plot <- renderPlot({
    if(is.null(tmp())){
      return(NULL)
    }

    t <- tmp()

    plot(t$LogFC,t$LogPV,col = t$Y_pc+1, pch=16, cex = t$Y_p + .1,
         xlim = c(-.5,.5),main = t$BR)


  })


  output$NetworkPlot <- renderPlot({

    if(is.null(tmp())){
      return(NULL)
    }

    t <- tmp()

    plot(t$g_obj, vertex.label.cex = 0.5, main = t$BR)

  })
  
  output$FeatureView1 <- renderPlot({
    
    if(is.null(tmp())){
      return(NULL)
    }
    
    if(is.null(FV())){
      return(NULL)
    }
    
    t <- tmp()
    FV <- FV()
    
    temp <- t$X[[FV]]
    #print(temp)
    p1 <- hist(temp, main = paste(c(t$Gene,'=',temp[t$In_gene]),sep=' '))


  })
  
  output$FeatureView2 <- renderPlot({
    
    if(is.null(tmp())){
      return(NULL)
    }
    
    if(is.null(FV())){
      return(NULL)
    }
    
    t <- tmp()
    FV <- FV()
    
    temp <- t$X[[FV]]
    #print(temp)
    
    pal = colorRampPalette(c("blue", "red"))
    
    SO = findInterval(temp, sort(temp))
    
    plot(t$PC1,t$PC2, col=pal(length(temp))[SO], pch=16, cex = t$Y_p + .1, main = t$Feat)

    legend("topright", col=pal(2), pch=19,
           legend=range(temp))
  })

  output$RT1 <- renderText({
    
    if(is.null(tmp())){
      return(NULL)
    }
    
    t <- tmp()
    
    txt <- paste(c('Score = ',t$Score), sep = '')
    print(txt)
    
  })  
  
  output$RT2 <- renderText({
    
    if(is.null(tmp())){
      return(NULL)
    }
    
    t <- tmp()
    
    txt <- paste(c('Rank = ',t$Rank), sep = '')
    print(txt)
    
  })  



})