library(shiny)


Dat <- readRDS('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/AD_UnionModelResults.rds')
GeneList <- unique(Dat$GeneSymb[Dat$Y_union>0])

Features <- c('DE','Module-Network','DE-Values','Global-Network')

NetDat <- read.csv('/Users/sumitmukherjee/Documents/DriverPrediction/MiscFiles/AggregateMods_March2018.csv',stringsAsFactors = F)

BR <- unique(NetDat$brainRegion)


#selectInput('featType', 'FeatureType', Features)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel('Data and feature selection'),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput('FeatureList', 'Features:',
                choices = Features), 
    selectInput('Genes', 'Genes:',
                choices = GeneList),
    selectInput('BR', 'Brain Region:',
                choices = BR), 
    uiOutput("VizFeature"),
    textOutput('RT1'),
    textOutput('RT2')
    
    
  ),
  
  # Let the user pick a particular feature from the feature list 
  mainPanel(
    #plotOutput('PCA_plot'),
    #plotOutput('Butterfly_plot'),
    #plotOutput('NetworkPlot')
    
    tabsetPanel(type = "tabs",
                tabPanel("PCA Plot", plotOutput('PCA_plot')),
                tabPanel("DE Butterfly plot", plotOutput('Butterfly_plot')),
                tabPanel("Co-expression module", plotOutput('NetworkPlot')),
                tabPanel("FeatureView", plotOutput('FeatureView1'),
                         plotOutput('FeatureView2'))
    )

  )
))