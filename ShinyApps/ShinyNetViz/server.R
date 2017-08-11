library(shiny)

#library(datasets)
library(ggplot2)

source('/Users/sumitmukherjee/Documents/SageStuff/SumitUtilFcns.R')
source('/Users/sumitmukherjee/Documents/SageStuff/CreateMagmaFiles.R')
source('/Users/sumitmukherjee/Documents/SageStuff/GenGraphViz.R')

MAGMA_pval <- 
  readRDS('/Users/sumitmukherjee/Documents/SageStuff/GWAS_files/DLPFC_results_magma/Magma_AllMethods.Rda')


#head(TempEnr)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  formulaText <- reactive(input$dataset)

  # Return the formula text for printing as a caption
  output$caption <- renderText({
    formulaText()
  })
  
  tmp <- reactive({
    if (is.null(formulaText())){
      return(NULL)
    }
    GeneSetName <- formulaText()
    for (i in 1:length(Types)){
      tmp <- data.frame(EnrList[[Types[i]]][[GeneSetName]])
      
      if (i == 1){
        ResFrame <- tmp
      } else {
        ResFrame <- rbind(ResFrame, tmp)
      }
    }
    EnrByCat <- AssignEnrichedModule(ResFrame)
    
    
    #9. For each module in the graph, obtain the MAGMA p-value
    Mod_pval <- rep(1,length(names(V(Net))))
    
    In1 <- which(names(V(Net)) %in% MAGMA_pval$SET)
    Mod_pval[In1] <- MAGMA_pval$P
    
    LegendNames <- unique(EnrByCat$CatList)
    Labels <- rep(1,length(Mod_pval))
    for (i in 1:length(Mod_pval)){
      Labels[i] <- which( LegendNames %in% EnrByCat$CatList[i])
    }
    
    pval_cut <- 0.01 
    
    tmp2 <- list()
    tmp2$Labels <- Labels
    tmp2$Mod_pval <- Mod_pval
    tmp2$LegendNames <- LegendNames
    tmp2$pval_Inv <- which(Mod_pval > pval_cut)
    tmp2$Labels[tmp2$pval_Inv] <- max(Labels) + 1
    tmp2$Mod_pval[tmp2$pval_Inv] <- 1
    tmp2$LegendNames[max(Labels)+1] <- 'irrelevant'
    
    tmp2
    })
  
  
  output$myTable <- renderDataTable({
    if(is.null(tmp())){
      return(NULL)
    } else {
      t <- tmp()
      tab <- list()
      tab$LegendNo <- c(1:length(t$LegendNames))
      tab$LegendNames <- t$LegendNames
      tab <- data.frame(tab, stringsAsFactors = F)
    }
    tab
  })
  
  output$mpgPlot <- renderPlot({
    if(!is.null(tmp())){
      t <- tmp()
      Names <- names(V(Net))
      Names[t$pval_Inv] <- ''
      
      PlotGraphGivenLabels(Net, l, t$Labels, -log10(t$Mod_pval),
                           sizeMin = .5, VertName = Names)
    }
  })
})