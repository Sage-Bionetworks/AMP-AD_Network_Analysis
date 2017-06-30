library(shiny)

#library(datasets)
library(ggplot2)

source('/Users/sumitmukherjee/Documents/SageStuff/SumitUtilFcns.R')
source('/Users/sumitmukherjee/Documents/SageStuff/CreateMagmaFiles.R')
 
ModType <- 'consensus'
ModNameEnr <- 'consensus37DLPFC'


#head(TempEnr)

# Define server logic required to plot various variables against mpg
shinyServer(function(input, output) {
  
  # Compute the forumla text in a reactive expression since it is 
  # shared by the output$caption and output$mpgPlot expressions
  formulaText <- reactive(input$dataset)
  formulaText2 <- reactive(input$ModName)
  
  # Return the formula text for printing as a caption
  output$caption <- renderText({
    formulaText()
  })
  
  tmp <- reactive({
    if (is.null(formulaText())|is.null(formulaText2())){
      return(NULL)
    }
      
    TempEnr <- TempEnr[[formulaText()]]
    ModNameEnr <- formulaText2()
  In <- which(TempEnr$ModuleNameFull %in% ModNameEnr)
  tmp2 <- list()
  Category<- TempEnr$category[In]
  fisherOR <- TempEnr$fisherOR[In]
  Log10_pval <- -log10(TempEnr$fisherPval[In])
  I = sort(Log10_pval, decreasing = T, index.return = T)
  I <- I$ix
  tmp2$Category<- Category[I]
  tmp2$fisherOR <- fisherOR[I]
  tmp2$Log10_pval<- Log10_pval[I]
  tmp2 <- data.frame(tmp2, stringsAsFactors = F)
  tmp2
  #Log10_pval <- Log10_pval[I]
  #Log10_pval
  })
  
output$myTable <- renderDataTable({
  if(is.null(tmp())){
    return(NULL)
  } else {
    t <- tmp()
  }
  t
})
    
#  Generate a plot of the requested variable against mpg and only
#  include outliers if requested
  # output$mpgPlot <- renderPlot({
  #   if(!is.null(tmp())){
  #     dat <- tmp()
  #   p <- ggplot(data = dat, aes(x = Category, y = Log10_pval)) +
  #   geom_bar(stat = 'identity') +
  #   theme(axis.text.x = element_text(angle = 90, size = 5)) +
  #   scale_x_discrete(limits = dat$Category)
  #   p}
  # })
  output$mpgPlot <- renderPlot({
    if(!is.null(tmp())){
      t <- tmp()
      barplot(t$Log10_pval, xlab = 'genesets',
              ylab = 'Minus_log10(Pval)')
    }
  })
})