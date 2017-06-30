library(shiny)

GSE_names <- readRDS('/Users/sumitmukherjee/Documents/SageStuff/ConservedGenesets.rds')
GSE_names <- GSE_names$NamesInt


# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("Gene set enrichment visualization"),
  
  # Sidebar with controls to select the variable to plot against mpg
  # and to specify whether outliers should be included
  sidebarPanel(
    selectInput("dataset", "Choose a dataset:",
                choices = GSE_names) 
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    h3(textOutput("caption")),
    
    plotOutput("mpgPlot"), 
    
    tabsetPanel(
      tabPanel('Legend', 
               dataTableOutput('myTable'))
    )
    
  )
))