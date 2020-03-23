#==============================================================================
# The UI to the corona.R code, to download and plot the latest 
# Covid-19 trajectories. 
#==============================================================================
library(shiny)
library(shinyIncubator)
library(rhandsontable)
#Load libraries
library(tidyverse)
library(lubridate)
#==============================================================================


# Define UI for function that draws a population growth trajectory
shinyUI(fluidPage(
  # Application title
  titlePanel("Covid-19 infection rates"),
  
  # Sidebar with input for the countries to plot: 
  sidebarLayout(
    sidebarPanel(

   	selectizeInput("countries", "Countries", choices=NULL, multiple =T, 
   		selected = c("Italy", "Canada", "China", "US","Switzerland"),
   		options = list( placeholder = 'Type a name, e.g. US')),

      # Button
      downloadButton("downloadData", "Countries in data set")
     
    ),
    
    # Plot the infection rates
    mainPanel(
      plotOutput("mainplot")
    )
  ),

  uiOutput("tag1"),
  uiOutput("tag2")

))