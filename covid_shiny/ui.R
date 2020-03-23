#==============================================================================
# The UI to the corona.R code, to download and plot the latest 
# Covid-19 trajectories. 
#==============================================================================
library(shiny)
library(shinyIncubator)
library(rhandsontable)
#==============================================================================


# Define UI for function that draws a population growth trajectory
shinyUI(fluidPage(
  # Application title
  titlePanel("Covid-19 infection rates"),
  
  # Sidebar with input for the countries to plot: 
  sidebarLayout(
    sidebarPanel(

   	selectizeInput(
        "countries"
        , "Enter countries to plot (Default: China, US, Italy, Switzerland, Canada).
        	Include a space between entries"
        , choices = NULL
        , multiple = TRUE
        , options = list(create = TRUE)
      )

      # Button
      downloadButton("downloadData", "Download"),
     
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("infection_ratest")
      )
    )
  ),

  uiOutput("tag1"),
  uiOutput("tag2")

))