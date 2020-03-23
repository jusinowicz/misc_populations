library(shiny)
#==============================================================================
# This is R code to download and plot the latest Covid-19 trajectories. 
#==============================================================================
#Load libraries
library(tidyverse)
library(lubridate)
#==============================================================================

#==============================================================================
# Plot the infection rates. 
#==============================================================================
infection_rates = function (countries,get_countries=F) {

	#variables from the UI: 
	#countries = c("Italy", "Canada", "China", "US","Switzerland")
	
	#Download the latest data set
	cv1=read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
	colnames(cv1)[1] = "State"  #Rename the headings
	colnames(cv1)[2] = "Country"

	d2 = dim(cv1)[2] #Get the length of the max time series
	#Reshape the data set into a time series:
	cv1_ts = cv1 %>%
		gather(key = day, value = N, colnames(cv1)[5:d2]) 
	cv1_ts$day = mdy(cv1_ts$day)

	#Pull out the countries specified in countries
	cv1_use = cv1_ts %>% 
		filter(Country %in%  countries)

	#Sum across states in a country for totals
	cv1_cr = cv1_use %>%
		group_by(Country,day )%>% 
	  	summarise( N = sum(N))

	#Find the day at which reported cases reaches 150
	cv1_cr150 = subset(cv1_cr, N > 150 & N< 60000) %>%
		mutate(day_from = as.integer(day))

	#Rescale the time series so that they all start on day_150	
	for (n in 1:length(countries)){
		cv1_cr150$day_from[cv1_cr150$Country == countries[n]]  = 
			subset(cv1_cr150, Country == countries[n])$ day_from - 
				min(subset(cv1_cr150, Country == countries[n])$ day_from)		
	}


	#all of the countries in the list: 
	if(get_countries == T ) { 
		return(	all_countries = unique(cv1_ts$Country) )
	}

	#cv1_cr150

	ggplot(cv1_cr150, aes( x=day_from, y = N, color= Country) ) +
		geom_line()+
		geom_point()+
		xlab("Days since case number 150")+
	  	ylab("Reported cases")


}

#==============================================================================
# Shiny server logic
#==============================================================================
# Tell the server how to run this: 

shinyServer( function(input, output,session) {
	
	#Initialize countries:
	countries = c("Italy", "Canada", "China", "US","Switzerland")

	#Get the full country list: 
	country_list =  infection_rates (countries,get_countries=T)
	#Update the full list with this command. This makes plotting interactive, 
	#because it will pull the names of countries from the actual data set.
	updateSelectizeInput(session, 'countries',
              choices = country_list,
              selected = c("Italy", "Canada", "China", "US","Switzerland"),
              server = TRUE
  	)

	#The main plot. This takes the updated user-input country choices and 
	#passes them to the main function above. 
	output$mainplot = renderPlot({
		countries = input$countries
		 infection_rates(input$countries,get_countries=F)
  	})

	#For a complete list of all of the countries in the data set: 
	datasetInput = reactive( infection_rates (input$countries,get_countries=T))

	output$downloadData =	 downloadHandler(
    	filename= "country_list.csv",
    	
    	content = function(file) {
      	write.csv(datasetInput(), file, row.names = FALSE)
		}
	)

	#Just tag things with my own info
	url3 = a("CSSE at Johns Hopkins University", href = "https://systems.jhu.edu/" )
	url1 = a("R code", href="https://github.com/jusinowicz/misc_populations/tree/master/covid_shiny")
	url2 = a("me", href="http://jacobusinowicz.com/")


    output$tag1 <- renderUI({
      tagList("Date come from", url3)
	})
    output$tag2 <- renderUI({
      tagList("Go to Jacob Usinowicz's github page for:", url1)
	})
	output$tag3 <- renderUI({
      tagList("More about", url2)
	})

})