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
infection_rates = function (countries) {

	#variables from the UI: 
	countries = c("Italy", "Canada", "China", "US","Switzerland")
	
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

	#plot with ggplot
	cv1_cr150 %>%
		ggplot(aes( x=day_from, y = N, color= Country) ) +
		geom_line()+
		geom_point()+
		xlab("Days since case number 150")+
	  	ylab("Reported cases")

	#all of the countries in the list: 
	all_countries = unique(cv1_ts$Country)

}

#==============================================================================
# Shiny server logic
#==============================================================================
# Tell the server how to run this: 

shinyServer( function(input, output) {
	
	#Basic plot of the infection rates
	output$infection_rates = renderPlot({
		infection_rates(input$countries)
		})

	#Download the possible countries to plot: 
	datasetInput = reactive( infection_rates (countries))

	output$downloadData =	 downloadHandler(
    	filename= "country_list.csv",
    	
    	content = function(file) {
      	write.csv(datasetInput(), file, row.names = FALSE)
		}
	)



	url1 = a("R code", href="https://github.com/jusinowicz/misc_populations/tree/master/covid_shiny")
	url2 = a("me", href="http://jacobusinowicz.com/")

    output$tag1 <- renderUI({
      tagList("Go to Jacob Usinowicz's github page for:", url1)
	})
	output$tag2 <- renderUI({
      tagList("More about", url2)
	})
})