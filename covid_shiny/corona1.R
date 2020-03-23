library(tidyverse)
library(lubridate)
countries = c("Italy", "Canada", "China", "US","Switzerland")
#cv1=read_csv("time_series_19-covid-Confirmed.csv")
cv1=read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
colnames(cv1)[1] = "PS"
colnames(cv1)[2] = "CR"
#cv1=as.data.frame(cv1)
# colnames(cv1)[5:61]=mdy(colnames(cv1)[5:61])
# rename_at(cv1, cv1[,5:61], mdy(colnames(cv1)[5:61])
d2 = dim(cv1)[2]
cv1_ts = cv1 %>%
	gather(key = day, value = N, colnames(cv1)[5:d2]) 
cv1_ts$day = mdy(cv1_ts$day)

cv1_use = cv1_ts %>% 
	filter(CR %in%  countries)

cv1_cr = cv1_use %>%
	group_by(CR,day )%>% 
  	summarise( N = sum(N))

cv1_cr150 = subset(cv1_cr, N > 150 & N< 60000) %>%
	mutate(day_from = as.integer(day))
for (n in 1:length(countries)){
	cv1_cr150$day_from[cv1_cr150$CR == countries[n]]  = 
		subset(cv1_cr150, CR == countries[n])$ day_from - 
			min(subset(cv1_cr150, CR == countries[n])$ day_from)		
}

cv1_cr150 %>%
	ggplot(aes( x=day_from, y = N, color= CR) ) +
	geom_line()+
	geom_point()