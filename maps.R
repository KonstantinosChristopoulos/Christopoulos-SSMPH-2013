#Figure 1
library(usmap)
library(ggplot2)
#import data
d<-read.csv(file=".../maps.csv",
colClasses=c("fips"="character"))
d$siP<-d$si/100

# empirical suicide rates
d1<-data.frame(fips=d$fips, sui=d$sui_rate)
plot_usmap(
  regions = "states",  data = d1, values = "sui", color= "black") +
   scale_fill_continuous(name = "Suicide mortality per 100k", low = "white", high = "red", label=scales::comma)+ theme(legend.position = "right")
 
#Stringency Index 
d2<-data.frame(fips=d$fips, si=d$siP)
plot_usmap(
  regions = "states",  data = d2, values = "si", color= "black") +
   scale_fill_continuous(name = "Stringency Index %", low = "white", high = "orange", label=scales::comma)+ theme(legend.position = "right")
