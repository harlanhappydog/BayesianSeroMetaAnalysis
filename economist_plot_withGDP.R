# Code to replicate approximately figure from the Economist

theta_md <- -5.23
theta1_md <- 0.48
theta2_md <- -0.20


predictIFR <- function(age65, GDPppp){

Z1_star <- (log(age65)-mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older))
Z2_star <- (log(GDPppp)-mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp))

100*( 1 - exp(-exp(theta_md + theta1_md*(Z1_star) 
				+  theta2_md*(Z2_star)  )) )	
}


library(countrycode)
library(WDI)
wdi <- WDI(country = 'all',
           indicator = c('wdi_gdppc_ppp' = 'NY.GDP.PCAP.PP.CD',
                         'wdi_pop_over_65' = 'SP.POP.65UP.TO.ZS'))
                     
wdi$iso3c <- countrycode(wdi$iso2c, "iso2c", "iso3c")
wdi_all <- wdi 
wdi$iso2c <- NULL
wdi$country <- NULL

# Only latest observation
wdi <- wdi[order(wdi$year), ]
wdi <- wdi[!is.na(wdi$iso3c), ]  


for(i in setdiff(colnames(wdi), c("year", "iso3c"))){
  wdi[, i] <- ave(wdi[, i], wdi$iso3c, 
                  FUN = function(x){
                    if(max(which(!is.na(x))) == -Inf){
                      NA
                    } else {
                      x[max(which(!is.na(x)))]
                    }
                  })
}

# Collapse rows to one per country (multiple happens when data comes from different years)
for(i in setdiff(colnames(wdi), "iso3c")){
  wdi[, i] <- ave(wdi[, i], wdi$iso3c, FUN = function(x) mean(x, na.rm = T))
}
wdi <- unique(wdi)
wdi_1<-unique(merge(wdi, wdi_all[,c("iso3c","country")], by = "iso3c", all.x=TRUE, all.y=FALSE))

cov_data <-wdi_1
cov_data[,"Country.Name"]<-cov_data[,"country"]

cov_data[,"ifr"] <- apply(cov_data[,c("wdi_pop_over_65", "wdi_gdppc_ppp")] , 1, function(z) predictIFR(z[1],z[2]))


cov_data[, "selectnames"]<-NA

cov_data[cov_data[,"Country.Name"]=="Singapore", "selectnames"]<-"Singapore"
cov_data[cov_data[,"Country.Name"]=="Qatar", "selectnames"]<-"Qatar"
cov_data[cov_data[,"Country.Name"]=="Saudi Arabia", "selectnames"]<-"Saudi Arabia"
cov_data[cov_data[,"Country.Name"]=="Nigeria", "selectnames"]<-"Nigeria"
cov_data[cov_data[,"Country.Name"]=="India", "selectnames"]<-"India"
cov_data[cov_data[,"Country.Name"]=="Peru", "selectnames"]<-"Peru"
cov_data[cov_data[,"Country.Name"]=="South Africa", "selectnames"]<-"S. Africa"
cov_data[cov_data[,"Country.Name"]=="United States", "selectnames"]<-"US"
cov_data[cov_data[,"Country.Name"]=="United Kingdom", "selectnames"]<-"UK"
cov_data[cov_data[,"Country.Name"]=="Italy", "selectnames"]<-"Italy"
cov_data[cov_data[,"Country.Name"]=="Japan", "selectnames"]<-"Japan"
cov_data[cov_data[,"Country.Name"]=="Spain", "selectnames"]<-"Spain" 
 

  annotations <- data.frame(
        xpos = c(-Inf,1200,Inf,Inf),
        ypos =  c(-Inf, 1.5,-Inf,Inf),
        annotateText = c("","Size = Share of \n population over 65"
                        ,"","")) #<- adjust

cov_data[, "aborder"] <-  !is.na(cov_data[, "selectnames"])
ggplot(cov_data, aes(x = wdi_gdppc_ppp, y = ifr, label = selectnames ))  +  geom_point(aes(size =(wdi_pop_over_65)), alpha = 0.5, col = "cornflowerblue") +   geom_point(aes(size =(wdi_pop_over_65), col = border), fill="blue", alpha = 0.5) + 
  scale_color_manual(values = c( "#E7B800","#00AFBB", "#FC4E07"))+
  scale_size(range = c(0.25, 10)) +  scale_x_continuous(trans = "log10", breaks=c(1000,5000,10000,50000,100000), labels=c("1,000", "5,000", "10,000", "50,000", "100,000")) + xlab("GDP per capita at purchasing-power parity ($), log scale") +    ylab("") +   geom_text(size = 4, colour = "black", vjust = -1) +  ylim(0,1.5) + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_line(colour = "grey"),
    panel.border = element_blank(),
    panel.background = element_blank(), plot.title = element_text(hjust = 0.0, vjust=2.12), legend.position="none")  + 
        geom_text(data = annotations, aes(x=xpos,y=ypos,label=annotateText))+ggtitle("Expected IFR for a population with a similar age and GDP (%)") 
  