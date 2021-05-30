##########################################
####  load required libraries and data:
library("rjags")
library("hBayesDM")
library("RCurl")
library("coda")
library("xtable")

csvfile <- getURL("https://raw.githubusercontent.com/harlanhappydog/BayesianSeroMetaAnalysis/main/IFRdata_withGDP_etc.csv")
IFRdata <- read.csv(text = csvfile)

MCMCn <- 2000000

##########################################
####  functions:
summary_md_ci <- function(xx, samps){
c(  md=summary(samps)$quantiles[xx, "50%"], 
	lower=HDIofMCMC(unlist(((samps[[1]])[,xx])), credMass = 0.95)[1], 
	higher=HDIofMCMC(unlist(((samps[[1]])[,xx])), credMass = 0.95)[2]
	)}

	
##########################################
####  model in JAGS:	
metaIFR <- "model {

# Priors:

	icloglog_theta0 ~ dbeta(0.3, 30); 
	icloglog_beta  ~ dbeta(1, 3);
	theta0 <- log(-log(1-icloglog_theta0));
	beta  <- log(-log(1-icloglog_beta));

	inv.var_sig   <- (1/sigma)^2 ;
	inv.var_tau   <- (1/tau)^2 ;
	sigma     ~ dnorm(0, 1/10) T(0,);
	tau     ~ dnorm(0, 1/10) T(0,);
	theta1 ~ dnorm(0, 1/10);
	theta2 ~ dnorm(0, 1/10);

# Likelihood:
	
for(k in 1:K){
	cc[k] ~ dbin(ir[k], tests[k]);
    censor.index[k] ~ dinterval(deaths[k], c(deaths_lower[k], deaths_upper[k]))
	deaths[k] ~ dbin(ifr[k]*ir[k], pop[k]);

	cloglog(ir[k]) <- cloglog_ir[k];
	cloglog(ifr[k]) <- cloglog_ifr[k];
	
	cloglog_ir[k] ~ dnorm(beta, inv.var_sig);
	cloglog_ifr[k] ~ dnorm(theta0 + theta1*Z1[k] + theta2*Z2[k], inv.var_tau);
	
  }
  
# Summary
g_IFR_9p = theta0 + theta1*(age65_values[1]) + theta2*(gdp_values[1]);  
g_IFR_16p = theta0 + theta1*(age65_values[2])+ theta2*(gdp_values[2]);    
g_IFR_20p = theta0 + theta1*(age65_values[3])+ theta2*(gdp_values[3]);  

IFR_9p <- 1 - exp(-exp(g_IFR_9p))  
IFR_16p <- 1 - exp(-exp(g_IFR_16p))  
IFR_20p <- 1 - exp(-exp(g_IFR_20p))  

epsilon ~ dnorm(0,1);  
predictIFR_9p <- 1 - exp(-exp(g_IFR_9p + tau*epsilon))  
predictIFR_16p <- 1 - exp(-exp(g_IFR_16p + tau*epsilon))  
predictIFR_20p <- 1 - exp(-exp(g_IFR_20p + tau*epsilon))  

}
"


##########################################
### Fit model ###
IFRdata <- IFRdata[!IFRdata[,"Author"]%in%c("Hallal et al."),]

IFRdata <- IFRdata[!IFRdata[,"Author"]%in%c("Appa et al.", "Bendavid et al.", "McLaughlin et al.", "Rosenberg et al."),]

print(xtable(IFRdata[,c("Author", "Location", "window", "IRinterval")], digits =c(rep(0,3), 2,2)), include.rownames=FALSE)

print(xtable(IFRdata[,c("Location", "Population",  "deaths14_lower" ,  "deaths14_upper", "total_tests", "total_cases", "aged_65_older", "GDPppp_st")], digits = 0), include.rownames=FALSE)
	
IFRdata[IFRdata[,"Location"]=="Delhi, India","Population"]<- 19800000	
K <- length(IFRdata$total_tests)
K
	
# GDP PPP per capita 
# World = 17,811.3;     
# US = 65,297.5
# EU = 47,827.7	

gdp_values <- c(
(log(17811.3) - mean(log(IFRdata$GDPppp_st)))/sd(log(IFRdata$GDPppp_st)),
(log(65297.5) - mean(log(IFRdata$GDPppp_st)))/sd(log(IFRdata$GDPppp_st)),
(log(47827.7) - mean(log(IFRdata$GDPppp_st)))/sd(log(IFRdata$GDPppp_st)))

# Age over 65
# World = 9    
# US = 16   
# EU = 20	

age65_values <-  c(
(log(9) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)),
(log(16) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)),
(log(20) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)))


jags.modelIFR <- jags.model(textConnection(metaIFR), 
	data = list(
    K = K,
    tests = IFRdata$total_tests,
    cc = IFRdata$total_cases,    
    pop = IFRdata$Population,
    deaths_lower = IFRdata$deaths14_lower-1,
    deaths_upper = IFRdata$deaths14_upper,
    deaths = rep(NA, K),
    Z1 = c(scale(log(IFRdata$aged_65_older))),
    Z2 = c(scale(log(IFRdata$GDPppp_st))),
    censor.index = rep(1, K),
    age65_values = age65_values,
    gdp_values = gdp_values
    ), 
    n.chains = 5, 
    n.adapt = 5000,
    inits = list(deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean))))
    

params <- c( 
		"theta0", "theta1", "theta2", "ir", "ifr", "tau", "sigma",
		"IFR_9p", "predictIFR_9p", "IFR_16p", "predictIFR_16p", "IFR_20p", "predictIFR_20p")

sampsIFR <- coda.samples(jags.modelIFR, params, 
				n.iter = MCMCn,  thin = 100, n.adapt = round(MCMCn*0.20))

saveRDS(sampsIFR, "sampsIFR_analysis.rds")

# How is MCMC mixing?

pgd <- purrr::possibly(coda::gelman.diag, list(mpsrf = NA), quiet = FALSE)
diagnostic <- pgd(sampsIFR, autoburnin = FALSE, multivariate = TRUE)$mpsrf
diagnostic

# Obtain posterior medians, HDP cred./pred. intervals:

modelparams1<-c("theta0", "theta1", "theta2", "tau", "sigma")
results1 <- data.frame(modelparams1, t(apply(cbind(modelparams1), 1, 
	function(xx){(round(summary_md_ci(xx, sampsIFR),2))})))

results1

modelparams2<-c("IFR_9p", "predictIFR_9p", "IFR_16p",
				 "predictIFR_16p", "IFR_20p", "predictIFR_20p")
results2 <- data.frame(modelparams2, t(apply(cbind(modelparams2), 1, 
	function(xx){ (round(100*summary_md_ci(xx, sampsIFR),2))})))

results2


# Obtain "study-specific" estimates for IFR and IR:

study_names <- IFRdata[,"Location"]

# IR variables
IRvars <- cbind(apply(cbind(1:K), 1, function(ii) paste("ir[",ii,"]", sep="")))
meta_IR14 <- data.frame(study_names, t(
apply(IRvars, 1, 	function(xx){round(100*summary_md_ci(xx, sampsIFR), 2)})
))

# IFR variables
IFRvars <- cbind(apply(cbind(1:K), 1, function(ii) paste("ifr[",ii,"]", sep="")))
meta_IFR14 <- data.frame(study_names, t(
apply(IFRvars, 1,   function(xx){round(100*summary_md_ci(xx, sampsIFR), 2)})
))

# Obtain "fitted values":

theta0_md <- results1[results1[,"modelparams1"]=="theta0", "md"]
theta1_md <- results1[results1[,"modelparams1"]=="theta1", "md"]
theta2_md <- results1[results1[,"modelparams1"]=="theta2", "md"]

predictIFR <- function(age65, GDPppp){
Z1_star <- (log(age65)-mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older))
Z2_star <- (log(GDPppp)-mean(log(IFRdata$GDPppp_st)))/sd(log(IFRdata$GDPppp_st))
100*( 1 - exp(-exp(theta0_md + theta1_md*(Z1_star) + theta2_md*(Z2_star)  )) )	
}


meta_IFR14[,"IFR_fitted"] <- apply(IFRdata[,c("aged_65_older", "GDPppp_st")],1, function(z) predictIFR(z[1],z[2]))
meta_IFR14


##########################################
# Plotting
##########################################

library("ggplot2")
library("gridExtra")

####  functions:
	cloglog <- function(x){log(-log(1-x))}
	icloglog <- function(x){1 - exp(-exp(x))}
	
	
####  organize results:
meta_IFR <- meta_IFR14[order(meta_IFR14[,"IFR_fitted"], IFRdata$Location),]
meta_IFR


 world_wide<-results2[results2[,"modelparams2"]=="IFR_9p",-1] 
 	
 usa_wide<-results2[results2[,"modelparams2"]=="IFR_16p",-1]

 eu_wide<-results2[results2[,"modelparams2"]=="IFR_20p",-1]

# exp(-1.313337*sd(log(IFRdata$GDPppp))+mean(log(IFRdata$GDPppp)))

meta_IFR_plus <- rbind(meta_IFR, 
data.frame(study_names="World (65yo=9%, GDP=17.8K)", rbind (c(world_wide, IFR_fitted =NA))),
data.frame(study_names="USA (65yo=16%, GDP=65.3K)", rbind (c(usa_wide, IFR_fitted =NA))),
data.frame(study_names="EU (65yo=20%, GDP=47.8K)", rbind (c(eu_wide, IFR_fitted =NA)))
)


colnames(meta_IFR_plus) <- c("Study","IFR","lower","upper", "fitted")
meta_IFR_plus<-droplevels(meta_IFR_plus)
meta_IFR_plus[,"fitted"]<-unlist(meta_IFR_plus[,"fitted"])

for(kk in 1:4){
meta_IFR_plus[,kk]<-unlist(meta_IFR_plus[,kk])
}

meta_IFR_plus <- meta_IFR_plus[order(meta_IFR_plus[,"fitted"]),]




a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=5))+ ylab("Infection Fatality Rate (%)")+ xlab("")


b <- a +geom_point(data= meta_IFR_plus , aes(x= Study ,y= fitted), col="black", pch = 4, cex=2)+ geom_pointrange(size = 0.5, shape=c(rep(20, K), 15, 15, 15))

#+geom_pointrange(data= meta_IFR_plus , aes(x= Study, y=IFR,ymax= CFR_adj,ymin= PFR_excess), col="blue", size=0.25, alpha=0.5) 

#+geom_point(data= meta_IFR_plus , aes(x= Study ,y= PFR_excess), col="blue", pch = "|", cex=2)+geom_point(data= meta_IFR_plus , aes(x= Study ,y= CFR_adj), col="blue", pch = "|", cex=2)

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b + coord_flip(ylim=c(0, 2.5))  + scale_x_discrete(limits = rev(unique(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

meta_IFR_plus[,"IFR_labs"]<-paste0(sprintf("%.2f",round(unlist(meta_IFR_plus[,"IFR"]),2)), " [",sprintf("%.2f",round(unlist(meta_IFR_plus[,"lower"]),2)),"," ,sprintf("%.2f",round(unlist(meta_IFR_plus[,"upper"]),2)),"]", sep="")

preds <- meta_IFR_plus[meta_IFR_plus [,"Study"] %in% 
				c("World (65yo=9%, GDP=18.4K)", 
                "USA (65yo=16%, GDP=65.3K)", 
                "EU (65yo=20%, GDP=47.8K)"), "IFR_labs"]

#preds <- gsub("]", "}", preds, fixed=TRUE)
#preds <- gsub("[", "{", preds, fixed=TRUE)

meta_IFR_plus[meta_IFR_plus [,"Study"] %in% c("World (65yo=9%, GDP=18.4K)", 
                "USA (65yo=16%, GDP=65.3K)", 
                "EU (65yo=20%, GDP=47.8K)"), "IFR_labs"] <- preds


IFRplot_sero <- e + geom_text(aes(label=(meta_IFR_plus [,"IFR_labs"])), hjust=-0.22+ unlist(meta_IFR_plus[,"IFR"])/2,vjust= -.75 +c(0.2,rep(0,length(unlist(meta_IFR_plus[,"IFR"]))-1)) , cex = c(rep(3, K), 3.5,3.5,3.5) )
IFRplot_sero
##########################################


##########################################
# IR variables
IRvars <- cbind(apply(cbind(1:K), 1, function(ii) paste("ir[",ii,"]", sep="")))

meta_IR <- meta_IR14[order(meta_IFR14[,"IFR_fitted"], IFRdata$Location),]

eu_wideIR <- usa_wideIR <- world_wideIR <- c(md=NA,lower=NA,higher=NA)
meta_IR_plus <- rbind(meta_IR, 
data.frame(study_names="World (65yo=9%, GDP=18.4K)", rbind (c(world_wideIR))),
data.frame(study_names="USA (65yo=16%, GDP=65.3K)", rbind (c(usa_wideIR))),
data.frame(study_names="EU (65yo=20%, GDP=47.8K)", rbind (c(eu_wideIR)))
)


colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus <-droplevels(meta_IR_plus)

for(kk in 1:4){
meta_IR_plus[,kk]<-unlist(meta_IR_plus[,kk])
}

a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=5)) + ylab("Infection Rate (%)")+ xlab("")

b <- a + geom_pointrange(size = c(rep(0.5, K), 0.75, 0.75, 0.75), shape=c(rep(20, K), 15, 15, 15))
#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip(ylim=c(0, 50)) + scale_x_discrete(limits = rev(unique(meta_IR_plus$Study)))


#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

meta_IR_plus[,"IR_labs"]<-paste0(sprintf("%.1f",round(meta_IR_plus[,"IR"],1)), " [",sprintf("%.1f",round(meta_IR_plus[,"lower"],1)),"," ,sprintf("%.1f",round(meta_IR_plus[,"upper"],1)),"]", sep="")

IRplot_sero <- e + geom_text(aes(label=(meta_IR_plus [,"IR_labs"])), hjust=-0.03 + meta_IR_plus[,"IR"]/90,vjust=-.75, cex=3) 

IRplot_sero 


##########################################
IFRplot_sero_noyaxis<- IFRplot_sero + theme(axis.title.y = element_blank(),axis.text.y = element_blank())
IRplot_sero_noyaxis<- IRplot_sero + theme(axis.title.y = element_blank(),axis.text.y = element_blank())
## Figure 2:
sero_plot <- grid.arrange(IFRplot_sero, IRplot_sero_noyaxis, ncol=2,widths=c(12,4))

# note: requires large plotting area
sero_plot

##########################################
##########################################
##########################################

# Code to replicate approximately figure from the Economist

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

cov_data[,"ifr_fitted"] <- apply(cov_data[,c("wdi_pop_over_65", "wdi_gdppc_ppp")] , 1, function(z) predictIFR(z[1],z[2]))


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
cov_data[cov_data[,"Country.Name"]=="Afghanistan", "selectnames"]<-"Afghanistan"  
cov_data[cov_data[,"Country.Name"]=="Burundi", "selectnames"]<-"Burundi"  
cov_data[cov_data[,"Country.Name"]=="Bulgaria", "selectnames"]<-"Bulgaria"  
cov_data[cov_data[,"Country.Name"]=="Ukraine", "selectnames"]<-"Ukraine"  
cov_data[cov_data[,"Country.Name"]=="Nicaragua", "selectnames"]<-"Nicaragua"  

  annotations <- data.frame(
        xpos = c(-Inf,1200,Inf,Inf),
        ypos =  c(-Inf, 1.5,-Inf,Inf),
        annotateText = c("","Size = Share of \n population over 65"
                        ,"","")) #<- adjust

cov_data[, "aborder"] <-  !is.na(cov_data[, "selectnames"])


cov_data<-data.frame(rbind(cov_data[cov_data[, "aborder"]==FALSE,], cov_data[cov_data[, "aborder"]==TRUE,]))
ggplot(cov_data, aes(x = wdi_gdppc_ppp, y = ifr_fitted, label = selectnames)) + theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_line(colour = "lightgrey", size=0.1),
    panel.grid.minor = element_line(colour = "lightgrey", size=0.1),
    panel.border = element_blank(), 
    panel.background = element_blank(), plot.title = element_text(hjust = 0.0, vjust=2.12), legend.position="none") +  geom_point(aes(size =(wdi_pop_over_65)), alpha = 0.2, col = "cornflowerblue") +   geom_point(aes(size =(wdi_pop_over_65), col = aborder), fill="blue", alpha = 0.2) + 
  scale_color_manual(values = c( "cornflowerblue","darkblue"))+
  scale_size(range = c(0.25, 10)) +  scale_x_continuous(trans = "log10", breaks=c(1000,5000,10000,50000,100000), labels=c("1,000", "5,000", "10,000", "50,000", "100,000")) + xlab("GDP per capita at purchasing-power parity (US$), log scale") +    ylab("IFR (%)") +   geom_text(size = 4, colour = "black", vjust = -1) +  ylim(0,1.75)   + 
        geom_text(data = annotations, aes(x=xpos,y=ypos,label=annotateText))+ggtitle("Expected IFR for a population with a similar age structure and GDP as:") 
  
