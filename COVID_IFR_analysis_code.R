##########################################
####  load required libraries and data:
library("rjags")
library("hBayesDM")
library("RCurl")
library("coda")
library("xtable")
library("rriskDistributions")
library("stringr")
library("ggplot2")
library("gridExtra")
library("plyr")

MCMCn <- 2000000

##########################################
####  functions:

cloglog <- function(x){log(-log(1-x))}
icloglog <- function(x){1 - exp(-exp(x))}
	
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
	icloglog_beta  ~ dbeta(1, 3);   # note that in the medRxiv preprint (doi: https://doi.org/10.1101/2021.05.12.21256975;) this is described erroneously as dbeta(1, 30);
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


# with alternative priors:

##########################################
####  model in JAGS:	
metaIFR_alt <- "model {

# Priors:

	icloglog_theta0 ~ dunif(0, 1); 
	icloglog_beta  ~ dunif(0, 1);
	theta0 <- log(-log(1-icloglog_theta0));
	beta  <- log(-log(1-icloglog_beta));

	inv.var_sig   <- (1/sigma)^2 ;
	inv.var_tau   <- (1/tau)^2 ;
	sigma     ~ dnorm(0, 1/100) T(0,);
	tau     ~ dnorm(0, 1/100) T(0,);
	theta1 ~ dnorm(0, 1/100);
	theta2 ~ dnorm(0, 1/100);

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


# loading the data:

sero <- read.csv("s_data_52.csv")
econ2 <- read.csv("REGION_ECONOM_11082021015737664.csv")


# merging the OECD GDP data with our seroprevalence dataset:

for(rr  in sero[sero[,c("oecd_region")]%in%c(1), "Region"]){

	my_region <- as.character(rr)

	print("     ")
	print(my_region)
	print("     ")

	print(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"Year"]%in%max(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,"Year"])  & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,])

	sero[sero[,c("Region")]%in%c(rr), "OECD_GDP"] <- mean(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"Year"]%in%max(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,"Year"])  & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,"Value"])

}


for(rr  in sero[sero[,c("oecd_region")]%in%c(2), "Region"]){

	region_names <- as.character(rr)

	my_region <- unique(econ2[grep(region_names,(econ2$Region)),"Region"])
	my_region <- c(as.character(my_region))

	print("     ")
	print(my_region)
	print("     ")


	print(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"Year"]%in%max(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,"Year"])  & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,])

	sero[sero[,c("Region")]%in%c(rr), "OECD_GDP"] <- mean(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"Year"]%in%max(econ2[econ2[,"Region"]%in%c(my_region) & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,"Year"]) & econ2[,"MEAS"]%in% c("PC_USD_PPP") ,"Value"])

}

sero[,"GDP_OECD_original"] <- sero[,c("GDPppp")]
sero[!is.na(sero[,"OECD_GDP"]), "OECD_GDP_source"] <- 1
sero[is.na(sero[,"OECD_GDP"]), "OECD_GDP_source"] <- 0
sero[is.na(sero[,"OECD_GDP"]),"OECD_GDP"] <- sero[is.na(sero[,"OECD_GDP"]),"GDPppp"]
sero_data <- sero
sero_data[,c("GDPppp")] <- sero_data[,c("OECD_GDP")]
sero[sero[,"Excluded"]==0 ,c("Author", "Location", "OECD_GDP_source")]

##########

##############################
##############################
# Chen et al.  - based  analysis
##############################
##############################

dim(sero_data[sero_data[,"Dataset_Chen"]==1,])
dim(sero_data[sero_data[,"Dataset_Chen"]==1 & sero_data[,"Excluded_screening"]==1,])
table(sero_data[sero_data[,"Dataset_Chen"]==1 & sero_data[,"Excluded_screening"]==1,"exclusion_reason"])
dim(sero_data[sero_data[,"Dataset_Chen"]==1 & sero_data[,"Excluded_data"]==1,])
table(sero_data[sero_data[,"Dataset_Chen"]==1 & sero_data[,"Excluded_data"]==1,"exclusion_reason"])
dim(sero_data[sero_data[,"Dataset_Chen"]==1 & sero_data[,"Excluded"]==1,])


IFRdata_chen <- sero_data[sero_data[,"Dataset_Chen"]==1,]


####  Making Table 1 and Table 2 ####

look0 <- (IFRdata_chen[,c("Author", "Location", "Excluded", "Notes_exclusion", "Location", "Population",  "date_start", "date_end",   "IR_lower",  "IR_upper", "deaths14_lower", "deaths14_upper", "aged_65_older", "GDPppp", "Country", "Prevalence.Estimate.Name", "Dataset_both", "death_data_quality")])

look1<-(look0[look0[,"Excluded"]==0,])


#look0 <- (IFRdata_chen[,c("Author", "Location", "Excluded", "Notes_exclusion", "Location", "Population",  "date_start", "date_end",   "IR_lower",  "IR_upper", "deaths14_lower", "deaths14_upper", "aged_65_older", "GDPppp", "Country", "Prevalence.Estimate.Name", "death_data_quality")])

#look1<-(look0[look0[,"Excluded"]==0 & look0[,"death_data_quality"]==1,])

excluded_studies<- look0[look0[,"Excluded"]==1,c("Author", "Location", "Notes_exclusion")]

print(xtable(excluded_studies[order(excluded_studies[,"Author"]),]), include.rownames=FALSE)

look1[,"IR_lower"] <- as.numeric(as.character(look1[,"IR_lower"] ))
look1[,"IR_upper"] <- as.numeric(as.character(look1[,"IR_upper"]))

for(x in 1:dim(look1)[1]){
mytol<-0.001
IR_CI <-(c(unlist(look1[x,c("IR_lower", "IR_upper")]/100)))
ab_param <- get.beta.par(p = c(0.025,0.975), q = IR_CI, plot = FALSE, tol = mytol)

while(is.na(ab_param[1])){
mytol <- mytol+0.001
ab_param <- get.beta.par(p = c(0.025,0.975), q = IR_CI, plot = FALSE, tol = mytol)
}
look1[x,"total_tests"] <- round(ab_param)[1]+round(ab_param)[2]+1
look1[x,"total_cases"] <- round(ab_param)[1]
}


look2 <- look1

look2[,"Authors"] <- look2[,"Author"]

look2[,"study_names"] <-make.unique(as.character(look2[ ,"Location"]))

look2[,"date_start_lab"] <- apply(cbind(1:dim(look2)[1]),1,  function(qq) {paste(c(str_split(as.character(look2[qq,"date_start"]), "-")[[1]][2:3]), collapse="/")})


look2[,"date_end_lab"] <- apply(cbind(1:dim(look2)[1]),1,  function(qq) {paste(c(str_split(as.character(look2[qq,"date_end"]), "-")[[1]][2:3]), collapse="/")})

look2[,"window"] <- apply(cbind(1:dim(look2)[1]), 1, function(qq) {paste(look2[qq,c("date_start_lab", "date_end_lab")], collapse=" - ")})

look2[,"IRinterval"] <-apply(cbind(1:dim(look2)[1]),1, function(qq) {paste("[",paste(sprintf("%.2f", look2[qq,c("IR_lower", "IR_upper")]), collapse=", "),"]", sep="")})


look2[,"author_with_citation"] <- apply(cbind(1:dim(look2)[1]), 1, function(qq) { paste(look2[qq,c("Authors")], "...", look2[qq,"Prevalence.Estimate.Name"], "...", collapse="")})

look2[,"Dataset_both"][look2[,"Dataset_both"]==1]<-"yes"
look2[,"Dataset_both"][look2[,"Dataset_both"]=="0"]<-""

#Table 1
print(xtable(look2[order(look2[,"Authors"]),][,c("author_with_citation", "study_names", "window", "IRinterval", "Dataset_both")], digits =c(rep(0,3), 2,2,1)), include.rownames=FALSE)

#Table 2
print(xtable(look2[order(look2[,"Authors"]),][,c("study_names", "Population",  "deaths14_lower" ,  "deaths14_upper", "total_tests", "total_cases", "aged_65_older", "GDPppp")], digits = 0), include.rownames=FALSE)


IFRdata <- look2

K <- length(IFRdata$total_tests)
K

# For sensitivity analysis
 IFRdata <- look2[look2[,"death_data_quality"]==1,]
 K <- length(IFRdata$total_tests)
 K

#########

# GDP PPP per capita 
# World = 17,811.3;     
# US = 65,297.5
# EU = 47,827.7	

gdp_values <- c(
(log(17811.3) - mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp)),
(log(65297.5) - mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp)),
(log(47827.7) - mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp)))

# Age over 65
# World = 9    
# US = 16   
# EU = 20	

age65_values <-  c(
(log(9) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)),
(log(16) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)),
(log(20) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)))

# Initial values for MCMC random seeds

	inits1<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(123), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits2<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(456), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits3<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(789), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits4<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(999), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits5<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(111), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))



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
    Z2 = c(scale(log(IFRdata$GDPppp))),
    censor.index = rep(1, K),
    age65_values = age65_values,
    gdp_values = gdp_values
    ), 
    n.chains = 5, 
    n.adapt = 5000,
    inits = list(inits1, inits2, inits3, inits4, inits5))
    

params <- c( 
		"theta0", "theta1", "theta2", "ir", "ifr", "tau", "sigma",
		"IFR_9p", "predictIFR_9p", "IFR_16p", "predictIFR_16p", "IFR_20p", "predictIFR_20p")

sampsIFR <- coda.samples(jags.modelIFR, params, 
				n.iter = MCMCn,  thin = 100, n.adapt = round(MCMCn*0.20))

#saveRDS(sampsIFR, "chen_1.rds")

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


results1_chen <- results1
results2_chen <- results2


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
Z2_star <- (log(GDPppp)-mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp))
100*( 1 - exp(-exp(theta0_md + theta1_md*(Z1_star) + theta2_md*(Z2_star)  )) )	
}


meta_IFR14[,"IFR_fitted"] <- apply(IFRdata[,c("aged_65_older", "GDPppp")],1, function(z) predictIFR(z[1],z[2]))
meta_IFR14


meta_IR14

##########################################
# Plotting
##########################################

library("ggplot2")
library("gridExtra")

####  functions:
	cloglog <- function(x){log(-log(1-x))}
	icloglog <- function(x){1 - exp(-exp(x))}
	
	
####  organize results:
meta_IFR <- meta_IFR14[order(meta_IFR14[,"IFR_fitted"], meta_IFR14[,"study_names"]),]
meta_IFR


####  organize results:
meta_IR <- meta_IR14[order(meta_IFR14[,"IFR_fitted"], meta_IFR14[,"study_names"]),]
meta_IR

 world_wide<-results2[results2[,"modelparams2"]=="IFR_9p",-1] 
 	
 usa_wide<-results2[results2[,"modelparams2"]=="IFR_16p",-1]

 eu_wide<-results2[results2[,"modelparams2"]=="IFR_20p",-1]

# exp(-1.313337*sd(log(IFRdata$GDPppp))+mean(log(IFRdata$GDPppp)))

meta_IFR_plus <- rbind(meta_IFR, 
data.frame(study_names="World (65yo=9%, GDP=17.8k)", rbind (c(world_wide, IFR_fitted =NA))),
data.frame(study_names="USA (65yo=16%, GDP=65.3k)", rbind (c(usa_wide, IFR_fitted =NA))),
data.frame(study_names="EU (65yo=20%, GDP=47.8k)", rbind (c(eu_wide, IFR_fitted =NA)))
)

colnames(meta_IFR_plus) <- c("Study","IFR","lower","upper", "fitted")
meta_IFR_plus<-droplevels(meta_IFR_plus)
meta_IFR_plus[,"fitted"]<-unlist(meta_IFR_plus[,"fitted"])

for(kk in 1:4){
meta_IFR_plus[,kk]<-unlist(meta_IFR_plus[,kk])
}

meta_IFR_plus <- meta_IFR_plus[order(meta_IFR_plus[,"fitted"]),]


a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=5))+ ylab("Infection Fatality Rate (%)")+ xlab("")


b <- a +geom_point(data= meta_IFR_plus , aes(x= Study ,y= fitted), col="black", pch = 4, cex=2)+ geom_pointrange(size = 0.5, shape=c(rep(20, K), 4, 4, 4))

#+geom_pointrange(data= meta_IFR_plus , aes(x= Study, y=IFR,ymax= CFR_adj,ymin= PFR_excess), col="blue", size=0.25, alpha=0.5) 

#+geom_point(data= meta_IFR_plus , aes(x= Study ,y= PFR_excess), col="blue", pch = "|", cex=2)+geom_point(data= meta_IFR_plus , aes(x= Study ,y= CFR_adj), col="blue", pch = "|", cex=2)

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b + coord_flip(ylim=c(0, 2.5))  + scale_x_discrete(limits = rev(unique(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

meta_IFR_plus[,"IFR_labs"]<-paste0(sprintf("%.2f",round(unlist(meta_IFR_plus[,"IFR"]),2)), " [",sprintf("%.2f",round(unlist(meta_IFR_plus[,"lower"]),2)),"," ,sprintf("%.2f",round(unlist(meta_IFR_plus[,"upper"]),2)),"]", sep="")

preds <- meta_IFR_plus[meta_IFR_plus [,"Study"] %in% 
				c("World (65yo=9%, GDP=18.4k)", 
                "USA (65yo=16%, GDP=65.3k)", 
                "EU (65yo=20%, GDP=47.8k)"), "IFR_labs"]

#preds <- gsub("]", "}", preds, fixed=TRUE)
#preds <- gsub("[", "{", preds, fixed=TRUE)

meta_IFR_plus[meta_IFR_plus [,"Study"] %in% c("World (65yo=9%, GDP=18.4k)", 
                "USA (65yo=16%, GDP=65.3k)", 
                "EU (65yo=20%, GDP=47.8k)"), "IFR_labs"] <- preds


IFRplot_chen <- e + geom_text(aes(label=(meta_IFR_plus [,"IFR_labs"])), hjust = -0.22+ unlist(meta_IFR_plus[,"IFR"])/2, vjust = -.75 +c(0.2,rep(0,length(unlist(meta_IFR_plus[,"IFR"]))-1)) , cex = c(rep(2.5, K), 3.25,3.25,3.25) )
IFRplot_chen
##########################################

## IR plot
####  organize results:

meta_IR_plus <- rbind(meta_IR)
colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus<-droplevels(meta_IR_plus)

a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=5))+ ylab("Infection Rate (%)")+ xlab("")

b <- a + geom_pointrange(size = 0.5, shape=c(rep(20, K)))

cc <- b + coord_flip(ylim=c(0, 30))  + scale_x_discrete(limits = rev(unique(meta_IR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()
e

##############################
##############################
# Serotracker  - based  analysis
##############################
##############################

dim(sero_data[sero_data[,"Dataset_Serotracker"]==1,])
dim(sero_data[sero_data[,"Dataset_Serotracker"]==1 & sero_data[,"Excluded_screening"]==1,])
table(sero_data[sero_data[,"Dataset_Serotracker"]==1 & sero_data[,"Excluded_screening"]==1,"exclusion_reason"])
dim(sero_data[sero_data[,"Dataset_Serotracker"]==1 & sero_data[,"Excluded_data"]==1,])
table(sero_data[sero_data[,"Dataset_Serotracker"]==1 & sero_data[,"Excluded_data"]==1,"exclusion_reason"])
dim(sero_data[sero_data[,"Dataset_Serotracker"]==1 & sero_data[,"Excluded"]==1,])

IFRdata_sero <- sero_data[sero_data[,"Dataset_Serotracker"]==1,]


####  Making Table 1 and Table 2 ####

look0 <- (IFRdata_sero[,c("Author", "Location", "Excluded", "Notes_exclusion", "Location", "Population",  "date_start", "date_end",   "IR_lower",  "IR_upper", "deaths14_lower", "deaths14_upper", "aged_65_older", "GDPppp", "Country", "Prevalence.Estimate.Name", "Dataset_both", "death_data_quality")])

look1<-(look0[look0[,"Excluded"]==0,])

#look0 <- (IFRdata_sero[,c("Author", "Location", "Excluded", "Notes_exclusion", "Location", "Population",  "date_start", "date_end",   "IR_lower",  "IR_upper", "deaths14_lower", "deaths14_upper", "aged_65_older", "GDPppp", "Country", "Prevalence.Estimate.Name", "death_data_quality")])

#look1<-(look0[look0[,"Excluded"]==0 & look0[,"death_data_quality"]==1,])


excluded_studies<- look0[look0[,"Excluded"]==1,c("Author", "Location", "Notes_exclusion")]

print(xtable(excluded_studies[order(excluded_studies[,"Author"]),]), include.rownames=FALSE)

look1[,"IR_lower"] <- as.numeric(as.character(look1[,"IR_lower"] ))
look1[,"IR_upper"] <- as.numeric(as.character(look1[,"IR_upper"]))

for(x in 1:dim(look1)[1]){
mytol<-0.001
IR_CI <-(c(unlist(look1[x,c("IR_lower", "IR_upper")]/100)))
ab_param <- get.beta.par(p = c(0.025,0.975), q = IR_CI, plot = FALSE, tol = mytol)

while(is.na(ab_param[1])){
mytol <- mytol+0.001
ab_param <- get.beta.par(p = c(0.025,0.975), q = IR_CI, plot = FALSE, tol = mytol)
}
look1[x,"total_tests"] <- round(ab_param)[1]+round(ab_param)[2]+1
look1[x,"total_cases"] <- round(ab_param)[1]
}


look2 <- look1

look2[,"Authors"] <- look2[,"Author"]

look2[,"study_names"] <-make.unique(as.character(look2[ ,"Location"]))

look2[,"date_start_lab"] <- apply(cbind(1:dim(look2)[1]),1,  function(qq) {paste(c(str_split(as.character(look2[qq,"date_start"]), "-")[[1]][2:3]), collapse="/")})


look2[,"date_end_lab"] <- apply(cbind(1:dim(look2)[1]),1,  function(qq) {paste(c(str_split(as.character(look2[qq,"date_end"]), "-")[[1]][2:3]), collapse="/")})

look2[,"window"] <- apply(cbind(1:dim(look2)[1]), 1, function(qq) {paste(look2[qq,c("date_start_lab", "date_end_lab")], collapse=" - ")})

look2[,"IRinterval"] <-apply(cbind(1:dim(look2)[1]),1, function(qq) {paste("[",paste(sprintf("%.2f", look2[qq,c("IR_lower", "IR_upper")]), collapse=", "),"]", sep="")})



look2[,"author_with_citation"] <- apply(cbind(1:dim(look2)[1]), 1, function(qq) { paste(look2[qq,c("Authors")], "...", look2[qq,"Prevalence.Estimate.Name"], "...", collapse="")})

look2[,"Dataset_both"][look2[,"Dataset_both"]==1]<-"yes"
look2[,"Dataset_both"][look2[,"Dataset_both"]=="0"]<-""

#Table 1
print(xtable(look2[order(look2[,"Authors"]),][,c("author_with_citation", "study_names", "window", "IRinterval", "Dataset_both")], digits =c(rep(0,3), 2,2,1)), include.rownames=FALSE)


#Table 2
print(xtable(look2[order(look2[,"Authors"]),][,c("study_names", "Population",  "deaths14_lower" ,  "deaths14_upper", "total_tests", "total_cases", "aged_65_older", "GDPppp")], digits = 0), include.rownames=FALSE)





IFRdata <- look2
K <- length(IFRdata$total_tests)
K

# For sensitivity analysis
 IFRdata <- look2[look2[,"death_data_quality"]==1,]
 K <- length(IFRdata$total_tests)
 K

#########

# GDP PPP per capita 
# World = 17,811.3;     
# US = 65,297.5
# EU = 47,827.7	

gdp_values <- c(
(log(17811.3) - mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp)),
(log(65297.5) - mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp)),
(log(47827.7) - mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp)))

# Age over 65
# World = 9    
# US = 16   
# EU = 20	

age65_values <-  c(
(log(9) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)),
(log(16) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)),
(log(20) - mean(log(IFRdata$aged_65_older)))/sd(log(IFRdata$aged_65_older)))

# Initial values for MCMC random seeds

	inits1<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(123), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits2<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(456), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits3<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(789), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits4<- list(.RNG.name = "base::Super-Duper", 
                  .RNG.seed = c(999), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))
	inits5<- list(.RNG.name = "base::Wichmann-Hill", 
                  .RNG.seed = c(111), deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean)))



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
    Z2 = c(scale(log(IFRdata$GDPppp))),
    censor.index = rep(1, K),
    age65_values = age65_values,
    gdp_values = gdp_values
    ), 
    n.chains = 5, 
    n.adapt = 5000,
    inits = list(inits1, inits2, inits3, inits4, inits5))
    

params <- c( 
		"theta0", "theta1", "theta2", "ir", "ifr", "tau", "sigma",
		"IFR_9p", "predictIFR_9p", "IFR_16p", "predictIFR_16p", "IFR_20p", "predictIFR_20p")

sampsIFR <- coda.samples(jags.modelIFR, params, 
				n.iter = MCMCn,  thin = 100, n.adapt = round(MCMCn*0.20))


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

results1_sero <- results1
results2_sero <- results2

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
Z2_star <- (log(GDPppp)-mean(log(IFRdata$GDPppp)))/sd(log(IFRdata$GDPppp))
100*( 1 - exp(-exp(theta0_md + theta1_md*(Z1_star) + theta2_md*(Z2_star)  )) )	
}


meta_IFR14[,"IFR_fitted"] <- apply(IFRdata[,c("aged_65_older", "GDPppp")],1, function(z) predictIFR(z[1],z[2]))
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
meta_IFR <- meta_IFR14[order(meta_IFR14[,"IFR_fitted"], meta_IFR14[,"study_names"]),]
meta_IFR

####  organize results:
meta_IR <- meta_IR14[order(meta_IFR14[,"IFR_fitted"], meta_IFR14[,"study_names"]),]
meta_IR


 world_wide<-results2[results2[,"modelparams2"]=="IFR_9p",-1] 
 	
 usa_wide<-results2[results2[,"modelparams2"]=="IFR_16p",-1]

 eu_wide<-results2[results2[,"modelparams2"]=="IFR_20p",-1]

# exp(-1.313337*sd(log(IFRdata$GDPppp))+mean(log(IFRdata$GDPppp)))

meta_IFR_plus <- rbind(meta_IFR, 
data.frame(study_names="World (65yo=9%, GDP=17.8k)", rbind (c(world_wide, IFR_fitted =NA))),
data.frame(study_names="USA (65yo=16%, GDP=65.3k)", rbind (c(usa_wide, IFR_fitted =NA))),
data.frame(study_names="EU (65yo=20%, GDP=47.8k)", rbind (c(eu_wide, IFR_fitted =NA)))
)

colnames(meta_IFR_plus) <- c("Study","IFR","lower","upper", "fitted")
meta_IFR_plus<-droplevels(meta_IFR_plus)
meta_IFR_plus[,"fitted"]<-unlist(meta_IFR_plus[,"fitted"])

for(kk in 1:4){
meta_IFR_plus[,kk]<-unlist(meta_IFR_plus[,kk])
}

meta_IFR_plus <- meta_IFR_plus[order(meta_IFR_plus[,"fitted"]),]



a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=5))+ ylab("Infection Fatality Rate (%)")+ xlab("")

b <- a +geom_point(data= meta_IFR_plus , aes(x= Study ,y= fitted), col="black", pch = 4, cex=2)+ geom_pointrange(size = 0.5, shape=c(rep(20, K), 4, 4, 4))

#+geom_pointrange(data= meta_IFR_plus , aes(x= Study, y=IFR,ymax= CFR_adj,ymin= PFR_excess), col="blue", size=0.25, alpha=0.5) 

#+geom_point(data= meta_IFR_plus , aes(x= Study ,y= PFR_excess), col="blue", pch = "|", cex=2)+geom_point(data= meta_IFR_plus , aes(x= Study ,y= CFR_adj), col="blue", pch = "|", cex=2)

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b + coord_flip(ylim=c(0, 4.5))  + scale_x_discrete(limits = rev(unique(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

meta_IFR_plus[,"IFR_labs"]<-paste0(sprintf("%.2f",round(unlist(meta_IFR_plus[,"IFR"]),2)), " [",sprintf("%.2f",round(unlist(meta_IFR_plus[,"lower"]),2)),"," ,sprintf("%.2f",round(unlist(meta_IFR_plus[,"upper"]),2)),"]", sep="")

preds <- meta_IFR_plus[meta_IFR_plus [,"Study"] %in% 
				c("World (65yo=9%, GDP=18.4k)", 
                "USA (65yo=16%, GDP=65.3k)", 
                "EU (65yo=20%, GDP=47.8k)"), "IFR_labs"]

#preds <- gsub("]", "}", preds, fixed=TRUE)
#preds <- gsub("[", "{", preds, fixed=TRUE)

meta_IFR_plus[meta_IFR_plus [,"Study"] %in% c("World (65yo=9%, GDP=18.4k)", 
                "USA (65yo=16%, GDP=65.3k)", 
                "EU (65yo=20%, GDP=47.8k)"), "IFR_labs"] <- preds


IFRplot_sero <- e + geom_text(aes(label=(meta_IFR_plus [,"IFR_labs"])), hjust = -0.22+ unlist(meta_IFR_plus[,"IFR"])/2, vjust = -.75 +c(0.2,rep(0,length(unlist(meta_IFR_plus[,"IFR"]))-1)) , cex = c(rep(2.5, K), 3.25,3.25,3.25) )
IFRplot_sero
##########################################

##########################################

## IR plot
####  organize results:

meta_IR_plus <- rbind(meta_IR)
colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus<-droplevels(meta_IR_plus)

a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=5))+ ylab("Infection Rate (%)")+ xlab("")

b <- a + geom_pointrange(size = 0.5, shape=c(rep(20, K)))

cc <- b + coord_flip(ylim=c(0, 80))  + scale_x_discrete(limits = rev(unique(meta_IR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

e
