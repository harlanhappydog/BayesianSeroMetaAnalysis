##########################################
####  load required libraries and data:
library("rjags")
library("hBayesDM")
library("RCurl")

csvfile <- getURL("https://raw.githubusercontent.com/harlanhappydog/BayesianSeroMetaAnalysis/main/IFRdata_withGDP.csv")
IFRdata <- read.csv(text = csvfile)

MCMCn <- 1000000

##########################################
####  functions:
summary_md_ci <- function(xx, samps){
c(  md=summary(samps)$quantiles[xx, "50%"], 
	lower=HDIofMCMC(unlist(((samps[[1]])[,xx])))[1], 
	higher=HDIofMCMC(unlist(((samps[[1]])[,xx])))[2]
	)}
	
##########################################
####  model in JAGS:	
metaIFR <- "model {

# Priors:

	icloglog_theta ~ dbeta(0.3, 30); 
	icloglog_beta  ~ dbeta(1, 3);
	theta <- log(-log(1-icloglog_theta));
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
	cloglog_ifr[k] ~ dnorm(theta + theta1*Z1[k] + theta2*Z2[k], inv.var_tau);
	
  }
  
# Summary

g_IFR_9p = theta + theta1*(-1.041098) +  theta2*(-1.313337);  
g_IFR_16p = theta + theta1*(0.1890347) + theta2*(0.4797427);  
g_IFR_20p = theta + theta1*(0.6661173) + theta2*(0.03933778);  

IFR_ <-  1 - exp(-exp(theta))
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
K <- length(IFRdata$total_tests)

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
    censor.index = rep(1, K)
    ), 
    n.chains = 5, 
    n.adapt = 5000,
    inits = list(deaths = round(apply(cbind(IFRdata$deaths14_lower, IFRdata$deaths14_upper), 1, mean))))

params <- c("IFR_9p", "predictIFR_9p", "IFR_16p", 
		"predictIFR_16p", "IFR_20p", "predictIFR_20p", 
		"icloglog_theta", "theta", "theta1", "theta2", "ir", "ifr", "tau", "sigma")
sampsIFR <- coda.samples(jags.modelIFR, params, 
				n.iter = MCMCn,  thin = 100, n.adapt = round(MCMCn*0.10))


# obtain posterior medians, HDP cred./pred. intervals:
modelparams1<-c("theta", "theta1", "theta2", "tau", "sigma")
results1 <- data.frame(modelparams1, t(apply(cbind(modelparams1), 1, 
	function(xx){(round(summary_md_ci(xx, sampsIFR),2))})))

modelparams2<-c("IFR_9p", "predictIFR_9p", "IFR_16p",
				 "predictIFR_16p", "IFR_20p", "predictIFR_20p")
results2 <- data.frame(modelparams2, t(apply(cbind(modelparams2), 1, 
	function(xx){ (round(100*summary_md_ci(xx, sampsIFR),2))})))

results1
results2

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


##########################################
# Plotting
##########################################

library("ggplot2")
library("gridExtra")

####  functions:
	cloglog <- function(x){log(-log(1-x))}
	icloglog <- function(x){1 - exp(-exp(x))}
	
	
####  organize results:
meta_IFR <- meta_IFR14[order(IFRdata$aged_65_older, IFRdata$Location),]

 world_wide<-results2[results2[,"modelparams2"]=="predictIFR_9p",-1] 
 	
 usa_wide<-results2[results2[,"modelparams2"]=="predictIFR_16p",-1]

 eu_wide<-results2[results2[,"modelparams2"]=="predictIFR_20p",-1]

# exp(-1.313337*sd(log(IFRdata$GDPppp))+mean(log(IFRdata$GDPppp)))

meta_IFR_plus <- rbind(meta_IFR, 
data.frame(study_names="World (65yo=9%, GDP(ppp)=18.4K)", rbind (c(world_wide))),
data.frame(study_names="USA (65yo=16%, GDP(ppp)=65.3K)", rbind (c(usa_wide))),
data.frame(study_names="EU (65yo=20%, GDP(ppp)=47.8K)", rbind (c(eu_wide)))
)

colnames(meta_IFR_plus) <- c("Study","IFR","lower","upper")
meta_IFR_plus<-droplevels(meta_IFR_plus)

for(kk in 1:4){
meta_IFR_plus[,kk]<-unlist(meta_IFR_plus[,kk])
}


a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=5))+ ylab("Infection Fatality Rate (%)")+ xlab("")


b <- a + geom_pointrange(size = c(rep(0.5, K), 0.75, 0.75, 0.75), shape=c(rep(20, K), 15, 15, 15))


#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b + coord_flip(ylim=c(0, 2.5))  + scale_x_discrete(limits = rev(unique(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

meta_IFR_plus[,"IFR_labs"]<-paste0(sprintf("%.2f",round(unlist(meta_IFR_plus[,"IFR"]),2)), " [",sprintf("%.2f",round(unlist(meta_IFR_plus[,"lower"]),2)),"," ,sprintf("%.2f",round(unlist(meta_IFR_plus[,"upper"]),2)),"]", sep="")

preds <- meta_IFR_plus[meta_IFR_plus [,"Study"] %in% 
				c("World (65yo=9%, GDP(ppp)=18.4K)", 
                "USA (65yo=16%, GDP(ppp)=65.3K)", 
                "EU (65yo=20%, GDP(ppp)=47.8K)"), "IFR_labs"]

preds <- gsub("]", "}", preds, fixed=TRUE)
preds <- gsub("[", "{", preds, fixed=TRUE)

meta_IFR_plus[meta_IFR_plus [,"Study"] %in% c("World (65yo=9%, GDP(ppp)=18.4K)", 
                "USA (65yo=16%, GDP(ppp)=65.3K)", 
                "EU (65yo=20%, GDP(ppp)=47.8K)"), "IFR_labs"] <- preds


IFRplot_sero <- e + geom_text(aes(label=(meta_IFR_plus [,"IFR_labs"])), hjust=-0.22+ unlist(meta_IFR_plus[,"IFR"])/2,vjust= -.75 +c(0.2,rep(0,length(unlist(meta_IFR_plus[,"IFR"]))-1)) , cex = c(rep(3, K), 3.5,3.5,3.5) )
IFRplot_sero
##########################################


##########################################
# IR variables
IRvars <- cbind(apply(cbind(1:K), 1, function(ii) paste("ir[",ii,"]", sep="")))

meta_IR <- meta_IR14[order(IFRdata$aged_65_older, IFRdata$Location),]

eu_wideIR <- usa_wideIR <- world_wideIR <- c(md=NA,lower=NA,higher=NA)
meta_IR_plus <- rbind(meta_IR, 
data.frame(study_names="World (65yo=9%, GDP(ppp)=18.4K)", rbind (c(world_wideIR))),
data.frame(study_names="USA (65yo=16%, GDP(ppp)=65.3K)", rbind (c(usa_wideIR))),
data.frame(study_names="EU (65yo=20%, GDP(ppp)=47.8K)", rbind (c(eu_wideIR)))
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
