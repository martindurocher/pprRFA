#load in packages
rm(list = ls())
library(mgcv)
library(evd)
source('~/github/pprRFA/R/libpprRFA.R')

library(nsRFA)
data(hydroSIMN)

attrib <- data.frame(parameters[,1],scale(log(parameters[,3:6])))
names(attrib) <- names(parameters[,c(1,3:6)])

mydata <- new_rfa_data(annualflows,attrib)
	
#Calculate estimate of the return level of 10 years
# at logarithm scale by bootstrap and its variance
# NOTE that nboot is low to limit computing time.
z <- at_site_boot_log(mydata, ret = 10, nboot = 30, 
	control = list(maxit=500))
mydata <- set_response(mydata,z)
		
#PPR initial fitting
fit <- zppr(mydata, nterms = 2, criteria = 'wls')
print(fit)
par(mfrow = c(2,2))
plot(fit)
		
#Improved fitting using gam function
gfit <- ppr2gam(fit, k = 5)
print(gfit)
plot(gfit)
terms_plot(gfit)

# PPR fitting without using mydata  
fit <- zppr(z[,-1], mydata$x[,-1], nterms = 2, criteria = 'ols')
print(fit)
gfit <- ppr2gam(fit, k = 5)
print(gfit)
plot(gfit)
terms_plot(gfit)
