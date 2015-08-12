#
# objective: Personal program for testing
#

rm(list = ls())
set.seed(1)

library(mgcv)
library(evd)
source('~/github/pprRFA/R/libpprRFA.R')

library(nsRFA)
data(hydroSIMN)

USE_CACHE <- TRUE
CACHE_FILE <- '~/github/pprRFA/inst/cache.Rdata'
 
attrib <- data.frame(parameters[,1],scale(log(parameters[,3:6])))
names(attrib) <- names(parameters[,c(1,3:6)])

mydata <- new_rfa_data(annualflows,attrib)
	
#Calculate estimate of the return level of 10 years
# at logarithm scale by bootstrap and its variance
# NOTE that nboot is low to limit computing time.
if(USE_CACHE) {
	load(CACHE_FILE)
} else {
	z <- at_site_boot_log(mydata, ret = 10, nboot = 50, 
		control = list(maxit=500))
	save(z, file = CACHE_FILE)
}

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

#####################################
# PPR fitting without using rfa object 
X <- mydata$x[,-1]

y <- z[,2]
tmp <- zppr(y, X, nterms = 1)

y <- z[,-1]
tmp <- zppr(y, X, nterms = 2, criteria = 'ols')

gfit <- ppr2gam(tmp, k = 5, subs = seq(20))
print(gfit)
print(predict(gfit))
print(predict(gfit,X[21:30,], se.fit = TRUE))

