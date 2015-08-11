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

ldata <- get_split_an(mydata)
lmom <- t(sapply(ldata,Lmoments))

# Modeling the mean.
mydata <- set_response(mydata, 
	data.frame(no = mydata$x[,1], z = log(lmom[,'l1']) , varz =1)
	)
	
fit <- zppr(mydata, nterms = 2, criteria = 'wls')
gfit_l1 <- ppr2gam(fit, k = 5)

# Modeling the lcv.
mydata <- set_response(mydata, 
	data.frame(no = mydata$x[,1], z = lmom[,'lcv'] , varz =1)
	)
	
fit <- zppr(mydata, nterms = 1, criteria = 'wls')
gfit_lcv <- ppr2gam(fit, k = 5)

# Modeling the lca.
mydata <- set_response(mydata, 
	data.frame(no = mydata$x[,1], z = lmom[,'lca'] , varz =1)
	)
	
fit <- zppr(mydata, nterms = 2, criteria = 'wls')
gfit_lca <- ppr2gam(fit, k = 5)

# predict return level according to index flood
pred <- cbind(exp(predict(gfit_l1)),predict(gfit_lcv),
	predict(gfit_lca))
	
pars <- par.GEV(pred[,1],pred[,1]*pred[,2],pred[,3])
ret10 <- invF.GEV(.9,pars$xi,pars$alfa,pars$k)
