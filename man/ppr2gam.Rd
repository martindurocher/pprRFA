\name{ppr2gam}
\alias{ppr2gam}
\alias{plot.gppr}
\alias{terms_plot.gppr}
\alias{print.gppr}

\title{
Fit a PPR model using gam
}
\description{
The  ouput of the zppr method is used as an initial solution. A 
gam model (mgcv packages) is fitted for given direction alpha and 
the alpha are optimized by the optim function. 
}

\usage{
ppr2gam(self, k = 5, basis = "cr", m = 2, fx = FALSE, gls_tol = 1e-6, 
    gls_maxit = 10,...)

# Residuals graphics
plot(fit,...) 

# Graphics of the smooth function (see: plot.gam)
terms_plot(fit,...)

# Predicted values (see: predict.gam, residuals.gam)
predict(fit,...)
residuals(fit)

#Printing the summary of the gam object + alphas 
print(fit,...)


}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{self}{ Data object, output from new_rfa_data.}

  \item{fit}{Output from the zppr method.}
  
  \item{gls_tol}{Tolerance values for the covergeance of the model variance.}
  
  \item{gls_maxit}{Maximum number of iteration for estimating the model variance.}

  \item{\dots}{Others settings for the smooth functions (see. gam, s).}
}
\details{

See zppr for more details on the PPR model. The final output is a gam model 
with direction alpha. 

If the gls are used, the optimization is iterated until model variance 
has converged.

}
\value{

\item{gam}{A gam model for optimal alpha.}
\item{alpha}{PPR directions.}
\item{w}{Weights used in the final regression.}

}
\references{

Wood, S. (2006). Generalized additive models: an introduction with R 
(Vol. 66). Chapman & Hall/CRC.

Roosen, C. B., & Hastie, T. J. (1994). Automatic smoothing spline 
projection pursuit. Journal of Computational and Graphical Statistics, 3
(3), 235–248.

Yu, Y., & Ruppert, D. (2002). Penalized spline estimation for partially 
linear single-index models. Journal of the American Statistical 
Association, 97, 1042–1054. doi:10.1198/016214502388618861

}

\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{ppr}, link{gam}, \link{s}, \link{optim}
}
\examples{
####################################
## Regression of the return level ##
####################################

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
fit <- zppr(mydata, nterms = 2, criteria = 'gls')
print(fit)
plot(fit)
		
#Improved fitting using gam function
gfit <- ppr2gam(fit, k = 5)
print(gfit)
plot(gfit)
terms_plot(gfit)
	
#########################################
## Prediction by the index flood model ##
#########################################

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

# predict return level according to index flood assumption
pred <- cbind(exp(predict(gfit_l1)),predict(gfit_lcv),
	predict(gfit_lca))
	
pars <- par.GEV(pred[,1],pred[,1]*pred[,2],pred[,3])
print(invF.GEV(.9,pars$xi,pars$alfa,pars$k))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
