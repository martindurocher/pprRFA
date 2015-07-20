new_rfa_data <- function(annual,attrib)
{
	# verify format of the 
	y <- as.data.frame(annual)
	if( ncol(y) != 3)
		return(-1)
	else
		names(y) <- c('no','year','ann')
		
	# verify that the stations number matching between
	# annual and attrib
	tmp <- unique(annual[,1])
	if( sum(tmp == attrib[,1]) != length(attrib[,1]) )
		return(-2)
		
	# calculate the record length of each stations
	nsite <- sapply(split(y[,3],y[,1]),length)
		
	self <- list(x=attrib, y = annual, n = length(tmp), 
		response = NULL, nsite= nsite, 
		p = ncol(attrib) - 1)
	
	class(self) <- 'rfa'
	return(self)
}

get_split_an <- function(x) UseMethod('get_split_an',x)
get_split_an.rfa <- function(self) split(self$y[,3],self$y[,1])

at_site_boot_log <- function(x,...) UseMethod('at_site_boot_log',x)
at_site_boot_log.rfa <- function(self, ret = 2, nboot=500, ...)
{

	ylist <- get_split_an(self)
	m <- length(ylist)
	ans <- data.frame(z = rep(0,m), var = rep(0,m))
	
	yboot <- rep(NA,nboot)
	
	bar <- txtProgressBar(1,m*nboot)
	for(i in seq(m))
	{
		yi <- ylist[[i]]
		
		for(k in seq(nboot))
		{
			setTxtProgressBar(bar,(i-1)*nboot+k)
			
			yi_sub <- sample(yi,length(yi),replace=TRUE)
			
			fit <- fgev(yi_sub, std.err = FALSE, ...)
			
			yboot[k] <- log(qgev(1-1/ret,fit$estimate[1],
				fit$estimate[2],fit$estimate[3]))
		}
		
		ans[i,1] <- self$x[i,1]
		ans[i,2] <- mean(yboot)
		ans[i,3] <- var(yboot)	
	}
	
	return(ans)
}

set_response <- function(x,...) UseMethod('set_response',x)
set_response.rfa <- function(self,z)
{
	#adding or change a response value from rfa data
	z <- as.data.frame(z)
	names(z)<-c('no','z','var')
	
	if( sum(z$no != self$x[,1]) > 0 )
		return(NA)
	else
	{
		self$response <- z
		return(self)
	}
}

zppr <- function(x,...) UseMethod('zppr',x)
zppr.rfa <- function(self, w=NULL, criteria = 'wls', nterms = 1, ...)
{
	#wrapper for basic ppr functions in stats packages
	if(criteria == 'ols' )
	{
		fit <- ppr(self$response[,2]~., data=self$x[,-1], 
			nterms = nterms,...)
	}
	else if(criteria == 'wls')
	{	
		if(is.null(w))
			w <- self$nsite
		
		fit <- ppr(self$response[,2]~., data=self$x[,-1], nterms = nterms,
			weights = w,...)
	}
	else if(criteria == 'gls')
	{
		fit <- ppr(self$response[,2]~., data=self$x[,-1], 
			nterms = nterms, ...)
			
		nloop <- 100
		sig2 <- var(residuals(fit))
		
		w <- 1/(1 + self$response[,3]/sig2)
		
		fit <- ppr(self$response[,2]~., data=self$x[,-1], nterms = nterms,
			weights = w,...)
			
	}
	
	sig2 <- var(residuals(fit))
	
	ans <- list(response = self$response[,2], fit = fit, alpha = fit$alpha,
		w = w, criteria = criteria, x = self$x[,-1], beta = fit$beta,
		sample_var = self$response[,3], sig2 = sig2, nterms = nterms,
		p = self$p )
		
	class(ans) <- append(class(ans),'zppr')	
	
	return(ans)
}

print.zppr <- function(self)
{
	# print a basic summary of the ppr fit
	alpha <- as.matrix(self$alpha)
	betas <- self$beta
	
	zerr <- residuals(self$fit)
	zhat <- predict(self$fit)
	
	sigma <- sqrt(self$sig2)
	r2 <- 1-var(zerr)/var(self$response)
	
	cat('\n#####################\n')
	cat('Summary for PPR model\n')
	cat('\nAlpha:\n')
	print(round(alpha,4))
	
	cat('\nBeta:\n')
	print(round(betas,4))
	
	cat('\nSD: \n',round(sigma,4) ,'\n')
	cat('\nR2: \n',round(r2,4) ,'\n')	
}

plot.zppr <- function(self)
{
	# plot basic graphics of the residuals
	zerr <- residuals(self$fit)
	zhat <- predict(self$fit)
	
	qqnorm(scale(zerr), ylab = 'Residuals Quantiles', main = '')
	abline(0,1)	
	
	plot(zhat,zerr/sd(zerr), ylab = 'Standard Residuals', 
		xlab = 'fitted')
	abline(h = 0, lty =3)

}

alpha2mat <- function(alpha,nterms)
{
	# transform a vector of alphas in a matrix with unitary columns
	# used by ppr2gam
	ans <- matrix(alpha, ncol = nterms)

	for(i in seq(nterms))
		ans[,i] <- ans[,i]/sqrt(sum(ans[,i]^2))
		
	return(ans)
}

ppr_alpha <- function(alpha,arg, opti)
{
	#optimize a gam model for given alphas
	# used by ppr2gam
	amat <- alpha2mat(alpha,arg$nterms)
	nu <- arg$x %*% amat
	yx <- data.frame(arg$y,nu)
	names(yx) <- c('Y',paste('NU',seq(arg$nterms),sep=''))
	
	if(is.null(arg$w))
		fit <- gam(as.formula(arg$cmd), data = yx)
	else
		fit <- gam(as.formula(arg$cmd), data = yx, weights = arg$w)
	
	# opti == TRUE is used when optimizing alpha
	if(opti)
		return(fit$gcv.ubre)
	else
		return(fit)
}

ppr2gam <- function(self, k = 5, basis = 'cr', m = 2, fx = FALSE, 
	gls_tol = 1e-6, gls_maxit = 10,...)
{
	# create string for the model formula
	cmd <- "Y~"
	sopt <- paste(",bs='", basis, "',fx=",fx,",m=",m,")",sep='')

	for(i in seq(self$nterms))
	{
		if(i == 1)
			cmd <- paste(cmd,"s(NU1, k=",k,	sopt, sep='') 
		else
			cmd <- paste(cmd,"+s(NU", i,", k=",k,sopt,sep='')
	}

	w <- self$w
	a <- as.vector(self$alpha)

	#Update the initial fitting of the PPR with the gam settings
	if(self$criteria == 'gls')
	{
		sig2 <- self$sig2
		sig2_old <- Inf
		iter <- 0
	
		while(abs(sig2-sig2_old) > gls_tol)
		{
			arg <- list(nterms = self$nterms , y = self$response, 
				x = as.matrix(self$x), w = w, cmd = cmd)

			sol <- optim(a, ppr_alpha, arg = arg, opti = TRUE, ...)
			a <- sol$par
	
			# Fit a gam object with the optimal alphas
			fit <- ppr_alpha(a, arg = arg, opti = FALSE,...)
			
			sig2_old <- sig2 	
			sig2 <- fit$sig2
			w <- 1 / (1 + self$sample_var/sig2)
			iter <- iter + 1
		
			if(iter > gls_maxit)
			{
				cat('\nWarning GLS maximum iterations reached\n')
				break
			}
		} # gls while
	}
	else
	{
		arg <- list(nterms = self$nterms , y = self$response, 
				x = as.matrix(self$x), w = w, cmd = cmd)
		
		sol <- optim(a, ppr_alpha, arg = arg, opti = TRUE, ...)
		a <- sol$par
				
		fit <- ppr_alpha(a, arg = arg, opti = FALSE, ...)
	}
	
	amat <- as.data.frame(alpha2mat(a,arg$nterms))
	rownames(amat) <- rownames(amat)
	colnames(amat) <- colnames(amat)
	
	ans <- list(gam = fit, w = w, alpha = amat)
			
	class(ans) <- append('list','gppr')
	
	return(ans)
}

terms_plot <- function(x,...) UseMethod('terms_plot',x)
terms_plot.gppr <- function(self, show_resid = TRUE, pch = 1, ...)
{
	plot(self$gam, residuals=show_resid, pch = 1,...)
}

plot.gppr <- function(self,...)
{
	serr <- scale(residuals(self$gam))
	qqnorm(serr, xlab = 'Normal quantiles',
		ylab = 'Residuals quantiles', main = '')
	abline(0,1)
	
	plot(predict(self$gam),serr,
		xlab = "Predictions", ylab = 'Standard residuals' )
	abline(h=0)
}

print.gppr <- function(self)
{
	print(summary(self$gam))
	cat('\nAlpha:\n')
	print(self$alpha)
}

residuals.gppr <- function(self,...)
{
	residuals.gam(self$gam,...)
}

predict.gppr <- function(self,...)
{
	predict.gam(self$gam,...)
}
