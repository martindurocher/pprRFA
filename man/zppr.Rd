\name{zppr}
\alias{zppr}
\alias{new_rfa_data}
\alias{set_response}

\title{
Wrapper for Projection Pursuit Regression with RFA data
}

\description{
Wrapping function for the ppr methods in stats packages. Mainly used as initial 
steps for the ppr2gam method 
}

\usage{

# To create a rfa data object
new_rfa_data <- function(annual,attrib)

# To add or change the hydrologial variables to fit
set_response.rfa <- function(self,z)

# To fit a PPR model
zppr(self, w=NULL, criteria = 'wls', nterms = 1, ...)
}

\arguments{

\item{annual}{Dataset of the hydrological variable, First row must be a 
  identification number (or string) for a site. Second row must be 
  a year. Third row must be the annual maximum values.
}

\item{attrib}{Dataset of the site characteristics. The first row must be
  the same identification number as in annual. 
}  

\item{z}{
Dataset of the hydrological variables to predict (response). The first row
is the site identification, the second row is the response variable 
and the third row is its sampling variance.
}

\item{self}{
 Data object, output from new_rfa_data 
}

\item{w}{Weights for the sites. Default is the length of the time series
calculated from annual}

\item{criteria}{Least-squares criteria used to fit the PPR model.
  Choice are : Ordinary = 'ols', "weigthed = " wls, generalized = 'gls'} 

\item{nterms}{Number of terms in the PPR model} 

\item{\dots}{Additional parameter for the ppr method. See the documentation
of respective function for more delails }
}

\details{

A projection pursuit regression (PPR) is a model of the form

$ mu + sum_p left[ g_p (alpha_p X) right] $

where $alpha_p$ are direction vectors, X is a design matrix of covariates
and $g_p$ are smooth functions.

An iterative procedure is used to estimate jointly $mu$, $alpha_p$ and
smooth functions $g_p$ by least-squares. See  Friedman &
Stuetzle (1981)

If generalized least-squares are used (criteria == gls), the total error 
is characterized by a sampling error and a model error. 
See Tasker & Stedinger (1989). The resulting covariance matrix isgit

$sigma^2 I + Sigma$

where $sigma^2$ is the variance of the model error and $Sigma$ is 
the covariance matrix (diagonal) of the sampling error. First a PPR 
model with is fitted using ordinary least-squares (criteria == ols) to 
estimate the $sigma^2$ by the empirical residual variance. Afterwards, 
the gls weights are calculated. This process is not iterated. An improved 
solution is later obtained by the ppr2gam method. 

}
\value{
%% put
}
\references{

  Durocher, M., Chebana, F., & Ouarda, T. B. M. J. (2015). 
  A Nonlinear Approach to Regional Flood Frequency Analysis Using 
  Projection Pursuit Regression. Journal of Hydrometeorology (In press). 
  doi:10.1175/JHM-D-14-0227.1 

  Friedman, J. H., & Stuetzle, W. (1981). Projection Pursuit Regression. 
  Journal of the American Statistical Association, 76, 817–823. 
  doi:10.1080/01621459.1981.10477729

  Hwang, J.-N., Lay, S.-R., Maechler, M., Martin, R. D., & Schimert, J. 
  (1994). Regression modeling in back-propagation and projection 
  pursuit learning. Neural Networks, IEEE Transactions on, 5(3), 
  342–353. doi:10.1109/72.286906

  Tasker, G., & Stedinger, J. (1989). An operational GLS model for 
  hydrologic regression. Journal of Hydrology, 111(1), 361–375. 
  doi:10.1016/0022-1694(89)90268-0

}
\author{
Martin Durocher <martin.durocher@ete.inrs.ca>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{ppr}, \link{ppr2gam}
}

\examples{

#See ppr2gam method

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
