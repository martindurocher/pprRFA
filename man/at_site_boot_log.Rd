\name{at_site_boot_log}
\alias{at_site_boot_log}

\title{
Perform frequency analysis
}
\description{
Perform at-site frequency analysis using maximum likelihood and 
generalized extreme value distributions. The function estimate directly 
the desired return level at the logarithm scale with its variance using 
bootstrap. 
}

\usage{
at_site_boot_log(self, ret = 2, nboot=500, ...)
}

\arguments{

  \item{self}{Data object, output from new_rfa_data.}
  
  \item{ret}{Return periods in year.}
  
  \item{nboot}{Number of bootstrap to perform.}
  
  \item{\dots}{More setting for the underlying fgev method. See 
  documentation in respective evd packages).}
}

\details{

}

\value{

\item{z}{Dataset with three rows: identification number, return level, 
variance.}
  

}
\references{

Coles, S. (2001). An introduction to statistical modeling of extreme 
values. Springer Verlag.

Efron, B., & Tibshirani, R. (1993). An introduction to the bootstrap 
(Vol. 57). Chapman & Hall/CRC.

}
\author{

}
\note{

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \link{fgev}
}
\examples{

# see ppr2gam

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
