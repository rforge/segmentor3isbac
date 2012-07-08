\name{Segmentor}
\alias{Segmentor}
\alias{print.Segmentor}
\alias{Segmentor.default}
\title{Segmentor}
\description{
  The functions are used for change-point problems. Given a loss function (Poisson, Normal homoscedastic, Negative Binomial or Normal Heteroscedastic (with given constant mean)), the function Segmentor estimates the optimal segmentation with respect to the log-likelihood. The Segmentor gives estimates of the breakpoint locations as well as the loss function parameter of each segment.
}
\usage{
Segmentor(data=numeric(), file=character(), model=1, Kmax = 15, param = 0.15)
}
\arguments{

\item{data}{A vector of observations to be segmented. Must have no missing values.}
\item{file}{Instead of data, the location of a file containing the data to segment. This file must have one numeric value per line and have no missing values. The filename MUST consist of the whole path to the file.}
\item{model}{Emission distribution of observed data, 1: poisson (default), 2: normal-homoscedastic, 3: negative binomial or 4: normal-heteroscedastic}
\item{Kmax}{The maximum number of segments wanted for the data. The Segmentor will find all optimal segmentations in 1 to Kmax segments.}
\item{param}{Needed only for the Negative Binomial and the Normal Heteroscedastic distributions: for model 3, the inverse of the overdispersion parameter. By default set to 0.15 as a heuristically good approximate of the parameter for NGS data. For model 4, the mean of the signal. If not entered, the function uses the empirical mean of the data.}
}

\value{
\item{model}{Emission distribution (Poisson, Normal Homoscedastic, Negative binomial or Normal Heteroscedastic)} 
\item{breaks}{Matrix of size Kmax*Kmax of estimated change-point locations for each optimal segmentation in 1 to Kmax segments.}
\item{parameters}{Matrix of size Kmax*Kmax which elements are the estimated parameters for each segment of the optimal segmentation. If model is Poisson or Normal, the parameter corresponds to the mean of the signal in each segment. If model is Negative binomial, the parameter corresponds to the success-probability of the signal in each segment. If model is normal heteroscedastic, the parameter is the variance assuming known mean.} 
\item{likelihood}{Vector of size Kmax of resulting negative log-likelihood for each optimal segmentation.}
}



\references{
  Pruned dynamic programming for optimal multiple change-point detection: \url{http://arxiv.org/abs/1004.0887}

  A Generic Implementation of the Pruned Dynamic Programing Algorithm: \url{http://arxiv.org/abs/1204.5564}
}

\details{
  \tabular{ll}{
    Package: \tab Segmentor3IsBack\cr
    Type: \tab Package\cr
    Version: \tab 1.1\cr
    Date: \tab 2012-05-10\cr
    License: \tab GPL (>= 2)\cr
  }
}

\author{
Alice Cleynen, Michel Koskas and Guillem Rigaill

Maintainer: Who to complain to <alice.cleynen@agroparistech.fr>
}

\keyword{ R-package }
\keyword{ PDPA }
\keyword{ Segmentation }


\examples{
require(Segmentor3IsBack);
N=2000; x=c(rpois(N,2.0),rpois(2*N,2.2),rpois(N,1.9));
res=Segmentor(data=x,Kmax=3);  # Finds the optimal segmentation in 1, 2 and 3 segments with respect to the Poisson model

## alternate method, use file containing data 
N=1500; y=c(rnorm(N,0,1),rnorm(2*N,2.2,1),rnorm(N,1.7,1),rnorm(N,0.3,1));
write.table(y,file="mytest.txt",col.names=FALSE, row.names=FALSE)
##res1=Segmentor(file="~/mytest.txt",model=2,Kmax=5); 

}


