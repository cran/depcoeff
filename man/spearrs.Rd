% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/packdepcoeff.R
\name{spearrs}
\alias{spearrs}
\title{Spearman regression coefficient for split domains}
\usage{
spearrs(x,y,splitp=NULL)
}
\arguments{
\item{x}{datamatrix of regressor variables with two columns,}

\item{y}{data vector of the target variable}

\item{splitp}{vector of length 2 of the splitting points,
If p1 is the first component of this vector, then the point splits the domain of the
first regressor into a left region of fraction p1 of data items and a right region
of the remaining data items. The same is done for the second regressor. As the
result we obtain 4 subregions of the regressor domain. default=c(0.5,0.5)}
}
\value{
list of Kendall regression coefficients for the 4 split regions
and the total coefficient together with the corresponding optimal directions.
  direction ++ means that y increases whenever both regressors increases
  direction +- means that y increases whenever the first regressor increases and the
other regressor decreases..etc.
}
\description{
The function spearrs evaluates the multivariate Spearman regression coefficient
for two regressors and split regressor region.
It describes how well the target variable can be fit in each split region
by a function which is increasing w.r.t. some regressors and decreasing
w.r.t. the other regressors.
}
\examples{
library(MASS)
data<- gilgais
spearrs(data[,1:2],data[,3],splitp=c(0.4,0.6))
}
\references{
Eckhard Liebscher (2019). A copula-based dependence measure for regression analysis. submitted
}
