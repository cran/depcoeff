% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/packdepcoeff.R
\name{spearr}
\alias{spearr}
\title{Spearman regression coefficient}
\usage{
spearr(x,y,direction=NULL,out=0)
}
\arguments{
\item{x}{data matrix of regressor variables}

\item{y}{data vector of the target variable}

\item{direction}{vector of length d (d is number of regressors),
value 1 refers to regressors leading to increasing y whenever this regressor increases,
value -1 refers to regressors leading to decreasing y whenever this regressor increases.
If direction=NULL, then all coefficients are computed.}

\item{out}{value 1: full output, value 0: reduced output, only coefficients
that are largest in absolute value}
}
\value{
list of Spearman regression coefficients for several directions
}
\description{
The function spearr evaluates the multivariate Spearman regression coefficient.
It describes how well the target variable y can be fit by a function of regressor variables
which is increasing w.r.t. some regressors and decreasing w.r.t. the other
regressors.
}
\examples{
library(MASS)
data <- gilgais
spearr(data[,1:3],data[,4],out=1)
}
\references{
Eckhard Liebscher (2019). A copula-based dependence measure for regression analysis. submitted
}
