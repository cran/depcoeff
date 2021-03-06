% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/packdepcoeff.R
\name{zetapm}
\alias{zetapm}
\title{Zeta dependence coefficient of piecewise monotonicity}
\usage{
zetapm(x,y,amin=0.25,method="all",methodF=1,parp=1.5,parH=0.5)
}
\arguments{
\item{x, y}{data vectors of the two variables whose dependence is analysed.}

\item{amin}{minimum fraction of sample items to be used for one split region}

\item{method}{vector of chosen special coefficients:
Spearman...Spearman coefficient
footrule...Spearman's footrule
power...power coefficient
Huber...Huber function coefficient,
"all" refers to all coefficients}

\item{methodF}{value 1,2 or 3 refers to several methods for computation of the
distribution function values, 1 is the default value.}

\item{parp}{parameter of the power function (default 1.5). The parameter has
to be positive.}

\item{parH}{parameter of the Huber function (default 0.5). Valid values for parH
are between 0 and 1.}
}
\value{
list of zeta dependence coefficients (plusminus coefficient and minusplus one)
of piecewise monotonicity of two random variables containing the
following elements or a subset of it in this order:
Spearman coefficient, footrule, power coefficient, Huber function coefficient.
position1 and position2 indicate the number of the sample items where the optimized
split point is located
}
\description{
zetapm is a function to evaluate the zeta dependence coefficients of piecewise
monotonicity of two random variables x and y which is based on the copula. The regressor
domain (domain of x) is split into two parts. The function searches
for the optimal splitting point to obtain maximum
depedence. The main part of the function is coded as C++ procedure
}
\details{
Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
for the estimators of values \eqn{F(X_{i})} of the distribution function:
  methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
  methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
  methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
The values of the distribution function of \eqn{Y} are treated analogously.
}
\examples{
library(MASS)
data<- gilgais
zetapm(data[,1],data[,2])
}
\references{
Eckhard Liebscher (2017). Copula-based dependence measures for piecewise
monotonicity. Dependence Modeling 5 (2017), 198-220
}
