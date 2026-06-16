### R package depcoeff

#' @import copula

#### function psi
psi<- function(x,id,par=-1){
  h<- abs(x)
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if (id==1){return(x^2)}    # Spearman coefficient
  if (id==2){return(h)}      # Spearman's footrule
  if (id==3){return(h^par)}  # power coefficient
  if (id==4){                # Huber function
    h0<- ifelse(h<par,0.5*h*h,par*(h-0.5*par))
    return(h0)}
}
#### derivative of function psi
psider<- function(x,id,par=-1){
  h<- abs(x)
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if (id==1){return(2.0*x)}   # Spearman-coefficient
  if (id==2){return(sign(x))}  # Spearman's footrule
  if (id==3){return(sign(x)*par*h^(par-1))}  # power coefficient
  if (id==4){                # Huber function
    h0<- ifelse(h<par,x,par*sign(x))
    return(h0)}
}
#### normalisation factor psi_bar
psif<- function(id,par=-1){
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if ((id==4)&(par>=1.0)) {par<- 0.5}
  if (id==1){return(0.166666666667)}
  if (id==2){return(0.333333333333)}
  if (id==3){return(2.0/(par*(par+3.0)+2.0))}
  if (id==4){return(par*(2.0-par)*(par*(par-2.0)+2.0)/12.0)}
}
#### normalisation factor psi_breve for xi coefficient
psif2<- function(id,par=-1){
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if ((id==4)&(par>=1.0)) {par<- 0.5}
  if (id==3){return(1.0/(par+1.0))}
  if (id==4){return(par*(par*(par-3.0)+3.0)/6.0)}
}
### minimal value of zeta for the specific coefficients
zetamin<- function(id,par=-1){
  if (par<=0) {
    if (id==3) {par<- 1.5}
    if (id==4) {par<- 0.5}}
  if ((id==4)&(par>=1.0)) {par<- 0.5}
  if (id==1){return(-1.0)}
  if (id==2){return(-0.5)}
  if (id==3){return(-par/2.0)}
  if (id==4){return((par^3-2.0*par^2+2.0)/(par^3-4.0*par^2+6.0*par-4.0))}
}

#' Zeta dependence coefficient
#'
#' zetac is a function to evaluate the zeta dependence coefficient (one interval)
#' of two random variables x and y which is based on the copula. Four specific coefficients
#' are available: the Spearman coefficient, Spearman's footrule, the power coefficient
#' and the Huber function coefficient.
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage zetac(x,y,method="Spearman",methodF=1,parH=0.5,parp=1.5,outsd=TRUE,eps=0.95)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param method list of names of the coefficients: "Spearman" stands for the
#' Spearman coefficient, "footrule" means Spearman's footrule, "power" stands
#' for the power function coefficient, "Huber" means the Huber function
#' coefficient. If "all" is assigned to method then all methods are used.
#' @param methodF value 1,2 or 3 refers to several methods for computation of
#' the distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for
#' parH are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @param outsd logical. If TRUE, the estimated standard deviation and confidence intervals
#' are evaluated.  If FALSE, the evaluation is suppressed
#' @param eps real value, confidence level.
#' @return A list of the following vectors:
#'   \describe{
#'     \item{Spearman}{zeta Spearman coefficient, standard error, confidence interval}
#'     \item{footrule}{zeta footrule coefficient, standard error, confidence interval}
#'     \item{power}{zeta coefficient with power function, standard error, confidence interval}
#'     \item{Huber}{zeta coefficient with Huber function, standard error, confidence interval}
#'   }
#' The zeta dependence coefficient of two random variables is bounded
#' by 1. The higher the value the stronger is the dependence.
#' @references Eckhard Liebscher (2014). Copula-based dependence measures. Dependence Modeling 2 (2014), 49-64
#' @importFrom stats var qnorm
#' @examples
#' library(MASS)
#' data<- gilgais
#' zetac(data[,1],data[,2],method="all")
#' @export

zetac<- function(x,y,method="Spearman",methodF=1,parH=0.5,parp=1.5,outsd=TRUE,eps=0.95){
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")} #check the parameters
  if (n<4) {stop("not enough sample items")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if ((eps<=0.5)|(eps>=1)){ stop("wrong epsilon")}
  if (!methodF %in% 1:3) methodF<- 1

    ### compute normed difference of ranks
  if (methodF==1) {
    u<- rank(x)/n
    v<- rank(y)/n
  }else{
    if (methodF==2){
      u<- rank(x)/(n+1)
      v<- rank(y)/(n+1)
    } else {
      u<- rank(x)/sqrt(n^2-1)
      v<- rank(y)/sqrt(n^2-1)
  }}
  uu<- u-v

  if (method=="all") {method<- c("Spearman","footrule","power","Huber")}

  if (outsd){ z<- qnorm((1+eps)/2)}   #quantile
  ### Spearman coefficient
  if ("Spearman" %in% method) {
    h<- psi(uu,1)
    s<- 1.0-sum(h)/(n*psif(1))  #coefficient
    if (outsd){
      h2<- rep(0.0,n)
      for (j in 1:n){
        h2<- h2+psider(u[j]-v[j],1)*((u<=u[j])-u[j]-(v<=v[j])+v[j])
      }
      sd<- sqrt(var(h+(h2/n))/n)/psif(1)     #standard deviation
      outl<- list(Spearman=s,Spearman_se=sd,Spearman_conf=c(s-sd*z,s+sd*z))
    }else{
      outl<- list(Spearman=s)
    }
  } else {outl<- list()}
  ### Spearman's footrule
  if ("footrule" %in% method) {
    h<- psi(uu,2)
    s<- 1.0-sum(h)/(n*psif(2))  #coefficient
    if (outsd){
      h2<- rep(0.0,n)
      for (j in 1:n){
        h2<- h2+psider(u[j]-v[j],2)*((u<=u[j])-u[j]-(v<=v[j])+v[j])
      }
      sd<- sqrt(var(h+(h2/n))/n)/psif(2)     #standard deviation
      outl<- append(outl,list(footrule=s,footrule_se=sd,footrule_conf=c(s-sd*z,s+sd*z)))
    }else{
      outl<- append(outl,list(footrule=s))
  }}
  ### power function coefficient
  if ("power" %in% method) {
    h<- psi(uu,3,parp)
    s<- 1.0-sum(h)/(n*psif(3,parp))  #coefficient
    if (outsd){
      h2<- rep(0.0,n)
      for (j in 1:n){
        h2<- h2+psider(u[j]-v[j],3,parp)*((u<=u[j])-u[j]-(v<=v[j])+v[j])
      }
      sd<- sqrt(var(h+(h2/n))/n)/psif(3,parp)     #standard deviation
      outl<- append(outl,list(power=s,power_se=sd,power_conf=c(s-sd*z,s+sd*z)))
    }else{
      outl<- append(outl,list(power=s))
  }}
  ### Huber function coefficient
  if ("Huber" %in% method) {
    h<- psi(uu,4,parH)
    s<- 1.0-sum(h)/(n*psif(4,parH))  #coefficient
    if (outsd){
      h2<- rep(0.0,n)
      for (j in 1:n){
        h2<- h2+psider(u[j]-v[j],4,parH)*((u<=u[j])-u[j]-(v<=v[j])+v[j])
      }
      sd<- sqrt(var(h+(h2/n))/n)/psif(4,parH)     #standard deviation
      outl<- append(outl,list(Huber=s,Huber_se=sd,Huber_conf=c(s-sd*z,s+sd*z)))
    }else{
      outl<- append(outl,list(Huber=s))
  }}

  return(outl)  ## output: list of computed coefficients,
}

#' Zeta dependence coefficient of piecewise monotonicity
#'
#' zetapm is a function to evaluate the zeta dependence coefficients of piecewise
#' monotonicity of two random variables x and y which is based on the copula. The regressor
#' domain (domain of x) is split into two parts. The function searches
#' for the optimal splitting point to obtain maximum
#' depedence. The main part of the function is coded as C++ procedure
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage zetapm(x,y,amin=0.25,method="all",methodF=1,parp=1.5,parH=0.5)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param amin minimum fraction of sample items to be used for one split region
#' @param method vector of chosen special coefficients:
#'   Spearman...Spearman coefficient
#'   footrule...Spearman's footrule
#'   power...power coefficient
#'   Huber...Huber function coefficient,
#'   "all" refers to all coefficients
#' @param methodF value 1,2 or 3 refers to several methods for computation of the
#' distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for parH
#' are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @return A list of zeta dependence coefficients of piecewise monotonicity of two
#' random variables containing the following elements:
#'   \describe{
#'     \item{plusminuscoeff}{zeta plusminus coefficients: Spearman, footrule, power, Huber coefficients if chosen}
#'     \item{position1}{positions of the optimal split points within the x-value (zeta plusminus coefficients)}
#'     \item{minuspluscoeff}{zeta minuspluscoeff coefficients: Spearman, footrule, power, Huber coefficients if chosen}
#'     \item{position2}{positions of the optimal split points within the x-value (zeta minusplus coefficients)}
#'     }
#' @references Eckhard Liebscher (2017). Copula-based dependence measures for piecewise
#' monotonicity. Dependence Modeling 5 (2017), 198-220
#' @examples
#' x<- seq(0,1.0,by=0.01)
#' eps<- rnorm(length(x),sd=0.05)
#' y<- x-2*(x>=0.4)*(x-0.4)+eps
#' zetapm(x,y,amin=0.2,method="all",methodF=1,parp=1.5,parH=0.5)
#' @export

zetapm<- function(x,y,amin=0.25,method="all",methodF=1,parp=1.5,parH=0.5){
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")} #check the parameters
  if ((parp<=0.0)|(parH<=0.0)|(parH>=1.0)){ stop("parameter(s) not valid")}
  if (!methodF %in% 1:3) methodF<- 1
  if (n<8) {  stop("not enough sample items")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if ((amin>=0.5)|(amin<=0.0)) {amin<- 0.25}  ## amin--> default 0.25
  method0<- c("Spearman","footrule","power","Huber")    #vector of all methods
  if (length(method)==1) {if (method=="all") {method<- method0}}
  mv<- which(method == method0)  # vector with numbers of chosen methods

  u<- rank(x)
  u1<-u2<-v1<-v2<- vector()
  na<- max(ceiling(amin*n+0.5),4)  #number of items in region of minimum length
  ne<- n-na  #length of u1,u2
  # data for the left-side part
  u1<- u[u<=ne]  # x ranks of the left part
  v1<- rank(y[u<=ne])  # y ranks of the left part
  v3<- ne+1.0-v1   # y ranks for the minus-plus coefficient

  # data for the right-side part
  u2<- u[u>na]-na  # x ranks of the right part
  v2<- rank(y[u>na])  # y ranks of the leftright part
  v4<- ne+1.0-v2  # y ranks for the minus-plus coefficient

  ## plus-minus-coefficient
  h1<- coeffpml(u1,v1,u2,v2,amin,parp,parH,n,na,methodF)

   # data for the left-side part
  u1<- u[u<=ne]  # x ranks of the left part

  # data for the right-side part
  u2<- u[u>na]-na  # x ranks of the right part

  ### minus-plus coefficient
  h2<- coeffpml(u1,v3,u2,v4,amin,parp,parH,n,na,methodF)

  ###output of coefficients: Spearman, Spearman's footrule, power function, Huber function
  outl<- list(h1[mv],h1[(mv+4)],h2[mv],h2[(mv+4)])
  names(outl)<- c("plusminuscoeff","position1","minuspluscoeff","position2")
  return(outl)
}

#' Zeta coefficient of piecewise monotonicity with split domain
#'
#' The function zetaci evaluates the coefficient of piecewise monotonicity of variables x and y where the x-domain
#' is split into a fixed number of intervals.
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage zetaci(x,y,a,method="Spearman",methodF=1,parH=0.5,parp=1.5)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param a vector of fractions \eqn{a_{i},0<a_{i}<a_{i+1}<1} for the splitting. A fraction
#' of \eqn{a_{1},a_{2}-a_{1},a_{3}-a{2}}... of data points are in the corresponding split region.
#' The number of split regions is equal to the length of \eqn{a} plus 1.
#' @param method value (default "Spearman")
#' @param methodF value 1,2 or 3 refers to several methods for computation of
#' the distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for
#' parH are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @return A list of zeta dependence coefficients of piecewise monotonicity of two
#' random variables containing the following elements:
#'   \describe{
#'     \item{Spearman}{zeta Spearman plusminus/minusplus coefficient}
#'     \item{footrule}{zeta Spearman's footrule plusminus/minusplus coefficient}
#'     \item{power}{zeta power plusminus/minusplus coefficient}
#'     \item{Huber}{zeta Huber plusminus/minusplus coefficient}
#'   }
#' @references Eckhard Liebscher (2017). Copula-based dependence measures for piecewise
#' monotonicity. Dependence Modeling 5 (2017), 198-220
#' @examples
#' library(MASS)
#' x<- seq(0,1.0,by=0.01)
#' eps<- rnorm(length(x),sd=0.05)
#' y<- x-2*(x>=0.4)*(x-0.4)+4*(x>=0.75)*(x-0.75)+eps
#' zetaci(x,y, a=c(0.25, 0.5, 0.75))
#' @export


zetaci<- function(x,y,a,method="Spearman",methodF=1,parH=0.5,parp=1.5){
  if ((parp<=0.0)|(parH<=0.0)|(parH>=1.0)|(min(a)<=0.0)|(max(a)>=1.0)){
    stop("parameter(s) not valid")}  #check the parameters
  if (length(method)==1) { if (method=="all") {method<- c("Spearman","footrule","power","Huber")}}
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if (n<4) {  stop("not enough sample items")}
  if (!methodF %in% 1:3) methodF<- 1

  m<- length(a)+1
  a[m]<- 1.0
  A<- vector()   # vector of cumulative numbers of sample items with X_i in intervals I_1...I_k
  A[1]<- as.integer(floor(a[1]*n))
  tv<- (A[1]<3)
  if (m>2){
  for (j in 2:(m-1)) {
    A[j]<- as.integer(floor(a[j]*n))
    if (A[j]-A[j-1]<3) {tv<- TRUE}
  }}
  A[m]<- n
  if (A[m]-A[m-1]<3) {tv<- TRUE}   ## a minimum of 3 sample items must be in every interval
  if (tv) {
    stop("vector a not valid")}
    outl<- list()
  u<- rank(x)
  ss1<- ss2<- sf1<- sf2<- sp1<- sp2<- sh1<- sh2<- 0.0
  signi<- 1  # sign of Y in the interval
  for (k in 1:m) {
    ### compute normed difference of ranks
    if (k==1) {
      AA<- 0
      u1<- u[u<=A[1]]
      v1<- y[u<=A[1]]   # sample of Y-values with X_i in I_k
    } else {
      AA<- A[k-1]
      u1<- u[(u>A[k-1])&(u<=A[k])]
      v1<- signi*y[(u>A[k-1])&(u<=A[k])]  # sample of Y-or (-Y)-values with X_i in I_k
    }
    v2<- rank(v1)
    if (methodF==1) {
      uu<- (u1-AA-v2)/(A[k]-AA)  # normed difference of the ranks
      vv<- (u1-A[k]-1+v2)/(A[k]-AA)
    }else{
      if (methodF==2){
        uu<- (u1-AA-v2)/(A[k]-AA+1.0)  # normed difference of the ranks
        vv<- (u1-A[k]-1+v2)/(A[k]-AA+1.0)
      } else {
        uu<- (u1-AA-v2)/sqrt((A[k]-AA)^2-1.0)  # normed difference of the ranks
        vv<- (u1-A[k]-1+v2)/sqrt((A[k]-AA)^2-1.0)
      }}
    ### Spearman coefficient
    if ("Spearman" %in% method) {
      ss1<- ss1+sum(psi(uu,1))
      ss2<- ss2+sum(psi(vv,1))}
    ### Spearman's footrule
    if ("footrule" %in% method) {
      sf1<- sf1+sum(psi(uu,2))
      sf2<- sf2+sum(psi(vv,2))}
    ### power function coefficient
    if ("power" %in% method) {
      sp1<- sp1+sum(psi(uu,3,parp))
      sp2<- sp2+sum(psi(vv,3,parp))}
    ### Huber function coefficient
    if ("Huber" %in% method) {
      for (i in 1:(length(uu))){
        sh1<- sh1+ psi(uu[i],4,parH)
        sh2<- sh2+ psi(vv[i],4,parH)
      }
    }
    signi<- -signi
  }
  ## evaluation of the coefficients
  if ("Spearman" %in% method) {
    ss1<- 1.0-(ss1/n/psif(1))
    ss2<- 1.0-(ss2/n/psif(1))
    outl<- append(outl,list(Spearman=c(ss1,ss2)))}
  if ("footrule" %in% method) {
    sf1<- 1.0-(sf1/n/psif(2))
    sf2<- 1.0-(sf2/n/psif(2))
    outl<- append(outl,list(footrule=c(sf1,sf2)))}
  if ("power" %in% method) {
    sp1<- 1.0-(sp1/n/psif(3,parp))
    sp2<- 1.0-(sp2/n/psif(3,parp))
    outl<- append(outl,list(power=c(sp1,sp2)))}
  if ("Huber" %in% method) {
    sh1<- 1.0-(sh1/n/psif(4,parH))
    sh2<- 1.0-(sh2/n/psif(4,parH))
    outl<- append(outl,list(Huber=c(sh1,sh2)))}
  return(outl)
}

#' Xi dependence coefficient
#'
#' xic is a function to evaluate the xi dependence coefficient (one interval)
#' of two random variables x and y which is based on the copula. Two specific coefficients
#' are available: the power coefficient and the Huber function coefficient.
#'
#' Let \eqn{X_{1},\ldots ,X_{n}} be the sample of the \eqn{X} variable. Formulas
#' for the estimators of values \eqn{F(X_{i})} of the distribution function:
#'   methodF = 1 \eqn{\rightarrow \hat{F}(X_{i})=\frac{1}{n}\textrm{rank}(X_{i})}
#'   methodF = 2 \eqn{\rightarrow \hat{F}^{1}(X_{i})=\frac{1}{n+1}\textrm{rank}(X_{i})}
#'   methodF = 3 \eqn{\rightarrow \hat{F}^{2}(X_{i})=\frac{1}{\sqrt{n^{2}-1}}\textrm{rank}(X_{i})}
#' The values of the distribution function of \eqn{Y} are treated analogously.
#' @usage xic(x,y,method="power",methodF=1,parH=0.5,parp=1.5,outsd=TRUE,eps=0.95)
#' @param x,y data vectors of the two variables whose dependence is analysed.
#' @param method list of names of the coefficients: "power" stands
#' for the power function coefficient, "Huber" means the Huber function
#' coefficient. If "all" is assigned to method then all methods are used.
#' @param methodF value 1,2 or 3 refers to several methods for computation of
#' the distribution function values, 1 is the default value.
#' @param parH parameter of the Huber function (default 0.5). Valid values for
#' parH are between 0 and 1.
#' @param parp parameter of the power function (default 1.5). The parameter has
#' to be positive.
#' @param outsd logical. If TRUE, the estimated standard deviation and confidence intervals
#' are evaluated.  If FALSE, the evaluation is suppressed
#' @param eps real value, confidence level.
#' @return A list of the following vectors:
#'   \describe{
#'     \item{power}{zeta coefficient with power function, standard error, confidence interval}
#'     \item{Huber}{zeta coefficient with Huber function, standard error, confidence interval}
#'   }
#' The zeta dependence coefficient of two random variables is bounded
#' by 1. The higher the value the stronger is the dependence.
#' @references Eckhard Liebscher (2014). Copula-based dependence measures. Dependence Modeling 2 (2014), 49-64
#' @importFrom stats var qnorm
#' @examples
#' library(MASS)
#' data<- gilgais
#' xic(data[,1],data[,2])
#' @export

xic<- function(x,y,method="power",methodF=1,parH=0.5,parp=1.5,outsd=TRUE,eps=0.95){
  n<- length(x)
  if (n!=length(y)){ stop("number of sample items in x and y are different")} #check the parameters
  if (n<4) {stop("not enough sample items")}
  if ((!is.vector(x))|(!is.vector(y))){ stop("data format wrong")}
  if ((eps<=0.5)|(eps>=1)){ stop("wrong epsilon")}
  if (!methodF %in% 1:3) methodF<- 1

  ### compute normed difference of ranks
  if (methodF==1) {
    u<- rank(x)/n
    v<- rank(y)/n
  }else{
    if (methodF==2){
      u<- rank(x)/(n+1)
      v<- rank(y)/(n+1)
    } else {
      u<- rank(x)/sqrt(n^2-1)
      v<- rank(y)/sqrt(n^2-1)
    }}
  uu<- u-v
  uuu<- u+v-1

  if (method=="all") {method<- c("power","Huber")}

  if (outsd){ z<- qnorm((1+eps)/2)}   #quantile
  ### power function coefficient
  if ("power" %in% method) {
    h<- psi(uuu,3,parp)-psi(uu,3,parp)
    s<- min(1.0,sum(h)/(n*psif2(3,parp)))  #coefficient
    if (outsd){
      h2<- rep(0.0,n)
      for (j in 1:n){
        h2<- h2+(psider(u[j]+v[j]-1,3,parp)*((u<=u[j])-u[j]+(v<=v[j])-v[j])
                 -psider(u[j]-v[j],3,parp)*((u<=u[j])-u[j]-(v<=v[j])+v[j]))
      }
      sd<- sqrt(var(h+(h2/n))/n)/psif2(3,parp)     #standard deviation
      outl<- list(power=s,power_se=sd,power_conf=c(s-sd*z,s+sd*z))
    }else{
      outl<- list(power=s)
    }
  } else {outl<- list()}
  ### Huber function coefficient
  if ("Huber" %in% method) {
    h<- psi(uuu,4,parH)-psi(uu,4,parH)
    s<- min(1.0,sum(h)/(n*psif2(4,parH)))  #coefficient
    if (outsd){
      h2<- rep(0.0,n)
      for (j in 1:n){
        h2<- h2+(psider(u[j]+v[j]-1,4,parH)*((u<=u[j])-u[j]+(v<=v[j])-v[j])
                 -psider(u[j]-v[j],4,parH)*((u<=u[j])-u[j]-(v<=v[j])+v[j]))
      }
      sd<- sqrt(var(h+(h2/n))/n)/psif2(4,parH)     #standard deviation
      outl<- append(outl,list(Huber=s,Huber_se=sd,Huber_conf=c(s-sd*z,s+sd*z)))
    }else{
      outl<- append(outl,list(Huber=s))
    }}

  return(outl)  ## output: list of computed coefficients,
}

#' Spearman regression coefficient
#'
#' The function spearr evaluates the multivariate Spearman regression coefficient.
#' It describes how well the target variable y can be fit by a function of regressor variables
#' which is increasing w.r.t. some regressors and decreasing w.r.t. the other
#' regressors.
#'
#' @usage spearr(x,y,direction=NULL,out=0)
#' @param x data matrix of regressor variables
#' @param y data vector of the target variable
#' @param out value 1: full output, value 0: reduced output, only coefficient
#' that is largest in absolute value
#' @param direction vector of length d (d is number of regressors),
#' value 1 refers to regressors leading to increasing y whenever this regressor increases,
#' value -1 refers to regressors leading to decreasing y whenever this regressor increases.
#' If direction=NULL, then all coefficients are computed.
#' @return A list containing
#'   \describe{
#'     \item{dcoeff}{Spearman regression coefficient}
#'     \item{dir}{direction vector}
#'     }
#' @references Eckhard Liebscher (2021). On a multivariate version of Spearman's
#' correlation coefficient for regression: Properties and Applications. Asian Journal
#' of Statistical Sciences 1, No. 2, 123-150.
#' @examples
#' library(MASS)
#' data <- gilgais
#' spearr(data[,1:3],data[,4],out=1)
#' @export

spearr<-function(x,y,direction=NULL,out=0){
  n<- length(y)
  x<- as.matrix(x)
  d<- ncol(x)

  if (n!=nrow(x)){ stop("number of sample items in x and y are different")}
  if ((!is.vector(y))|(d<2)){ stop("data format wrong")}
  if (n<4) {  stop("not enough sample items")}

  mxc<- -10.0 #maximum coefficient
  if (is.null(direction)) {
    ka<-(2^(d-1)-1)
  } else {
    if (any(direction==0)|(length(direction)!=d)){ stop("parameter direction not valid")}
    ka<- 0
    out<- 0
  }

  ret<- list()
  x<- cbind(x,-x)
  u<- apply(x,2,function(cc) (rank(cc)-0.5)/n)  #columnwise copula-transform
  for (k in ka:0){
    if (ka>0){  #all directions
      direction<- 2*as.numeric(c(1,intToBits(k)[1:(d-1)]))-1   #vector of signs, first element always plus 1
      if (k<ka)  {
        sel<- (direction<0)
        im<- (1:d) + d * sel   #select the columns of u
        uu<- u[,im]
      } else {
        uu<- u[,1:d]
      }
    } else {  # one direction  case ka=0
      sel<- (direction<0)
      im<- (1:d) + d * sel    #select the columns of u
      uu<- u[,im]
    }

    if (d>2) {z<- apply(uu,1,function(rw) (-prod(1.0-rw)+prod(rw))) #negative phi-function (rowwise)
    } else{ z<- apply(uu,1,function(rw) sum(rw)-1.0)}

    v<- 2.0*rank(y)-1.0-n  # transform the y-values
    Bn<- sum(z*v)                       #numerator in the estimator
    An<- sum(sort(z)*(2*(1:n)-1.0-n)) #denominator in the estimator
    if (out==1){   # full output
      if (An==0.0) {ret[[k+1]]<- list(dcoeff="not defined",dir=direction)
      }else{ret[[k+1]]<- list(dcoeff=Bn/An,dir=direction)}  #An,Bn,
    }else{
      if (An!=0.0){  ##reduced output
        h<- Bn/An
        if (abs(h)>mxc) {
          mxc<- abs(h)
          ret<- list(dcoeff=h,dir=direction)
        }
      }
    }
  }
  return(ret)
}

#' Spearman regression coefficient for split domains
#'
#' The function spearrs evaluates the multivariate Spearman regression coefficient
#' for two regressors and split regressor region.
#' It describes how well the response variable can be fitted in each split region
#' by a function which increases or decreases as the regressors increase.
#'
#' @usage spearrs(x,y,splitp=NULL)
#' @param x data matrix of regressor variables with two columns,
#' @param y data vector of the target variable
#' @param splitp vector of length 2 of the splitting points,
#' If p1 is the first component of this vector, then the point splits the domain of the
#' first regressor into a left region of fraction p1 of data items and a right region
#' of the remaining data items. The same is done for the second regressor. As the
#' result we obtain 4 subregions of the regressor domain. default=c(0.5,0.5)
#' @return A list of Spearman regression coefficients for the 4 split regions
#' (11=left lower region, 12=left upper region...) with components:
#'   \describe{
#'     \item{dcoeff++}{split coefficient ++}
#'     \item{dcoeff+-}{split coefficient +-}
#'     \item{totalcoeff}{total coefficient}
#'     \item{directions}{optimal directions for the several split regions}
#'   }
#' direction ++ means that y increases whenever both regressors increases
#' direction +- means that y increases whenever the first regressor increases and the
#' other regressor decreases..etc.
#' @references Eckhard Liebscher (2021). On a multivariate version of Spearman's
#' correlation coefficient for regression: Properties and Applications. Asian Journal
#' of Statistical Sciences 1, No. 2, 123-150.
#' @examples
#' library(MASS)
#' data<- gilgais
#' spearrs(data[,1:2],data[,3],splitp=c(0.4,0.6))
#' @export

spearrs<-function(x,y,splitp=NULL){
  n<- length(y)
  x<- as.matrix(x)

  if (nrow(x)!=n){ stop("number of sample items in x and y are different")} #check the parameters
  if (!is.vector(y)){ stop("data format wrong")}
  if (n<4) {  stop("not enough sample items")}

  cg<- 0.0  # total coefficient
  if (length(x[1,])>2) {x<- x[,1:2]} #first columns are selected in the case d>2

  if (is.null(splitp)) {
    splitp<- c(0.5,0.5)}  # default case: Split points in the middle
  if ((min(splitp)<=0.0)|(max(splitp)>=1.0)){ stop("parameter(s) not valid")}

  ret<- list()
  dca<- dcb<- array(0.0,dim=c(2,2,2))
  dc<- vector(length=4)
  u<- apply(x,2,function(cc) (rank(cc)-0.5)/n)  #columnwise copula transform

  us1<- (u[,1]<=splitp[1])  #logical vector for Split region 1 (left)
  us2<- (u[,2]<=splitp[2])  #logical vector for Split region 2 (right)
  for (k in 1:2){
    k0<- 1
    for (l1 in 1:2){  #Split for u[,1]
      l01<- (l1==1)
      for (l2 in 1:2){  #Split for u[,2]
        l02<- (l2==1)
        u0<- subset(u,(us1==l01)&(us2==l02)) #choose the subset
        y0<- subset(y,(us1==l01)&(us2==l02))
        n1<- length(u0[,1])  #number of sample items in the subset
        if (k==1) {   #k=1 --> "++" coefficient
          u0<- apply(u0,2,function(cc) (rank(cc)-0.5)/n1)
        } else{       #k=2 --> "+-" coefficient
          u0[,1]<- (rank(u0[,1])-0.5)/n1
          u0[,2]<- (rank(-u0[,2])-0.5)/n1
        }
        z<- apply(u0,1,function(rw) sum(rw)-1.0) #apply phi-function rowwise

        v<- 2.0*rank(y0)-1.0-n1
        Bn<- sum(z*v)  #numerator of the estimator
        An<- sum(sort(z)*(2*(1:n1)-1.0-n1)) #denominator
        if (k==1) {
          dc[k0]<- ifelse((An!=0.0),Bn/An,0.0)
        } else{
          if (An!=0.0) {h<- Bn/An}else{h<- 0.0}
          ret[[k0]]<- dc[k0]
          names(ret)[[k0]]<- paste0("dcoeff++",toString(l1),l2)
          ret[[k0+1]]<- h
          names(ret)[[k0+1]]<- paste0("dcoeff+-",toString(l1),l2)
        }
        dca[k,l1,l2]<- ifelse((An>0.0),n1*An,0.0)
        dcb[k,l1,l2]<- n1*Bn
        k0<- k0+2
      }}
  }
  mc<- 0.0
  for (t1 in 1:2){
    for (t2 in 1:2){
      for (t3 in 1:2){
        for (t4 in 1:2){
          h<- (abs(dcb[t1,1,1])+abs(dcb[t2,1,2])+abs(dcb[t3,2,1])+abs(dcb[t4,2,2]))/
            (dca[t1,1,1]+dca[t2,1,2]+dca[t3,2,1]+dca[t4,2,2])
          if (h>mc) {
            mc<- h
            if (t1==1){ ic<- ifelse((dcb[t1,1,1]>0),"++","--")
            } else { ic<- ifelse((dcb[t1,1,1]>0),"+-","-+")}
            if (t2==1){ ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"+-","-+"))}
            if (t3==1){ ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"+-","-+"))}
            if (t4==1){ ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"+-","-+"))}
            #cat(ic,"\n")
          }
        }}}}
  ret[[k0]]<- mc     #output total coefficient
  names(ret)[[k0]]<- "totalcoeff"
  ret[[k0+1]]<- ic
  names(ret)[[k0+1]]<- "directions"
  return(ret)
}

#' Kendall regression coefficient
#'
#' The function kendr evaluates the Kendall regression coefficient.
#' It describes how well the target variable y can be fitted by a function of the regressor variables
#' which increases or decreases as the regressors increase.
#'
#' @usage kendr(x,y,direction=NULL,out=0,outsd=TRUE,eps=0.95)
#' @param x data matrix of regressor variables
#' @param y data vector of the response variable
#' @param out value 1: full output, value 0: reduced output of one coefficient
#' that is largest in absolute value
#' @param direction vector of length d (d is number of regressors),
#' value 1 refers to regressors leading to increasing y whenever this regressor increases,
#' value -1 refers to regressors leading to decreasing y whenever this regressor increases.
#' If direction=NULL, then all coefficients are computed.
#' @param outsd logical. If TRUE, the estimated standard deviation and confidence intervals
#' are evaluated.  If FALSE, the evaluation is suppressed
#' @param eps real value, confidence level.
#' @return A list of the coefficients for several directions with components:
#'   \describe{
#'     \item{dcoeff}{Kendall regression coefficient}
#'     \item{dir}{direction vector}
#'     \item{Pxx}{fraction of \eqn{X \le \breve{X}} where \eqn{\breve{X}} has the same distribution as \eqn{X}}
#'     \item{sd}{standard error}
#'     \item{conf}{confidence interval}
#'   }
#' @references Eckhard Liebscher (2021). Kendall regression coefficient. Computational Statistics and Data Analysis 157 (2021). 107140
#' @importFrom stats var qnorm
#' @examples
#' library(MASS)
#' data <- gilgais
#' kendr(data[,1:3],data[,4],out=1)
#' @export

kendr<-function(x,y,direction=NULL,out=0,outsd=TRUE,eps=0.95){
  n<- length(y)
  x<- as.matrix(x)
  d<- ncol(x)

  if (nrow(x)!=n){ stop("number of sample items in x and y are different")} #check the parameters
  if ((!is.vector(y))|(d<2)){ stop("data format wrong")}
  if (n<4) {  stop("not enough sample items")}
  if ((eps<=0.5)|(eps>=1)){ stop("wrong epsilon")}

  mxc<- -10.0 #maximum coefficient
  if (is.null(direction)) {
    ka<-(2^(d-1)-1)
  } else {
    if (any(direction==0)|(length(direction)!=d)){ stop("parameter direction not valid")}
    ka<- 0
    out<- 0
  }
  ret<- list()

  x<- cbind(x,-x)
  for (k in ka:0){
    if (ka>0){  #all directions
      direction<- 2*as.integer(c(1,intToBits(k)[1:(d-1)]))-1   #vector of signs, first element always plus 1
      if (k<ka)  { # direction does not consist solely of ones
        sel<- (direction<0)
        im<- (1:d) + d * sel   #select the columns of x
        u<- x[,im]
      }else{
        u<- x[,1:d]
      }
    } else {  #  one direction
      sel<- (direction<0)
      im<- (1:d) + d * sel    #select the columns of x
      u<- x[,im]
    }
    Fuu<- F.n(u,u,offset = 0, smoothing = "none")
    v<- cbind(u,y)
    Fvv<- F.n(v,v,offset = 0, smoothing = "none")

    A<-(sum(Fuu)-1.0)  #denominator
    B<-(sum(Fvv)-1.0)  #numerator

    if (A<=0.0) {
      if (out==1){  ret[[k+1]]<- list(dcoeff="Pxx=0, coefficient not defined",dir=direction)
      } else {   stop("Pxx=0, coefficient not defined")
      }
    }else{
      kendrc<- 2*B/A-1.0
      ff<- TRUE
      if (out==1){  # output of all coefficients
        ret[[k+1]]<- list(dcoeff=kendrc,dir=direction,Pxx=A/(n-1.0))  # output Kendall coefficient
      } else {       # reduced output of the largest (absolute) value
        if (kendrc<0) {
          direction<- -direction
          kendrc<- -kendrc}
        if (kendrc>mxc) {
          mxc<- kendrc
          ret<- list(dcoeff=kendrc,dir=direction,Pxx=A/(n-1.0))
        } else { ff<-FALSE }
      }

      if ((outsd)&ff) { # evaluation of the standard deviation/confidence interval
        vv<- Fvv+F.n(-v,-v,offset = 0, smoothing = "none")-
          (Fuu+F.n(-u,-u,offset = 0, smoothing = "none"))*(1.0+kendrc)/2.0

        sd<- 2*(n-1)*sqrt(var(vv)/n)/A  #standard deviation
        z<- qnorm((1+eps)/2)
        if (out==1){
          ret[[k+1]][["sd"]]<-sd
          ret[[k+1]][["conf"]]<- c(kendrc-sd*z,kendrc+sd*z)
        }else{
          ret[["sd"]]<- sd
          ret[["conf"]]<- c(kendrc-sd*z,kendrc+sd*z)
        }
      } # end outsd
    }
  } # end k-loop
  return(ret)
}


#' Kendall regression coefficient for split domains
#'
#' The function kendrs evaluates the multivariate Kendall regression coefficient
#' for two regressors and split regressor region.
#' It describes how well the response variable can be fitted in each split region
#' by a function which increases or decreases as the regressors increase.
#'
#' @usage kendrs(x,y,splitp=NULL)
#' @param x data matrix of regressor variables with two columns,
#' @param y data vector of the response variable
#' @param splitp vector of length 2 of the splitting points,
#' If p1 is the first component of this vector, then the point splits the domain of the
#' first regressor into a left region of fraction p1 of data items and a right region
#' of the remaining data items. The same is done for the second regressor. As the
#' result we obtain 4 subregions of the regressor domain. default=c(0.5,0.5)
#' @return A list of Kendall regression coefficients for the 4 split regions
#' (11=left lower region, 12=left upper region...) with components
#'   \describe{
#'     \item{dcoeff++}{split coefficient ++}
#'     \item{dcoeff+-}{split coefficient +-}
#'     \item{totalcoeff}{total coefficient}
#'     \item{directions}{optimal directions}
#'   }
#' direction ++ means that y increases whenever both regressors increases
#' direction +- means that y increases whenever the first regressor increases and the
#' other regressor decreases..etc.
#' @references Eckhard Liebscher (2021). Kendall regression coefficient. Computational Statistics and Data Analysis 157 (2021). 107140
#' @examples
#' library(MASS)
#' data<- gilgais
#' kendrs(data[,1:2],data[,3],splitp=c(0.4,0.6))
#' @export

kendrs<-function(x,y,splitp=NULL){
  n<- length(y)
  cg<- 0.0  # total coefficient
  if (length(x[1,])>2) {x<- x[,1:2]} #the first two columns are chosen for d>2

  if (is.null(splitp)) {
    splitp<- c(0.5,0.5)}  # Split points in the middle
  ret<- list()
  dca<- dcb<- array(0.0,dim=c(2,2,2))
  dc<- vector(length=4)
  s1<- stats::quantile(x[,1],splitp[1])
  s2<- stats::quantile(x[,2],splitp[2])
  us1<- (x[,1]<=s1)  #logical Vector for split region 1
  us2<- (x[,2]<=s2)  #logical Vector for split region 2
  for (k in 1:2){
    k0<- 1
    for (l1 in 1:2){  #Split for u[,1]
      l01<- (l1==1)
      for (l2 in 1:2){  #Split for u[,2]
        l02<- (l2==1)
        x0<- as.matrix(subset(x,(us1==l01)&(us2==l02))) #choose the subsets
        y0<- subset(y,(us1==l01)&(us2==l02))
        n1<- length(x0[,1])  #number of items
        if (k==2) {   #k=1 --> "++" coefficient, k=2 --> "+-" coefficient
          x0[,2]<- -x0[,2]
        }
        v<- as.matrix(cbind(x0,y0))
        B<- (sum(copula::F.n(x0,x0,offset = 0, smoothing = "none"))-1.0)/(n1-1.0) #denominator
        A<- (sum(copula::F.n(v,v,offset = 0, smoothing = "none"))-1.0)/(n1-1.0) #numerator
        if (k==1) {
          dc[k0]<- ifelse((B!=0.0),2*A/B-1.0,0.0)
        } else{
          if (B!=0.0) {h<- 2*A/B-1.0}else{h<- 0.0}
          ret[[k0]]<- dc[k0]
          names(ret)[[k0]]<- paste0("dcoeff++",toString(l1),l2)
          ret[[k0+1]]<- h
          names(ret)[[k0+1]]<- paste0("dcoeff+-",toString(l1),l2)
        }
        dca[k,l1,l2]<- n1*A
        dcb[k,l1,l2]<- ifelse((B>0.0),n1*B,0.0)
        k0<- k0+2
      }}
  }
  mc<- 0.0
  for (t1 in 1:2){
    for (t2 in 1:2){
      for (t3 in 1:2){
        for (t4 in 1:2){
          h<- -1.0+2.0*(abs(dca[t1,1,1])+abs(dca[t2,1,2])+abs(dca[t3,2,1])+abs(dca[t4,2,2]))/
            (dcb[t1,1,1]+dcb[t2,1,2]+dcb[t3,2,1]+dcb[t4,2,2])
          if (h>mc) {
            mc<- h
            if (t1==1){ ic<- ifelse((dcb[t1,1,1]>0),"++","--")
            } else { ic<- ifelse((dcb[t1,1,1]>0),"+-","-+")}
            if (t2==1){ ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t2,1,2]>0),"+-","-+"))}
            if (t3==1){ ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t3,2,1]>0),"+-","-+"))}
            if (t4==1){ ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"++","--"))
            } else { ic<- paste(ic,ifelse((dcb[t4,2,2]>0),"+-","-+"))}
          }
        }}}}
  ret[[k0]]<- mc
  names(ret)[[k0]]<- "gencoeff"
  ret[[k0+1]]<- ic
  names(ret)[[k0+1]]<- "directions"
  return(ret)
}


#' Multivariate Kendall regression coefficient
#'
#' The function kendrm evaluates the multivariate Kendall regression coefficient.
#' It describes how well the response vector y can be fitted by a function of the regressor variables
#' which increases or decreases as the regressors increases.
#'
#' @usage kendrm(x,y,direction=NULL,out=0,outsd=TRUE,eps=0.95)
#' @param x data matrix of regressor variables
#' @param y data matrix of the response vector
#' @param out value 1: full output, value 0: reduced output of one coefficient
#' that is largest in absolute value, value 2: only the largest coefficient
#' @param direction vector of length d (d is number of regressors),
#' value 1 refers to regressors leading to increasing y whenever this regressor increases,
#' value -1 refers to regressors leading to decreasing y whenever this regressor increases.
#' If direction=NULL, then all coefficients are computed.
#' @param outsd logical. If TRUE, the estimated standard deviation and confidence intervals
#' are evaluated.  If FALSE, the evaluation is suppressed
#' @param eps real value, confidence level.
#' @return A list of the coefficients for several directions with components:
#'   \describe{
#'     \item{dcoeff}{Kendall regression coefficient}
#'     \item{dir}{direction vector}
#'     \item{Pxx}{fraction of \eqn{X \le \breve{X}} where \eqn{\breve{X}} has the same distribution as \eqn{X}}
#'     \item{sd}{standard deviation}
#'     \item{conf}{confidence interval}
#'   }
#' @references Eckhard Liebscher (2026). Kendall regression coefficient revisited.
#' @importFrom stats var qnorm
#' @examples
#' library(faraway)
#' y<- fat[,c("brozek","siri")]
#' x<- fat[,c("weight","adipos","abdom","hip")]
#' kendrm(x,y,direction=NULL,out=1,outsd=FALSE)
#' kendrm(x,y,direction=NULL,out=0,outsd=TRUE)
#' @export

kendrm<-function(x,y,direction=NULL,out=0,outsd=TRUE,eps=0.95){
  n<- nrow(x)
  d<- ncol(x)
  if (n!=nrow(y)){ stop("number of sample items in x and y are different")}
  if (n<4) {  stop("not enough sample items")}
  if ((eps<=0.5)|(eps>=1)){ stop("wrong epsilon")}
  if (ncol(y)<2) { stop("too few columns in matrix y")}

  x<- as.matrix(x)
  y<- as.matrix(y)
  mxc<- -10.0 #maximum coefficient
  if (is.null(direction)) {
    ka<-(2^d-1)
  } else {
    if (any(direction==0)|length(direction)!=d){stop("parameter direction not valid")}
    ka<- 0
    out<- 2}
  ret<- list()

  x<- cbind(x,-x)
  for (k in ka:0){
    if (ka>0){  #all directions
      direction<- 2*as.numeric(intToBits(k)[1:d])-1   #vector of signs, first element always plus 1
      if (k<ka)  { # direction does not consist solely of ones
        sel<- (direction<0)
        im<- (1:d) + d * sel   #select the columns of x
        u<- x[,im]
      }else{
        u<- x[,1:d]
      }
    } else {  #  one direction
      sel<- (direction<0)
      im<- (1:d) + d * sel
      u<- x[,im]
    }
    v<- cbind(u,y)
    Fuu<- F.n(u,u,offset = 0, smoothing = "none")
    Fvv<- F.n(v,v,offset = 0, smoothing = "none")
    Fyy<- F.n(y,y,offset = 0, smoothing = "none")

    PX<-(sum(Fuu)-1.0)/(n-1.0)
    PXY<-(sum(Fvv)-1.0)/(n-1.0)
    PY<-(sum(Fyy)-1.0)/(n-1.0)
    if ((PX<=0.0)|(PY>=1.0)) {
      if (out==1){
        ret[[k+1]]<- list(dcoeff="not defined",dir=direction)
      }else{  stop("Pxx=0, coefficient not defined")}
    }else{
      kendr<- (PXY-PX*PY)/(PX*(1.0-PY))    #Kendall coefficient
      ff<- TRUE
      if (out==1){  # output of all coefficients
        ret[[k+1]]<- list(dcoeff=kendr,dir=direction,Pxx=PX)
      }else{  ## reduced output of the largest absolute value
        if (out==2) {
          akendr<- kendr
        } else { akendr<-abs(kendr)}
        if (akendr>mxc) {
          mxc<- akendr
          ret<- list(dcoeff=kendr,dir=direction,Pxx=PX)
        } else {  ff<- FALSE}
      }

      if ((outsd)&ff) {
        vv<- Fvv+F.n(-v,-v,offset = 0, smoothing = "none")-
          (Fuu+F.n(-u,-u,offset = 0, smoothing = "none"))*
          (PY*(1.0-kendr)+kendr)-
          PX*((Fyy+F.n(-y,-y,offset = 0, smoothing = "none"))*
                (1.0-kendr)+2*kendr)
        sd<- sqrt(var(vv)/n)/(PX*(1.0-PY))
        z<- qnorm((1+eps)/2)
        if (out==1){
          ret[[k+1]][["sd"]]<-sd
          ret[[k+1]][["conf"]]<- c(kendr-sd*z,kendr+sd*z)
        }else{
          ret[["sd"]]<- sd
          ret[["conf"]]<- c(kendr-sd*z,kendr+sd*z)
        }
      } # end outsd
    }
  } # end k-loop
  return(ret)
}


#' Multivariate Kendall's tau
#'
#' The function kendtaum evaluates the multivariate Kendall's tau coefficient.
#' It describes the dependence of the variables in the data matrix
#'
#' @usage kendtaum(x,outsd=TRUE,eps=0.95)
#' @param x data matrix of regressor variables
#' @param outsd logical. If TRUE, the estimated standard deviation and confidence intervals
#' are evaluated.  If FALSE, the evaluation is suppressed
#' @param eps real value, confidence level.
#' @return A list of the coefficients for several directions with components:
#'   \describe{
#'     \item{dcoeff}{multivariate Kendall's tau}
#'     \item{sd}{standard error}
#'     \item{conf}{confidence interval}
#'   }
#' @references Schmid, F.; Schmidt, R.; Blumentritt, T.; Gaißer, S.; Ruppert, M. (2010).
#' Copula-based measures of multivariate association. in F. Durante, W. Härdle,
#' P. Jaworski, T. Rychlik (eds.) Copula Theory and Its Applications. Springer Berlin, 2010.
#' @importFrom stats var qnorm
#' @examples
#' library(MASS)
#' data <- gilgais
#' kendtaum(data[,1:4],outsd=TRUE)
#' @export

kendtaum<-function(x,outsd=TRUE,eps=0.95){
  n<- nrow(x)
  d<- ncol(x)
  x<- as.matrix(x)
  if (n<4) {  stop("not enough sample items")}
  if ((eps<=0.5)|(eps>=1)){ stop("wrong epsilon")}
  if (d==1){stop("x has only one column")}

  Fxx<- F.n(x,x,offset = 0, smoothing = "none")
  d2<- 2^(d-1)
  ret<- list(dcoeff=((2.0*d2*(sum(Fxx)-1.0)/(n-1.0))-1.0)/(d2-1.0))
  if (outsd){
    vv<- Fxx+F.n(-x,-x,offset = 0, smoothing = "none")
    z<- qnorm((1+eps)/2)
    ret[["sd"]]<- 2.0*d2*sqrt(var(vv)/n)/(d2-1.0)
    ret[["conf"]]<- c(ret[[1]]-ret[[2]]*z,ret[[1]]+ret[[2]]*z)
  }
  return(ret)
}
