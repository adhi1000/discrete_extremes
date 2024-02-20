

# ===============
# Libraries
library(gsl)



# =====================
# AUXILIARY

print_loading = function(it,it.max,step=2/100,newline=F){
  if(it%%ceiling(it.max*step)==0){
    if(newline){
      cat("\r",round((it/it.max)*100),'%',"\n")      
    } else {
      cat("\r",round((it/it.max)*100),'%')
    }
  }
}


# BIC 
bic = function(nllh,nbr.par,nbr.obs) nbr.par*log(nbr.obs)+2*nllh


# Compute confidence intervals from hessian
# Based on asymptotic normality of the mle 
# input: 
# - par: the mean of the mle  
# - hessian: the hessian at par
get.par.conf.int = function(par,hessian,alpha=0.95,pm=F){
  I = solve(hessian) # inverse fisher information matrix
  var = diag(I) 
  u = par+qnorm(1-(1-alpha)/2)*sqrt(var)
  l = par+qnorm((1-alpha)/2)*sqrt(var)
  if(pm){
    return(u-par) 
  } else{
    return(data.frame(par=par, l=l, u=u))
  }
}
# Alternative:
# input:
# x: a function of some parameters
# grad: the gradient of x w.r.t. to the parameters
# hessian: the hessian at the parameters
get.conf.int = function(x,grad,hessian,pm=F){
  I = solve(hessian)
  var = grad.stat %*% I %*% grad.stat
  u = par+qnorm(1-(1-alpha)/2)*sqrt(var)
  l = par+qnorm((1-alpha)/2)*sqrt(var)
  if(pm){
    return(u-par) 
  } else{
    return(data.frame(x=x, u=u, l=l))
  }
}

# ===========================
# Create a stepfun object from a cdf
# The step function f(x) is 0 when x<a, takes the values of the cdf for x\in[a,b), and is 1 for x>= b
from.pdist.to.stepfun = function(pdist,a,b){
  x = a:b
  y = c(0,sapply(a:(b-1),pdist),1)
  return(stepfun(x,y,f=0)) # f=0: right-continuous
}
stepFun.test = from.pdist.to.stepfun(pdist=function(x) (x+1)/6,a=0,b=2)
stopifnot(stepFun.test(-0.001)==0  & stepFun.test(0)>0 & stepFun.test(1)<1 & stepFun.test(2)==1)


# ===========================
# DISTRIBUTIONS
# ============================


# ============================

# Generalized Pareto distribution

# is x in the support 
in.gpd.sup = function(x,xi,si,mu){
  if(xi>=0){
    return(x>=mu) 
  } else
    return(x>=mu && x<= mu-si/xi)
}

# cdf 
p.gpd = function(x,xi,si,mu=0){
  if(x<mu){
    return(0)
  } else if(in.gpd.sup(x,xi=xi,si=si,mu=mu)){
    return(1-(1+xi*(x-mu)/si)^(-1/xi))
  } else {
    return(1)
  }
}

# density
d.gpd = function(x,xi,si,mu=0){
  if(in.gpd.sup(x,xi=xi,si=si,mu=mu)){
    if(xi==0){
      return(exp(-(x-mu)/si)/si)
    } else{
      return((1+xi*(x-mu)/si)^(-1/xi-1)/si)
    }
  } else 
    return(0)
}


# log density
log.d.gpd = function(x,xi,si,mu=0){
  if(in.gpd.sup(x,xi=xi,si=si,mu=mu)){
    if(xi==0){
      return(-(x-mu)/si-log(si) )
    } else{
      return((-1/xi-1)*log(1+xi*(x-mu)/si) -log(si))
    }
  } else 
    return(0)
}
stopifnot(round(log(d.gpd(30,4,5,6))-log.d.gpd(30,4,5,6),6)==0)
stopifnot(round(log(d.gpd(30,0,5,6))-log.d.gpd(30,0,5,6),6)==0)


# quantile
q.gpd = function(x,xi,si,mu=0){
  r = rep(NA,length(x))
  is.in.range = x<=1 & x>=0
  r[is.in.range] = (1/xi)*si*((1-x[is.in.range])^(-xi)-1)+mu
 return(r)
}
if(F){ # depreciated
q.gpd = function(x,xi,si,mu=0){
  if(x>1 || x<0){
    return(NA)
  } else{
    return((1/xi)*si*((1-x)^(-xi)-1)+mu)
  }
}
}
stopifnot(round(p.gpd(q.gpd(0.4,2.1,1.1,4),2.1,1.1,4),6)==0.4)

# sample
r.gpd = function(n,xi,si,mu=0) sapply(runif(n),q.gpd,xi=xi,si=si,mu=mu)

# Discrete generalized Pareto 
# A discrete GPD rv is defined as X = floor(Y) such that Y is GPD
# It has survival function Pr(X>=k)=Pr(Y>=k) and cdf Pr(X<=k)=Pr(Y<=k+1)
pm.dgpd = function(x,xi,si,mu=0) p.gpd(x+1,xi=xi,si=si,mu=mu)-p.gpd(x,xi=xi,si=si,mu=mu)
p.dgpd = function(x,xi,si,mu=0) p.gpd(x+1,xi=xi,si=si,mu=mu)
q.dgpd = function(x,xi,si,mu=0) q.gpd(x,xi=xi,si=si,mu=mu)-1
r.dgpd = function(n,xi,si,mu=0) floor(r.gpd(n,xi=xi,si=si,mu=mu))

if(F){ # wrong, use "pmf.rigth.trunc.dist" instead
# right-truncated GPD 
# truncation point: m
# discrete right-truncated GPD: X = floor(Y) s.t. Y is a truncated GPD. Thus its maximal possible value is m-1, not m
p.tgpd = function(x,xi,si,m,mu=0) ifelse(x>=m,1,p.gpd(x,xi=xi,si=si,mu=mu)/p.gpd(m,xi=xi,si=si,mu=mu))
pm.tgpd = function(x,xi,si,m,mu=0) ifelse(x>=m,0,p.tgpd(x+1,xi=xi,si=si,m=m,mu=mu)-p.tgpd(x,xi=xi,si=si,m=m,mu=mu))
q.tgpd = function(x,xi,si,m,mu=0)  q.gpd(p.gpd(m,xi=xi,si=si,mu=mu)*x,xi=xi,si=si,mu=mu)
stopifnot(round(p.tgpd(q.tgpd(0.3,2,3,4),2,3,4)-0.3,6)==0)
stopifnot(pm.tgpd(1,xi=1,si=1,m=4)==pm.dgpd(1,xi=1,si=1)) # wrong
}


# right truncated discrete distribution
# truncation point: m
pmf.rigth.trunc.dist = function(x,pmf,cdf,lower,upper){
  stopifnot(length(x)==1)
  stopifnot(is.finite(lower)&&is.finite(upper))  
  if(x<lower || x>upper) return(0)
  if(x>=lower && x<upper) return(pmf(x))
  if(x==upper) return(1-cdf(upper-1))
}

stopifnot(sum(sapply(c(0,1,2,3,4,5),function(x) pmf.rigth.trunc.dist(x,pmf=function(x) dpois(x,lambda=1),cdf=function(x) ppois(x,lambda=1),lower=0,upper=4)))==1)
stopifnot(sum(sapply(c(-1,0,1,2,3,4),function(x) pmf.rigth.trunc.dist(x,pmf=function(x) dpois(x,lambda=1),cdf=function(x) ppois(x,lambda=1),lower=0,upper=4)))==1)
stopifnot(all.equal(sum(sapply(c(0,1,2,3),function(x) pmf.rigth.trunc.dist(x,pmf=function(x) dpois(x,lambda=1),cdf=function(x) ppois(x,lambda=1),lower=0,upper=4))),ppois(3,lambda=1)))
stopifnot(sapply(c(0,1,2,3),function(x) pmf.rigth.trunc.dist(x,pmf=function(x) dpois(x,lambda=1),cdf=function(x) ppois(x,lambda=1),lower=0,upper=4))==sapply(c(0,1,2,3),dpois,lambda=1))
stopifnot(pmf.rigth.trunc.dist(1,pmf=function(x) pm.dgpd(x,xi=1,si=1),cdf=function(x) p.dgpd(x,xi=1,si=1),lower=0,upper=4)==pm.dgpd(1,xi=1,si=1))




# ============================

# Zipf-Mandelbrot distribution

# cdf
p.zm = function(x,s,q) ifelse(x>=0, 1-gsl::hzeta(s,x+q+1)/gsl::hzeta(s,q), 0)


# pmf
pm.zm = function(x,s,q){
  ifelse((x>=0) && ((1+1/s>0) || x<=-q), (x+q)^(-s)/gsl::hzeta(s,q), 0)
}
stopifnot(round(sum(sapply(0:1000, function(x) pm.zm(x,s=3.2,q=2.1))),5)==1)
stopifnot(round(sum(sapply(0:10, function(x) pm.zm(x,s=3.2,q=2.1)))- p.zm(10,s=3.2,q=2.1),5)==0)

# The normalizing constant is the Hurwitz-Zeta function:
# gsl::hzeta(s,q) = sum_{i=0}^\infty (i+q)^{-s}
pm.zm.norm = function(s,q) gsl::hzeta(s,q)
stopifnot(round(sum( 1/(rep(0:100000)+2)^3 )-gsl::hzeta(s=3,q=2))==0)

# unnormalized pmf
pm.zm.un = function(x,s,q) ifelse((x>=0)&& ((1+1/s>0) || x<=-q), (x+sign(1/(s-1))*q)^(-s), 0)

# continuous Zipf-Mandelbrot distribution
# extend zm to a continuous distribution czm with rv Y s.t. floor(Y) = X is zm
# remark: pr(Y<=k)=pr(X<=k-1)
p.czm = function(x,s,q) ifelse(x>=0, 1-gsl::hzeta(s,x+q)/gsl::hzeta(s,q), 0)
q.czm = function(x,s,q,upper=99999) uniroot( function(y) p.czm(y,s=s,q=q)-x, interval=c(0,upper))$root
stopifnot(round(p.czm(q.czm(0.4,2.1,1.1),2.1,1.1)-0.4,4)==0) 
stopifnot(round(p.czm(q.czm(0.9999,4.1,3.1),4.1,3.1)-0.9999,6)==0)

if(F){#  wrong, use "pmf.rigth.trunc.dist" instead
# right-truncated Zipf-Mandelbrot 
p.tczm = function(x,s=s,q=q,m) ifelse(x>=m,1,p.czm(x,s=s,q=q)/p.czm(m,s=s,q=q))
pm.tczm = function(x,s=s,q=q,m) ifelse(x>=m,0,p.tczm(x+1,s=s,q=q,m=m)-p.tczm(x,s=s,q=q,m=m))
q.tczm = function(x,s=s,q=q,m)  q.czm(p.czm(m,s=s,q=q)*x,s=s,q=q)
stopifnot(round(p.tczm(q.tczm(0.3,2,3,4),2,3,4)-0.3,6)==0)
pm.tzm = function(x,s,q,m) ifelse(x>=0 && x<=m-1 && ((1+1/s>0) || x<=-q), (x+sign(1/(s-1))*q)^(-s)/(sum((q+0:(m-1))^(-s))), 0)
pm.tzm.un = function(x,s,q,m) ifelse(x>=0 && x<=m-1 && ((1+1/s>0) || m-1<=-q), (x+sign(1/(s-1))*q)^(-s), 0)
pm.tzm.norm = function(s,q,m) ifelse( ((1+1/s>0) || m-1<=-q), sum((sign(1/(s-1))*q+0:(m-1))^(-s)),NA)
stopifnot(round(pm.tzm(0,s=3.4,q=3,m=3)+pm.tzm(1,s=3.4,q=3,m=3)+pm.tzm(2,s=3.4,q=3,m=3),6)==1)
}




# Generalized Zipf-Mandelbrot distribution

# cdf
p.gzm = function(x,xi,si){
    if(x<0) return(0)
    
    if(xi>0){
        return(p.zm(x,s=1+1/xi,q=si/xi))
    } else if(xi==0){
        return(exp(-x/si))
    } else if(xi<0){
        e = floor(si/abs(xi))
        if(x>=e) return(1)
        f = function(x) (1+xi*x/si)^(-(1+1/xi))
        S = sum(sapply(0:x,f))
        return(S/(S+sum(sapply((x+1):e,f))))
    }
}

# pmf
pm.gzm = function(x,xi,si){
    if(x<0) return(0)
    if(xi>0){
        return(pm.zm(x,s=1+1/xi,q=si/xi))
    } else if(xi==0){
        return(exp(-x/si)*(1-exp(-1/si)))
    } else if(xi<0){
        e = floor(si/abs(xi))
        if(x>e) return(0)
        f = function(x) (1+xi*x/si)^(-(1+1/xi))
        S = sum(sapply(0:e,f))
        return(f(x)/S)
    }
}
stopifnot(round(p.gzm(3,xi=-0.2,si=10)-sum(sapply(0:3,pm.gzm,xi=-0.2,si=10)),6)==0)
stopifnot(round(p.gzm(100,xi=-0.2,si=10)-sum(sapply(0:100,pm.gzm,xi=-0.2,si=10)),6)==0)


# pmf unnormalized
pm.gzm.un = function(x,xi,si){
    if(x<0) return(0)
    if(xi>0){
        return(pm.zm.un(x,s=1+1/xi,q=si/xi))
    } else if(xi==0){
        return(exp(-x/si)*(1-exp(-1/si)))
    } else if(xi<0){
        e = floor(si/abs(xi))
        if(x>e) return(0)
        f = function(x) (1+xi*x/si)^(-(1+1/xi))
        return(f(x))
    }
}

# normalization
pm.gzm.norm = function(xi,si){
    if(xi>0){
        return(pm.zm.norm(s=1+1/xi,q=si/xi))
    } else if(xi==0){
        return(1)
    } else if(xi<0){
        e = floor(si/abs(xi))
        f = function(x) (1+xi*x/si)^(-(1+1/xi))
        return(sum(sapply(0:e,f)))
    }
}
stopifnot(round(pm.gzm.un(0.2,3,2)/pm.gzm.norm(3,2)-pm.gzm(0.2,3,2),6)==0)
stopifnot(round(pm.gzm.un(0.2,-0.3,5)/pm.gzm.norm(-0.3,5)-pm.gzm(0.2,-0.3,5),6)==0)


# continuous Zipf-Mandelbrot distribution
# extend zm to a continuous distribution czm with rv Y s.t. floor(Y) = X is zm
# remark: pr(Y<=k)=pr(X<=k-1)
p.cgzm = function(x,xi,si) ifelse(x>=0,p.gzm(x-1,xi=xi,si=si),0)
q.cgzm = function(x,xi,si,upper=99999) uniroot( function(y) p.cgzm(y,xi=xi,si=si)-x, interval=c(0,upper))$root
q.cgzm = function(x,xi,si) q.cgzm(x,s=1+1/xi,si=si/xi)


q.cgzm = function(x,xi,si,upper=99999) uniroot( function(y) p.cgzm(y,xi=xi,si=si)-x, interval=c(0,upper))$root
stopifnot(round(q.cgzm(p.cgzm(33.2,4.1,3.1),4.1,3.1)-33.2,6)==0)
stopifnot(round(q.cgzm(p.cgzm(2.2,-0.5,3.1),-0.5,3.1)-2.2,6)==0)

if(F){
stopifnot(round(p.cgzm(q.cgzm(0.4,2.1,1.1),2.1,1.1)-0.4,4)==0) 
stopifnot(round(p.cgzm(q.cgzm(0.9999,4.1,3.1),4.1,3.1)-0.9999,6)==0)
}

# Gradient of probability of an Hurwitz Zeta to be >= x
# w.r.t. (s,q)
get.grad.sprob.zm = function(x,pu,s,q){
    f = function(par){
        pu*(1-p.zm(x-1,s=par[1],q=par[2]))
    }
    return(numDeriv::grad(f,c(s,q)))
}

if(F){ # comparision density DGPD and ZM
v = seq(0,10,length.out=100)
xi.sim=0.5;si.sim=0.3
plot(sapply(v,pm.gzm,si=si.sim,xi=xi.sim),sapply(v,pm.dgpd,si=si.sim,xi=xi.sim));abline(0,1)
}
     


# Create a quantile-quantile plot 
# input:
# - data: observations
# - q.fct: a quantile function
# - only quantiles corresponding to probabilites between pmin and pmax are plotted
# - only a percentage (prop) of the observations will be used between a and b to avoid producing heavy images

QQplot.old2 = function(data,q.fct,pmin=0,pmax=1,a=0,b=0,prop=1,ci=F,xlab="exp. quantile",ylab="th. quantile",
                       type=1,alpha=0.95,log=F,B){
  l = function(x){
    if(log){
      return(log(x))
    } else {
      return(x)
    }
  }
  
  n=length(data)
  t=seq(1,n)/(n+1) 
  if(prop<1){
    w = (t<a) | (t>b)
    s = sample.int(length(w[w==F]),round(prop*length(w[w==F])))
    t = c(t[w],t[w==F][s])
  }
  t=t[t>=pmin]
  t=t[t<=pmax]
  q.exp = quantile(data,t,type=type) 
  q.th = sapply(t,q.fct)
  plot(l(q.exp),l(q.th),xlab=xlab,ylab=ylab)
  abline(0,1)
  
  
  if(ci){ # plot confidence interval
    # simulate from the model B times
    # at each iteration, compute quantiles at points t and save the vector in a row of mat
    # compute confidence interval for each quantile, using the columns of mat
    mat = matrix(NA,B,length(t)) 
    for(i in 1:B){
      x.sim = q.fct(runif(n))
      q.sim = quantile(x.sim,t,type=type) 
      mat[i,] = q.sim
    }
    stopifnot(sum(is.na(mat))==0)
    q.u = apply(mat,2,function(x) quantile(x,1-(1-alpha)/2))
    q.l = apply(mat,2,function(x) quantile(x,(1-alpha)/2))
    
    ind.u = sort.int(q.u, index.return=TRUE)$ix 
    ind.l = sort.int(q.l, index.return=TRUE)$ix
    lines(l(q.u[ind.u]),l(q.th[ind.u]))
    lines(l(q.l[ind.l]),l(q.th[ind.l]))
  }
}