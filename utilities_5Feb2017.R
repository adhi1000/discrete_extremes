library(plyr)

# ==============================
# TRANSFORMATION IN PSEUDO ANGULAR COMPONENTS
pseudo.angular.coordinates = function(x){
  if(is.matrix(x)==F){
    coord = list(NA)
    coord$r = sum(x)
    coord$w = x/sum(x)
    return(coord)
  }
  else{
    n=dim(x)[1]
    d=dim(x)[2]
    r = apply(x,1,sum)
    w = matrix(NA,n,d)
    for(i in 1:n){w[i,] = x[i,]/r[i]}
    coord = list(NA)
    coord$r = r
    coord$w = w
    return(coord)
  }
}
x.test = cbind(c(1,2,3),c(4,5,6))
coord=pseudo.angular.coordinates(x.test)
stopifnot(coord$r==c(5,7,9))
stopifnot(coord$w==cbind(c(1,2,3)/coord$r,c(4,5,6)/coord$r))
x.test = c(1,3)
coord=pseudo.angular.coordinates(x.test)
stopifnot(coord$r==4)
stopifnot(coord$w==c(1/4,3/4))
# ==============================

# ===============================
# REJECTION SAMPLING
rejection.sampling = function(r.g,f.Mg){
  r = 0
  upper.lim = 10000000
  x = NA
  while(r<upper.lim){
    u = runif(1) 
    x = r.g(n=1)
    r = r+1
    if(u<f.Mg(x)){r=upper.lim+1}
  }
  if(r==upper.lim){warning("Attention: upper limit reached in rej sampling")}
  return(x)
}
# r.g: auxiliary function from which to sample
# f.Mg: ratio target density f over cte * auxiliary

########################################
# difference négative des logs (neg log returns)
ndifflog = function(v){
  return(-diff(log(v)))
}
# test
stopifnot(ndifflog(c(1:3)) == c(-(log(2)-log(1)), -(log(3)-log(2))))

########################################
# create a ts object from a list, indexed at 1 for the first observation, 
# with a given start and frequency
createTs = function(alist,start,frequency){
  ats = ts(alist/alist[1],start=start,frequency=frequency)
  return(ats)
}

########################################
# Plot n time series. 
# Inputs: a matrix of the ts in the columns, a title, a path, a list of colors
plotMyTs = function(tsMatrix,colors,xlab,ylab,path){ 
  stopifnot(length(colors)==dim(tsMatrix)[2])
  pdf(pointsize=13,width=11,height=6.5,file=path,family="Bookman")
  plot(tsMatrix[,1],col=col[1],ylim=range(na.omit(tsMatrix)),
       xlab=xlab,ylab=ylab,type="l")
  for(i in c(2:dim(tsMatrix)[2])){
    lines(tsMatrix[,i],col=col[i],lty=2)
  }
  dev.off()
}
########################################
# divise un nombre n = s*x+r*y tel que dist(x-y) est minimisé
# Remarque: utile pour mfrow
#getTwoDivisors = function(nbr){
#if( isPrime(nbr) ){ return(getTwoDivisors(nbr+1))}
#else {
# factors = findfactors(nbr)
# div = c( product(factors[1]) , product(factors[2:length(factors)]) )
#  for(i in c(2:(length(factor
#}
#}




########################################
# test if the number is prime
isPrime = function(number) {
  if(number==1){return(TRUE)}
  if(number==2){return(TRUE)}
  for (i in c(2:(number-1))){
    if (number %% i == 0){return(FALSE)}
  }
  return(TRUE)
}
stopifnot(isPrime(1))
stopifnot(isPrime(2))
stopifnot(isPrime(3))
stopifnot(isPrime(4)==FALSE)
stopifnot(isPrime(37))
stopifnot(isPrime(38)==FALSE)
########################################
# factor decomposition
findfactors <- function(n) {
  d <- c()
  div <- 2; nxt <- 3; rest <- n
  while( rest != 1 ) {
    while( rest%%div == 0 ) {
      d <- c(d, div)
      rest <- floor(rest / div)
    }
    div <- nxt
    nxt <- nxt + 2
  }
  d
}
# test
stopifnot(findfactors(3)==3)
stopifnot(findfactors(4)==c(2,2))
stopifnot(findfactors(18)==c(2,3,3))

##############################

# return all pairs of a vector in a matrix
library("sets")
getAllPairs = function(my.set){
  stopifnot(my.set==unique(my.set))
  my.set=sort(my.set)
  my.pairs = as.list( set_combn(my.set,2) )
  l = length(my.pairs)
  my.pairs.matrix = matrix(NA,l,2)
  for(i in c(1:l)){
    my.pairs.matrix[i,1] = as.list( as.list(my.pairs)[[i]] )[[1]]
    my.pairs.matrix[i,2] = as.list( as.list(my.pairs)[[i]] )[[2]]
  }
  return(my.pairs.matrix)
}

print(getAllPairs(c(2,1,5)))
# test
stopifnot(getAllPairs(c(2,1,5))==matrix(rbind(c(1,2),c(1,5),c(2,5)),3,2))
stopifnot(getAllPairs(c(1,2))==rbind(c(1,2)))  
#



# ========================================
# AIC 
aic = function(nllh,nbr.par){
  stopifnot(nbr.par==round(nbr.par))
  return(2*nbr.par+2*nllh)
}
# BIC 
bic = function(nllh,nbr.par,nbr.obs){
  stopifnot(nbr.par==round(nbr.par))
  return(nbr.par*log(nbr.obs)+2*nllh)
}


########################################
binarystuff = function(modelno, numberOfVars){
  mapply(function(i) modelno %/% 2^i - 2*modelno %/% 2^(i+1), (numberOfVars-1):0)
}

f.matrix=function(subids){
  fammembers=length(subids)
  n2=2^fammembers
  binmatrix1=matrix(rep(0,times=n2*fammembers),nrow=n2,ncol=fammembers)  ## binary matrix...all possible "simple" models (no interaction)
  for(i in 1:n2){
    binmatrix1[i,]=c(binarystuff(i-1,fammembers))  
  }
  size=c()
  for(j in 0:fammembers){size[j+1]=choose(fammembers,j)}
  colnames(binmatrix1)=subids
  binmatrix1
}

# ============================
# ATT: tout vérifier! semble que coincide pas avec grad(c(1,0))*Hessian^(-1)*grad(c(1,0))

# 95% confidence interval from hessian
# PRO: the covariance of the estimate is (asymptotically) 
# the inverse of the hessian.
from.varcov.to.ci = function(par, varcov, alpha=0.05,pm=pm){
  return(from.var.to.ci(par=par,var=diag(varcov),alpha=alpha,pm=pm))
}
from.var.to.ci = function(par, var, alpha=0.05,pm=F){
  stopifnot(length(par)==length(var))
  prop_sigma<-sqrt(var)
  lower<-par+qnorm(alpha/2)*prop_sigma
  upper<-par+qnorm(1-alpha/2)*prop_sigma
  if(pm==T){
    return((upper-lower)/2)
  } else{
    interval = NA
    interval<-data.frame(value=par, lower=lower,upper=upper)
    return(interval)
  }
}


from.hessian.to.ci = function(par, hessian, alpha=0.05,pm=F){
  inv_fisher_info<-solve(hessian)
  return(from.varcov.to.ci(par=par,varcov=inv_fisher_info,alpha=alpha,pm=pm))
}
stopifnot(round(from.hessian.to.ci(par=0.3,hessian=300,alpha=0.1,pm=T)-qnorm(0.95,sd=1/sqrt(300)),6)==0)

from.hessian.to.varcov = function(hessian){
  inv_fisher_info<-solve(hessian)
  return(inv_fisher_info)
}
#deviance" = (-2)*log(likelihood)

from.hessian.to.ci.2 = function(par,hessian){
  inv_fisher_info<-solve(hessian)
  prop_sigma<-sqrt(diag(inv_fisher_info))
  prop_sigma<-diag(prop_sigma)
  upper<- par+1.96*prop_sigma
  lower<- par-1.96*prop_sigma
  interval<-data.frame(value=par, upper=upper, lower=lower)
  return(interval)
}


# =====================================
# return the rows of a matrix which are also the rows of a second matrix
# ATT! not equals, but set equality...
row.of.mat1.is.a.row.of.mat2 = function(mat1,mat2){
  f = function(y){return(apply(mat1,1,setequal,y=y))}
  f0 = function(y){return(apply(mat2,1,setequal,y=y))}
  dim1 = dim(mat1)[1]
  if(is.null(dim1)){
    return(any(f0(mat1)))
  } else if(dim1==1){
    return((any(apply(apply(mat2,1,is.element,el=mat1),2,all))))
  }
  else{
    return(apply(apply(mat2,1,f),1,any))
  }
}
m1.test = rbind(c(1,2),c(4,5),c(1,5))
m2.test = rbind(c(4,5),c(2,5),c(1,2))
stopifnot(row.of.mat1.is.a.row.of.mat2(m1.test,m2.test)==c(T,T,F))
stopifnot(row.of.mat1.is.a.row.of.mat2(m2.test,m1.test)==c(T,F,T))
m1.test = rbind(c(1,2,3),c(4,5,6))
m2.test = rbind(c(4,5,6),c(2,2,4))
stopifnot(row.of.mat1.is.a.row.of.mat2(m1.test,m2.test)==c(F,T))
stopifnot(row.of.mat1.is.a.row.of.mat2(m2.test,m1.test)==c(T,F))
m1.test = c(1,2)
m3.test = c(1,5)
m2.test = rbind(c(4,5),c(2,5),c(1,2))
stopifnot(row.of.mat1.is.a.row.of.mat2(m1.test,m2.test)==T)
stopifnot(row.of.mat1.is.a.row.of.mat2(m3.test,m2.test)==F)

# ==========
# which rows of a matrix is the same set (or the same) than a given row

get.pos.row.in.mat = function(row,mat,set.equality=T){
  d=dim(mat)[1]
  if(set.equality){
  f = function(y) apply(mat,1,setequal,y=y)
  r = (1:d)[f(row)]
  ifelse(r==is.integer(0),NULL,r)
  return(r)
  } else{
    f = function(x) all(x==row)
    r = which(apply(mat,1,f))
   if(length(r)==0){ return(NULL)} else {return(r)}
}
}
stopifnot(get.pos.row.in.mat(c(1,2),rbind(c(1,2),c(4,5),c(1,5)))==1)
stopifnot(get.pos.row.in.mat(c(2,1),rbind(c(1,2),c(4,5),c(1,5)))==1)
stopifnot(get.pos.row.in.mat(c(1,2),rbind(c(4,5),c(1,5),c(1,2)))==3)
stopifnot(get.pos.row.in.mat(c(1,2),rbind(c(1,2),c(4,5),c(1,5),c(1,2)))==c(1,4))
stopifnot(get.pos.row.in.mat(c(1,2),rbind(c(4,5),c(1,5)))==integer(0))
stopifnot(get.pos.row.in.mat(c(1,2),rbind(c(1,2),c(4,5),c(1,5)),set.equality=F)==1)
stopifnot(is.null(get.pos.row.in.mat(c(2,1),rbind(c(1,2),c(4,5),c(1,5)),set.equality=F)))
stopifnot(is.null(get.pos.row.in.mat(c(2,1),rbind(c(4,5),c(1,5),c(1,2)),set.equality=F)))

# ===========
# Check if a package is loaded
is.loaded = function(mypkg){return(is.element(paste0('package:',mypkg), search()))}
# ==========
# do a qqplot
do.a.qqplot =  function(x,th.quantile.fn,log=F,line=F,...){
  k=(1:length(x))/(length(x)+1)
  if(log==F){
    plot(th.quantile.fn(k),quantile(x,k),... ,xlab='Theoretical Quantiles',ylab='Sample Quantiles')
    if(line==T) abline(0,1,col="blue")
  } else{
    plot(log(th.quantile.fn(k)),log(quantile(x,k)),log="xy",... ,xlab='Theoretical Quantiles',ylab='Sample Quantiles')
    if(line==T) abline(0,1,col="blue",untf=T)
  }
}
# not sure log qq plot is correct...
# 
# ===========
# do a tail qqplot
do.a.tail.qqplot = function(x,qfn,lower.q,line=F){
  l=length(x)
  plot(qfn((1:l)/(l+1))[quantile(x,(1:l)/(l+1))>lower.q],quantile(x,(1:l)/(l+1))[quantile(x,(1:l)/(l+1))>lower.q],
       xlab='Theoretical Quantiles',ylab='Sample Quantiles')
  if(line){abline(0,1,col='red',lwd=2)}
}

# =========
# Qqplot

# Create a quantile-quantile plot 
# input:
# - data: observations
# - q.fct: a quantile function
# - only quantiles corresponding to probabilites between pmin and pmax are plotted
# - only a percentage (prop) of the observations will be used between a and b to avoid producing heavy images

# Att: q.fct must be vectorized fct (otherwise set q.fct.is.vec=F)

QQplot = function(data,q.fct,pmin=0,pmax=1,a=0,b=0,prop=1,ci=F,xlab="exp. quantiles",ylab="th. quantiles",
                  type=1,alpha=0.95,log=F,B,cex.lab=1,lty=1,q.fct.is.vec=T,xaxt="s",yaxt="s",pcimin=0,pcimax=1){
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
  plot(l(q.exp),l(q.th),xlab=xlab,ylab=ylab,cex.lab=cex.lab,xaxt=xaxt,yaxt=yaxt)
  abline(0,1)
  
  
  if(ci){ # plot confidence interval
    # simulate from the model B times
    # at each iteration, compute quantiles at points t and save the vector in a row of mat
    # compute confidence interval for each quantile, using the columns of mat
    mat = matrix(NA,B,length(t)) 
    for(i in 1:B){
        if(q.fct.is.vec){
            x.sim = q.fct(runif(n))      
        } else { 
      x.sim = sapply(runif(n),q.fct)
      }
      q.sim = quantile(x.sim,t,type=type) 
      mat[i,] = q.sim
    }
    stopifnot(sum(is.na(mat))==0)
    q.u = apply(mat,2,function(x) quantile(x,1-(1-alpha)/2))
    q.l = apply(mat,2,function(x) quantile(x,(1-alpha)/2))
    
    ind.u = sort.int(q.u, index.return=TRUE)$ix 
    ind.l = sort.int(q.l, index.return=TRUE)$ix
    
    # cut the confidence interval 
    stopifnot(pcimax>pcimin)
    up.bd.ci = pcimax*max(q.u)
    low.bd.ci = pcimin*max(q.u)
    w = (q.u[ind.u]<=up.bd.ci)& (q.l[ind.l]<=up.bd.ci)& (q.u[ind.u]>=low.bd.ci)&(q.l[ind.l]>=low.bd.ci)
    
    # plot ci
    lines(l(q.u[ind.u][w]),l(q.th[ind.u][w]),lty=lty)
    lines(l(q.l[ind.l][w]),l(q.th[ind.l][w]),lty=lty)
  }
}

# depreciated
QQplot.old = function(data,q.fct,pmin=0,pmax=1,log=F,disc=F,type=6,light.a=0,light.b=0,light.prop=1){
  n=length(data)
  if(disc==F){
    t=seq(1,n)/(n+1)
  } else{
    table = table(data)
    cumsum = c(0,cumsum(table))
    nbr.val = length(table)
    t=rep(NA,nbr.val)
    for(i in 1:nbr.val){
      t[i] = (1+cumsum[i])/(n+1)
    } 
  }
  if(light.a!=0 | light.b!=0){
    
    w = (t<light.a) | (t>light.b)
    s = sample.int(length(w[w==F]),round(light.prop*length(w[w==F])))
    t = c(t[w],t[w==F][s])
  }
  
  t=t[t>pmin]
  t=t[t<pmax]
  q.exp = quantile(data,t,type=type) 
  q.th = sapply(t,q.fct)
  if(log==F){
    if((pmin==0)&&(pmax==1)){
      plot(q.exp,q.th)
    }else {
      warning("ATT")
      plot(q.exp,q.th,xlim=c(quantile(data,pmin),quantile(data,pmax)),ylim=c(quantile(data,pmin),quantile(data,pmax)))
    }
  } else{
    warning("ATT")
    plot(exp(q.exp),exp(q.th),xlim=c(quantile(data,pmin),quantile(data,pmax)),ylim=c(quantile(data,pmin),quantile(data,pmax))) #,log="xy" doesn't work...
  }
  abline(0,1)
}


# =========
# Box test plot
plot.box.test = function(x,lag.max,h.line=NULL,...){
  v = rep(NA,lag.max)
  for(i in 1:lag.max){
    v[i] = Box.test(x,lag=i)$p.value
  }
  plot(v,...)
  if(is.null(h.line)==F){abline(h=h.line,col='red')}
}
# =============================
# sample uniformly between the integers 1 and int.max
r.unif.discrete= function(n,int.max){
  return(floor(runif(n)*nbr.max+1))
}
# test   hist(r.unif.discrete(1000,6),breaks=100)
# =============================


# =============================
# Convert a list of list into a numerical matrix (deeply)
# (inner list must be of the same length)
list.of.list.as.matrix = function(listoflist,byrow){
  l=length(listoflist)
  len = rep(NA,l)
  for(i in 1:l){
    len[i]=length(listoflist[[i]])
  }
  stopifnot(length(unique(len))==1)
  n=l
  m=unique(len)
  if(byrow==F){
    n=unique(len)
    m=l
  }
  M=matrix(NA,n,m)
  for(i in 1:n){
    for(j in 1:m){
      if(byrow){
        M[i,j]=as.numeric(listoflist[[i]][[j]])
      }
      else{
        M[i,j]=as.numeric(listoflist[[j]][[i]])
      }
    }
  }
  stopifnot(sum(is.na(M))==0)
  return(M)
}
listoflist.test = list(
  list('dep'=0.4,'asy1'=0.7,'asy2'=0.3),
  list('dep'=0.7,'asy1'=0.6,'asy2'=0.4),
  list('dep'=0.3,'asy1'=0.7,'asy2'=0.6),
  list('dep'=0.2,'asy1'=0.3,'asy2'=0.2)
)
v.test = c(0.4,0.7,0.3,0.7,0.6,0.4,0.3,0.7,0.6,0.2,0.3,0.2)
m.test.truth.1 = matrix(v.test,byrow=T,4,3)
m.test.truth.2 = matrix(v.test,byrow=F,3,4)
m.test.1=list.of.list.as.matrix(listoflist.test,byrow=T)
m.test.2=list.of.list.as.matrix(listoflist.test,byrow=F)
stopifnot(m.test.truth.1==m.test.1)
stopifnot(m.test.truth.2==m.test.2)
# =============================
# print a loading indication in percent
print.loading = function(it,it.max,step=2/100){
  if(it%%ceiling(it.max*step)==0){
    cat(round((it/it.max)*100),'%','\n')
  }
}

# ===========================
# Estimation of the tail dependence coeff
# P(X>u|Y>u)
tail.dep = function(x,y,threshold){
  n=length(x)
  stopifnot(n==length(y))
  a1=sum((x>threshold)*(y>threshold))
  a2=sum(y>threshold)
  if(a2==0){
    result=0
  }
  else{
    result=a1/a2
  }
  stopifnot(is.numeric(result))
  return(result)
}
x.test=c(1,2,4,0,2)
y.test=c(0,3,1,5,4)
stopifnot(tail.dep(x.test,y.test,1.8)==2/3)
# ===========================
# Extremogram



extremogram = function(x,u,nbr.lag,ret=F,do.plot=T,permutationBounds=F,
                       boundPrecision=100,
                       xlim=c(2.1,nbr.lag),xlab='Lag',ylab='Extremogram'
                       ,...){
  n=length(x)
  stopifnot(length(nbr.lag)<=n)
  value.per.lag = rep(NA,nbr.lag)
  value.per.lag[1]=1
  for(h in 2:nbr.lag){
    x1=x[(h+1):n]
    x2=x[1:(n-h)]
    value.per.lag[h]=tail.dep(x1,x2,u)
  }
  if(do.plot){
    if(permutationBounds){
      K=boundPrecision
      x.perm=matrix(NA,K,nbr.lag)
      for(k in 1:K){
        x.perm[k,] = extremogram(sample(x),u,nbr.lag=nbr.lag,do.plot=F,ret=T)
        #print.loading(it=k,it.max=K,step=2/100)
      }
      f=function(x) quantile(x,0.95)
      ic.95 = apply(x.perm,2,f)
    }
    plot(value.per.lag,type="h",
         xlab=xlab,ylab=ylab,xlim=xlim,...)
    stopifnot(length(ic.95)==nbr.lag)
    bb=round(xlim)
    if(permutationBounds) lines(bb[1]:bb[2],ic.95[bb[1]:bb[2]],col='dodgerblue4')
  }
  if(ret) return(value.per.lag)
}
x=rnorm(500)
extremogram(x,u=quantile(x,0.95),nbr.lag=30,permutationBounds=T)
dev.off()



# ===========================
# Estimation of the tail dependence coeff
# P(X>u|Y>u)
asy.tail.dep = function(x,y,x.threshold,y.threshold){
  n=length(x)
  stopifnot(n==length(y))
  a1=sum((x>x.threshold)*(y>y.threshold))
  a2=sum(y>y.threshold)
  if(a2==0){
    result=0
  }
  else{
    result=a1/a2
  }
  stopifnot(is.numeric(result))
  return(result)
}
x.test=c(1,2,4,0,2)
y.test=c(0,3,1,5,4)
stopifnot(asy.tail.dep(x.test,y.test,0.5,2.1)==2/3)
stopifnot(asy.tail.dep(x.test,y.test,2.1,0.5)==1/4)
# ======================
which.min.matrix = function(mat){
  stopifnot(any(is.na(mat))==F)
  inds = which(mat == min(mat), arr.ind=TRUE)
  return(inds)
}
m.test = matrix(c(1,2,-3,3,5,6,2,-11,3),3,3)
stopifnot(which.min.matrix(m.test)==c(2,3))
# ==========================
library(MASS) 
library(RColorBrewer)
scatterplotProfile = function(X,...){
  stopifnot(dim(X)[2]==2)
  k <- 11
  my.cols <- rev(brewer.pal(k, "RdYlBu"))
  ## compute 2D kernel density, see MASS book, pp. 130-131
  z <- kde2d(X[,1], X[,2], n=50)  
  plot(X, xlab="X label", ylab="Y label", pch=19, cex=.4,...)
  contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
  abline(h=mean(X[,2]), v=mean(X[,1]), lwd=2)
  legend("topleft", paste("R=", round(cor(X)[1,2],2)), bty="n")
}

# ===================================
# Iteration in Binary
count.in.binary = function(n){
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  stopifnot(is.wholenumber(n)&&n>=2)
  c=rep(NA,2^n)
  for(i in 1:n){
    v1 = rep(0,2^(i-1))
    v2 = rep(1,2^(i-1))
    v=NA
    for(j in 1:(2^(n-1)/length(v1))){
      v =  c(v,v1,v2)
    }
    v=v[2:length(v)]
    
    c=cbind(c,v)
  }
  c= c[,2:dim(c)[2]]
  return(c)
}
stopifnot(as.matrix(count.in.binary(2))==matrix(c(0,0,1,0,0,1,1,1),4,2,byrow=T))
stopifnot(dim(count.in.binary(7))[1]==2^7)

# ==================================
inf.sum = function(f,from=1,maxit=1000,tol=0.0001){
  s=0
  for(i in from:maxit){
    t = f(i)
    s = s + t 
    if(abs(t/s)<tol){
      return(s)
    }
  }
  warning("In inf.sum: max iteration")
  return(s)
}
stopifnot(round(inf.sum(function(n) 2^(-n), from=0, maxit=1000,tol=0.0001),3)==2)

# ===============
# inverse a cdf numerically
# get numerical quantile
# (low,up): range on which the distribution lies
q.from.cdf = function(p,cdf,low=-99999,up=99999){
  f=function(q) cdf(q)-p
  return(uniroot(f,interval=c(low,up))$root)
}
# remark: up and low denotes the upper/lower endpoint of the distribution.
# test
stopifnot(round(pexp(q.from.cdf(0.3,pexp,low=0)),5)==0.3)
stopifnot(round(pexp(q.from.cdf(0,pexp,low=0)),5)==0)
stopifnot(round(pexp(q.from.cdf(0.999,pexp,low=0)),5)==0.999)

#======================
transform.sample.dist1.to.dist2 = function(x,pdist1,qdist2,fct.condition.on.x=NULL){
  if(is.null(fct.condition.on.x)){
    return(sapply(x,function(x) qdist2(pdist1(x)))) 
  } else{
    return(sapply(x,function(x) ifelse(fct.condition.on.x(x),qdist2(pdist1(x)),x))) 
  }
}
stopifnot(round(transform.sample.dist1.to.dist2(transform.sample.dist1.to.dist2(3.1,pdist1=pnorm,qdist2=qexp),pdist1=pexp,qdist2=qnorm)-3.1,6)==0)
stopifnot(round(transform.sample.dist1.to.dist2(2,pdist1=pnorm,qdist2=qexp,fct.condition.on.x= function(x) x!=2))==2)
qqnorm(transform.sample.dist1.to.dist2(runif(100),pdist1=punif, qdist2=qnorm));abline(0,1)

# get the cdf or the quantile fct of a distribution conditioned to be above a threshold u
pdist.cond.geq.u = function(pdist,u) function(x){ 1-(1-pdist(x))/(1-pdist(u))}
qdist.cond.geq.u = function(qdist,pdist,u) function(x){ qdist(1-(1-x)*(1-pdist(u)))}
f.test =pdist.cond.geq.u(pdist=pnorm,u=0.2)
f2.test = qdist.cond.geq.u(qdist=qnorm,pdist=pnorm,u=0.2)
stopifnot(round(f2.test(f.test(0.3))-0.3,6)==0)





# ====================
# Transform the columns of a matrix to uniform using their empirical cumulative function
# If remove1=T: remove entire rows with 1s
# If replace1=T: replace 1s with second largest value

transform.to.unif.using.ecdf = function(x,fct.condition.on.x=function(x) is.na(x)==F,remove1=F,replace1=F){
  
  d = dim(x)
  f = function(x) ifelse(fct.condition.on.x(x),ecdf(x[fct.condition.on.x(x)])(x),x)
  u = apply(x,2,f)
  count = 0
  
  if(remove1){ # remove rows containing 1s
    for(i in 1:d[2]){
      w = which(u[,i]==1)
      count = count + length(w)
      u = u[-w,]
      stopifnot(sum(u[,i]==1)==0) 
    }
    cat(count," rows containing 1s removed","\n")
  } else {
    if(replace1){
      # replace it by the largest value < 1 
      for(i in 1:d[2]){
        w = which(u[,i]==1)
        count = count + length(w)
        max.val.2 = max(u[,i][u[,i]<1])
        u[w,i] = max.val.2
        stopifnot(sum(u[,i]==1)==0) 
      }
      cat(count," 1s replaced by largest column value < 1","\n")
    }
  }
  return(u)
}

# Apply the ecdf function on both the positive and negative part of each
# column of a matrix to transform them to uniform
ecdf.on.pos.and.neg.marg = function(x,remove1=F,replace1=F){
  
  u = NULL
  for(i in 1:dim(x)[2]){
    v = x[,i] 
  w0 = v[v==0]
  w1 = w2 = NULL
  if(any(v>0))  w1 = ecdf(v[v>0])(v[v>0])
  if(any(v<0))  w2 = -ecdf(-v[v<0])(-v[v<0])
  w = c(w0,w1,w2)
  u = cbind(u,w)
  }
  
  count = 0
  if(remove1){ # remove rows containing +/- 1
    for(i in 1:dim(x)[2]){
      w = which(abs(u[,i])==1)
      count = count + length(w)
      u = u[-w,]
      stopifnot(sum(u[,i]==1)==0) 
    }
    cat(count," rows containing 1s removed","\n")
  } else {
    if(replace1){
      # replace +/- 1 by the largest value < 1 
      for(i in 1:dim(x)[2]){
        for(sign in c(-1,1)){
        w = which(u[,i]==sign)
        count = count + length(w)
        m.val.2 = sign*max(sign*u[,i][abs(u[,i])<1])
        u[w,i] = m.val.2
        }
        stopifnot(sum(u[,i]==1)==0) 
      }
      cat(count," 1s replaced by largest column value < 1","\n")
    }
  }
 return(u)
}

# ==========================

probMassPlot = function(table, probability.mass, xylim=NULL){
  stopifnot(is.table(table))
  s=sum(table)
  v=names(table)
  n=length(v)
  prob.th = rep(NA, n)
  prob.exp = rep(NA, n)
  for(i in 1: n){
    prob.th[i] = as.numeric(probability.mass(as.numeric(v[i])))
    prob.exp[i] = as.numeric(table[v[i]]/s)
  }
  plot(prob.th,prob.exp,xlab="prob th.", ylab="prob exp.")
  abline(0,1,col="blue")
}

PmPmPlot = function(table,pm,xylim=NULL){
  stopifnot(is.table(table))
  s=sum(table)
  v=names(table)
  n=length(v)
  y = rep(NA, n)
  for(i in 1:n){
    prob.th = as.numeric(pm(as.numeric(v[i])))
    prob.exp = as.numeric(table[v[i]]/s)
    y[i] = (prob.th-prob.exp) #/prob.th
  }
  plot(1:n,y,xlab="", ylab="(prob th. - prob exp.)/prob th.")
  abline(0,0,col="blue")
  cat("sum abs rel errors: ",sum(abs(y)),"\n")
}
PmPmPlot(table=table(rpois(1000,2)),pm=function(x) dpois(x,2))


discCdfPlot = function(table, cdf, xylim=NULL,scale=NULL){
  stopifnot(is.table(table))
  n0 = as.numeric(names(table[1]))
  n1 = as.numeric(names(table[length(table)]))
  prob.th = rep(NA, n1-n0+1)
  prob.exp = rep(NA, n1-n0+1)
  N = sum(table)
  
  cumsum = 0
  for(i in n0:n1){
    cumsum = cumsum + table[as.character(i)]
    prob.th[n0+i-1] = cdf(i)
    prob.exp[n0+i-1] = cumsum/N
  }
  if(is.null(scale)){
    plot(prob.th,prob.exp,xlab="prob th.", ylab="prob exp.")
    abline(0,1,col="blue")
  } else if(scale=="exp"){
    plot(-prob.th+1,-prob.exp+1,xlab="prob th.",ylab="prob exp.",log="xy")
    abline(0,1,col="blue")
  } else if(scale=="log"){
    plot(prob.th,prob.exp,xlab="prob th.",ylab="prob exp.",log="xy")
    abline(0,1,col="blue")
  }
}

# ===========================
# Create a stepfun object from a cdf
# The step function f(x) is 0 when x<a, takes the values of the cdf for x\in[a,b), and is 1 for x>= b
from.pdist.to.stepfun = function(pdist,a,b){
  stopifnot(is.finite(pdist(b)) & is.finite(pdist(a)))
  x = a:b
  y = c(0,sapply(a:(b-1),pdist),1)
  return(stepfun(x,y,f=0)) # f=0: right-continuous
}
stepFun.test = from.pdist.to.stepfun(pdist=function(x) (x+1)/6,a=0,b=2)
stopifnot(stepFun.test(-0.001)==0  & stepFun.test(0)>0 & stepFun.test(1)<1 & stepFun.test(2)==1)

# ===========================


# depreciated:
# Create a stepfun object from a cdf. 
# The step fct is 0 under a and 1 between b-1 and b
from.pdist.to.stepfun.old = function(pdist,a,b,silent=F){
  # stopifnot(pdist(a-1)==0)
  x = a:b
  y = c(0,sapply(a:(b-1),pdist),1)
  if(silent==F) cat("In from.pdist.to.stepfun: approx: ", 1-pdist(b)," between ",b-1," and ",b,"\n") 
  return(stepfun(x,y,f=0))
}
stopifnot(identical(from.pdist.to.stepfun(pdist=function(x) x/6,a=1,b=3)(-1:6),stepfun(1:3,c(0,1/6,2/6,1))(-1:6)))
fn.test = from.pdist.to.stepfun(pdist=function(x) ifelse(x<=6,x/6,1),a=0,b=6)
plot(fn.test)
stopifnot(fn.test(-50)==0 & fn.test(6)==1 & fn.test(50)==1 & round(fn.test(5)-0.8333333,6)==0)


# ===========================
# Count the number of characters (without counting spaces)
library(stringi)
nbr.char.without.space = function(x) stringi::stri_length(x) - stringi::stri_count_fixed(x, " ")
stopifnot(nbr.char.without.space("hello to you é")==11)


# ====================
# Truncated Cdf and Quantile
# If X has cdf, then its truncated cdf on [l,u] is: {F(x)-F(l)}/{F(u)-F(l)}
# and its quantile functions is: F^{-1}[{F(u)-F(l)}x+F(l)]
p.trunc.dist = function(p.dist,low=-Inf,up=Inf){
  function(x) (p.dist(x)-ifelse(is.infinite(low),0,p.dist(low)))/(ifelse(is.infinite(up),1,p.dist(up))-ifelse(is.infinite(low),0,p.dist(low))) 
}
q.trunc.dist = function(p.dist,q.dist,low=-Inf,up=Inf){
  function(x) q.dist( ( ifelse(is.infinite(up),1,p.dist(up))-ifelse(is.infinite(low),0,p.dist(low)))*x+ifelse(is.infinite(low),0,p.dist(low))) 
}
stopifnot(p.trunc.dist(p.dist=pnorm,low=0)(2.3)==(pnorm(2.3)-pnorm(0))/(1-pnorm(0)))
stopifnot(q.trunc.dist(p.dist=pnorm,q.dist=qnorm,low=0)(0.2)==qnorm(0.2*(1-pnorm(0))+pnorm(0)))


# ========================
# Transform a rv uniformly distributed to a given distribution truncated from below
from.unif.to.trunc.dist = function(x,pdist,qdist,low){
  
  qdist2 = q.trunc.dist(p.dist=pdist,q.dist=qdist, low=low) 
  
  y = transform.sample.dist1.to.dist2(x,pdist1=punif,qdist2=qdist2)
  return(y)
}


# ===========================
# Return which rows of M equal v
which.rows.equals = function(M,v) which(apply(M,1,function(x) all(x==v)))
stopifnot(which.rows.equals(M=matrix(c(3,2,3,4,5,6,7,7),ncol=2,nrow=4),v=c(3,7))==3)
stopifnot(which.rows.equals(M=matrix(c(3,2,3,7,7,6,7,7),ncol=2,nrow=4),v=c(3,7))==c(1,3))

# =======================
# Get probability that the rows of a matrix exceed a vector of thresholds
get.prob.exced = function(x,b,thresh=rep(0,dim(x)[2]),abs=T,strict=T){
  stopifnot(length(b)==length(thresh))
  stopifnot(length(b)==dim(x)[2])
  a = function(x) if(abs) return(abs(x)) else return(x)
  f = function(x) if(strict) all((a(x)>thresh)==b) else all((a(x)>=thresh)==b)
  v = apply(x,1,f)
  return(sum(v)/length(v))
}  
stopifnot(get.prob.exced(x=matrix(c(1,2,3,-3,2,3),byrow=T,nrow=2),b=c(0,0,1),thresh=c(5,5,5)/2,abs=T)==1/2)
stopifnot(get.prob.exced(x=matrix(c(1,2,3,-3,2,3),byrow=T,nrow=2),b=c(0,0,1),thresh=c(5,5,5)/2,abs=F)==1)
stopifnot(get.prob.exced(x=matrix(c(1,2),byrow=T,nrow=1),b=c(1,1),thresh=c(1,1),abs=F,strict=T)==0)  
stopifnot(get.prob.exced(x=matrix(c(1,2),byrow=T,nrow=1),b=c(1,1),thresh=c(1,1),abs=F,strict=F)==1)  

get.all.prob.exced = function(x,thresh=rep(0,dim(x)[2]),abs=T,strict=T,probs=F){
  stopifnot(length(thresh)==dim(x)[2])
  a = function(x) if(abs){ return(abs(x))} else{ return(x) }
  f = function(x) if(strict) as.numeric(a(x)>thresh) else as.numeric(a(x)>=thresh) 
  x.binary = t(apply(x,1,f))  
  t.binary = plyr::count(df=as.data.frame(x.binary))
  if(probs==F){
    return(t.binary)
  } else {
    t.binary[,"freq"]=t.binary[,"freq"]/dim(x)[1]
    colnames(t.binary)[dim(x)[2]+1]="probs"
    return(t.binary)
  }
} 
stopifnot(get.all.prob.exced(x=matrix(c(1,2,3,-3,2,3),byrow=T,nrow=3),thresh=c(3,3),abs=T,strict=F)[,"freq"]==c(1,1,1))
stopifnot(get.all.prob.exced(x=matrix(c(1,2,3,-3,2,3),byrow=T,nrow=3),thresh=c(3,3),abs=T,strict=T)[,"freq"]==3)

get.prob.exced.biv = function(x,y){
  # get indices where obs fall in 3 different quadrants: uu, ud, du
  uu = (x>0)*(y>0)==1
  ud = (x>0)*(y==0)==1
  du = (x==0)*(y>0)==1
  dd = (x==0)*(y==0)==1
  
  # get frequency
  nuu = sum(uu)
  nud = sum(ud)
  ndu = sum(du) 
  ndd = sum(du) 
  
  s= nuu+nud+ndu+ndd
  return(list(puu=nuu/s,pud=nud/s,pdu=ndu/s))
}
stopifnot(all(unlist(get.prob.exced.biv(c(1,2,0,0),c(0,1,0,1)))==unlist(list(puu=0.25,pud=0.25,pdu=0.25))))

# Get the rows of a matrix that exceed a vector of thresholds
get.exced.values = function(x,b,thresh,abs=T,strict=T){
  stopifnot(length(b)==length(thresh))
  stopifnot(length(b)==dim(x)[2])
  a = function(x) if(abs) return(abs(x)) else return(x)
  f = function(x) if(strict) all((a(x)>thresh)==b) else all((a(x)>=thresh)==b)
  v = apply(x,1,f)
  return(x[v,])
}
stopifnot(get.exced.values(x=matrix(c(1,2,3,-3,2,3),byrow=T,nrow=2),b=c(0,0,1),thresh=c(5,5,5)/2,abs=T)==c(1,2,3))

# ==================
# Obtain a correlation matrix from a vector
# Input: a vector corresponding to the upper triangular matrix read by row,
# i.e., m12, m13, m14,..., m23, m24,...

from.vector.to.corr.mat = function(x){
  m = diag(1)
  k = 1
  l = length(x)
  
  # find size of the matrix:
  # solve d(d-1)/2 = l 
  d = 1/2 * (1 + sqrt(1 + 8*l))
  
  stopifnot(d==round(d))
  m = diag(d)
  
  for(i in 1:d){ 
    for(j in 1:d){
      if(i<j){
        m[i,j] = m[j,i] = x[k]
        k = k+1
      }
    }
  }
  stopifnot(isSymmetric(m))
  return(m)
}

# inverse

from.corr.mat.to.vector = function(m) m[upper.tri(m)] 

m.test= matrix(c(1,0.4,0.7,0.4,1,-0.8,0.7,-0.8,1),nrow=3)
stopifnot(from.vector.to.corr.mat(m.test[upper.tri(m.test)])==m.test)
stopifnot(from.corr.mat.to.vector(c(2,3,4))==c(2,3,4))

# Compute return levels
# alpha: probability of being below the return level, e.g. 0.999
# pu: probability to exceed the threshold u
# q.exceed.dist: distribution of X-u|X>u
get.return.level = function(alpha,u,pu,q.exceed.dist){
  return(u+q.exceed.dist(1-(1-alpha)/pu))
}

# Compute return level of a GPD and its gradient w.r.t. (si,xi)
# alpha: return level = quantile alpha
# u: threshold; pu: prob exceed this threshold
# si, xi: scale and shape parameters of the GPD
get.return.level.gpd = function(alpha,u,pu,si,xi){
  z = (1-alpha)/pu
  return( u+(si/xi)*(z^(-xi)-1))
}
get.grad.return.level.gpd = function(alpha,u,pu,si,xi){
  z = (1-alpha)/pu
  return( c( (z^(-xi)-1)/xi, (si/xi^2)*(1-z^(-xi)*(1+xi*log(z))) ) ) # (df/dsi, df/dxi)
}

# Gradient of probability of a GPD to be >= x
# w.r.t. (xi,si)
get.grad.sprob.gpd = function(x,si,xi,pu){
  v1 = pu*(1+xi*x/si)^(-1/xi)*(-x*xi+(si+x*xi)*log(1+xi*x/si)) / (xi^2*(si+xi*x)) # d/dxi
  v2 = pu*x*(1+xi*x/si)^(-1/xi-1)/si^2 # d/dsi
  return(c(v1,v2)) 
}
# Gradient of probability mass of a discrete GPD in k
# w.r.t. (xi,si)
get.grad.pm.gpd = function(x,si,xi,pu){
  return(get.grad.sprob.gpd(x,si,xi,pu)-get.grad.sprob.gpd(x+1,si,xi,pu))
}


# =============
# Remove all rows from a matrix that are not finite (eg NA, NaN, +/- Inf)
remove.non.finite.rows = function(m) m[!rowSums(!is.finite(m)),]
stopifnot(all(remove.non.finite.rows(matrix(c(1,2,NA,NaN,1,Inf,-1,1,9,3),5))==matrix(c(2,1,-1,3),2)))

# ============
# EXAMPLE: EMPIRICAL QUANTILE INVERSE ECDF

x.test = rnorm(1000)
p.emp.pos = function(e,x) ecdf(x[x>0])(e)
q.emp.pos = function(e,x) as.vector(quantile(x[x>0],probs=e,type=6))
p.emp.pos(e=q.emp.pos(e=0.3,x=x.test),x=x.test)
a=p.emp.pos(e=x.test,x=x.test)
plot(q.emp.pos(e=a[a>0],x=x.test),x.test[x.test>0])

# =======================

two.dim.integral.loop = function(f,low1,up1,low2,up2){
  integrate(function(y) { 
    sapply(y, function(y) {
      integrate(function(x) f(x,y), low1,up1)$value
    })
  }, low2, up2)
}
p.test =0.3
test.f = function(x,y) exp(-(x^2-2*p.test*x*y+y^2)/(2*(1-p.test^2)))/(2*sqrt(1-p.test^2)*pi)
stopifnot(round(two.dim.integral.loop(test.f,low1=-Inf,up1=Inf,low2=-Inf,up2=Inf)$value,4)==1)

# Packages for integration 
# pracma::quad2d(test.f,xa=-20,xb=20,ya=-20,yb=20) # unreliable

# =============
# CHI-SQUARE TEST
lambda.test=100
z.test= rpois(100,lambda=lambda.test)
tz.test = setNames(tabulate(z.test), 1:max(z.test))
pp.test=sapply(as.numeric(names(tz.test)),dpois,lambda=lambda.test)
chisq.test(x=tz.test,p=pp.test/sum(pp.test))
# remark: does not work when state space too large

# ====================
# Normalize a covariance matrix
from.cov.to.corr = function(cov){
  D = diag(sqrt(diag(cov)))
  return(solve(D) %*% cov %*% solve(D))
}  
S.test = rbind(c(1.0,  1.0,  8.1), c(1.0, 16.0, 18.0), c(8.1, 18.0, 81.0))
U.test = rbind(c(1.0,  0.25,  0.9), c(0.25, 1, 0.5), c(0.9, 0.5, 1))
stopifnot(all(round(from.cov.to.corr(S.test)-U.test,6)==0))

# ======================
# Find the position / indices of the largest / maximal values in a matrix
get.indices.of.largest.val.in.mat = function(mat,decreasing=T){
    mat.od <- order(mat, decreasing = decreasing)
    mat.od.arr <- cbind(mat.od%%nrow(mat), mat.od%/%nrow(mat)+1)
    mat.od.arr[,2][mat.od.arr[,1]==0] <- mat.od.arr[,2][mat.od.arr[,1]==0] - 1
    mat.od.arr[,1][mat.od.arr[,1]==0] <- nrow(mat)
    return(mat.od.arr)
}
# ======================
# Log likelihood of a multivariate Normal distribution
# n is the sample size
# K is the precision matrix
# SIgma.emp is the empirical covariance matrix x^T x/ n
loglik.mvnorm = function(K, Sigma.emp, n){
    d = dim(K)[1]
    stopifnot(dim(Sigma.emp)[1]==d)
    return( -n/2 * ( d*log(2*pi) - log(det(K)) + sum(diag(K %*% Sigma.emp))))
}


# Print a table in latex
print.table.latex = function(table){
    text = ""
    for(i in 1:dim(table)[1]){
        for(j in 1:dim(table)[2]){
            text = paste(text, table[i,j])
            if(j<dim(table)[2]){
                text = paste(text, "&")
            }
        }
        text = paste(text, "\\\\","\n")
        
    }
    cat(text)
}

# print fixed number of decimals
prt.rd = function(x,rd) sprintf(paste0("%.",as.character(rd),"f"),x)
# build confidence intervals in latex
ci.in.latex= function(par,par.d,par.u,rd) as.character(paste0("$", prt.rd(par,rd), "_{[",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"]}","$"))
ci.in.latex.pm = function(par,pm,rd) as.character(paste0("$", prt.rd(par,rd), "_{ \\pm ",prt.rd(pm,rd),"}","$"))


# ================
# Conditional distribution of a multivariate Normal / Student distribution
# f(x)=(1+ (x-mu) Sigma (x-mu) / df)^(-(df+d)/2)
# mu: not the mean but a shift parameter mu such that X -> X+mu
# sigma: matrix Sigma (covariance matrix is proportional to Sigma) 
# input: for the moment, sigma has to be a correlation matrix
condMvtnorm = function(sigma,mu=rep(0,dim(sigma)[1]),df, dependent.ind, given.ind, X.given){
    
    stopifnot(diag(sigma)==1)
    stopifnot(length(mu)==dim(sigma)[2])
    stopifnot(length(given.ind)==length(X.given))
    stopifnot(length(given.ind)+length(dependent.ind)==length(mu))
    
    p1 = length(X.given)
    
    sigma.21 = sigma[dependent.ind,given.ind]
    sigma.12 = sigma[given.ind,dependent.ind]
    sigma.11 = sigma[given.ind,given.ind]
    sigma.22 = sigma[dependent.ind,dependent.ind]
    mu.1 = mu[given.ind]
    mu.2 = mu[dependent.ind]
    
    sigma.11.inv = solve(sigma.11)
    
    A = sigma.21 %*% sigma.11.inv
    mu.2g1 = mu.2 + A %*% (X.given-mu.1)
    sigma.22g1 = sigma.22 - A %*% sigma.12
    d1 = as.vector( t(X.given - mu.1) %*% sigma.11.inv %*% (X.given - mu.1) )
    lambda = (df+d1)/(df+p1)
    if(df==Inf) lambda=1
    return(list(condMu=mu.2g1,condVar=lambda * sigma.22g1,condDf=df+p1))
    
}
# test
sigma.test = matrix(c(1,0.4,0.4,1),byrow=T,nrow=2)
mu.test = c(0.2,-0.7)
df.test = 4

r1 = mvtnorm::dmvt(c(0.2,-0.1),df=df.test,sigma=sigma.test,delta=mu.test,log=F)
L= condMvtnorm(sigma=sigma.test,mu=mu.test,df=df.test, dependent.ind=2, given.ind=1, X.given=0.2)
r2 = dt(0.2-mu.test[1],df=df.test)*dt((-0.1-L$condMean)/sqrt(L$condVar),df=L$condDf)/sqrt(L$condVar) 
stopifnot(round(r1-r2,6)==0)
# check same as conMVNorm package when df = Inf
A = condMvtnorm(sigma=sigma.test,mu=mu.test,df=Inf, dependent.ind=2, given.ind=1, X.given=0.2)
A.condMVNorm = condMVNorm::condMVN(mean=mu.test,sigma=sigma.test,dependent.ind=2,given.ind=1,X.given=0.2)
stopifnot( round(sum(abs( c(A$condMean-A.condMVNorm$condMean,A$condVar-A.condMVNorm$condVar))),6)==0 && A$condDf==Inf)
# Remark
a1 = dt(0.2,ncp=0.3,df=3) # non-central Student (mu in condMvtnorm is not ncp)
a2 = dt(0.2-0.3,df=3) # non-standardize Student (mu in condMvtnorm is the shift mu in (x-mu)/si )
stopifnot(a1!=a2)
# check if same as conditional distribution of bivariate copula in package copula
# (involved test)
p.test=0.6; df.test=11.2
v.test = c(0.7,-0.3)
x.test = qt(v.test[1],df=df.test)
y.test = qt(v.test[2],df=df.test)
Dist = condMvtnorm(sigma=matrix(c(1,p.test,p.test,1),nrow=2,byrow=T),df=df.test, dependent.ind=2, given.ind=1, X.given=x.test)
r = pt((y.test-Dist$condMu)/sqrt(Dist$condVar),df=Dist$condDf)
# todo: 
if(F) stopifnot(round(copula::cCopula(v.test,copula::tCopula(p.test,df=df.test))[2]-r,6)==0)


# Bivariate cdf of the multivariate student
if(F){ #todo
p.biv.mvt = function(sigma,mu=rep(0,dim(sigma)[1]),df)
    stopifnot(diag(sigma)==1) # sigma must be a correlation matrix 
 condMvtnorm(sigma=sigma,mu=mu,df=df, dependent.ind=2, given.ind=1, X.given=x)

# only available for df integer valued
mvtnorm::pmvt(c(1,2),sigma=matrix(c(1,0.3,0.3,1),nrow=2,byrow=T),df=2.2)[1]

}


# ==============
    
                             