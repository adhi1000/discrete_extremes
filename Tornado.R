# ================
# FIT GPD, D-GPD, GZD to earthquake outcomes.


# FIT DIFFERENT MODELS TO EXCEEDANCES
# COMPARE PARAMETER AND QUANTiLE ESTIMATION


myRFolder = "~/Downloads/discrete_extremes"
setwd(myRFolder)
source("Functions.R")
source("utilities_5Feb2017.R")

# Graphical settings
pointsize = 18
pointsize.large = 12


# ===========

# IMPORT DATAScreen Shot 2017-06-24 at 12.37.09
# Article: Tipett, Lepore and Cohen, Extreme Tornado Outbreaks, Science, 216
# Unit: Tornado outbreaks are sequences of six or more tornadoes rated F1 and greater on the Fujita scale, 
# or rated EF1 and greater on the Enhanced Fujita scale, that occur in close succession

# Outbreak definition:
#1. Only tornadoes ratedE/F1 and greater
#2. Arrange tornadoes by start time
#3. Mark time gaps of more than 6h
#4. Keep clusters with more than 6 tornadoes  (defined as 12 in the paper!)
#5. Count out breaks/yr and tornadoes/outbreak


# Findings
# "increases in the annual mean and variance of the number of tornadoes per outbreak (6)"
# "Linear trends in the percentiles of the number of torna- does per outbreak (Fig. 1a) are positive, statistically signifi-
# cant and increase exponentially faster with percentile probability (Fig. 1b)"
# Figure 1a:  Annual 20th, 40th, 60th and 80th percentiles of the number of E/F1+ tornadoes per outbreak (6 or more E/F1+ tornadoes),
# 1954-2015 (solid lines), and quantile regression fits to 1965- 2015 assuming linear growth in time (dashed lines).

# We refer to outbreaks with 12 or more E/F1+ tor- nadoes as “extreme outbreaks” (2). 
# There were 435 extreme outbreaks 1965-2015, 
#no statistically significant trends in the annual number of extreme outbreaks (P = 0.66) (Fig. 2a), 
# and no statistically significant autocorrelation in the num- bers of tornadoes per extreme outbreak (fig. S2c).

# The GP distributions found here have shape parameter around 0.3
#The percentiles of the number of tornadoes per extreme outbreak (Fig. 2b) also have upward trends that are statisti- cally significant (above the 30th percentile) 

# EXTREME EASTERN U.S. WINTER OF 2015 NOT SYMPTOMATIC OF CLIMATE CHANGE, Laurie Trenary,


tornado = read.table(file="outbreakdata1950_2015.txt",sep=",")
dim(tornado)
tornado



# v: tornadoes per outbreak
v = as.vector(unlist(t(tornado[,2:dim(tornado)[2]])))
v = v[is.na(v)==F]
stopifnot(length(tornado[tornado>=12 & is.na(tornado)==F])-dim(tornado)[1]==length(v[v>=12])) 
plot(table(v))

v.t = NULL
for(i in 1:dim(tornado)[1]){
  v.t = c(v.t, rep(tornado[i,1], sum(is.na(tornado[i,2:dim(tornado)[2]])==F)))
}
stopifnot(length(v)==length(v.t))

# q: outbreaks/yr
q = apply(tornado[,2:dim(tornado)[2]],1,function(x) sum(x[is.na(x)==F]))
stopifnot(sum(tornado[1,][is.na(tornado[1,])==F])-1950==q[1])
plot(q)
plot(table(q))

initial.date = 1965
stopifnot(length(v[v>=12 & v.t>=initial.date])==435) # check corresponds to article




models = c("dgpd","zm","gpd2","tdgpd","tgpd","gpd","tgpd2")
models = c("dgpd","zm","gpd2","tdgpd","tgpd2")
models = c("dgpd","tgpd2","tdgpd")


models.qqplot = c("tdgpd")

# For each model, compute the following quantities:
# - par: maximum likelihood estimator 
# - val: log-likelihood at the mle 
# - par.pm: value pm such that [par-pm, par+pm] is a confidence interval for the mle (par.pm)
# - ks: Kolmogorov-Smirnov goodness-of-fit test for discrete data
# - latexName: name of the model in a latex tabular


lang="torn"
dataset = ""
start.year = 1965
data = v[v.t>=start.year]

qt = NA
u = 6 # in Tippett article: 6 in fig 1, 12 in fig 2 and table 1
u = 12 # in Tippett: u=12
s = data[data>=u]-u
t = table(s)
cov.t = (v.t[v.t>=start.year]-min(v.t[v.t>=start.year]))/max(v.t[v.t>=start.year]-min(v.t[v.t>=start.year]))
cov.t = cov.t[data>=u]
stopifnot(length(cov.t)==length(s))
length(s) # 435 in Tippett

# Kolgomorov Smirnov setting
number.of.bin = 5 # to compute KS test with temporal trend: nbr of bins to divide KS 
simulated.p.value = T
B.sim = 200 # nbr of MC simulation for p value
ks.bootstrap = T # if compute average of ks test on bootstrapped data 
B.bootstrap = 200

Clauset=F

B = 200 

# Test: same nbr of extreme outcomes per year
len.year = v.t[length(v.t)]-initial.date+1
stopifnot(length(table(cov.t))==len.year)
chisq.test(table(cov.t))
chisq.test(table(cov.t),p=rep(1/len.year,len.year)) # ATT: rejected... different in article

y00=as.numeric(names(table(cov.t)))[sample.int(n=len.year,size=length(s),replace=T)]
chisq.test(cov.t,y00)


par = val = cri = par.pm = ks = ks.t = ks.boot = ks.boot.t = latexName = B.simulations = list(NA)
for(name in models) par[[name]] = list(NA) # initialization

# number of tied observations with less than 20 ties
round(sum(t[t<=20])/sum(t),3)

pdf(file=paste0("images/","torn","_plotTable.pdf"),pointsize=pointsize,width=7,height=7)
par(pty="s")
plot(table(data[data>=u]),xlab="Number of extreme tornadoes per outbreak",ylab="Frequency",xaxt="n")
axis(side = 1, at = c(1,10,100,180))
dev.off()

pdf(file=paste0("images/","torn","large","_plotTable.pdf"),pointsize=pointsize.large,width=7,height=7)
par(pty="s")
plot(table(data[data>=u]),xlab="Number of extreme tornadoes per outbreak",ylab="Frequency",xaxt="n")
axis(side = 1, at = c(1,10,100,180))
dev.off()

# Temporal Discrete generalized Pareto 

name = "tdgpd"
if(any(name==models)){
  
  nll = function(par){
    sum = 0
    for(i in 1:length(s)){
      sum =  sum + log(pm.dgpd(s[i],xi=par[1],si=par[2]+par[3]*cov.t[i])) 
    }
    return(-sum)
  }
  
  o = optim(par=c(1,1,0),fn=nll,hessian=T)  
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value,nbr.par=length(o$par),nbr.obs=length(s))
  par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
  
  
  # Compute KS test
  for(iter in c(1,2)){
    if(iter==1) nbr = 1 
    if(iter==2) nbr = number.of.bin 
    
    seg = seq(from=0,to=1,length.out=nbr+1)
    ks.bin = seg.nbr.obs = rep(NA,length(seg)-1)
    seg[length(seg)]=seg[length(seg)]+0.001
    
    for(i in 1:(length(seg)-1)){
      w = (cov.t>=seg[i])&(cov.t<seg[i+1])# select variable with similar time covariate
      stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3]),a=0,b=max(s)+1)
      ks.bin[i]  = dgof::ks.test(x=s[w],y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
      seg.nbr.obs[i] = sum(w)
    }
    if(F) seg.nbr.obs
    if(F) ks.bin
    
    if(iter==1) ks[[name]] = ks.bin
    if(iter==2) ks.t[[name]] = min(ks.bin)
  }
  
  if(ks.bootstrap){ 
    print("KS computed using bootstrap")
    
    # Compute KS test for bootstrap
    for(iter in c(1,2)){
      if(iter==1) nbr = 1 
      if(iter==2) nbr = number.of.bin 
      
      seg = seq(from=0,to=1,length.out=nbr+1)
      ks.bin = seg.nbr.obs = rep(NA,length(seg)-1)
      seg[length(seg)]=seg[length(seg)]+0.001
      
      for(i in 1:(length(seg)-1)){
        print_loading(it=i,it.max=(length(seg)-1))
        
        w = (cov.t>=seg[i])&(cov.t<seg[i+1])# select variable with similar time covariate
        
        v = rep(NA,B.bootstrap)
        for(j in 1:B.bootstrap){
          stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3]),a=0,b=max(s)+1)
          v[j] = dgof::ks.test(x=sample(s[w],replace=T),y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
        }
        ks.bin[i] = mean(v)
        seg.nbr.obs[i] = sum(w)
      }
      if(F) seg.nbr.obs
      if(F) ks.bin
      
      if(iter==1) ks.boot[[name]] = ks.bin
      if(iter==2) ks.boot.t[[name]] = min(ks.bin)
    }
  }
  
  
  
  if(Clauset){ # KS: Method of A. Clauset and C. R. Shalizi and M. E. J. Newman
    
    nbr = 1
    B.Clauset=200
    
    seg = seq(from=0,to=1,length.out=nbr+1)
    ks.bin = seg.nbr.obs = rep(NA,length(seg)-1)
    seg[length(seg)]=seg[length(seg)]+0.001
    
    for(i in 1:(length(seg)-1)){
      print_loading(it=i,it.max=(length(seg)-1))
      w = (cov.t>=seg[i])&(cov.t<seg[i+1])# select variable with similar time covariate
      
      v = rep(NA,B.Clauset)
      
      stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3]),a=0,b=max(s)+1)
      ks.bench  = dgof::ks.test(x=s[w],y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
      
      for(j in 1:B.Clauset){
        # print(j)
        # simulate from fitted model
        x.sub = r.dgpd(length(s[w]),xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3])
        t.sub = table(x.sub)
        
        nll.sub = function(par) -sum(log(sapply(as.numeric(names(t.sub)),pm.dgpd,xi=par[1],si=par[2]))*t.sub) 
        o.sub = optim(par=c(1,1),fn=nll.sub,hessian=T)  
        
        stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o.sub$par[1],si=o.sub$par[2]),a=0,b=max(x.sub)+1)
        v[j] =dgof::ks.test(x=x.sub,y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
        
      }
      ks.bin[i] = 1-sum(v>=ks.bench)/length(v)
    }
    ks.bin
    min(ks.bin)
  }
  
  
  stopifnot(length(cov.t)==length(s))
  
  # Generate nbr of extreme outcomes in each year randomly
  cov.t.sampled = as.numeric(names(table(cov.t)))[sample.int(n=len.year,size=length(s),replace=T)]
  B.simulations.mat = matrix(NA,B,length(s))
  for(i in 1:length(s)){
    B.simulations.mat[,i] = r.dgpd(B,xi=o$par[1],si=o$par[2]+o$par[3]*cov.t.sampled[i])
  }
  stopifnot(sum(is.na(B.simulations.mat))==0)
  B.simulations[[name]] = B.simulations.mat
  
  
  
  
  #pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
  #par(pty="s")
  #QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,B.simulations.mat=B.simulations.mat,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
  # todo: include temporal dep
  #dev.off()
  if(F){
    pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
    # todo: include temporal dep
    dev.off()
    
    pdf(file=paste0("images/qqplot99_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,pmax=0.99,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
    # todo: include temporal dep
    dev.off()
  }
  latexName[[name]] = "$\\text{D-GPD}_{t}$"
}

# Discrete generalized Pareto 

name = "dgpd"
if(any(name==models)){
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.dgpd,xi=par[1],si=par[2]))*t) 
  o = optim(par=c(1,1),fn=nll,hessian=T)  
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value,nbr.par=length(o$par),nbr.obs=length(s))
  par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
  stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s)+1)
  ks[[name]] = dgof::ks.test(x=s,y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
  
  if(ks.bootstrap){
    v = rep(NA,B.bootstrap)
    for(i in 1:B.bootstrap){
      v[i] = dgof::ks.test(x=sample(s,replace=T),y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
    }
    ks.boot[[name]] = mean(v)
  }
  
  
  if(F){
    pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
    dev.off()
    
    pdf(file=paste0("images/qqplot99_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,pmax=0.99,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
    dev.off()
    if(F) QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=1,ci=T,pmax=0.9,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
  }
  latexName[[name]] = "D-GPD"
  
}


# Zipf-Mandelbrot
name = "zm"
if(any(name==models)){
  # reparametrization: (s,q)=(1+1/xi,si/xi)
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.zm.un,s=1+1/par[1],q=par[2]/par[1]))*t) + length(s)*log(pm.zm.norm(s=1+1/par[1],q=par[2]/par[1]))
  
  o = NULL; o = optim(par=c(1,1),fn=nll,hessian=T)  
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value,nbr.par=length(o$par),nbr.obs=length(s))
  par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
  stepFun = from.pdist.to.stepfun(pdist=function(x) p.czm(x,s=1+1/o$par[1],q=o$par[2]/o$par[1]),a=0,b=max(s)+1)
  ks[[name]] = dgof::ks.test(x=s+1,y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
  
  if(F){ # QQ-plots
    QQplot(s, function(x) q.czm(x, s=1+1/o$par[1],q=o$par[2]/o$par[1]),a=0,b=0.9,prop=0.05,B=B,lty=2)
    QQplot(s, function(x) q.czm(x, s=1+1/o$par[1],q=o$par[2]/o$par[1]),a=0,b=0.9,prop=0.05,pmax=0.99,B=B,lty=2)
  }
  
  latexName[[name]] = "ZM"
}

# GPD 
name = "gpd"
if(any(name==models)){
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),d.gpd,xi=par[1],si=par[2]))*t) 
  o = try(optim(par=c(1,1),fn=nll,hessian=T),silent=T)
  if(is.character(o)==F){ 
    o$value
    
    val[[name]] = o$value
    nll_discretized = function(par) -sum(log(sapply(as.numeric(names(t)),pm.dgpd,xi=o$par[1],si=o$par[2]))*t) 
    cri[[name]] = bic(nllh=nll_discretized(o$par),nbr.par=length(o$par),nbr.obs=length(s))
    par[[name]] =  get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
    # fitted distribution is for continuous rv; following test assumes that rounding the rv gives the data
    # equivalent to test that the data -1/2 are floored function of the fitted distribution 
    stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s))
    ks[[name]] = dgof::ks.test(x=s-1/2,y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
    
    if(F){
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,B=B,lty=2)
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.99,B=B,lty=2)
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.9,B=B,lty=2)
    }
  }
  latexName[[name]] = "$\\text{GPD}_{\\delta=0}$"
  
  
  
}


name = "tgpd"
if(any(name==models)){
  
  nll = function(par){
    sum = 0
    for(i in 1:length(s)){
      sum =  sum + log(d.gpd(s[i],xi=par[1],si=par[2]+par[3]*cov.t[i])) 
    }
    return(-sum)
  }
  o = try(optim(par=c(1,1,0),fn=nll,hessian=T),silent=T)
  
  if(is.character(o)==F){ 
    o$value
    
    nll_discretized = function(par){
      sum = 0
      for(i in 1:length(s)){
        sum =  sum + log(pm.dgpd(s[i],xi=par[1],si=par[2]+par[3]*cov.t[i])) 
      }
      return(-sum)
    }
    
    val[[name]] = nll_discretized(o$par)
    cri[[name]] = bic(nllh=nll_discretized(o$par),nbr.par=length(o$par),nbr.obs=length(s))
    par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
    # fitted distribution is for continuous rv; following test assumes that rounding the rv gives the data
    # equivalent to test that the data -1/2 are floored function of the fitted distribution 
    
    if(F){
      # todo: add temporal 
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,B=B,lty=2)
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.99,B=B,lty=2)
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.9,B=B,lty=2)
    }
  }
  
  latexName[[name]] = "$\\text{GPD}_{\\delta=0,t}$"
  
  
  # Compute KS test
  for(iter in c(1,2)){
    if(iter==1) nbr = 1 
    if(iter==2) nbr = number.of.bin 
    
    seg = seq(from=0,to=1,length.out=nbr+1)
    ks.bin = seg.nbr.obs = rep(NA,length(seg)-1)
    seg[length(seg)]=seg[length(seg)]+0.001
    
    for(i in 1:(length(seg)-1)){
      w = (cov.t>=seg[i])&(cov.t<seg[i+1])# select variable with similar time covariate
      stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3]),a=0,b=max(s)+1)
      ks.bin[i]  = dgof::ks.test(x=s[w],y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
      seg.nbr.obs[i] = sum(w)
    }
    if(F) seg.nbr.obs
    if(F) ks.bin
    
    if(iter==1) ks[[name]] = ks.bin
    if(iter==2) ks.t[[name]] = min(ks.bin)
  }
  
  if(ks.bootstrap){
    v = rep(NA,B.bootstrap)
    for(i in 1:B.bootstrap){
      v[i] = dgof::ks.test(x=sample(s,replace=T),y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
    }
    ks.boot[[name]] = mean(v)
  }
  
}


# GPD -1/2
name = "gpd2"
if(any(name==models)){
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),d.gpd,xi=par[1],si=par[2],mu=-1/2))*t) 
  o = try(optim(par=c(1,1),fn=nll,hessian=T),silent=T)
  if(is.character(o)==F){
    o$value
    
    nll_discretized = function(par) -sum(log(sapply(as.numeric(names(t)),pm.dgpd,xi=o$par[1],si=o$par[2]))*t) 
    val[[name]] = nll_discretized(o$par)    
    cri[[name]] = bic(nllh=nll_discretized(o$par),nbr.par=length(o$par),nbr.obs=length(s))
    par[[name]] =  get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
    
    # fitted distribution is for continuous rv; following test assumes that rounding the rv gives the data
    # equivalent to test that the data are floored function of the fitted distribution but with mu=0 instead of mu=-1/2
    stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s))
    ks[[name]] = dgof::ks.test(x=s,y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
    if(dataset=="X15"){
      pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
      par(pty="s")
      QQplot(s-1/2, function(x) q.gpd(x,xi=o$par[1],si=o$par[2],mu=-1/2),a=0,b=0.3,prop=0.1,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
      dev.off()
      
      if(F) QQplot(s-1/2, function(x) q.gpd(x,xi=o$par[1],si=o$par[2],mu=-1/2),pmax=0.9,a=0,b=0.3,prop=0.1,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
    }
  }
  latexName[[name]] = "$ \\text{GPD}_{\\delta=-\\frac{1}{2}} $"
}




# GPD -1/2
name = "tgpd2"
if(any(name==models)){
  
  nll = function(par){
    sum = 0
    for(i in 1:length(s)){
      sum =  sum + log(d.gpd(s[i],xi=par[1],si=par[2]+par[3]*cov.t[i],mu=-1/2)) 
    }
    return(-sum)
  }
  o = try(optim(par=c(1,1,0),fn=nll,hessian=T),silent=T)
  
  if(is.character(o)==F){
    
    
    nll_discretized = function(par){
      sum = 0
      for(i in 1:length(s)){
        sum =  sum + log(pm.dgpd(s[i],xi=par[1],si=par[2]+par[3]*cov.t[i])) 
      }
      return(-sum)
    }
    
    val[[name]] = nll_discretized(o$par)
    cri[[name]] = bic(nllh=nll_discretized(o$par),nbr.par=length(o$par),nbr.obs=length(s))
    par[[name]] =  get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
    
    
    if(dataset=="X15"){
      pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
      par(pty="s")
      QQplot(s-1/2, function(x) q.gpd(x,xi=o$par[1],si=o$par[2],mu=-1/2),a=0,b=0.3,prop=0.1,B=B,lty=2,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
      dev.off()
    }
  }
  
  # Compute KS test
  for(iter in c(1,2)){
    if(iter==1) nbr = 1 
    if(iter==2) nbr = number.of.bin 
    
    seg = seq(from=0,to=1,length.out=nbr+1)
    ks.bin = seg.nbr.obs = rep(NA,length(seg)-1)
    seg[length(seg)]=seg[length(seg)]+0.001
    
    for(i in 1:(length(seg)-1)){
      w = (cov.t>=seg[i])&(cov.t<seg[i+1])# select variable with similar time covariate
      stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3]),a=0,b=max(s)+1)
      ks.bin[i]  = dgof::ks.test(x=s[w],y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
      seg.nbr.obs[i] = sum(w)
    }
    if(F) seg.nbr.obs
    if(F) ks.bin
    
    if(iter==1) ks[[name]] = ks.bin
    if(iter==2) ks.t[[name]] = min(ks.bin)
  }
  
  if(ks.bootstrap){ 
    
    # Compute KS test for bootstrap
    for(iter in c(1,2)){
      if(iter==1) nbr = 1 
      if(iter==2) nbr = number.of.bin 
      
      seg = seq(from=0,to=1,length.out=nbr+1)
      ks.bin = seg.nbr.obs = rep(NA,length(seg)-1)
      seg[length(seg)]=seg[length(seg)]+0.001
      
      for(i in 1:(length(seg)-1)){
        print_loading(it=i,it.max=(length(seg)-1))
        
        w = (cov.t>=seg[i])&(cov.t<seg[i+1])# select variable with similar time covariate
        
        v = rep(NA,B.bootstrap)
        for(j in 1:B.bootstrap){
          stepFun = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]+mean(cov.t[w])*o$par[3]),a=0,b=max(s)+1)
          v[j] = dgof::ks.test(x=sample(s[w],replace=T),y=stepFun,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
        }
        ks.bin[i] = mean(v)
        seg.nbr.obs[i] = sum(w)
      }
      if(F) seg.nbr.obs
      if(F) ks.bin
      
      if(iter==1) ks.boot[[name]] = ks.bin
      if(iter==2) ks.boot.t[[name]] = min(ks.bin)
    }
  }
  
  
  latexName[[name]] = "$ \\text{GPD}_{t,\\delta=-\\frac{1}{2}} $"
}


# Test if temporal trend
if(is.element("dgpd",models)&is.element("tdgpd",models)){
  D= 2 *(val[["dgpd"]] - val[["tdgpd"]] ) 
  p.value = 1-pchisq(D,df=3-2) 
  if(p.value<0.05) print("Null hyp stationarity rejected")
}


# QQ Plots from simulations
for(name in models.qqplot){
  if(is.null(B.simulations[[name]])==F){
    n=length(s)
    seq=seq(1,n)/(n+1) 
    q.exp = quantile(s,seq,type=1) 
    q.mean = q.u = q.l = rep(NA,length(seq))
    for(k in 1:n){
      v = apply(B.simulations[[name]],1,quantile,probs=seq[k],type=1)   
      q.mean[k] = mean(v)
      q.u[k] = quantile(v,0.975)
      q.l[k] = quantile(v,0.025)
    }
    
    pdf(file=paste0("images/qqplotSim_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    plot(q.exp,q.mean,xlab="Empirical Quantiles",ylab="Theoretical Quantiles")
    abline(0,1)
    lines(q.exp,q.u,lty=2)
    lines(q.exp,q.l,lty=2)
    dev.off()
  }
}


# =======
# EXPORT RESUlTS IN A LATEX TABLE

# print fixed number of decimals
prt.rd = function(x,rd) sprintf(paste0("%.",as.character(rd),"f"),x)
# build confidence intervals in latex
ci.in.latex= function(par,par.d,par.u,rd) as.character(paste0("$", prt.rd(par,rd), "_{[",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"]}","$"))


txt = "Model &  BIC & $ \\xi$ & $ \\sigma_0 $  & $ \\sigma_1 $   \\\\ \n"

for(name in models){
  if(all(is.na(par[[name]]))==F){
    txt = paste(txt, latexName[[name]], "& $",
                #                prt.rd(ks[[name]],rd=2), "$ & $", 
                prt.rd(cri[[name]],rd=1), "$ &", 
                ci.in.latex(par[[name]][1,1],par[[name]][1,2],par[[name]][1,3],rd=2), "&", 
                ci.in.latex(par[[name]][2,1],par[[name]][2,2],par[[name]][2,3],rd=2), "&", 
                ci.in.latex(par[[name]][3,1],par[[name]][3,2],par[[name]][3,3],rd=2), 
                "\\\\","\n")
  } else{
    txt = paste(txt,latexName[[name]], "& NA &  & &\\\\ \n")
  }
}  

cat(txt)
cat(paste0("Threshold = ",u,", quantile:",qt,", size of exceed. = ",length(s)),", size init data = ",length(data),"\n")

# other summary 
round(unlist(ks),2)
round(unlist(ks.t),2)
round(unlist(ks.boot),2)
round(unlist(ks.boot.t),2)

