
myRFolder = "~/Downloads/discrete_extremes"
setwd(myRFolder)
source("Functions.R")
source("utilities_5Feb2017.R")


# Import data
raw.data.FR = read.csv("Insee_AccouchMultiples_1995_2014.csv")
raw.data.US = read.csv("USMultipleBirths_1995_2014.csv",sep=" ")


# Process data
data.US = as.table(apply(raw.data.US[,c("single","twin","triplet","quadruplet","quintupletOrMore")],2,sum))
data.FR.2 = as.table(apply(raw.data.FR[,c("totBirths","twin","triplet","quadruplet","quintupletOrMore")],2,sum))
data.FR = c(data.FR.2[1]-sum(data.FR.2[-1]), data.FR.2[-1])
names(data.FR)=c("single","twin","triplet","quadruplet","quintupletOrMore")
stopifnot(sum(data.FR)==sum(raw.data.FR[,"totBirths"]))
data.FR.red = c(data.FR[1:3],sum(data.FR[4:5]))
names(data.FR.red)=c("single","twin","triplet","quadrupletOrMore")

# Display data
print(data.US)
print(data.FR)


# ====================
# ====================

# MODEL COMPARISON 

# ====================
# ====================


# ====================

# USER SETTINGS

# choice of models
name.list = c("gpd_disc","zm.rp","nb","pois")

# choice of data set

dataset="FR"
dataset="US"


# 100*(1-alpha)% confidence interval 
alpha=0.1

# graphical settings
pointsize = 18

# ====================


# select data set
data = NULL
if(dataset=="US"){
  data = data.US
} else if (dataset=="FR"){
  data = data.FR
}

# threshold 
u= 3 # try u=2 and u=1
m = 5-u # censored point: quintuplets and more

# select exceedances above u
t = data[u:length(data)]
v0 = s0 = NULL 
for(i in 1:length(t)){
  v0 = c(v0,i-1)
  s0 = c(s0,rep(i-1,t[i]))
}
names(t)= v0
s = s0 
stopifnot(length(s)==sum(t))

# Probability to be great or equal to the threshold u
pu = sum(t)/sum(data)

routine = NULL

# target: try to estimate Pr(X-u>= qe)
qe = NULL
if(dataset=="US"){
  qe = m
} else if (dataset=="FR"){
  qe = m-1
}

# set lists
val  = cri = par  = pe = pe.pm = chisq = list(NA)
for(name in name.list) par[[name]] = list(NA)

# Fit a range of models

# D-GPD
name = "gpd_disc"
if(any(name==name.list)){
  
  pm.dgpd.right.trunc = function(x,xi,si) pmf.rigth.trunc.dist(x,pmf=function(x) pm.dgpd(x,xi=xi,si=si),cdf=function(x) p.dgpd(x,xi=xi,si=si),lower=0,upper=m)
  nll = function(par) -sum(log(sapply(as.numeric(names(t)), pm.dgpd.right.trunc,xi=par[1],si=par[2]))*t) 
  # check if correct
  
  o = NULL; o = optim(par=c(1,1),fn=nll,hessian=T)  
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value, nbr.par=length(o$par),nbr.obs=sum(t))
  par[[name]] = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=F)
  
  # chi-squared test
  p.val = chisq.test(x=t,p=sapply(0:m,function(x) pm.dgpd.right.trunc(x,xi=o$par[1],si=o$par[2])))$p.value
  chisq[[name]] = round(p.val,2)
  
  if(F){
    pdf(file=paste0("images/","_",dataset,"_",routine,"_",name,".pdf"),pointsize=13,width=7,height=7,family="Bookman")
    par(pty="s")
    QQplot(s[s<=m],function(x) q.tgpd(x,xi=o$par[1],si=o$par[2],m=m-1),a=0,b=0.999,prop=0.05)
    dev.off()
    # warning: q.tgpd may be incorrect
  }
  
  # estimate and confidence interval for prob(X>=qe) 
  pr = pu*(1-p.dgpd(qe-1,xi=o$par[1],si=o$par[2]))
  grad.pr = get.grad.sprob.gpd(qe,pu=pu,xi=o$par[1],si=o$par[2])
  var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
  pr.pm = as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
  pu*sum(t[(qe+1):length(t)])/sum(t); pr-pr.pm; pr+pr.pm
  
  pe[[name]]= pr
  pe.pm[[name]] = pr.pm
  
}


# Reparametrized Zipf-Mandelbrot distribution
name = "zm.rp"
if(is.element(name,name.list)){
  pm.zm.right.trunc = function(x,s,q) pmf.rigth.trunc.dist(x,pmf=function(x) pm.zm(x,s=s,q=q),cdf=function(x) p.zm(x,s=s,q=q),lower=0,upper=m)
  nll = function(par) -sum(log(sapply(as.numeric(names(t)), pm.zm.right.trunc,s=1+1/par[1],q=par[2]/par[1]))*t) 
  
  
  o = optim(par=c(1,1),fn=nll,hessian=T)
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value, nbr.par=length(o$par),nbr.obs=sum(t))
  par[[name]] = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=F)
  p.val = chisq.test(x=t,p=sapply(0:m,function(x) pm.zm.right.trunc(x,s=1+1/o$par[1],q=o$par[2]/o$par[1])))$p.value
  chisq[[name]] = round(p.val,2)
  
  #  # prob(X>=qe)
  pr = pu*(1-p.zm(qe-1,s=1/o$par[1]+1,q=o$par[2]/o$par[1])) 
  grad.pr = numDeriv::grad(f= function(y) pu*(1-p.zm(qe-1,s=1/y[1]+1,q=y[2]/y[1])),c(o$par[1],o$par[2]))
  
  var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
  pr.pm= as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
  pu*sum(t[(qe+1):length(t)])/sum(t); pr-pr.pm; pr+pr.pm
  
  pe[[name]] = pr
  pe.pm[[name]] = pr.pm
}


# NB
name = "nb"
if(any(name==name.list)){
  
  dnb.right.trunc = function(x,size,prob) pmf.rigth.trunc.dist(x,pmf=function(x) dnbinom(x,size=size,prob=prob),cdf=function(x) pnbinom(x,size=size,prob=prob),lower=0,upper=m)
  nll = function(par) -sum(log(sapply(as.numeric(names(t)), dnb.right.trunc,size=par[1],prob=par[2]))*t) 
  
  o = NULL; o = optim(par=c(1,0.01),fn=nll,hessian=T)  
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value, nbr.par=length(o$par),nbr.obs=sum(t))
  par[[name]] = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=F)
  
  p.val = chisq.test(x=t,p=sapply(0:m,function(x) dnb.right.trunc(x,size=o$par[1],prob=o$par[2])))$p.value
  chisq[[name]] = round(p.val,2)
  
  #  # prob(X>=qe)
  pr = pu*(1-pnbinom(qe-1,size=o$par[1],prob=o$par[2])) 
  grad.pr = numDeriv::grad(f= function(par) pu*(1-pnbinom(qe-1,size=par[1],prob=par[2])),c(o$par[1],o$par[2]))
  
  var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
  pr.pm= as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
  pu*sum(t[(qe+1):length(t)])/sum(t); pr-pr.pm; pr+pr.pm
  
  pe[[name]] = pr
  pe.pm[[name]] = pr.pm
  
}



# Poisson
name = "pois"
if(any(name==name.list)){
  
  dpois.right.trunc = function(x,lambda) pmf.rigth.trunc.dist(x,pmf=function(x) dpois(x,lambda=lambda),cdf=function(x) ppois(x,lambda=lambda),lower=0,upper=m)
  nll = function(par) -sum(log(sapply(as.numeric(names(t)), dpois.right.trunc,lambda=par))*t) 
  
  o = NULL; o = optim(par=c(1),fn=nll,hessian=T)  
  val[[name]] = o$value
  cri[[name]] = bic(nllh=o$value, nbr.par=length(o$par),nbr.obs=sum(t))
  par[[name]] = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=F)
  
  p.val = chisq.test(x=t,p=sapply(0:m,function(x) dpois.right.trunc(x,lambda=o$par)))$p.value
  chisq[[name]] = round(p.val,2)
  
  #  # prob(X>=qe)
  pr = pu*(1-ppois(qe-1,lambda=o$par)) 
  grad.pr = numDeriv::grad(f= function(par) pu*(1-ppois(qe-1,lambda=par)),o$par)
  
  var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
  pr.pm= as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
  pu*sum(t[(qe+1):length(t)])/sum(t); pr-pr.pm; pr+pr.pm
  
  pe[[name]] = pr
  pe.pm[[name]] = pr.pm
  
}

# Display results

prt.rd = function(x,rd) sprintf(paste0("%.",as.character(rd),"f"),x)
ci.in.latex.1 = function(par,par.d,par.u,rd) as.character(paste0("$", prt.rd(par,rd), "_{[",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"]}","$"))
ci.in.latex.2 = function(par.d,par.u,rd) as.character(paste0("$ [",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"] $"))

for(name in name.list){
  
  rd = 2
  cat(paste(name, " & $",prt.rd(val[[name]],1),"$ &", ci.in.latex.1(par[[name]][1,1],par[[name]][1,2],par[[name]][1,3],rd), "&", 
            ci.in.latex.1(par[[name]][2,1],par[[name]][2,2],par[[name]][2,3],rd), "&",
            "$",prt.rd(chisq[[name]],rd),"$ \\\\","\n"))
  
}

cat(paste0("Threshold = ",u,", prob above: ",qe,", size of exceed. = ",length(s)),", size init data = ",length(data),"\n")

print(unlist(cri))





# ====================
# EXTREME QUANTILE ESTIMATION

set.seed(199)

B=500 

# data
data = data.FR
data = data.US


name.list = c("gpd_disc","nb","pois")
name.list = c("gpd_disc","zm.rp","nb")
name.list = c("gpd_disc","nb","pois")
name.list = c("zm.rp")
name.list = c("gpd_disc")
name.list = c("gpd_disc","zm.rp","nb","pois")
# threshold 

u= 2
m = 5-u # censored point on exceedances
qe = 5-u

# exceedances
t = data
names(t)= c("1","2","3","4","5")

pe.true = t[5]/sum(t)

# nbr of obs the data used in estimation
K = round(sum(t)/1000) # try with 1000 # TODO: selecton not well done ->>_-- wrong here

nbr.exced = nbr.4 = nbr.5 = rep(NA,B)
pe = pe.pm = val = list(NA)

for(i in 1:B){
  print_loading(it=i,it.max=B)
  
  x.sub = sample(1:length(t),size=K,replace=T,prob=t)   # sample K observations from data distribution
  pu = length(x.sub[x.sub>=u]) / length(x.sub) # probability of being above the threshold
  s.sub = x.sub[x.sub>=u]-u
  t.sub = table(s.sub)
  plot(t.sub)  
  nbr.exced[i] = length(s.sub)
  nbr.4[i] = length(s.sub[s.sub>=4-u])
  nbr.5[i] = length(s.sub[s.sub>=5-u]) # number of observations in the target set
  
  
  # D-GPD
  name = "gpd_disc"
  if(is.element(name,name.list)){
    
    pm.dgpd.right.trunc = function(x,xi,si) pmf.rigth.trunc.dist(x,pmf=function(x) pm.dgpd(x,xi=xi,si=si),cdf=function(x) p.dgpd(x,xi=xi,si=si),lower=0,upper=m) 
    nll = function(par) -sum(log(sapply(as.numeric(names(t.sub)), pm.dgpd.right.trunc,xi=par[1],si=par[2]))*t.sub) 
    
    o = NULL; o = optim(par=c(1,1),fn=nll,hessian=T)  
    val[[name]][i] = o$value
    
    if(F) p.val = chisq.test(x=t.sub,p=sapply(0:m,function(x) pm.dgpd.right.trunc(x,xi=o$par[1],si=o$par[2])))$p.value
    
    # estimate and confidence interval for prob(X>=qe) 
    pr = pu*(1-p.dgpd(qe-1,xi=o$par[1],si=o$par[2]))
    
    grad.pr = numDeriv::grad(f= function(par) pu*(1-p.dgpd(qe-1,xi=par[1],si=par[2])),c(o$par[1],o$par[2]))
    if(F) grad.pr = get.grad.sprob.gpd(qe-1,pu=pu,xi=o$par[1],si=o$par[2])
    
    var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
    pr.pm = as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
    t[5]/sum(t); pr-pr.pm; pr+pr.pm
    
    pe[[name]][i]= pr
    pe.pm[[name]][i] = pr.pm
  }
  
  # Reparametrized Zipf-Mandelbrot distribution
  name = "zm.rp"
  if(is.element(name,name.list)){
    
    pm.zm.right.trunc = function(x,s,q) pmf.rigth.trunc.dist(x,pmf=function(x) pm.zm(x,s=s,q=q),cdf=function(x) p.zm(x,s=s,q=q),lower=0,upper=m)
    nll = function(par) -sum(log(sapply(as.numeric(names(t.sub)), pm.zm.right.trunc,s=1+1/par[1],q=par[2]/par[1]))*t.sub) 
    
    
    # catch numerical error when computing hessian
    o = try( optim(par=c(1,1),fn=nll,hessian=T) ,silent=T)
    
    if(is.character(o)){ # if couldn't compute the hessian
      warning(o[1]); 
      o = try( optim(par=c(1,0.1),fn=nll,hessian=F), silent=T)
      val[[name]][i] = o$value
      #  # prob(X>=qe)
      pr = pu*(1-p.zm(qe-1,s=1/o$par[1]+1,q=o$par[2]/o$par[1]))  # TODO problem here
      pe[[name]][i] = pr
      pe.pm[[name]][i] = NA
      
    } else{
      val[[name]][i] = o$value
      
      #  # prob(X>=qe)
      pr = pu*(1-p.zm(qe-1,s=1/o$par[1]+1,q=o$par[2]/o$par[1])) 
      pe[[name]][i] = pr
      
      grad.pr = numDeriv::grad(f= function(y) pu*(1-p.zm(qe-1,s=1/y[1]+1,q=y[2]/y[1])),c(o$par[1],o$par[2]))
      var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
      pr.pm= as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
      t[5]/sum(t);  pr-pr.pm; pr+pr.pm
      pe.pm[[name]][i] = pr.pm
      
      if(F) p.val = chisq.test(x=t.sub,p=sapply(0:m,function(x) pm.zm.right.trunc(x,s=1+1/o$par[1],q=o$par[2]/o$par[1])))$p.value
    }
    #}
  }
  
  # Poisson
  name = "pois"
  if(any(name==name.list)){
    
    dpois.right.trunc = function(x,lambda) pmf.rigth.trunc.dist(x,pmf=function(x) dpois(x,lambda=lambda),cdf=function(x) ppois(x,lambda=lambda),lower=0,upper=m)
    nll = function(par) -sum(log(sapply(as.numeric(names(t.sub)), dpois.right.trunc,lambda=par))*t.sub) 
    
    o = NULL; o = optim(par=c(1),fn=nll,hessian=T)  
    val[[name]][i] = o$value
    
    if(F) p.val = chisq.test(x=t.sub,p=sapply(0:m,function(x) dpois.right.trunc(x,lambda=o$par)))$p.value
    
    #  # prob(X>=qe)
    pr = pu*(1-ppois(qe-1,lambda=o$par)) 
    grad.pr = numDeriv::grad(f= function(par) pu*(1-ppois(qe-1,lambda=par)),o$par)
    
    var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
    pr.pm= as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
    t[5]/sum(t); pr-pr.pm; pr+pr.pm
    
    pe[[name]][i] = pr
    pe.pm[[name]][i] = pr.pm
    
  }
  
  # NB
  name = "nb"
  if(any(name==name.list)){
    
    dnb.right.trunc = function(x,size,prob) pmf.rigth.trunc.dist(x,pmf=function(x) dnbinom(x,size=size,prob=prob),cdf=function(x) pnbinom(x,size=size,prob=prob),lower=0,upper=m)
    nll = function(par) -sum(log(sapply(as.numeric(names(t.sub)), dnb.right.trunc,size=par[1],prob=par[2]))*t.sub) 
    
    # catch numerical error when computing hessian
    o = try( optim(par=c(1,0.01),fn=nll,hessian=T)  ,silent=T)
    if(is.character(o)){ # if couldn't compute the hessian
      warning(o[1]); 
      o =  optim(par=c(1,0.01),fn=nll,hessian=F)
      val[[name]][i] = o$value
      pr = pu*(1-pnbinom(qe-1,size=o$par[1],prob=o$par[2])) 
      pe[[name]][i] = pr
      pe.pm[[name]][i] = NA
      
    } else {
      #  # prob(X>=qe)
      val[[name]][i] = o$value
      pr = pu*(1-pnbinom(qe-1,size=o$par[1],prob=o$par[2])) 
      grad.pr = numDeriv::grad(f= function(par) pu*(1-pnbinom(qe-1,size=par[1],prob=par[2])),c(o$par[1],o$par[2]))
      
      var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
      pr.pm= as.numeric(from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T))
      pu*sum(t[(qe+1):length(t)])/sum(t); pr-pr.pm; pr+pr.pm
      
      pe[[name]][i] = pr
      pe.pm[[name]][i] = pr.pm
    }
    
  }
  
  
}

mean(val$gpd_disc)
mean(val$nb)
mean(val$zm.rp)
mean(val$pois)
mean(nbr.exced)
mean(nbr.4)
mean(nbr.5)

pe.m = coverage = len = true.len = nbr.na = nbr.na.ci = list(NA)
pe.true 
for(name in name.list){
  
  w = is.na(pe[[name]])==F
  w.ci = is.na(pe.pm[[name]])==F
  nbr.na[[name]] = sum(length(w)-sum(w)) # nbr NA
  nbr.na.ci[[name]] =  sum(length(w.ci)-sum(w.ci)) # nbr NA in CI
  pe.m[[name]] = mean(pe[[name]][w])
  
  coverage[[name]] = sum( (pe[[name]][w.ci]-pe.pm[[name]][w.ci] <= pe.true)&(pe[[name]][w.ci]+pe.pm[[name]][w.ci] >= pe.true))/sum(w.ci)
  len[[name]] = mean(2*pe.pm[[name]][w.ci])
  true.len[[name]] = quantile(pe[[name]][w],1-alpha/2) - quantile(pe[[name]][w],alpha/2)
}

unlist(nbr.na)
unlist(nbr.na.ci)
round(unlist(pe.m),7)
round(unlist(coverage),2)
round(unlist(len),7)
round(unlist(true.len),7)



