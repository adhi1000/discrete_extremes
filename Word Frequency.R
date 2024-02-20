

# ==================================
# USER SETTINGS

# set directory


myRFolder = "~/Downloads/discrete_extremes"
setwd(myRFolder)
source("Functions.R")
source("utilities_5Feb2017.R")


# Choose data set code: 0,9,10 or 15
experiment = 15#0 #15, 0

# Quality of confidence intervals
B = 2000 # final version: 2000

# sim p.value ks test

sim.pvalue = NULL
ks.bootstrap = T
simulated.p.value=T
B.sim = 200 # nbr of simulations for p-value
B.bootstrap = 200 # nbr of resampling to compute KS



Clauset=F
# graphical settings
pointsize = 18
pointsize.large = 12



# ==================================

# Import data

if(exists("data.eng")==F) data.eng = read.table(file="British_word_freq.txt")
if(exists("data.fr")==F) data.fr = read.table(file="French_word_freq.txt")
x9 = round(sort(data.fr[,"X9_freqfilms2"],decreasing=T)) # Rounded French word frequencies in a collection of movie subtitles
x10 = round(sort(data.fr[,"X10_freqlivres"],decreasing=T)) # Rounded French word frequencies in a collection of books
x0 = data.eng[,"Freq"] # English word frequencies in a corpus
x15 = sort(data.fr[,"X15_nblettres"],decreasing=T) # length of French words

if(F){ # rounding
  x9 = 100*sort(data.fr[,"X9_freqfilms2"],decreasing=T) #  French word frequencies in a collection of movie subtitles
  x10 =100* sort(data.fr[,"X10_freqlivres"],decreasing=T) # French word frequencies in a collection of movie subtitles
}

x9.name = data.fr[,"X1_ortho"][sort(data.fr[,"X9_freqfilms2"],decreasing=T,index.return=T)$ix]
x10.name = data.fr[,"X1_ortho"][sort(data.fr[,"X10_freqlivres"],decreasing=T,index.return=T)$ix]
# Choose data set

if(experiment == 0){
  lang="eng"
  dataset = ""
  data = x0
  qt = NA
  if(F) u = min(data)+1 
  u= 16
  
  # u= 11 : pvalue DGPD:  0.09
  # u= 13 : pvalue DGPD:  0,18
  # u= 15 : pvalue DGPD: 0.34
  # u = 16 : pvalue DGPD: 0.56
  # u= 18 : pvalue DGPD: 0.56 (att ZM: 0.59, numerical diff)
  # u= 20 : pvalue DGPD: 0.66
  # u=60 : pvalue DGPD: 0.59. (att: 0.65 for gpd, mistake to correct)
  
  s = data[data>=u]-u
  t = table(s)
  
  sim.pvalue = F
  
  if(F){ 
    u.limit=306
    length(data[data>u.limit]) # at least 300 observations
    
    # Mean residual plot
    ismev::mrl.plot(data-1/2, umin = min(data), umax = u.limit,conf = 0.95, nint = 100) 
    ismev::mrl.plot(data-1/2, umin = min(data), umax = 60,conf = 0.95, nint = 100) # u=16: min point st linear
    
    # Fit over range of threshold
    ismev::gpd.fitrange(data, umin=min(data), umax=u.limit, nint = 50, show = FALSE) # u=60: min st all included
    ismev::gpd.fitrange(data, umin=min(data), umax=60, nint = 50, show = FALSE) #if u<50, then u=30
    
  }
  
}

if(experiment == 9){
  lang="fr"
  dataset = "X9" # word freq in subtitles
  qt = 0.95
  data = x9
  u = round(quantile(data,qt)) # select large threshold
  s=data[data>=u]-u  # get exceedances above the threshold
  t=table(s)
}

if(experiment == 10){
  lang="fr"
  dataset = "X10" # word freq in books
  qt = 0.95
  data = x10
  u = round(quantile(data,qt)) # select large threshold
  s=data[data>=u]-u  # get exceedances above the threshold
  t=table(s)
}

if(experiment == 15){
  lang="fr"
  dataset = "X15" # word lengths
  routine ="tail"
  qt = 0.98
  data = x15
  u = round(quantile(data,qt)) # select large threshold
  
  s=data[data>=u]-u  # get exceedances above the threshold
  t=table(s)
  
}


# Plot frequency table

if(lang=="eng"){
  pdf(file=paste0("images/",lang,"_",dataset,"_figure_1_right.pdf"),pointsize=pointsize,width=7,height=7)
  par(pty="s")
  plot(table(data[data>=u]),xlab="Word Frequency",ylab="Frequency",log="x",xlog=T,xaxt="n")
  axis(side = 1, at = c(1,10,100,1000,10000))
  dev.off()
}

if(experiment =="15"){
  article.table.reference = "_plotTable"
  if(lang=="fr" & dataset=="X15") article.table.reference = "_Figure_1_right"
  pdf(file=paste0("images/",lang,"_",dataset,article.table.reference,".pdf"),pointsize=pointsize,width=7,height=7)
  par(pty="s")
  plot(table(data[data>=u]),xlab="Word length frequency",ylab="Frequency",xaxt="n")
  axis(side = 1, at = pretty(as.numeric(names(table(data[data>=u])))))
  dev.off()
  
  sim.pvalue = T
}

# Fit the following models to the exceedances s by maximum likelihood:
# - gpd: generalized Pareto distribution (shape parameter xi, scale parameter sim and location parameter mu=0)
# - gpd2: generalized Pareto distribution with location mu=-1/2
# - dgpd: discrete generalized Pareto distribution
# - tdgpd: truncated generalized Pareto distribution (additional parameter: endpoint m)
# - zm: Zipf-Mandelbrot law (parameters s=1+1/xi, q= si/xi)

models = c("dgpd","zm","gpd2","gpd","nb")


# For each model, compute the following quantities:
# - par: maximum likelihood estimator 
# - val: log-likelihood at the mle 
# - par.pm: value pm such that [par-pm, par+pm] is a confidence interval for the mle (par.pm)
# - ks: Kolmogorov-Smirnov goodness-of-fit test for discrete data
# - latexName: name of the model in a latex tabular

par = val = cri = par.pm = ks = ks.sim =ks.clauset= ks.boot = stepFun = latexName = list(NA)
for(name in models) par[[name]] = list(NA) # initialization
for(name in models) ks.boot[[name]] = rep(NA,B.bootstrap)


# number of tied observations with less than 20 ties
round(sum(t[t<=20])/sum(t),3)

# Discrete generalized Pareto 

name = "dgpd"
if(any(name==models)){
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.dgpd,xi=par[1],si=par[2]))*t) 
  o = optim(par=c(1,1),fn=nll,hessian=T)  
  val[[name]] = o$value
  par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
  stepFun[[name]] = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s)+1)
  ks[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=F)$p.value # Kolgomorov Smirnov test
  if(sim.pvalue) ks.sim[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=T,B=B.sim)$p.value # Kolgomorov Smirnov test
  
  
  if(Clauset){ # KS: Method of A. Clauset and C. R. Shalizi and M. E. J. Newman
    
    B.Clauset=100
    v = rep(NA,B.Clauset)
    
    stepFun[[name]] = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s)+1)
    ks.bench = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=simulated.p.value,B=B.sim)$p.value # Kolgomorov Smirnov test
    
    for(j in 1:B.Clauset){
      print_loading(it=j,it.max=B.Clauset)
      # simulate from fitted model
      x.sub = r.dgpd(length(s),xi=o$par[1],si=o$par[2])
      t.sub = table(x.sub)
      
      nll.sub = function(par) -sum(log(sapply(as.numeric(names(t.sub)),pm.dgpd,xi=par[1],si=par[2]))*t.sub) 
      o.sub = optim(par=c(1,1),fn=nll.sub,hessian=T)  
      
      stepFunClauset = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o.sub$par[1],si=o.sub$par[2]),a=0,b=max(x.sub)+1)
      v[j] =dgof::ks.test(x=x.sub,y=stepFunClauset,alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
    }
    ks.clauset[[name]] = sum(ks.bench>=v)/length(v)
  }
  
  
  pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
  par(pty="s")
  QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
  dev.off()
  
  pdf(file=paste0("images/qqplot99_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
  par(pty="s")
  QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.75,prop=0.05,ci=T,pmax=0.98,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
  dev.off()
  
  if(dataset=="X15"){
    pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.3,prop=0.1,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
    abline(1,1)
    dev.off()
  }
  
  
  
  if(F){# Check if right-truncated discrete GPD distributon with a third parameter m is a better model
    nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.tgpd,xi=par[1],si=par[2],m=max(s)+par[3]))*t) 
    m0=max(s)*10; m0=10;
    o.trunc = optim(par=c(1,1,m0),fn=nll,hessian=T)
    o.trunc$par
    bic(o.trunc$value,3,nbr.obs=length(s))<bic(o$value,2,nbr.obs=length(s)) # improvement?
    # For experiment 0: no. For 9: yes (m small). For 10: no
  }  
  latexName[[name]] = "D-GPD"
}

# Zipf-Mandelbrot
name = "zm"
if(any(name==models)){
  # reparametrization: (s,q)=(1+1/xi,si/xi)
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.zm.un,s=1+1/par[1],q=par[2]/par[1]))*t) + length(s)*log(pm.zm.norm(s=1+1/par[1],q=par[2]/par[1]))
  
  o = NULL; o = try(optim(par=c(1,1),fn=nll,hessian=T),silent=T)
  if(is.character(o)){o = optim(par=c(1,1),fn=nll,hessian=F)}
  val[[name]] = o$value
  if(is.null(o$hessian)==F) par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
  stepFun[[name]] = from.pdist.to.stepfun(pdist=function(x) p.czm(x+1,s=1+1/o$par[1],q=o$par[2]/o$par[1]),a=0,b=max(s)+1)
  ks[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided")$p.value # Kolgomorov Smirnov test
  if(sim.pvalue) ks.sim[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=T,B=B.sim)$p.value # Kolgomorov Smirnov test
  
  
  if(F){ # QQ-plots
    QQplot(s, function(x) q.czm(x, s=1+1/o$par[1],q=o$par[2]),a=0,b=0.9,prop=0.05,B=B,lty=2)
    QQplot(s, function(x) q.czm(x, s=1+1/o$par[1],q=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.99,B=B,lty=2)
  }
  
  if(F){ # Check if right-truncated Zipf-Mandelbrot distributon with a third parameter m is a better model
    nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.tzm.un,s=1+1/par[1],q=par[2]/par[1],m=max(s)+par[3]))*t) + length(s)*log(pm.tzm.norm(s=1+1/par[1],q=par[2]/par[1],m=max(s)+par[3]))
    m0 = max(s)*10; m0=10
    o.trunc = optim(par=c(1,1,m0),fn=nll,hessian=T)
    o.trunc$par
    o.trunc$value<= o$value # lower likelihood?
    bic(o.trunc$value,3,nbr.obs=length(s))<=bic(o$value,2,nbr.obs=length(s)) # improvement?
  }
  latexName[[name]] = "ZM"
}

# Negative binomial

name = "nb"
if(any(name==models)){
  
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),dnbinom,size=par[1],prob=par[2]))*t) 
  o = try(optim(par=c(1,0.01),fn=nll,hessian=T) ,silent=T)
  if(is.character(o)) o = optim(par=c(1,0.01),fn=nll,hessian=F)
  
  val[[name]] = o$value
  if(is.null(o$hessian)==F) par[[name]] = get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
  stepFun[[name]] = from.pdist.to.stepfun(pdist=function(x) pnbinom(x,size=o$par[1],prob=o$par[2]),a=0,b=max(s)+1)
  ks[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=F)$p.value # Kolgomorov Smirnov test
  if(sim.pvalue) ks.sim[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=T,B=B.sim)$p.value # Kolgomorov Smirnov test
  
  
  
  #pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
  #par(pty="s")
  #QQplot(s,function(x) qnbinom(x,size=o$par[1],prob=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
  #dev.off()
  
  #pdf(file=paste0("images/qqplot99_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
  #par(pty="s")
  #QQplot(s,function(x) qnbinom(x,size=o$par[1],prob=o$par[2]),a=0,b=0.9,prop=0.05,ci=T,pmax=0.99,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
  #dev.off()
  
  #  if(dataset=="X15"){
  #    pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
  #    par(pty="s")
  #    QQplot(s, function(x) qnbinom(x,size=o$par[1],prob=o$par[2]),a=0,b=0.3,prop=0.1,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
  #    abline(1,1)
  #    dev.off()
  #  }
  
  latexName[[name]] = "NB"
}



# GPD 
name = "gpd"
if(any(name==models)){
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),d.gpd,xi=par[1],si=par[2]))*t) 
  o = try(optim(par=c(1,1),fn=nll,hessian=T),silent=T)
  if(is.character(o)==F){ 
    o$value
    par[[name]] =  get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
    
    # fitted distribution is for continuous rv; following test assumes that rounding the rv gives the data
    # equivalent to test that the data -1/2 are floored function of the fitted distribution 
    stepFun[[name]] = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s))
    ks[[name]] = dgof::ks.test(x=s-1/2,y=stepFun[[name]],alternative="two.sided")$p.value # Kolgomorov Smirnov test
    if(sim.pvalue) ks.sim[[name]] = dgof::ks.test(x=s-1/2,y=stepFun[[name]],alternative="two.sided",simulate.p.value=T,B=B.sim)$p.value # Kolgomorov Smirnov test
    # att: double check, seems not done correctly. 
    
    if(F){
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,B=B,lty=2)
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.99,B=B,lty=2)
      QQplot(s, function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),a=0,b=0.9,prop=0.05,pmax=0.9,B=B,lty=2)
    }
  }
  latexName[[name]] = "$\\text{GPD}_{\\delta=0}$"
  
}

# GPD -1/2
name = "gpd2"
if(any(name==models)){
  nll = function(par) -sum(log(sapply(as.numeric(names(t)),d.gpd,xi=par[1],si=par[2],mu=-1/2))*t) 
  o = try(optim(par=c(1,1),fn=nll,hessian=T),silent=T)
  if(is.character(o)==F){
    o$value
    par[[name]] =  get.par.conf.int(par=o$par,hessian=o$hessian,alpha=0.9)
    
    # fitted distribution is for continuous rv; following test assumes that rounding the rv gives the data
    # equivalent to test that the data are floored function of the fitted distribution but with mu=0 instead of mu=-1/2
    stepFun[[name]] = from.pdist.to.stepfun(pdist=function(x) p.dgpd(x,xi=o$par[1],si=o$par[2]),a=0,b=max(s)+1)
    ks[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided")$p.value # Kolgomorov Smirnov test
    if(sim.pvalue) ks.sim[[name]] = dgof::ks.test(x=s,y=stepFun[[name]],alternative="two.sided",simulate.p.value=T,B=B.sim)$p.value # Kolgomorov Smirnov test
    
    if(dataset=="X15"){
      pdf(file=paste0("images/qqplot_",lang,"_",dataset,"_",name,".pdf"),pointsize=pointsize,width=7,height=7)
      par(pty="s")
      QQplot(s-1/2, function(x) q.gpd(x,xi=o$par[1],si=o$par[2],mu=-1/2),a=0,b=0.3,prop=0.1,B=B,lty=2,xlab="Empirical Quantiles",ylab="Fitted Quantiles")
      abline(1,1)
      dev.off()
    }
  }
  latexName[[name]] = "$ \\text{GPD}_{\\delta=-\\frac{1}{2}} $"
}


# =======
# EXPORT RESUlTS IN A LATEX TABLE


# print fixed number of decimals
prt.rd = function(x,rd) sprintf(paste0("%.",as.character(rd),"f"),x)
# build confidence intervals in latex
ci.in.latex= function(par,par.d,par.u,rd) as.character(paste0("$", prt.rd(par,rd), "_{[",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"]}","$"))


txt = "Model & KS &  nll & $ \\xi$ & $ \\sigma $    \\\\ \n"
if(experiment==0) txt = paste(txt,"British", "& & & & \\\\ \n")
if(experiment==9) txt = paste(txt,"French (in subtitles)", "& & & & \\\\ \n")
if(experiment==10) txt = paste(txt,"French (in books)", "& & & & \\\\ \n")
if(experiment==15) txt = paste(txt,"French", "& & & & \\\\ \n")



for(name in models){
  if(all(!is.na(par[[name]]))){
    txt = paste(txt, latexName[[name]], "& $",
                prt.rd(ks[[name]],rd=2), "$ & $", 
                prt.rd(val[[name]],rd=1), "$ &", 
                ci.in.latex(par[[name]][1,1],par[[name]][1,2],par[[name]][1,3],rd=2), "&", 
                ci.in.latex(par[[name]][2,1],par[[name]][2,2],par[[name]][2,3],rd=2), 
                "\\\\","\n")
  } else{
    txt = paste(txt,latexName[[name]], "& NA & & & &\\\\ \n")
  }
}  

set.seed(191)
if(ks.bootstrap){ # Compute KS test for bootstrap
  
  for(j in 1:B.bootstrap){
    print_loading(it=j,it.max=B.bootstrap)
    s.sub = sample(s,replace=T)
    
    for(name in models){
      if(is.null(stepFun[[name]])==F){
        ks.boot[[name]][j] = dgof::ks.test(x=s.sub,y=stepFun[[name]],alternative="two.sided",simulate.p.value=simulated.p.value, B=B.sim)$p.value # Kolgomorov Smirnov test
      }
    }
  }
}
ks.boot.m = list(NULL)
for(name in models){
  ks.boot.m[[name]]= mean(ks.boot[[name]])
}
unlist(ks.boot.m)

# Table 2 (first part)
cat("Table 2 (first part) in Discrete Extremes: \n")
if(T){
  cat(txt)
  cat(paste0("Threshold = ",u,", quantile:",qt,", size of exceed. = ",length(s)),", size init data = ",length(data),"\n")
  
  round(unlist(val),2)
  round(unlist(ks),2)
  round(unlist(ks.sim),2)
  print(round(unlist(ks.boot.m),2))
  
  # compute BIC
  round(sapply(unlist(val),bic,nbr.par=2,nbr.obs=sum(t)),1)
  
}

