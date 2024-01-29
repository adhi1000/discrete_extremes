# ================
# SIMULATE DATA FROM DISCRETE DISTRIBUTIONS IN MDA for xi=0, xi>0, xi<0
# FIT DIFFERENT MODELS TO EXCEEDANCES
# COMPARE PARAMETER AND QUANTiLE ESTIMATION


myRFolder = "~/Downloads/discrete_extremes"
setwd(myRFolder)
source("Functions.R")
source("utilities_5Feb2017.R")

# Graphical settings
pointsize = 18

# ===========

# POSSIBLE MODELS

if(F) name.list = c("gpd.disc","hz","zmN","hz.rp","zm.rp","ismev.u0" ,"gpd.u0","ismev.u2" ,"gpd.u2","gpd.cont")


# ==============
# SETTINGS

# Number of total observations 
sampleSize = 8000 # 8000

# Quantile level at which the threshold is taken
q.u = 0.95

# Quantile level to estimate
p.target = 0.9999

# Nbr time experiment repeated
B = 500 

# Models to fit
name.list = c("gpd.disc","hz","zmN","hz.rp","zm.rp","ismev.u0" ,"gpd.u0","ismev.u2" ,"gpd.u2")
name.list = c("ismev.on.cont", "gpd.disc","ismev.u2" ,"gpd.u2","ismev.u0","gpd.u0","emp")
name.list = c("headerContModels","ismev.on.cont", "headerDiscrModels", "gpd.disc","gzm.rp","ismev.u2" ,"ismev.u0","gpd.u2","nb")



# Experiment
mda.list = c("fre","gum","wei")

mda.list = c("gum")
mda.list = c("wei")
mda.list = c("fre")

# Other settings

par(pty="s")
rd = 2
do.pdf = T
alpha = 0.1 # confidence interval at (1-alpha)*100 %


# =================
# SIMULATION FROM Y and its discretization X


for(mda in mda.list){
    
    
    
    if(mda=="fre"){
        # =================
        # FRECHET DOMAIN OF ATTRACTION
        
        
        # distribution of Y
        experiment = 2
        
        if(experiment ==1){ 
            df = 4
            
            rdist = function(n) rt(n,df=df)
            pdist = function(x) pt(x,df=df)  
            qdist = function(x) qt(x,df=df)  
            
            # true xi
            true.xi= 1/df 
        } else if(experiment == 2){
        
        a= 2 # in biometrika: 2, in thesis 2
        b = 0.75  # in biometrika: 0.75, in thesis: 1
        
        rdist = function(n) actuar::rinvgamma(n, shape=a, scale = b)
        pdist = function(x) actuar::pinvgamma(x, shape=a, scale = b)
        qdist = function(x) actuar::qinvgamma(x, shape=a, scale = b) 
        
        # true xi
        true.xi= 1/a
        }
        
        # round function
        lambda = 1
        h = 1
        
        
    } else if(mda=="gum"){
        # =================
        # GUMBEL DOMAIN OF ATTRACTION
        
        
        # distribution of Y
        sd = 3
        rdist = function(n) rnorm(n,mean=0,sd=sd)
        pdist = function(x) pnorm(x,mean=0,sd=sd)
        qdist = function(x)  qnorm(x,mean=0,sd=sd)  
        
        
        # true xi: 
        true.xi = 0
        
        # round function
        lambda = 1
        h = 1/2
        
        
    } else if(mda=="wei"){
        
        # =================
        # WEIBULL DOMAIN OF ATTRACTION
        
        
        # distribution of Y
        if(F){ # previous
            c=17;  d=0; shape1 = 1; shape2 = 10
        }
        
        c=14.7  # 7.1: dgpd fails, 7.9: zm better, 7.5: zm slightly better (NA's ci)
        d=0
        shape1 = 1
        shape2 = 2
        
        
        rdist = function(n) rbeta(n,shape1=shape1,shape2=shape2)*c+d
        pdist = function(x) pbeta((x-d)/c,shape1=shape1,shape2=shape2)
        qdist = function(x) c*qbeta(x,shape1=shape1,shape2=shape2)+d
        stopifnot(round(pdist(qdist(0.3))-0.3,6)==0)   
        
        # endpoint: e1=si/|xi|
        # true xi: (-1/shape1?) not -1/shape2?
        true.xi = -1/shape2 # todo check!
        
        # round function
        lambda = 1
        h = 0
        
        
    }
    
    # ============
    # TARGET 
    
    rd.function = function(x) floor(lambda*x+1-h)
    
    # given target
    q.target = rd.function(qdist(p.target))
    q.target.cont = (q.target-1+h)/lambda
    
    # probability to estimate
    true.p = 1-pdist(q.target.cont) # Prob(X>= q.target)
    
    
    
    # =================
    # RUN EXPERIMENTS: estimation of extreme quantile
    
    # parameter initialization
    val = xi = si = xi.pm = si.pm = xi.pmb = si.pmb = si = ks = tg = pe = pe.pm = u.sim = len.sim = len.in.set = max.obs= latexName = nbr.tied.to.less.than.20= list(NA)
    
    
    for(i in 1:B){
      print_loading(it=i,it.max=B)
      
        
        # Simulate
        
        y = rdist(sampleSize) 
        x = rd.function(y)
        
        
        # threshold
        u = floor(as.numeric(quantile(x,q.u,type=6))) # floored because quantile does not always return an integer
        
        # threshold for continuous data y 
        u.cont = (u-1+h)/lambda
        
        # Select exceedances
        
        
        x.u = x[x>=u]-u
        pu = length(x.u)/length(x)
        s = x.u
        t = table(s)
    
        # IMPROVEMENT: u should be selected to make sure "pu" is close to the target "pu"
        
        
        if(do.pdf){
            pdf(file=paste0("images/","freqTable","_",mda,".pdf"),pointsize=pointsize,width=7,height=7)
            par(pty="s")
            plot(t,ylab="Frequency",xlab="Simulated observations",xaxt = "n")
            axis(side = 1, at = pretty(as.numeric(names(t))))
            dev.off()
        }
        
        if(F){ # Visualize
            length(x.u)
            plot(t)
            plot(x)
            plot(y)
        }
        
        # Save threshold and nbr of exceedances
        u.sim[[i]] = u
        len.sim[[i]] = length(x.u)
        len.in.set[[i]] = length(x.u[x.u>q.target])
        max.obs[[i]] = max(x)
        
      
        nbr.tied.to.less.than.20[[i]] = sum(table(x.u)[table(x.u)<=20])/length(x.u)
        
        
        # GPD discrete
        name = "gpd.disc"
        if(is.element(name,name.list)){
            nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.dgpd,xi=par[1],si=par[2]))*t) 
            o = try( optim(par=c(0.5,1),fn=nll,hessian=T) ,silent=F)
            
            
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                
                val[[name]][i] =  o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
            
                
                # extreme set
                pr =  pu*(1-p.gpd(q.target-u,si=o$par[2],xi=o$par[1])) # pr(X>=q)=1-Pr(X<q)=1-Pr(X<=q-1)=1-F_GPD(q)
                grad.pr = get.grad.sprob.gpd(q.target-u,si=o$par[2],xi=o$par[1],pu=pu) # +1 because F_DPGD(k)=F_GPD(k+1)
                var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                
               # chi test
                if(F){
                    step.fct = from.pdist.to.stepfun(pdist=function(x) p.gpd(x+1,xi=o$par[1],si=o$par[2]),a=0,b=max(s)+10)
                    ks[[name]][i] = dgof::ks.test(x=s,y=step.fct,alternative="two.sided")$p.value # Kolgomorov Smirnov test
                }
                
                # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
                    dev.off()
                }
            }
            latexName[[name]] = "D-GPD"
            }
        
        # Exp discrete
        name = "exp.disc"
        if(is.element(name,name.list)){
            nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.exp,si=par))*t) 
            o = try( optim(par=1,fn=nll,hessian=T) ,silent=F)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                
                val[[name]][i] =  o$value
                si[[name]][i] = o$par
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                si.pm[[name]][i] = v
                
                # chi test
                if(F){
                    step.fct = from.pdist.to.stepfun(pdist=function(x) p.gpd(x+1,xi=o$par[1],si=o$par[2]),a=0,b=max(s)+10)
                    ks[[name]][i] = dgof::ks.test(x=s,y=step.fct,alternative="two.sided")$p.value # Kolgomorov Smirnov test
                }
                
                # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.exp(x,si=o$par))
                    dev.off()
                }
            }
        }
        
        
        # ISMEV gpd threshold 0
        
        name = "ismev.u0"
        o = NULL;
        if(is.element(name,name.list)){
            o = try(ismev::gpd.fit(s,threshold=-0.0000001,show=F) ,silent=T)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                val[[name]][i] =  o$nllh
                xi[[name]][i] = o$mle[2]  # reversed
                si[[name]][i] = o$mle[1]
                v = from.varcov.to.ci(par=o$mle,varcov=o$cov,pm=T)
                xi.pm[[name]][i] = v[2]  # reversed
                si.pm[[name]][i] = v[1]
                # ptarget[[name]][i] = pu*(1-p.gpd(target-u,a=1/o$mle[2],si=o$mle[1]))
                
                # extreme set
                pr =  pu*(1-p.gpd(q.target-u,si=o$mle[1],xi=o$mle[2]))
                grad.pr = get.grad.sprob.gpd(q.target-u,si=o$mle[1],xi=o$mle[2],pu=pu)
                var.pr = rev(grad.pr) %*% o$cov %*% rev(grad.pr)
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                
                # chi test
                #step.fct = from.pdist.to.stepfun(pdist=function(x) p.gpd(x,a=1/o$mle[2],si=o$mle[1]),a=0,b=max(s)+10)
                # todo: continuous version
                #ks[[name]][i] = dgof::ks.test(x=s,y=step.fct,alternative="two.sided")$p.value # Kolgomorov Smirnov test
                
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.gpd(x,xi=o$mle[2],si=o$mle[1]))
                    dev.off()
                }
            }
            latexName[[name]] = "$\\text{GPD}_{\\delta=0}$"
        }
        
        # ISMEV gpd threshold -1/2
        
        name = "ismev.u2"
        #s.m = s+1/2
        
      
        if(is.element(name,name.list)){
            o = try(ismev::gpd.fit(s,threshold=-0.5,show=F) ,silent=T)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                
                val[[name]][i] =   o$nllh
                xi[[name]][i] = o$mle[2] # reversed
                si[[name]][i] = o$mle[1]
                v = from.varcov.to.ci(par=o$mle,varcov=o$cov,pm=T)
                xi.pm[[name]][i] = v[2]  # reversed
                si.pm[[name]][i] = v[1]
                
                
                # extreme set
                pr =  pu*(1-p.gpd(q.target-u+1/2,si=o$mle[1],xi=o$mle[2])) # +1/2 because gpd model fitted to s+1/2
                grad.pr = get.grad.sprob.gpd(q.target-u+1/2,si=o$mle[1],xi=o$mle[2],pu=pu)
                var.pr = rev(grad.pr) %*% o$cov %*% rev(grad.pr)
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                
                if(F){
                # extreme set 2nd version!!! without shift 1/2
                pr =  pu*(1-p.gpd(q.target-u,si=o$mle[1],xi=o$mle[2])) 
                grad.pr = get.grad.sprob.gpd(q.target-u,si=o$mle[1],xi=o$mle[2],pu=pu)
                var.pr = rev(grad.pr) %*% o$cov %*% rev(grad.pr)
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                }
              
                
                
                # chi test
                #step.fct = from.pdist.to.stepfun(pdist=function(x) p.gpd(x-u0,a=1/o$mle[2],si=o$mle[1]),a=0,b=max(s)+10)
                # todo: continuous version
                #ks[[name]][i] = dgof::ks.test(x=s,y=step.fct,alternative="two.sided")$p.value # Kolgomorov Smirnov test
                
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.gpd(x,xi=o$mle[2],si=o$mle[1]))
                    dev.off()
                }
                
            }
            latexName[[name]] = "$ \\text{GPD}_{\\delta=-\\frac{1}{2}} $"
        }
        
        # GPD  continuous threshold 0
        
        name = "gpd.u0"
        if(is.element(name,name.list)){
            if(i==0)  name.list = c(name.list,name)
            nll = function(par) -sum(sapply(as.numeric(names(t)),log.d.gpd,xi=par[1],si=par[2])*t) 
            o = NULL
            o = try(optim(par=c(0.5,1),fn=nll,hessian=T)  ,silent=T)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                val[[name]][i] =  o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                # return level
                rl = get.return.level.gpd(alpha=p.target,u=u,pu=pu,si=o$par[2],xi=o$par[1])
                grad.rl = get.grad.return.level.gpd(alpha=p.target,u=u,pu=pu,si=o$par[2],xi=o$par[1])
                var.rl = rev(grad.rl) %*% solve(o$hessian) %*% rev(grad.rl)
                pe[[name]][i] = rl
                pe.pm[[name]][i] =from.var.to.ci(par=rl, var=var.rl, alpha=alpha,pm=T)
                
                
                # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
                    dev.off()
                }
            }
        }
        # GPD continuous threshold -1/2
      
        name = "gpd.u2"
        if(is.element(name,name.list)){
            s.m = s+1/2; t.m = table(s.m)
            nll = function(par) -sum(sapply(as.numeric(names(t.m)),log.d.gpd,xi=par[1],si=par[2])*t.m) 
            o = try( optim(par=c(0.5,1),fn=nll,hessian=T) ,silent=T)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
              
           
                val[[name]][i] =  o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                #tg[[name]][i] =  get.return.level(qt=q.target,u=u,pu=pu,q.exceed.dist= function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
                
                # extreme set
                pr =  pu*(1-p.gpd(q.target-u+1/2,si=o$par[2],xi=o$par[1])) 
                grad.pr = get.grad.sprob.gpd(q.target-u+1/2,si=o$par[2],xi=o$par[1],pu=pu)
                var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                
                
              
                
                
                # Export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
                    dev.off()
                }
            }
            latexName[[name]] = name
        }
        
        
        # Reparametrized Zipf Mandelbrot distribution (N=Inf)
        
        name = "zm.rp"
        if(is.element(name,name.list)){
            nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.zm.un,s=1+1/par[1],q=par[2]/par[1]))*t) + length(s)*log(pm.zm.norm(s=1+1/par[1],q=par[2]/par[1]))
            
            o = try( optim(par=c(0.5,1),fn=nll,hessian=T) ,silent=T)
            #if(is.character(o)){ 
            #    o = try( optim(par=c(1,1),fn=nll,hessian=F) ,silent=T)    
            #    val[[name]][i] = o$value
            #  
            
            if(is.character(o)){     
                val[[name]][i] = NA
            } else{
                val[[name]][i] = o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.czm(x,s=1+1/o$par[1],q=o$par[2]/o$par[1],upper=max(s)*100) )
                    dev.off()
                }
                
                # return level
                rl = get.return.level(alpha=p.target,u=u,pu=pu,q.exceed.dist= function(x) q.czm(x,s=1+1/o$par[1],q=o$par[2]/o$par[1],upper=max(s)*100))
                f =  function(par) get.return.level(alpha=p.target,u=u,pu=pu,q.exceed.dist= function(x) q.czm(x,s=1+1/par[1],q=par[2]/par[1],upper=max(s)*100))
                grad.rl = try(numDeriv::grad(func=f,x=o$par),silent=T)
                if(is.character(grad.rl)){
                    rl.pm = NA
                } else {
                    var.rl = grad.rl %*% solve(o$hessian) %*% grad.rl
                    rl.pm = from.var.to.ci(par=rl, var=var.rl, alpha=alpha,pm=T)
                }
                pe[[name]][i] = rl
                pe.pm[[name]][i] = rl.pm
                
                
             
            }
            
        } 
        
        
        
        # Generalized Zipf Mandelbrot distribution (N=Inf) 
        
        name = "gzm.rp"
        if(is.element(name,name.list)){
            nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.gzm.un,xi=par[1],si=par[2]))*t) + length(s)*log(pm.gzm.norm(xi=par[1],si=par[2]))
            o = try( optim(par=c(0.5,1),fn=nll,hessian=T) ,silent=T)
        
            if(is.character(o)){     
                val[[name]][i] = NA
            } else{
                val[[name]][i] = o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.cgzm(x,xi=o$par[1],si=o$par[2],upper=max(s)*100) )
                    dev.off()
                }
                
                if(F) plot( sapply(runif(1000),q.cgzm,xi=o$par[1],si=o$par[2],upper=max(s)*100) )
                
                
                # extreme set
                pr =  pu*(1-p.gzm(q.target-u-1,si=o$par[2],xi=o$par[1])) # pr(X>=q) = 1-pr(X<q)=1-pr(X<=q-1)
                f =  function(par) pu*(1-p.gzm(q.target-u-1,si=par[2],xi=par[1]))
                grad.pr = numDeriv::grad(func=f,x=o$par)
                sigma.inv = solve(o$hessian)
                var.pr = grad.pr %*% sigma.inv %*% grad.pr
                
                
                if(F){
                if(var.pr<0 && matrixcalc::is.positive.definite(sigma.inv)==F){
                    print("In gzm.rp: Hessian inverse replaced by closest semi-definite positive matrix") 
                    sigma.inv.PD = as.matrix(Matrix::nearPD(sigma.inv)$mat)
                    var.pr = grad.pr %*% sigma.inv.PD %*% grad.pr
                    
                    sigma.PD = as.matrix(Matrix::nearPD(o$hessian)$mat)
                    sigma.inv.PD = solve(sigma.PD)
                    
                    sigma.inv.PD %*% o$hessian
                    sigma.inv %*% o$hessian
                    
                    
                }
                }
                
                
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
             
                
            if(F){    grad.pr = try(numDeriv::grad(func=f,x=o$par),silent=T)
                if(is.character(grad.pr)){
                    pr.pm = NA
                } else {
                    var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
                    pr.pm = from.var.to.ci(par=rl, var=var.pr, alpha=alpha,pm=T)
                }
                pe[[name]][i] = pr
                pe.pm[[name]][i] = pr.pm
            }
                
                
            }   
            }
        latexName[[name]] = "ZM"
        
        
        # Zipf Mandelbrot distribution (N=Inf)
        # (usually same as zm.rp)
        name = "zm"
        if(is.element(name,name.list)){
            nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.zm.un,s=par[1],q=par[2]))*t) + length(s)*log(pm.zm.norm(s=par[1],q=par[2]))
            o = NULL; #o = optim(par=c(1,1),fn=nll,method="Nelder-Mead",hessian=T)  
            o = try( optim(par=c(2,1),fn=nll,hessian=T) ,silent=T)
            if(is.character(o)){ 
                warning(o[1])
                o = try( optim(par=c(2,1),fn=nll,hessian=F) ,silent=T)    
                val[[name]][i] = o$value
                if(is.character(o)){     
                    val[[name]][i] = NA
                }
            } else{
                val[[name]][i] = o$value
                
                
                # todo: parameter transformation
                if(F){
                    xi[[name]][i] = o$par[1]
                    si[[name]][i] = o$par[2]
                    v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                    xi.pm[[name]][i] = v[1]
                    si.pm[[name]][i] = v[2]
                }
                
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.czm(x,s=o$par[1],q=o$par[2]) )
                    dev.off()
                }
            }
            
            pe[[name]][i]  = get.return.level(alpha=p.target,u=u,pu=pu,q.exceed.dist= function(x) q.czm(x,s=o$par[1],q=o$par[2]))
        }
        
        
        # Hurwitz zeta
        
        name = "hz"
        if(is.element(name,name.list)){
            nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.hzeta.un,s=1+1/par[1],q=par[2]))*t) + length(s)*log(pm.hzeta.norm(s=1+1/par[1],q=par[2]))
            o = try( optim(par=c(1,1),fn=nll,hessian=T) ,silent=T)
            if(is.character(o)){ 
                warning(o[1])
                val[[name]][i] = NA
            } else{
                val[[name]][i] = o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                
                
                # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s,function(x) q.from.cdf(p=x,cdf=function(x) p.hzeta(x,s=1/o$par[1]+1,q=o$par[2]),low=0,up=999999))
                    dev.off()
                }
            }
        }
        
        
        
        name = "zmN"
        if(is.element(name,name.list)){
            s1 = s+1
            t1 = table(s1)
            N=max(s1)*10
            nll = function(par) -sum(sapply(as.numeric(names(t1)),tolerance::dzipfman,log=F,s=1/par[1],b=par[2],N=N)*t1)
            o = try(  optim(par=c(-0.01,1),fn=nll,hessian=T),silent=T)
            if(is.character(o)){ 
                warning(o[1])
                val[[name]][i] = NA
            } else{
                val[[name]][i] = o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s1,function(x) tolerance::qzipfman(x, s=1/o$par[1], b = o$par[2], N = N))
                    dev.off()
                }
            }
        }
        
        
        name = "zmN.rp"
        if(is.element(name,name.list)){
            s1 = s+1
            t1 = table(s1)
            N=max(s1)*10
            nll = function(par) -sum(sapply(as.numeric(names(t1)),tolerance::dzipfman,log=T,s=1+1/par[1],b=par[2]/par[1],N=N)*t1)
            o =  try(  optim(par=c(0.5,1),fn=nll,hessian=T),silent=T)
            if(is.character(o)){ 
                warning(o[1])
                val[[name]][i] = NA
            } else{
                val[[name]][i] = o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(s1,function(x) tolerance::qzipfman(x, s=1/o$par[1], b = o$par[2]/o$par[1], N = N))
                    dev.off()
                }
            }
        }
        
        name = "emp"
        if(is.element(name,name.list)){
            pe[[name]][i] = quantile(x,p.target,type=6)
            val[[name]][i] = NA
        }
        
        
        # GPD on continuous data
        name = "ismev.on.cont"
        o = NULL;
        y.u = y[y>=u.cont]-u.cont
        pu.cont = length(y.u)/length(y)
        stopifnot(round(pu.cont-pu,6)==0)
        
        if(is.element(name,name.list)){
            o = try(ismev::gpd.fit(y.u,threshold=0,show=F) ,silent=T)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                val[[name]][i] =  o$nllh
                xi[[name]][i] = o$mle[2]  # reversed
                si[[name]][i] = o$mle[1]
                v = from.varcov.to.ci(par=o$mle,varcov=o$cov,pm=T)
                xi.pm[[name]][i] = v[2]  # reversed
                si.pm[[name]][i] = v[1]
                #ptarget[[name]][i] = pu*(1-p.gpd(target-u.cont,a=1/o$mle[2],si=o$mle[1]))
                
                # extreme set
                pr =  pu.cont*(1-p.gpd(q.target.cont-u.cont,si=o$mle[1],xi=o$mle[2])) 
                grad.pr = get.grad.sprob.gpd(q.target.cont-u.cont,si=o$mle[1],xi=o$mle[2],pu=pu.cont)
                var.pr = rev(grad.pr) %*% o$cov %*% rev(grad.pr)
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                
                
                # chi test
                #step.fct = from.pdist.to.stepfun(pdist=function(x) p.gpd(x,a=1/o$mle[2],si=o$mle[1]),a=0,b=max(s)+10)
                # todo: continuous version
                #ks[[name]][i] = dgof::ks.test(x=s,y=step.fct,alternative="two.sided")$p.value # Kolgomorov Smirnov test
                
                if(i==1 && do.pdf){
                    QQplot(y.u,function(x) q.gpd(x,xi=o$mle[2],si=o$mle[1]))
                }
            }
            latexName[[name]] = "GPD"
        }
        
        
        # GPD on continuous data
        name = "gpd.cont"
        o = NULL;
        y.u = y[y>=u.cont]-u.cont
        pu.cont = length(y.u)/length(y)
        stopifnot(round(pu.cont-pu,6)==0)
        t.cont = table(y.u)
        
        if(is.element(name,name.list)){
            nll = function(par) -sum(sapply(as.numeric(names(t.cont)),log.d.gpd,xi=par[1],si=par[2])*t.cont) 
            o = NULL
            o = try(optim(par=c(1,1),fn=nll,hessian=T)  ,silent=T)
            if(is.character(o)){ 
                warning(o[1]); pe[[name]][i] = NA
            } else{
                val[[name]][i] =  o$value
                xi[[name]][i] = o$par[1]
                si[[name]][i] = o$par[2]
                v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
                xi.pm[[name]][i] = v[1]
                si.pm[[name]][i] = v[2]
                
                # extreme set
                pr =  pu.cont*(1-p.gpd(q.target.cont-u.cont,si=o$par[2],xi=o$par[1])) 
                grad.pr = get.grad.sprob.gpd(q.target.cont-u.cont,si=o$par[2],xi=o$par[1],pu=pu.cont)
                var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
                pe[[name]][i] = pr
                pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
                
              # export qq plot
                if(i==1 && do.pdf){
                    pdf(file=paste0("images/",mda,"_",name,".pdf"),pointsize=13,width=7,height=7)
                    par(pty="s")
                    QQplot(y.u,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]))
                    dev.off()
                }
            }
            latexName[[name]] = "GPD not ismev"
        }
    
        
        
        # Negative binomial 
        name = "nb"
        if(is.element(name,name.list)){
          nll = function(par) -sum(log(sapply(as.numeric(names(t)),dnbinom,size=par[1],prob=par[2]))*t) 
          
          o = try( optim(par=c(1,0.1),fn=nll,hessian=T) ,silent=F)
          
          
          if(is.character(o)){ 
            warning(o[1]); pe[[name]][i] = NA
          } else{
            
            val[[name]][i] =  o$value
            xi[[name]][i] = o$par[1]
            si[[name]][i] = o$par[2]
            v = from.hessian.to.ci(par=o$par,hessian=o$hessian,pm=T)
            xi.pm[[name]][i] = v[1]
            si.pm[[name]][i] = v[2]
            
            
            # extreme set
            pr =  pu*(1-pnbinom(q.target-u,prob=o$par[2],size=o$par[1])) # pr(X>=q)=1-Pr(X<q)=1-Pr(X<=q-1)=1-F_GPD(q)
            #grad.pr = get.grad.sprob.gpd(q.target-u,si=o$par[2],size=o$par[1],pu=pu) # +1 because F_DPGD(k)=F_GPD(k+1)
            #var.pr = grad.pr %*% solve(o$hessian) %*% grad.pr
            pe[[name]][i] = pr
            #pe.pm[[name]][i] = from.var.to.ci(par=pr, var=var.pr, alpha=alpha,pm=T)
            
          }
          latexName[[name]] = "NB"
        }
        
        
    }    
    
    
    
    
    
    
    
    name = "headerContModels"
    if(is.element(name,name.list)){
    latexName[[name]] = "\\hline $Y-u' \\mid Y\\geq u' $ "
    }
    name = "headerDiscrModels"
    if(is.element(name,name.list)){
        latexName[[name]] = "\\hline  $X-u\\mid X\\geq u$"
    }
    
    prt.rd = function(x,rd) sprintf(paste0("%.",as.character(rd),"f"),x)
    ci.in.latex.1 = function(par,par.d,par.u,rd) as.character(paste0("$", prt.rd(par,rd), "_{[",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"]}","$"))
    ci.in.latex.2 = function(par.d,par.u,rd) as.character(paste0("$ [",prt.rd(par.d,rd),",",prt.rd(par.u,rd),"] $"))
    rd = 2
    
    # Might be useful
    prt.rd = function(x,rd) sprintf(paste0("%.",as.character(rd),"f"),x)
    
    # constants
    pw1 = 4
    k1 = 10^pw1
    
    
    # TODO: add likelihood for ZM and DGPD
    
    # =============
    # Create latex tabular 
    
    # Header
    text = ""
    text = paste0(text," & $ \\pr(X\\geq ", q.target," ) \\times 10^",pw1,", $&  &  $ \\xi $  & $ \\sigma$   \\\\  \n",
                  " \\hline\\hline \n",
                 "&  mean (coverage, length)   & true length    & mean (coverage, length) &  mean (length) \\\\ \n",
        "truth & $ \\textbf{", prt.rd(true.p*k1,rd),"} $ & &  $",prt.rd(true.xi,rd),"$ &  \\\\ \n")
  
    xi.m = si.m = xi.u = si.u = xi.d = si.d = xi.pm.m = si.pm.m = ks.m= pe.m = pe.u = pe.d = pe.true.len = NA
    val.m = conv= conv.ci =nbr.not.conv.ci=nbr.not.conv= list(NA)

    
    for(name in name.list){
        # number of NA's 
        w = is.na(pe[[name]])==F
        w.ci = is.na(pe.pm[[name]])==F
        
        if(all(w==F)){ 
            text =  paste(text,latexName[[name]],"&  & $ $ & \\\\  \n")
        } else{
        
        # percentage of successful convergence of the optimization
        conv[[name]] =  (sum(w)/length(w))*100
        conv.ci[[name]] =  (sum(w.ci)/length(w.ci))*100
        
        # nbr of succ
        nbr.not.conv[[name]] =  length(w)-sum(w)
        nbr.not.conv.ci[[name]] =   length(w)-sum(w.ci)
        
        
        # mean negative log likelihood
        val.m[[name]] = mean(as.vector(val[[name]][w]))
        
        # MLE
        
        # mean
        xi.m = mean(as.vector(xi[[name]][w]))
        si.m = mean(as.vector(si[[name]][w]))
        
        # mean length
        xi.length = mean(as.vector(2*xi.pm[[name]][w.ci]))
        si.length = mean(as.vector(2*si.pm[[name]][w.ci]))
        
        # coverage of xi
        xi.inside = (xi[[name]][w.ci]-xi.pm[[name]][w.ci] <= true.xi) & (true.xi <= xi[[name]][w.ci]+xi.pm[[name]][w.ci])
        xi.coverage = sum(xi.inside)/length(w)*100
        
        # Extreme set
        pe.m =  mean(as.vector(pe[[name]][w]))
        pe.length = mean(as.vector(na.omit(2*pe.pm[[name]][w.ci])))
        pe.inside = (pe[[name]][w.ci]-pe.pm[[name]][w.ci] <= true.p) & (true.p <= pe[[name]][w.ci]+pe.pm[[name]][w.ci])
        pe.coverage = sum(pe.inside)/length(w)*100
        pe.true.len = quantile(as.vector(pe[[name]][w]),1-alpha/2)-quantile(as.vector(pe[[name]][w]),alpha/2)
        
        
       text =  paste(text,latexName[[name]],
                  "& $", prt.rd(pe.m*k1,rd),"\\, (",prt.rd(pe.coverage,0),"\\%, \\, ",prt.rd(pe.length*k1,rd),")", 
                  "$ & $", prt.rd(pe.true.len*k1,rd),
                  "$ & $", prt.rd(xi.m,rd),"\\, (",prt.rd(xi.coverage,0),"\\%, \\, ",prt.rd(xi.length,rd),")",
                  "$ & $",prt.rd(si.m,rd),"\\, (",prt.rd(si.length,rd),")", 
                  "$ \\\\ \n")
    }
    }
    cat(text)
    
    print(unlist(val.m)) # log-likelihood
    
    # Converging methods:
    print(unlist(conv))# (check 100%)
    print(unlist(conv.ci))# (check 100%)
    
    print(unlist(nbr.not.conv))# (check 0%)
    print(unlist(nbr.not.conv.ci))# (check 0%)
    
    # Check proportion exceedances is close to target q.u (gap due to rounding threshold)
    print(q.u); print(1-mean(unlist(len.sim))/length(y))
  
    # proportion of tied observations
    round(mean(unlist(nbr.tied.to.less.than.20)),2)
    
    cat(paste("Info","\n","MDA=",mda,",",lambda="lambda, ", "sample size =",sampleSize,", q.u=",q.u,", B=",B,", qt target =",round(q.target,rd), 
               ", prob above target = ", p.target,"\n","mean size of exceed = ",round(mean(unlist(len.sim)),1),", mean thresh = ",round(mean(unlist(u.sim)),1),
               "\n", "true xi = ", round(true.xi,3), "mean maximal obs = ",round(mean(unlist(len.in.set))),1),"\n")
    
}


if(F){
# ========================
# POISSON DISTRIBUTION
set.seed(90)
do.pdf =T
n = 5000
rate = 1
x = rpois(n,rate)
u = qpois(0.95,rate)
s = x[x>=u]-u
t = table(s)
plot(s)
length(s);u
plot(t)


if(do.pdf){
    pdf(file=paste0("images/","Poiss",".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    #plot(t,xlab="",ylab="Frequency")
    plot(table(x),xlab="",ylab="Frequency")
    dev.off()
}
# GPD discrete
name = "gpd_disc"
nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.dgpd,xi=par[1],si=par[2]))*t) 
o = try( optim(par=c(0.5,1),fn=nll,hessian=T) ,silent=F)
si.dgpd=o$par[2]
# export qq plot
if(do.pdf){
    pdf(file=paste0("images/","poisson","_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$par[1],si=o$par[2]),xlab="Empirical Quantiles",ylab="Fitted Quantiles")
    abline(1,1)
    dev.off()
}

# GPD continuous
u0 = 0
name = "ismev"
o = try(ismev::gpd.fit(s,threshold=u0,show=F) ,silent=T)
si.gpd=o$mle[1]
if(do.pdf){
    pdf(file=paste0("images/","poisson","_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$mle[2],si=o$mle[1]),xlab="Empirical Quantiles",ylab="Fitted Quantiles")
    abline(1,1)
    dev.off()
}

# GPD continuous
u0 = -1/2
name = "ismev_u2"
o = try(ismev::gpd.fit(s,threshold=u0,show=F) ,silent=T)
if(do.pdf){
    pdf(file=paste0("images/","poisson","_",name,".pdf"),pointsize=pointsize,width=7,height=7)
    par(pty="s")
    QQplot(s,function(x) q.gpd(x,xi=o$mle[2],si=o$mle[1]),xlab="Empirical Quantiles",ylab="Fitted Quantiles")  
    abline(1,1)
    dev.off()
}


# Generalized Zipf Mandelbrot distribution 
name = "gzm.rp"
    nll = function(par) -sum(log(sapply(as.numeric(names(t)),pm.gzm.un,xi=par[1],si=par[2]))*t) + length(s)*log(pm.gzm.norm(xi=par[1],si=par[2]))
    o = try( optim(par=c(0.5,1),fn=nll,hessian=T) ,silent=T)
 si.zm = o$par[2]
        pdf(file=paste0("images/","poisson","_",name,".pdf"),pointsize=pointsize,width=7,height=7)
            par(pty="s")
            QQplot(s,function(x) q.cgzm(x,xi=o$par[1],si=o$par[2],upper=max(s)*100) )
            abline(1,1)
            dev.off()
       
            
# Info in latex:
print(n) # sample size
si.gpd; si.dgpd; si.zm # mle for sigma

# todo: improve q.cgzm 

}

