# simulation from the ME2 model:
# random intercepts in both the membership models and outcome models

library(doParallel)
library(MASS)

source("fecode.R")
source("mecode.R")
source("me2code.R")
source("sim_outcome.R")

#### true parameter values

tbss = c(-0.5,1, 1.5) # beta_ss_1
p=length(tbss)

tbsn=c(-0.3,0.8,1.3) # beta_sn
tbss0=c(-0.2,1,1) # beta_ss_0

tass=c(1.6,0.2,0.1) # alpha_ss
tasn=c(-0.1,-0.1,-0.2) # alpha_sn

nc = 30 # number of clusters each arm
msize = 25 # average cluster size
sdsize = 3 # standard deviation of the cluster size

tt = 2
trho =0.01 # ICC in the outcome models
taut=tt*trho # variance of u_i
tsat=tt*(1-trho) # variance of the error

tgat=0.8 # variance of v_i

nsim=10# number of simulations
nbp=10 # number of bootstrap samples

################# data simulation

rn=matrix(round(rnorm(nc*nsim,msize,sdsize)),nc,nsim) # cluster size
sv1=matrix(rnorm(nc*nsim,0,sqrt(tgat)),nc,nsim) # random intercepts in the membership models of the treatment group
sv0=matrix(rnorm(nc*nsim,0,sqrt(tgat)),nc,nsim) # random intercepts in the membership models of the control group

cl=makeCluster(detectCores())
registerDoParallel(cl)

## kk is the simulation index
## ii is the cluster index

sda1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{
library(MASS)
			ni=rn[ii,kk]
			x1=rbinom(ni,1,1/2)
			x2=rnorm(ni,0,1)
			sx1=cbind(1,x1,x2)	
             simry1(sx1,ni,tbss,tbsn,taut,tsat,tass,tasn,sv1[ii,kk],ii,kk)
}

sda0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc,.combine='rbind') %dopar%{
library(MASS)
             ni=rn[ii,kk]
			x1=rbinom(ni,1,1/2)
			x2=rnorm(ni,0,1)
			sx0=cbind(1,x1,x2)	
             simry0(sx0,ni,tbss0,taut,tsat,tass,tasn,sv0[ii,kk],ii,kk)
}

stopCluster(cl)

#### true SACE value calculation

tsace=rep(0,nsim)

for(kk in 1:nsim){
da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)
tsace[kk]=mean(da1$y[which(da1$ss==1)])-
mean(da0$y[which(da0$ss==1)]) 
}
msace=mean(tsace)

xnames=c("x0","x1","x2")

#############

cl=makeCluster(nc)
registerDoParallel(cl)

#### SACE estimation: fixed effects only 

fe=foreach(kk=1:nsim,.combine='cbind',.errorhandling ="remove") %dopar%{
da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)
ffit(da1,da0,xnames)}

### bootstrap

bfe=foreach(kk=1:nsim) %:%
foreach(icount(nbp),.combine='cbind',.errorhandling ="remove") %dopar% {

da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)

m1=nrow(da1)
m0=nrow(da0)

s1=sample(m1,m1,replace=T)
s0=sample(m0,m0,replace=T)

bda1=da1[s1,]
bda0=da0[s0,]

ffit(bda1,bda0,xnames) 
}

#### SACE estimation using the ME method

me=foreach(kk=1:nsim,.combine='cbind',.errorhandling ="remove") %dopar%{
source("mecode.R")
da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)
fit(da1,da0,xnames)}

### bootstrap

bme=foreach(kk=1:nsim) %:%
foreach(icount(nbp),.combine='cbind',.errorhandling ="remove") %dopar% {
source("mecode.R")
da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)

s1=sample(nc,nc,replace=T)
s0=sample(nc,nc,replace=T)

bda1=bda0=NULL

for(ii in 1:nc){
da1j=subset(da1,cid==s1[ii])
da1j$cid=ii
bda1=rbind(bda1,da1j)}

for(ii in 1:nc){
da0k=subset(da0,cid==s0[ii])
da0k$cid=ii
bda0=rbind(bda0,da0k)}

fit(bda1,bda0,xnames) 
}

#### SACE estimation using the ME2 method

me2=foreach(kk=1:nsim,.combine='cbind',.errorhandling ="remove") %dopar%{
da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)
rfit(da1,da0,xnames)}

### bootstrap

bme2=foreach(kk=1:nsim) %:%
foreach(icount(nbp),.combine='cbind',.errorhandling ="remove") %dopar% {

da0=subset(sda0,simid==kk)
da1=subset(sda1,simid==kk)

s1=sample(nc,nc,replace=T)
s0=sample(nc,nc,replace=T)

bda1=bda0=NULL

for(ii in 1:nc){
da1j=subset(da1,cid==s1[ii])
da1j$cid=ii
bda1=rbind(bda1,da1j)}

for(ii in 1:nc){
da0k=subset(da0,cid==s0[ii])
da0k$cid=ii
bda0=rbind(bda0,da0k)}

rfit(bda1,bda0,xnames) 
}

stopCluster(cl)

################################### |Bias|, MSE and Coverage ############

fesace=fe[1,]
mesace=me[1,]
me2sace=me2[1,]

bias=round(abs(c(mean(fesace),mean(mesace),mean(me2sace))-msace)*100,2)
mse=round(c(mean((fesace-msace)^2),mean((mesace-msace)^2),mean((me2sace-msace)^2))*100,2)

fcover=mcover=m2cover=rep(0,nsim)
pb=c(0.025,0.975)

for(kk in 1:nsim){

qfe=quantile(bfe[[kk]][1,],prob=pb)
fcover[kk]=((qfe[1]<=msace)&(qfe[2]>=msace))

qme=quantile(bme[[kk]][1,],prob=pb)
mcover[kk]=((qme[1]<=msace)&(qme[2]>=msace))

qme2=quantile(bme2[[kk]][1,],prob=pb)
m2cover[kk]=((qme2[1]<=msace)&(qme2[2]>=msace))
}

cover=c(mean(fcover),mean(mcover),mean(m2cover))

rbind(bias,mse,cover*100)


