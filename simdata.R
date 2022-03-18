library(MASS)

csmfun = function(taut, tsat, n) {
	csm = matrix(taut,n,n)+diag(tsat,n) 
	return(csm)
}

simn=function(msize,sdsize,nc,nsim){
	rn=matrix(round(rnorm(nc*nsim,msize,sdsize)),nc,nsim)
	return(rn)
}

simx=function(rn){
nsim=ncol(rn)
nc=nrow(rn)
sx=NULL	
	for(k in 1:nsim){

		for(i in 1:nc){
			ni=rn[i,k]
			x1=rbinom(ni,1,1/2)
			x2=rnorm(ni,0,1)
			xi=data.frame(x0=1,x1=x1,x2=x2,id=1:ni,cid=i,simid=k)	
			sx=rbind(sx,xi)}}	
	return(sx)
	}
	
simy1=function(sx,kk,ii,tbss,tbsn,taut,tsat,tass,tasn){
library(MASS)
xk=subset(sx,simid==kk)
			xi=as.matrix(xk[which(xk$cid==ii),c("x0","x1","x2")])
			ni=nrow(xi)
			covi=csmfun(taut,tsat,ni)
						
			myiss=xi%*%tbss
			yiss=mvrnorm(1,myiss,covi)		
			myisn=xi%*%tbsn
			yisn=mvrnorm(1,myisn,covi)					
		lss=xi%*%tass
        lsn=xi%*%tasn
        ess=exp(lss)
        esn=exp(lsn)
        pnni=1/(1+ess+esn)
        pssi=ess*pnni
        psni=esn*pnni
        sta=matrix(0,ni,3)
        for(j in 1:ni){
        sta[j,]=rmultinom(n=1,size=1,prob=c(pssi[j],psni[j],pnni[j]))
		}
		yi=rowSums(cbind(yiss,yisn)*sta[,1:2])
		yind=rowSums(sta[,1:2])
		dai=data.frame(y=yi,yss=yiss,ysn=yisn,pss=pssi,psn=psni,pnn=pnni,
		ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=ii,simid=kk,x0=1,
		x1=xi[,"x1"],x2=xi[,"x2"],yind=yind)		

return(dai)
}

	#colnames(sda)=c("y","yss","ysn",
	#"pss","psn","pnn","ss","sn","nn","id","cid","simid","x0","x1","x2",
	#"survivor")

##########

simy0=function(sx0,kk,ii,tbss0,taut,tsat,tass,tasn){

		xk=subset(sx0,simid==kk)

	
			xi=as.matrix(xk[which(xk$cid==ii),c("x0","x1","x2")])
			ni=nrow(xi)
			covi=csmfun(taut,tsat,ni)
						
			myi=xi%*%tbss0
			yi=mvrnorm(1,myi,covi)			

		lss=xi%*%tass
        lsn=xi%*%tasn
        ess=exp(lss)
        esn=exp(lsn)
        pnni=1/(1+ess+esn)
        pssi=ess*pnni
        psni=esn*pnni
        
        sta=matrix(0,ni,3)
        for(j in 1:ni){
        sta[j,]=rmultinom(n=1,size=1,prob=c(pssi[j],psni[j],pnni[j]))
		}
		
		dai=data.frame(y=yi*sta[,1],yss=yi,pss=pssi,psn=psni,pnn=pnni,
		ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=ii,simid=kk,x0=1,
		x1=xi[,"x1"],x2=xi[,"x2"],yind=sta[,1])		

return(dai)
}

library(doParallel)

tbss = c(-0.5,1, 1.5)
p=length(tbss)

tbsn=c(-0.3,0.8,1.3)
tbss0=c(-0.2,1,1)

tass=c(1,2,1)
tasn=c(-0.5,-1.5,-1)

################################

nc = c(30,60)
lnc = length(nc)
msize = c(25, 50)
lms = length(msize)

sdsize = 3

tt = 2
trho = c(0.01,0.05,0.1)
lr = length(trho)

################################\n

nsim = 200


for (j in 1:lnc) {

	for (k in 1:lms) {

		rn = simn(msize[k], sdsize, nc[j], nsim)

		for (i in 1:lr) {

cl=makeCluster(detectCores())
registerDoParallel(cl)

sx=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc[j],.combine='rbind') %dopar%{

			ni=rn[ii,kk]
			x1=rbinom(ni,1,1/2)
			x2=rnorm(ni,0,1)
			data.frame(x0=1,x1=x1,x2=x2,id=1:ni,cid=ii,simid=kk)	
			}

sx0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc[j],.combine='rbind') %dopar%{

			ni=rn[ii,kk]
			x1=rbinom(ni,1,1/2)
			x2=rnorm(ni,0,1)
			data.frame(x0=1,x1=x1,x2=x2,id=1:ni,cid=ii,simid=kk)	
			}


taut=tt*trho[i]
tsat=tt*(1-trho[i])

sda1=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc[j],.combine='rbind') %dopar%{

simy1(sx,kk,ii,tbss,tbsn,taut,tsat,tass,tasn)

}

sda0=foreach(kk=1:nsim,.combine='rbind') %:%
foreach(ii=1:nc[j],.combine='rbind') %dopar%{

simy0(sx0,kk,ii,tbss0,taut,tsat,tass,tasn)

}

stopCluster(cl)


			simo = list(sda1=sda1,sda0=sda0)

			save(simo, file = paste("nc_", nc[j], "m_", msize[k], "rho_", i, ".RData", 
				sep = ""))

		}
	}
}

po = list(trho = trho, tbss = tbss, tbsn=tbsn,tt=tt,
tbss0=tbss0,lr = lr, nsim = nsim, nc = nc, p=p,tass=tass,tasn=tasn, 
	lnc = lnc, msize = msize, lms = lms, taut = tt * trho, tsat=tt*(1-trho))

save(po, file = paste("par.RData"))
