library(MASS)

csmfun = function(taut, tsat, ni) {
	csm = matrix(taut,ni,ni)+diag(tsat,ni) 
	return(csm)
}

### outcome simulation in the treatment group
### no random intercepts in the membership models
	
## ii is the cluster index
## kk is the simulation index

simy1=function(xi,ni,tbss,tbsn,taut,tsat,tass,tasn,ii,kk){
	covi=csmfun(taut,tsat,ni)
						
	myiss=xi%*%tbss
	yiss=mvrnorm(1,myiss,covi)	# y_ss_1	
	myisn=xi%*%tbsn
	yisn=mvrnorm(1,myisn,covi)	# y_sn 
				
	lss=xi%*%tass
        lsn=xi%*%tasn
        ess=exp(lss)
        esn=exp(lsn)
        pnni=1/(1+ess+esn)      # p_nn
        pssi=ess*pnni           # p_ss
        psni=esn*pnni           # p_sn

        sta=matrix(0,ni,3)
        for(j in 1:ni){
        sta[j,]=rmultinom(n=1,size=1,prob=c(pssi[j],psni[j],pnni[j]))
		}
	yi=rowSums(cbind(yiss,yisn)*sta[,1:2])
	yind=rowSums(sta[,1:2])  # index of observed outcomes
	dai=data.frame(y=yi,yss=yiss,ysn=yisn,pss=pssi,psn=psni,pnn=pnni,
	ss=sta[,1],sn=sta[,2],nn=sta[,3],id=1:ni,cid=ii,simid=kk,x0=1,
	x1=xi[,"x1"],x2=xi[,"x2"],yind=yind)		

return(dai)
}

### outcome simulation in the control group
### no random intercepts in the membership models

simy0=function(xi,ni,tbss0,taut,tsat,tass,tasn,ii,kk){

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

### outcome simulation in the treatment group
### Random intercepts are in both the membership models and outcome models

simry1=function(xi,ni,tbss,tbsn,taut,tsat,tass,tasn,v1i,ii,kk){

	covi=csmfun(taut,tsat,ni)
						
	myiss=xi%*%tbss
	yiss=mvrnorm(1,myiss,covi)		
	myisn=xi%*%tbsn
	yisn=mvrnorm(1,myisn,covi)	
				
	lss=xi%*%tass+v1i
        lsn=xi%*%tasn+v1i
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

### outcome simulation in the control group
### Random intercepts are in both the membership models and outcome models

simry0=function(xi,ni,tbss0,taut,tsat,tass,tasn,v0i,ii,kk){

	covi=csmfun(taut,tsat,ni)
						
	myi=xi%*%tbss0
	yi=mvrnorm(1,myi,covi)			

	lss=xi%*%tass+v0i
        lsn=xi%*%tasn+v0i
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