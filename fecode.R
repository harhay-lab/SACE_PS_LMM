pmod=function(xi,alphass,alphasn){

ess=exp(xi%*%alphass)
esn=exp(xi%*%alphasn)
dsn=(1+ess+esn)
pnn=1/dsn

pss=ess*pnn
psn=esn*pnn	

sdess=pss*(1-pss)
sdesn=psn*(1-psn)

eta0=esn/(esn+1)

return(list(pss=pss,psn=psn,pnn=pnn,
sdesn=sdesn,sdess=sdess,eta0=eta0))	
}

feta1=function(x11,y11,betass1,betasn,sat,pss1,psn1){

sa=sqrt(sat)

mss=x11%*%betass1
msn=x11%*%betasn

nume=dnorm(y11,mean=mss,sd=sa)*pss1
dsn=dnorm(y11,mean=msn,sd=sa)*psn1
eta1=nume/(nume+dsn)
return(eta1)

}

ftr=function(x1,tx1,x11,tx11,sind1,y11,m1,m11,
alphass,alphasn,betass1,betasn,sat,p){

######### E step

pmod1=pmod(x1,alphass,alphasn)
pss=pmod1$pss
psn=pmod1$psn
sdss=pmod1$sdess
sdsn=pmod1$sdesn
pnn=pmod1$pnn

pss1=pss[sind1]
psn1=psn[sind1]

eta=feta1(x11,y11,betass1,betasn,sat,pss1,psn1)

########### survivors 

#### beta

meta=matrix(eta,nrow=p,ncol=m11,byrow=T)
etax11=tx11*meta
eta1x11=tx11-etax11

bssp=etax11%*%x11
bssy=etax11%*%y11

bsnp=eta1x11%*%x11
bsny=eta1x11%*%y11

##### alpha

aess=rowSums(etax11)
aesn=rowSums(eta1x11)

############### all

mpss=matrix(pss,nrow=p,ncol=m1,byrow=T)
msdss=matrix(sdss,nrow=p,ncol=m1,byrow=T)

dass1=-rowSums(tx1*mpss)
dass2=-(tx1*msdss)%*%x1

mpsn=matrix(psn,nrow=p,ncol=m1,byrow=T)
msdsn=matrix(sdsn,nrow=p,ncol=m1,byrow=T)
	
dasn1=-rowSums(tx1*mpsn)
dasn2=-(tx1*msdsn)%*%x1

nbetass1=lm.fit(bssp,bssy)$coef
nbetasn=lm.fit(bsnp,bsny)$coef

return(list(nbetass1=nbetass1,nbetasn=nbetasn,dass1=dass1,
dass2=dass2,dasn1=dasn1,dasn2=dasn2,
aess=aess,aesn=aesn))
}

fcon=function(x0,tx0,x01,tx01,x00,tx00,sind0,y01,m0,m01,
alphass,alphasn,betass0,sat,p){
	
m00=m0-m01

pmod0=pmod(x0,alphass,alphasn)
pss0=pmod0$pss
psn0=pmod0$psn
sdss0=pmod0$sdess
sdsn0=pmod0$sdesn
pnn0=pmod0$pnn

eta=pmod0$eta0[-sind0]
	
########### non-survivors

##### alpha

meta=matrix(eta,nrow=p,ncol=m00,byrow=T)
aesnc=rowSums(tx00*meta)

########### survivors 

#### beta

bssp0=tx01%*%x01
bssy0=tx01%*%y01

##### alpha

aessc=rowSums(tx01)

############### all

###### alpha

mpss0=matrix(pss0,nrow=p,ncol=m0,byrow=T)
msdss0=matrix(sdss0,nrow=p,ncol=m0,byrow=T)

dass1c=-rowSums(tx0*mpss0)
dass2c=-(tx0*msdss0)%*%x0

mpsn0=matrix(psn0,nrow=p,ncol=m0,byrow=T)
msdsn0=matrix(sdsn0,nrow=p,ncol=m0,byrow=T)
	
dasn1c=-rowSums(tx0*mpsn0)
dasn2c=-(tx0*msdsn0)%*%x0


nbetass0=lm.fit(bssp0,bssy0)$coef

return(list(nbetass0=nbetass0,dass1c=dass1c,dass2c=dass2c,
dasn1c=dasn1c,dasn2c=dasn2c,aessc=aessc,aesnc=aesnc))
}

#####

ffe=function(x1,tx1,x11,tx11,sind1,y11,m1,m11,
x0,tx0,x01,tx01,x00,tx00,sind0,y01,m0,m01,
alphass,alphasn,betass1,betasn,sat,p){

to=ftr(x1,tx1,x11,tx11,sind1,y11,m1,m11,
alphass,alphasn,betass1,betasn,sat,p)

co=fcon(x0,tx0,x01,tx01,x00,tx00,sind0,y01,m0,m01,
alphass,alphasn,betass0,sat,p)

nalphass=alphass-solve(to$dass2+co$dass2c,to$dass1+co$dass1c+to$aess+co$aessc)

nalphasn=alphasn-solve(to$dasn2+co$dasn2c,to$dasn1+co$dasn1c+to$aesn+co$aesnc)

return(list(nbetass1=to$nbetass1,nbetasn=to$nbetasn,
nbetass0=co$nbetass0,
nalphass=nalphass,nalphasn=nalphasn))
}

#################################

fv=function(x1,x11,sind1,y11,m11,x0,x01,sind0,y01,m01,
betass1,betasn,betass0,alphass,alphasn,sat){

######### E step

pmod1=pmod(x1,alphass,alphasn)
pss=pmod1$pss
psn=pmod1$psn
pss1=pss[sind1]
psn1=psn[sind1]

eta=feta1(x11,y11,betass1,betasn,sat,pss1,psn1)
eta1=1-eta

ty1ss=eta*(y11-x11%*%betass1)^2
ty1sn=eta1*(y11-x11%*%betasn)^2
den1=sum(ty1ss+ty1sn)

ty0ss=(y01-x01%*%betass0)^2
den0=sum(ty0ss)

nsat=(den1+den0)/(m11+m01)

return(nsat)
}

#######################

ffit=function(da1,da0,xnames){

x1=as.matrix(da1[,xnames]) # m_{1} x p
tx1=t(x1) # p x m_{1}
m1=nrow(x1)
p=ncol(x1)
y1=da1[,"y"]
sind1=which(da1$yind==1)	    
m11=length(sind1)
x11=x1[sind1,,drop=F] # m_{1,1} x p
tx11=t(x11)  # p x m_{1,1}
y11=y1[sind1]

x0=as.matrix(da0[,xnames]) # m_{0} x p
tx0=t(x0) # p x m_{0}
m0=nrow(x0)
y0=da0[,"y"]
sind0=which(da0$yind==1)
m01=length(sind0)

x01=x0[sind0,,drop=F] # m_{0,1} x p
tx01=t(x01)  # p x m_{0,1}
y01=y0[sind0]

x00=x0[-sind0,,drop=F] # m_{0,0} x p
tx00=t(x00)  # p x m_{0,0}

######
 
alphass=(1:p)/(100*p)
alphasn=(p:1)/(90*p)

s1=lm(y1~x1-1)
s0=lm(y0~x0-1)

betass1=s1$coef
betass0=s0$coef
betasn=(betass1+betass0)/2
sat=(summary(s1)$sigma^2+summary(s0)$sigma^2)/2

dif=1
nit=40
thr=1e-3
nuit=0

while(any(dif>thr)&(nuit<nit)){

fo=ffe(x1,tx1,x11,tx11,sind1,y11,m1,m11,
x0,tx0,x01,tx01,x00,tx00,sind0,y01,m0,m01,
alphass,alphasn,betass1,betasn,sat,p)

dif.betass1=betass1-fo$nbetass1
dif.betasn=betasn-fo$nbetasn
dif.betass0=betass0-fo$nbetass0
dif.beta=c(dif.betass1,dif.betasn,dif.betass0)

betass1=fo$nbetass1
betasn=fo$nbetasn
betass0=fo$nbetass0

dif.alphass=alphass-fo$nalphass
dif.alphasn=alphasn-fo$nalphasn
dif.alpha=c(dif.alphass,dif.alphasn)

alphass=fo$nalphass
alphasn=fo$nalphasn

vo=fv(x1,x11,sind1,y11,m11,x0,x01,sind0,y01,m01,
betass1,betasn,betass0,alphass,alphasn,sat)

dif.v=c(sat-vo)
sat=vo

dif=abs(c(dif.v,dif.alpha,dif.beta))
nuit=nuit+1
}

po0=pmod(x0,alphass,alphasn)
fpss0=po0$pss
fy0=x0%*%betass0

po1=pmod(x1,alphass,alphasn)
fpss1=po1$pss
fy1=x1%*%betass1

fsace=sum(fy1*fpss1)/sum(fpss1)-
sum(fy0*fpss0)/sum(fpss0)

pss=(sum(fpss0)+sum(fpss1))/(m1+m0)
psn=(sum(po0$psn)+sum(po1$psn))/(m1+m0)
pnn=1-pss-psn

return(c(betass1,betasn,betass0,
alphass,alphasn,sat,fsace,pss,psn,pnn))
}

sace=function(da,y.name,tr.name,xnames){

da$y=da[,y.name]
da$yind=1
da$yind[is.na(da$y)]=0

da0=da[which(da[,tr.name]==0),]
da1=da[which(da[,tr.name]==1),]

fobj=ffit(da1,da0,xnames)
return(fobj)
}
