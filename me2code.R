##### SACE estimation 
## with random intercepts in both the outcome models and the membership models

### da: data set
### y.name: name of the outcome variable
### tr.name: name of the treatment variable
### xnames: names of the covariates
### cluster.name: name of the cluster variable

rsace=function(da,y.name,tr.name,xnames,cluster.name){

da$y=da[,y.name]
da$yind=1  # create the index variable.
da$yind[is.na(da$y)]=0 # Index values of the truncated outcomes (NA's) are set 0

da0=da[which(da[,tr.name]==0),] # control arm
da1=da[which(da[,tr.name]==1),] # treatment arm

da0$cid=da0[,cluster.name]
ugp0=sort(unique(da0[,cluster.name]))
lugp0=length(ugp0)
for(i in 1:lugp0){
da0$cid[which(da0[,cluster.name]==ugp0[i])]=i
}

da1$cid=da1[,cluster.name]  # create the cluster variable "cid" for the treatment group
ugp1=sort(unique(da1[,cluster.name]))
lugp1=length(ugp1)
for(i in 1:lugp1){
da1$cid[which(da1[,cluster.name]==ugp1[i])]=i
}

fobj=rfit(da1,da0,xnames)
return(fobj)
}

#######################

rfit=function(da1,da0,xnames){

nms=100  # number of Monte Carlo simulations in the E-step calculation

### treatment group 

n1=length(unique(da1$cid)) # number of treated clusters
m1v=as.vector(table(da1$cid)) # cluster size
m1=sum(m1v) # sample size
m1v1=as.vector(table(da1$cid,da1$yind)[,2]) # number of observed outcomes by cluster
xm1=da1[,c("cid",xnames),drop=F]
p=ncol(xm1)-1
y1=da1[,"y"]
sind1=tapply(da1$yind,da1$cid,function(x){which(x==1)})	 # index of observed outcomes by cluster    

### control group 

n0=length(unique(da0$cid))
m0v=as.vector(table(da0$cid))
m0=sum(m0v)
m0v1=as.vector(table(da0$cid,da0$yind)[,2])
xm0=da0[,c("cid",xnames),drop=F]
y0=da0[,"y"]
sind0=tapply(da0$yind,da0$cid,function(x){which(x==1)})	 
	  
### starting values 
  
alphass=t(t((1:p)/p))/20
alphasn=-t(t((1:p)/p))/10

s1=lm(y1~as.matrix(xm1[,-1])-1)
s0=lm(y0~as.matrix(xm0[,-1])-1)

betass1=s1$coef
betass0=s0$coef
betasn=(betass1+betass0)/2

sat=(summary(s1)$sigma^2+summary(s0)$sigma^2)/2
taut=sat/5
gat=sat/10

dif=1
nit=35 # maximum number of iterations
thr=1e-3
nuit=0

while(any(dif>thr)&(nuit<nit)){

fo=rfe(xm1,sind1,y1,n1,m1v,m1v1,
alphass,alphasn,betass1,betasn,
xm0,sind0,y0,n0,m0v,m0v1,m01,
betass0,taut,sat,gat,p,nms,xnames)

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

vo=rv(xm1,sind1,y1,n1,m1v,m1v1,
xm0,sind0,y0,n0,m0v,m0v1,
betass1,betasn,betass0,alphass,alphasn,
taut,sat,gat,nms,m0,m1,xnames)

dif.v=c(sat-vo$nsat,taut-vo$ntaut,vo$ngat-gat)
sat=vo$nsat
taut=vo$ntaut
gat=vo$ngat

dif=abs(c(dif.v,dif.alpha,dif.beta))
nuit=nuit+1
}

return(obj=c(vo$sace,betass1,betasn,betass0,
alphass,alphasn,sat,taut,gat,vo$pss,vo$psn,vo$pnn))

}

####################

rpmod=function(xi,sindi,alphass,alphasn,gat,ni,nms){

svi=rnorm(nms,0,sqrt(gat))

pss=psn=pnn=sdess=sdesn=matrix(0,ni,nms)
fz1vi=fz0vi=rep(0,nms)

indi=rep(0,ni)
indi[sindi]=1
nindi=1-indi

for(k in 1:nms){

ess=exp(xi%*%alphass+svi[k])
esn=exp(xi%*%alphasn+svi[k])
dsn=(1+ess+esn)

pnn[,k]=pnnk=1/dsn
fz1vi[k]=prod(((1-pnnk)^indi)*(pnnk^nindi))

pss[,k]=pssk=ess*pnnk
fz0vi[k]=prod(((1-pssk)^nindi)*(pssk^indi))

psn[,k]=psnk=esn*pnnk	

sdess[,k]=pssk*(1-pssk)
sdesn[,k]=psnk*(1-psnk)
}

mpss=rowMeans(pss)
mpsn=rowMeans(psn)
mpnn=rowMeans(pnn)

fz1=mean(fz1vi)
if(fz1==0){ev1i2=0
msdess1=rowMeans(sdess)
msdesn1=rowMeans(sdesn)
} else{ev1i2=mean(svi^2*fz1vi)/fz1
msdess1=rowMeans(sdess*fz1vi)/fz1
msdesn1=rowMeans(sdesn*fz1vi)/fz1}

fz0=mean(fz0vi)
if(fz0==0){
ev0i2=0
msdess0=rowMeans(sdess)
msdesn0=rowMeans(sdesn)
} else{
ev0i2=mean(svi^2*fz0vi)/fz0
msdess0=rowMeans(sdess*fz0vi)/fz0
msdesn0=rowMeans(sdesn*fz0vi)/fz0}

meta0=mpsn/(mpsn+mpnn)

return(list(mpss=mpss,mpsn=mpsn,mpnn=mpnn,meta0=meta0,ev1i2=ev1i2,ev0i2=ev0i2,
msdesn1=msdesn1,msdess1=msdess1,msdesn0=msdesn0,msdess0=msdess0))	
}

reta1=function(x1i1,y1i,betass1,betasn,taut,sat,mpss1,mpsn1){

mss=x1i1%*%betass1
msn=x1i1%*%betasn
jointsd=sqrt(sat+taut)

nume=dnorm(y1i,mean=mss,sd=jointsd)*mpss1
dsn=dnorm(y1i,mean=msn,sd=jointsd)*mpsn1
eta1i=nume/(nume+dsn)
return(eta1i)
}

reu1=function(x1i1,y1i,betass1,betasn,taut,sat,mpss1,mpsn1,nms){

mss=x1i1%*%betass1
msn=x1i1%*%betasn

sdtau=sqrt(taut)
sdsa=sqrt(sat) 

sui=rnorm(nms,mean=0,sd=sdtau)
fyij=rep(0,nms)

for(k in 1:nms){
	uk=sui[k]
	d1=dnorm(y1i,mean=mss+uk,sd=sdsa)*mpss1
	d2=dnorm(y1i,mean=msn+uk,sd=sdsa)*mpsn1
	fij=d1+d2
	fyij[k]=prod(fij)}

fyi=mean(fyij)

if(fyi==0){
eui=eui2=vui=0
} else{
eui=mean(sui*fyij)/fyi
eui2=mean(sui^2*fyij)/fyi
vui=eui2-(eui)^2}

return(list(vui=vui,eui2=eui2,eui=eui))
}

reu0=function(x0i1,y0i,m0i1,betass0,taut,sat){

mss0=x0i1%*%betass0

coef=taut/(m0i1*taut+sat)
eui=coef*sum(y0i-mss0)
vui=coef*sat
eui2=vui+eui^2
return(list(vui=vui,eui2=eui2,eui=eui))
}

rtr=function(xm1,sind1,y1,n1,m1v,m1v1,
alphass,alphasn,betass1,betasn,taut,sat,gat,p,nms,xnames){
	
bssp=bsnp=matrix(0,p,p)
bssy=bsny=matrix(0,p,1)
dass1=dasn1=matrix(0,p,1)
dass2=dasn2=matrix(0,p,p)
aess=aesn=matrix(0,p,1)

for(i in 1:n1){
	
cind=which(xm1$cid==i)
sind1i=sind1[[i]]
	
x1i=as.matrix(xm1[cind,xnames,drop=F]) # m_{1,i} x p
tx1i=t(x1i) # p x m_{1,i}
m1i1=m1v1[i]
m1i=m1v[i]

######### E step

pmod1=rpmod(x1i,sind1i,alphass,alphasn,gat,m1i,nms)

mpss=pmod1$mpss
mpsn=pmod1$mpsn
msdss=pmod1$msdess1
msdsn=pmod1$msdesn1
mpnn=pmod1$mpnn

mpss1=mpss[sind1i]
mpsn1=mpsn[sind1i]

if(m1i1>0){

x1i1=x1i[sind1i,,drop=F] # m_{1,i,1} x p
tx1i1=t(x1i1)  # p x m_{1,i,1}

y1i=y1[cind][sind1i]

etai=reta1(x1i1,y1i,betass1,betasn,taut,sat,mpss1,mpsn1)

########### survivors 

euo=reu1(x1i1,y1i,betass1,betasn,taut,sat,mpss1,mpsn1,nms)
cyi=y1i-euo$eui

#### beta

metai=matrix(etai,nrow=p,ncol=m1i1,byrow=T)
etaix1i1=tx1i1*metai
etai1x1i1=tx1i1-etaix1i1

bssp=bssp+etaix1i1%*%x1i1
bssy=bssy+etaix1i1%*%cyi

bsnp=bsnp+etai1x1i1%*%x1i1
bsny=bsny+etai1x1i1%*%cyi

##### alpha

aess=aess+rowSums(etaix1i1)
aesn=aesn+rowSums(etai1x1i1)
} else{bssp=bssp
       bssy=bssy
       bsnp=bsnp
       bsny=bsny
       aess=aess
       aesn=aesn}

############### all

mpssm=matrix(mpss,nrow=p,ncol=m1i,byrow=T)
msdssm=matrix(msdss,nrow=p,ncol=m1i,byrow=T)

dass1=dass1-rowSums(tx1i*mpssm)
dass2=dass2-(tx1i*msdssm)%*%x1i

mpsnm=matrix(mpsn,nrow=p,ncol=m1i,byrow=T)
msdsnm=matrix(msdsn,nrow=p,ncol=m1i,byrow=T)
	
dasn1=dasn1-rowSums(tx1i*mpsnm)
dasn2=dasn2-(tx1i*msdsnm)%*%x1i

}

nbetass1=lm.fit(bssp,bssy)$coef
nbetasn=lm.fit(bsnp,bsny)$coef

return(list(nbetass1=nbetass1,nbetasn=nbetasn,dass1=dass1,
dass2=dass2,dasn1=dasn1,dasn2=dasn2,aess=aess,aesn=aesn))
}


rcon=function(xm0,sind0,y0,n0,m0v,m0v1,
alphass,alphasn,betass0,taut,sat,gat,p,nms,xnames){
	
bssp0=matrix(0,p,p)
bssy0=matrix(0,p,1)
dasn1c=dass1c=matrix(0,p,1)
dasn2c=dass2c=matrix(0,p,p)
aessc=aesnc=matrix(0,p,1)

for(i in 1:n0){

cind=which(xm0$cid==i)	
x0i=as.matrix(xm0[cind,xnames,drop=F]) # m_{0,i} x p

tx0i=t(x0i) # p x m_{0,i}
sind0i=sind0[[i]]

m0i1=m0v1[i]
m0i=m0v[i]
m0i0=m0i-m0i1

pmod0=rpmod(x0i,sind0i,alphass,alphasn,gat,m0i,nms)
mpss0=pmod0$mpss
mpsn0=pmod0$mpsn
msdss0=pmod0$msdess0
msdsn0=pmod0$msdesn0
mpnn0=pmod0$mpnn
	
########### non-survivors

if(m0i0>0){

##### alpha

if(m0i1>0){
x0i0=x0i[-sind0i,,drop=F] # m_{0,i,0} x p
etai=pmod0$meta0[-sind0i]} else {
x0i0=x0i
etai=pmod0$meta0
}

tx0i0=t(x0i0)  # p x m_{0,i,0} 
metai=matrix(etai,nrow=p,ncol=m0i0,byrow=T)
p0=rowSums(tx0i0*metai)
aesnc=aesnc+p0 } else {aesnc=aesnc}

########### survivors 

if(m0i1>0){

x0i1=x0i[sind0i,,drop=F] # m_{0,i,1} x p
tx0i1=t(x0i1)  # p x m_{0,i,1}

y0i=y0[cind][sind0i]

eu0=reu0(x0i1,y0i,m0i1,betass0,taut,sat)

cyi=y0i-eu0$eui

#### beta

bssp0=bssp0+tx0i1%*%x0i1
bssy0=bssy0+tx0i1%*%cyi

##### alpha

aessc=aessc+rowSums(tx0i1)
} else{bssp0=bssp0
       bssy0=bssy0
       aessc=aessc
}

############### all

###### alpha

mpss0m=matrix(mpss0,nrow=p,ncol=m0i,byrow=T)
msdss0m=matrix(msdss0,nrow=p,ncol=m0i,byrow=T)

dass1c=dass1c-rowSums(tx0i*mpss0m)
dass2c=dass2c-(tx0i*msdss0m)%*%x0i

mpsn0m=matrix(mpsn0,nrow=p,ncol=m0i,byrow=T)
msdsn0m=matrix(msdsn0,nrow=p,ncol=m0i,byrow=T)
	
dasn1c=dasn1c-rowSums(tx0i*mpsn0m)
dasn2c=dasn2c-(tx0i*msdsn0m)%*%x0i

}

nbetass0=lm.fit(bssp0,bssy0)$coef

return(list(nbetass0=nbetass0,dass1c=dass1c,dass2c=dass2c,
dasn1c=dasn1c,dasn2c=dasn2c,aessc=aessc,aesnc=aesnc))
}

#####

rfe=function(xm1,sind1,y1,n1,m1v,m1v1,
alphass,alphasn,betass1,betasn,
xm0,sind0,y0,n0,m0v,m0v1,m01,
betass0,taut,sat,gat,p,nms,xnames){

to=rtr(xm1,sind1,y1,n1,m1v,m1v1,
alphass,alphasn,betass1,betasn,taut,sat,gat,p,nms,xnames)

co=rcon(xm0,sind0,y0,n0,m0v,m0v1,
alphass,alphasn,betass0,taut,sat,gat,p,nms,xnames)

nalphass=alphass-solve(to$dass2+co$dass2c,to$dass1+co$dass1c+to$aess+co$aessc)

nalphasn=alphasn-solve(to$dasn2+co$dasn2c,to$dasn1+co$dasn1c+to$aesn+co$aesnc)

return(list(nbetass1=to$nbetass1,nbetasn=to$nbetasn,
nbetass0=co$nbetass0,
nalphass=nalphass,nalphasn=nalphasn))
}

#################################

rv=function(xm1,sind1,y1,n1,m1v,m1v1,
xm0,sind0,y0,n0,m0v,m0v1,
betass1,betasn,betass0,alphass,alphasn,
taut,sat,gat,nms,m0,m1,xnames){

den1=den0=0
eu02=ev02=rep(0,n0)
eu12=ev12=rep(0,n1)

m11=sum(m1v1)
m01=sum(m0v1)

pssi1=fyi1=psni1=0

for(i in 1:n1){

cind=which(xm1$cid==i)
	
x1i=as.matrix(xm1[cind,xnames,drop=F])# m_{1,i} x p
sind1i=sind1[[i]]
x1i1=x1i[sind1i,,drop=F]
y1i=y1[cind][sind1i]
m1i1=m1v1[i]

######### E step

pmod1=rpmod(x1i,sind1i,alphass,alphasn,gat,m1v[i],nms)
mpss=pmod1$mpss
mpsn=pmod1$mpsn
mpss1=mpss[sind1i]
mpsn1=mpsn[sind1i]
ev12[i]=pmod1$ev1i2

etai=reta1(x1i1,y1i,betass1,betasn,taut,sat,mpss1,mpsn1)
etai1=1-etai

euo=reu1(x1i1,y1i,betass1,betasn,taut,sat,mpss1,mpsn1,nms)
eu1i=euo$eui
eu12[i]=euo$eui2
ty1ss=etai*(y1i-x1i1%*%betass1-eu1i)^2
ty1sn=etai1*(y1i-x1i1%*%betasn-eu1i)^2
den1=den1+sum(ty1ss+ty1sn+euo$vui)

fyi1=fyi1+sum(mpss*(x1i%*%betass1+eu1i))
pssi1=pssi1+sum(mpss)
psni1=psni1+sum(mpsn)
}

pssi0=fyi0=psni0=0

for(i in 1:n0){

cind=which(xm0$cid==i)
	
x0i=as.matrix(xm0[cind,xnames,drop=F])
sind0i=sind0[[i]]
x0i1=x0i[sind0i,,drop=F]
y0i=y0[cind][sind0i]
m0i1=m0v1[i]

eu0=reu0(x0i1,y0i,m0i1,betass0,taut,sat)
eu02[i]=eu0$eui2
ty0ss=(y0i-x0i1%*%betass0-eu0$eui)^2
den0=den0+sum(ty0ss+eu0$vui)

pmod0=rpmod(x0i,sind0i,alphass,alphasn,gat,m0v[i],nms)
mpss0=pmod0$mpss
ev02[i]=pmod0$ev0i2

psni0=psni0+sum(pmod0$mpsn)

fyi0=fyi0+sum(mpss0*(x0i%*%betass0+eu0$eui))
pssi0=pssi0+sum(mpss0)
}

nsat=(den1+den0)/(m11+m01)

ntaut=(sum(eu12)+sum(eu02))/(n1+n0)

ngat=(sum(ev12)+sum(ev02))/(n1+n0)

sace=fyi1/pssi1-fyi0/pssi0
pss=(pssi0+pssi1)/(m1+m0)
psn=(psni0+psni1)/(m1+m0)
pnn=1-pss-psn
return(list(nsat=nsat,ntaut=ntaut,sace=sace,
ngat=ngat,pss=pss,psn=psn,pnn=pnn))
}

