rm(list=ls())

#AMS Function
amseval <- function (pred,actual,wgt,wfactor) {
	thresh<-sort(pred)
	
	tmat<-as.data.frame(cbind(pred,wgt*(actual),wgt*(1-actual)))
		
	tmato<-tmat[order(tmat[,1]),]
	tmato[,2]<-cumsum(tmato[,2])*wfactor
	tmato[,3]<-cumsum(tmato[,3])*wfactor
	
	tmato[,2]<-tmato[nrow(tmato),2]-tmato[,2]
	tmato[,3]<-tmato[nrow(tmato),3]-tmato[,3]
	
	tmato[which(tmato[,2]<0),2]<-0
	tmato[which(tmato[,3]<0),3]<-0

	tmato<-tmato[c(1:(nrow(tmato)-1)),]
	
	amsmat<-matrix(NA,nrow=length(pred)-1,ncol=2)
	
	amsmat[,1]<-tmato[,1]
	amsmat[,2]<-sqrt(2 * ((tmato[,2] + tmato[,3] + 10) * log(1 + tmato[,2] / (tmato[,3] + 10)) - tmato[,2]))
	
	amsmat
}


load(file="HB_dev_val_samps.RData")

#Discretize all numerical variables
catvars<-c('EventId','Label','PRI_jet_num')
numvars<-setdiff(names(devsamp),catvars)
numvars<-setdiff(numvars,'Weight')

devsamp[,'PRI_jet_num']<-as.factor(devsamp[,'PRI_jet_num'])
devsamp[,'Label']<-as.factor(devsamp[,'Label'])
devsamp[,'EventId']<-as.character(devsamp[,'EventId'])

valsamp[,'PRI_jet_num']<-as.factor(valsamp[,'PRI_jet_num'])
valsamp[,'Label']<-as.factor(valsamp[,'Label'])
valsamp[,'EventId']<-as.character(valsamp[,'EventId'])

	
devsamp$Weight<-as.numeric(as.character(devsamp$Weight))
valsamp$Weight<-as.numeric(as.character(valsamp$Weight))

devwgt<-devsamp$Weight
valwgt<-valsamp$Weight


for (i in numvars) {

	devsamp[,i]<-as.numeric(as.character(devsamp[,i]))
	valsamp[,i]<-as.numeric(as.character(valsamp[,i]))

	devsamp[which(devsamp[,i]==-999),i]<-0
	valsamp[which(valsamp[,i]==-999),i]<-0
	
}

devsamp<-devsamp[-1]
valsamp<-valsamp[-1]

devsamp<-devsamp[-which(names(devsamp)=='Weight')]
valsamp<-valsamp[-which(names(valsamp)=='Weight')]


#------BUNCH OF EXPERIMENTS
#LR model
devsamp$resp<-0
devsamp$resp[which(devsamp$Label=='s')]<-1

valsamp$resp<-0
valsamp$resp[which(valsamp$Label=='s')]<-1

devsamp<-devsamp[-which(names(devsamp)=='Label')]
valsamp<-valsamp[-which(names(valsamp)=='Label')]

dev1s<-which(devsamp$PRI_jet_num==0)
dev2s<-which(devsamp$PRI_jet_num==1)
dev3s<-which(devsamp$PRI_jet_num==2)
dev4s<-which(devsamp$PRI_jet_num==3)

val1s<-which(valsamp$PRI_jet_num==0)
val2s<-which(valsamp$PRI_jet_num==1)
val3s<-which(valsamp$PRI_jet_num==2)
val4s<-which(valsamp$PRI_jet_num==3)

s1var<-setdiff(names(devsamp),c('PRI_jet_num','DER_deltaeta_jet_jet','DER_mass_jet_jet','DER_prodeta_jet_jet','DER_lep_eta_centrality','PRI_jet_leading_pt','PRI_jet_leading_eta','PRI_jet_leading_phi','PRI_jet_subleading_pt','PRI_jet_subleading_eta','PRI_jet_subleading_phi','PRI_jet_all_pt'))
s2var<-setdiff(names(devsamp),c('PRI_jet_num','DER_deltaeta_jet_jet','DER_mass_jet_jet','DER_prodeta_jet_jet','DER_lep_eta_centrality','PRI_jet_subleading_pt','PRI_jet_subleading_eta','PRI_jet_subleading_phi'))
s3var<-setdiff(names(devsamp),c('PRI_jet_num'))
s4var<-setdiff(names(devsamp),c('PRI_jet_num'))



#-use this for 2_1_1 version
frm1<-as.formula(paste("resp~",paste(setdiff(s1var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s1var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))
frm2<-as.formula(paste("resp~",paste(setdiff(s2var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s2var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))
frm3<-as.formula(paste("resp~",paste(setdiff(s3var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s3var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))
frm4<-as.formula(paste("resp~",paste(setdiff(s4var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s4var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))


require(earth)
mm1<-earth(frm1,data=devsamp[dev1s,s1var],degree=1,glm=list(family=binomial,weights=devwgt[dev1s]),trace=1,fast.beta=0,thresh=1e-6,pmethod="none")
mm2<-earth(frm2,data=devsamp[dev2s,s2var],degree=1,glm=list(family=binomial,weights=devwgt[dev2s]),trace=1,fast.beta=0,thresh=1e-6,pmethod="none")
mm3<-earth(frm3,data=devsamp[dev3s,s3var],degree=1,glm=list(family=binomial,weights=devwgt[dev3s]),trace=1,fast.beta=0,thresh=1e-6,pmethod="none")
mm4<-earth(frm4,data=devsamp[dev4s,s4var],degree=1,glm=list(family=binomial,weights=devwgt[dev4s]),trace=1,fast.beta=0,thresh=1e-6,pmethod="none")


mod1v<-0*c(1:nrow(valsamp))
mod1v[val1s]<-predict(mm1,newdata=model.frame(frm1,valsamp[val1s,s1var]), type = "response")
mod1v[val2s]<-predict(mm2,newdata=model.frame(frm2,valsamp[val2s,s2var]), type = "response")
mod1v[val3s]<-predict(mm3,newdata=model.frame(frm3,valsamp[val3s,s3var]), type = "response")
mod1v[val4s]<-predict(mm4,newdata=model.frame(frm4,valsamp[val4s,s4var]), type = "response")


amsmat<-amseval(mod1v,valsamp$resp,valwgt,(1+nrow(devsamp)/nrow(valsamp)))

########
amsmax<-max(amsmat[,2],na.rm=T)
threshmax<-amsmat[which(amsmat[,2]==amsmax),1]

plot(amsmat[,1],amsmat[,2])

save(list=c("mm1","mm2","mm3","mm4","s1var","s2var","s3var","s4var","amsmat","amsmax","threshmax"),file="mod1Params_v_2_1_2.RData")

