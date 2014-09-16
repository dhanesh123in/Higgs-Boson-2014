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
	
}

#Missing Value Imputation for DER_mass_MMC
#	devsamp[which(devsamp[,'DER_mass_MMC']==-999),'DER_mass_MMC']<-0
#	valsamp[which(valsamp[,'DER_mass_MMC']==-999),'DER_mass_MMC']<-0


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

#Cases with missing DER_mass_MMC
dmiind<-which(devsamp[,'DER_mass_MMC']==-999)
vmiind<-which(valsamp[,'DER_mass_MMC']==-999)

#Just impute mass to zero
devsamp[dmiind,'DER_mass_MMC']<-0
valsamp[vmiind,'DER_mass_MMC']<-0

##Log-variables 
lvars<-c('DER_mass_MMC','DER_mass_transverse_met_lep','DER_mass_vis','DER_pt_h','DER_deltar_tau_lep','DER_pt_tot','DER_sum_pt','DER_pt_ratio_lep_tau','PRI_tau_pt','PRI_lep_pt','PRI_met','PRI_met_sumet','PRI_jet_leading_pt','PRI_jet_all_pt')

for (xn in intersect(lvars,s1var)) {
	devsamp[dev1s,xn]<-log(devsamp[dev1s,xn]+1)
	valsamp[val1s,xn]<-log(valsamp[val1s,xn]+1)
}

for (xn in intersect(lvars,s2var)) {
	devsamp[dev2s,xn]<-log(devsamp[dev2s,xn]+1)
	valsamp[val2s,xn]<-log(valsamp[val2s,xn]+1)
}

for (xn in intersect(lvars,s3var)) {
	devsamp[dev3s,xn]<-log(devsamp[dev3s,xn]+1)
	valsamp[val3s,xn]<-log(valsamp[val3s,xn]+1)
}

for (xn in intersect(lvars,s4var)) {
	devsamp[dev4s,xn]<-log(devsamp[dev4s,xn]+1)
	valsamp[val4s,xn]<-log(valsamp[val4s,xn]+1)
}



###8-Earth Models for PRI_jet_num's and 

require(h2o)
localH2O<-h2o.init()

dev1samp.hex<-as.h2o(localH2O,devsamp[dev1s,s1var])
dev2samp.hex<-as.h2o(localH2O,devsamp[dev2s,s2var])
dev3samp.hex<-as.h2o(localH2O,devsamp[dev3s,s3var])
dev4samp.hex<-as.h2o(localH2O,devsamp[dev4s,s4var])

val1samp.hex<-as.h2o(localH2O,valsamp[val1s,s1var])
val2samp.hex<-as.h2o(localH2O,valsamp[val2s,s2var])
val3samp.hex<-as.h2o(localH2O,valsamp[val3s,s3var])
val4samp.hex<-as.h2o(localH2O,valsamp[val4s,s4var])


h2o.mm1a<-h2o.gbm(x=1:(length(s1var)-1),y=length(s1var),data=dev1samp.hex,validation=val1samp.hex,balance.classes=T,n.trees=200)
h2o.mm2a<-h2o.gbm(x=1:(length(s2var)-1),y=length(s2var),data=dev2samp.hex,validation=val2samp.hex,balance.classes=T,n.trees=200)
h2o.mm3a<-h2o.gbm(x=1:(length(s3var)-1),y=length(s3var),data=dev3samp.hex,validation=val3samp.hex,balance.classes=T,n.trees=200)
h2o.mm4a<-h2o.gbm(x=1:(length(s4var)-1),y=length(s4var),data=dev4samp.hex,validation=val4samp.hex,balance.classes=T,n.trees=200)
#--2nd try 3.659 (original dataset), 3.717 with log transforms -- BEST


#h2o.mm1<-h2o.randomForest(x=1:(length(s1var)-1),y=length(s1var),data=dev1samp.hex,validation=val1samp.hex,balance.classes=T,stat.type="GINI")
#h2o.mm1<-h2o.randomForest(x=1:(length(s1var)-1),y=length(s1var),data=dev1samp.hex,validation=val1samp.hex,balance.classes=T,stat.type="GINI")

h2o.mm1b<-h2o.deeplearning(x=1:(length(s1var)-1),y=length(s1var),data=dev1samp.hex,validation=val1samp.hex,hidden=c(10,10,10,10),activation='Tanh',balance_classes=T,train_samples_per_iteration=-1,epochs=100)
h2o.mm2b<-h2o.deeplearning(x=1:(length(s2var)-1),y=length(s2var),data=dev2samp.hex,validation=val2samp.hex,hidden=c(10,10,10,10),activation='Tanh',balance_classes=T,train_samples_per_iteration=-1,epochs=100)
h2o.mm3b<-h2o.deeplearning(x=1:(length(s3var)-1),y=length(s3var),data=dev3samp.hex,validation=val3samp.hex,hidden=c(10,10,10,10),activation='Tanh',balance_classes=T,train_samples_per_iteration=-1,epochs=100)
h2o.mm4b<-h2o.deeplearning(x=1:(length(s4var)-1),y=length(s4var),data=dev4samp.hex,validation=val4samp.hex,hidden=c(10,10,10,10),activation='Tanh',balance_classes=T,train_samples_per_iteration=-1,epochs=100)

mod1va<-0*c(1:nrow(valsamp))
mod1va[val1s]<-as.data.frame(h2o.predict(h2o.mm1a,newdata=val1samp.hex))[,3]
mod1va[val2s]<-as.data.frame(h2o.predict(h2o.mm2a,newdata=val2samp.hex))[,3]
mod1va[val3s]<-as.data.frame(h2o.predict(h2o.mm3a,newdata=val3samp.hex))[,3]
mod1va[val4s]<-as.data.frame(h2o.predict(h2o.mm4a,newdata=val4samp.hex))[,3]

mod1vb<-0*c(1:nrow(valsamp))
mod1vb[val1s]<-as.data.frame(h2o.predict(h2o.mm1b,newdata=val1samp.hex))[,3]
mod1vb[val2s]<-as.data.frame(h2o.predict(h2o.mm2b,newdata=val2samp.hex))[,3]
mod1vb[val3s]<-as.data.frame(h2o.predict(h2o.mm3b,newdata=val3samp.hex))[,3]
mod1vb[val4s]<-as.data.frame(h2o.predict(h2o.mm4b,newdata=val4samp.hex))[,3]

Wl<-seq(0,1,0.05)
amsmax<-c()
threshmax<-c()
sprintf('--combining ML models--')
sprintf('Weight    AMS-max    Best Threshold')
########
for (i in c(1:length(Wl))) {
	mod1v<-Wl[i]*mod1va+(1-Wl[i])*mod1vb
	amsmat<-amseval(mod1v,valsamp$resp,valwgt,(1+nrow(devsamp)/nrow(valsamp)))
	amsmax[i]<-max(amsmat[,2],na.rm=T)
	threshmax[i]<-amsmat[which(amsmat[,2]==amsmax[i]),1]
	sprintf('%5.4f    %5.4f    %5.4f',Wl[i],amsmax[i],threshmax[i])
}

amsmax_m<-max(amsmax)
threshmax_m<-threshmax[which(amsmax==amsmax_m)]
wtmax<-Wl[which(amsmax==amsmax_m)]

save(list=c("h2o.mm1a","h2o.mm2a","h2o.mm3a","h2o.mm4a","h2o.mm1b","h2o.mm2b","h2o.mm3b","h2o.mm4b","s1var","s2var","s3var","s4var","lvars","amsmax_m","threshmax_m","wtmax"),file="mod1Params_v_5_0.RData")

h2o.saveModel(h2o.mm1a, dir="v_5_0", name="h2o.mm1a", save_cv = FALSE, force=FALSE)
h2o.saveModel(h2o.mm2a, dir="v_5_0", name="h2o.mm2a", save_cv = FALSE, force=FALSE)
h2o.saveModel(h2o.mm3a, dir="v_5_0", name="h2o.mm3a", save_cv = FALSE, force=FALSE)
h2o.saveModel(h2o.mm4a, dir="v_5_0", name="h2o.mm4a", save_cv = FALSE, force=FALSE)

h2o.saveModel(h2o.mm1b, dir="v_5_0", name="h2o.mm1b", save_cv = FALSE, force=FALSE)
h2o.saveModel(h2o.mm2b, dir="v_5_0", name="h2o.mm2b", save_cv = FALSE, force=FALSE)
h2o.saveModel(h2o.mm3b, dir="v_5_0", name="h2o.mm3b", save_cv = FALSE, force=FALSE)
h2o.saveModel(h2o.mm4b, dir="v_5_0", name="h2o.mm4b", save_cv = FALSE, force=FALSE)
