rm(list=ls())
require(earth)
load(file="mod1Params_v_5_0.RData")

tstdat<-read.csv2("test.csv",header=T,sep=',')

catvars<-c('EventId','PRI_jet_num')
numvars<-setdiff(names(tstdat),catvars)

tstdat[,'PRI_jet_num']<-as.factor(tstdat[,'PRI_jet_num'])
tstdat[,'EventId']<-as.character(tstdat[,'EventId'])

for (i in numvars) {

	tstdat[,i]<-as.numeric(as.character(tstdat[,i]))
	
}

tstevents<-tstdat['EventId']
tstdat<-tstdat[-1]

tst1s<-which(tstdat$PRI_jet_num==0)
tst2s<-which(tstdat$PRI_jet_num==1)
tst3s<-which(tstdat$PRI_jet_num==2)
tst4s<-which(tstdat$PRI_jet_num==3)


tmiind<-which(tstdat[,'DER_mass_MMC']==-999)
tstdat[tmiind,'DER_mass_MMC']<-0

##Log-variables 
lvars<-c('DER_mass_MMC','DER_mass_transverse_met_lep','DER_mass_vis','DER_pt_h','DER_deltar_tau_lep','DER_pt_tot','DER_sum_pt','DER_pt_ratio_lep_tau','PRI_tau_pt','PRI_lep_pt','PRI_met','PRI_met_sumet','PRI_jet_leading_pt','PRI_jet_all_pt')

for (xn in intersect(lvars,s1var)) {
	tstdat[tst1s,xn]<-log(tstdat[tst1s,xn]+1)
}

for (xn in intersect(lvars,s2var)) {
	tstdat[tst2s,xn]<-log(tstdat[tst2s,xn]+1)
}

for (xn in intersect(lvars,s3var)) {
	tstdat[tst3s,xn]<-log(tstdat[tst3s,xn]+1)
}

for (xn in intersect(lvars,s4var)) {
	tstdat[tst4s,xn]<-log(tstdat[tst4s,xn]+1)
}

require(h2o)
localH2O<-h2o.init()

tst1samp.hex<-as.h2o(localH2O,tstdat[tst1s,setdiff(s1var,'resp')])
tst2samp.hex<-as.h2o(localH2O,tstdat[tst2s,setdiff(s2var,'resp')])
tst3samp.hex<-as.h2o(localH2O,tstdat[tst3s,setdiff(s3var,'resp')])
tst4samp.hex<-as.h2o(localH2O,tstdat[tst4s,setdiff(s4var,'resp')])


mod1ta<-0*c(1:nrow(tstdat))
mod1ta[tst1s]<-as.data.frame(h2o.predict(h2o.mm1a,newdata=tst1samp.hex))[,3]
mod1ta[tst2s]<-as.data.frame(h2o.predict(h2o.mm2a,newdata=tst2samp.hex))[,3]
mod1ta[tst3s]<-as.data.frame(h2o.predict(h2o.mm3a,newdata=tst3samp.hex))[,3]
mod1ta[tst4s]<-as.data.frame(h2o.predict(h2o.mm4a,newdata=tst4samp.hex))[,3]

mod1tb<-0*c(1:nrow(tstdat))
mod1tb[tst1s]<-as.data.frame(h2o.predict(h2o.mm1b,newdata=tst1samp.hex))[,3]
mod1tb[tst2s]<-as.data.frame(h2o.predict(h2o.mm2b,newdata=tst2samp.hex))[,3]
mod1tb[tst3s]<-as.data.frame(h2o.predict(h2o.mm3b,newdata=tst3samp.hex))[,3]
mod1tb[tst4s]<-as.data.frame(h2o.predict(h2o.mm4b,newdata=tst4samp.hex))[,3]

mod1t<-wtmax*mod1ta+(1-wtmax)*mod1tb

pred1t<-rep('b',nrow(tstevents))

pred1t[which(mod1t>threshmax_m)]='s'

pred1t<-as.factor(as.character(pred1t))

submt<-as.data.frame(cbind(tstevents,order(mod1t),pred1t))
names(submt)<-c("EventId","RankOrder","Class")
write.csv(submt,"submission_h2o_v5_0.csv",row.names=FALSE,quote=F)

