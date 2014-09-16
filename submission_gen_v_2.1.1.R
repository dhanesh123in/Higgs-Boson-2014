rm(list=ls())
require(earth)
load(file="mod1Params_v_2_1_1.RData")

tstdat<-read.csv2("test.csv",header=T,sep=',')

catvars<-c('EventId','PRI_jet_num')
numvars<-setdiff(names(tstdat),catvars)

tstdat[,'PRI_jet_num']<-as.factor(tstdat[,'PRI_jet_num'])
tstdat[,'EventId']<-as.character(tstdat[,'EventId'])

for (i in numvars) {

	tstdat[,i]<-as.numeric(as.character(tstdat[,i]))

	tstdat[which(tstdat[,i]==-999),i]<-0
	
}

tstevents<-tstdat['EventId']
tstdat<-tstdat[-1]

tst1s<-which(tstdat$PRI_jet_num==0)
tst2s<-which(tstdat$PRI_jet_num==1)
tst3s<-which(tstdat$PRI_jet_num==2)
tst4s<-which(tstdat$PRI_jet_num==3)



frm1<-as.formula(paste("~",paste(setdiff(s1var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s1var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))
frm2<-as.formula(paste("~",paste(setdiff(s2var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s2var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))
frm3<-as.formula(paste("~",paste(setdiff(s3var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s3var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))
frm4<-as.formula(paste("~",paste(setdiff(s4var,'resp'),sep="",collapse="+"),"+",paste(paste("I(",setdiff(s4var,'resp'),"^2)",sep=""),sep="",collapse="+"),sep="",collapse="+"))


mod1t<-0*c(1:nrow(tstdat))
mod1t[tst1s]<-predict(mm1,newdata=model.frame(frm1,tstdat[tst1s,setdiff(s1var,'resp')]), type = "response")
mod1t[tst2s]<-predict(mm2,newdata=model.frame(frm2,tstdat[tst2s,setdiff(s2var,'resp')]), type = "response")
mod1t[tst3s]<-predict(mm3,newdata=model.frame(frm3,tstdat[tst3s,setdiff(s3var,'resp')]), type = "response")
mod1t[tst4s]<-predict(mm4,newdata=model.frame(frm4,tstdat[tst4s,setdiff(s4var,'resp')]), type = "response")



pred1t<-rep('b',nrow(tstevents))

pred1t[which(mod1t>threshmax)]='s'

pred1t<-as.factor(as.character(pred1t))

submt<-as.data.frame(cbind(tstevents,order(mod1t),pred1t))
names(submt)<-c("EventId","RankOrder","Class")
write.csv(submt,"submission_lr_v2_1_1.csv",row.names=FALSE,quote=F)

