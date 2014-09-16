rm(list=ls())
# install xgboost package, see R-package in root folder
require(xgboost)
require(methods)

testsize <- 550000

dtrain <- read.csv("data/training.csv", header=TRUE)
dtrain[33] <- dtrain[33] == "s"
label <- as.numeric(dtrain[[33]])
data <- as.matrix(dtrain[2:31])

lvars<-c('DER_mass_MMC','DER_mass_transverse_met_lep','DER_mass_vis','DER_pt_h','DER_deltar_tau_lep','DER_pt_tot','DER_sum_pt','DER_pt_ratio_lep_tau','PRI_tau_pt','PRI_lep_pt','PRI_met','PRI_met_sumet','PRI_jet_leading_pt','PRI_jet_all_pt')

for (xn in lvars) {
	data[which(data[,xn]!=-999),xn]<-log(data[which(data[,xn]!=-999),xn]+1)
}

set.seed(123);
tset<-which(runif(nrow(data))>0.2)
vset<-setdiff(c(1:nrow(data)),tset)


weightt <- as.numeric(dtrain[[32]])[tset] * testsize / length(label[tset])

sumwpos <- sum(weightt * (label[tset]==1.0))
sumwneg <- sum(weightt * (label[tset]==0.0))
print(paste("weight statistics: wpos=", sumwpos, "wneg=", sumwneg, "ratio=", sumwneg / sumwpos))

xgmat <- xgb.DMatrix(data[tset,], label = label[tset], weight = weightt, missing = -999.0)
param <- list("objective" = "binary:logitraw",
              "scale_pos_weight" = sumwneg / sumwpos,
              "bst:eta" = 0.1,
              "bst:max_depth" = 6,
              "eval_metric" = "auc",
              "eval_metric" = "ams@0.15",
              "silent" = 1,
              "nthread" = 16)
watchlist <- list("train" = xgmat)
nround = 200
print ("loading data end, start to boost trees")
bst = xgb.train(param, xgmat, nround, watchlist );
# save out model
xgb.save(bst, "higgs.model")
print ('finish training')


#Validation set
weightv <- as.numeric(dtrain[[32]])[vset] * testsize / length(label[vset])
sumwposv <- sum(weightv * (label[vset]==1.0))
sumwnegv <- sum(weightv * (label[vset]==0.0))
print(paste("weight statistics: wpos=", sumwposv, "wneg=", sumwnegv, "ratio=", sumwnegv / sumwposv))

xgmatv <- xgb.DMatrix(data[vset,], label = label[vset], weight = weightv, missing = -999.0)
ypredv <- predict(bst, xgmatv)

amseval <- function (pred,actual,wgt,wfactor) {
	thresh<-sort(pred)
	
	tmat<-as.data.frame(cbind(pred,wgt*(actual),wgt*(1-actual)))
		
	tmato<-tmat[order(tmat[,1]),]
	tmato[,2]<-cumsum(tmato[,2])*wfactor
	tmato[,3]<-cumsum(tmato[,3])*wfactor
	
	sm<-tmato[1,2]
	bm<-tmato[1,3]
	
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


amsmat<-amseval(ypredv,dtrain[vset,33],as.numeric(dtrain[[32]])[vset],length(label)/length(vset))


save(list=c("tset","vset","amsmat"),file="xgboost-4.RData")