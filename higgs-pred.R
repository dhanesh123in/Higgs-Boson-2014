# install xgboost package, see R-package in root folder
require(xgboost)
require(methods)

modelfile <- "higgs.model"
outfile <- "submission_xgboost-4.csv"
dtest <- read.csv("data/test.csv", header=TRUE)
data <- as.matrix(dtest[2:31])
idx <- dtest[[1]]

lvars<-c('DER_mass_MMC','DER_mass_transverse_met_lep','DER_mass_vis','DER_pt_h','DER_deltar_tau_lep','DER_pt_tot','DER_sum_pt','DER_pt_ratio_lep_tau','PRI_tau_pt','PRI_lep_pt','PRI_met','PRI_met_sumet','PRI_jet_leading_pt','PRI_jet_all_pt')

for (xn in lvars) {
	data[which(data[,xn]!=-999),xn]<-log(data[which(data[,xn]!=-999),xn]+1)
}

xgmat <- xgb.DMatrix(data, missing = -999.0)
bst <- xgb.load(modelfile=modelfile)
ypred <- predict(bst, xgmat)


load(file="xgboost-4.RData")
threshold=amsmat[which(amsmat[,2]==max(amsmat[,2])),1]
rorder <- rank(ypred, ties.method="first")
#threshold <- 0.15
# to be completed
#ntop <- length(rorder) - as.integer(threshold*length(rorder))
#plabel <- ifelse(rorder > ntop, "s", "b")

plabel <- ifelse(ypred>threshold, "s", "b")
outdata <- list("EventId" = idx,
                "RankOrder" = rorder,
                "Class" = plabel)
write.csv(outdata, file = outfile, quote=FALSE, row.names=FALSE)
