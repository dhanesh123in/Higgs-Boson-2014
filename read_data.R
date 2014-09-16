
traindat<-read.csv2("training.csv",header=T,sep=',')

set.seed(111)
ind <- sample(2, nrow(traindat), replace = TRUE, prob=c(0.8, 0.2))

devsamp<-traindat[ind==1,]
valsamp<-traindat[ind==2,]

save(list=c("ind","devsamp","valsamp"),file="HB_dev_val_samps.RData")

