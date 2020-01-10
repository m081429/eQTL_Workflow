args <- commandArgs(TRUE)
data<-read.table(args[1],sep="\t",head=T)
s1<-seq(1:(ncol(data)-1))
s2<-sample(s1, size=length(s1), replace = FALSE, prob = NULL)
s3<-s2+1
s3<-c(1,s3)
data1<-data[,s3]
write.table(data1,args[2],quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
