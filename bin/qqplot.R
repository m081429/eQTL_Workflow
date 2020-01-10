args <- commandArgs(TRUE)
setwd(args[1])
infile<-args[2]
outfile<-args[3]
pval<-read.table(infile,sep="\t",head=T)
png(filename = outfile, units="in", width=11, height=8.5, res=300)
vec<-sort(rev(as.numeric(as.vector(pval[,4]))))
qqplot(x=-log10(runif(length(vec))),-log10(vec),main="QQ PLOT PVAL",xlab="theory quantiles",ylab="data quantiles", plot.it=TRUE, datax=FALSE,distribution=dunif)
#qqline(vec,distribution=qunif(vec),col="red")
qqline(-log10(vec),distribution=function(p) qchisq(p),col="red")
dev.off()
