args <- commandArgs(TRUE)
library(MASS,lib.loc=args[1])
library(sfsmisc,lib.loc=args[1])
snp<-read.table(args[2],sep="\t",skip=1)
rownames(snp)<-snp[,1]
gene<-read.table(args[3],sep="\t",head=T)
data2<-read.table(args[4],sep="\t",head=T)
#colnames(data2)<-c("rsid","geneid")
rownames(gene)<-as.character(gene[,1])
k1<-c()
k1<-c("gene","snp","Intercept_Estimate","snp_Estimate","Intercept_Std.Error","snp_Std.Error","Intercept_t value","snp_t value","ROBUST_F_TEST_RLMFIT")
z<-1
while(z<=nrow(data2))
{
        #print(z)
        kk<-c()
        flm<-rlm(as.numeric(as.character(snp[as.character(data2[z,1]),2:ncol(snp)])) ~ as.numeric(as.character(gene[as.character(data2[z,2]),2:ncol(gene)])), family=gaussian(link="log"))
        kk<-c(as.character(data2[z,2]),as.character(data2[z,1]),as.vector(coefficients(summary(flm))),f.robftest(flm)$p.value)
        k1<-rbind(k1,kk)
        z<-z+1
}
write.table(k1,args[5],sep="\t",quote=F,row.names=F,col.names=F)

