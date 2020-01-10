# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R (Revolution R and Matlab are faster on PC).
#
# Set working directory:
# setwd('');
#
args <- commandArgs(TRUE);
library(methods);
library(MatrixEQTL,lib.loc=args[1]);
#source("/data5/bsi/bioinf_ext1/s113625.eQTL/EQTL_WORKFLOW/bin/tools/Matrix_eQTL_engine.R")
#source("/data2/bsi/tertiary/Wang_Chen_m092469/eQTL_mrna_seq_florida_alz/Res/AD_workflow/tmp/run1_Chen/processing/MatrixEQTL/R/Matrix_eQTL_engine.R");
#setwd(args[1])
#source(args[1]);
setwd(args[2]);
tempdir<-function(){args[2]}
SNP_file_name ="SNP.txt"; 
snps_location_file_name="snploc.txt";
# Gene expression file name
#expression_file_name = "/data4/bsi/bioinf_ext1/s112132.BreastCancer_eQTL/GE.txt";
expression_file_name =args[3];
gene_location_file_name = 'geneloc.txt';

# Covariates file name
# Set to character() for no covariates
#covariates_file_name = 'gCovariates.txt';
covariates_file_name =character();
if(file.exists("covariates.txt"))
{
	covariates_file_name ='covariates.txt';
	print("covariates exists");
}
# Output file name
#output_file_name = "result_main_eqtl.txt";
output_file_name = args[4];
output_file_name.cis=args[7];
# Only associations significant at this level will be saved
#pvOutputThreshold_cis = 1e-2;
#pvOutputThreshold_tra = 1e-2;
pvOutputThreshold_cis = as.numeric(as.character(args[5]));
pvOutputThreshold_tra = as.numeric(as.character(args[6]));
#print("test")
#print(pvOutputThreshold_cis)
#print("test1")
#print(pvOutputThreshold_tra)
# Error covariance matrix
# Set to character() for identity.
errorCovariance = character();
# errorCovariance = read.table('gerrorCovariance.txt');

#cisDist = 1e5;
#cisDist=100000
cisDist = as.numeric(args[8]);
#print(cisDist)
#exit
useModel = modelLINEAR; # modelANOVA or modelLINEAR
#useModel =args[7];
## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = '\t'; # the TAB character
snps$fileOmitCharacters = 'NA'; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 2000; # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = '\t'; # the TAB character
gene$fileOmitCharacters = 'NA'; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 2000; # read file in pieces of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = '\t'; # the TAB character
cvrt$fileOmitCharacters = 'NA'; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = 2000; # read file in one piece
if(length(covariates_file_name)>0) {
	cvrt$LoadFile(covariates_file_name);
}
#cvrt=character()
## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
errorCovariance = numeric();
Matrix_eQTL_main(snps=snps,gene=gene, cvrt=cvrt, output_file_name=output_file_name,output_file_name.cis=output_file_name.cis,pvOutputThreshold.cis=pvOutputThreshold_cis,pvOutputThreshold=pvOutputThreshold_tra,useModel=useModel,errorCovariance=errorCovariance,verbose=TRUE,snpspos=snpspos,genepos=genepos,cisDist=cisDist,pvalue.hist=FALSE);
#data<-read.table(output_file_name,sep="\t",head=T)
#s5<-table(data[,6])
#write.table(s5,"out.txt",sep="\t",col.names=FALSE,row.names=FALSE)
#data<-data[with(data, order(p.value)), ]
#write.table(data,"/data4/bsi/bioinf_ext1/s112132.BreastCancer_eQTL/PVAL.txt",sep="\t",col.names=TRUE,row.names=FALSE)
#data<-data[with(data, order(FDR)), ]
#write.table(data,"/data4/bsi/bioinf_ext1/s112132.BreastCancer_eQTL/FDR.txt",sep="\t",col.names=TRUE,row.names=FALSE)
q()
n
