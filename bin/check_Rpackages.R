args <- commandArgs(TRUE)
library(MASS,lib.loc=args[1])
library(sfsmisc,lib.loc=args[1])
library(MatrixEQTL,lib.loc=args[1])
#?Matrix_eQTL_engine
print("package installation sucessful")
