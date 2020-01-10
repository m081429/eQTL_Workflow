#!/bin/bash
#
# This script will install all anciliary tools for eqtl
# usage: install_tools.sh path_to_eqtl_workflow_install_without_a_slash_at_the_end
#
#DETAILED DESCRIPTION TO DOWNLOAD THE TOOLS

if [ "$#" -eq 2 ]
then


#extracting download path
PLINK=`cat $2 |grep PLINK=|cut -f2 -d '='`
ANNOVAR=`cat $2 |grep ANNOVAR=|cut -f2 -d '='`
MASS=`cat $2 |grep MASS=|cut -f2 -d '='`
SFS=`cat $2 |grep SFS=|cut -f2 -d '='`
MATRIX=`cat $2 |grep MATRIX=|cut -f2 -d '='`

testp=`$EQTL/EXTERNALTOOLS/PLINK/$PLINK_BASE/plink -h| grep -c 'Purcell'`
export EQTL=$1
if [ $testp -eq 0 ]
then
	#Create main tools directory
	mkdir -p $EQTL/EXTERNALTOOLS
	#Change directory
	cd  $EQTL/EXTERNALTOOLS

	#PLINK 
	#Create main plink directory
	mkdir $EQTL/EXTERNALTOOLS/PLINK
	#Change to plink main directory
	cd $EQTL/EXTERNALTOOLS/PLINK
	#Download the plink package
	wget $PLINK
	#Uncompresszip files
	PLINK_BASE=`basename $PLINK`
	unzip $PLINK_BASE
	rm $PLINK_BASE
	PLINK_BASE=`echo $PLINK_BASE|sed 's/.zip//g'`;
	echo $PLINK_BASE
	testp=`$EQTL/EXTERNALTOOLS/PLINK/$PLINK_BASE/plink -h| grep -c 'Purcell'`
	if [ $testp -eq 0 ]
	then
		echo "Error during plink installation"
		exit 1
	fi
else
	echo "Skipping plink installation as exists.If you want to reinstall , delete the PLINK directory and reinstall\n";	
fi


#ANNOVAR
#Create main ANNOVAR directory	
mkdir -p $EQTL/EXTERNALTOOLS/ANNOVAR

testp=`$perl $EQTL/EXTERNALTOOLS/ANNOVAR/annovar/annotate_variation.pl | grep -c 'Version'`
if [ $testp -eq 0 ]
then
	#Change to ANNOVAR main directory
	cd $EQTL/EXTERNALTOOLS/ANNOVAR

	#Download the ANNOVAR packge (link in the email)
	wget $ANNOVAR
	ANNOVAR_BASE=`basename $ANNOVAR`
	#Uncompress zip files
	tar -zxvf $ANNOVAR_BASE
	rm $ANNOVAR_BASE
	cd annovar
	echo `pwd`
	perl=`which perl`
	if [ ! -f $perl ]
		then
		echo "perl path $perl not found"
		exit 1
	fi
	testp=`$perl $EQTL/EXTERNALTOOLS/ANNOVAR/annovar/annotate_variation.pl | grep -c 'Version'`
	if [ $testp -eq 0 ]
	then
		echo "Error during ANNOVAR installation"
		exit 1
	fi
else
	echo "Skipping ANNOVAR installation as exists.If you want to reinstall, delete the ANNOVAR directory and reinstall\n";	
fi
	
#RPACKAGES
#Create main Rpackages directory
mkdir $EQTL/EXTERNALTOOLS/Rpackages
#Change to Rpackages main directory
cd $EQTL/EXTERNALTOOLS/Rpackages
wget $MATRIX
MATRIX=`basename $MATRIX`
wget $MASS
MASS=`basename $MASS`
wget $SFS
SFS=`basename $SFS`
#Create directory .rlib.
mkdir rlib
#Install R package in that rlib directory
R=`which R`
if [ ! -f $R ];
then
   echo "R is not installed"
   exit 1
fi
$R CMD INSTALL -l ./rlib/ $MATRIX
$R CMD INSTALL -l ./rlib/ $MASS
$R CMD INSTALL -l ./rlib/ $SFS
#test whether R packages installed or not
RSCRIPT=`which Rscript` 
$RSCRIPT $EQTL/bin/check_Rpackages.R ./rlib/

#<<'end_long_comment'
#end_long_comment
else

echo "usage: install_tools.sh <path_to_ezimputer_install_without_trailing_slash> <CONFIG FILE>"
fi
