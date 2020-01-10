#!/bin/bash
# If you installed all the tools using our own scripts, this script will produce a
# prototype of a config file as long as you have the right system tools
# installed. perl, python, java, and sun grid engine (SGE) or portable Batch System PBS.
#
# To use this script, 
# bash make_tool_info.pl full_path_for_ezimputer_without_a_slash_at_the_end
#
if [ "$#" -eq 1 ]
then
EQTL=$1
PLINK=`find ${EQTL}/EXTERNALTOOLS/ -name 'plink'`
echo "path_plink=${PLINK}"
PERL=`which perl`
echo "path_perl=${PERL}"
SH=`which sh`
echo "path_sh=${SH}"
QSUB=`which qsub`
echo "path_qsub=$QSUB"
RSCRIPT=`which Rscript`
echo "path_rscript=$RSCRIPT"
echo "path_rlib=${EQTL}/EXTERNALTOOLS/Rpackages/rlib"
echo "path_annovar=${EQTL}/EXTERNALTOOLS/ANNOVAR/annovar"

else
    echo "usage:  make_tool_info.sh path_to_eQTL_install_without_a_slash_at_the_end"
fi
