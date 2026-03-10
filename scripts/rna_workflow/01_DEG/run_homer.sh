#!/bin/bash
 
# This script requires an installation of HOMER with added path to a bash profile for loading.

geneFile=$1
outDir=$3
bg_genes=$2

mkdir -p ${outDir}




sbatch -t 72:00:00 --mem=4G --wrap="bash -c 'export PATH=/users/s/e/seyoun/tools/homer/5.1/bin:$PATH; echo $PATH;findMotifs.pl ${geneFile} human ${outDir} -bg ${bg_genes}'"

