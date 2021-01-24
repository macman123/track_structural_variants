@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @
@ ~~ ANALYSIS OF STRUCTURAL VARIATIONS AROUND A GENE ~~ @
@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

The scripts available here can be used to track structural variations around a gene of interest.
A detailed explanation of how the method works with a practical example can be found in:

	Role of the mobilome in the global dissemination of the carbapenem resistance gene blaNDM
	Mislav Acman, Ruobing Wang, Lucy van Dorp, Liam P. Shaw, Qi Wang, Nina Luhmann, Yuyao Yin, 
	Shijun Sun, Hongbin Chen, Hui Wang, Francois Balloux
	bioRxiv 2021.01.14.426698; doi: https://doi.org/10.1101/2021.01.14.426698

If you are using this code please cite the paper.


###########################
## Required Dependencies ##
###########################

Python3 with the following packages/modules: 
	os, pandas, and joblib (with Parallel and delayed)

BLAST tools (2.6.0 or above) (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	makeblastdb and blastn

R (v4.0.2+) with the following libraries:
	data.table, igraph, seqinr, data.tree, ape, ggtree, ggplot2, scales, cowplot, svglite

############################
## Required Preprocessing ##
############################

!! Flip all of the sequences/contigs you wish to analyse such that the gene of interest is always facing the same direction !!
!! If you wish to inspect the structural variations upstream of the gene then flip the gene sequence and all of the contigs 
   to face the opposite direction !!


#############################
## How to use the scripts? ##
#############################

Run the scripts within this directory in the following order:
1) BLAST_allVSall.sh
2) blast2net.py
3) track_variants.R

These are scripts, so you will have to change parts of the code to make them work. Usually these are file paths and variables.
Each script contains coments which are here to help one understand the code.
There is also a test dataset which you are free to run to familiarise yourself with the scripts or to troubleshoot. 



