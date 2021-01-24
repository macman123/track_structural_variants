#! /bin/bash

## CHANGE THESE PARAMETERS:
# directory containing sequences to analyse (don't forget the *):
contigs=test_data/*
# a gene of interest:
gene=test_gene.fna
# BLAST bin path
blast=~/bin/ncbi-blast-2.10.1+/bin


# The script is blasting pairs of sequences. Here all possible pairs of sequences are made:
echo "--> MAKE COMBINATIONS"
f=($contigs.fna)
mkdir -p results 
> results/BLAST_comb.txt
for ((i = 0; i < ${#f[@]}; i++)); do  
  for ((j = i + 1; j < ${#f[@]}; j++)); do 
    printf "${f[i]} ${f[j]}\n" >> results/BLAST_comb.txt; 
  done;
done


# Making blast DBs for pairwise comparisons:
echo "--> MAKE BLAST DBs"
for f in $contigs; do
  ${blast}/makeblastdb -in $f -dbtype nucl
done


# Finding position of the gene:
echo "--> Finding gene positions"
for f in $contigs.fna; do
  ${blast}/blastn -db $f -query $gene -max_hsps 1 -outfmt 6 -out $f.blastout
done
>results/GENE_positions.txt
cat ${contigs}*blastout >> results/GENE_positions.txt
rm -rf ${contigs}/*blastout


# BLASTING ALL VS ALL
mkdir -p results/blastout
echo "--> BLAST ALL VS ALL"
calc(){ awk "BEGIN { print $* }"; }
counter=1
lines=`wc -l < results/BLAST_comb.txt`
while read comb; do
  f1=`echo $comb | cut -d" " -f1`
  f2=`echo $comb | cut -d" " -f2`
  perc=`calc $counter/$lines*100`
  echo -ne "\r${perc}% completed"
  $blast/blastn -db $f1 -query $f2 -task dc-megablast -outfmt 6 -out results/blastout/$counter.blastout -num_threads 12
  counter=$((counter + 1))
done < results/BLAST_comb.txt

