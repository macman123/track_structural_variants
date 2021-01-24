#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 13:44:05 2019

@author: macman
"""

import os
import pandas as pd
from joblib import Parallel, delayed

# It might take a while for this script to run so be patient.

# This function selects blast hits that cover the entire gene of interest
# and outputs a network in pajek format. The network is basically a table
# where each row represents a pair of sequences/contigs. There are three
# columns: ID of 1st sequence; ID of 2nd sequence; BLAST hit length.
# In subsequent analysis this network is filtered using a splitting
# threshold and the structural variations are uncovered.

def filter_blast(file, GENE_pos):
    
    filename = os.fsdecode(file)
     
    if(os.path.getsize("results/blastout/" + filename) == 0):
        return
    blastout = pd.read_table("results/blastout/" + filename, header=None)

    seqID_1 = blastout.iloc[0][0] 
    seqID_2 = blastout.iloc[0][1]
    GENE_pos_seq_1 = GENE_pos.loc[GENE_pos[1] == seqID_1]
    GENE_pos_seq_2 = GENE_pos.loc[GENE_pos[1] == seqID_2]
     
 
    # NOTE: All sequences/contigs should be fliped so the gene of interest is
    # always facing the same direction. Therfore, there is no need to flip the 
    # blast results. If the qend is bigger than qstart for example it means 
    # that the insertion happened in the other direction and these sequences 
    # need to be treated  differently 
    blastout_filtered = blastout[(blastout[6] <= GENE_pos_seq_1.iloc[0][8]) & \
                                 (blastout[6] <= GENE_pos_seq_1.iloc[0][9]) & \
                                 (blastout[7] >= GENE_pos_seq_1.iloc[0][8]) & \
                                 (blastout[7] >= GENE_pos_seq_1.iloc[0][9]) & \
                                 (blastout[8] <= GENE_pos_seq_2.iloc[0][8]) & \
                                 (blastout[8] <= GENE_pos_seq_2.iloc[0][9]) & \
                                 (blastout[9] >= GENE_pos_seq_2.iloc[0][8]) & \
                                 (blastout[9] >= GENE_pos_seq_2.iloc[0][9])]

    #save filtered blast
    blastout_filtered.to_csv("results/blastout/" + file+".filtered.tsv",sep="\t", header=False, index=False)
    
    #writing a file
    if(len(blastout_filtered.index) == 0):
        return
    max_alignment_length=blastout_filtered[3]-blastout_filtered[5]
    return(seqID_1 + " " + seqID_2 + " " + str(max_alignment_length.iloc[0])+ "\n")


GENE_pos =  pd.read_table("results/GENE_positions.txt", header=None)

# currently requires 12 cores. Change n_jobs if you want more or less cores to be used
net = Parallel(n_jobs=12)(delayed(filter_blast)(file, GENE_pos) \
               for file in os.listdir("results/blastout"))
net = list(filter(None, net))
with open("results/blast_network.net","w+") as f:
    for item in net:
        f.write(item)


