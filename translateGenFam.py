#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

translateGenFam.py

Translate genFamClust index names to true gene names


@author: mfortz
"""


import sys


# GenFam.fasta file
GenFam_fasta = open(sys.argv[1], 'r').readlines()

# input file to genFamClust
sequence_fa = open(sys.argv[2], 'r').readlines()

# GenFam.syc  or GenFam.nnc file
scoreFile = open(sys.argv[3],'r').readlines()

# Output file
translatedScoreFile = open(sys.argv[4],'w')



geneName = {}

numGenes = len(sequence_fa)//2

for i in range(numGenes):
    indexName = GenFam_fasta[2*i]
    indexName = indexName.strip()[1:]
    
    trueName = sequence_fa[2*i]
    trueName = trueName.strip()[1:]
    
    geneName[indexName] = trueName
    


for l in scoreFile:
    x = l.split()
    g1 = geneName[x[0]]
    g2 = geneName[x[1]]
    nc = x[2]
    sc = x[3]
    
    translatedScoreFile.writelines([g1,'\t',g2,'\t',nc,'\t',sc,'\n'])
