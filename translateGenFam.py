#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

translateGenFam.py

Translate genFamClust index names to true gene names

@author: mfortz
"""

import sys

<<<<<<< HEAD


# hsf input file to genFamClust
syntenyFile = open(sys.argv[1], 'r').readlines()

# GenFam.syc  or GenFam.nnc file
scoreFile = open(sys.argv[2],'r').readlines()

# Output file
translatedScoreFile = open(sys.argv[3],'w')


# Derive true gene names from syntenyFile
#
geneName = {}


# let x.y.z count species.ctg.gene numbers
x = 0
y = 0
z = 0

# initialize names using first line
firstline = syntenyFile[0].split()
oldSpecies = firstline[0]
oldCtg = firstline[1]
geneName['0.0.0'] = firstline[2] 


doneSpecies = {}
doneSpecies[oldSpecies] = ''

# go through synteny file
for l in syntenyFile[1:]:
    g = l.split()
    trueName = g[2] 
    species = g[0]
    ctg = g[1]
    
    # update counters
    if (species != oldSpecies):
        # record last entry in case species name shows up again
        doneSpecies[oldSpecies] = {'x':x, 'y':y}
        oldSpecies = species
        oldCtg = ctg
        
        try: 
            prev = doneSpecies[species]
            x = prev['x']
            y = prev['y'] + 1
            z = 0
        except KeyError:
            x += 1
            y = 0
            z = 0
        
    elif (ctg != oldCtg):
        y += 1
        z = 0
        
        oldCtg = ctg
    else:     
        z += 1
=======
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
>>>>>>> bcc42b85afcbbac3bc8f283eca437a8354e25873
    
    # record translated name
    indexName = str(x)+'.'+str(y)+'.'+str(z)
    geneName[indexName] = trueName 
    
    print(indexName, oldSpecies, oldCtg)
    
for l in scoreFile:
    if l[0]!='#':
        x = l.split()
        g1 = geneName[x[0]]
        g2 = geneName[x[1]]
        nc = x[2]
    translatedScoreFile.writelines([g1,'\t',g2,'\t',nc,'\n'])
#        sc = x[3]
#    translatedScoreFile.writelines([g1,'\t',g2,'\t',nc,'\t',sc,'\n'])

