#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
makeSequence.py

Formats gene sequence for use as genfamclust input.

Arguments:

    Input:  
        1. File with two columns: species name, and path to fasta file 
                of the corresponding species.
           Fasta files must be in standard fasta format, but with one 
                line for each nucleotide sequence.
            
        2. ALL_GENE_file, containing position of genes on their contig

    Output: 
        3. Sequence File (*.fasta) :
            List with two columns, tab delimited. First columnn for gene name, 
            second column for full sequence.
        4. Synteny File (*.hsf):
            List with four columns, tab delimited. Format: 
            [species name, chromosome name, gene name, position on chromosome]
            
Created on Fri Jun  2 13:10:11 2017

@author: mfortz
"""

import sys, re

#NOTE
#currently reading fastaFile for one species. 
#need to read all species

fastaLocation = open(sys.argv[1],'r').readlines()   
allGeneFile = open(sys.argv[2],'r').readlines()

sequenceFile = open(sys.argv[3],'w')
syntenyFile = open(sys.argv[4], 'w')

#temp storage
fullSeqDict = {}
syntenyDict = {}

def makeSeqDict(fastaFile):
    #read fastaFile
    for l in fastaFile:
        #differentiate between lines with ctg name and lines with nucleotides
        if l[0] == '>':
            ctg = l.split()[0][1:]      #takes new ctg name
            fullSeqDict[ctg] = ''
        elif len(l)>2:
            #add nucleotide line to the sequence belonging to last ctg
            #ignore '\n' or '\t at the end of sequence
            fullSeqDict[ctg] = l.rstrip()                                           

def makeSynDict(currentSpecies):
    #read ALL_GENE_file to get preliminary synteny information     
    for l in allGeneFile:
        if l[0]!='#':
            g = l.split()        
            #take genes of currentSpecies only
            sp = g[0]
            
            if sp != currentSpecies:
                continue
            else:        
                #info needed for synteny file output
                ctg = g[1]            
                name = g[3]
                orient = g[4]
                st = int(g[5])            
            
                #split exons position into a list of (start,end) tuples
                #some lines have no data for exons position so need to check 
                if len(g) > 8:
                    pos = re.split('\D',g[8])
                    segments = []
                    for i in range(len(pos)//2):
                        start = int(pos[2*i])
                        end = int(pos[2*i+1])
                        segments.append((start,end))
                                    
                    #add to (or make) a list of genes for every contig
                    try:
                        syntenyDict[ctg].append([sp, ctg, name, st, segments, orient])
                    except KeyError:
                        syntenyDict[ctg] = [[sp, ctg, name, st, segments, orient]]
                                        
def rev_comp(seq):
    complement = {'a':'t','A':'T','c':'g','C':'G','g':'c','G':'C','t':'a',
                  'T':'A','n':'n','N':'N'}
    
    seqComplement = ''
    
    #add complement to start of string
    for x in seq:
        seqComplement = complement[x] + seqComplement
    
    return seqComplement
                            
for x in fastaLocation:
    if x[0]!='#':
        x = x.split()
        species = x[0]
        fastaFile = open(x[1], 'r')
        
        makeSeqDict(fastaFile)
        makeSynDict(species)
        
        fastaFile.close()       #free up memory
      
#sort genes on each ctg by start position
for c in syntenyDict:
    syntenyDict[c].sort(key=lambda line: line[3])
        
#make output files
for c in syntenyDict:
    for g in syntenyDict[c]:
        name = g[2]
        sp = g[0]
        st = str(g[3])
        orient = g[5]
        
        #find full gene sequence
        segments = g[4]
        geneSeq = ''
        for (start, end) in segments:
            geneSeq += fullSeqDict[c][start-1:end]
                
        #get the reversed complement of '-' oriented genes
        if orient == '-':
            geneSeq = rev_comp(geneSeq)
                        
        #write to output files
        sequenceFile.writelines(['>',name,'\n'])
        sequenceFile.write(geneSeq)
        sequenceFile.write('\n')
        
        syntenyFile.writelines([sp, '\t', c, '\t', name,'\t', st])
        syntenyFile.write('\n')


