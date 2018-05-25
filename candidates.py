#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: mfortz
"""

from __future__ import division

import sys    
import itertools    
from deriveGeneStructure import deriveGeneStructure

all_gene_file = open(sys.argv[1],'r').readlines()
data, geneOrder = deriveGeneStructure(all_gene_file)

#output file
candidatesOut = open(sys.argv[2], 'w+')
candidatesOut.write("#F1 \t F2_gene \n")

nonEmptySpecies = {}

### Shared neighbours  
if len(sys.argv) > 2:  
    num_species_spanned_by_f1 = int(sys.argv[3])
    if num_species_spanned_by_f1 < 2:
        print("Error. Max size must be >= 2.")
        sys.exit()
    
    ### Tandem duplicates
    min_species_spanned_by_f1 = int(sys.argv[4])
    if min_species_spanned_by_f1 < 2:
        print("Error. Max size must be >= 2.")
        sys.exit()

    ### Ratio observed synteny / size of smaller cluster
    ratioIn = float(sys.argv[5])
    if ratioIn > 1 or ratioIn < 0:
        print("Error. RatioIn must be in [0,1].")
        sys.exit()
    
else :
    ### Shared neighbours    
    num_species_spanned_by_f1 = 5
    ### Tandem duplicates
    min_species_spanned_by_f1 = 5
    ### Ratio observed synteny / size of smaller cluster
    ratioIn = 0.8 

########################################
# Tandem Duplicates
########################################

for fam in data.familiesDict:
    nonEmptySpecies[fam] = 0    
    currentFamily = data.familiesDict[fam]
    
    for species in data.speciesDict:
        if len(currentFamily.familyMembersDict[species]) > 0:
            nonEmptySpecies[fam] += 1

Result = {}

for gene in data.genesDict:
    currentGene = data.genesDict[gene]
    currentFamName = currentGene.family
    
    #look for f1 candidates, ie less than max size
    if nonEmptySpecies[currentFamName] >= min_species_spanned_by_f1:
        
        #load an ordered list of genes on the current contig and species
        currentFragment = geneOrder[currentGene.species][currentGene.ctg]   
        positionOnCtg = currentFragment.index(currentGene) 
               
        #adds the neighbour to the left or to the right
        #offset = 0 if left nbr, 1 if right nbr
        def addNbr(offset):  
            
            #takes care of errors arising from the gene being at the end of its contig
            if offset == 0:
                if positionOnCtg == 0:
                    return
                else:
                    nbr = currentFragment[positionOnCtg - 1]
            elif offset == 1:
                try:
                    nbr = currentFragment[positionOnCtg + 1]
                except IndexError:
                    return
    
            #aif f1 != f2 and size_f1 <= size_f2, add to the appropriate list
            if currentFamName != nbr.family and nonEmptySpecies[currentFamName] <= nonEmptySpecies[nbr.family]:
                try:
                    Result[(currentFamName,nbr.family)].append((currentGene.gene,nbr.gene))
                except KeyError:
                    Result[(currentFamName,nbr.family)] = []
                    Result[(currentFamName,nbr.family)].append((currentGene.gene,nbr.gene))
                
        #add both left and right neighbours to Result
        addNbr(0)
        addNbr(1)
        
"""Filter result for pairs such that the ratio of syntenic order to max_size is above threshold"""
doneMap = {}

for f1,f2 in Result:
    
    """Filter for repeated entries, ie instances of (f1,f2) and (f2,f1)"""
    try:
        doneMap[(f2,f1)]
        continue
    except KeyError:
        pass
    
    observedSynteny = len(Result[(f1,f2)])
    f1size = nonEmptySpecies[f1]
    ratio = observedSynteny/f1size
    
    if ratio >= ratioIn:
        doneMap[(f1,f2)] = None
        doneMap[(f2,f1)] = None        
        candidatesOut.writelines([f1,'\t',f2,'\n'])
        
##########################################
# Shared Neighbour
##########################################

def findFamilyInXSpecies(x,q):

    addressBook = []
    onlyXSpeciesDict = {}
    allButXplusQSpeciesDict = {}
        
    #Find all size 'x' subsets of species. We want each subset to be in order.
    speciesList = list(data.speciesDict.keys())
    speciesList.sort()
    addressBook = list(itertools.combinations(speciesList,x))    
    
    #Make an entry in the dictionaries for every x-subset
    for xSet in addressBook:
        onlyXSpeciesDict[xSet] = []
        allButXplusQSpeciesDict[xSet] = []
    
    #Go through all families to add them in the appropriate dictionary
    for f in data.familiesDict:
        currentFamily = data.familiesDict[f]   #family object

        nonemptySpecies = []    #gives species in which the family appears in
        missingSpecies = []     #gives species in which the family does not appear in

        #find a set of species which the family appears in, and another set which the family does not appear in
        for species in currentFamily.familyMembersDict:   
            if len(currentFamily.familyMembersDict[species]) > 0:
                nonemptySpecies.append(species)
            else:
                missingSpecies.append(species)
                
        #in order to use nonemptySpecies as a dictionary key, we will sort it and change to tuple
        nonemptySpecies.sort()
        nonemptySpecies = tuple(s for s in nonemptySpecies)
        
        #add the family to the appropriate dictionary and key
        if len(nonemptySpecies) == x:
             onlyXSpeciesDict[nonemptySpecies].append(f)

        elif len(missingSpecies) >= x and len(missingSpecies) <= x+q :
            #if more than x species are missing, we want to add the family to all possible spots
            for xSubset in addressBook:
                
                #check every xSubset that is a subset of missingSpecies
                isSubset = True
                
                for xS in xSubset:
                    if not xS in missingSpecies:
                        isSubset = False
                        break
                if isSubset:    
                    allButXplusQSpeciesDict[xSubset].append(f)
                
    return (addressBook, onlyXSpeciesDict, allButXplusQSpeciesDict)

'''Functions to add entries to result. 
*****Make sure currentGene's family is in F1 or F2******

Input: 
    xSubset = size x subset of species on addressBook
    flag= 0 if family is in only x species (i.e. F1, F3 configuration)
          1 if family is in all but x species (i.e. F2, F3 configuration)
Output: (none) Adds entry to Result.
'''

def addResult(xSubset, flag):
    
    #load an ordered list of genes on the current contig and species
    currentFragment = geneOrder[currentGene.species][currentGene.ctg]   
    positionOnCtg = currentFragment.index(currentGene)  
        
    #add xSubset as a dictionary key if necessary
    if not xSubset in Result:
        Result[xSubset] = {}
        
    #adds the neighbour to the left or to the right
    #offset = 0 if left nbr, 1 if right nbr
    def addNbr(offset):  
                
        #takes care of errors arising from the gene being at the end of its contig
        if offset == 0:
            if positionOnCtg == 0:
                return
            else:
                nbr = currentFragment[positionOnCtg - 1]
        elif offset == 1:
            try:
                nbr = currentFragment[positionOnCtg + 1]
            except IndexError:
                return

        #adds to the appropriate list
        try:
            Result[xSubset][nbr.family][flag].append((currentGene.gene, nbr.gene))
        except KeyError:
            Result[xSubset][nbr.family] = ([],[])
            Result[xSubset][nbr.family][flag].append((currentGene.gene, nbr.gene))
    
    #add both left and right neighbours to Result
    addNbr(0)
    addNbr(1)

'''Creates Result dictionaries
ResultDict is a dictionary with keys (x,q) mapping to the appropriate Result dictionary.

Result = A double indexed dictionary. 
        Format :  Result[D][F3] = (list0,list1)
        
        D is a tuple of x families, sorted in default string order
        F3 is the name of a family that is close to F1 or F2. 
        list0 is a list of ordered pairs (gF1,gF3), where gF1 is in F1 and gF3 is in F3
        list1 is a list of ordered pairs (gF2,gF3), where gF2 is in F2 and gF3 is in F3

'''

ResultDict = {}

for x in range(2,num_species_spanned_by_f1+1):
    q =2

    Result = {}

    #Load dictionaries containing families of interest
    addressBook, onlyXSpeciesDict, allButXplusQSpeciesDict =  findFamilyInXSpecies(x,q)

    for g in data.genesDict:
        currentGene = data.genesDict[g]

        #go through addressBook to find the set of species where the currentGene's family belongs
        ##this method is kinda slow since at worst we go through (n choose x) sets for every gene. Consider revising.
        for xSubset in addressBook:
            
            if currentGene.species in xSubset:   
                if currentGene.family in onlyXSpeciesDict[xSubset]:
                    addResult(xSubset,0)
                    break

            else:
                if currentGene.family in allButXplusQSpeciesDict[xSubset]:
                    addResult(xSubset,1)
                    #family may be in multiple allButX sets, so no need to break.

    #update complete results dict
    ResultDict[(x,q)] = Result
        

'''Filters Result

Requires Result dictiopnaries and files to write in.
Outputs a txt file. Each line contains F1, F3 genes which are close to each other and F2, F3 genes which are close.

'''
    
candidateList = []

for x in range(2,num_species_spanned_by_f1+1):
    q =2
    
    #Load dictionary and txt file
    Result = ResultDict[(x,q)]
    
    for d in Result:
        for f in Result[d]:
            if Result[d][f][0] and Result[d][f][1]:  #we have a match if both lists are non-empty
   
                #go through the cartesian product of the two lists
                for g1,g3a in Result[d][f][0]:
                    F1_gene = data.genesDict[g1]
                    F1 = F1_gene.family

                    for g2,g3b in Result[d][f][1]:
                        F2_gene = data.genesDict[g2]
                        F2 = F2_gene.family
                                                
                        candidateList.append((F1,F2))

candidateSet = set(candidateList)

for f1,f2 in candidateSet:
    candidatesOut.writelines([f1, '\t', f2, '\n'])
