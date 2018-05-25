#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

deriveGeneStructure.py

Created on Wed Jul 26 19:00:24 2017

@author: mfortz
"""

def deriveGeneStructure(all_gene_file):
    
    """
    copy of structures/main.py
    find out how to load object from another script
    """
    
    # might need to switch directory to structures

    from DataSet import DataSet
    from GeneFamily import GeneFamily
        
    """Create all Gene objects"""
    data = DataSet()
    
    for l in all_gene_file:
        if l[0]!='#':
            data.addGene(l.split())
        
    """Create a list of Species"""
    def generateSpeciesDict():
        for g in data.genesDict:
            currentSpecies = data.genesDict[g].species
            if currentSpecies not in data.speciesDict:
                data.addSpecies(currentSpecies)
                
    generateSpeciesDict()
        
    """Make Family Dictionary"""
    def generateFamilyDict():
        for g in data.genesDict:
            currentGene = data.genesDict[g]
            
            if currentGene.family not in data.familiesDict:
                currentFamily = GeneFamily(currentGene.family,data.speciesDict) 
                data.addFamily(currentFamily.familyName,currentFamily)
            else:
                currentFamily = data.familiesDict[currentGene.family]
            currentFamily.addToFamily(currentGene)
        
    generateFamilyDict()
            
    """Make a dictionary which gives information about gene order """

    geneOrder = {}
    for s in data.speciesDict:
        geneOrder[s] = {}
    
    for g in data.genesDict:
        currentGene = data.genesDict[g]
        try:
            geneOrder[currentGene.species][currentGene.ctg].append(currentGene)   #pass gene object
        except KeyError:
            geneOrder[currentGene.species][currentGene.ctg] = [currentGene]
                
    #sort genes in their contig
    for s in geneOrder:
        for c in geneOrder[s]:
            geneOrder[s][c].sort(key = lambda gene: int(gene.start))
        
    return data, geneOrder
    
    
