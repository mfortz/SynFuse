# -*- coding: utf-8 -*-
"""
DataSet class

Contains three dictionaries: "genesDict", "speciesDict" and "familiesDict". 
These dictionaries contain the appropriate objects as values and their names as index..

Contains functions to add a **new index and an appropriate object** to each dictionary. Note that "addFamily" takes a 
GeneFamily object as an argument because such an object requires the complete speciesList to be created beforehand. 

Contains functions to print all genes (and gene information) of a specific family or all families.

Contains a function to create a dictionary which is a subset of "genesDict". Takes a list of family names as argument.
                                     

Created on Sun May  7 11:53:11 2017
@author: mfortz
"""

from Gene import Gene
from GeneFamily import GeneFamily
from Species import Species


class DataSet(object):
    
    def __init__(self):
        self.familiesDict = {}
        self.speciesDict = {}
        self.genesDict = {}
        
        
    def addFamily(self,familyName,familyObject):
        self.familiesDict[familyName] = familyObject

               
    def addSpecies(self,species):
        species = Species(species)
        self.speciesDict[species.speciesName] = species
        
    def addGene(self,lineFromFile):
        gene = Gene(lineFromFile)
        self.genesDict[gene.gene] = gene    #not a list since gene is unique
        
       
    def subset(self,famNameList):
        subsetDict = {}
        for f in famNameList:
            subsetDict[f] = self.familiesDict[f]
        return subsetDict
    
    
    def printFamily(self,familyName):
        currentFamDict = self.familiesDict[familyName].familyMembersDict
        for species in currentFamDict:
            for gene in currentFamDict[species]:
                gene.printGene()
        print()
    
    def printAllFamilies(self):
        for f in self.familiesDict:
            self.printFamily(f)
        ###need to show empty species?






#        
    
        
        
    
        