# -*- coding: utf-8 -*-
"""
GeneFamily class

Initialize with a family name and a list of species.
An object of this class contains a familyName and a dictionary indexed by species. Each index contains genes 
belonging to the family and that species.
Member function "addToFamily" adds a gene at the appropriate species index.



Created on Sun May  7 11:36:08 2017
@author: mfortz
"""



class GeneFamily(object):
    
    #find way to not pass listofspecies every time
    def __init__(self,name,listOfSpecies):
        self.familyName = name  
        self.familyMembersDict = {}
        for s in listOfSpecies:
            self.familyMembersDict[s]=[]
            
    
    def addToFamily(self,gene):
        if gene.family == self.familyName:
            self.familyMembersDict[gene.species].append(gene)
            


            
        

    
    
    
        
        
    
        
        
        
    


    
    