# -*- coding: utf-8 -*-
"""
Gene class

Conatins information about a single gene from the data file. Allows for formatted printing of this information.
Note that some  information may be missing from the file. In that case, the appropriate member variable will contain 
the string "No data".


Created on Thu May  4 14:18:24 2017
@author: mfortz
"""


class Gene(object):
    
    def __init__(self,data = [None]*9):
        self.species = data[0]
        self.ctg = data[1]
        self.family = data[2]
        self.gene = data[3]
        self.orient = data[4]
        self.start = data[5]
        self.end = data[6]
        self.exons = data[7]
        try:
            self.exons_pos = data[8]
        except IndexError:
            self.exons_pos = 'No data'
  

    def printGene(self):
        print(self.species, end="\t")
        print(self.ctg, end="\t")
        print(self.family, end="\t")
        print(self.gene, end="\t")
        print(self.orient, end="\t")
        print(self.start, end="\t")
        print(self.end, end="\t")
        print(self.exons, end="\t")
        print(self.exons_pos, end="\t")
        print()  
        
        
        