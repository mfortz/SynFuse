#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
silhoutte.py

Calculates the silhouette of gene families.
    
Usage: 
    $ python silhouette.py [args]
    
See README for argument lists


Author: Mark Forteza, mforteza@sfu.ca
Last modified: 28 Jul 2017

"""

from __future__ import division

import sys, random
from deriveGeneStructure import deriveGeneStructure
#from silhouette_functions import clusterSilhouette

"""Helper functions"""

# Create dictionary that maps gene pairs to NC scores
def createScoreDict(scoreFile):
    scoreDict = {}
    
    for l in scoreFile:
        if l[0]!='#':
            x = l.split()
            g1 = x[0]
            g2 = x[1]
            NC = float(x[2])
            scoreDict[(g1,g2)] = NC
    
    return scoreDict

#Create dictionary that maps family pairs to edges between them
#Note that GenFam.syc contains two entries for each pair: (g1,g2) and (g2,g1). 
#These have the same scores. Therefore, key (f1,f2) is equal to (f2,f1), and 
#(f1,f1) image needs to be halved.
def computeNumEdges(scoreDict):
    numEdgesDict = {}
    
    for g1,g2 in scoreDict:
        f1 = data.genesDict[g1].family
        f2 = data.genesDict[g2].family
        
        try:
            numEdgesDict[(f1,f2)] += 1
        except KeyError:
            numEdgesDict[(f1,f2)] = 1
    
    return numEdgesDict


# Gives a list of family members for each family
def createMembersList(data):
    membersList = {}
    
    for famName in data.familiesDict:
        currentFam = data.familiesDict[famName]
        geneObjectList = []
        
        for s in data.speciesDict:
            geneObjectList += currentFam.familyMembersDict[s]
            
        membersList[famName] = [g.gene for g in geneObjectList]
        
    return membersList


"""Silhouette computations"""

# ave similarity of x to vertices in the same family
def geneCohesion(x,ownFam):
    global membersList, scoreDict, data
    
    familyList = membersList[ownFam]
    NCsum = 0
    
    if len(familyList) > 1:
        for g in familyList:
            if (x,g) in scoreDict:
                NC = scoreDict[(x,g)]
                NCsum += NC

        NCcohesion = NCsum/(len(familyList)-1)
    else:
        NCcohesion = 0.00
    
    return NCcohesion

# ave similarity of x to vertices in a different family
def geneSeparation(x,otherFam):
    global membersList, scoreDict, data
    
    familyList = membersList[otherFam]
    NCsum = 0
    
    for g in familyList:
        if (x,g) in scoreDict:
            NC = scoreDict[(x,g)]
            NCsum += NC
            
    NCseparation = NCsum/len(familyList)
    return NCseparation

# (cohesion - separation) / max(cohesion,separation)
# joinedFams=None if computing on a true cluster, else if computing for joined
# families joinedFams = (f1,f2)
def geneSilhouette(x,ownFam,joinedFams):
    global membersList, scoreDict, data

    NCcohesion = geneCohesion(x,ownFam)
    NCseparationList = []
    
    #compute separtion on all families except the one containing x.
    if not joinedFams:
        for otherFam in data.familiesDict:
            if ownFam != otherFam:
                NC = geneSeparation(x,otherFam)
                NCseparationList.append(NC)
    else:
        for otherFam in data.familiesDict:
            if joinedFams[0] != otherFam and joinedFams[1] != otherFam:
                NC = geneSeparation(x,otherFam)
                NCseparationList.append(NC)
        
    #use the max separation among families
    NCseparation = max(NCseparationList)
    
    #compute silhouette
    try: 
        NCsilhouette = (NCcohesion - NCseparation)/ max(NCcohesion,NCseparation)
    except ZeroDivisionError:
        NCsilhouette = 0

    return NCsilhouette


# Compute min,max,ave silhouette of a cluster
# joinedFams=None if computing on a true cluster, else if computing for joined
# families joinedFams = (f1,f2)
def clusterSilhouette(family, joinedFams=None):
    global membersList, scoreDict, data

    familyList = membersList[family]
    NCsil = []      #record silhouette of every gene for output   
    
    for g in familyList:
        NC = geneSilhouette(g,family,joinedFams)
        NCsil.append((g,NC))
    
    # MIN,MAX,AVE
    NCsilList = [score[1] for score in NCsil]
    
    NCsilMin = min(NCsilList)
    NCsilMax = max(NCsilList)
    NCsilAve = sum(NCsilList)/len(NCsilList)
    
    valuesOutNC = [NCsilMin, NCsilMax, NCsilAve]
    
    return NCsil, valuesOutNC

"""Computes silhouettes for each family and prints output"""
def computeAllFamSilhouette(membersList, numEdgesDict, ncOut):
    ncOut.write("#family \t num_vertices \t num_edges \t min,max,ave silhouette")
    ncOut.write("\t gene silhouette \n")

    for famName in data.familiesDict:
        numV = len(membersList[famName])
        
        #count number of edges inside a cluster. Divided by 2 because for
        #a,b in f, both nc(a,b) and nc(b,a) are recorded in GenFam.nnc
        try:
            numE = numEdgesDict[(famName, famName)]//2
        except KeyError:
            numE = 0
        
        #call silhouette calculations
        NCsil, valuesOutNC = clusterSilhouette(famName)

        #print header        
        ncOut.writelines([famName,'\t',str(numV),'\t',str(numE),'\t'])
        
        #write min,max,ave sil on a cluster and sil for every member
        for value in valuesOutNC:
            short = "%.6f" % value
            ncOut.write(str(short))
            ncOut.write('\t')  
                
        for (gene,sil) in NCsil:
            short = "%.6f" % sil
            ncOut.writelines([gene,'\t',str(short),'\t'])
        ncOut.write('\n')
        
    #close to ensure all outputs are printed
    ncOut.close() 

"""Join families, compute silhouette, and print output"""
def computeJoinedSilhouette(pairsList, famSilhouette, ncOut):
    # print header
    ncOut.write("#F1 \t F2 \t num_edges_between \t s1 \t s2 \t s12 \t s \t d \n")
    
    global membersList
    
    # function for joining families
    def joinFams(famName1, famName2):
        if famName1 == famName2:
            print("Error. Attempting to join same family.")
            return
        
        #combine families
        membersJ = membersList[famName1] + membersList[famName2]
        #members1 = membersList[famName1]
        #members2 = membersList[famName2]
        #membersJ = members1 + members2
        
        #create new membersList with the f1,f2 joined and the originals deleted
        revisedMembersList = membersList.copy()
        revisedMembersList["joined"] = membersJ
        del revisedMembersList[f1]
        del revisedMembersList[f2]
    
        return revisedMembersList

    #Compute joined silhouettes
    for pair in pairsList:
        f1 = pair[0]
        f2 = pair[1]
        s1 = famSilhouette[f1]
        s2 = famSilhouette[f2]
        
        #copy the true membersList and make new one with f1,f2 joined. 
        #This method is used because membersList is a global variable on
        #silhouette computation funtions
        trueMembersList = membersList.copy()
        membersList = joinFams(f1,f2)
    
        #count number of edges between f1 and f2
        try:
            numEdgesBetween = numEdgesDict[(f1,f2)]
        except KeyError:
            numEdgesBetween = 0
    
        #call silhouette computations
        NCsil, valuesOutNC = clusterSilhouette("joined",(f1,f2))
        
        #reassign the true membersList
        membersList = trueMembersList.copy()
        
        #compute scaled average( sil(f1),sil(f2) )
        #scale by size of respective families
        sizef1 = len(membersList[f1])
        sizef2 = len(membersList[f2])
        silSumf1 = s1 * sizef1
        silSumf2 = s2 * sizef2
        s= (silSumf1 + silSumf2) / (sizef1 + sizef2)
        
        #average silhouette of the cluster after joining
        s12 = valuesOutNC[2]
        
        #change in silhouette 
        d = s12 - s
        
        #print output
        values = [s1, s2, s12, s, d]
        ncOut.writelines([f1,'\t',f2,'\t',str(numEdgesBetween),'\t'])
    
        for v in values:
            short = "%.6f" % v
            ncOut.write(str(short))
            ncOut.write('\t')
        ncOut.write('\n')

    #close to ensure all outputs are printed
    ncOut.close()

"""Shuffle edge weights"""

#intercluster edges shuffled (not only weights)
#insideShuffled = 0 or 1
def shuffleWeights(scoreDict,insideShuffled):
    
    intraClusterEdges = []
    intraClusterWeights = []
    #interClusterEdges = []
    interClusterWeights = []
    
    global data   
        
    #record processed gene pairs
    doneDict = {}
    
    for x,y in scoreDict:
        g1 = data.genesDict[x]
        g2 = data.genesDict[y]
        f1 = g1.family
        f2 = g2.family
        
        if f1 == f2:
            try:
                doneDict[(x,y)]   #KeyError if (x,y) has not been processed
                continue
            except KeyError:
                intraClusterEdges.append((x,y))
                intraClusterWeights.append(scoreDict[(x,y)])
                #avoid duplicate scores, i.e. score(x,y) and score(y,x)
                doneDict[(x,y)] = None
                doneDict[(y,x)] = None
        else:
            try:
                doneDict[(x,y)]   #KeyError if (x,y) has not been processed
                continue
            except KeyError:
                #interClusterEdges.append((x,y))
                interClusterWeights.append(scoreDict[(x,y)])
                #avoid duplicate scores, i.e. score(x,y) and score(y,x)
                doneDict[(x,y)] = None
                doneDict[(y,x)] = None
    
    #shuffle edges
    shuffledScoreDict = {}
    
    """New way """
    geneList = list(data.genesDict.keys())    
    
    for weight in interClusterWeights:
        sameFam = True        
        
        while sameFam:
            x,y = random.sample(geneList,2)
            g1 = data.genesDict[x]
            g2 = data.genesDict[y]
            if g1.family != g2.family:
                sameFam = False
                
        #want score(x,y) = score(y,x)
        shuffledScoreDict[(x,y)] = weight
        shuffledScoreDict[(y,x)] = weight
        
    #shuffle intracluster weights
    if insideShuffled == 1:
        #random.shuffle(intraClusterEdges)
        random.shuffle(intraClusterWeights)        
        
        for i in range(len(intraClusterEdges)):
            x,y = intraClusterEdges[i]
            weight = intraClusterWeights[i]
            #want score(x,y) = score(y,x)
            shuffledScoreDict[(x,y)] = weight
            shuffledScoreDict[(y,x)] = weight
    else:
        for x,y in intraClusterEdges:
            shuffledScoreDict[(x,y)] = scoreDict[(x,y)]
            shuffledScoreDict[(y,x)] = scoreDict[(y,x)]
        
    return shuffledScoreDict

###############################################################
'''Main'''

'''Parameters'''

NC_SCORE_FILE        = sys.argv[1]
GENE_FILE            = sys.argv[2]
CANDIDATES_FILE      = sys.argv[3]
NB_SAMPLES           = int(sys.argv[4])
SILHOUETTE_ALL       = sys.argv[5]
SILHOUETTE_JOINED    = sys.argv[6]
SILHOUETTE_SH_ALL    = sys.argv[7]
SILHOUETTE_SH_JOINED = sys.argv[8]

'''all_families'''

print("Computing initial clustering silhouette...")

#input edge weights
scoreFile = open(NC_SCORE_FILE, 'r').readlines()

# load gene structure
all_gene_file = open(GENE_FILE,'r').readlines()
data, geneOrder = deriveGeneStructure(all_gene_file)
    
#output all family silhouette 
silhouetteOut = open(SILHOUETTE_ALL,'w')

membersList = createMembersList(data)
scoreDict = createScoreDict(scoreFile)
numEdgesDict = computeNumEdges(scoreDict)
computeAllFamSilhouette(membersList, numEdgesDict, silhouetteOut)

'''joined'''

print("Computing joined cluster silhouettes...")

#input list of family pairs (eg candidates, null pairs)
pairsFile = open(CANDIDATES_FILE,'r').readlines()
#input computed all fam silhouette
silhouetteIn = open(SILHOUETTE_ALL,'r').readlines()    
#output all family silhouette 
joinedOut = open(SILHOUETTE_JOINED,'w')

#collect family pairs in a list
pairsList = []
for l in pairsFile:
    if l[0]!='#':
        x = l.split()
        pairsList.append((x[0],x[1]))

#dictionary that maps family name to ave nc silhouette
famSilhouette = {}
for l in silhouetteIn:
    if l[0]!='#':
        x = l.split()
        fam = x[0]
        ncSil = float(x[5])   #ave nc silhouette
        famSilhouette[fam] = ncSil

computeJoinedSilhouette(pairsList,famSilhouette,joinedOut)

'''shuffled'''

print("Shuffling edges...")

#output all family silhouette 
silhouetteOutShuf = open(SILHOUETTE_SH_ALL,'w')
# Shuffle weights 
insideShuffled = 1 # Shuffle intra cluster weights? 
scoreDict= shuffleWeights(scoreDict,insideShuffled)

print("Computing silhouette of shuffled clusters...")
# Compute silhouettes of shuffled clusters
computeAllFamSilhouette(membersList, numEdgesDict, silhouetteOutShuf)

print("Computing silhouette of joined shuffled clusters...")
# Reload the just computed all fam silhouette file
silhouetteInShuf = open(SILHOUETTE_SH_ALL,'r').readlines()

# Find all fam pairs with an edge between
pairsList = list(set(numEdgesDict.keys()))

# Eliminate intracluster entries
tempList = []
for x,y in pairsList:
    if x != y:
        tempList.append((x,y))
pairsList = list(tempList)

#eliminate entry (y,x) if (x,y) is present
tempList = []
doneMap = {}
for x,y in pairsList:
    try:
        doneMap[(x,y)]
    except KeyError:
        doneMap[(x,y)] = None
        doneMap[(y,x)] = None
        tempList.append((x,y))
pairsList = list(tempList)
        
#dictionary that maps family name to ave silhouette
famSilhouetteShuf = {}
for l in silhouetteInShuf:
    if l[0]!='#':
        x = l.split()
        fam = x[0]
        ncSil = float(x[5])   #ave nc silhouette
        famSilhouetteShuf[fam] = ncSil

#try not to hard code this    
pairsList = random.sample(pairsList,NB_SAMPLES)

#output joined family silhouette
joinedOutShuf = open(SILHOUETTE_SH_JOINED,'w')
#compute silhouette
computeJoinedSilhouette(pairsList,famSilhouetteShuf,joinedOutShuf)
