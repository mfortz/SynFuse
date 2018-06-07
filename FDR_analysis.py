# -*- coding: utf-8 -*-
"""
computeFDR.py

Created on Thu Jul 27 15:15:55 2017

@author: mfortz
"""

from __future__ import division
import sys

#input observed (candidate) and null (shuffled) scores
observedIn = open(sys.argv[1],'r').readlines()[1:]
nullIn     = open(sys.argv[2],'r').readlines()[1:]
#set FDR
setFDR = float(sys.argv[3])
#output accepted pairs
outputFile = open(sys.argv[4],'w')


'''extract silhouette difference scores from input files'''
observedScores = []
for l in observedIn:
    if l[0]!='#':
        x = l.split()
        d = float(x[7])
        observedScores.append(d)
    
nullScores = []
for l in nullIn:
    if l[0]!='#':
        x = l.split()
        d = float(x[7])
        nullScores.append(d)

#data set sizes
b = len(observedScores)
n = len(nullScores)

#sort score lists from highest to lowest
observedScores.sort(reverse = True)
nullScores.sort(reverse = True)


'''main FDR computation'''
def computeFDR(current_t, previousScores_n, previousScores_b):
    #take scores above threshold
    trimmed_n = []
    trimmed_b = []
    
    #take null scores above threshold
    for score in previousScores_n:
        if score >= current_t:
            trimmed_n.append(score)
        else:
            break  #all scores after is <t since list is reverse sorted
    
    #take observed scores above threshold
    for score in previousScores_b:
        if score >= current_t:
            trimmed_b.append(score)
        else:
            break  #all scores after is <t since list is reverse sorted     
            
    '''compute FDR, ratio of observed over null scores above threshold'''
    s_n = len(trimmed_n)
    s_b = len(trimmed_b)
    
    if s_b == 0:
        print("Unable to reach set FDR.")
        return -1, trimmed_n, trimmed_b

    #scompute FDR, scale by size of data sets
    currentFDR = (s_n/s_b) * (b/n)  
    return currentFDR, trimmed_n, trimmed_b


''' compute t that satisfies setFDR'''

#score lists copies that will be trimmed at each call of computeFDR(..)
trimmed_n = [z for z in nullScores]
trimmed_b = [z for z in observedScores]

#initialize values
t = -0.1   
FDR = 1.0

#increase threshold until target FDR is achieved
while FDR > setFDR:
    t += 0.001
    FDR, trimmed_n, trimmed_b = computeFDR(t, trimmed_n, trimmed_b)
    
fdrOut = "%.4f" % (FDR*100)
tOut = "%.4f" % t

#print computed FDR and t
print("Computed treshold t = " + tOut + " with FDR = " + fdrOut +"%")
print("Observed "+ str(len(trimmed_b)) + " pairs above threshold.")


'''record accepted pairs'''
acceptedPairs = []

#find accepted family pairs
for l in observedIn:
    if l[0]!='#':
        x = l.split()
        d = float(x[7])
        if d >= (t):
            pair = (x[0],x[1],"%.4f" % d)
            acceptedPairs.append(pair)

#print to output file
outputFile.writelines(['#','F1','\t','F2','\t','d','\n'])

for p in acceptedPairs:
    outputFile.writelines([p[0],'\t',p[1],'\t',p[2],'\n'])

        
