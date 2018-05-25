# -*- coding: utf-8 -*-
"""
computeFDR.py

Created on Thu Jul 27 15:15:55 2017

@author: mfortz
"""

from __future__ import division
import sys

#input observed (candidate) and null (shuffled) scores
observedIn = open(sys.argv[1],'r').readlines()
nullIn     = open(sys.argv[2],'r').readlines()
#set FDR
setFDR = float(sys.argv[3])
#output accepted pairs
outputFile = open(sys.argv[4],'w')

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

observedScores.sort(reverse = True)
nullScores.sort(reverse = True)

def computeFDR(current_t, previousScores_n, previousScores_b):
    #take scores above threshold
    trimmed_n = []
    trimmed_b = []
    
    for score in previousScores_n:
        if score >= current_t:
            trimmed_n.append(score)
        else:
            break  #all scores after is <t since list is reverse sorted
 
    for score in previousScores_b:
        if score >= current_t:
            trimmed_b.append(score)
        else:
            break  #all scores after is <t since list is reverse sorted     
            
    s_n = len(trimmed_n)
    s_b = len(trimmed_b)
    
    if s_b == 0:
        print("Unable to reach set FDR.")
        return -1, trimmed_n, trimmed_b

    currentFDR = (s_n/s_b) * (b/n)
    return currentFDR, trimmed_n, trimmed_b

# compute t that satisfies setFDR
trimmed_n = [z for z in nullScores]
trimmed_b = [z for z in observedScores]

t = -0.1   #initialize values
FDR = 1.0

while FDR > setFDR:
    t += 0.001
    FDR, trimmed_n, trimmed_b = computeFDR(t, trimmed_n, trimmed_b)
    
fdrOut = "%.4f" % (FDR*100)
tOut = "%.4f" % (t-0.001)

#print computed FDR and t
print("Computed treshold t = " + tOut + " with FDR = " + fdrOut +"%")
print("Observed "+ str(len(trimmed_b)) + " pairs above threshold.")

acceptedPairs = []

#find accepted family pairs
for l in observedIn:
    if l[0]!='#':
        x = l.split()
        d = float(x[7])
        if d <= (t-0.001):
            pair = (x[0],x[1])
            acceptedPairs.append(pair)

for p in acceptedPairs:
    outputFile.writelines([p[0],'\t',p[1],'\n'])

        
