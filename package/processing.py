#!/usr/bin/python

import math
import random
import sys

MAX_AFF = 50000.0
THRESHOLD = 500.0

"""
Module for processing input and output data
Possibly including R functions in the future
"""

def checkAffinity(affinity):
	if affinity == 0:
		affinity = 1000000
	if affinity < 1.000001:
		affinity = float(random.randint(1,10))

	return affinity

def scalePrediction(affinity):
    score = 1 - math.log(checkAffinity(affinity))/math.log(MAX_AFF)
    if score < 0:
        score = 0

    return score

##def rank(affinity):

def findMax(array):
    max = 1.0
    for idx,val in enumerate(array):
        val = float(val)
        if val > max:
        	max = val
        
    return max

def convertPBL50(list):
    result = []
    for i in range(len(list)):
        pBL = -float(list[i])
        convert_value = math.pow(10,pBL)*math.pow(10,9)
        result.append(convert_value)
    
    return result
    
def checkRepeat(current,list):
	exist = False
	index_prev = 0
	for index, value in enumerate(list):
		if current == value:
			exist = True
			index_prev = index

	return exist, index_prev

def getDataRepeat(idx,idx_prev,PepData,list1,list2,list3):
    meas = float(list1[idx])
    meas_prev = float(list1[idx_prev])
    meas_new = 0.5*(meas+meas_prev)
    pred = float(list2[idx])
    ##max_affinity = findMax(list2)
    if meas_new < THRESHOLD:
        PepData.meas_bi[idx] = 1
    else:
        PepData.meas_bi[idx] = 0

    PepData.meas_contin[idx] = scalePrediction(meas_new)
    PepData.meas_contin[idx_prev] = scalePrediction(meas_new)
    PepData.predict[idx] = scalePrediction(pred)
    PepData.predict_rank[idx] = float(list3[idx])/100

def getData(idx,PepData,list1,list2,list3): 
    meas = float(list1[idx])
    pred = float(list2[idx])
    ##max_affinity = findMax(list2)
    if meas < THRESHOLD:
        PepData.meas_bi[idx] = 1
    else:
        PepData.meas_bi[idx] = 0
    
    PepData.meas_nm[idx] = meas
    PepData.meas_contin[idx] = scalePrediction(meas)
    PepData.predict[idx] = scalePrediction(pred)
    PepData.predict_rank[idx] = float(list3[idx])/100

def twoColsMean(list1, list2):
    outlist = [None]*len(list1) 
    if not (len(list1) == len(list2)):
        print("two columns not same length!!")
        raise
    else:
        for i in range(len(list1)):
            outlist[i] = math.exp((math.log(float(list1[i])) 
                + math.log(float(list2[i])))/2)

    return outlist

def reorder(reflist, targetlist):
    orderedlist = [None]*len(reflist)
    repeat = set()
    uniq = []
    multiple = {}
    for idx, x in enumerate(reflist):
        if x in repeat:
            multiple[x] += 1
            reflist[idx] = int(x) + 50*multiple[x]
        else:
            multiple[x] = 0
            repeat.add(x)
            uniq.append(x)
            
    d = dict(zip(reflist,targetlist))
    for key in sorted(d):
        i = int(key) - 1
        orderedlist[i] = d[key]

    return orderedlist
