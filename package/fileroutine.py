#!/usr/bin/python

import subprocess
import os
import re
import sys

import processing

def getColumns(inFile, header=True, delimiter="\t"):
    """
    Get columns of data from inFile. The order of the rows is respected
    """
    cols = {}
    indexToName = {}
    filetype = 0
    for lineNum, line in enumerate(inFile):
        if not line.startswith("##"):
            if lineNum == 0:
            	fields = line.split()
                if len(fields) == 1:
                    filetype = 1 ##netMHC.xls output
                else:
                    filetype = 0 ##sequence input
                    cols, indexToName = getHeadings(line,header)
                            	
            elif filetype:
            	cols, indexToName = getHeadings(line,header)
            	filetype = 0 ##already read header, switch to normal read

            else:	
                cells = line.split()
                i = 0
                for cell in cells:
                    cell = cell.strip()
                    cols[indexToName[i]] += [cell]
                    i += 1

    return cols, indexToName


def getHeadings(line, header):
	indexToName = {}
	cols = {}
	i = 0
	if header:
		headings = line.split()
		for heading in headings:
			heading = heading.strip()
			indexToName[i] = heading
			cols[heading] = []
			i += 1
	else:
		fields = line.split()
		for field in fields:
			indexToName[i] = i
			cols[i] = [field]
			i += 1

	return cols, indexToName

##Find allele name and initiate a list
def findAllele(indexFile):
	AlleleName = {}
	for line in indexFile:
		cols = line.split()
		##item = re.sub('[*:]','',cols[0])
		item = cols[0].strip()
		AlleleName[item] = 0
		##AlleleName[item] = 1
	
	return AlleleName

##Process the input sequence file
def InputSeq(AlleleName, Dir):
	print("Converting Sequence File....")
	for key, val in AlleleName.iteritems():
		data = []
		data.append([])
		data.append([])

		if os.path.getsize(os.path.join(Dir.input_dir,"HLA-%s.txt") % key):
			AlleleName[key] = 1
			with open(os.path.join(Dir.input_dir,"HLA-%s.txt") % key,'r') as f:
				for line in f:
					cols = line.split()
					data[0].append(cols[1].strip())
					data[1].append(cols[2].strip())
		else:
			print("File empty; no sequence converted for allele:%s" % (key))
			AlleleName[key] = 0
    		
		if AlleleName[key]:
			with open(os.path.join(Dir.tmp_dir,"HLA-%s.txt") % key,'w') as outFile:
				print("Writing File:%s" % (key))
				for i in range(len(data[0])):
					pep_flag, index_prev = processing.checkRepeat(data[0][i],data[0][:i])
					if not pep_flag:
						outFile.write("%s\t%.3f\n" % (data[0][i], float(data[1][i])))
					else:
						pep_flag = False

def writeRdata(PepData, dataout):
	dataout.write("peptide" + "\t"
		+ "meas_nm" + "\t" 
		+ "meas_bi" + "\t" 
		+ "meas_contin" + "\t" 
		+ "predict" + "\t" 
		+ "predict_rank" + "\n")
        				
	for i in range(len(PepData.name)):
		dataout.write("%s\t%.3f\t%d\t%.3f\t%.3f\t%.3f\n" % (PepData.name[i], PepData.meas_nm[i],
    		PepData.meas_bi[i], PepData.meas_contin[i], 
    		PepData.predict[i], PepData.predict_rank[i]))


