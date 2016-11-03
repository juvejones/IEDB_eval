#!/usr/bin/python

import subprocess,re,os,sys,math

import fileroutine
import processing
from container import PepData

"""
Module for running RESTApi on IEDB
"""

def runPredict(methods, allele_list, Dir):
	for method in methods:
		print("RUN {0}----->".format(method))
		if not os.path.exists(str(method)):
			os.mkdir(str(method))

		output_path = os.path.join(str(method),Dir.output_dir)
		if not os.path.exists(output_path):
			os.mkdir(output_path)

		dataout_path = os.path.join(str(method),Dir.data_dir)
		if not os.path.exists(dataout_path):
			os.mkdir(dataout_path)
		summaryf = open(dataout_path+'/summary.txt','w')

		for key, value in allele_list.iteritems():
		    print(key)
		    allele_rename = re.sub('[*:]','',key)
		    summaryf.write("###%s\n" % allele_rename)
		    ## input sequence is in two-column format, converted by fileroutine.InputSeq
		    ## and stored in tmp_dir
		    input_path = Dir.tmp_dir 
		    output = output_path+'/HLA-{0}.txt'.format(allele_rename)
			
		    with open(input_path+'/HLA-{0}.txt'.format(key),'r') as file:
		        lines = file.readlines()
		        num_lines = sum(1 for line in lines)
		        splitlen = 50

		        if num_lines > 100:
		            outf = open(output,'w')
		            batch_lines = []
		            j = 0
		            for l_start in range(0, num_lines, splitlen):
		                batch_lines.append([])
		                batch_lines[j] = lines[l_start:l_start+splitlen]
		                lastline = batch_lines[j][-1]

		                i = 0
		                sequences = "%3E"
		                for line in batch_lines[j]:
		                	i += 1
		                	cols = line.split()
		                	sequences = sequences + "peptide" + str(i) + "%0A" + str(cols[0])
		                	if not line == lastline:
		                		sequences = sequences + "%0A%3E"	

		                allele = "HLA-" + key
		                data = "method="+str(method)+"&sequence_text="+str(sequences)+"&allele="+str(allele)+"&length=9"
		                subprocess.call([
		                	'curl',
		                	'-d',
		                	data,
		                	'http://tools-api.iedb.org/tools_api/mhci/',
		                	'-o',
		                	'temp_output'
		                	])
		                temp = open('temp_output','r')
		                if j == 0:
		                	outputData = temp.readlines()
		                else:
		                	header = temp.readline()
		                	outputData = temp.readlines()
		                temp.close()
		                outf.write("".join(outputData))
		                j += 1

		            outf.close()
		        else:
		            i = 0
		            sequences = "%3E"
		            lastline = lines[-1]
		            for line in lines:
		                i += 1
		                cols = line.split()
		                sequences = sequences + "peptide" + str(i) + "%0A" + str(cols[0])
		                if not line == lastline:
		                	sequences = sequences + "%0A%3E"	

		            allele = "HLA-" + key
		            data = "method="+str(method)+"&sequence_text="+str(sequences)+"&allele="+str(allele)+"&length=9"
		            subprocess.call([
		            	'curl',
		            	'-d',
		            	data,
		            	'http://tools-api.iedb.org/tools_api/mhci/',
		            	'-o',
		            	output
		            	])
		                
		        file.seek(0)
		        cols, indexToName = fileroutine.getColumns(file, header=False)

		    if not os.path.isfile(output):
		    	print("%s output not found for HLA-%s" % (method, allele_rename))

		    else:
		    	with open(output,'r') as outf:
		    		firstline = outf.readline()
		    		if re.match("Invalid allele", firstline):
		    		    print("allele not available for {0}\n".format(allele_rename) )
		    		    continue 
		    		elif re.match("<!DOCTYPE",firstline):
		    		    print("network error for {0}\n".format(allele_rename) )
		    		    continue
		    		else:
		    		    print("Calculation Finished!")
		    		    print("Now Process Statistics")
		    		    outf.seek(0)
		    		    cols_pred, indexToName_pred = fileroutine.getColumns(outf, header=True)

		    	peptides = cols[indexToName[0]]
		    	if allele_rename == "B2705-Flower":
		    		meas_cols = processing.convertPBL50(cols[indexToName[1]])
		    	else:
		    		meas_cols = cols[indexToName[1]]
		    	if method == "consensus":
		    		pred_cols = processing.twoColsMean(cols_pred[indexToName_pred[7]], 
		    			cols_pred[indexToName_pred[9]])
		    		pred_rank = cols_pred[indexToName_pred[6]]
		    	else:
		    		pred_cols = cols_pred[indexToName_pred[6]]
		    		pred_rank = cols_pred[indexToName_pred[7]]
		    	
		    	indxcols = cols_pred[indexToName_pred[1]]
		    	pred_cols = processing.reorder(indxcols, pred_cols)
		    	pred_rank = processing.reorder(indxcols, pred_rank)

		    with PepData(len(peptides)) as peptidedata:
		    	for index,val in enumerate(peptides):
		    		peptidedata.name[index] = val
		    		processing.getData(index, peptidedata, 
		    			meas_cols, pred_cols, pred_rank)
		    
		    	with open(dataout_path+'/HLA-{0}.txt'.format(allele_rename),'w') as dataoutfile:
		    		fileroutine.writeRdata(peptidedata, dataoutfile)

		    	fileroutine.writeRdata(peptidedata, summaryf)   

        summaryf.close()			

