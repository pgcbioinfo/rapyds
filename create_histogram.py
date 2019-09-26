#!/usr/bin/python

"""
RApyDS
Restriction Site Associated DNA Python-Digested Simulation 

create_histogram.py
"""


from __future__ import print_function
import shutil, os, errno, sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import json
import os.path

def create_histogram_image(xbins,output_path,enzyme,sequence,input_json, seq_name):
	with open(input_json) as jsonfile:
		data = json.load(jsonfile)
		seq_len = data[enzyme][0]
		cut_sites = data[enzyme][1:]
		num_bins = np.arange(0,seq_len+1,float(seq_len)/xbins)
		fig1, fig = plt.subplots()
		n, bins, patches = fig.hist(cut_sites, num_bins, facecolor='blue', alpha=0.5,  edgecolor='black', linewidth=0.75)
		fig.set_xticks(bins[::5])
		plt.xlabel("Base Pairs")
		plt.ylabel("Number of Loci")
		plt.title(seq_name+" - "+enzyme)
		plt.savefig(output_path+'/'+sequence+"_"+enzyme+'.png')

def create_histograms(bins,output_dir,input_path):
	## clean up and making directories
	output_path = output_dir + "/images/"
	if(os.path.exists(output_path) == True):
		shutil.rmtree(output_path)
	os.makedirs(output_path)

	if (os.path.exists(input_path+"/genome_names.txt") == False):
		print("No genome_names.txt file in "+input_path)
		raise SystemExit

	if (os.path.exists(input_path+"/RE.txt") == False):
		print("No RE.txt file in "+input_path)
		raise SystemExit

	renzymes = []
	with open(input_path+"/RE.txt") as f:
		for content in f:
			line = content.strip().rstrip()
			renzymes.append(line)

	with open(input_path+"/genome_names.txt") as f:
		for content in f:
			line = content.strip().rstrip()
			sequence_name = line.split(" ")[0]

			for enzyme in renzymes:
				jsonfile = input_path+"/cut"+sequence_name+".json"
				create_histogram_image(bins,output_dir+"/images/",enzyme,sequence_name,jsonfile,line)


# create_histogram_image(20,'output',sys.argv[2], sys.argv[1])