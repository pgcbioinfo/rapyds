#!/usr/bin/python

"""
RApyDS
Restriction site-associated DNA from Python-implemented Digestion Simulations
https://github.com/pgcbioinfo/rapyds

create_histogram.py
"""

from __future__ import print_function
import shutil, os, errno, sys
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import json
import os.path
import math
import argparse

def create_density_histogram_image(xbins,output_dir,enzyme,sequence,input_json, seq_name):
	file_open = open(output_dir+"/histogram.txt","a+")
	file_open.write("Sequence Name: "+seq_name+"\n")
	file_open.write("Enzyme: "+enzyme+"\n")
	with open(input_json) as jsonfile:
		data = json.load(jsonfile)
		seq_len = data[enzyme][0]
		cut_sites = data[enzyme][1:]

		num_bins = np.arange(0,seq_len+1,float(seq_len)/xbins)
		# if(args.nbin):
		# # 	num_bins = np.arange(0,seq_len+1,float(seq_len)/xbins)
		# # else:
		# num_bins = [x for x in range(0,seq_len+1,xbins)]
		# num_bins.append(xbins*((seq_len/xbins)+1))

		fig1, fig = plt.subplots()
		n, bins, patches = fig.hist(cut_sites, num_bins, facecolor='blue', alpha=0.5,  edgecolor='black', linewidth=0.75)
		fig.set_xticks(bins[::5])
		fig.ticklabel_format(style='plain',useOffset=False)
		plt.xlabel("Base Pairs")
		plt.ylabel("Number of Loci")
		plt.title(seq_name+" - "+enzyme)
		plt.legend(["Bin size: %d" % bins[1]])
		file_name = output_dir+"/images/"+sequence+"_"+enzyme+'.png'
		# file_name = file_name.replace(' ','-')
		plt.savefig(file_name)
		plt.cla()
		plt.close(fig1)

		# printing of histogram bins
		for i in range(0,len(num_bins)-1):
			if(i == 0):
				file_open.write("%d - %d \t %d\n" % (0,bins[i+1],n[i]))
			else:
				file_open.write("%d - %d \t %d\n" % (bins[i]+1,bins[i+1],n[i]))
		
	file_open.write("\n")
	file_open.close()

def create_density_histograms(bins,output_dir,input_path):
	

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

	file_open = open(output_dir+"/histogram.txt","w+")
	file_open.write("=======================\n")
	file_open.write("LOCI DENSITY\n")
	file_open.write("=======================\n")
	file_open.close()

	with open(input_path+"/genome_names.txt") as f:
		for content in f:
			line = content.strip().rstrip()
			sequence_name = line.split(" ")[0]

			for enzyme in renzymes:
				jsonfile = input_path+"/cut"+sequence_name+".json"
				create_density_histogram_image(bins,output_dir,enzyme,sequence_name,jsonfile,line)
	print("Finished creating histograms")


if __name__ == '__main__':
	global args
	parser = argparse.ArgumentParser(description='RApyDS loci density histogram plotter script')
	settings = parser.add_mutually_exclusive_group(required=True)
	settings.add_argument('-nbin', help='number of bins', action='store_true')
	settings.add_argument('-binsize', help='size of bins', action='store_true')
	parser.add_argument('bin', type=int, help='number of bins or bin size')
	parser.add_argument('input_dir', default='output', help='path of RApyDS output directory containing the json files')
	args = parser.parse_args()

	create_density_histograms(args.bin,args.input_dir,args.input_dir)
	

# create_histogram_image(20,'output',sys.argv[2], sys.argv[1])