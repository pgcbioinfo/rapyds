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

def create_histogram_image(bins,path,enzyme,output_json):
	with open(output_json) as jsonfile:
		data = json.load(jsonfile)
		cut_sites = data[enzyme]
		num_bins = bins
		n, bins, patches = plt.hist(cut_sites, num_bins, facecolor='blue', alpha=0.5,  edgecolor='black', linewidth=0.75)
		plt.xlabel(enzyme)
		plt.savefig(path+'/output_'+enzyme+'.png')

create_histogram_image(20,'output',sys.argv[2], sys.argv[1])