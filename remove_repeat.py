#!/usr/bin/python

"""
RApyDS
Restriction site-associated DNA from Python-implemented Digestion Simulations
https://github.com/pgcbioinfo/rapyds

remove_repeat.py
"""

from __future__ import print_function
import sys, os
from operator import itemgetter

def remove_XAs(enzyme, genome_name,bwa_dir):
	if(os.path.isdir(bwa_dir) == False and os.path.exists(bwa_dir+"/aligned_pairs_"+genome_name+"_"+enzyme+".sam")==False) :
		print("Run bwa_aln.sh first")
		raise SystemExit

	with open(bwa_dir+"/aligned_pairs_"+genome_name+"_"+enzyme+".sam","r") as orig:
		try:
			fragments = []
			XAs = []
			total = 0
			unique = 0
			repeat = 0
			unique_frag_num = []
			repeat_frag_num = []
			done_frag_num = []
			# Check through lines of sam file
			for line in orig:
				# Get the LN part of sam file
				if("LN:" in line):
					total = int(line.split("\t")[2][3:])
				if("Frag_" not in line):
					continue

				line = line.strip().rstrip().split("\t")
				frag_name_split = line[0].split("_")
				frag_number = frag_name_split[1]
				## all fragments's start and end are in this list fragments
				fragments.append([int(frag_name_split[2]),int(frag_name_split[3])])

				
				## get all XAs
				if(len(line) > 21):
					list_XA = line[21].split(":")[2].split(";")[:-1]
					for alt_hits in list_XA:
						alt = alt_hits.split(",")
						## each XA's location is in XAs list
						XAs.append(int(alt[1]))
				for x in line:
					if x=="XT:A:U" in line[11]:
						if frag_number not in done_frag_num:
							unique_frag_num.append(frag_number)
							done_frag_num.append(frag_number)
						else:
							unique += 1
					elif x=="XT:A:R" in line[11]:
						if frag_number not in done_frag_num:
							repeat_frag_num.append(frag_number)
							done_frag_num.append(frag_number)
						else:
							if frag_number in repeat_frag_num:
								repeat += 1
							else:
								unique +=1

		except Exception as e:
			print("Error occured in remove_repeat {}".format(e))	
	set_XAs = set(XAs) ## remove duplicates XAs
	fragments = sorted(fragments, key=itemgetter(0)) ## sort fragment's in increasing order
	ctr = 0
	for XA in set_XAs:
		if(XA < 0):
			continue
		for i in range(len(fragments)):
			if(fragments[i][0] > XA):
				break
			if(fragments[i][0] <= XA and XA <= fragments[i][1]):
				ctr+=1
				break

	return ctr,unique,repeat

if __name__ == '__main__':
	pass
