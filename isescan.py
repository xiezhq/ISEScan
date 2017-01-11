#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ISEScan version
version = '1.1'

from sys import argv
import os
from sys import platform as _platform

import isPredict
######## OPERATIONAL SYSTEM IMPROVEMENTS & CORRECTIONS     #
#if _platform == "linux" or _platform == "linux2":	   #
#	## if they was using this on linux		   #
#elif _platform == "darwin":				   #
#	## MAC OS X					   #
#elif _platform == "win32":				   #
#	## WINDOWS					   #
#							   #
############################################################


#def main():
#	print ("Welcome! Here an args list:\n")
#	print ("--version shows the version.\n")
#	print ("seqfile Sequence file in fasta format.\n")
#	print ("path2proteome directory where proteome (each line corresponds to a protein sequence database translated from a genome) files will be placed.\n")
#	print ("path2hmm directory where the results of hmm search will be placed")

#	if sys.argv[1] == "--version": ##action here
	
#	if sys.argv[1] == "seqfile": ## action here
		
#	if sys.argv[1] == "path2proteome": ## action here
		
#	if sys.argv[1] == "path2hmm": ## action here
	
	
def isPredictSingle(args):
	seqfile = args['seqfile']
	path2proteome = args['path2proteome']
	path2hmm = args['path2hmm']
	filelist = seqfile.strip()+'.list'
	with open(filelist, 'w') as fp:
		fp.write(seqfile+'\n')

	isPredict.isPredict(filelist, path2proteome, path2hmm)
	os.remove(filelist)

if __name__ == "__main__":
	# Parse command line arguments
	descriptStr = 'Search IS Profile HMMs against gene database. A typical invocation would be:\n\
			python3 isPredictSingle.py seqfile proteome hmm'


	args4isPredictSingle = {
					'seqfile': args.seqfile,
					'path2proteome': args.path2proteome,
					'path2hmm': args.path2hmm,
				}

	isPredictSingle(args4isPredictSingle)
