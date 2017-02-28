#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ISEScan version
version = '1.2'

import argparse
import os

import isPredict

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
			python3 isescan.py seqfile proteome hmm'
	parser = argparse.ArgumentParser(prog='isescan', description = descriptStr)

	parser.add_argument('--version', action='version', version='%(prog)s' + ' ' + version)


	helpStr = 'sequence file in fasta format'
	parser.add_argument('seqfile', help = helpStr)

	helpStr = 'directory where proteome (each line corresponds to a protein sequence database translated from a genome) files will be placed'
	parser.add_argument('path2proteome', help = helpStr)

	helpStr = 'directory where the results of hmmsearch will be placed'
	parser.add_argument('path2hmm', help = helpStr)

	args = parser.parse_args()

	args4isPredictSingle = {
					'seqfile': args.seqfile,
					'path2proteome': args.path2proteome,
					'path2hmm': args.path2hmm,
				}

	isPredictSingle(args4isPredictSingle)
