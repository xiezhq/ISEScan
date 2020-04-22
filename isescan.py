#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ISEScan version
version = '1.7.2.1'

import argparse
import os
import sys

import isPredict

def isPredictSingle(args):
	seqfile = args['seqfile'].strip()
	path2proteome = args['path2proteome']
	path2hmm = args['path2hmm']
	seqfilename = os.path.basename(seqfile)
	org = os.path.basename(os.path.dirname(seqfile))
	filelist = org + '_' + seqfilename + '.list'
	with open(filelist, 'w') as fp:
		fp.write(seqfile+'\n')

	isPredict.isPredict(filelist, path2proteome, path2hmm, args['removeShortIS'], args['translateGenome'],
			args['nthread'])
	os.remove(filelist)

if __name__ == "__main__":
	import textwrap

	# Parse command line arguments
	descriptStr = '''\
			ISEScan is a python pipeline to identify Insertion Sequence elements (both complete and incomplete IS elements) in genom. A typical invocation would be:
			python3 isescan.py seqfile proteome hmm

			- If you want isescan to report only complete IS elements, you need to set command line option --removeShortIS.'''
	parser = argparse.ArgumentParser(prog='isescan', description = textwrap.dedent(descriptStr), 
			formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('--version', action='version', version='%(prog)s' + ' ' + version)

	helpStr= 'remove incomplete (partial) IS elements which include IS element with length < 400 or single copy IS element without perfect TIR.'
	parser.add_argument('--removeShortIS', action='store_true', help = helpStr)

	helpStr= 'use the annotated protein sequences in NCBI GenBank file (.gbk which must be in the same folder with genome sequence file), instead of the protein sequences predicted/translated by FragGeneScan. (Experimental feature!)'
	parser.add_argument('--no-FragGeneScan', action='store_false', help = helpStr)

	helpStr = 'sequence file in fasta format'
	parser.add_argument('seqfile', help = helpStr)

	helpStr = 'directory where proteome (each line corresponds to a protein sequence database translated from a genome) files will be placed'
	parser.add_argument('path2proteome', help = helpStr)

	helpStr = 'directory where the results of hmmsearch will be placed'
	parser.add_argument('path2hmm', help = helpStr)

	parser.add_argument('--nthread', required = False, type = int, default = 1, 
			help = 'number of CPU cores used for FragGeneScan and hmmer. By default one will be used.')

	args = parser.parse_args()

	args4isPredictSingle = {
					'seqfile': args.seqfile,
					'path2proteome': args.path2proteome,
					'path2hmm': args.path2hmm,
					'removeShortIS' : args.removeShortIS,
					'translateGenome' : args.no_FragGeneScan,
					'nthread': args.nthread,
				}

	isPredictSingle(args4isPredictSingle)
