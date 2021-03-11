#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ISEScan version
version = '1.7.2.2'

import argparse
import os
import sys

import isPredict

def isPredictSingle(args):
	seqfile = args['seqfile']
	output = args['output']
	seqfilename = os.path.basename(seqfile)
	org = os.path.basename(os.path.dirname(seqfile))
	filelist = org + '_' + seqfilename + '.list'
	with open(filelist, 'w') as fp:
		fp.write(seqfile+'\n')

	isPredict.isPredict(filelist, args['output'], args['removeShortIS'], args['translateGenome'],
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

	parser.add_argument(
		'--removeShortIS', 
		action='store_true', 
		help = "Remove incomplete (partial) IS elements which include IS element with length < 400 or single copy IS element without perfect TIR.",
		)

	parser.add_argument(
		'--no-FragGeneScan', 
		action='store_false', 
		help = "Use the annotated protein sequences in NCBI GenBank file (.gbk which must be in the same folder with genome sequence file), instead of the protein sequences predicted/translated by FragGeneScan. (Experimental feature!)",
		)

	parser.add_argument(
		'--seqfile',
		required = True,
		default='',
		help = "Sequence file in fasta format, '' by default",
                )

	parser.add_argument(
		'--output',
		required = True,
		default='results',
		help = "Output directory, 'results' by default",
                )

	parser.add_argument(
		'--nthread', 
		required = False, 
		type = int, 
		default = 1, 
		help = 'Number of CPU cores used for FragGeneScan and hmmer, 1 by default.')

	args = parser.parse_args()

	args4isPredictSingle = {
					'removeShortIS' : args.removeShortIS,
					'translateGenome' : args.no_FragGeneScan,
					'seqfile': args.seqfile.strip(),
					'output': args.output.strip(),
					'nthread': args.nthread,
				}

	isPredictSingle(args4isPredictSingle)
