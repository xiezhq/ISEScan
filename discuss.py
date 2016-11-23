#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import subprocess
import shlex
import os
import datetime
import itertools
import concurrent.futures

import tools

import fastcluster
import numpy
import scipy.spatial.distance
import scipy.cluster.hierarchy
import matplotlib.pyplot as plt

# number of processes to start for computation
nproc = 16
#nproc = 32

#silva = '/home/data/insertion_sequence/silva/SILVA_123.1_SSURef_Nr99_tax_silva_trunc.fasta'
silva = '/home/data/insertion_sequence/silva/SILVA_123.1_SSURef_tax_silva_trunc.fasta'
#silva = '/home/data/insertion_sequence/silva/SILVA_123.1_SSUParc_tax_silva_trunc.fasta'
#ltp = '/home/data/insertion_sequence/living_tree_project/LTPs123_SSU.compressed.fasta'
#silva = ltp

#dir2out4org = '/home/data/insertion_sequence/results/prediction'
#dir4blastout = '/home/data/insertion_sequence/results/blastout'
dir4blastout = '/tmp'
#perc_ident4mapRNA2genome = 99 # identity of alignment between silva 16S rRNA gene and genome
perc_ident4mapRNA2genome = 100 # identity of alignment between silva 16S rRNA gene and genome
coverage2mapRNA2genome = 1.0 # 100%, coverage of RNA gene when map rRNA to genome

# maxlen4taxpath = 6 in taxmap_slv_ssu_ref_123.1.txt.noEukaryota, 
# maxlen4taxpath = 14 in taxmap_slv_ssu_ref_123.1.txt
maxlen4taxpath = 6
# maxlen4taxname = 40 in taxmap_slv_ssu_ref_123.1.txt.noEukaryota
maxlen4taxname = 40

file4rnaMapTax = 'org2dna2rna2tax.list'

# Switch the granularity of IS type between family and familyCluster
IStype2family = True # use family as the minimal granularity of IS type
#IStype2family = False # use familyCluster as the minimal granularity of IS type


# isFile: .out file created by pred.py, listing all ISs in a genome sequence
# islist: [is, ..., is]
# is: {'seqid':seqid, 'family':family, 'familycluster':familyCluster, 'start':start, 'end':end, 'ncopy4is':ncopy4is}
def rdISfile(isFile):
	fp = open(isFile, 'r')
	islist = []
	for line in fp:
		if line[:5] == 'seqID' or line[0] == '#':
			continue
		element = {}
		line = line.strip()
		items = line.split(maxsplit=7)
		element['seqid'] = items[0]
		element['family'] = items[1]
		element['familyCluster'] = items[2]
		element['start'] = int(items[3])
		element['end'] = int(items[4])
		element['ncopy4is'] = int(items[6])
		element['evalue'] = float(line.rsplit(maxsplit=2)[1])
		islist.append(element)
	fp.close()
	return islist

def checkCopies(mislist):
	mncopy = {} # {org: {fileid:ncopy, ...}, ...}
	for org, islist4org in mislist.items():
		for fileid, islist in islist4org.items():
			ncopy = {'copies':0, 'single':0}
			for element in islist:
				if element['ncopy4is'] > 1:
					ncopy['copies'] += 1
				elif element['ncopy4is'] == 1:
					ncopy['single'] += 1
				else:
					e = 'Error: copy number of IS element ({},{}) is 0 in {}'.format(
							element['start'],element['end'], element['seqid'])
					raise RuntimeError(e)
			if org not in mncopy.keys():
				mncopy[org] = {}
			mncopy[org][fileid] = ncopy
	copies, single = 0, 0
	nseq = 0
	norg = 0
	for org in sorted(mncopy.keys()):
		copies4org, single4org = 0, 0
		norg += 1
		ncopy4org = mncopy[org]
		for fileid in sorted(ncopy4org.keys()):
			ncopy = ncopy4org[fileid]
			nseq += 1
			copies += ncopy['copies']
			single += ncopy['single']
			copies4org += ncopy['copies']
			single4org += ncopy['single']
			#print(fileid, ncopy['copies'], ncopy['single'])
		print('{} copies {} single {} copies% {:.2f}'.format(
			org, copies4org, single4org, copies4org/(copies4org+single4org)*100))

	print('Total norg {} nseq {} copies {} single {} copies% {:.2f}'.format(
		norg, nseq, copies, single, copies/(copies+single)*100))

# Return vecArr
# vecArr: [[n, ..., n], ...], len(orgSorted) * len(familySorted)
# n: number of ISs belonging to the specific family in the specific org if boolSwitch=False else 1
#
# family4org: {org:{family:n, ...}, ...}
# n: number of the specific family available in the specific organism
# familysSorted: [family, ...]
def vectorArray(orgSorted, familys4org, familysSorted, boolSwitch=False):
	vecArr = []
	for org in orgSorted:
		vec4org = []
		familys = familys4org[org].keys()
		for name in familysSorted:
			if name in familys:
				if boolSwitch == True:
					vec4org.append(True)
				else:
					vec4org.append(familys4org[org][name])
			else:
				if boolSwitch == True:
					vec4org.append(False)
				else:
					vec4org.append(0)
		vecArr.append(vec4org)
	return vecArr

# u,v: boolean vectors,
#	vectors with 27 states (if using family as dimension of vector), u and v are two different
#	organisms in our project.
# Return dist: 1 or 0, distance between two vectors, u and v.
#	dist = 0 if any common IS (family or familyCluster or IS type) are available in both else 1.
#	u and v are different by defenition if u==v==0, namely, dist(u,v)=1
def singleMatch(u, v):
	dist = 1 # u and v are completely different by default
	for x,y in zip(u, v):
		if x and y: # if both u and v are True
			dist = 0
			break
	return dist

# compute the distance between any two organisms
def distBetweenOrgsByIS(orgSorted, family4org, familySorted, boolSwitch=True):
	orgVect = vectorArray(orgSorted, family4org, familySorted, boolSwitch=boolSwitch)
	print(len(orgSorted), len(familySorted))
	#for item in orgVect:
	#	print(len(orgVect), len(item), item)

	# get the condensed distance matrix (dissimilarity matrix) from contigency table (orgVect)
	if boolSwitch == True:
		X = numpy.array(orgVect, bool)
		#value4metric = 'dice'
		value4metric = singleMatch
	else:
		X = numpy.array(orgVect, int)
		value4metric = 'euclidean'
	print('orgVect: {}\n{}'.format(X.shape, X))
	distMatrix = scipy.spatial.distance.pdist(X, metric=value4metric)
	return distMatrix

def clusterOrgByIS(orgSorted, family4org, familySorted, boolSwitch=True):
	orgVect = vectorArray(orgSorted, family4org, familySorted, boolSwitch=boolSwitch)
	print(len(orgSorted), len(familySorted))
	#for item in orgVect:
	#	print(len(orgVect), len(item), item)

	# get the condensed distance matrix (dissimilarity matrix) from contigency table (orgVect)
	if boolSwitch == True:
		X = numpy.array(orgVect, bool)
		#value4metric = 'dice'
		value4metric = singleMatch
	else:
		X = numpy.array(orgVect, int)
		value4metric = 'euclidean'
	print('orgVect: {}\n{}'.format(X.shape, X))
	distMatrix = scipy.spatial.distance.pdist(X, metric=value4metric)

	#method4dist = 'single' # nearest 
	method4dist = 'average' # UPGMA

	hclusters = fastcluster.linkage(distMatrix, method=method4dist, preserve_input=True)
	#hclusters = fastcluster.linkage(orgVect, method=method4dist, metric=value4metric, preserve_input=False)
	#del distMatrix 

	cophenet = scipy.cluster.hierarchy.cophenet(hclusters, distMatrix)
	print('cophenetCorrelation = {}'.format(cophenet[0]))

	tree = scipy.cluster.hierarchy.to_tree(hclusters, rd=False)
	#newick = tools.linkageTree2newick(tree, orgSorted)
	#fp = open('newick.txt', 'w')
	#fp.write(newick+'\n')
	#fp.close()

	#plt.figure()
	#dn = scipy.cluster.hierarchy.dendrogram(hclusters)
	#dn = scipy.cluster.hierarchy.dendrogram(hclusters, p=30, truncate_mode='lastp')
	#dn = scipy.cluster.hierarchy.dendrogram(hclusters, p=30, truncate_mode='level')
	#plt.show()

	return (distMatrix, tree)

def outputTaxByOrg(file, orgs):
	fp = open(file, 'w')
	for org in sorted(orgs.keys()):
		taxp = orgs[org]
		fp.write('{:<85}'.format(org))
		for p in taxp:
			fp.write(' {:>10}'.format(p))
		fp.write('\n')
	fp.close()

def outputAccid(file, rnaMap):
	fp = open(file, 'w')
	print('# organism qseqid qstart qend sseqid sstart send slen perc_scoverage perc_ident length', file=fp)
	n = 0
	for org in sorted(rnaMap.keys()):
		rnas = rnaMap[org]
		nseqWithRNA = len(rnas)
		if nseqWithRNA > 1:
			n += 1
			print('multiple genome sequences have 16S rRNA gene in {}: {}'.format(org, rnas.keys()))

		for fileid,rna in rnas.items():
			#print('{:<85} {:<10} {:>8} {:>8} {:<40} {:>8} {:>8} {:>4} {:>6.2f} {:>6.2f} {:>4}'.format(
			print('{:<85} {:<50} {:>8} {:>8} {:<40} {:>8} {:>8} {:>4} {:>6.2f} {:>6.2f} {:>4}'.format(
				org, 
				#fileid, rna['qstart'], rna['qend'],
				rna['qseqid'], rna['qstart'], rna['qend'],
				rna['sseqid'], rna['sstart'], rna['send'], rna['slen'], (rna['length']/rna['slen'])*100,
				rna['pident'], rna['length'],
				), file=fp)
	fp.close()
	print('number of organisms with 16S rRNA gene in multiple genome sequences:', n)

# return rnaMap
# rnaMap: {org:rna4seq, ...}
# rna4seq: {fileid:rna, ...}, rna for each sequence in the specific organism
# rna: {'qseqid':qseqid, 'sseqid':sseqid, 'pident':pident, 'length':length, 'qlen':qlen, 'slen':slen}
# qseqid: query sequence ID, eg. 'gi|556503834|ref|NC_000913.3|'
# sseqid: subject sequence ID, eg. 'CP009685.3470134.3471689', 'CP000946.4686640.4688173'
rnaMap = {}
def id4Genome2rna(file4rnaAccid):
	rnaMap = {}
	with open(file4rnaAccid) as fp:
		for line in fp:
			if line[0] == '#':
				continue
			items = line.split()
			org, qseqid, qstart, qend, sseqid, sstart, send, slen, perc_scoverage, perc_ident, length = items
			qlen = 0
			rna = {'qseqid':qseqid, 'sseqid':sseqid, 'pident':float(perc_ident), 
					'length':int(length), 'qlen':qlen, 'slen':int(slen)}
			if org not in rnaMap.keys():
				rnaMap[org] = {}
			rnaMap[org][qseqid] = rna
	return rnaMap

def doBlast(args):
	query, blastdb, blastout, strand, task, perc_ident = args
	#tools.doBlastn(query, blastdb, blastout, strand, task, perc_ident) 
	min4coverage = coverage2mapRNA2genome
	hitpairs = tools.getBlastout(blastout, min4coverage)
	return hitpairs

"""
Do multiple sequence alignment and output MSA alignment in Phylip format.
fastafile: fasta format file holding multiple sequences

clustalo -i SILVA_123.1_SSURef_tax_silva_trunc.fasta.ise -t RNA \
		--distmat-out=SILVA_123.1_SSURef_tax_silva_trunc.fasta.ise.distMatrix \
		--outfmt=phy -o SILVA_123.1_SSURef_tax_silva_trunc.fasta.ise.phy \
		--threads=32 --full
"""
def doMsa(fastaFile, distFile, algnFile):
	cmd = '/u/zhiqxie/informatics/inst/clustal-omega-1.2.1/bin/clustalo'
	cmdline = [cmd, '-i', fastaFile, '-t', 'RNA', '--distmat-out=distFile', '-o', algnFile, '--outfmt=phy', '--threads=32', '--full']
	callcmd = shlex.split(' '.join(cmdline))
	subprocess.check_output(callcmd, shell=False, universal_newlines=False, stderr=subprocess.STDOUT)


# map genome identifiers to 16S rRNA identifiers by blasting genome sequence against 16S rRNA sequence database
def genome2rna(dnaFiles, silva, dir4blastout, perc_ident4mapRNA2genome, coverage2mapRNA2genome):
	# map SILVA 16S rRNA gene to genome sequence to get 16S rRNA gene in genome sequence, 
	# #and then calculate distance matrix and tree (phynogenetic tree) based on 16S rRNA sequences 
	# #from genome sequences
	# dnaFiles: [(file, org), ..., (file, org)]
	# rnaMap: {org:rna4seq, ...}
	# rna4seq: {fileid:rna, ...}, rna for each sequence in the specific organism
	# rna: {'qseqid':qseqid, 'sseqid':sseqid, 'pident':pident, 'length':length, 'qlen':qlen, 'slen':slen}
	# qseqid: query sequence ID, eg. 'gi|556503834|ref|NC_000913.3|'
	# sseqid: subject sequence ID, eg. 'CP009685.3470134.3471689', 'CP000946.4686640.4688173'
	rnaMap = {}

	'''
	for file_org in dnaFiles:
		file, org = file_org
		fileid = os.path.basename(file).split('.', maxsplit=1)[0]
		query = file
		blastdb = silva
		blastout = os.path.join(dir4blastout, org, fileid+'.'+os.path.basename(blastdb)+'.out')
		tools.makedir(os.path.dirname(blastout))
		perc_ident = perc_ident4mapRNA2genome
		tools.doBlastn(query, blastdb, blastout, strand='both', task='megablast', perc_ident) 
		min4coverage = coverage2mapRNA2genome
		hitpairs = tools.getBlastout(blastout, min4coverage)
		if len(hitpairs) == 0:
			print('No exact 16S rRNA gene (with same length and 100% identity) for', org, fileid) 
		else:
			hitpairs.sort(key = lambda x: x['length'], reverse = True)
			hitpairs.sort(key = lambda x: x['pident'], reverse = True)
			if org not in rnaMap.keys():
				rnaMap[org] = {}
			rnaMap[org][fileid] = hitpairs[0]
	'''
	args4concurrent = []
	for file_org in dnaFiles:
		file, org = file_org
		fileid = os.path.basename(file).split('.', maxsplit=1)[0]
		query = file
		blastdb = silva
		blastout = os.path.join(dir4blastout, org, fileid+'.'+os.path.basename(blastdb)+'.out')
		tools.makedir(os.path.dirname(blastout))
		strand = 'both'
		task = 'megablast'
		perc_ident = perc_ident4mapRNA2genome
		args4concurrent.append((query, blastdb, blastout, strand, task, perc_ident))
	if len(args4concurrent) > nproc:
		nprocess = nproc
	else:
		nprocess = len(args4concurrent)
	with concurrent.futures.ProcessPoolExecutor(max_workers = nprocess) as executor:
		for args, hitpairs in zip(args4concurrent, executor.map(doBlast, args4concurrent)):
			file = args[0] # query sequence
			org = os.path.dirname(file).rsplit('/', maxsplit=1)[-1]
			fileid = os.path.basename(file).split('.', maxsplit=1)[0]
			if len(hitpairs) == 0:
				print('No exact 16S rRNA gene (with same length and 100% identity) for', org, fileid) 
			else:
				hitpairs.sort(key = lambda x: x['length'], reverse = True)
				hitpairs.sort(key = lambda x: x['pident'], reverse = True)
				if org not in rnaMap.keys():
					rnaMap[org] = {}
				rnaMap[org][fileid] = hitpairs[0]
	return rnaMap


# seqMeta: {'dnaType':type, 'dnaLen':dnaLen}, type can be 0,1,2 to represent phase or plasmid or chromosome DNA, 
#	 respectively
def getMeta(fnafile):
	seqMeta = {}
	headlines = []
	with open(fnafile, 'r') as fp:
		for line in fp:
			line = line.strip()
			if line[0] != '>':
				continue
			headlines.append(line)
	for line in headlines:
		if 'phage' in line.lower():
			seqMeta['dnaType'] = 0
		elif 'plasmid' in line.lower():
			seqMeta['dnaType'] = 1
		else:
			seqMeta['dnaType'] = 2
	return seqMeta

# tax: {seqid:taxorg, ...}
# taxorg: {'taxpath':taxpath, 'org':org, 'taxid':taxid}
# taxpath: [rank1, rank2, ...], e.g.: 
#	['Bacteria', 'Spirochaetae', 'Spirochaetes', 'Spirochaetales', 'Spirochaetaceae', 'Treponema 2']
#	['Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Rickettsiales', 'Mitochondria']
#	['Bacteria', 'Cyanobacteria', 'Chloroplast']
#	['Archaea', 'Thaumarchaeota', 'Marine Group I']
# org: orgnism name, e.g. 'Treponema berlinense', 'Chlorella vulgaris', 'Puccinia triticina'
def getSilvaTax(silvaTax):
	tax = {}
	with open(silvaTax) as fp:
		for line in fp:
			line = line.strip()
			if line[:16] == 'primaryAccession':
				continue
			items = line.split('\t')
			seqid = '.'.join(items[:3])
			taxpath = items[3].strip(' ;').split(';')
			taxorg = {'taxpath':taxpath, 'org':items[4], 'taxid':int(items[5])}
			tax[seqid] = taxorg
	return tax

# check lengths of taxonomic paths and the lengths of taxon names
def checkTaxLen(tax):
	taxpathlist = [tax4seq['taxpath'] for tax4seq in tax.values()]
	taxpathlist.sort(key = lambda x: len(x))
	print('shortest taxonomy path: {} {}, longest taxonomy path: {} {}'.format(
		len(taxpathlist[0]), taxpathlist[0], len(taxpathlist[-1]), taxpathlist[-1]))
	taxnamelist = []
	for tax4seq in tax.values():
		taxnamelist.extend(tax4seq['taxpath']) 
	taxnamelist = list(set(taxnamelist))
	taxnamelist.sort(key = lambda x: len(x))
	print('longest length of taxnomy name in taxonomy path:', len(taxnamelist[-1]), taxnamelist[-1])

# Update rnaMap with attaching taxonomy info, rnaMap is updated as below
# rnaMapTax: {org:rnas, ...}
# rnas: {fileid:dna2rna, ...}, mapping DNA to rna for each DNA sequence in the specific organism
# dna2rna: {'tax4seq':tax4seq, 'dnaType':type, 'qseqid':qseqid, 'sseqid':sseqid, ...}
# tax4seq: {'taxpath':taxpath, 'org':org, 'taxid':taxid}
# taxpath: [rank1, rank2, ...], e.g. ['Archaea', 'Thaumarchaeota', 'Marine Group I']
def org2tax(rnaMap, tax):
	rnaMapTax = {}
	for org,rnas in rnaMap.items():
		rnasTax = {}
		for fileid,dna2rna in rnas.items():
			if dna2rna['sseqid'] in tax.keys():
				tax4seq = tax[dna2rna['sseqid']]
				dna2rna['tax4seq'] = tax4seq
				rnasTax[fileid] = dna2rna
			else:
				e = 'No taxon (taxonomic path) for rna gene {} in {}'.format(dna2rna['sseqid'], org)
				raise RuntimeError(e)
		if len(rnasTax) > 0:
			rnaMapTax[org] = rnasTax
	return rnaMapTax

# check if multiple different taxons are assigned (mapped) to one organism.
# Return updated rnaMap including only organisms with only one taxon assigned to one or more chromosome DNA
# in the organism.
def checkOrg4multipleTax(rnaMap):
	# rnaMapNew: rnaMap including only organisms with only one taxon assigned to one or more chromosome DNA
	rnaMapNew = {}
	for org,rnas in rnaMap.items():
		ntax4chrom = 0
		taxids = []
		for dna2rna in rnas.values():
			if 'Eukaryota' in dna2rna['tax4seq']['taxpath']:
				print('Eukaryota found:', org, dna2rna['dnaType'], dna2rna['qseqid'])
			if dna2rna['dnaType'] != 2:
				continue
			# get taxid of chromosome DNA
			taxid = dna2rna['tax4seq']['taxid']
			if taxid not in taxids:
				taxids.append(taxid)
		ntaxids = len(taxids)
		if ntaxids < 1:
			print('no taxon is assigned to chromosome DNA in', org)
		elif ntaxids == 1:
			rnaMapNew[org] = rnas
		else:
			print('{} different taxons are assigned to chromosome DNAs in {}'.format(ntaxids, org))
			for fileid,dna2rna in rnas.items():
				print(org, fileid, dna2rna['dnaType'], dna2rna['tax4seq']['taxpath'])
	return rnaMapNew

# rnaMap: {org:rnas, ...}
# rnas: {fileid:dna2rna, ...}, mapping DNA to rna for each DNA sequence in the specific organism
# dna2rna: {'tax4seq':tax4seq, 'dnaType':type, 'qseqid':qseqid, 'sseqid':sseqid, 'pident':pident, 'length':length, 'qlen':qlen, 'slen':slen}
# tax4seq: {'taxpath':taxpath, 'org':org, 'taxid':taxid}
# taxpath: [rank1, rank2, ...], e.g. ['Archaea', 'Thaumarchaeota', 'Marine Group I']
def outputOrg2dna2rna2tax(outfile, rnaMap):
	#'organism qseqid qstart qend sseqid sstart send slen perc_scoverage perc_ident length dna taxid taxpath'
	fmtstr4org2dna2rna = '{:<85} {:<10} {:<40} {:>4} {:>6.2f} {:>6.2f} {:>4}' + ' {:>1} {:>6}'
	fmtstr4taxpath = ' '.join(['{:<'+str(maxlen4taxname)+'}'] * maxlen4taxpath)
	fmtstr = fmtstr4org2dna2rna + ' ' + fmtstr4taxpath
	with open(outfile, 'w') as fp:
		for org in sorted(rnaMap.keys()):
			rnas = rnaMap[org]
			for fileid,dna2rna in rnas.items():
				taxpath = dna2rna['tax4seq']['taxpath']
				fillnames = ['' for i in range(maxlen4taxpath - len(taxpath))]
				taxnames = []
				taxnames.extend(taxpath)
				taxnames.extend(fillnames)
				line = fmtstr.format(
						org, fileid, dna2rna['sseqid'], dna2rna['slen'], 
						(dna2rna['length']/dna2rna['slen'])*100,
						dna2rna['pident'], dna2rna['length'],
						dna2rna['dnaType'], dna2rna['tax4seq']['taxid'],
						taxnames[0], taxnames[1], taxnames[2], taxnames[3], 
						taxnames[4], taxnames[5],
						)
				fp.write(line + '\n')

# add taxonomy to the content of the old file and then write the updated content to the new file
def tax2issum(sumfileOld, sumfileNew, rnaMap):
	fmtstr = ' ' + ' '.join(['{:<'+str(maxlen4taxname)+'}'] * maxlen4taxpath) + '\n'
	with open(sumfileOld, 'r') as infp, open(sumfileNew, 'w') as outfp:
		for line in infp:
			if line[:8] == 'organism' or line[:5] == 'total':
				taxnames = ['rank'+str(i) for i in range(maxlen4taxpath)]
			else:
				org = line.split(maxsplit=1)[0]
				# organism with taxonomy info
				if org in rnaMap.keys():
					# it must be guaranteed that only one sequence in rnaMap[org], eg one chromosome DNA
					dna2rna = next(iter(rnaMap[org].values()))
					taxpath = dna2rna['tax4seq']['taxpath']
				# organism without taxonomy info
				else:
					taxpath = []
				fillnames = ['' for i in range(maxlen4taxpath - len(taxpath))]
				taxnames = []
				taxnames.extend(taxpath)
				taxnames.extend(fillnames)
			taxnamesStr = fmtstr.format(taxnames[0], taxnames[1], taxnames[2], taxnames[3], 
					taxnames[4], taxnames[5])
			outfp.write(line.replace('\n', taxnamesStr))

# Return org2istax
# org2istax: {org:istax, ...}
# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
# isfamily: {familyname:(nIS,norg), ...}
# taxpath: [rank0, rank1, ...]
def getsumtax(sumtaxfile):
	org2istax = {}
	with open(sumtaxfile, 'r') as fp:
		for line in fp:
			line = line.strip()
			if line[:5] == 'total':
				continue
			# items: organism, nIS, %genome, bps4is, dnaLen4is, dnaLen,
			#	ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage,
			#	family1, family1_%, family1_bps, family1_s,
			#	...,
			#	family27, family27_%, family27_bps, family27_s,
			#	rank0, rank1, ..., rank5
			# Note: in is.sum.tax file, we can have 'ISNCY_s rank0' content in each line, 
			#	the last letter of 'ISNCY_s' ('s') is the 1707th column of the line.
			#	We need to be careful of the taxpath content in each line because some
			#	ranknames from silva is a string with space, for example, 'Acidobacteriaceae (Subgroup 1)'.
			#	We therefore process the content before 1707th column and the content after that column
			#	in each line.
			items = line[:1707].split()
			if items[0] == 'organism':
				familynames = items[12:120:4]
				continue
			org, nIS, bps4is, dnaLen4is, dnaLen = items[0], int(items[1]), int(items[3]), int(items[4]), int(items[5])
			if org not in org2istax.keys():
				org2istax[org] = {}
			org2istax[org]['nIS'] = nIS
			#org2istax[org]['bps4is'] = bps4is
			org2istax[org]['dnaLen4is'] = dnaLen4is
			org2istax[org]['dnaLen'] = dnaLen
			# process items[12:12+108], namely, family1 through family27
			# Note: number of items for 27 families is 108 = 4*27
			org2istax[org]['is'] = {}
			for familyname,family_nIS,family_norg in zip(familynames, items[12:120:4], items[15:120:4]):
				if int(family_nIS) > 0:
					org2istax[org]['is'][familyname] = [int(family_nIS), int(family_norg)]
			# process the taxpath content in  range of line[1708:] in each line
			# Note: in is.sum.tax file, taxpath content starts at the 1709th column in each line, namely, 
			# line[1708:].  And, there are 6 blocks corresponding 6 taxonomy ranks in each line for taxpath, 
			# namely, rank0-rank5 each of which occupies 41 (40+1:rank_name + space_between_two_neighboring_rnak_names)
			# columns. We need process each rank block separately as there are some special rank names 
			# which have space in the rank name, e.g. 'Acidobacteriaceae (Subgroup 1)'.
			# In addition, org2istax[org]['taxpath'] may be [] when an organism is not mapped by SILVA rna gene.
			str4taxpath = line[1708:].strip()
			taxpath = [str4taxpath[i:i+41].strip() for i in range(0, len(str4taxpath), 41)]
			org2istax[org]['taxpath'] = taxpath
			
	return org2istax


# mislist: {org:{fileid:islist, ...}, ...}
# islist: [is, ..., is]
# is: {'seqid':seqid, 'family':family, 'familyCluster':familyCluster, 
#		'start':start, 'end':end, 'ncopy4is':ncopy4is}
#
# seqmeta4org: {org:meta4org, ...}
# meta4org: {fileid:metainfo, ...}
# metainfo: {'dnaType':dnaType, 'dnaLen':dnaLen}
#
# sum4org: {org:s, ...}
# s: {'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is, 
#		'family':family4org}
# family4org: {familyname:(nIS,norg), ...}
# sum4seq: {'nIS':nIS, 'family':family, 'metainfo':metainfo}
# family: {familyname:nIS, ...}
def getsum4org(mistlist, seqmeta4org):
	sum4org = {}
	for org,meta4org in seqmeta4org.items():
		sum4org[org] = {}
		sum4org[org]['dnaLen'] = sum([metainfo['dnaLen'] for metainfo in meta4org.values()])
		sum4org[org]['dnaLen4is'] = 0
		sum4org[org]['nIS'] = 0
		sum4org[org]['family'] = {}
		if org not in mistlist.keys():
			continue

		seqs = {}
		for fileid,islist in mistlist[org].items():
			sum4seq = {}
			family = {}
			for x in islist:
				if IStype2family == True:
					familyname = x['family']
				else:
					familyname = x['familyCluster']
				if familyname not in family.keys():
					family[familyname] = 0
				family[familyname] += 1
			sum4seq['family'] = family
			sum4seq['metainfo'] = meta4org[fileid]
			seqs[fileid] = sum4seq
		family4org = {} # for org
		for fileid,sum4seq in seqs.items():
			for familyname,nIS in sum4seq['family'].items():
				if familyname not in family4org.keys():
					family4org[familyname] = [0,1]
				family4org[familyname][0] += nIS

			if sum([v[0] for v in family4org.values()]) > 0:
				sum4org[org]['dnaLen4is'] += sum4seq['metainfo']['dnaLen']
		#sum4org[org]['seqs'] = seqs
		sum4org[org]['family'] = family4org
		sum4org[org]['nIS'] = sum([v[0] for v in family4org.values()])
	return sum4org

def getseqmeta(dnaFiles):
	seqmeta4org = {}
	for x in dnaFiles:
		fnafile, org = x
		fileid = os.path.basename(fnafile).split('.', 1)[0]
		seqs = tools.getFastaFull(fnafile)
		# seqs: [(header, seq), ..., (header, seq)]
		if org not in seqmeta4org.keys():
			seqmeta4org[org] = {}
		for seqi in seqs:
			# header: >header, namely, header line with '>' removed
			header, seq = seqi
			seqid = header.split(maxsplit=1)[0]
			metainfo = {}
			metainfo['dnaLen'] = len(seq)
			if 'phage' in header.lower():
				metainfo['dnaType'] = 0
			elif 'plasmid' in header.lower():
				metainfo['dnaType'] = 1
			else:
				metainfo['dnaType'] = 2
			#seqmeta4org[org][seqid] = metainfo
		
		if len(seqs) > 0:
			# use the metainfo of the last seq as the metainfo of fnafile
			seqmeta4org[org][fileid] = metainfo
		else:
			seqmeta4org[org][fileid] = {}
	return seqmeta4org

# Return org2istax
# org2istax: {org:istax, ...}
# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
# isfamily: {familyname:(nIS,norg), ...}
# taxpath: [rank0, rank1, ...]
#
# rnaMap: must have only one taxonomy clade assigned, e.g. clade assigned to one chromosome DNA in organism.
def getsumtax4org(sum4org, rnaMap):
	org2istax = {}
	'''
	familynames = set()
	for s in sum4org.values():
		familynames.update(set(s['family'].keys()))
	'''
	for org,s in sum4org.items():
		if org not in org2istax.keys():
			org2istax[org] = {}
		org2istax[org]['dnaLen'] = s['dnaLen']
		org2istax[org]['dnaLen4is'] = s['dnaLen4is']
		org2istax[org]['nIS'] = s['nIS']
		org2istax[org]['is'] = s['family']

		'''
		# fill the the left familys into org2istax[org]['is'] in order that it contains all familys
		names2fill = familynames - set(s['family'].keys())
		for name in names2fill:
			org2istax[org]['is'][name] = [0,0]
		'''

		if org in rnaMap.keys():
			# it must be guaranteed that only one taxonomy clade is assigned to org, 
			# e.g. clade assigned to one chromosome DNA in organism.
			dna2rna = next(iter(rnaMap[org].values()))
			org2istax[org]['taxpath'] = dna2rna['tax4seq']['taxpath']
		else:
			org2istax[org]['taxpath'] = []
	return org2istax

# construct string for each genome size bin against the number of ISs per genome
def binlines(ks, avgs, stds, ng):
	lines = []
	for k,avg,std,n in zip(ks, avgs, stds, ng):
		line = '{:>8} {:>10.2f} {:>10.2f} {:>4}'.format(k, avg, std, n)
		lines.append(line)
	return lines

# calculate the distributions of number of ISs per genome against genome size bin
def isVSgenomeSizeBin(org2istax):
	# istaxlist: [(org,istax), ...]
	# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen}
	istaxlist = list(org2istax.items())
	istaxlist.sort(key = lambda x: x[1]['dnaLen'])

	#step = 1 # 1 bp
	#step = 5000 # 5 kb
	#step = 10000 # 10 kb
	#step = 100000 # 100 kb
	#step = 500000 # 500 kb
	step = 1000000 # 1 Mb
	ks = []
	bins = []
	if step > 1:
		keyfunc = lambda x: x[1]['dnaLen']//step
	else:
		keyfunc = lambda x: x[0]
	for k,g in itertools.groupby(istaxlist, key=keyfunc):
		g = list(g)
		if step > 1:
			ks.append(k)
		else:
			ks.append(g[0][1]['dnaLen'])
		bins.append(g)
	ng = []
	for g in bins:
		ng.append(len(g))

	# average of number of ISs per genome (within the specific genome size bin) vs bin of genome size,
	# the genome size bins are like 0-0.5Mb, 0.5-1.0Mb, 1.0-1.5Mb, ...
	avgs = []
	stds = []
	median = []
	for g,n in zip(bins,ng):
		'''
		nIS = [x[1]['nIS'] for x in g]
		avgs.append(numpy.mean(nIS))
		stds.append(numpy.std(nIS))
		median.append(numpy.median(nIS))
		'''
		bps4is = [x[1]['bps4is'] for x in g]
		avgs.append(numpy.mean(bps4is))
		stds.append(numpy.std(bps4is))
		median.append(numpy.median(bps4is))

	headline = '{:>8} {:>10} {:>10} {:>4}'.format('bin', 'mean', 'stdev', 'ng')
	lines = binlines(ks, avgs, stds, ng)
	with open('bin2avg.'+str(step), 'w') as fp:
		fp.write(headline+'\n'+'\n'.join(lines)+'\n')

	# median value of number of ISs per genome (within the specific genome size bin) vs bin of genome size,
	# the genome size bins are like 0-0.5Mb, 0.5-1.0Mb, 1.0-1.5Mb, ...
	'''
	for g in bins:
		q, r = divmod(len(g), 2)
		if r != 0:
			m = g[q][1]['nIS']
		# When the number of data points is even, the median is interpolated 
		# by taking the average of the two middle values:
		else:
			m = (g[q-1][1]['nIS']+g[q][1]['nIS'])/2
	'''
	lines = []
	headline = '{:>8} {:>10} {:>10} {:>4}'.format('bin', 'median', 'stdev', 'ng')
	lines = binlines(ks, median, stds, ng)
	with open('bin2median.'+str(step), 'w') as fp:
		fp.write(headline+'\n'+'\n'.join(lines)+'\n')

	# average of number of ISs per genome (within the specific genome size bin) vs bin of genome size,
	# the genome size bins are like 0-0.5Mb, 0.5-1.0Mb, 1.0-2.0Mb, 2.0-4.0Mb, 4.0-8.0Mb, ...

	# median value of number of ISs per genome (within the specific genome size bin) vs bin of genome size,
	# the genome size bins are like 0-0.5Mb, 0.5-1.0Mb, 1.0-2.0Mb, 2.0-4.0Mb, 4.0-8.0Mb, ...


# org2istax: {org:istax, ...}
# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
# isfamily: {familyname:(nIS,norg), ...}
# taxpath: [rank0, rank1, ...]
#
# rank: int, the first rank taxonomic levels employed in taxon-based analysis,
#	for example, rank=2 means classified genomes based on the rank0;rank1 
#	(e.g. Bacteria;Cyanobacteria, Bacteria;Proteobacteria, ...)
#	to do summarization.
def analysis2org2istax(org2istax, rank, ofile):
	# taxsum: {taxon:sum, ...}
	# 	{rank0:sum, ...} means {'Bacteria':sum, 'Archaea':sum}
	#	{(rank0,rank1):sum, ...} means {('Bacteria','Cyanobacteria'):sum, ('Bacteria,Proteobacteria'):sum, ..., 
	#		('Archaea,Crenarchaeota'):sum, ('Archaea,Euryarchaeota'):sum, ...}
	# sum: {'norg':norg, 'norgis':norgis, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is,
	#		'orgs':orgs, 'orgs4is':orgs4is, 'familys':familys}
	# norg: number of organisms assigned to the specific taxon
	# norgis: number of organisms (with IS identified) assigned to the specific taxon
	# nIS: number of ISs occuring
	# dnaLen: total length of all DNAs in an organism, namely, size of a genome
	# dnaLen4is: total length of all DNAs with IS in an organism
	# orgs: [org, ...], names of organisms
	# orgs4is: [org, ...], names of organisms with IS identified
	# familys: {family1:[nIS, norg], ...}
	# norg: number of organisms containing the specific IS family
	taxsum = {}
	for org,istax in org2istax.items():
		isfamily = istax['is']
		taxpath = istax['taxpath']

		# only consider the first rank, namely, domain e.g. Bacteria and Archaea
		#taxon = tuple(taxpath[:rank])

		# only consider the first two ranks (rank0,rank1), 
		# namely, phylum e.g. Bacteria;Acetothermia and Archaea;Crenarchaeota
		# Note: taxon = () when taxpath = []
		taxon = tuple(taxpath[:rank])

		if taxon not in taxsum.keys():
			# initialization for each taxon
			taxsum[taxon] = {'norg':0, 'norgis':0, 'nIS':0, 'dnaLen':0, 'dnaLen4is':0,
					'orgs':[], 'orgs4is':[], 'familys':{}}
		taxsum[taxon]['nIS'] += istax['nIS']
		taxsum[taxon]['norg'] += 1
		taxsum[taxon]['orgs'].append(org)
		if istax['nIS'] > 0:
			taxsum[taxon]['norgis'] += 1
			taxsum[taxon]['orgs4is'].append(org)
		taxsum[taxon]['dnaLen'] += istax['dnaLen']
		taxsum[taxon]['dnaLen4is'] += istax['dnaLen4is']
		for family in isfamily.keys():
			# initilization for each family in the specific taxon
			if family not in taxsum[taxon]['familys'].keys():
				taxsum[taxon]['familys'][family] = [0,0]
			taxsum[taxon]['familys'][family][0] += isfamily[family][0]
			taxsum[taxon]['familys'][family][1] += isfamily[family][1]
	# Header line in output file: taxon nOrg nIS dnaLen familys
	# familys: family1_nIS family1_norg ... family27_nIS family27_norg
	#fmt4str1 = '{:<9} {:>4} {:>6} {:>12}'
	# Header line in output file: nOrg nIS dnaLen familys taxon
	# familys: family1_nIS family1_norg ... family27_nIS family27_norg
	fmt4str1 = '{:>4} {:>6} {:>6} {:>12} {:>12}'
	# max length of familyCluster name is 64, so we need here 64+3 and 64+4 for familyCluster_nIS
	#	and familyCluster_nOrg, respectively
	#fmt4family = '{:>18} {:>19}' # for family
	fmt4family = '{:>67} {:>68}' # for familyCluster
	fmt4taxon = '{:<}'
	familynames = set()
	for s in taxsum.values():
		familynames.update(set(s['familys'].keys()))
	familynamesSorted = sorted(familynames)
	with open(ofile, 'w') as fp:
		lines = []
		#str1 = fmt4str1.format('taxon', 'nOrg', 'nIS', 'dnaLen')
		#str1 = fmt4str1.format('nOrg', 'nOrgIS', 'nIS', 'dnaLen')
		str1 = fmt4str1.format('nOrg', 'nOrgIS', 'nIS', 'dnaLen', 'dnaLen4is')
		str4familys = []
		for family in familynamesSorted:
			str4family = fmt4family.format(family+'_nIS', family+'_norg')
			str4familys.append(str4family)
		#line = str1 + ' ' + ' '.join(str4familys) + '\n'
		str4taxon = fmt4taxon.format('taxon')
		line = str1 + ' ' + ' '.join(str4familys) + ' ' + str4taxon + '\n'
		lines.append(line)
		for taxon in sorted(taxsum.keys()):
			# taxon nIS familys
			# familys: family1_nIS family1_norg family2_nIS family2_norg ...
			str4taxpath = ';'.join(taxon)
			#str1 = fmt4str1.format(str4taxpath, taxsum[taxon]['norg'], 
			str1 = fmt4str1.format(taxsum[taxon]['norg'], taxsum[taxon]['norgis'],
					taxsum[taxon]['nIS'], taxsum[taxon]['dnaLen'], taxsum[taxon]['dnaLen4is'])
			str4familys = []
			for family in familynamesSorted:
				if family not in taxsum[taxon]['familys'].keys():
					nISnorg = (0, 0)
				else:
					nISnorg = taxsum[taxon]['familys'][family]
				str4family = fmt4family.format(nISnorg[0], nISnorg[1])
				str4familys.append(str4family)
			#line = str1 + ' ' + ' '.join(str4familys) + '\n'
			str4taxon = fmt4taxon.format(str4taxpath)
			line = str1 + ' ' + ' '.join(str4familys) + ' ' + str4taxon + '\n'
			lines.append(line)
		fp.write(''.join(lines))
	return taxsum

# Find potential organism pairs involved in HGT events where two organisms have at least one common
# IS type but are assigned to the different taxonomy clades.
def findOrgpairsInDiffTax(ijorgnames, taxsum):
	orgs4difftax = []
	for orgs in ijorgnames:
		sametax = False
		for k,v in taxsum.items():
			if orgs[0] in v['orgs4is'] and orgs[1] in v['orgs4is']:
				sametax = True
				break
		if sametax == False:
			orgs4difftax.append(orgs)
	return orgs4difftax

# Summarize IS and taxon for each IS family
# org2istax: {org:istax, ...}
# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen}
# isfamily: {familyname:(nIS,norg), ...}
# taxpath: [rank0, rank1, ...]
def getsum4istax(org2istax, rank, ofile):
	# sum4istax: {family:sum4family, ...}
	# sum4family: {'nIS':nIS, 'nCluster':nCluster, 'nOrg':nOrg, 
	#		'taxon': {taxon, ...},
	#		}
	# taxon: set, unique taxpath with the specific IS family
	sum4istax = {}
	for istax in org2istax.values():
		# Do not do summarization for organisms without IS
		if istax['nIS'] == 0:
			continue
		# Do not count organisms without taxonomy assignment.
		#if len(istax['taxpath']) == 0:
		#	continue
		# If only consider the first two ranks (rank0,rank1), 
		# namely, phylum e.g. Bacteria;Acetothermia and Archaea;Crenarchaeota
		# Note: taxon = () when taxpath = []
		taxpath = istax['taxpath']
		taxon = tuple(taxpath[:rank])
		for familyname,family in istax['is'].items():
			# nIS of the specific family in the specific organism
			if istax['is'][familyname][0] == 0:
				continue
			if familyname not in sum4istax.keys():
				sum4istax[familyname] = {}
				sum4istax[familyname]['nIS'] = 0
				sum4istax[familyname]['nOrg'] = 0
				sum4istax[familyname]['taxon'] = set()
			sum4istax[familyname]['nIS'] += istax['is'][familyname][0]
			sum4istax[familyname]['nOrg'] += istax['is'][familyname][1]
			sum4istax[familyname]['taxon'].add(taxon)

	# fmtstr: familyname, nIS, nOrg, nTax, nIS/nOrg
	#fmtstr4title = '{:<11} {:>6} {:>4} {:>4} {:>8}'
	#fmtstr = '{:<11} {:>6} {:>4} {:>4} {:>8.2f}'
	# The longest familyCluster is 64, then we need here 64 istead of 11.
	fmtstr4title = '{:<64} {:>6} {:>4} {:>4} {:>8}'
	fmtstr = '{:<64} {:>6} {:>4} {:>4} {:>8.2f}'
	rows = []
	for family in sorted(sum4istax.keys()):
		sum4family = sum4istax[family]
		row = fmtstr.format(family, 
				sum4family['nIS'], sum4family['nOrg'], len(sum4family['taxon']),
				sum4family['nIS']/sum4family['nOrg'])
		rows.append(row)
	lines = '\n'.join(rows)
		
	headline = fmtstr4title.format('family', 'nIS', 'nOrg', 'nTax', 'nIS/nOrg')
	with open(ofile, 'w') as fp:
		fp.write(headline + '\n' + lines + '\n')
	return sum4istax

# Do analysis based on IS element prediction in prokaryotic genomes
def discuss(args):
	dnaListFile4orgs = args['fileList']
	print('Discussion begins at', datetime.datetime.now().ctime())

	# Get the list of fasta files to process
	#
	# dnaFiles: [(file, org), ..., (file, org)]
	dnaFiles = tools.rdDNAlist(dnaListFile4orgs)

	# Get IS list from fileid.out for each fasta file, organisms without IS are not included.
	#
	# mislist: {org: {fileid:islist, ...}, ...}
	# islist: [is, ..., is]
	# is: {'seqid':seqid, 'family':family, 'familyCluster':familyCluster, 
	#		'start':start, 'end':end, 'ncopy4is':ncopy4is}
	# fullOrg: {org:fileids, ...}, full list of organisms from input file list, e.g. fna.list
	# fileids: [fileid, ...], all file IDs belonging to one organism (species)
	mislist = {}
	fullOrg = {}
	dir2out4org = args['dir2prediction']
	for fileOrg in dnaFiles:
		file, org = fileOrg
		fileName = os.path.basename(file)
		#fileid = fileName.split('.', 1)[0]
		fileid = fileName
		if org not in fullOrg.keys():
			fullOrg[org] = []
		fullOrg[org].append(fileid)

		isFile = os.path.join(dir2out4org, org, fileid+'.out')
		if os.path.isfile(isFile) and os.stat(isFile).st_size > 0:
			islist = rdISfile(isFile)
			if len(islist) < 1:
				continue
			if org not in mislist.keys():
				mislist[org] = {}
			# The organism with ISs is retained.
			mislist[org][fileid] = islist
		else:
			print('Not valid file containing prediected IS elements', isFile)
	
	# check ration of multi-copy IS vs single-copy IS per genome (species)
	#checkCopies(mislist)

	# Cluster species with IS elements

	family = set() # familyNames availalbe in all species in data set
	family4org = {}
	# family4org: {org:{family:n, ...}, ...}
	# n: number of the specific family available in the specific organism
	for org, islist4org in mislist.items():
		family4org[org] = {}
		for islist in islist4org.values():
			for element in islist:
				if IStype2family == True:
					familyName = element['family']
				else:
					familyName = element['familyCluster']
				family.add(familyName)
				if familyName not in family4org[org].keys():
					family4org[org][familyName] = 0
				family4org[org][familyName] += 1
	# family only includes familyName appearing in the examined species (orgs),
	# len(family) <= numberOfFamilyInISfinder (27 in current isfinder)
	# len(familyCluster) <= numberOfHMM (496 in current isfinder with the threshold of 30% in cdhit)
	familySorted = sorted(family)

	# cluster organisms by occurence of IS family or familyCluster
	org4is4sorted = sorted(mislist.keys())

	boolSwitch = True
	distMatrix = distBetweenOrgsByIS(org4is4sorted, family4org, familySorted, boolSwitch=boolSwitch)

	'''
	# tree: returned by scipy.cluster.hierarchy.to_tree()
	distMatrix, tree = clusterOrgByIS(org4is4sorted, family4org, familySorted, boolSwitch=boolSwitch)
	# dengrogram of hierachical clustering
	plt.figure()
	plt.plot(sorted(distMatrix), 'bo',)
	plt.xlabel('Pair of species')
	plt.ylabel('dissimilarity')
	#plt.axis([0,2700000, 0,1.1])
	plt.show()
	'''

	print('Number of organism:', len(fullOrg))
	print('Number of organism with IS:', len(org4is4sorted))
	print('Number of organism (with IS) pairs:', len(distMatrix))
	print(distMatrix[:100])

	"""
	distMatirx is a one-dimensional array, it is considered a condensed matrix of pairwise 
	dissimilarities in the format which is returned by scipy.spatial.distance.pdist. It
	contains the flattened, upper-triangular part of a pairwise dissimilarity matrix. That
	is, if there are N data points and the matrix d contains the dissimilarity between the
	i-th and j-th observation at position d[i][j] , the vector X has length N*(N-1)/2
	and is ordered as follows: 
	d[0][1], d[0][2], ..., d[0][n-1], d[1][2], ..., d[n-2][n-1] = 
	X[0], X[1], ..., X[n-2], X[n-1], ..., X[n*(n-1)/2-1]
	"""

	n = len(org4is4sorted) # number of organisms with IS element identified
	iorg = set()
	# iorg: {i, ...}, set of index of organism sharing IS type (dist=0) with any other organism
	ijorg = set() 
	# ijorg: {(i,j), ...}, set of pair of indice of organisms with IS type match (dist=0)
	c = 0 
	for ij,d in zip(itertools.combinations(range(n),2), distMatrix):
		c += 1
		if d < 1:
			i, j = ij # ij: (i, j), tuple
			iorg.update(ij)
			ijorg.add(ij)

	print('potential organisms involved in HGT:', len(sorted(iorg)))
	print('potential organism pairs involved in HGT:', len(sorted(ijorg)))


	# rnaMap: {org:rnas, ...}
	# rnas: {fileid:dna2rna, ...}, mapping DNA to rna for each DNA sequence in the specific organism
	# dna2rna: {'qseqid':qseqid, 'sseqid':sseqid, 'pident':pident, 'length':length, 'qlen':qlen, 'slen':slen}
	# qseqid: query sequence ID, eg. 'gi|556503834|ref|NC_000913.3|'
	# sseqid: subject sequence ID, eg. 'CP009685.3470134.3471689', 'CP000946.4686640.4688173'
	'''
	file4rnaAccid = '16sRNAacc.list'
	rnaMapRef = genome2rna(dnaFiles, silva, dir4blastout, perc_ident4mapRNA2genome, coverage2mapRNA2genome)
	outputAccid(file4rnaAccid, rnaMapRef)
	'''
	file4rnaAccid = '16sRNAacc.list.Ref.exactGene'
	rnaMapRef = id4Genome2rna(file4rnaAccid)
	# dnaFiles: [(file, org), ..., (file, org)]
	# genomesMeta = {org:meta4seqs, ...}
	# meta4seqs: {fileid:seqMeta}
	# seqMeta: {'dnaType', type}, type can be 0,1,2 to represent phage or plasmid or chromosome DNA, respectively
	genomesMeta = {}
	for file_org in dnaFiles:
		fnafile, org = file_org
		fileid = os.path.basename(fnafile).rsplit('.', 1)[0]
		seqMeta = getMeta(fnafile)
		if org not in genomesMeta.keys():
			genomesMeta[org] = {}
		genomesMeta[org][fileid] = seqMeta

	# make a copy of rnaMapRef, rnaMap, and add meta info to rnaMap
	orgs4genomeMeta = genomesMeta.keys()
	rnaMap = {}
	for org,rna4seq in rnaMapRef.items():
		if org not in orgs4genomeMeta:
			continue
		rnaMap[org] = rna4seq
		for fileid in rna4seq.keys():
			rnaMap[org][fileid]['dnaType'] = genomesMeta[org][fileid]['dnaType']

	# rnaMap: {org:rnas, ...}
	# rnas: {fileid:dna2rna, ...}, mapping DNA to rna for each DNA sequence in the specific organism
	# dna2rna: {'dnaType':type, 'qseqid':qseqid, 'sseqid':sseqid, ...}
	norgWithMultipleRNA = 0
	for org in sorted(rnaMap.keys()):
		rnas = rnaMap[org]
		nseqWithRNA = len(rnas)
		if nseqWithRNA > 1:
			norgWithMultipleRNA += 1
			print('multiple genome sequences have 16S rRNA gene in {}: {}'.format(
				org, [(fileid, dna2rna['dnaType']) for fileid,dna2rna in rnas.items()]))
	print('number of organisms with multiple genome sequences mapped to 16S rRNA:', norgWithMultipleRNA)

	# Keep in rnaMap only the dna2rna for chromosome DNA sequence.
	rnaMap4chrom = {}
	norg4mchroms = 0
	for org,rnas in rnaMap.items():
		rnas4chrom = {}
		for fileid,dna2rna in rnas.items():
			if dna2rna['dnaType'] == 2:
				rnas4chrom[fileid] = dna2rna
		nchrom = len(rnas4chrom)
		if nchrom < 1:
			print('no chromosome DNA occurs in organism:', org)
		else:
			if nchrom > 1:
				print('multiple chromosome DNAs ({}) in {}'.format(rnas4chrom.keys(), org))
				norg4mchroms += 1
			rnaMap4chrom[org] = rnas4chrom

	orgs4rna = set(rnaMap.keys())	
	orgs4NoRNA = set(fullOrg.keys()) - orgs4rna
	print('number of organisms without 16S rRNA gene:', len(orgs4NoRNA))
	for org in sorted(orgs4NoRNA):
		print('No 16S rRNA gene found for', org)
	accids4rna = set()
	for rnas in rnaMap.values():
		for rna in rnas.values():
			accids4rna.add(rna['sseqid'])
	print('{} genomes (1 or more DNA sequences) mapped to {} 16S rRNA genes from silva REF database'.format(
		len(orgs4rna), len(accids4rna)))
	print('number of organisms with chromosome DNA mapped to rna gene:', len(rnaMap4chrom))
	org4nonChromRNA = orgs4rna - set(rnaMap4chrom.keys())
	print('number of organisms with non-chromosome DNA mapped to rna gene:', len(org4nonChromRNA), org4nonChromRNA)
	print('number of organisms with multiple chromosome DNAs to rna genes:', norg4mchroms)

	print('Discuss running before getSilvaTax()', datetime.datetime.now().ctime())
	# Assign taxmonomy path (hierachical rank) to each genome (organism)
	# Note: the deepest rank assigned in silva is genus which is the upper level of species.
	# Note: organism 'Puccinia triticina' occurs in multiple taxonomy paths in Eukaryota, Bacteria and Archaea, 
	#	it is interesting.
	silvaTax = '/home/data/insertion_sequence/silva/taxonomy/taxmap_slv_ssu_ref_123.1.txt'
	#silvaTax = '/home/data/insertion_sequence/silva/taxonomy/taxmap_slv_ssu_ref_123.1.txt.noEukaryota'
	# tax: {seqid:tax4seq, ...}
	# tax4seq: {'taxpath':taxpath, 'org':org, 'taxid':taxid}
	# taxpath: [rank1, rank2, ...], e.g.: 
	#	['Bacteria', 'Spirochaetae', 'Spirochaetes', 'Spirochaetales', 'Spirochaetaceae', 'Treponema 2']
	#	['Bacteria', 'Proteobacteria', 'Alphaproteobacteria', 'Rickettsiales', 'Mitochondria']
	#	['Bacteria', 'Cyanobacteria', 'Chloroplast']
	#	['Archaea', 'Thaumarchaeota', 'Marine Group I']
	# org: orgnism name, e.g. 'Treponema berlinense', 'Chlorella vulgaris', 'Puccinia triticina'
	tax = getSilvaTax(silvaTax)
	print('Discuss running after getSilvaTax()', datetime.datetime.now().ctime())

	#checkTaxLen(tax)


	# Update rnaMap by attaching taxonomy info to each DNA sequence in each genome
	# rnaMapTax: {org:rnas, ...}
	# rnas: {fileid:dna2rna, ...}, mapping DNA to rna for each DNA sequence in the specific organism
	# dna2rna: {'tax4seq':tax4seq, 'dnaType':type, 'qseqid':qseqid, 'sseqid':sseqid, ...}
	# tax4seq: {'taxpath':taxpath, 'org':org, 'taxid':taxid}
	# taxpath: [rank1, rank2, ...], e.g. ['Archaea', 'Thaumarchaeota', 'Marine Group I']
	rnaMapTax = org2tax(rnaMap, tax)

	# check if multiple different taxons are assigned (mapped) to one organism.
	# Return updated rnaMap including only organisms with only one taxon assigned to one or more chromosome DNA
	# in the organism.
	rnaMap4oneTax2chrom = checkOrg4multipleTax(rnaMapTax)

	# output organisms with each dna (chromosome, plasmid, phage) mapped to rna gene mapped to taxon (taxpath, taxid)
	#outputOrg2dna2rna2tax(file4rnaMapTax, rnaMap4oneTax2chrom)

	'''
	#---
	# add taxonomy to is.sum-like file
	# note: each organism must has only one taxonomy, 
	#	e.g. the taxonomy mapped to one chromosome DNA in the organism.
	sumfileOld = args['isSum']
	sumfileNew = os.path.basename(sumfileOld) + '.tax'
	#tax2issum(sumfileOld, sumfileNew, rnaMap4oneTax2chrom)

	# Get summarization from is.sum like file
	# org2istax: {org:istax, ...}
	# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
	# isfamily: {familyname:(nIS,norg), ...}
	# taxpath: [rank0, rank1, ...]
	sumtaxfile = sumfileNew
	org2istax4sumfile = getsumtax(sumtaxfile)
	# remove organisms is not in fullOrg
	# fullOrg: {org:fileids, ...}, full list of organisms from input file list, e.g. fna.list
	# fileids: [fileid, ...], all file IDs belonging to one organism (species)
	org2istax4sumfile = {k:v for k,v in org2istax4sumfile.items() if k in fullOrg.keys()}
	org2istax = org2istax4sumfile
	#---
	'''

	#===
	# Get summarization from mislist which is retrieved from fileid.out file for the fasta file of each DNA sequence, 
	# organisms without IS are not included.

	# dnaFiles: [(file, org), ..., (file, org)]
	seqmeta4org = getseqmeta(dnaFiles)
	# seqmeta4org: {org:meta4org, ...}
	# meta4org: {seqid:metainfo, ...}
	# metainfo: {'dnaType':dnaType, 'dnaLen':dnaLen}

	# summarize IS elements in each organism
	#
	# mislist: {org:{fileid:islist, ...}, ...}
	# islist: [is, ..., is]
	# is: {'seqid':seqid, 'family':family, 'familyCluster':familyCluster, 
	#		'start':start, 'end':end, 'ncopy4is':ncopy4is}
	sum4org = getsum4org(mislist, seqmeta4org)
	# sum4org: {org:s, ...}, organisms without IS are included now.
	# s: {'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is, 
	#		'family':family4org}
	# family4org: {familyname:(nIS,norg), ...}
	# seqs: {seqid:sum4seq, ...}
	# sum4seq: {'nIS':nIS, 'family':family, 'metainfo':metainfo}
	# family: {familyname:nIS, ...}

	org2istax = getsumtax4org(sum4org, rnaMap4oneTax2chrom)
	# org2istax: {org:istax, ...}
	# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
	# isfamily: {familyname:[nIS,norg], ...}
	#===


	'''
	# To ensure that data from is.sum.tax and data from the combination of mislist and rnaMap
	# are exact same and using different data source and methods will not produce the different HGT events.
	for org1,org2 in zip(sorted(org2istax4sumfile.keys()), sorted(org2istax)):
		istax1 = org2istax4sumfile[org1]
		istax2 = org2istax[org2]
		for k1,k2 in zip(sorted(istax1.keys()), sorted(istax2.keys())):
			if istax1[k1] != istax2[k2]:
				print('hello diff istax', org1, org2, k1,k2,istax1[k1],istax2[k2]) 
				print([(k,istax1['is'][k]) for k in sorted(istax1['is'].keys())])
				print([(k,istax2['is'][k]) for k in sorted(istax2['is'].keys())])
				return 0
	'''

	# remove organisms without assigned taxonomy clade name or IS identified
	org2istaxcore = {k:v for k,v in org2istax.items() if len(v['taxpath'])>0 and v['nIS']>0}
	orgnames4core = set(org2istaxcore.keys())
	print('hello norg2istax and norg2istaxcore', len(org2istax), len(org2istaxcore))

	#isVSgenomeSizeBin(org2istax)

	
	# Summarize organisms based on taxonomy assignment and output results of analysis on org2istax
	#
	# taxsum: {taxon:sum, ...}
	# 	{rank0:sum, ...} means {'':sum, 'Bacteria':sum, 'Archaea':sum}
	#	{(rank0,rank1):sum, ...} means {'':sum, ('Bacteria','Cyanobacteria'):sum, ('Bacteria,Proteobacteria'):sum, ..., 
	#		('Archaea,Crenarchaeota'):sum, ('Archaea,Euryarchaeota'):sum, ...}
	# sum: {'norg':norg, 'norgis':norgis, 'nIS':nIS, 'dnaLen':dnaLen, 
	#		'orgs':orgs, 'orgs4is':orgs4is, 'familys':familys}
	# norg: number of organisms assigned to the specific taxon
	# norgis: number of organisms (with IS identified) assigned to the specific taxon
	# nIS: number of ISs occuring
	# dnaLen: total length of all DNAs in an organism, namely, size of a genome
	# orgs: [org, ...], names of organisms
	# orgs4is: [org, ...], names of organisms with IS identified
	# familys: {family1:[nIS, norg], ...}
	#	 norg: number of organisms containing the specific IS family
	rank = args['rank']
	ratio = args['ratio']
	ofile = 'analysis2org2istax.tax' + str(rank) + '.' + str(ratio)
	taxsum = analysis2org2istax(org2istax, rank, ofile)
	ofile = 'analysis2org2istax.taxcore' + str(rank) + '.' + str(ratio)
	taxsumcore = analysis2org2istax(org2istaxcore, rank, ofile)


	# Find potential organism pairs in HGT events
	#
	# org4is4sorted: [org, ...], sorted organism names with IS identified
	# iorg: {i, ...}, set of index of organism sharing IS type (dist=0) with any other organism
	# ijorg: {(i,j), ...}, set of pair of indice of organisms with IS type match (dist=0)
	# ijorgnames: [(orgi, orgj), ...], orgi occurs in org4is4sorted earlier than orgj,
	#	orgi and orgj are organism names sharing at least one IS type (dist=0) with each other.
	ijorgnames = [(org4is4sorted[ij[0]], org4is4sorted[ij[1]]) for ij in ijorg]
	ijorgnames = {p for p in ijorgnames if set(p) <= orgnames4core}
	orgs4difftax = findOrgpairsInDiffTax(ijorgnames, taxsumcore)
	# orgs4difftax: [(orgi, orgj), ...], 
	#	orgi and orgj share at least one IS type but assigned to the different taxonomy clades.
	norgs4difftax = len(orgs4difftax)
	npairs4difftax = 0
	for vpair in itertools.combinations(taxsumcore.values(), 2):
		# vpair[0] and vpair[1] are from different clades.
		npairs4difftax += len(vpair[0]['orgs4is']) * len(vpair[1]['orgs4is'])
	print('{:.0f} ({}/{}) potential HGT organism pairs (in different clades) share common IS, rank={} clades={}'.format(
		100*norgs4difftax/npairs4difftax, norgs4difftax, npairs4difftax, rank, len(taxsumcore)))

	# get the list of organisms involved in potential HGT events
	orgnames4orgpairs = set()
	for p in orgs4difftax:
		orgi, orgj = p
		orgnames4orgpairs.update(p)
	print('number of organisms involved in potential HGT:', len(orgnames4orgpairs))

	# get list of IS types identified in each org in orgpairs
	#
	# family4org: {org:{family:n, ...}, ...}
	# n: number of the specific family available in the specific organism
	family4orgs = {}
	# family4orgs: {org:familys, ...}, different from family4org because family4org may contain familys with n=0
	#	but family4orgs contain only familys with n > 0
	# familys: {family, ...}, set
	family4orgs = {}
	for org,familys in family4org.items():
		# remove org not in orgnames4orgpairs
		if org not in orgnames4orgpairs:
			continue
		# for familys in each org, remove fake family with 0 IS identified in the org
		family4orgs[org] = {family for family,n in familys.items() if n > 0}

	# find common ISs shared by a pair of two organisms from different taxonomy clades
	orgpairs4commonIS = {}
	# orgpairs4commonIS: {orgpair:commonISs, ...}
	# orgpair: (orgi, orgj), pair of organism names
	# commonISs: {family, ...}, IS type (e.g. family or familyCluster) shared by two oragnisms
	for p in sorted(orgs4difftax):
		orgi, orgj = p
		orgpairs4commonIS[p] = family4orgs[orgi] & family4orgs[orgj]

	#---
	# Observe the distribution of ratios for all ISs in all taxonomy clades
	# Print: IS, clade, ratio
	ratioCutoff = {}
	# ratioCutoff: {(taxon,common):info4ratio, ...}
	# info4ratio: {'ratio':ratio, 'norgcommon':norgcommon, 'norg':norg, ...}
	for taxon,sum in taxsumcore.items():
		norg = taxsum[taxon]['norg']
		for common,v in taxsumcore[taxon]['familys'].items():
			norgcommon = v[1]
			ratioCutoff[(taxon,common)] = {}
			ratioCutoff[(taxon,common)]['ratio'] = norgcommon/norg
			ratioCutoff[(taxon,common)]['norgcommon'] = norgcommon
			ratioCutoff[(taxon,common)]['norg'] = norg
	dir2hgt = '/home/data/insertion_sequence'
	ofile = os.path.join(dir2hgt, 'ratio.list'+str(rank))
	lines = []
	headerline = '{:>6} {:>10} {:>4} {:<140} {:<15}'.format(
			'ratio', 'norgcommon', 'norg', 'taxpath', 'common')
	lines.append(headerline)
	for k in sorted(ratioCutoff.keys(), key = lambda x: ratioCutoff[x]['ratio']):
		v = ratioCutoff[k]
		line = '{:>6.2f} {:>10} {:>4} {:<140} {:<15}'.format(
				v['ratio']*100, v['norgcommon'], v['norg'], 
				';'.join(k[0]), k[1])
		lines.append(line)
	with open(ofile, 'w') as fp:
		fp.write('\n'.join(lines)+'\n')
	#---


	# hunt for HGT organism pairs in orgpairs4commonIS
	#
	# org2istax: {org:istax, ...}
	# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
	# isfamily: {familyname:(nIS,norg), ...}
	# taxpath: [rank0, rank1, ...]
	#
	# taxsum: {taxon:sum, ...}
	# sum: {'norg':norg, 'norgis':norgis, 'nIS':nIS, 'dnaLen':dnaLen, 
	#		'orgs':orgs, 'orgs4is':orgs4is, 'familys':familys}
	# norg: number of organisms assigned to the specific taxon
	# norgis: number of organisms (with IS identified) assigned to the specific taxon
	# nIS: number of ISs occuring
	# dnaLen: total length of all DNAs in an organism, namely, size of a genome
	# orgs: [org, ...], names of organisms
	# orgs4is: [org, ...], names of organisms with IS identified
	# familys: {family1:[nIS, norg], ...}
	#	 norg: number of organisms containing the specific IS family
	#
	# taxon = tuple(taxpath[:rank])
	hgts = {}
	# htgs: {orgpair:htg, ...}
	# orgpair: (orgi, orgj)
	# hgt: {common:event, ...}
	# event: {'ratioi':ratioi, 'ratioj':ratioj, ...}
	# IStype: name of IStype, e.g. familyname or familyClusterName
	ratio = args['ratio']
	ratioCutoff = float(ratio) # norgWithCommon/norg
	norgCutoff = 3 # minimal value of number of organisms in the clade involved in HGT
	#for p,commons in orgpairs4commonIS.items():
	for p in sorted(orgpairs4commonIS.keys()):
		commons = orgpairs4commonIS[p]
		orgi, orgj = p
		# in the taxonomy clade which the organism is assigned to
		taxoni = tuple(org2istaxcore[orgi]['taxpath'][:rank])
		taxonj = tuple(org2istaxcore[orgj]['taxpath'][:rank])
		# number of organisms (with or without IS )
		norgi = taxsum[taxoni]['norg']
		norgj = taxsum[taxonj]['norg']
		# find number of organisms with same IStype in same clade for each common IStype in orgpairs
		for common in sorted(commons):
			# number of organisms with the specific common IS in the specific clade
			norgcommoni = taxsumcore[taxoni]['familys'][common][1]
			norgcommonj = taxsumcore[taxonj]['familys'][common][1]
			# norgWithCommon/norg
			ratioi = norgcommoni/norgi
			ratioj = norgcommonj/norgj
			hgt = {}
			# hgt: {common:event, ...}
			# event: {'ratioi':ratioi, 'ratioj':ratioj, ...}
			if (ratioi<ratioCutoff and norgi>=norgCutoff) or (ratioj<ratioCutoff and norgj>=norgCutoff):
				event = {}
				event['ratioi'] = ratioi
				event['ratioj'] = ratioj
				hgt[common] = event
		if len(hgt) > 0:
			hgts[p] = hgt
	
	# Output HGT events
	# the longest taxpath contains 140 characters.
	len4taxpath = 125
	str4fmt4taxon = '{:<' + str(len4taxpath) + '} '

	#  max length of familyCluster name is 64
	#fmt4title = str4fmt4taxon*2 + '{:<85} {:<85} {:<64} {:>10} {:>7} {:>4} {:>10} {:>7} {:>4} {:>6} {:>6} {:>10} {:>6} {:>7} {:>5} {:>10} {:>6} {:>7} {:>5}'
	#fmt = str4fmt4taxon*2 + '{:<85} {:<85} {:<64} {:>10} {:>7} {:>4} {:>10} {:>7} {:>4} {:>6.2f} {:>6.2f} {:>10} {:>6} {:>7} {:>5} {:>10} {:>6} {:>7} {:>5}'

	# Print IS acceptor preceding IS donor in each line as it is more convenient to group events by acceptor and then remove
	# the redundant events with same acceptor (combination of the organism and the transfered IS) but different donor.
	fmt4title = '{:<64} ' + str4fmt4taxon + '{:<85} {:>10} {:>7} {:>4} {:>6} {:>10} {:>6} {:>7} {:>5} ' + str4fmt4taxon + '{:<85} {:>10} {:>7} {:>4} {:>6} {:>10} {:>6} {:>7} {:>5}'
	fmt = '{:<64} ' + str4fmt4taxon + '{:<85} {:>10} {:>7} {:>4} {:>6.2f} {:>10} {:>6} {:>7} {:>5} ' + str4fmt4taxon + '{:<85} {:>10} {:>7} {:>4} {:>6.2f} {:>10} {:>6} {:>7} {:>5}'
	# 'hgtIS', 
	# 'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',
	# 'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',

	lines = []
	headerline = fmt4title.format(
			# IS type in HGT event, shared by the pair of organisms
			'hgtIS', 

			# taxpaths of the clades assigned to the pair of organisms
			'taxpath', 
			# pair of organisms related to HGT
			'organism', 
			# number of organisms with the specific common IS in the specific clade
			'norg4HGTIS', 
			# number of organisms with IS in the specific clade
			'norg4IS', 
			# number of organisms (with or without IS ) in the specific clade
			'norg', 
			# ratio of the organisms with the specific common IS among all organisms in the specific
			# clade assigned to the specific organism
			'ratio',
			# number of the specific common ISs (HGT IStype) in the specific organism
			'nhgtIS4org',
			# number of the specific common ISs in the specific clade
			'nhgtIS',
			# number of ISs in the specific organism
			'nIS4org',
			# number of ISs in the specific clade
			'nIS',

			'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',
			)
	# taxsum: {taxon:sum, ...}
	# sum: {'norg':norg, 'norgis':norgis, 'nIS':nIS, 'dnaLen':dnaLen, 
	#		'orgs':orgs, 'orgs4is':orgs4is, 'familys':familys}
	# familys: {family1:[nIS, norg], ...}
	#	 norg: number of organisms containing the specific IS family
	# org2istax: {org:istax, ...}
	# istax: {'is':isfamily, 'taxpath':taxpath, 'nIS':nIS, 'dnaLen':dnaLen, 'dnaLen4is':dnaLen4is}
	# isfamily: {familyname:(nIS,norg), ...}
	orgset = set()
	hgtISset = set()
	hgtEventSet = set()
	lines.append(headerline)
	hgtevents = []
	# hgtevents: [(orgA, orgD, common), ...]
	# orgA,orgD,common: acceptor organism, donor organism, transfered IS
	ratiodic = {}
	# ratiodic: {(org, common):ratio}
	for p in hgts.keys():
		orgset.update(p)
		orgi, orgj = p
		taxoni = tuple(org2istaxcore[orgi]['taxpath'][:rank])
		taxonj = tuple(org2istaxcore[orgj]['taxpath'][:rank])
		norg4isi = taxsumcore[taxoni]['norgis']
		norg4isj = taxsumcore[taxonj]['norgis']
		norgi = taxsum[taxoni]['norg']
		norgj = taxsum[taxonj]['norg']
		nISi = taxsumcore[taxoni]['nIS']
		nISj = taxsumcore[taxonj]['nIS']

		hgt = hgts[p]
		for common in sorted(hgt.keys()):
			hgtISset.add(common)
			hgtEventSet.add((p, common))
			nIScommoni, norgcommoni = taxsumcore[taxoni]['familys'][common]
			nIScommonj, norgcommonj = taxsumcore[taxonj]['familys'][common]
			nIScommon4orgi = org2istaxcore[orgi]['is'][common][0]
			nIScommon4orgj = org2istaxcore[orgj]['is'][common][0]
			nIS4orgi = org2istaxcore[orgi]['nIS']
			nIS4orgj = org2istaxcore[orgj]['nIS']

			hgtcommon = hgt[common]

			# define IS acceptor and donor
			if hgtcommon['ratioi'] < hgtcommon['ratioj']:
				taxonA, orgA, norgcommonA, norg4isA, norgA = taxoni, orgi, norgcommoni, norg4isi, norgi
				ratioA = hgtcommon['ratioi']
				nIScommon4orgA, nIScommonA, nIS4orgA, nISA = nIScommon4orgi, nIScommoni, nIS4orgi, nISi
				taxonD, orgD, norgcommonD, norg4isD, norgD = taxonj, orgj, norgcommonj, norg4isj, norgj
				ratioD = hgtcommon['ratioj']
				nIScommon4orgD, nIScommonD, nIS4orgD, nISD = nIScommon4orgj, nIScommonj, nIS4orgj, nISj
			else:
				taxonA, orgA, norgcommonA, norg4isA, norgA = taxonj, orgj, norgcommonj, norg4isj, norgj
				ratioA = hgtcommon['ratioj']
				nIScommon4orgA, nIScommonA, nIS4orgA, nISA = nIScommon4orgj, nIScommonj, nIS4orgj, nISj
				taxonD, orgD, norgcommonD, norg4isD, norgD = taxoni, orgi, norgcommoni, norg4isi, norgi
				ratioD = hgtcommon['ratioi']
				nIScommon4orgD, nIScommonD, nIS4orgD, nISD = nIScommon4orgi, nIScommoni, nIS4orgi, nISi


			# 'hgtIS', 
			# 'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',
			# 'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',
			line = fmt.format(
				common,
				';'.join(taxonA), orgA, norgcommonA, norg4isA, norgA, ratioA, nIScommon4orgA, nIScommonA, nIS4orgA, nISA,
				';'.join(taxonD), orgD, norgcommonD, norg4isD, norgD, ratioD, nIScommon4orgD, nIScommonD, nIS4orgD, nISD,
				)
			lines.append(line)

			# store acceptor and donor info for each HGT event
			hgtevents.append((orgA, orgD, common))
			ratiodic[(orgA,common)] = ratioA
			ratiodic[(orgD,common)] = ratioD
	print('HGT numbers: orgs={} hgtIS={} orgpairs={} hgtEvents={}'.format(
		len(orgset), len(hgtISset), len(hgts.keys()), len(hgtEventSet)))

	dir2hgt = '/home/data/insertion_sequence'
	ofile = os.path.join(dir2hgt, 'HGTs.list'+str(rank)+'.'+str(ratio))
	with open(ofile, 'w') as fp:
		fp.write('\n'.join(lines)+'\n')

	# Get HGT events with each (orgA,common) mapped to only one orgD with maximal ratioD
	#
	# sort and group all events by (orgA,common), and then retain only event (orgA,common,orgD) with the 
	# maximal ratioD in each events group
	events = {}
	# events:{(orgA,common):orgD, ...}
	import operator
	for k,g in itertools.groupby(sorted(hgtevents, key=operator.itemgetter(0,2)), key=operator.itemgetter(0,2)):
	# hgtevents: [(orgA, orgD, common), ...]
		g = list(g)
		# just get it if only one orgD could be found for the specific acceptor event (orgA,common)
		if len(g) == 1:
			events[k] = g[0][1]
			continue

		# retain only orgD with maximal ratioD in a group of redundant events (orgDs)
		# There are probably more than one events with maxRatioD, we need to pick the most reasonable orgD
		# from multiple orgDs.
		eventgroupsByRatioD = itertools.groupby(sorted(g, key = lambda x: ratiodic[(x[1],x[2])], reverse=True), 
								key = lambda x: ratiodic[(x[1],x[2])])
		eventgroup4maxRatioD = list(next(eventgroupsByRatioD)[1])
		# next(eventgroups4maxRatioD) is like (k,g).
		events[k] = eventgroup4maxRatioD[0][1]
		continue
		if len(eventgroup4maxRatioD) == 1:
			events[k] = eventgroup4maxRatioD[0][1]
			continue

		# We now need to pick the most probable orgD from multiple orgDs with maxRatioD.
		# Option1: 
		#	1) pick orgD4common with the most number of common ISs in its genome.
		#	2) pick orgD4is with the most number of ISs in its genome if multiple orgD4commonISs exist.
		#	3) simply pick the first orgD if still multiple orgD4is exist.
		eventgroupsByCommon = itertools.groupby(
				sorted(eventgroup4maxRatioD, key=lambda x: org2istaxcore[x[1]]['is'][x[2]][0], reverse=True),
				key = lambda x: org2istaxcore[x[1]]['is'][x[2]][0])
		eventgroupByCommon = list(next(eventgroupsByCommon)[1])
		if len(eventgroupByCommon) == 1:
			events[k] = eventgroupByCommon[0][1]
			continue

		# pick orgD4is with the most number of ISs in its genome if multiple orgD4commonISs exist
		eventgroupsByIS = itertools.groupby(
				sorted(eventgroupByCommon, key=lambda x: org2istaxcore[x[1]]['nIS'], reverse=True), 
				key = lambda x: org2istaxcore[x[1]]['nIS'])
		eventgroupByIS = list(next(eventgroupsByIS)[1])
		events[k] = eventgroupByIS[0] # simply pick the first orgD if one or more orgD4is exist
		#if len(eventgroupByIS) == 1:
		#	events[k] = eventgroupByIS[0][1]
		#	continue

	# Output HGT events with each (orgA,common) mapped to only one orgD with maximal ratioD
	#
	# ofile = os.path.join(dir2hgt, 'HGTs.list'+str(rank)+'.'+str(ratio))
	ofile4hgtsUniq = ofile+'.uniq'
	lines = []
	lines.append(headerline)

	orgset = set()
	hgtISset = set()
	orgpairset = set()
	for k in sorted(events.keys()):
		orgA, common = k
		orgD = events[k]
		orgset.update((orgA,orgD))
		hgtISset.add(common)
		orgpairset.add((orgA,orgD))

		taxonA = tuple(org2istaxcore[orgA]['taxpath'][:rank])
		nIScommonA, norgcommonA = taxsumcore[taxonA]['familys'][common]
		norg4isA = taxsumcore[taxonA]['norgis']
		norgA = taxsumcore[taxonA]['norg']
		ratioA = ratiodic[(orgA, common)]
		nIScommon4orgA = org2istaxcore[orgA]['is'][common][0]
		nIS4orgA = org2istaxcore[orgA]['nIS']
		nISA = taxsumcore[taxonA]['nIS']

		taxonD = tuple(org2istaxcore[orgD]['taxpath'][:rank])
		nIScommonD, norgcommonD = taxsumcore[taxonD]['familys'][common]
		norg4isD = taxsumcore[taxonD]['norgis']
		norgD = taxsumcore[taxonD]['norg']
		ratioD = ratiodic[(orgD, common)]
		nIScommon4orgD = org2istaxcore[orgD]['is'][common][0]
		nIS4orgD = org2istaxcore[orgD]['nIS']
		nISD = taxsumcore[taxonD]['nIS']

		# 'hgtIS', 
		# 'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',
		# 'taxpath', 'organism', 'norg4HGTIS', 'norg4IS', 'norg', 'ratio', 'nhgtIS4org', 'nhgtIS', 'nIS4org', 'nIS',
		line = fmt.format(
			common,
			';'.join(taxonA), orgA, norgcommonA, norg4isA, norgA, ratioA, nIScommon4orgA, nIScommonA, nIS4orgA, nISA,
			';'.join(taxonD), orgD, norgcommonD, norg4isD, norgD, ratioD, nIScommon4orgD, nIScommonD, nIS4orgD, nISD,
			)
		lines.append(line)
	print('HGT (uniq) numbers: orgs={} hgtIS={} orgpairs={} hgtEvents={}'.format(
		len(orgset), len(hgtISset), len(orgpairset), len(events)))
	with open(ofile4hgtsUniq, 'w') as fp:
		fp.write('\n'.join(lines)+'\n')


	# Summarize IS and taxon for each IS family
	#
	# sum4istax: {family:sum4family, ...}
	# sum4family: {'nIS':nIS, 'nCluster':nCluster, 'nOrg':nOrg, 
	#		'taxon': {taxon, ...},
	#		}
	# taxon: set, unique taxpath with the specific IS family identified
	ratio = args['ratio']
	ofile = 'sum4family.tax' + str(rank) + '.' + str(ratio)
	sum4istax = getsum4istax(org2istax, rank, ofile)
	ofile = 'sum4family.taxcore' + str(rank) + '.' + str(ratio)
	sum4istaxcore = getsum4istax(org2istaxcore, rank, ofile)

	print('Discuss running finishes at', datetime.datetime.now().ctime())


if __name__ == "__main__":
	descriptionStr = 'Do analysis for discussion, such as horizontonal gene transfer'
	parser = argparse.ArgumentParser(description = descriptionStr)

	helpStr = 'input file containing NCBI genome fasta files, one file per line, e.g. bacteria.fna.list'
	parser.add_argument('fileList', help = helpStr)

	helpStr = 'directory holding the results of IS prediction, one organism per sub-directory, e.g. /home/data/insertion_sequence/results4MGEScan-IS/prediction'
	parser.add_argument('dir2prediction', help = helpStr)

	helpStr = 'input file summarizing IS distribution of all families in all organisms, e.g. ../mgescan-is/is.sum'
	parser.add_argument('isSum', help = helpStr)

	helpStr = 'number of taxonomic rank levels used to classified genomes in summarization result, the value can be [1,2,3,4,5,6]'
	parser.add_argument('rank', help = helpStr)

	helpStr = 'number of the organisms with the shared IS VS number of all organisms in the specific clade which the organism belongs to, the value can be [0.1,0.2,0.3,0.4,0.5], HGT event is defined when caculated_value < ratio'
	parser.add_argument('ratio', help = helpStr)

	args = parser.parse_args()

	args4discuss = {
			'fileList':		args.fileList,
			'dir2prediction':	args.dir2prediction,
			'isSum':		args.isSum,
			'rank':			int(args.rank),
			'ratio':		args.ratio
			}
	discuss(args4discuss)
