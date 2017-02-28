#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#import re
import os
import sys
import argparse
import datetime
import concurrent.futures
import subprocess, shlex

import tools
import constants

# HMM files will be placed in rootpath.
#rootpath = '/data2/zhiqxie/aclame/buildhmm'

# for classifyTpase()
BLASTDB = ''

# clustal-omega
clustalo_cmd = "/u/zhiqxie/informatics/inst/clustal-omega-1.2.1/bin/clustalo"
# hmmer
hmmbuild_cmd = "/u/zhiqxie/informatics/inst/hmmer-3.1b2/bin/hmmbuild"

# cd-hit command to cluster protein sequences at 1.0 identity
CDHIT = '/u/zhiqxie/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit'

# common optins to run cd-hit
option4cdhit = '-g 1 -G 0 -aS 0.8  -d 0 -M 0'

# psi-cd-hit clusters proteins at very low threshold, it also cluster long DNA sequences, through blastp, blastn and metablast. 
# psi-cd-hit is a Perl script, which runs similar incremental algorithm like cd-hit, but using BLAST to calculate similarities. 
# example:
# cd-hit -i nr100 -o nr90 -c 0.9 -n 5 -g 1 -G 0 -aS 0.8  -d 0 -M 0 > nr90.log
thresholds4cdhit90 = ('nr90', '-c 0.9 -n 5')
# cd-hit -i nr90 -o nr60 -c 0.6 -n 4 -g 1 -G 0 -aS 0.8  -d 0 -M 0 > nr60.log
thresholds4cdhit60 = ('nr60', '-c 0.6 -n 4')
# psi-cd-hit.pl -i nr60 -o nr30 -c 0.3 -aS 0.8 -G 0 -g 1

# clstr_rev.pl nr90.clstr nr60.clstr > nr9060.clstr
# clstr_rev.pl nr9060.clstr nr30.clstr > nr906030.clstr
PSICDHIT = '/u/zhiqxie/informatics/inst/cd-hit-v4.6.4-2015-0603/psi-cd-hit/psi-cd-hit.pl'
CLSTR_REV = '/u/zhiqxie/informatics/inst/cd-hit-v4.6.4-2015-0603/clstr_rev.pl'

# switch to nr100, nr90, nr60 hierachical clustering
#LOW_IDENT = False
# switch to nr100, nr90, nr60, nr30 hierachical clustering
LOW_IDENT = True


def build_msa(seqfilename, msafilename):
	infile = "--infile=" + seqfilename
	seqformat = "--infmt=fa"
	seqtype = "--seqtype=Protein"
	outfile = "--outfile=" + msafilename
	msaformat = "--outfmt=st"
	overwrite = "--force"
	ncpu = "--threads=8"
	do_msa = [clustalo_cmd, infile, seqformat, seqtype, outfile, msaformat, overwrite, ncpu]

	return subprocess.call(do_msa, shell=False, universal_newlines=False)

def build_msa_v2(args):
	seqfilename, msafilename = args
	infile = "--infile=" + seqfilename
	seqformat = "--infmt=fa"
	seqtype = "--seqtype=Protein"
	outfile = "--outfile=" + msafilename
	msaformat = "--outfmt=st"
	overwrite = "--force"
	ncpu = "--threads=" + str(constants.nthread)
	do_msa = [clustalo_cmd, infile, seqformat, seqtype, outfile, msaformat, overwrite, ncpu]

	return subprocess.call(do_msa, shell=False, universal_newlines=False)

def is_hmmbuild(msafile, hmmfile):
	seqtype = "--amino"
	msa_format = "--informat stockholm"
	ncpu = '--cpu ' + str(constants.nthread)
	cmd_line = '{0} {1} {2} {3} {4} {5}'.format(hmmbuild_cmd, seqtype, msa_format, ncpu, hmmfile, msafile)
	#cmd_line = '{0} {1} {2} {3} {4}'.format(hmmbuild_cmd, seqtype, msa_format, hmmfile, msafile)
	do_hmmbuild = shlex.split(cmd_line)

	return subprocess.call(do_hmmbuild, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)

def buildMSA(args2concurrent4seq, nproc, fp_log):
	print("\nBegin to build multiple sequence alignments.")
	fp_log.write("Begin to build multiple sequence alignments.\n" + datetime.datetime.now().ctime() + '\n\n')

	if len(args2concurrent4seq) < nproc:
		nproc = len(args2concurrent4seq)
	with concurrent.futures.ProcessPoolExecutor(max_workers = nproc) as executor:
		for arg, outs in zip(args2concurrent4seq, executor.map(build_msa_v2, args2concurrent4seq)):
			seq_filename = os.path.basename(arg[0])
			if outs == 0:
				fp_log.write('Building multiple sequence alignment for ' + seq_filename + ', return ' + str(outs) + '\n')
			else:
				fp_log.write('Building multiple sequence alignment for ' + seq_filename + ', return error!\n')
	print("\nFinish building multiple sequence alignments.")
	fp_log.write("Finish building multiple sequence alignments.\n" + datetime.datetime.now().ctime() + '\n\n')

def buildHMM(args2concurrent4seq, hmms_file, fp_log):
	fp_hmms_file = open(hmms_file, 'w')
	print("\nBegin to build profile HMM models.")
	fp_log.write("Begin to build profile HMM models.\n" + datetime.datetime.now().ctime() + '\n\n')
	# args2concurrent4seq: [(seq_pathfilename, msa_pathfilename), ..., (seq_pathfilename, msa_pathfilename)]
	for item in args2concurrent4seq:
		# msa_file, hmm_file: familyPath/fileName
		# example: IS1/IS1_ISMhu11.long.pep.all.msa, IS1/IS1_ISMhu11.long.pep.all.hmm
		msa_pathfilename = item[1]
		hmm_pathfilename = msa_pathfilename.replace('msa', 'hmm')

		outs = is_hmmbuild(msa_pathfilename, hmm_pathfilename)
		if outs == 0:
			fp_log.write('Building profile HMM model for ' + msa_pathfilename + ', return ' + str(outs) + '\n')
		else:
			fp_log.write('Building profile HMM model for ' + msa_pathfilename + ', return error!\n')

		# write all clusters of HMM models into one hmm file
		fp_hmm_file = open(hmm_pathfilename, 'r')
		fp_hmms_file.write(fp_hmm_file.read())
		fp_hmm_file.close()
		print('Append profile HMM model from {} into {}'.format(hmm_pathfilename, hmms_file))

	print("\nFinish building profile HMM models.")
	fp_log.write("Finish building profile HMM models.\n" + datetime.datetime.now().ctime() + '\n\n')
	fp_hmms_file.close()
	print('Finish appending profile HMM models into', hmms_file)


def readClustersByCDhit(clusterFile):
	clusters = []
	fp = open(clusterFile, 'r')
	i = 0
	cluster = []
	for line in fp:
		line = line.strip()
		if line[0] == '>':
			if len(cluster) > 0:
				clusters.append((i, cluster))
				cluster = []
				i += 1
			continue

		index2 = line.find('...')
		left = line[:index2]
		index1 = left.rfind('>')
		isName = left[index1+1:]
		cluster.append(isName)
	clusters.append((i, cluster))
	fp.close()
	return clusters


def callCDhit4hierarch(familyPepFile, lowIdent=False):
	cdhit = CDHIT
	option = option4cdhit

	# remove identical sequence in same family
	# cmd = ~/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit options
	# options = -i IS5/IS5.long.pep.all -o IS5/IS5.long.pep.all.nr100 -c 1.0 -n 5 -g 1 -G 0 -aS 0.8 -d 0 -M 0
	input = familyPepFile
	output = '.'.join([familyPepFile, 'nr100'])
	options = ' '.join(['-i', input, '-o', output, '-c 1.0 -n 5', option])
	cmd_line = ' '.join([cdhit, options])
	do_cluster = shlex.split(cmd_line)
	subprocess.check_call(do_cluster, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)

	# cluster the nonredundant sequences in same family
	# cmd = ~/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit options
	# options = -i IS5/IS5.long.pep.all.nr100 -o IS5/IS5.long.pep.all.nr100.nr90 -c 0.9 -n 5 -g 1 -G 0 -aS 0.8 -d 0 -M 0
	input = output
	output = '.'.join([input, thresholds4cdhit90[0]])
	options = ' '.join(['-i', input, '-o', output, thresholds4cdhit90[1], option])
	cmd_line = ' '.join([cdhit, options])
	do_cluster = shlex.split(cmd_line)
	subprocess.check_call(do_cluster, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)

	# cluster the nonredundant sequences in same family
	# cmd = ~/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit options
	# options = -i IS5/IS5.long.pep.all.nr100.nr90 -o IS5/IS5.long.pep.all.nr100.nr90.nr60 -c 0.6 -n 4 -g 1 -G 0 -aS 0.8  -d 0 -M 0
	input = output
	output = '.'.join([input, thresholds4cdhit60[0]])
	options = ' '.join(['-i', input, '-o', output, thresholds4cdhit60[1], option])
	cmd_line = ' '.join([cdhit, options])
	do_cluster = shlex.split(cmd_line)
	subprocess.check_call(do_cluster, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)

	# clstr_rev.pl nr90.clstr nr60.clstr > nr9060.clstr
	clstrFile90 = '.'.join([input, 'clstr'])
	clstrFile60 = '.'.join([output, 'clstr'])
	clstrFile9060 = '.'.join([clstrFile90, os.path.basename(clstrFile60)])
	cmd_line = ' '.join([CLSTR_REV, clstrFile90, clstrFile60])
	do_cmd = shlex.split(cmd_line)
	with open(clstrFile9060, 'w') as fp:
		subprocess.check_call(do_cmd, shell=False, universal_newlines=False, stdout=fp)

	if lowIdent == False:
		# return cluster fileName
		return clstrFile9060

	# psi-cd-hit.pl -i nr60 -o nr30 -c 0.3 -aS 0.8 -G 0 -g 1 # accurate but slow mode
	# psi-cd-hit.pl -i nr60 -o nr30 -c 0.3 -aS 0.8 -G 0 -g 0 # fast
	# Note: psi-cd-hit.pl will fail to read input if filename of input is very long.
	input = output
	output = '.'.join([input, 'nr30'])
	options = ' '.join(['-i', os.path.basename(input), '-o', output, '-c 0.3 -aS 0.8 -G 0 -g 1'])
	#options = ' '.join(['-i', os.path.basename(input), '-o', output, '-c 0.3 -aS 0.8 -G 0 -g 0'])
	cmd_line = ' '.join([PSICDHIT, options])
	do_cluster = shlex.split(cmd_line)
	subprocess.check_call(do_cluster, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL, cwd=os.path.dirname(input))

	# clstr_rev.pl nr9060.clstr nr30.clstr > nr906030.clstr
	clstrFile30 = '.'.join([output, 'clstr'])
	clstrFile906030 = '.'.join([clstrFile9060, os.path.basename(clstrFile30)])
	cmd_line = ' '.join([CLSTR_REV, clstrFile9060, clstrFile30])
	do_cmd = shlex.split(cmd_line)
	with open(clstrFile906030, 'w') as fp:
		subprocess.check_call(do_cmd, shell=False, universal_newlines=False, stdout=fp)

	# return cluster fileName
	return clstrFile906030

def callCDhit(familyPepFile):
	cdhit = cdhit

	# remove identical sequence in same family
	# cmd = ~/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit options
	# options = -i IS5/IS5.long.pep.all -o IS5/IS5.long.pep.all.nr100 -c 1.0 -n 5 -M 0 -d 0 -g 1
	nr100File = '.'.join([familyPepFile, 'nr100'])
	options = ' '.join(['-i', familyPepFile, '-o', nr100File, '-c 1.0 -n 5 -M 0 -d 0 -g 1'])
	cmd_line = ' '.join([cdhit, options])
	do_cluster = shlex.split(cmd_line)
	if subprocess.call(do_cluster, shell=False, universal_newlines=False) != 0:
		e = ' '.join(['Fail to run', cmd_line])
		raise RuntimeError(e)

	# cluster the nonredundant sequences in same family
	# cmd = ~/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit options
	# options = -i IS5/IS5.long.pep.all.nr100 -o IS5/IS5.long.pep.all.nr100.nr60 -c 0.6 -n 4 -M 0 -d 0 -g 1
	nrFile = '.'.join([nr100File, thresholds4cdhit[0]])
	options = ' '.join(['-i', nr100File, '-o', nrFile, thresholds4cdhit[1], '-M 0 -d 0 -g 1'])
	cmd_line = ' '.join([cdhit, options])
	do_cluster = shlex.split(cmd_line)
	if subprocess.call(do_cluster, shell=False, universal_newlines=False) != 0:
		e = ' '.join(['Fail to run', cmd_line])
		raise RuntimeError(e)

	# return cluster fileName
	return '.'.join([nrFile, 'clstr'])

def doClusterByGroup(familyFeatures): 
	clusters = []
	# clusters: [cluster, ..., cluster]
	# cluster: (clusterName, g)
	# g: [isFeatures, ..., isFeatures]
	for k, g in itertools.groupby(
			sorted(familyFeatures, key = lambda x: x['group']),
			key = lambda x: x['group']):
		cluster = list(g)
		# clusterName: familyName_groupName
		# group name is just '' or blank spaces
		if k.strip() == '':
			k = 'Null'
		elif k.strip() == '.':
			k = 'Dot'
		clusterName = cluster[0]['familyName'] + '_' + k
		clusters.append((clusterName, cluster))

	return clusters
		
def doClusterByElement(familyFeatures):
	clusters = []
	# clusters: [cluster, ..., cluster]
	# cluster: (clusterName, g)
	# g: [isFeatures, ..., isFeatures]
	for k, g in itertools.groupby(
			sorted(familyFeatures, key = lambda x: x['isName']),
			key = lambda x: x['isName']):
		cluster = list(g)
		# clusterName: familyName_isName
		clusterName = cluster[0]['familyName'] + '_' + k
		clusters.append((clusterName, cluster))

	return clusters


def doClusterBySimilarity(family):
	familyName, familyFeatures, familyTpaseSeqFile = family
	familyPepFile = familyTpaseSeqFile

	# example: ISH3.faa.nr100.nr90.clstr.ISH3.faa.nr100.nr90.nr60.clstr.ISH3.faa.nr100.nr90.nr60.nr30.clstr
	dir, basename = os.path.split(familyPepFile)
	nr90 = '.'.join([basename, 'nr100.nr90.clstr'])
	nr9060 = '.'.join([basename, 'nr100.nr90.nr60.clstr'])
	nr906030 = '.'.join([basename, 'nr100.nr90.nr60.nr30.clstr'])

	if LOW_IDENT == False:
		clusterFile = os.path.join(dir, nr90+'.'+nr9060)
	else:
		clusterFile = os.path.join(dir, '.'.join([nr90, nr9060, nr906030]))

	if os.path.isfile(clusterFile) and os.stat(clusterFile).st_size > 0:
		print('Skip (cd-hit and psi-cd-hit) clustering {}'.format(familyPepFile))
	else:
		clusterFile = callCDhit4hierarch(familyPepFile, lowIdent=LOW_IDENT)

	clustersByCDhit = readClustersByCDhit(clusterFile)
	clusters = []
	isNames = {isFeatures['isName']:isFeatures for isFeatures in family[1]}
	for g in clustersByCDhit:
		index4cluster, memberNames4cluster = g
		clusterName = familyName + '_' + str(index4cluster)
		cluster = []
		#print('hello', len(set(isNames.keys())), len(set(memberNames4cluster)))
		commonNames = set(isNames.keys()) & set(memberNames4cluster)
		#print('hello1 ncommonNames={} for cluster {}'.format(len(commonNames), clusterName))
		for isName in commonNames:
			isFeatures = isNames[isName]
			isFeatures['clusterName'] = clusterName
			cluster.append(isFeatures)
		clusters.append((clusterName, cluster))
		
	return clusters

# Return all clusters of IS elements where clustering is based on transposase similarity.
# Return mCluster
# mCluster: [cluster, ...]
# cluster: (clusterName, members)
# members: [isFeatures, ...]
# isFeatures: {'isName': isName, ..., 'seq': seq}
#
# mfamilyFeatures: [(familyName, familyFeatures, familyTpaseSeqFile), ...]
# familyFeatures: [isFeatures, ...]
# isFeatures: {'isName': isName, ..., 'seq': seq}
# seq: character string
#
def getCluster(mfamilyFeatures):
	mCluster = []
	for family in mfamilyFeatures:
		# clusters: [(clusterName, cluster), ...]
		# cluster: [isFeatures, ..., isFeatures]
		mCluster.extend(doClusterBySimilarity(family))
	return mCluster

# mCluster: [cluster, ..., cluster]
# cluster: (clusterName, members)
# members: [isFeatures, ..., isFeatures]
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
def outputClusters(mCluster, fp):
	print('\nCluster list in families', file=fp)
	print('-' * 31, file=fp)
	print('{:<23} {:>8}'.format('clusterName', 'nMembers'), file=fp)
	print('-' * 31, file=fp)
	nCluster= 0
	nTotal = 0
	for cluster in mCluster:
		nCluster += 1
		nMembers = len(cluster[1])
		nTotal += nMembers
		print('{:<23} {:>8}'.format(cluster[0], nMembers), file=fp)
		#if nMembers > 1:
		#	e = '{:<23} {:>8}'.format(cluster[0], nMembers)
		#	raise RuntimeError(e)
	print('Number_of_clusters: {}, Number_of_members_in_all_clusters: {}'.format(nCluster, nTotal), file=fp)
	# all members in each cluster
	print('\n{0:<23} {1:<11}'.format('clusterName', 'memberName'), file=fp)
	print('-' * 34)
	for cluster in mCluster:
		for member in cluster[1]:
			print('{0:<23} {1:<11}'.format(cluster[0], member['isName']), file=fp)


# Return mfamilyFeatures
# mfamilyFeatures: [(familyName, familyFeatures, familyTpaseSeqFile), ...]
# familyFeatures: [isFeatures, ...]
# isFeatures: {'isName': isName, ..., 'tpase': seq}
# seq: character string
def getFamilyFeatures(tpaseSeqFile):
	mfamilyFeatures = []
	# seqs: [(header, seq), ..., (header, seq)]
	# seq: character string
	seqs = tools.getFastaFull(tpaseSeqFile)
	fileName = os.path.basename(tpaseSeqFile)
	familyName = fileName
	familyFeatures = []
	for item in seqs:
		header, seq = item
		isFeatures = {}
		isFeatures['isName'] = header.split(maxsplit=1)[0]
		isFeatures['tpase'] = seq
		familyFeatures.append(isFeatures)
	mfamilyFeatures.append((familyName, familyFeatures, tpaseSeqFile))
	return mfamilyFeatures
	
# format one sequence into fasta
def writeTpaseSeq2fileOnStream(member):
	isName = member['isName']
	seq = member['seq']
	fastaSeq = '\n'.join(tools.chunkstring(seq, constants.fastaLineWidth))
	headline = '>' + isName
	return [headline, fastaSeq]

# Add classification info such as family/group name into each transposase in mCluster
# Return mClusterNew
# mClusterNew: [cluster, ...]
# cluster: [clusterName, members]
# members: [member, ...]
# member: {'familyName':familyName, 'isName': isName, ..., 'tpase': pepSeq}
#
# blastdb: blast database created by makeblastdb
def classifyTpase(mCluster, blastdb=''):
	# Get the representative tpase of each cluster and combine all representative tapses into one fasta file/string
	tpaseFasta = []
	for cluster in mCluster:
		clusterName, members = cluster
		# Get sequence of the representative tpase of each cluster
		# represent: {'isName': isName, ..., 'tpase': seq}
		represent = members[0]
		tpaseSeqId = represent['isName']
		seq = represent['tpase']
		fastaSeq = '\n'.join(tools.chunkstring(seq, constants.fastaLineWidth))
		#headline = '>' + '_'.join([clusterName,tpaseSeqId])
		headline = '>' + tpaseSeqId
		tpaseFasta.extend([headline, fastaSeq])
	tpaseFile = '\n'.join(tpaseFasta)

	query = tpaseFile
	evalue = constants.min4evalue 
	print('Start blastping {} against {} at {}'.format('representative_tpases', blastdb, 
		datetime.datetime.now().ctime()))
	blastOut, err = tools.doBlastpOnStream(query, blastdb, task='blastp', e_value=evalue, 
			nthreads=constants.nthread)
	if len(err) > 0:
		e = 'Blastp {} against {}: {}'.format('representative_tpases', blastdb, err)
		raise RuntimeError(e)

	print('Start reading blastp output for cluster {} at {}'.format(
		clusterName, datetime.datetime.now().ctime()))
	# Read blast result
	# hits: [hit, ...]
	# hit: {'qseqid':qseqid, 'sseqid':sseqid, 'pident':pident, 'length':length, ..., 
	#		'evalue':evalue, 'nident':nident, 'qlen':qlen, 'slen':slen}
	# qseqid: e.g. 'protein:vir:93615'
	# sseqid: e.g. 'IS3|IS51|ISEC25|'
	hits = tools.getBlastpResultOnStream(blastOut)
	print('Finish reading output of blastp at {}'.format(datetime.datetime.now().ctime()))
	print('Number of clusters with IS family assignment:', len(hits))
	print('Number of clusters without IS family assignment:', len(mCluster)-len(hits))
	return hits

def addClassification2cluster(mCluster, hits):
	# Assign family information to each cluster
	mClusterNew = []
	for cluster in mCluster:
		clusterName, members = cluster
		represent = members[0]
		tpaseSeqId = represent['isName']
		familyname = 'new'
		for hit in hits:
			if tpaseSeqId == hit['qseqid']:
				# hit['sseqid']: e.g. 'IS3|IS51|ISEC25|'
				familyname = hit['sseqid'].split('|')[0]
				break
		membersNew = []
		for member in members:
			member['familyName'] = familyname
			membersNew.append(member)

		# Change cluasterName in according to the future use by the following procedures such as pred.py:
		# aclameFile_clusterIndex => ISfamily_clusterIndex
		# e.g. aclame_proteins_all_0.4.fasta.refined.final_0 => IS1_0
		clusterName = '_'.join([familyname, clusterName.rsplit('_',maxsplit=1)[-1]])

		mClusterNew.append([clusterName, membersNew])
	return mClusterNew

# The smallest tpase domain is assumed to be 50-amino-acid long.
# The smallest tpase protein (single domain) is assumed to be 50-amino-acid long.
#MIN4TPASE_DOMAIN = 50
#MIN4TPASE_DOMAIN = 100
def refineTpaseData(seqfile):
	seqs = tools.getFastaFull(seqfile)
	fasta = []
	for item in seqs:
		header, seq = item
		len4seq = len(seq)
		'''
		if len4seq < MIN4TPASE_DOMAIN:
			continue
		elif len4seq > 3000:
			print(len4seq, header)
		'''
		headerL = header.lower()
		#
		# Transposase related words and words with spelling mistakes:
		# IS (IS family), insert(insertion sequence, insertion element),
		#	ransposas(transposase) and transpoase, ranspositio(transposition), 
		#
		# Bad words which are not related to transposase: ransposo(transposon),
		#
		# Negative words excluding transposase:
		# helper (helper protein), accessory (accessory protein), binding(ATP-binding, atp binding protein),
		#	resolvase, recombinase, repressor(repressor protein), inactiv(inactivated derivative)
		#
		if ((' IS' in header or 
			'insert' in headerL or 
			'ransposas' in headerL or 'transpoase' in headerL or 
			'ranspositio' in headerL 
			# or 'ransposo' in headerL
			)
				and ('helper' not in headerL 
					and 'accessory' not in headerL 
					and 'binding' not in headerL
					and 'resolvase' not in headerL
					and 'recombinase' not in headerL
					and 'repressor' not in headerL # remove proteins with 'repressor' to get version1 of transposase data set
					and 'inactiv' not in headerL # remove more proteins to get verstion2 of transposase data set
					)):

			# get the fake tpase IDs from the curated file with negative tpase list file, *.negative
			with open('/data2/zhiqxie/aclame/aclame_proteins_all_0.4.fasta.refined.negative', 'r') as fp:
				seqids = []
				for line in fp:
					line = line.strip()
					seqid = line[1:].split(maxsplit=1)[0]
					seqids.append(seqid)

			seqid = header.split(maxsplit=1)[0]
			if seqid in seqids:
				continue
			fasta.append(tools.fastaFormat(header, seq))

			# re-examine candidates related to ' IS'(IS family) or 'ransposo'(transposon)
			# After examined all candiates related to transposon, it confirmed that 'transposon'
			#	was not directly related to transposase and it usually related to 'resolvase'.
			if (not ('insert' in headerL or
					'ransposas' in headerL or 'transpoase' in headerL or
					'ranspositio' in headerL)): 
				print('>'+header)
			'''
			if len4seq < MIN4TPASE_DOMAIN:
				print('short tpase',header)
			'''

	#seqfileRefined = '.'.join([seqfile, 'refined', str(MIN4TPASE_DOMAIN)])
	seqfileRefined = '.'.join([seqfile, 'refined.final'])
	#seqfileRefined = '.'.join([seqfile, 'refined.test'])
	with open(seqfileRefined, 'w') as fp:
		fp.write('\n'.join(fasta)+'\n')

def buildTpaseHMM(args):
	tpaseSeqFile = args['tpaseSeqFile']
	rootpath = os.path.join(os.path.dirname(tpaseSeqFile), 'buildhmm')

	# for pre-processing aclame fasta dat to retrieve the real tpase sequences from aclame data set
	#refineTpaseData(tpaseSeqFile)
	#return 0

	logFile = '.'.join([args['program'], 'log'])
	fp_log = open(logFile, 'w')
	cmdline = args['cmdline']
	print(cmdline, datetime.datetime.now().ctime())
	fp_log.write(cmdline + ' ' + datetime.datetime.now().ctime()+'\n')

	# get tranposase sequences from fasta file
	seqs = tools.getFastaFull(tpaseSeqFile)
	# seqs: seqs: [(header, seq), ..., (header, seq)]
	# seq: character string
	mfamilyFeatures = getFamilyFeatures(tpaseSeqFile)

	print('Start clustering transposases in {} at {}'.format(tpaseSeqFile, datetime.datetime.now().ctime()))
	# Based on Tpase peptide sequence, do cluster analysis on IS elements within same family
	# in order to divide the diverse family into multiple subfamilies or groups. 

	# mfamilyFeatures: [(familyName, familyFeatures, familyTpaseSeqFile), ...]
	# familyFeatures: [isFeatures, ...]
	# isFeatures: {'isName': isName, ..., 'seq': seq}
	# seq: character string
	#
	# mCluster: [cluster, ..., cluster]
	# cluster: (clusterName, members)
	# members: [isFeatures, ..., isFeatures]
	# isFeatures: {'isName': isName, 'clusterName':clusterName, ..., 'seq': seq}
	#
	if LOW_IDENT == False:
		print('Clustering each family at sequence identity 60%')
	elif LOW_IDENT == True:
		print('Clustering each family at sequence identity 30%')
	mCluster = getCluster(mfamilyFeatures)

	print('Finish clustering transposases in {} into {} clusters at {}'.format(
		tpaseSeqFile, len(mCluster), datetime.datetime.now().ctime()))

	# Query the representative member (transposase) of each cluster against IS faimly/group name database to 
	# assign IS classification to each cluster of the transposases in original data set 
	# (transposase sequences).
	# mCluster: [cluster, ..., cluster]
	# cluster: (clusterName, members)
	# members: [member, ...]
	# member: {'isName': isName, ..., 'tpase': pepSeq}
	blastdb = BLASTDB
	if len(blastdb) > 0:
		hits = classifyTpase(mCluster, blastdb)
	else: # no classification database is available
		hits = []
	mClusterNew = addClassification2cluster(mCluster, hits)
	# mClusterNew: [cluster, ...]
	# cluster: (clusterName, members)
	# members: [member, ...]
	# member: {'familyName':familyName, 'isName': isName, ..., 'tpase': pepSeq}
	print('Finish classifying transposases in {} clusters at {}'.format(
		len(mClusterNew), datetime.datetime.now().ctime()))

	args2concurrent4seq = []

	DOphmmer = False
	DOhmmsearch = False
	fasta4singleMemberClusters = []
	for cluster in mClusterNew:
		clusterName, members = cluster
		# the first member of cluster
		familyname = members[0]['familyName']

		familyNameOnDisk = familyname.replace('/', '_')
		family_path = os.path.join(rootpath, familyNameOnDisk)

		clusterNameOnDisk = clusterName.replace('/', '_')

		seq_filename = '.'.join([clusterNameOnDisk, 'faa'])
		msa_filename = '.'.join([clusterNameOnDisk, 'faa.msa'])
		seq_pathfilename = os.path.join(family_path, seq_filename)
		msa_pathfilename = os.path.join(family_path, msa_filename)

		# output multiple sequences in fasta format into a file, one file each family or cluster
		cluster_fasta = []
		nseq = 0
		for element in members:
			if 'tpase' in element.keys():
				pepSeq = element['tpase']
			else:
				continue
			if len(pepSeq) == 0:
				continue

			familyname = element['familyName']
			groupname = '' # we do not have subgroup name
			elementname = element['isName']
			header = '|'.join((clusterName, familyname, groupname, elementname))
			element_fasta = tools.fasta_format(header, pepSeq)

			cluster_fasta.append(element_fasta)
			nseq += 1

		if nseq > 1:
			# cluster with two or more members
			args2concurrent4seq.append((seq_pathfilename, msa_pathfilename))
			if DOhmmsearch == False:
				DOhmmsearch = True
		elif nseq == 1:
			# single-member cluster
			# 1) multiple sequence alignement can not be done
			# 2) phmmer instead of hmmsearch will be used for searching a single query
			# against a sequence database. phmmer works essentially just like hmmsearch
			# does, except we provide a query sequence instead of a query profile HMM.

			fasta4singleMemberClusters.append(cluster_fasta[0])
			if DOphmmer == False:
				DOphmmer = True

		#if not os.path.isfile(seq_pathfilename):
		tools.makedir(family_path)
		with open(seq_pathfilename, 'w') as fp_seq_pathfilename:
			fp_seq_pathfilename.write('\n'.join(cluster_fasta))
		print('Write peptide sequences from cluster {} of family {} into {}'.format(
				clusterName, familyname, seq_pathfilename))

	# sequence file holding multiple peptide sequences used as input for phmmer
	clustersSeqFile4phmmer = os.path.join(rootpath, '.'.join(['clusters', 'single', 'faa']))
	with open(clustersSeqFile4phmmer, 'w') as fp_clustersSeqFile4phmmer:
		fp_clustersSeqFile4phmmer.write('\n'.join(fasta4singleMemberClusters))
	print('Write peptide sequences from single-member clusters into', clustersSeqFile4phmmer)

	# Build multiple alignment and HMM models
	hmms_file = os.path.join(rootpath, '.'.join(['clusters', 'faa', 'hmm']))
	print('Building multiple alignment and HMM models.')
	buildMSA(args2concurrent4seq, constants.nproc, fp_log)
	buildHMM(args2concurrent4seq, hmms_file, fp_log)

	fp_log.close()
	print('Log into {} at {}'.format(logFile, datetime.datetime.now().ctime()))


if __name__ == "__main__":
	# Parse command line arguments
	descriptStr = "Build profile hidden markov model for each insertion sequence cluster (family) from transposase sequences. A typical invocation would be: python3 buildProfileHMM.py tpase.fasta"
	parser = argparse.ArgumentParser(description = descriptStr)

	helpStr = 'file holding the sequences of transposases, fasta format'
	parser.add_argument('tpaseSeqFile', help = helpStr)

	args = parser.parse_args()
	args4buildTpaseHMM = {
				'cmdline': ' '.join(sys.argv),
				'program': sys.argv[0],
				'tpaseSeqFile': args.tpaseSeqFile,
				}

	buildTpaseHMM(args4buildTpaseHMM)

