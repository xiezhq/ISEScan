import subprocess, shlex
import string
import operator
import re
import constants
import tools
import concurrent.futures
import os.path
import itertools
import ssw_wrap
import sys
import datetime


# Return mInput4ssw
## mInput4ssw: [input4orfHits, ..., input4orfHits]
## input4orfHits: [input4orfhit, ..., input4orfhit]
# mInput4ssw: [input4orfhit, ..., input4orfhit]
# input4orfhit: (familyName, orfStr, seq1, seq2, minScore, minLen)
# familyName: the profile HMM model of which is used to search protein database
# orfStr: string, orf == ('0000941.1', 20, 303, '+') -> orfStr == '0000941.1_20_303_+'
#
# mOrfHits: {seqid: orfHits, ..., seqid: orfHits}
# orfHits: [orfhit, ..., orfhit]
# orfhit: (orf, familyName, best_1_domain_E-value, full_sequence_E-value, overlap_number)
# orf: (seqid, begin, end, strand)
# mDna: {seqid: (org, fileid, sequence), ..., seqid: (org, fileid, sequence)}
# sequence: string
#
def prepare4ssw2findIRbyDNAbyFar(mOrfHits, mDna):
	mInput4ssw = []
	mboundary = {}
	for seqid in mOrfHits:
		if len(mOrfHits[seqid]) == 0:
			continue
		DNAlen = len(mDna[seqid][-1])
		for orfHit in mOrfHits[seqid]:
			familyName = orfHit[1]

			maxDist4ter2orf = constants.maxDist4ter2orf
			minDist4ter2orf = constants.minDist4ter2orf

			orf = orfHit[0]
			orfLeftPos, orfRightPos, strand = orf[1:]
			start1, end1, start2, end2 = pseudoSeqBoundary_v4(orfLeftPos, orfRightPos, 
									maxDist4ter2orf, minDist4ter2orf)
			# Set start as the first bp if distance between left end and the nearest ORF is 
			# less than maxDist4ter2orf.
			if start1 < 1:
				start1 = 1
			if end2 > DNAlen:
				end2 = DNAlen

			# The strand does not matter when extract two terminal sequences to align, namely,
			# which sequence is the first sequence in pairwise alignement does not make sense.
			lSeq = mDna[seqid][-1][start1-1: end1]
			rSeq = tools.complementDNA(mDna[seqid][-1][start2-1: end2], '1')[::-1]

			minScore = 0.0

			# familyName example: 'IS200/IS605_4|IS200/IS605|IS1341|ISCARN12|', 
			#	'IS30_0', 'IS110_25|IS110||ISLIN1|'
			min4tir = constants.min4tir
			#min4tir = constants.optmin4tir
			if '|' in familyName:
				familyCluster = familyName.split('|',1)[0]
			else:
				familyCluster = familyName
			family, cluster = familyCluster.rsplit('_', 1)
			if familyCluster == 'IS200/IS605_8':
				minLen = min4tir[familyCluster]
			else:
				minLen = min4tir[family]

			# transformat orf to a string, for example,
			# orf == ('gi|256374160|ref|NC_013093.1|', 20, 303, '+') 
			#	-> orfStr == 'gi|256374160|ref|NC_013093.1|_20_303_+'
			# So, orfStr can be used an ID like isName in other called sub-functions
			orfStr = '_'.join((orf[0], str(orf[1]), str(orf[2]), orf[3]))
			mInput4ssw.append((familyName, orfStr, lSeq, rSeq, minScore, minLen))
			mboundary[orfStr] = (start1, end1, start2, end2)
		#mInput4ssw.append(input4orfHits)
	return (mInput4ssw, mboundary)

# check the intersection between the current upsteam and downstream tir search windows
# and the neighboring tpase ORFs, and then probably shrink the tir search windows to
# avoid the insersection with the neighboring tpase ORFs.
#
# orfhitsNeighbors: {orf: orfhitneighbors, ..., orfhitneighbors}
# orfhitneighbors: [before, orfhit, after], before/after is None or orfhit
# orfhit: (orf, familyName, best_1_domain_E-value, full_sequence_E-value, ncopy4tpase, raworfhits)
# orf: (seqid, begin, end, strand)
def tirwindowIntersectORF(start1, end1, start2, end2, 
		orfhit,
		orfhitsNeighbors, minDist4ter2orf):
	orf = orfhit[0]
	orfBegin, orfEnd = orf[1:3]
	before = orfhitsNeighbors[orf][0]
	after = orfhitsNeighbors[orf][1]
	start1New, end1New, start2New, end2New = start1, end1, start2, end2
	if before != None and orf[3] == before[0][3] and start1 <= before[0][2]:
		start1New = before[0][2] + 1
		#end1New = end1 - (start1-start1New)
		print('hello, before', before)
		print('shrink the boundary of tir search region around orfhit {}, start1 ({}) to {}'.format(
			orfhit, start1,  start1New))

	if after != None and orf[3] == after[0][3] and end2 >= after[0][1]:
		end2New = after[0][1] - 1
		#start2New = start2 - (end2-end2New)
		print('hello, after', after)
		print('shrink the boundary of tir search region around orfhit {}, end2 ({}) to {}'.format(
			orfhit, end2,  end2New))
	return (start1New, end1New, start2New, end2New)

# mDNA: {accid: (org, fileid, sequence), ..., accid: (org, fileid, sequence)}
#
# morfhits: {seqid: orfhits, ...}
# orfhits: [orfhit, ...]
# orfhit: (orf, familyName, best_1_domain_E-value, full_sequence_E-value, ncopy4tpase, raworfhits)
# orf: (seqid, begin, end, strand)
# raworfhits: {'orfhits4tpase':orfhits4tpase}
# orfhits4tpase: [] or [orfhit4tpase, ...]
# orfhit4tpase: (orf, familyName, best_1_domain_E-value, full_sequence_E-value, ncopy4tpase)
#
# orfhitsNeighbors: {orf: orfhitneighbors, ..., orfhitneighbors}
# orfhitneighbors: [before, orfhit, after], before/after is None or orfhit
#
def prepare4ssw2findIRbyDNAbyFar4orfhits(morfhits, mDna, maxDist4ter2orf, minDist4ter2orf,
		morfhitsNeighbors):
	mInput4ssw = []
	mboundary = {}
	for seqid, orfhits in morfhits.items():
		if len(orfhits) == 0:
			continue
		orfhitsNeighbors = morfhitsNeighbors[seqid]
		DNAlen = len(mDna[seqid][-1])
		for orfhit in orfhits:
			orfHit = orfhit
			familyName = orfHit[1]
			orf = orfhit[0]
			orfBegin, orfEnd = orf[1:3]
			orfLen = orfEnd - orfBegin + 1

			# familyName example: 'IS200/IS605_4|IS200/IS605|IS1341|ISCARN12|', 
			#	'IS30_0', 'IS110_25|IS110||ISLIN1|'
			if '|' in familyName:
				familyCluster = familyName.split('|',1)[0]
			else:
				familyCluster = familyName
			family, cluster = familyCluster.rsplit('_', 1)

			#if familyCluster == 'IS200/IS605_8':
			#	minMax4tir = constants.minMax4tir[familyCluster]
			#else:
			#	minMax4tir = constants.minMax4tir[family]
			minMax4tir = constants.minMax4tir[family]

			if constants.useOPTtir == True:
				minLen = minMax4tir[2]
			else:
				minLen = minMax4tir[0]

			ncopy4tpase = orfHit[4]
			virtualORF = False
			if ncopy4tpase > 1:
				virtualORF = True
				# consider the identified aligned-region as a virtual ORF
				qstart, qend = orfBegin, orfEnd

				# if aligned region spans two or more Tpases, eg. composite transposon,
				# the current hit in the aligned region is trimmed to the tpase ORF.
				if constants.splitAlign2orf == True:
					if ((before != None and qstart <= before[0][2]) or 
							(after != None and qend >= after[0][1])):
						print('split ncopy4is={} qstart={} qend={} orf={} before={} after={}'.format(
							ncopy4is, qstart, qend, orf, before, after))
						virtualORF = False
			if virtualORF == True:
				dist = minMax4tir[1]
				#dist = minMax4tir[1] + minMax4dr[1]
				#start1, end1, start2, end2 = qstart, orfBegin, orfEnd, qend
				maxDist = 0
				minDist = -dist
				start1, end1, start2, end2 = pseudoSeqBoundary_v4(qstart, qend, maxDist, minDist)
			else:
				start1, end1, start2, end2 = pseudoSeqBoundary_v4(orfBegin, orfEnd, 
									maxDist4ter2orf, minDist4ter2orf)
			'''
			if ncopy4tpase <= 1:
				start1, end1, start2, end2 = tirwindowIntersectORF(
						start1, end1, start2, end2, 
						orfhit, orfhitsNeighbors, minDist4ter2orf)
			'''

			'''
			if end1 >= start2:
				end1 = int((end1+start2)/2)
				start2 = end1 + 1
			'''

			# check DNA termini to assure that tir is within DNA termini
			if start1 < 1:
				start1 = 1
				if start1 > end1:
					end1 = 1
			if end2 > DNAlen:
				end2 = DNAlen
				if start2 > end2:
					start2 = DNAlen

			# check the sequencial order of tir search boundaries and set the tir search window of 1 bp
			# to prevent SSW alignment from returning fake tir.
			if not (start1 < end1 < start2 < end2):

				#e = 'Error, invalid tir search window (org={} fastafile={} seq={}): {}-{} {}-{} around ORF {}'.format(
				#		mDna[seqid][0], mDna[seqid][1], seqid, start1, end1, start2, end2, orf)
				#raise RuntimeError(e)
				e = 'No tir will be found in the invalid tir search window (org={} fastafile={}): {}-{} {}-{} around ORF {}'.format(
						mDna[seqid][0], mDna[seqid][1], start1, end1, start2, end2, orfhit)
				print(e)

				# Set a tir search windows the 1-bp-long termini of ORF, 
				# which will later prevent alignment algorithm from returning any tir.
				start1 = end1 = orf[1]
				start2 = end2 = orf[2]

			# The strand does not matter when extracting two terminal sequences to align, namely,
			# which sequence is the first sequence in pairwise alignement does not make sense.
			lSeq = mDna[seqid][-1][start1-1: end1]
			rSeq = tools.complementDNA(mDna[seqid][-1][start2-1: end2], '1')[::-1]

			minScore = 0.0

			# transformat orf to a string, for example,
			# orf == ('gi|256374160|ref|NC_013093.1|', 20, 303, '+') 
			#	-> orfStr == 'gi|256374160|ref|NC_013093.1|_20_303_+'
			# So, orfStr can be used an ID like isName in other called sub-functions
			orfStr = '_'.join((orf[0], str(orf[1]), str(orf[2]), orf[3]))
			mInput4ssw.append((familyName, orfStr, lSeq, rSeq, minScore, minLen))
			mboundary[orfStr] = (start1, end1, start2, end2)
	return (mInput4ssw, mboundary)

		
# Prepare sequences for alignment to search for multipe-copy full-length IS elements.
def prepare4ssw2findIScopyByDNA(hits, dna):
	mInput4ssw = []
	mboundary = {}
	if len(hits) == 0:
		return (mInput4ssw, mboundary)

	hitPairs = itertools.combinations(hits, 2)

	'''
	for hitpair in hitPairs:
		input4ssw, boundary = prepare4ssw2findIScopyByDNA4hitPair((hitpair, dna))
		mInput4ssw.append(input4ssw)
		mboundary[boundary[0]] = boundary[1]
	'''
	nthread = constants.nthread
	# for +/+ or -/- strands comparison of each hit pair
	with concurrent.futures.ThreadPoolExecutor(max_workers = nthread) as executor:
		future2pairs = {executor.submit(prepare4ssw2findIScopyByDNA4hitPair, (pair, dna)): pair for pair in hitPairs}
		for future in concurrent.futures.as_completed(future2pairs):
			pair = future2pairs[future]
			try:
				input4ssw, boundary = future.result()
			except Exception as exc:
				print('{} generated an exception: {}'.format(pair, exc))
			else:
				# for +/+ alignment
				# seq2 is the normal strand of the 2nd sequence.
				#
				mInput4ssw.append(input4ssw)

				# input4ssw: [tirSeqs, orfStr, seq1, seq2, minScore, minLen]
				# for +/- alignment
				# seq2 is the reversed complementary sequence of the 2nd sequence.
				# seq2 = tools.complementDNA(dna[start2-1: end2], '1')[::-1]
				#
				input4sswCopy = input4ssw[:]
				input4sswCopy[3] = tools.complementDNA(input4sswCopy[3], '1')[::-1]
				mInput4ssw.append(input4sswCopy)

				mboundary[boundary[0]] = boundary[1]

	return (mInput4ssw, mboundary)

# Search ORF copy, then we can check if TIRs in two copies are same or not. 
# Attach tir sequences info to the returned data.
#def prepare4ssw2findIScopyByDNA4hitPair((hitPair, dna)):
def prepare4ssw2findIScopyByDNA4hitPair(input):
	hitPair, dna = input

	'''
	# Option1:
	# 	1) Search full IS element copy, then we can check if TIRs in two copies are same or not.
	# Use the boundary of the first TIR of hit as the boundary of full-length IS element 
	# if multiple TIRs are available, else use the boundary of the ORF of hit.
	if len(hitPair[0]['tirs']) > 0 and len(hitPair[0]['tirs'][0]) > 0:
		start1 = hitPair[0]['tirs'][0][-6]
		end1 = hitPair[0]['tirs'][0][-3]
		tirSeqs = list(hitPair[0]['tirs'][0][-2:])
	else:
		start1, end1 = hitPair[0]['orf'][1:3]
		tirSeqs = ['', '']
	if len(hitPair[1]['tirs']) > 0 and len(hitPair[1]['tirs'][0]) > 0:
		start2 = hitPair[1]['tirs'][0][-6]
		end2 = hitPair[1]['tirs'][0][-3]
		tirSeqs.extend(hitPair[1]['tirs'][0][-2:])
	else:
		start2, end2 = hitPair[1]['orf'][1:3]
		tirSeqs.extend(['',''])
	'''
	# Option2:
	# 	2) Search ORF copy, then we can check if TIRs in two copies are same or not. 
	start1, end1 = hitPair[0]['orf'][1:3]
	start2, end2 = hitPair[1]['orf'][1:3]
	if len(hitPair[0]['tirs']) > 0 and len(hitPair[0]['tirs'][0]) > 0:
		tirSeqs = list(hitPair[0]['tirs'][0][-2:])
	else:
		tirSeqs = ['', '']
	if len(hitPair[1]['tirs']) > 0 and len(hitPair[1]['tirs'][0]) > 0:
		tirSeqs.extend(hitPair[1]['tirs'][0][-2:])
	else:
		tirSeqs.extend(['',''])

	# Use orf of the hit as IS element ID
	orf1 = hitPair[0]['orf']
	orf2 = hitPair[1]['orf']
	#
	# OptionA: IS element is strand-specific DNA element.
	#	1) compare only '+' strands, +/+ or -/-
	seq1 = dna[start1-1: end1]
	seq2 = dna[start2-1: end2]
	#
	# OptionB: IS element is strand-nonspecific DNA element
	#	2) compare both '+' strands and '-' strands, namely, (+/+ and +/-) or (-/- and +/+)
	#seq1 = dna[start1-1: end1]
	#seq21 = dna[start2-1: end2]
	#seq22 = tools.complementDNA(seq21, '1')[::-1]

	# Option0:
	#	minLen = min(seq1, seq2)/2
	#minScore, minLen = 0.0, int(min(len(seq1), len(seq2))/2)
	#
	# Option1:
	#	minLen = min(seq1, seq2) * sim4iso
	minScore, minLen = 0.0, int(min(len(seq1),len(seq2)) * constants.sim4iso)
	#
	# Option2:
	#	minLen = min(orf1Len, orf2Len) * sim4iso
	#orf1Len = hitPair[0]['orf'][2] - hitPair[0]['orf'][1] + 1
	#orf2Len = hitPair[1]['orf'][2] - hitPair[1]['orf'][1] + 1
	#minScore, minLen = 0.0, int(min(orf1Len, orf2Len) * constants.sim4iso)

	# transformat orf to a string, for example,
	# orf == ('0000941.1', 20, 303, '+') -> orfStr == '0000941.1_20_303_+'
	# So, orfStr can be used an ID like isName in other called sub-functions
	orfStr1 = '_'.join((orf1[0], str(orf1[1]), str(orf1[2]), orf1[3]))
	orfStr2 = '_'.join((orf2[0], str(orf2[1]), str(orf2[2]), orf2[3]))
	orfStr = orfStr1 + '|' + orfStr2

	# transformat tir aligned sequences to a string, for example,
	# tir1 for orfstr1
	# CCCTTTGCTTCGCAAAGGCCCTCTC
	# |||*||||||||| |||||||||||
	# CCCCTTGCTTCGC-AAGGCCCTCTC
	# tir2 for orfstr2
	# CCCTTTGCTTCGC
	# |||*|||||||||
	# CCCCTTGCTTCGC
	#
	# tirSeqs = 'CCCTTTGCTTCGCAAAGGCCCTCTC|CCCCTTGCTTCGC-AAGGCCCTCTC|CCCTTTGCTTCGC|CCCCTTGCTTCGC'
	# or tirSeqs = '||CCCTTTGCTTCGC|CCCCTTGCTTCGC'
	# or tirSeqs = 'CCCTTTGCTTCGCAAAGGCCCTCTC|CCCCTTGCTTCGC-AAGGCCCTCTC||'
	# So, tirSeqs can be used to find full IS copies with identical tir sequences.
	#tirSeqs = '|'.join(tirSeqs)
	#
	# tirSeqs = 'CCCTTTGCTTCGCAAAGGCCCTCTC|CCCCTTGCTTCGCAAGGCCCTCTC|CCCTTTGCTTCGC|CCCCTTGCTTCGC'
	# or tirSeqs = '||CCCTTTGCTTCGC|CCCCTTGCTTCGC'
	# or tirSeqs = 'CCCTTTGCTTCGCAAAGGCCCTCTC|CCCCTTGCTTCGCAAGGCCCTCTC||'
	tirSeqs = '|'.join(tirSeqs).replace('-', '')
	
	# Data items correspondence bwteen TIR searching and IS copy searching. 
	# familyName, isName = tirSeqs, orfStr, in order that isName, namely, orfStr, is an unique ID

	input4ssw = [tirSeqs, orfStr, seq1, seq2, minScore, minLen]
	boundary = ((orf1, orf2), (start1, end1, start2, end2))
	return (input4ssw, boundary)

# Based on comparing two full IS elements with TIRs, find the copies of IS element in same DNA sequence.
def prepare4ssw2findIScopyByDNA4hitPairByTIR(input):
	hitPair, dna = input
	# Use the boundary of the first TIR of hit as the boundary of full-length IS element if TIR is available,
	# else use the boundary of the ORF of hit.
	if len(hitPair[0]['tirs']) > 0 and len(hitPair[0]['tirs'][0]) > 0:
		start1 = hitPair[0]['tirs'][0][-6]
		end1 = hitPair[0]['tirs'][0][-3]
	else:
		start1, end1 = hitPair[0]['orf'][1:3]
	if len(hitPair[1]['tirs']) > 0 and len(hitPair[1]['tirs'][0]) > 0:
		start2 = hitPair[1]['tirs'][0][-6]
		end2 = hitPair[1]['tirs'][0][-3]
	else:
		start2, end2 = hitPair[1]['orf'][1:3]

	# Use orf of the hit as IS element ID
	orf1 = hitPair[0]['orf']
	orf2 = hitPair[1]['orf']
	if orf1[3] == '+':
		seq1 = dna[start1-1: end1]
	else:
		seq1 = tools.complementDNA(dna[start1-1: end1], '1')[::-1]
	if orf2[3] == '+':
		seq2 = dna[start2-1: end2]
	else:
		seq2 = tools.complementDNA(dna[start2-1: end2], '1')[::-1]

	# Option1:
	#	minLen = min(seq1, seq2)/2
	#minScore, minLen = 0.0, int(min(len(seq1), len(seq2))/2)
	#
	# Option2:
	#	minLen = min(seq1, seq2) * sim4iso
	#minScore, minLen = 0.0, int(min(len(seq1),len(seq2)) * constants.sim4iso)
	#
	# Option3:
	#	minLen = min(orf1Len, orf2Len) * sim4iso
	orf1Len = hitPair[0]['orf'][2] - hitPair[0]['orf'][1] + 1
	orf2Len = hitPair[1]['orf'][2] - hitPair[1]['orf'][1] + 1
	minScore, minLen = 0.0, int(min(orf1Len, orf2Len) * constants.sim4iso)

	# transformat orf to a string, for example,
	# orf == ('0000941.1', 20, 303, '+') -> orfStr == '0000941.1_20_303_+'
	# So, orfStr can be used an ID like isName in other called sub-functions
	orfStr1 = '_'.join((orf1[0], str(orf1[1]), str(orf1[2]), orf1[3]))
	orfStr2 = '_'.join((orf2[0], str(orf2[1]), str(orf2[2]), orf2[3]))
	orfStr = orfStr1 + '|' + orfStr2
	# familyName, isName = orfStr, orfStr, in order that isName, namely, orfStr, is an unique ID

	input4ssw = (orfStr, orfStr, seq1, seq2, minScore, minLen)
	boundary = ((orf1, orf2), (start1, end1, start2, end2))
	return (input4ssw, boundary)


# Note: finding IR requires sharp values for gap and mismatch comparing with finding homology sequence
#	as the two short sequences of IR usually have much higher identity than common homology sequences
# Return: filters = [filter, .., filter]
# filter: (gap, gapextend, match, mismatch)
def buildFilter4ssw(gap, gapextend, match, mismatch):
	# gapopen: default 2 (refer to consistants.py)
	gapRange = range(gap - 1, gap + 12, 1)
	#gapRange = range(gap - 1, gap + 8, 1)
	#gapRange = range(gap - 1, gap + 4, 1)
	#gapRange = range(gap - 0, gap + 2, 1)

	# gapextend: default 1
	gapextendRange = range(gapextend - 0, gapextend + 12, 1)
	#gapextendRange = range(gapextend - 0, gapextend + 8, 1)
	#gapextendRange = range(gapextend - 0, gapextend + 4, 1)
	#gapextendRange = range(gapextend - 0, gapextend + 2, 1)

	# match: default 2
	matchRange = range(match - 1, match + 12, 1)
	#matchRange = range(match - 1, match + 8, 1)
	#matchRange = range(match - 1, match + 4, 1)
	#matchRange = range(match - 0, match + 2, 1)

	# mismatch: default 2
	#mismatchRange = range(mismatch - 1, mismatch + 13, 1)
	mismatchRange = range(mismatch - 1, mismatch + 8, 1)
	#mismatchRange = range(mismatch - 1, mismatch + 4, 1)
	#mismatchRange = range(mismatch - 0, mismatch + 2, 1)

	filters = []
	for gap in gapRange:
		for gapextend in gapextendRange:
			for match in matchRange:
				if (gap + gapextend) < match:
					continue
				for mismatch in mismatchRange:
					if mismatch*2 < match:
						continue
					filter = (gap, gapextend, match, mismatch)
					filters.append(filter)
	return filters


# Return group of TIRs with the greatest irScore.
# TIRs: [(TIR1, filter1), ..., (TIRn, filtern)]
#
# elementTIR: [(TIR1, filter1), ..., (TIRn, filtern)], TIR information for one IS element
# TIR: [familyName, isName, ir]
def keepBestTIR_v3(elementTIR):
	for score, g in itertools.groupby(
			sorted(elementTIR, key = lambda x: tools.irScore(x[0][2]), reverse = True), 
			key = lambda x: tools.irScore(x[0][2])):
		return list(g)

# Return the group of TIRs: [(TIR1, filter1), ..., (TIRn, filtern)]
# list(item): [(TIR1, filter1), ..., (TIRn, filtern)]
#
# elementTIR: [(TIR1, filter1), ..., (TIRn, filtern)], TIR information for one IS element
# TIR: [familyName, isName, ir]
def keepBestTIR_v2(elementTIR):
	# sort by nGaps/irLen all possible TIRs for one IS
	for nGaps, g in itertools.groupby(
			sorted(elementTIR, key = lambda x: x[0][2][3]/x[0][2][2]), 
			key = lambda x: x[0][2][3]/x[0][2][2]):
		#return list(g)

		# Sorted and grouped by irId/irLen, then get the TIR with greatest irId/irLen.
		# list(g): [(TIR1, filter1), ..., (TIRn, filtern)]
		for irId, item in itertools.groupby(
				sorted(g, key = lambda x: x[0][2][1]/x[0][2][2], reverse = True), 
				key = lambda x: x[0][2][1]/x[0][2][2]):
			# Return the best TIRs with least nGaps/irLen and greatest irId/irLen
			return list(item)

# Return the group of TIRs with least nGaps: [(TIR1, filter1), ..., (TIRn, filtern)]
# list(item): [(TIR1, filter1), ..., (TIRn, filtern)]
# Note: TIR1 and TIRn have same nGaps and irId but may have different sequence alignments.
# 	For example, ir1 and ir2 have the same nGaps and irId but might have different irLen and 
#	start1/end1/start2/end2 and seq1/seq2.
#
# elementTIR: iterator returned by itertools.groupby, TIR information for one IS element
def keepBestTIR(elementTIR):
	# sort by nGaps all possible TIRs for one IS
	for nGaps, g in itertools.groupby(
			sorted(elementTIR, key = lambda x: x[0][2][3]), 
			key = lambda x: x[0][2][3]):
		return list(g)


# Return bestTIR which is the best matched TIR with the best irScore comparing with other matched TIRs 
# if there is any other matched TIR
# Note: the bestTIR may be found under multiple filter conditions
# 
# TIRfilters: [(TIR_1, filter_1), ..., (TIR_i, filter_j), ...]
# TIR: [familyName, isName, ir]
#
# TIRs: [element1TIR, ..., elementnTIR]
# elementTIR: [(TIR1, filter1), ..., (TIRn, filtern)]
# TIR: [familyName, isName, ir]
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
#
def checkTIRseq(TIRfilters):
	bestTIR = []
	#isNames = []
	elements = []
	for k, g in itertools.groupby(
			# sort and group by isName
			sorted(TIRfilters, key = lambda x: x[0][1]), 
			lambda x: x[0][1]):
		#isNames.append(k)
		elements.append(list(g))
	for elementTIR in elements:

		# keep the best TIRs which have same irScore but might have different alignments under different filters
		#
		elementTIR = keepBestTIR_v3(elementTIR)

		# elementTIR: [(TIR1, filter1), ..., (TIRn, filtern)]
		# Note:
		#	TIR1 and TIRn have same irScore but may have different sequence alignments. For example, 
		#	ir1 and ir2 have the same nGaps and irId and irLen but might 
		#	have different seq1/seq2.
		#
		# get all ISs, each of which keeps only the best TIR group
		if len(elementTIR[0][0][2]) > 0:
			bestTIR.append(elementTIR)

	return bestTIR

# Return bestTIR which is the best matched TIR with least gaps comparing with other matched TIRs if there is any other matched TIR
# Note: the bestTIR may be found under multiple filter conditions
# 
# TIRfilters: [(TIR_1, filter_1), ..., (TIR_i, filter_j), ...]
# TIR: [familyName, isName, ir]
#
# TIRs: [element1TIR, ..., elementnTIR]
# elementTIR: [familyName, isName, ir]
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
#
def checkTIRseq_v1(TIRfilters):
	bestTIR = []
	#isNames = []
	elements = []
	for k, g in itertools.groupby(sorted(TIRfilters, key = lambda x: x[0][1]), lambda x: x[0][1]):
		#isNames.append(k)
		elements.append(list(g))
	for elementTIR in elements:

		# keep the best TIRs with same nGaps and irId and different alignment paths under different filters
		#
		# Given one IS, group by nGaps all possible TIRs and then get best TIR group (same sequence under one or more filters)
		elementTIR = keepBestTIR(elementTIR)
		#

		# elementTIR: [(TIR1, filter1), ..., (TIRn, filtern)]
		# Note: TIR1 and TIRn have same nGaps and irId and irLen but may have different sequence alignments.
		# 	For example, ir1 and ir2 have the same nGaps and irId and irLen but might have different start1/end1/start2/end2 and seq1/seq2.
		#
		# get all ISs, each of which keeps only the best TIR group
		bestTIR.append(elementTIR)

	return bestTIR


# filterPerformance: [(filter, perf, TIRs), ..., (filter, perf, TIRs)]
# filter: (gapopen, gapextend, match, mismatch)
# perf: (nmatch, ndismatch, ndiscard, nIS)
# TIRs: [TIR, ..., TIR]
# TIR: [familyName, isName, ir]
def outputPerformanceBySSW(filterPerformance):
	fmtWaterPerfStr = '{:>7} {:>9} {:>5} {:>8} {:>9} {:>12} {:>9}'
	print(fmtWaterPerfStr.format('gapOpen', 'gapExtend', 'match', 'mismatch', 'matchedIS', 'notMatchedIS', 'discardIS'))
	print('-' * 50)
	for item in filterPerformance:
		gap, gapextend, match, mismatch = item[0]
		matched, notMatched, discard, nIS = item[1]
		print(fmtWaterPerfStr.format(gap, gapextend, match, mismatch, matched, notMatched, discard))

# filterPer: {filter: [nIS, {isName, ..., isName}], ..., filter: [nIS, {isName, ..., isName}]}
def outputPerf_v2(filterPerf):
	fmtPerfStr = '{:>7} {:>9} {:>5} {:>8} {:<}'
	print(fmtPerfStr.format('gapOpen', 'gapExtend', 'match', 'mismatch', 'matchedISwithBestTIR'))
	print('-' * 50)
	for filter, value in sorted(filterPerf.items(), key = lambda x: x[1][0], reverse = True):
		gap, gapextend, match, mismatch = filter
		print(fmtPerfStr.format(gap, gapextend, match, mismatch, value[0]))

# filterPerf: {filter: nIS, ..., filter: nIS}
def outputPerf(filterPerf):
	fmtPerfStr = '{:>7} {:>9} {:>5} {:>8} {:<}'
	print(fmtPerfStr.format('gapOpen', 'gapExtend', 'match', 'mismatch', 'matchedISwithBestTIR'))
	print('-' * 50)
	for filter, num in sorted(filterPerf.items(), key = operator.itemgetter(1), reverse = True):
		gap, gapextend, match, mismatch = filter
		print(fmtPerfStr.format(gap, gapextend, match, mismatch, num))

# find the IS elements which have TIR found under filters other than best filter
#
# filterPer: {filter: [nIS, {isName, ..., isName}], ..., filter: [nIS, {isName, ..., isName}]}
# filter: (gapopen, gapextend, match, mismatch)
#
# bestTIRfilters: [element1TIRgroup, ..., elementnTIRgroup]
# elementTIRgroup: [(TIR, filter), ..., (TIR, filter)]
#
def TIRbyNonbestfilter_v2(filterPerf, bestTIRfilters):
	# get matched TIR info under the best filter
	isNamesByBestFilter = (sorted(filterPerf.values(), key = operator.itemgetter(0)))[-1][1]
	
	# get best TIR info under one or more filters
	isNamesByFilters = {item[0][0][1] for item in bestTIRfilters}
	diffNames = sorted(isNamesByFilters - isNamesByBestFilter)
	print('matched ISs under any filter other than the best filter: {}={}-{}\n{}'.format(	len(diffNames), 
												len(isNamesByFilters), 
												len(isNamesByBestFilter), 
												diffNames))
	print('best TIRs found by filters other than best filter')
	print('output ONE of the filters producing the best TIR of the IS element')
	print('-' * 50)
	for item in bestTIRfilters:
		# item[0][0][1]: isName
		if item[0][0][1] in diffNames:
			# get filters
			#for TIRfilter in item:
			#	gap, gapextend, matrixFile = TIRfilter[1]
			#	match, mismatch = tools.resolveMatrixFileName(matrixFile)
			filters = [TIRfilter[1] for TIRfilter in item]

			# output all filters producing the best TIR of the IS element
			#print('{} {} {}:'.format(item[0][0][0], item[0][0][1], filters))
			#
			# output one of the filters producing the best TIR of the IS element
			print('{} {} {}:'.format(item[0][0][0], item[0][0][1], filters[0]))
			start1, end1, start2, end2, seq1, seq2 = item[0][0][2][-6:]			
			print('{:>6} {} {:<6}\n{:6} {} {:6}\n{:>6} {} {:<6}'.format(
				start1, seq1, end1, 
				' ', tools.buildMatchLine(seq1, seq2), ' ', 
				end2, seq2, start2))

# find the IS elements which have TIR found under filters other than best filter
#
# filterPerformance: [(filter, perf, TIRs), ..., (filter, perf, TIRs)]
# filter: (gapopen, gapextend, match, mismatch)
# perf: (nmatch, ndismatch, ndiscard, nIS)
# TIRs: [TIR, ..., TIR]
# TIR: [familyName, isName, ir]
#
# bestTIRfilters: [element1TIRgroup, ..., elementnTIRgroup]
# elementTIRgroup: [(TIR, filter), ..., (TIR, filter)]
#
def TIRbyNonbestfilter(filterPerformance, bestTIRfilters):
	# get matched TIR info under the best filter
	isNamesByBestFilter = {item[1] for item in filterPerformance[0][2]}
	
	# get best TIR info under one or more filters
	isNamesByFilters = {item[0][0][1] for item in bestTIRfilters}
	diffNames = sorted(isNamesByFilters - isNamesByBestFilter)
	print('matched ISs under any filter other than the best filter: {}={}-{}\n{}'.format(	len(diffNames), 
												len(isNamesByFilters), 
												len(isNamesByBestFilter), 
												diffNames))
	print('matched TIR found by filters other than best filter')
	print('output the first filter producing the best TIR of the IS element')
	print('-' * 50)
	for item in bestTIRfilters:
		# item[0][0][1]: isName
		if item[0][0][1] in diffNames:
			# get filters
			#for TIRfilter in item:
			#	gap, gapextend, matrixFile = TIRfilter[1]
			#	match, mismatch = tools.resolveMatrixFileName(matrixFile)
			filters = [TIRfilter[1] for TIRfilter in item]

			# Option1: output all filters producing the best TIR of the IS element
			print('{} {} {}:'.format(item[0][0][0], item[0][0][1], filters))
			# Option2: output the first filter producing the best TIR of the IS element
			#print('{} {} {}:'.format(item[0][0][0], item[0][0][1], filters[0]))

			start1, end1, start2, end2, seq1, seq2 = item[0][0][2][-6:]			
			matchLine = tools.buildMatchLine(seq1, seq2)
			print('{:>6} {} {:<6}\n{:6} {} {:6}\n{:>6} {} {:<6}'.format(
				start1, seq1, end1,
				' ', matchLine, ' ',
				end2, seq2, start2))


# Combine all matched bestTIRs (bestMatchedTIRfilters) and other bestTIRs which are from bestTIRfilters
# but not found in bestMatchedTIRfilters, return the combined TIRfilters.
# bestTIRfilters: [element1TIRgroup, ..., elementnTIRgroup]
# elementTIRgroup: [(TIR, filter), ..., (TIR, filter)]
# filter: (gapopen, gapextend, match, mismatch)
# TIR: [familyName, isName, ir]
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
def combineBestTIRfilters(bestMatchedTIRfilters, bestTIRfilters):
	TIRfilters = []
	for item in bestTIRfilters:
		isName = item[0][1]
		for itemMatched in bestMatchedTIRfilters:
			isNameMatched = itemMatched[0][1]
			if isNameMatched == isName:
				break
		else:
			TIRfilters.append(item)
	return bestMatchedTIRfilters + TIRfilters


# Calculate TIRs in all IS elements under the specific filter condition
#
#def getPerformanceByFilterBySSW(args2concurrent):
#	mfamilyFeatures, mInput4ssw, filter = args2concurrent
def getPerformanceByFilterBySSW(mfamilyFeatures, mInput4ssw, filter):
	bestIRs = findIRbySSW(mInput4ssw, filter)

	"""
	# Shorten each ir to the reasonable length based on constants.maxLenIR 
	IRs = []
	# bestIRs: [IR, ..., IR]
	# IR: [familyName, isName, ir]
	for IR in bestIRs:
		ir = tools.shortenIR(IR[2])
		IRs.append([IR[0], IR[1], ir])
	"""
	IRs = bestIRs

	# compare IRs found by SSW with IRs found in IS dataset
	#perf, TIRs = compareIRbyISfinder(IRs, mfamilyFeatures)
	#
	# Compare the predicted IRs with IRs from isfinder, and return the IRs matched by records in ISfinder
	perf, matchedTIRs = compareIRbyISfinder_v2(IRs, mfamilyFeatures)

	return (perf, matchedTIRs, IRs)

# Calculate performance for each filter
# Return perf: {filter: [nIS, {isNames, ..., isNames}]}
# nIS: len({isNames, ..., isNames}), number of IS element with the best TIR which has the best irScore 
#	comparing with other TIRs if other TIRs available.
#
# TIRfilters: [element1TIRgroup, ..., elementnTIRgroup]
# elementTIRgroup: [(TIR, filter), ..., (TIR, filter)]
# filter: (gapopen, gapextend, match, mismatch)
# TIR: [familyName, isName, ir]
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
def calculatePerf_v2(TIRfilters):
	perf = {}
	for item in TIRfilters:
		isName = item[0][0][1]
		for tirfilter in item:
			if tirfilter[1] in perf:

				perf[tirfilter[1]][0] += 1
				perf[tirfilter[1]][1].add(isName)
			else:
				value = [1, {isName}]
				perf[tirfilter[1]] = value
	return perf

# Calculate performance for each filter
# Return perf: {filter: nIS}
# nIS: number of IS element with the best TIR which has the best irScore comparing with other TIRs if other TIRs available.
# TIRfilters: [element1TIRgroup, ..., elementnTIRgroup]
# elementTIRgroup: [(TIR, filter), ..., (TIR, filter)]
# filter: (gapopen, gapextend, match, mismatch)
# TIR: [familyName, isName, ir]
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
def calculatePerf(TIRfilters):
	perf = {}
	for item in TIRfilters:
		for tirfilter in item:
			if tirfilter[1] in perf:

				perf[tirfilter[1]] += 1
			else:
				perf[tirfilter[1]] = 1
	return perf


# Return mTIR where only unique tirs among the best TIRs for each IS will be kept, 
# Note:	there are probably multiple different best TIRs (with identical irScore)
#	found by local alignment algorithm.
# mTIR: {isName: (familyName, isName, tirs), ..., isName: (familyName, isName, tirs)}
# tirs: [tir, ..., tir]
# tir: [irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2], ir without score attached
# 
# args for function:
# bestTIRfilters: [element1_TIRgroup, ..., elementn_TIRgroup],
#       best matched TIRs which can be found under one specific filter or multiple different filters
# elementTIRgroup: [(TIR1, filter1), ..., (TIRn, filtern)]
# TIR: [familyName, isName, ir]
# ir: (score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2)
def independentTIR(bestTIRfilters):
	mTIR = {}
	# isTirfilters: multiple tirs for one IS element
	for isTirfilters in bestTIRfilters:
		tirs = []
		familyName, isName = isTirfilters[0][0][:2]
		# check multiple tirs for each IS element
		# tirfilter: one tir
		for tirfilter in isTirfilters:
			# keep alignment score for each ir
			#ir = tirfilter[0][2]
			# remove alignment score for each ir
			tir = tirfilter[0][2][1:]
			tirs.append(tuple(tir))

		# remove duplicate ir for each IS
		tirs = set(tirs)

		mTIR[isName] = (familyName, isName, tirs)
	return mTIR

def independentTIRwithScore(bestTIRfilters):
	mTIR = {}
	# isTirfilters: multiple tirs for one IS element
	for isTirfilters in bestTIRfilters:
		tirs = []
		familyName, isName = isTirfilters[0][0][:2]
		# check multiple tirs for each IS element
		# tirfilter: one tir
		for tirfilter in isTirfilters:
			# keep alignment score for each ir
			tir = tirfilter[0][2]
			tirs.append(tuple(tir))

		# remove duplicate ir for each IS
		tirs = set(tirs)

		mTIR[isName] = (familyName, isName, tirs)
	return mTIR

# Convert boundary numbering in short sequence segment to boundary numbering in original full sequence
# mTIR: {isName: (familyName, isName, tirs), ..., isName: (familyName, isName, tirs)}
# tirs: [tir, ..., tir]
# tir: (irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2), ir without score attached
# 	or
#	(score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2), ir with score attached
def restoreBoundary4tir(mTIR, mboundary):
	new_mTIR = {}
	for isName in mTIR:
		# actual boundary
		bdstart1, bdend2 = mboundary[isName][0], mboundary[isName][-1]

		# pure boundary is used to define tir boundary found by SSW
		# Note: SSW consume two short sequences and define boundaries of two sequences
		#	as 1, Length1 and 1, Length2. When to identify TIR, we use Length1 = Length2 
		p_bdstart1 = 1
		p_bdstart2 = 1

		tirs = (mTIR[isName])[-1]
		new_tirs = []
		for tir in tirs:
			p_start1, p_end1, p_start2, p_end2 = tir[-6:-2]
			# start1 = bdstart1 + move1 and move1 = p_start1 - p_bdstart1
			start1 = bdstart1 + (p_start1 - p_bdstart1)
			# end1 - start1 = p_end1 - p_start1
			end1 = start1 + (p_end1 - p_start1)
			# end2 = bdend2 - move2 and move2 = p_start2 - p_bdstart2
			end2 = bdend2 - (p_start2 - p_bdstart2)
			# end2 - start2 = p_end2 - p_start2
			start2 = end2 - (p_end2 - p_start2)

			if len(tir) == 10:
				# with alignment score 
				new_tir = (tir[0], tir[1], tir[2], tir[3], 
						start1, end1, start2, end2, 
						tir[-2], tir[-1])
			else:
				# without alignment score
				new_tir = (tir[0], tir[1], tir[2], 
						start1, end1, start2, end2, 
						tir[-2], tir[-1])
			new_tirs.append(new_tir)

		new_tirs.sort(key = tools.irScore, reverse = True)
		new_mTIR[isName] = (mTIR[isName][0], mTIR[isName][1], new_tirs)

	return new_mTIR


# Define the boundary of two sub-sequences to be concatenated into one pseudo sequence.
# Note: IR does not always starts from terminus though it is true for most IS elements. 
#
# orf: (left, right), maximal and maximal distances from the left terminus genome sequence to the most left 
#	and the most right ORFs, respectively.
# distIR2Orf: (leftMax, leftMin, rightmin, rightMax), the longest and the shortest distances from the potential IS element 
#	to the most left ORF and those from the most right ORF to the potential IS element, respectively,
#	where the IRs lie out of ORFs. The negative value of left or right means that IR is within ORF.
# 
def pseudoSeqBoundary(orf, distIR2Orf):
	bdStart1 = orf[0] - distIR2Orf[0]
	bdEnd1 = orf[0] - distIR2Orf[1]
	bdStart2 = orf[1] + distIR2Orf[2]
	bdEnd2 = orf[1] + distIR2Orf[3]

	return (bdStart1, bdEnd1, bdStart2, bdEnd2)

# orfLeftPos, orfRightPos: the most left and right positions of ORFs
# irLong: allowed maximal length of IR
def pseudoSeqBoundary_v3(orfLeftPos, orfRightPos, irLong):
	bdEnd1 = orfLeftPos - 1
	bdStart1 = orfLeftPos - irLong
	bdStart2 = orfRightPos + 1
	bdEnd2 = orfRightPos + irLong
	return (bdStart1, bdEnd1, bdStart2, bdEnd2)

# orfLeftPos, orfRightPos: the most left and right positions of ORFs.
# maxDist4ter2orf, minDist4ter2orf: allowed maximum and minimum distance between ends of 
#	IS element and the nearest ORFs in the same IS element.
def pseudoSeqBoundary_v4(orfLeftPos, orfRightPos, maxDist4ter2orf, minDist4ter2orf):
	bdEnd1 = orfLeftPos - minDist4ter2orf
	bdStart1 = orfLeftPos - maxDist4ter2orf
	bdStart2 = orfRightPos + minDist4ter2orf
	bdEnd2 = orfRightPos + maxDist4ter2orf
	return (bdStart1, bdEnd1, bdStart2, bdEnd2)


# Return ir
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# seq1, seq2: 'ATCG', only aligned sequence without gap included
# align: PyAlignRes object for alignment result description, which is returend by ssw_wrap.Aligner.align()
# cigarPair: [(4, 'M'), (2, 'I'), (8, 'M'), (1, 'D'), (10, 'M'), (6, 'S')]
#       Note: for details of cigar sting, Please check the document "The SAM Format Specification",
#       http://samtools.github.io/hts-specs/SAMv1.pdf, particularly the "An example" section and
#       "CIGAR: CIGAR string" section.
def getIRbySSWnoGap(seq1, seq2, align, cigarPair):
	score = align.score
	start1, end1 = align.ref_begin + 1, align.ref_end + 1
	# seq2 is the reverse complementary sequence of right end of full seq
	#start2, end2 = align.query_end + 1, align.query_begin + 1
	start2, end2 = align.query_begin + 1, align.query_end + 1

	nGaps = 0
	irId = 0
	# length of seq2
	irLen = align.query_end - align.query_begin + 1
	index1, index2 = 0, 0
	# process only the aligned sequence area defined by 
	# seq1[align.ref_begin: align.ref_end+1] and seq2[align.query_begin: align.query_end+1]
	seq1 = seq1[align.ref_begin: align.ref_end+1]
	seq2 = seq2[align.query_begin: align.query_end+1]
	for pair in cigarPair:
		if pair[1] == 'I':
			nGaps += pair[0]
			index2 += pair[0]
		elif pair[1] == 'D':
			nGaps += pair[0]
			irLen += pair[0]
			index1 += pair[0]
		elif pair[1] == 'M':
			s1 = seq1[index1: index1+pair[0]]
			s2 = seq2[index2: index2+pair[0]]
			for c1, c2 in zip(s1, s2):
				if c1 == c2:
					irId += 1
			index1 += pair[0]
			index2 += pair[0]
	return [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]

# Get TIR from alignment result
# alignment: tuple returned by tools.buildAlignment()
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
shortestAlignment = 2
def  getIRbySSW(alignment):
	header, line1, line2, line3 = alignment
	# get the alignement
	line1, line2, line3 = line1[20:-9], line2[20:], line3[20:-9] 
	irLen = len(line1)
	if irLen < shortestAlignment: # alignment must be 2 bp or longer.
		return []
	nGap1 = line1.count('-')
	nGap2 = line3.count('-')
	ir = [	header['score'], # score
		line2.count('|'), # irId
		#header['end2'] - header['begin2'] + 1 + nGap2, # irLen
		irLen, # irLen
		nGap1 + nGap2, # nGaps
		header['begin1'], # start1
		header['end1'], # end1
		header['begin2'], # start2
		header['end2'], # end2
		line1, # seq1 with gap
		line3] # seq2 with gap
	return ir


# Return ir
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# args: (input4IS, filter)
# input4IS: (familyName, isName, seq1, seq2, minScore, minLen)
# filter: (gapopen, gapextend, match, mismatch)
# 
#def findIR4elementBySSW(input4IS, filter):
def findIR4elementBySSW(args):
	input4IS, filter = args
	familyName, isName, seq1, seq2, minScore, minLen = input4IS 
	gapopen, gapextend, match, mismatch = filter
	# set minScore based on such rule that an IR must have at least two consecutive matches 
	minScore = match * 2
	#minScore = match * 5
	#minScore = match * 10
	#minScore = match * 15
	#minScore = match * 20

	# If sequence is an empty string, then no alignment is done.
	if len(seq1) < 1 or len(seq2) < 1:
		#print('Warning: {} {} {} no alignment can be done with seq1= {} and seq2={}!'.format(familyName, isName, filter, seq1, seq2))
		return []

	# For SSW, the gap_open is defined as the total penalty when opening a gap
	gapopen = gapopen + gapextend
	'''
	ssw = ssw_wrap.Aligner(	str(seq1), 
				match = int(match), mismatch = int(mismatch), 
				gap_open = int(gapopen), gap_extend = int(gapextend), 
				report_secondary = False, report_cigar = True)
	align = ssw.align(str(seq2), min_score = float(minScore), min_len = int(minLen))
	'''

	ssw = ssw_wrap.Aligner(	seq1, 
				match = match, mismatch = mismatch, 
				gap_open = gapopen, gap_extend = gapextend, 
				report_secondary = False, report_cigar = True)
	align = ssw.align(seq2, min_score = minScore, min_len = minLen)

	if align:
		#cigarPair = tools.parseCigarString(align.cigar_string)
		alignment = tools.buildAlignment(seq1, seq2, align, align.cigar_string)
		#ir = getIRbySSWnoGap(seq1, seq2, align, cigarPair)
		ir = getIRbySSW(alignment)
		'''
		header = alignment[0]
		if header['conflict'] == True:
			line1, line2, line3 = alignment[1:]
			print('Alignment conflict: {} {} {}\n {}\n {}\n {}'.format(
				familyName, isName, filter, line1, line2, line3))
		'''
	else:
		ir = []
		#print('Warning: {} {} {} no alignment returned by SSW!'.format(familyName, isName, filter))
	
	return ir


# Find best IR for each element, and return all best IRs with one IR per element
#
# Input:
# mInput4ssw: [input4IS1, .., input4ISn]
# input4ISn: (familyName, isName, seq1, seq2, minScore, minLen)
# minScore, minLen: the minimal score and length for the alignment to be reported
# filter: (gapopen, gapextend, match, mismatch)
# 
# Return:
# mBestIR: [TIR, ...]
# TIR: [familyName, isName, ir]
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# seq1, seq2: inverted repeat sequences
def findIRbySSW(mInput4ssw, filter):
	mBestIR = []

	'''
	args2concurrent = []
	for input4IS in mInput4ssw:
		args2concurrent.append((input4IS, filter))

	if len(args2concurrent) > constants.nproc:
		nprocess = constants.nproc
	else:
		nprocess = len(args2concurrent)
	with concurrent.futures.ProcessPoolExecutor(max_workers = nprocess) as executor:
		for args, ir in zip(args2concurrent, executor.map(findIR4elementBySSW, args2concurrent)):
		# args: (input4IS, filter)
		# input4IS: (familyName, isName, seq1, seq2, minScore, minLen)
			#familyName, isName = args[0][:2]
			#mBestIR.append([familyName, isName, ir])
			mBestIR.append([args[0][0], args[0][1], ir])
	'''
	for input4IS in mInput4ssw:
		# input4IS: (familyName, isName, seq1, seq2, minScore, minLen)
		ir = findIR4elementBySSW((input4IS, filter))

		#if len(ir) > 0 and (ir[3] > 0 or ir[1]/ir[2] < constants.irSim4singleCopy):
		#	ir = []

		mBestIR.append([input4IS[0], input4IS[1], ir])

	return mBestIR


# get the sequence in which all letters are upper case
# Note: the IRs of some elements do not start from terminus of IS element, we need retrieve non-terminal IR sequence
#	from such special IS elements. In such case, the terminal sequence is often represented as lowercase letter like 
#	atcgu instead of ATCGU.
# Note: some IRs are longer than 50 which is the maximal length of element['lSeq'], then we need retrieve long sequence
#	from element['isSeq'], namely, IS_SEQ record of isfinder
def getIRstartFromEnd(endseq):
	initialIRseq = ''
	if endseq == '':
		 initialIRseq = ''
	elif endseq[0] in 'ATCGU':
		initialIRseq = endseq
	# if the initial codes are lowercase letter, it means IR do not start from terminus of IS element.
	# however, the initial codes being uppercase does not mean IR DO start from terminus of IS element, e.g. ISSod25.
	else:
		for i, c in enumerate(endseq):
			if c.isupper():
				initialIRseq = endseq[i:]
				break
	return initialIRseq
	
# Find the lowest index of endseq in isSeq, return 0 if endseq not in isSeq
# Note: assume that IR is at termini of IS if LEFT END not in IS_SEQ
def getIRstart(endSeq, isSeq):
	endseq = getIRstartFromEnd(endSeq)
	'''
	if endseq == '':
		return -1
	else:
		
		index = isSeq.find(endseq)
		if index == -1:
			index = 0
		return index
	'''
	index = isSeq.find(endseq)
	if index == -1:
		index = 0
	return index

# Comparre the predicted IR with IR records in isfinder
# return 1 if matched, 0 if not matched, else -1
# ir: [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# element['irLen']: -1 or >=0, no information about IR length is available in isfinder when -1, else information is available
# Note:	'matched' just means irId from predicted IR is equal to that from isfinder and
#	Option1: 
#		irLen from predicted IR is equal or shorter than that from isfinder.
#	Option2:
#		irLen from predicted IR is equal to that from isfinder.
#	Option3:
#		irLen from predicted IR can be lessnthan /euqal to/ greater than that from isfinder, namely, unlimited.
#
#	It is difficult to compare the exact sequences of IRs from the predicted and isfinder
#	as there is no exact IR sequence claimed in IS element record from isfinder.
#	Only elements with element['irLen'] >= 0 is used for comparison.
#
def matchIR_v2(ir, element):
	if len(ir) == 0:
		if element['irLen'][-1] == 0:
			return 1
		else:
			return 0

	num = 0

	# Option0:
	#	1) The initial bases of both lSeq and rSeq of the predicted ir are located in TIR annotated in ISfinder
	#
	#irSeq1 = (ir[-2].replace('-', '')).upper()
	#irSeq2 = (ir[-1].replace('-', '')).upper()
	# The initial 100 bp of both sequences to align are noise sequences to simulate the real world situation.
	# Sequences: lSeq = lNoise + TIR1 + lNoise, rSeq = rNoise + TIR2 + rNoise
	# Lengths: lSeqLen = rSeqLen = 100 + element['irLen'][-1] + 100
	# Both lNoise and rNoise are 100 bp long. 
	# TIR1 and TIR2 are terminal sequences of full IS element sequence, the lengths of both TIR1 and TIR2 
	# are defined by element['irLen'][-1].
	# Refer to prepare4ssw2findIRbyISfinder() for more details.
	#
	if 100 < ir[-6] <= element['irLen'][-1]+100 and 100 < ir[-4] <= element['irLen'][-1]+100:

	#irId, irLen = ir[1], ir[2]

	# Option1:
	#	1) predicted irId should be equal to that from ISfinder
	#	and
	#	2) predicted irLen should be equal or than that from ISfinder.
	#if irLen <= element['irLen'][-1] and irId == element['irLen'][0]:
	
	# Option2:
	#	1) both predicted irId and irLen should be equal to those from ISfinder.
	#if irLen == element['irLen'][-1] and irId == element['irLen'][0]:
	
	# Option3:
	#	1) irLen from predicted IR can be lessnthan /euqal to/ greater than that from isfinder, namely, unlimited.
	#if irId == element['irLen'][0]:

	# Option4:
	#	1) both irId and irId/irLen are equal to or greater than that those from isfinder
	# Refer to IS element ISEc9 which is active IS. 
	# We relax criteria because no TIR can be found for ISEc9 under Option3 condition.
	#if irId >= element['irLen'][0] and (not (irId/irLen < element['irLen'][0]/element['irLen'][-1])):

	# Option5:
	#	1) irId is equal to or greater than that one from isfinder.
	#	and
	#	2) number of matches in core regions is greater than those in other regions, 
	#		namely, irIdCore > (irId - irIdCore), which is same as irIdCore/irId > 50%.
	#	and
	#	3) irId/irLen is greater than a cutoff like 40%/45%/50%.
	#if irId >= element['irLen'][0]:
	#	irIdCore = tools.getIrIdCore(ir[-2], ir[-1])
	#	#irSeqLen = ir[-5] - ir[-6] + 1
	# cutoff: 50%, 50%
	#if irId >= element['irLen'][0] and irIdCore > (irId - irIdCore) and irId > (irSeqLen-irId):
	# cutoff: 50%, 40%
	#if irId >= element['irLen'][0] and irIdCore > (irId - irIdCore) and irId/irLen > 0.4:
	#if irId >= element['irLen'][0] and irIdCore > (irId - irIdCore) and irId/irSeqLen > 0.4:
		num = 1
	else:
		num = 0
	return num

# Comparre the predicted IR with IR records in isfinder
# return 1 if matched, 0 if not matched, else -1
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# element['irLen']: -1 or >=0, no information about IR length is available in isfinder when -1, else information is available
# Note1:	'matched' just means irLen (and irId) from the predicted IR is equal or shorter than those from the IR 
#		in ISfinder. It is difficult to compare the exact sequences of IRs from the predicted and isfinder
#		as there is no exact IR sequence claimed in IS element record from isfinder.
# Note2:	improvement, only elements with element['irLen'] >= 0 is used for comparison.
# Note3:	improvement, roughly compare the sequences from element['isSeq'] and ir when check if two IRs are matched.
def matchIR(ir, element, index):
	if len(ir) == 0:
		if element['irLen'][-1] == 0:
			return 1
		else:
			return 0
	irId, irLen, irSeq1, irSeq2 = ir[1], ir[2], ir[-2], ir[-1]
	#index = getIRstart(element['lSeq'], element['isSeq'])
	#if irLen == element['irLen'][-1] and irId == element['irLen'][0]:
	#
	if irLen <= element['irLen'][-1] and irId == element['irLen'][0]:
		if index == -1:
			#print('Warning: mismatch between LEFT END and IS_SEQ', element['isName'])
			return -1
		isIRseq = element['isSeq'][index:]
		irSeq = (irSeq1.replace('-', '')).upper()
		if  irSeq in isIRseq:
			# Note: ISSod25, IR obviously does not start from terminus of IS_SEQ though the codes in LEFT END are all uppercase.
			return 1
		else:
			return 0
	else:
		return 0

# Compare the predicted IR with IR from isfinder, and return the IR matched with ISfinder
# Return: ((nmatch, nnotmatch, ndiscard, nIS), TIR)
# Note: no sequence but only length of IR is available in ISfinder
# Note: two IRs being matched means that the predicted IR and ISfinder IR have same irIds.
#
# mfamilyFeatures: [(familyName, familyFeatures), ..., (familyName, familyFeatures)]
# familyFeatures: [isFeatures, ..., isFeatures]
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
#
# IRs, TIR: [IR, ..., IR]
# IR: [familyName, isName, ir]
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# seq1, seq2: inverted repeat sequence and its comlementary sequence, both sequence run from 5' to 3'
#
def compareIRbyISfinder_v2(IRs, mfamilyFeatures):
	TIR = []
	# nmatch: how many elements the IRs of which are matched between prediction and isfinder
	nmatch = 0
	# nnotmatch: how many elements the IRs of which are not matched between prediction and isfinder
	nnotmatch = 0
	# ndiscard: how many elements are not used for comparison because of no clear or correct IR information available in ISfinder
	ndiscard = 0
	# how many elements in total
	nIS = 0
	for IR in IRs:
		for family in mfamilyFeatures:
			if family[0] != IR[0]:
				continue
			for i, element in enumerate(family[1]):
				if element['isName'] != IR[1]:
					continue
				nIS += 1

				# No IR is found in IS element by prediction.
				if len(IR[2]) == 0:
					continue

				# No TIR information is available in ISfinder.
				if element['irLen'][-1] == -1:
					ndiscard += 1
					print('Warning: discard IS without TIR info {} {}'.format(family[0], element['isName']))
					continue

				num = matchIR_v2(IR[2], element)

				if num == 1:
					nmatch += 1
					TIR.append(IR)

				elif num == 0:
					nnotmatch += 1
				else:
					e = 'Erorr: num must be 1 or 0'
					raise RuntimeError(e)
				break
			break
	return ((nmatch, nnotmatch, ndiscard, nIS), TIR)

# Compare the predicted IR with IR from isfinder, and return the IR matched with ISfinder
# Note: no sequence but only length of IR is available in ISfinder
# Return: (nmatch, ndismatch, ndiscard, nIS, TIR)
# TIR: [IR, ..., IR]
# IR: [familyName, isName, ir]
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# seq1, seq2: inverted repeat sequence and its comlementary sequence, both sequence run from 5' to 3'
def compareIRbyISfinder(IRs, mfamilyFeatures):
	TIR = []
	# nmatch: how many elements the IRs of which are matched between prediction and isfinder
	nmatch = 0
	# nnotmatch: how many elements the IRs of which are not matched between prediction and isfinder
	nnotmatch = 0
	# ndiscard: how many elements are not used for comparison because of no clear or correct IR information available in ISfinder
	ndiscard = 0
	# how many elements in total
	nIS = 0
	for IR in IRs:
		for family in mfamilyFeatures:
			if family[0] != IR[0]:
				continue
			for i, element in enumerate(family[1]):
				if element['isName'] != IR[1]:
					continue
				nIS += 1
				if len(IR[2]) > 0:
					irSeq1, irSeq2 = IR[2][-2], IR[2][-1]
				else:
					irSeq1 = irSeq2 = ''

				# left location of expected IR in isfinder
				index = getIRstart(element['lSeq'], element['isSeq'])

				# left location of IR found by prediction
				irSeq = (irSeq1.replace('-', '')).upper()
				irIndex = element['isSeq'].find(irSeq)

				if element['irLen'][-1] == -1:
					ndiscard += 1
					print('Warning: discard IS without TIR info (isfinder {} {}) {} {} ({} {}) {}'.format(element['irLen'], index+1,
						family[0], element['isName'], len(irSeq1), irIndex+1, irSeq1.upper()))
					continue
				num = matchIR(IR[2], element, index)
				if num == -1:
					ndiscard += 1
					print('Warning: discard IS with LEFT END not in IS_SEQ (isfinder {} {}) {} {} ({} {}) {}'.format(element['irLen'], index+1,
						family[0], element['isName'], len(irSeq1), irIndex+1, irSeq1.upper()))
					continue
				elif num == 1:
					nmatch += 1
					TIR.append(IR)
					'''
					if len(irSeq1) != element['irLen'][-1]:
						print('Lucky: matched with different lengths of IR (isfinder {} {}) {} {} ({} {})\n{}\n{}'.format(element['irLen'], index+1,
								family[0], element['isName'], len(irSeq1), irIndex+1, irSeq1.upper(), irSeq2.upper()))
					else:
						print('Lucky: matched (isfinder {} {}) {} {} ({} {})\n{}\n{}'.format(element['irLen'], index+1,
								family[0], element['isName'], len(irSeq1), irIndex+1, irSeq1.upper(), irSeq2.upper()))
					'''
				else:
					nnotmatch += 1
					'''
					print('Warning: not matched (isfinder {} {}) {} {} ({} {})\n{}\n{}'.format(element['irLen'], index+1,
						family[0], element['isName'], len(irSeq1), irIndex+1, irSeq1.upper(), irSeq2.upper()))
					'''
	return ((nmatch, nnotmatch, ndiscard, nIS), TIR)


# return True if element is better than best_element, else return False (equal or worse)
def compare_element(element, best_element):
	nseq1, nseq2 = 0, 0
	npep1, npep2 = 0, 0
	for item in element[6:-2]:
		if 'IS_SEQ' in item[0][0]:
			nseq1 += 1
			continue
		if 'IS_PEP' in item[0][0]:
			npep1 += 1
	for item in best_element[6:-2]:
		if 'IS_SEQ' in item[0][0]:
			nseq2 += 1
			continue
		if 'IS_PEP' in item[0][0]:
			npep2 += 1
	#print(element[0][0][0], npep1, nseq1, best_element[0][0][0], npep2, nseq2)
	if npep1 > npep2:
		return True
	elif npep1 == npep2 and nseq1 > nseq2:
		return True
	else:
		return False

# select one element with the most complete sequence information for each species (origin)
def best_element_per_origin(family):
	family_refined = []
	# origin_element = {origin_name: full_element_data}
	origin_element = {}
	for element in family:
		origin_name = element[1][1][0]
		if origin_name in origin_element:
			if compare_element(element, origin_element[origin_name]):
				origin_element[origin_name] = element
		else:
			origin_element[origin_name] = element
	for best in origin_element:
		family_refined.append(origin_element[best])

	return family_refined

# select one element with the most complete sequence information for each Group
def best_element_per_group(family):
	family_refined = []
	group_element = {}
	for element in family:
		group_name = element[0][2][3].strip()
		if group_name == '' or group_name == '-':
			family_refined.append(element)
		elif group_name in group_element:
			if compare_element(element, group_element[group_name]):
				group_element[group_name] = element
		else:
			group_element[group_name] = element
	
	for best in group_element:
		family_refined.append(group_element[best])

	return family_refined


def translate_genome_dna_v1(dna, output_path, seq_type, train_model):
#./run_FragGeneScan.pl -genome=./example/NC_000913.fna -out=./example/NC_000913.test  -complete=1  -train=complete
	file_name_index = dna.rfind('/')
	output_file = output_path + dna[file_name_index:]
	gene_translate_cmd = "/u/zhiqxie/informatics/inst/FragGeneScan1.19/run_FragGeneScan.pl"
	input = "-genome=" + dna
	output = "-out=" + output_file
	seq_type = "-complete=" + seq_type
	train_model = "-train=" + train_model
	#nthread = '-thread=8'
	nthread = '-thread=2'
	cmd_line = '{0} {1} {2} {3} {4} {5}'.format(gene_translate_cmd, input, output, seq_type, train_model, nthread)
	do_FragGeneScan = shlex.split(cmd_line)
	p = subprocess.Popen(do_FragGeneScan, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	try:
		outs, errs = p.communicate(timeout = 1500)
	except subprocess.TimeoutExpired:
		p.kill()
		outs, errs = p.communicate()
	return (outs, errs)

def translate_genome_dna_v2(dna, output_file, seq_type, train_model):
	#./run_FragGeneScan.pl -genome=./example/NC_000913.fna -out=./example/NC_000913.test  -complete=1  -train=complete

	gene_translate_cmd = "/u/zhiqxie/informatics/inst/FragGeneScan1.19/run_FragGeneScan.pl"
	input = "-genome=" + dna
	output = "-out=" + output_file
	seq_type = "-complete=" + seq_type
	train_model = "-train=" + train_model
	#nthread = '-thread=8'
	nthread = '-thread=1'
	cmd_line = '{0} {1} {2} {3} {4} {5}'.format(gene_translate_cmd, input, output, seq_type, train_model, nthread)
	do_FragGeneScan = shlex.split(cmd_line)

	#return subprocess.check_call(do_FragGeneScan, shell=False, universal_newlines=False)
	return subprocess.call(do_FragGeneScan, shell=False, universal_newlines=False)
	#return subprocess.check_output(do_FragGeneScan, stderr=subprocess.STDOUT, shell=False, universal_newlines=False)
	'''
	p = subprocess.Popen(do_FragGeneScan, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	try:
		outs, errs = p.communicate(timeout = 1500)
	except subprocess.TimeoutExpired:
		p.kill()
		outs, errs = p.communicate()
	return (outs, errs)
	'''

def translate_genome_dna_v3(args):
	dna, output_file, seq_type, train_model = args
	#./run_FragGeneScan.pl -genome=./example/NC_000913.fna -out=./example/NC_000913.test  -complete=1  -train=complete -thread==4

	gene_translate_cmd = constants.FragGeneScan
	input = "-genome=" + dna
	output = "-out=" + output_file
	seq_type = "-complete=" + seq_type
	train_model = "-train=" + train_model
	#nthread = '-thread=16'
	#nthread = '-thread=8'
	#nthread = '-thread=4'
	nthread = '-thread='+str(constants.nthread)
	cmd_line = '{0} {1} {2} {3} {4} {5}'.format(gene_translate_cmd, input, output, seq_type, train_model, nthread)
	#cmd_line = '{0} {1} {2} {3} {4}'.format(gene_translate_cmd, input, output, seq_type, train_model)
	do_FragGeneScan = shlex.split(cmd_line)

	#return subprocess.call(do_FragGeneScan, shell=False, universal_newlines=False)
	return subprocess.call(cmd_line, shell=True, universal_newlines=False)

def is_hmmsearch(hmm, database, output):
	hmmsearch_cmd = "/u/zhiqxie/informatics/inst/hmmer-3.1b2/bin/hmmsearch"
	options = "--tblout " + output + " " + "--max"
	cmd_line = '{0} {1} {2} {3}'.format(hmmsearch_cmd, options, hmm, database)
	do_hmmsearch = shlex.split(cmd_line)

	#return subprocess.check_call(do_hmmsearch, shell=False, universal_newlines=False)
	return subprocess.call(do_hmmsearch, shell=False, universal_newlines=False)
	#return subprocess.check_output(do_hmmsearch, stderr=subprocess.STDOUT, shell=False, universal_newlines=False)
	'''
	p = subprocess.Popen(do_hmmsearch, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	try:
		outs, errs = p.communicate(timeout = 150)
	except subprocess.TimeoutExpired:
		p.kill()
		outs, errs = p.communicate()
	return (outs, errs)
	'''

def is_hmmsearch_v2(args):
	hmm, database, output = args
	hmmsearch_cmd = constants.hmmsearch
	#dir, output = os.path.split(output)
	#nthread = constants.nthread
	nthread = constants.nproc
	options = ' '.join(["--tblout", output, "--max --noali", "--cpu", str(constants.nthread)])
	cmd_line = ' '.join([hmmsearch_cmd, options, hmm, database])
	do_hmmsearch = shlex.split(cmd_line)

	#return subprocess.call(cmd_line, shell=True, universal_newlines=False, stdout=subprocess.DEVNULL)
	return subprocess.call(do_hmmsearch, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)
	#return subprocess.check_call(do_hmmsearch, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)

# run phmmer as:
# phmmer --tblout phmmerHitsFile --max --noali --cpu nthread seqFile databaseFile
# seqFile: profile HMM models file, created by hmmbuild in HMMer package
# databaseFile: proteome file in which multiple protein amino acid sequence are placed in FASTA format
def is_phmmer(args):
	seqFile, database, output = args
	phmmer_cmd = constants.phmmer
	#dir, output = os.path.split(output)
	#nthread = constants.nthread
	nthread = constants.nproc
	options = ' '.join(["--tblout", output, "--max --noali", "--cpu", str(constants.nthread)])
	cmd_line = ' '.join([phmmer_cmd, options, seqFile, database])
	do_search = shlex.split(cmd_line)

	#return subprocess.call(cmd_line, shell=True, universal_newlines=False, stdout=subprocess.DEVNULL)
	return subprocess.call(do_search, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)
	#return subprocess.check_call(do_search, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)


if __name__ == "__main__":
	pass
