#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import time, random
import os
import sys
import argparse
import datetime
import operator
import concurrent.futures
import subprocess, shlex

import is_source_data_parse, isfinder_page 
import tools
import constants

isfinderGenomePath = '/u/zhiqxie/xie/is/isfinder_genome'
isfinderAnnotationPath = isfinderGenomePath

rootpath = constants.rootpath

PRINT_ISFINDER_FEATURES = True
#PRINT_ISFINDER_FEATURES = False

# shortest IR registered in ISfinder database
shortestIR = 2

# cd-hit command to cluster protein sequences at 1.0 identity
cdhit = '/u/zhiqxie/informatics/inst/cd-hit-v4.6.4-2015-0603/cd-hit'

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
psicdhit = '/u/zhiqxie/informatics/inst/cd-hit-v4.6.4-2015-0603/psi-cd-hit/psi-cd-hit.pl'
clstr_rev = '/u/zhiqxie/informatics/inst/cd-hit-v4.6.4-2015-0603/clstr_rev.pl'

# switch to nr100, nr90, nr60 hierachical clustering
#lowIdent = False
# switch to nr100, nr90, nr60, nr30 hierachical clustering
lowIdent = True

# IS elements in ISfinder database, which have multiple IS_PEP records and the longest IS_PEP is not
# the Transposase.
# Note: family IS200/IS605 has 3 groups, IS200 which usually has one short tpase ORF, 
#	IS1341 which usually has one long accessory gene ORF, 
#	IS606 which usually has one short tpase ORF and one long accessory gene ORF.
#	The next version of profile HMM should pick the shortest tpase ORF (or both tpase and accessory gene ORFs)
#	for each member	of IS606 group when we choose tpase peptide to build profile hidden markov models.
shortTpases = {
		'ISKPN25': 4, # ISL3. The tpase is IS_PEP4 among four IS_PEP records.
		'ISSM4': 5, # ISL3. The tpase is IS_PEP5 among five IS_PEP records.
		'IS5564': 2, # IS481. The tpase is IS_PEP2 among two IS_PEP records.
		'ISCGL1': 2, # IS481. The tpase is IS_PEP2 among two IS_PEP records.
		#'IS606': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		#'IS1253A': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		#'IS1253B': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		#'IS605': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		#'ISHP608': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		#'ISCLTE2': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		#'ISTAC1': 1, # IS200/IS605. The tpase is IS_PEP among two IS_PEP records.
		}

is_family_list_file_url = 'https://www-is.biotoul.fr/is/IS_infos/is_family.html'
is_family_list_file = 'is_family.html'
prokaryote_family_list = [
			'IS1',
			'IS1595',
			'IS3',
			'IS481',
			'IS4',
			'IS701',
			'ISH3',
			'IS1634',
			'IS5',
			'IS1182',
			'IS6',
			'IS21',
			'IS30',
			'IS66',
			'IS91',
			'IS110',
			'IS200/IS605',
			'IS607',
			'IS256',
			'IS630',
			'IS982',
			'IS1380',
			'ISAs1',
			'ISL3',
			'ISAzo13',
			'ISKra4',
			'ISNCY',
			#'ISH6', # new family found by xie on 1/4/2016
			#'Tn3', # composite transposon
			#'Tn554',
			]


def check_synonyms_iso_family_group(element0):

	if(element0[1] == ['Synonyms', 'Iso', 'Family', 'Group']):
		return False
	return True

def check_origin_hosts(element1):
	if(element1[0] == ['Origin', 'Hosts']):
		return False
	return True

def check_accessionNumber_len(element2):
	if(element2[0] == ['Accession Number', 'Length', 'IR', 'DR']):
		return False
	return True

def check_transActivity_orfLen(element3):
	if((element3[0][0][:13] == 'Transposition') and (element3[1][0][:3] == 'ORF')):
		return False
	return True

def check_ends(element4):
	if(element4[0][0][:8] == 'Left End' and element4[1][0][:9] == 'Right End'):
		return False
	return True

def check_insertionSites(element5):
	if(element5[0][0][:15] == 'Insertion Sites'):
		return False
	return True


def check_isSeq(element6):
	if(element6[0][0][:6] == 'IS_SEQ'):
		return False
	return True

def check_isPep(element7, element_3):
	#print(element7[0][0][:6], element_3[0][0][:6])
	if(element7[0][0][:6] == element_3[0][0][:6]):
		return False
	return True

def check_comments(element_2):
	if(element_2[0][0][:8] == 'Comments'):
		return False
	return True

def check_references(element_1):
	if(element_1[0][0][:10] == 'References'):
		return False
	return True

def check_is_element_data_integrity(element):
	n_data_item = len(element)
	if n_data_item < 10:
		print("Error! Please check IS data items missing:", n_data_item, "< 10 data items:")
		return True
	if check_synonyms_iso_family_group(element[0]):
		print("Error! Please check IS synonyms_iso_family_group integrity.")
		return True
	if check_origin_hosts(element[1]):
		print("Error! Please check IS origin_hosts.")
		return True
	if check_accessionNumber_len(element[2]):
		print("Error! Please check IS 'Accession Number', 'Length', 'IR', 'DR'")
		return True
	if check_transActivity_orfLen(element[3]):
		print("Error! Please check Transposition and ORF")
		return True
	if check_ends(element[4]):
		print("Error! Please check IS Left End and Right End")
		return True
	if check_insertionSites(element[5]):
		print("Error! Please check IS Insertion Sites")
		return True
	if check_isSeq(element[6]):
		print("Error! Please check IS IS_SEQ")
		return True
	if check_isPep(element[7], element[-3]):
		print("Error! Please check IS IS_PEP* and Comments")
		return True
	if check_comments(element[-2]):
		print("Error! Please check IS Comments")
		return True
	if check_references(element[-1]):
		print("Error! Please check IS References")
		return True

	return False


def is_seq_missing(element):
	for item in element:
		if 'IS_SEQ' in item[0][0]:
			return False
	return True

def is_pep_missing(element):
	for item in element:
		if 'IS_PEP' in item[0][0]:
			return False
	return True

# Check if there is any type of data item missed in IS element raw data.
def find_element_data_item_missing(element):
	data_missing_flag = False

	if not ('Synonyms' in element[0][1] and 'Iso' in element[0][1] 
	and 'Family' in element[0][1] and 'Group' in element[0][1]):
		print('Warning: ', lement[0][2][2], element[0][0][0], "No 'Synonyms', 'Iso', 'Family', 'Group' records found")
		data_missing_flag = True
	if not ('Origin' in element[1][0] and 'Hosts' in element[1][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'Origin', 'Hosts' records found")
		data_missing_flag = True
	if not ('Accession Number' in element[2][0] and 'Length' in element[2][0]
	and 'IR' in element[2][0] and 'DR' in element[2][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'Accession Number', 'Length', 'IR', 'DR' records found")
		data_missing_flag = True
	if not ('Transposition' in element[3][0][0] and 'ORF' in element[3][1][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'Transposition', 'ORF' records found")
		data_missing_flag = True
	if not ('Left End' in element[4][0][0] and 'Right End' in element[4][1][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'Left End', 'Right End' records found")
		data_missing_flag = True
	if not ('Insertion Sites' in element[5][0][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'Insertion Sites' records found")
		data_missing_flag = True
	if not ('IS_SEQ' in element[6][0][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'IS_SEQ' records found")
		data_missing_flag = True

	"""
	if not ('IS_PEP' in element[7][0][0]):
		print(element[0][2][2], element[0][0][0], "No 'IS_PEP' records found")
		data_missing_flag = True

	npep = element[-3][0][0][6:].rstrip(': ') 
	if not ('IS_PEP' in element[-3][0][0][:6] and npep.isdecimal()):
		#print(element[0][2][2], element[0][0][0], "No 'IS_PEPn' records found")
		data_missing_flag = True
	elif(int(npep) > 3):
		print(element[0][2][2], element[0][0][0], npep, "IS_PEP records found")
	"""
	if is_pep_missing(element[6:]):
		print('Warning:', element[0][2][2], element[0][0][0], "No 'IS_PEP' records found")
		data_missing_flag = True

	if not ('Comments' in element[-2][0][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'Comments' records found")
		data_missing_flag = True
	if not ('References' in element[-1][0][0]):
		print('Warning: ', element[0][2][2], element[0][0][0], "No 'References' records found")
		data_missing_flag = True

	return data_missing_flag


def is_family_missing(element):
	return False
def is_origin_missing(element):
	return False
	pass
def is_accession_number_missing(element):
	return False
	pass
def is_orf_missing(element):
	return False
	pass
def is_end_missing(element):
	return False
	pass
def is_insertion_sites_missing(element):
	return False
	pass

# Clean IS_SEQ record
# Return cleaned IS_SEQ record
def clean_seq(is_name, is_seq):
	cleaned_seq = ''
	for nuc in is_seq[1][0]:
		if nuc.isalpha():
			cleaned_seq += nuc
		elif not nuc.isspace():
			print('Warning: {} remove strange base code ({})'.format(is_name, nuc))

	is_seq[1][0] = cleaned_seq
	
	return is_seq

# Clean IS_PEP* records
# Return cleaned IS_PEP* record
def clean_pep(is_name, is_pep):
	cleaned_pep = '' 
	# remove any character in ' *()'
	#map = str.maketrans('', '', ' *()')
	# remove any character in '()'
	#map = str.maketrans('', '', '()')
	#is_pep[1][0] = is_pep[1][0].translate(map) 
	for a in is_pep[1][0]:
		if a.isalpha():
			cleaned_pep += a
		elif not a.isspace():
			print('Warning: {} remove strange amino acid code ({})'.format(is_name, a))

	is_pep[1][0] = cleaned_pep
	return is_pep

# Function: clean known problem for each data type, e.g., removing unknown codes (' ', '*', ')' etc.) 
# in nucleic acid and peptide sequence.
# Return: cleaned element
#         [[], [], [], [], [], [], [['IS_SEQ'], ['...']], [[IS_PEP], ['...']], [[IS_PEPn], ['...']], [['Comments:'], ['...']], [['References:'], ['...']]]
def clean_element(element):
	cleaned_element = element
	# The normal IS element will have 10 or more data records in each IS element
	'''
	if is_family_missing(element):
		cleaned_element.append([])
	else:
		cleaned_element.append(element[0])
	if is_origin_missing(element):
		cleaned_element.append([])
	else:
		cleaned_element.append(element[1])
	if is_accession_number_missing(element):
		cleaned_element.append([])
	else:
		cleaned_element.append(element[2])
	if is_orf_missing(element):
		cleaned_element.append([])
	else:
		cleaned_element.append(element[3])
	if is_end_missing(element):
		cleaned_element.append([])
	else:
		cleaned_element.append(element[4])
	if is_insertion_sites_missing(element):
		cleaned_element.append([])
	else:
		cleaned_element.append(element[5])
	'''
	if not is_seq_missing(element[6:-2]):
		index_pep = 7
		cleaned_element[6] = clean_seq(element[0][0][0], element[6])
		if cleaned_element[6][1][0] == '':
			print('Warning:', element[0][2][2], element[0][0][0], 'empty IS_SEQ record')
	else:
		index_pep = 6
	if not is_pep_missing(element[index_pep:-2]):
		# one or more peptide sequences are available
		for i, is_pep in enumerate(element[index_pep:-2]):
			cleaned_element[i + index_pep] = clean_pep(element[0][0][0], is_pep)
		if cleaned_element[i + index_pep][1][0] == '':
			print('Warning:', element[0][2][2], element[0][0][0], 'empty IS_SEQ record')

	return cleaned_element


def is_group_stat(family):
	group = {}
	for element in family:
		ISname = element[0][0][0]
		key = element[0][2][3]
		if key in group:
			#group[key] += 1
			group[key].append(ISname)
		else:
			#group[key] = 1
			group[key] = [ISname]
	return group

def is_origin_stat(family):
	pass

def is_seq_stat(family):
	seq = {}
	#seq['1'] = 0
	seq['1'] = []
	seq['len'] = []
	for element in family:
		for item in element[6:-2]:
			if 'IS_SEQ' in item[0][0]:
				#seq['1'] += 1
				seq['1'].append(element[0][0][0])
				seq['len'].append((element[0][0][0], len(item[1][0])))
				break
		else:
			seq['len'].append((element[0][0][0], 0))
	if len(family) - len(seq['1']) > 0:
		#seq['0'] = len(family) - seq['1']
		seq['0'] = [element[0][0][0] for element in family if element[0][0][0] not in seq['1']]
	seq['len'] = sorted(seq['len'], key = operator.itemgetter(1))
	return seq


# count how many IS_PEP records in an IS element and register the longest one as the transposase.
def is_pep_stat(family):
	pep = {}
	pep['len'] = []
	for element in family:
		ISname = element[0][0][0]
		npep = 0
		pepLen = []
		for item in element[6:-2]:
			if 'IS_PEP' in item[0][0]:
				npep += 1
				#pep['len'].append((ISname, len(item[1][0])))
				# register the lengths of all PEPs
				pepLen.append(len(item[1][0]))
		if npep == 0:
			pep['len'].append((ISname, 0))
		else:
			pep['len'].append((ISname, max(pepLen)))
			#pep['len'].append((ISname, min(pepLen)))
		key = str(npep)
		if key in pep:
			#pep[key] += 1
			pep[key].append(ISname)
		else:
			#pep[key] = 1
			pep[key] = [ISname]
		if npep == 0 or npep > 1:
			print(element[0][2][2], ISname, npep, "IS_PEP records")
		# register the longest PEP
	pep['len'] = sorted(pep['len'], key = operator.itemgetter(1))
	return pep


def is_seq_pep_stat(family):
	seq_pep = {'NO_SEQ_PEP': 0, 'NO_SEQ': 0, 'NO_PEP': 0, 'SEQ_PEP': 0}
	for element in family:
		seq_flag = False
		pep_flag = False
		for item in element[6:-2]:
			if 'IS_SEQ' in item[0][0]:
				seq_flag = True
				continue
			if 'IS_PEP' in item[0][0]:
				pep_flag = True
				break
		if seq_flag == False and pep_flag == False:
			seq_pep['NO_SEQ_PEP'] += 1
		elif seq_flag == False and pep_flag == True:
			seq_pep['NO_SEQ'] += 1
		elif seq_flag == True and pep_flag == False:
			seq_pep['NO_PEP'] += 1
	seq_pep['SEQ_PEP'] = len(family) - seq_pep['NO_SEQ_PEP'] - seq_pep['NO_SEQ'] - seq_pep['NO_PEP']

	return seq_pep

def isfinder_statistics(families):
	n_IS = 0
	n_family = 0
	stat = []
	for family in families:
		group = is_group_stat(family)
		origin = is_origin_stat(family)
		seq = is_seq_stat(family)
		pep = is_pep_stat(family)
		seq_pep = is_seq_pep_stat(family) 

		family_name = family[0][0][2][2]

		#nIS = sum(group.values())
		nIS = sum([len(groupIS) for groupIS in group.values()])
		newstat = ({'family_name': family_name,
				family_name: nIS, 
				'Group': group, 
				'IS_SEQ': seq, 'IS_PEP': pep, 'IS_SEQ_PEP': seq_pep})
		stat.append(newstat)

	return stat


# Remove imperfect IS elements from isfinder 
# - remove IS which has no IS_SEQ or IS_PEP records, namely, keep only IS with both IS_SEQ and IS_PEP records.
# - remove IS which has no ORF with the complete boundary coordinates on IS sequence.
def refineIsfinder(families):
	mfamily = []
	for family in families:
		familyName = family[0][0][2][2].strip()
		newFamily = []
		for element in family:
			isName = element[0][0][0].strip()
			if element[6][0][0][:6] != 'IS_SEQ':
				print('Warning: remove IS element without IS_SEQ record', element[6], familyName, isName)
				continue
			if element[7][0][0][:6] != 'IS_PEP':
				print('Warning: remove IS element without IS_PEP records', element[7], familyName, isName)
				continue
			orfsStr = element[3][1][1].strip()
			if tools.hasNumbers(orfsStr) != True or tools.hasBrackets(orfsStr) != True or len(getOrfs(orfsStr)) == 0:
				print('Warning: reomove IS element withouth known ORF coordinates info', element[3][1], familyName, isName)
				continue
			newFamily.append(element)
		mfamily.append(newFamily)

	return mfamily


# Output statistics of isfinder database
def output_isfinder_stat(isfinder_stat, fp):
	m_family_minmax_len = []
	#tb1_title = '# member distribution per group in family'
	#tb1_fmt = '{:11} {:6}'
	tb2_title = '# IS length distribution in family'
	tb2_fmt = '{:11} {:6}'
	for family_stat in isfinder_stat:
		family_name = family_stat['family_name']
		seq = family_stat['IS_SEQ']
		pep = family_stat['IS_PEP']
		for element in seq['len']:
			if element[1] == 0:
				continue
			else:
				break
		seq_shortest_name, seq_shortest_len = element[0], element[1]
		for element in pep['len']:
			if element[1] == 0:
				continue
			else:
				break
		pep_shortest_name, pep_shortest_len = element[0], element[1]
		family_minmax_len = (family_name, 
						(seq_shortest_name, seq_shortest_len, seq['len'][-1][0], seq['len'][-1][1]), 
						(pep_shortest_name, pep_shortest_len, pep['len'][-1][0], pep['len'][-1][1]))
		m_family_minmax_len.append(family_minmax_len)

		print(tb2_title, family_name, file = fp)
		for element in family_stat['IS_SEQ']['len']:
			print(tb2_fmt.format(element[0], element[1]), file = fp)

	# Features for each family
	print('\n# The features of each IS family', file = fp)
	featuresFmtStr = '{:11} {:>7} {:>6} {:>5} {:>5} {:>8} {:>6} {:11} {:>11} {:11} {:>10} {:11} {:>11} {:11} {:>10}'
	print(featuresFmtStr.format(
		'Family', 'Members', 'Groups', 
		'NoSeq', 'NoPep', 'NoSeqPep', 'pepSeq', 
		'IS', 'ShortestSeq', 'IS', 'LongestSeq', 
		'IS', 'ShortestPep', 'IS', 'LongestPep'), 
		file = fp)
	for family, family1 in zip(isfinder_stat, m_family_minmax_len):
		familyName = family['family_name']
		if familyName != family1[0]:
			e = 'Family name error: {} {}'.format(familyName, family[0])
			raise RuntimeError(e)
		noSeqPep = family['IS_SEQ_PEP']
		print(featuresFmtStr.format(
			familyName, family[familyName], len(family['Group']), 
			noSeqPep['NO_SEQ'], noSeqPep['NO_PEP'], noSeqPep['NO_SEQ_PEP'], noSeqPep['SEQ_PEP'], 
			family1[1][0], family1[1][1], family1[1][2], family1[1][3], 
			family1[2][0], family1[2][1], family1[2][2], family1[2][3]), 
			file = fp)


# Return isNames with irLen >= shortestIR 
def getIsNamesTIR(mfamilyFeatures):
	isNames = set()
	for family in mfamilyFeatures:
		for element in family[1]:
			if element['irLen'][-1] < shortestIR:
				continue
			isNames.add(element['isName'])
	return isNames

# find out in which ISs TIRs art not found by the specific TIR finding alignment algorithm
def notMatched(isNamesTIR, mTIRbyAlign):
	return isNamesTIR - mTIRbyAlign.keys()


# Convert the data structure of mfamilyFeatures from list to dictionary
# Input
# mfamilyFeatures: [(familyName, familyFeatures), ..., (familyName, familyFeatures)]
# familyFeatures: [isFeatures, ..., isFeatures]
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
# Output
# mfamilyFeaturesDic: {'familyName': familyDic, ..., 'familyName': familyDic}
# familyDic: {'isName': isFeatures, ..., 'isName': isFeatures}
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
def convertFromList2Dic(mfamilyFeatures):
	mfamilyFeaturesDic = {}
	for family in mfamilyFeatures:
		familyDic = {}
		for isFeatures in family[1]:
			familyDic[isFeatures['isName']] = isFeatures
		mfamilyFeaturesDic[family[0]] = familyDic
	return mfamilyFeaturesDic

# Output TIRs for all IS elements with TIR information available in isfinder database
# mTIR: {isName: (familyName, isName, tirs), ..., isName: (familyName, isName, tirs)}
# tirs: [tir, ..., tir]
# tir: (irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2), ir without score attached
def addTIR2features(mfamilyFeaturesDic, mTIR):
	for item in mTIR.values():
		familyName, isName = item[:2]

		# mfamilyFeaturesDic[familyName][isName], isFeatures, which is a dictionary.
		#
		# The 'tir' feature is created for the IS element with irLen >= shortestIR in ISfinder
		#  while no 'tir' feature is created for IS element with isLen < shortestIR in ISfinder.
		# Option1:
		# 	1) We simply choose the first tir as the TIR of IS element if multiple tirs found 
		#		by alignment algorithm.
		#if len(item[3]) == 0:
		#	mfamilyFeaturesDic[familyName][isName]['tir'] = ()
		#else:
		#	mfamilyFeaturesDic[familyName][isName]['tir'] = item[2][0]
		# Option2:
		#	1) We keep all best tirs for TIR of IS element.
		mfamilyFeaturesDic[familyName][isName]['tir'] = item[2]

#def outputTIR(mfamilyFeaturesDic, mTIR):
def outputTIR(mTIR):
	fp = open('ISfinderTIR.list', 'w')
	print('TIR sequence for each IS element with irLen >= ', shortestIR, 'in ISfinder', file=fp)
	print('Note: TIRs of different IS elements might be identified by alignment algorithm under different filters', file=fp)
	print('-' * 50, file=fp)
	for item in sorted(mTIR.values(), key = operator.itemgetter(0, 1)):
		# Option1: output only the first tir among all tirs of each IS element
		print('{:<15} {:<13} {}\n {:28} {}\n {:28} {}\n {:28} {}'.format(
			item[0], item[1], item[2][0][:-2], 
			' ', item[2][0][-2], 
			' ', tools.buildMatchLine(item[2][0][-2], item[2][0][-1]), 
			' ', item[2][0][-1]), file=fp)
		'''
		#
		# Option2: output all tirs of each IS element
		familyName, isName, tirs = item
		line1 = '{:<15} {:<13}'.format(familyName, isName)
		ntir = 0
		for tir in tirs:
			#line2 = ' {:28} {}'.format(' ', tir[:7])

			if ntir < 1:
				line2 = ' {:<3} {:<3} {:<3} {:>7} {:>7} {:>7} {:>7}'.format(
						tir[0], tir[1], tir[2], tir[3], tir[4], tir[5], tir[6])
			else:
				line2 = '\n {:28} {:<3} {:<3} {:<3} {:>7} {:>7} {:>7} {:>7}'.format(
						' ', tir[0], tir[1], tir[2], tir[3], tir[4], tir[5], tir[6])

			line3 = '\n {:28} {}'.format(' ', tir[7])
			line4 = '\n {:28} {}'.format(' ', tools.buildMatchLine(tir[7], tir[8]))
			line5 = '\n {:28} {}'.format(' ', tir[8])
			#print(line2 + line3 + line4 + line5)
			line1 += line2 + line3 + line4 + line5
			ntir += 1
		print(line1, file=fp)
		if ntir > 1:
			print('hello, multiple TIRs found', file=fp)
		'''
	'''
	for familyName, familyDic in mfamilyFeaturesDic.items():
		for isName, features in familyDic.items():
			if 'tir' not in features.keys():
				continue
			print('{:<15} {:<13} {}'.format(familyName, isName, features['tir']))
	'''

def build_msa(seqfilename, msafilename):
	clustalo_cmd = "/u/zhiqxie/informatics/inst/clustal-omega-1.2.1/bin/clustalo"
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
	clustalo_cmd = "/u/zhiqxie/informatics/inst/clustal-omega-1.2.1/bin/clustalo"
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
	hmmbuild_cmd = "/u/zhiqxie/informatics/inst/hmmer-3.1b2/bin/hmmbuild"
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


# update list of family names based on the isfinder annotations refered by path
# familyNames: list, a list of family names
# path: a chracter string, path refered to the directory where isfinder genome annotations are available
# return names: set, containning all family names appearing in isfinder genome annotations
def updateFamilyList(familyNames, path):
	filePattern = os.path.join(path, '*.csv')
	names = set(familyNames)
	import glob
	for file in glob.iglob(filePattern):
		seqAnnotations = tools.isfinder_IS_in_genome(file)
		names.update([item[1] for item in seqAnnotations])

	return names

# Return non-Insertion Sequence family set without IS family
# Note: we assume family names begining with 'IS' are insertion sequence families.
# familyNames: {'name1', ..., 'namen'}
# names: {'name1', ..., 'namen'}
def filterOutIS(familyNames):
	names = set()
	for family in familyNames:
		if family[:2].upper() == 'IS':
			continue
		names.add(family)
	return names

# search digits from the right side of character sequence
def rGetNum(s):
	startR = ''
	# look for the digits in s till non-digital character is met
	for ch in reversed(s):
		if ch.isdigit():
			startR += ch
		elif startR == '':
			continue
		else:
			break
	start = startR[::-1]
	if start.isdigit():
		return int(start)
	else:
		return None

# search for digits from the left side of character sequence
def lGetNum(s):
	end = ''
	# look for the starting digits in s till non-digital character is met
	for ch in s:
		if ch.isdigit():
			end += ch
		elif end == '':
			continue
		else:
			break
	if end.isdigit():
		return int(end)
	else:
		return None

# Find ORF related character strings like 'ORFn1', 'ORFn2'.
# Return the probable number of ORFs in IS, which is max(n1, ..., ni)
def getOrfByComments(comments):
	orfs = [int(chrs[3:]) for chrs in re.findall('orf\d+', comments.lower())]
	if len(orfs) > 0:
		return max(orfs)
	else:
		return 0


# Convert ORFs record into orfs with digital start and end.  # Format of ORF record examples: 
# s = '91 (56-331) 174(229-753) 232(56-753)[F]', 
# orfs = [[56-331], [229-753], [56-753]]
# ISRle39d: s = '250(83-834', 
# '413(30->958)'
def getOrfs(s):
	orfs = []
	segments = s.split('-')
	for s1, s2 in zip(segments, segments[1:]):
		start, end = rGetNum(s1.rstrip()), lGetNum(s2.lstrip())
		if start is None or end is None:
			continue
		orfs.append([start, end])
	return orfs

# convert ORFs record into orfs with digital start and end
def getOrfs1(s):
	if ')' in s:
		# s = '91 (56-331) 174(229-753) 232(56-753)[F]'
		# orfStr = ['56-331', '229-753', '56-753']
		orfsStr = [item[:item.find(')')] for item in s.split('(') if '-' in item]
	else:
		# s = '250(83-834', IS6, ISRle39d.
		# orfsStr = ['83-834']
		orfsStr = [item for item in s.split('(') if '-' in item]
		

	# get a flat list
	# orfs = [56, 331, 229, 753, 56, 753]
	#orfs = [int(num) for item in orfsStr for num in item.split('-') if num.isdigit()]

	# get a list of list
	# orfs = [[56, 331], [229, 753], [56, 753]]
	orfs = [[int(num) for num in item.split('-') if num.isdigit()] for item in orfsStr]
	#print('orfsStr', orfsStr)
	
	return orfs

# convert Insertion Sites record into insertSites with all possible sites
# Insertion Sites example:
# s = 'CAGCTGCTGGTGCAGCAGGC ( ) CGTTACCAGTTCTTCGATGT'
# s = 'AGACGGCACT ( CGATAATT ) CCAACGCTGT'
# s = 'AGT ( ) NNNNNNNNNN'
# s = 'ATTTACTTAT ( GACATTAAA ) AGTAACTTTT CACCGGCATA ( CTCTGCGAC ) ATCGTATAAC ATCAGATGAA ( TAGCCTGTC ) TTACAGTAAA ATCTAATAAC ( TGTTTCTTT ) TTGTTTTTAT NNNNNNNNTT ( GTGTTTTTC ) ATNNNNNNNN'
# ISSod25: s = 'GCACGAATAAAATGTT GAAC ( ) CAAG TCCAACAAGCGATTTA'
#
def getInsertSites(s):
	if s.count('(') == 1:
		p = re.compile(r'.+\(\s*[A-Z,a-z]*\s*\).+')
	else:
		p = re.compile(r'[A-Z,a-z]+\s*\(\s*[A-Z,a-z]*\s*\)\s*[A-Z,a-z]+')
	insertSites = [site.replace(' ', '') for site in p.findall(s)]
	return insertSites

def specifyPep(isPep, id):
	return isPep[id-1]

# Get all IS features including IS_SEQ and IS_PEP
# Return: mfamilyFeatures
# mfamilyFeatures: [(familyName, familyFeatures), ..., (familyName, familyFeatures)]
# familyFeatures: [isFeatures, ..., isFeatures]
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
def getFeatures(families, fp):
	mfamilyFeatures = []
	for family in families:
		# get family name from the first member of family
		familyName = family[0][0][2][2].strip().upper()
		familyFeatures = []
		for element in family:
			features = {}
			isName = element[0][0][0].strip().upper()
			group = element[0][2][3].strip().upper()
			origin = element[1][1][0].strip().upper()
			features['isName'] = isName
			features['group'] = group
			features['origin'] = origin
			features['familyName'] = familyName

			# length of IS.
			#
			# Note for IS finder submission:
			# LENGTH:
			# Should only be included if the ends of the element can be defined either by transposition 
			# or by comparison with other elements. It should not include the direct target repeats 
			# (green triangles). In case of any doubt, ND should be added to this field.
			#
			isLen = element[2][1][1].strip().strip('><?').strip()
			if isLen.isdecimal():
				isLen = int(isLen)
			else:
				print('Warning: unknown lenght of IS element #{}# {} {}'.format(isLen, familyName, isName))
				isLen = -1
			features['isLen'] = isLen

			# lenght of IR
			#
			# Note for IS finder submission:
			# INVERTED REPEATS: 
			# The length of the inverted repeats (red), if any, should be included. If there are 
			# no missmatches then a single number is sufficient. If there are a few differences then 
			# these should be included as #identical/#total. If no inverted repeats could be found, 
			# insert "0" (zero). If the information is not available (for example a partial sequence 
			# or a truncated element) a "-" (dash) should be inserted. 
			#
			ir = element[2][1][2].strip()
			if '/' in ir:
				irLen = ir.split('/', 1)
				irLen = [int(irLen[0].strip()), int(irLen[1].strip())]
			elif ir.isdecimal():
				irLen = [int(ir)]
			else:
				if ir != '-':
					print('Warning: unknown number for the length of IR #{}# {} {}'.format(ir, familyName, isName))
				irLen = [-1]
			features['irLen'] = irLen

			# length of DR
			#
			# Note for IS finder submission:
			# DIRECT REPEATS: 
			# Direct repeats (green triangles) should be defined as follows: There can be NO mismatches. 
			# Great care must be taken whenever symmetrical or palindromic terminal repeats are present. 
			# In case of doubt, this should be clearly stated in the comments field. If no direct repeats 
			# could be found, insert "0" (zero). If the information is not available (for example a partial 
			# sequence or a truncated element) a "-" (dash) should be inserted. 
			#
			dr = element[2][1][3].strip()
			if dr.isdecimal():
				drLen = int(dr)
			else:
				drLen = -1
				if dr != '-':
					print('Warning: unknown number for the length of DR #{}# {} {}'.format(dr, familyName, isName))
			features['drLen'] = drLen

			# activity of transposition.
			#
			# Note for IS finder submission:
			# TRANSPOSITION: 
			# Y = yes; ND = not determined; I = inactive element.
			#
			act = element[3][0][1].strip().lower()
			if not (act == 'y' or act == 'i'):
				act = '-'
				print('Warning: unknown transposition activity of IS element #{}# {} {}'.format(act, familyName, isName))
			features['act'] = act


			# ORFs
			#
			# Note for IS finder submission:
			# ORF: 
			# When more than one consecutive reading frame is present, they should be entered in their order 
			# of occurrence from left to right. Length is in amino acids and does not include the stop codon. 
			# Start and end are in nucleotides and stand for the first nucleotide of the start codon and the 
			# terminal nucleotide of the stop codon. In cases where transposase is synthesised by programmed 
			# translational frameshifting, the "start" of the distal frame should be taken as the first 
			# nucleotide following a stop codon. 
			#
			orfsStr = element[3][1][1].strip().upper()
			# There are many irregular records which are not compatible with format recommended by ISfinder
			# submission rules above.
			# We can check if there is a number in raw data record, and then check if there is a pair of 
			# round brackets, in order that we can get the 'start' and 'end' location of ORFs in IS element.
			# 
			if tools.hasNumbers(orfsStr) and tools.hasBrackets(orfsStr):
				# [] or [[start, end], ..., [start, end]]
				#print(orfsStr, familyName, isName)
				orfs = getOrfs(orfsStr)
			else:
				orfs = []
			if len(orfs) == 0:
				print('Warning: irregular or empty ORFs #{}# {} {}'.format(orfsStr, familyName, isName), file=fp)
			features['orf'] = orfs
			

			# sequences of Left End and Right End
			#
			# Note for IS finder submission:
			# SEQUENCE: 
			# By convention, the left end of the IS (IRL - red) is defined as that proximal to the transposase 
			# promoter. The sequence should be entered without flanking direct target repeats (DR- green) with 
			# the orientation of transcription/translation of the transposase gene (yellow rectangle) from left 
			# to right. 
			#
			# Note that the 50 terminal base pairs of the DNA sequence are automatically copied into the 
			# appropriate fields in the correct relative orientations.
			#
			lEnd, rEnd = element[4][0][1].strip(), element[4][1][1].strip()

			lEnd = tools.cleanDNA(lEnd)
			#if lEnd.isalpha() and len(lEnd) == 50:
			if lEnd.isalpha():
				lSeq = lEnd
				if len(lEnd) < 50:
					print('Warning: incomplete Left End (<50) #{}# {} {}'.format(lEnd, familyName, isName))
			else:
				lSeq = ''
				print('Warning: irregular sequence code in Left End #{}# {} {}'.format(lEnd, familyName, isName))
			features['lSeq'] = lSeq
			rEnd = tools.cleanDNA(rEnd)
			#if rEnd.isalpha() and len(rEnd) == 50:
			if rEnd.isalpha():
				rSeq = rEnd
				if len(rEnd) < 50:
					print('Warning: incomplete Right End (<50) #{}# {} {}'.format(rEnd, familyName, isName))
			else:
				rSeq = ''
				print('Warning: irregular sequence code in Right End #{}# {} {}'.format(rEnd, familyName, isName))
			features['rSeq'] = rSeq

			# Insertion sites
			#
			# Note for IS finder submission:
			# INSERTION SITES: 
			# When defined by transposition events (i.e. knowledge of the target sequence before and after 
			# insertion) or where the perfect direct repeat is unambiguous, the insertion site should be 
			# reported using the convention: N10 ( DR ) N10 
			# e.g. GGATCGAATG ( TGC ) GATACCCTGC 
			#
			s = element[5][0][1].strip()
			if s == '-' or s == '' or s.isspace():
				print('Warning: no INSERTION SITES', familyName, isName)
			elif '(' in s and ')' in s:
				insertSites = getInsertSites(s)
				if len(insertSites) == 0:
					print('Warning: unknown Insert Sites #{}# {} {}', s, familyName, isName)
			else:
				print('Warning: unknown Insert Sites #{}# {} {}', s, familyName, isName)


			# IS sequence (without DRs)
			#
			# Note for IS finder submission:
			# SEQUENCE: 
			# By convention, the left end of the IS (IRL) is defined as that proximal to the transposase 
			# promoter. The sequence should be entered without flanking direct target repeats (DR) with 
			# the orientation of transcription/translation of the transposase gene from left 
			# to right. 
			#
			if element[6][0][0][:6].strip().upper() == 'IS_SEQ':
				isSeq = element[6][1][0]
			else:
				isSeq = ''
			features['isSeq'] = isSeq

			# IS peptide sequence
			# 
			# Note for IS finder submissions:
			# IS PEP: 
			# When more than one consecutive reading frame is present, they should be entered in their order of 
			# occurrence from left to right. 
			# The amino acid sequence of each translated reading frame should be included using the single 
			# letter code. 
			#
			isPep = []

			# IS_PEP record is available.
			if element[7][0][0][:6].strip().upper() == 'IS_PEP':
				for pep in element[7:-2]:
					isPep.append(pep[1][0])
				# count number of peptides in IS elements
				if pep[0][0][6].isdigit():
					nOrfByPep = int(pep[0][0][6: pep[0][0].rfind(':')])
				else:
					nOrfByPep = 1

			# ignore the IS element without peptide sequence
			if len(isPep) < 1:
				continue

			features['isPep'] = isPep
			# number of IS_PEP records
			features['nIsPep'] = len(isPep)

			features['nOrfByPep'] = nOrfByPep

			# Collect peptide sequences for the transposase of IS element. 
			# The tpase is usually the lonegest peptide in IS_PEP records but there are special cases where tpase
			# is not the longest peptide among IS_PEPs.
			# Note: family IS200/IS605 has 3 groups, IS200 which usually has one short tpase ORF,
			#       IS1341 which usually has one long accessory gene ORF except ISHHU11 and ISNAMO21 (both have 
			#	one short passenger gene and one long accessory gene),
			#       IS606 which usually has one short tpase ORF and one long accessory gene ORF.
			#
			'''
			if features['familyName'] == 'IS200/IS605':
				if (features['group'] == 'IS1341' or features['group'] == 'IS200') and len(isPep) != 1:
					print('Warning: {} {} {} has {} peptides in ISfinder'.format(
						features['familyName'], features['group'], isName, len(isPep)))
				elif features['group'] == 'IS606' and len(isPep) != 2:
					print('Warning: {} {} {} has {} peptides in ISfinder'.format(
						features['familyName'], features['group'], isName, len(isPep)))
				elif len(isPep) > 1:
					print('Warning: {} {} {} has {} peptides in ISfinder'.format(
						features['familyName'], features['group'], isName, len(isPep)))

				if features['group'] == 'IS1341':
					#features['tpase'] = ''
					features['tpase'] = max(isPep, key=len)
				else: # IS200, IS606, others without groupname assigned
					#features['tpase'] = min(isPep, key=len)
					features['tpase'] = ''.join(isPep)
			elif isName in shortTpases:
			'''
			if isName in shortTpases:

				features['tpase'] = specifyPep(isPep, shortTpases[isName])
			else: # tpase is the longest peptide.
				features['tpase'] = max(isPep, key=len)

			# Comments: we try to retrieve valueable information from Comments record.
			# 
			# Note for IS finder submission:
			# COMMENTS: 
			# The following information will be welcome: 
			# • For ISs from genome sequences, the genome co-ordinates. 
			# • For unusual species, some indication of their behaviour, pathogenicity, particular environmental niche. 
			# • Some idea of the degree of relatedness to other key members of the family. 
			# • Association with other ISs or particular genes. 
			# • Unusual structural or organisational features (e.g.internal repeated elements) 
			# • Types of transposition product (e.g. cointegrate, direct insertions...) 
			# • Copy number 
			# • Occurrence (i.e. strains tested which proved negative) 
			#
			comments = element[-2][1][0].strip()
			nOrfByComment = getOrfByComments(comments)
			features['nOrfByComment'] = nOrfByComment

			if nOrfByComment == 0:
				if len(orfs) !=  nOrfByPep:
					print('Warning: number of ORFs may be unclear, nOrfs={} nOrfByPep={} nOrfByComment={} {} {}'.format(len(orfs), nOrfByPep, nOrfByComment, familyName, isName))
			elif max(len(orfs), nOrfByPep, nOrfByComment) != min(len(orfs), nOrfByPep, nOrfByComment):
				print('Warning: number of ORFs may be unclear, nOrfs={} nOrfByPep={} nOrfByComment={} {} {}'.format(len(orfs), nOrfByPep, nOrfByComment, familyName, isName))

			familyFeatures.append(features)
		mfamilyFeatures.append((familyName, familyFeatures))

	return mfamilyFeatures

# Remove IS with incorrect isLen/irLen.
def cleanFeatures(mfamilyFeatures):
	newmfamilyFeatures = []
	for family in mfamilyFeatures:
		newFamily = []
		familyName = family[0]
		for element in family[1]:
			isName = element['isName']
			isLen, irId, irLen = element['isLen'], element['irLen'][0], element['irLen'][-1]
			if 0 < isLen < irLen:
				print('Warning: remove IS element with weird isLen {}, irLen {}', isLen, irLen, 
												familyName, isName)
				continue
			newFamily.append(element)
		newmfamilyFeatures.append(newFamily)
	return mfamilyFeatures

# - exchange values of isLen and irLen when isLen < irLen
def cleanFeatures1(mfamilyFeatures):
	for family in mfamilyFeatures:
		for element in family[1]:
			isLen, irId, irLen = element['isLen'], element['irLen'][0], element['irLen'][-1]
			if 0 < isLen < irLen:
				if len(element['irLen']) > 1:
					print('Error: incorrect isLen {}, isId {}, isLen {} {} {}'.format(isLen, irId, irLen,
						family[0], element['isName']), file=sys.stderr)
					exit(0)
				else:
					print('Warning: isLen {} and irLen {} is exchanged {} {}'.format(isLen, irLen, 
						family[0], element['isName']))
					element['isLen'], element['irLen'][0] = irLen, isLen
	return mfamilyFeatures

def long_pep_v1(is_pep):
	maxlen = 0
	longpep = ''
	for item in is_pep:
		if 'IS_PEP' in item[0][0]:
			length = len(item[1][0])
			if length > maxlen:
				maxlen = length
				longpep = item[1][0]
	return longpep.translate(str.maketrans(' ', '\n'))

def long_pep_v2(is_pep):
	maxlen = 0
	longpep = ''
	for item in is_pep:
		if 'IS_PEP' in item[0][0]:
			length = len(item[1][0])
			if length > maxlen:
				maxlen = length
				longpep = item[1][0]
	return longpep

def longPep(isPep):
	maxLen = 0
	longPep = ''
	for pep in isPep:
		pepLen = len(pep)
		if pepLen > maxLen:
			maxLen = pepLen
			longPep = pep
	return longPep

def getLongPep1(mfamilyFeatures):
	pep4msa = []
	for family in mfamilyFeatures:
		familyName = family[0]
		familyPep = []
		for element in family[1]:
			elementPep = longPep(element['isPep'])
			if elementPep != '':
				familyPep.append([familyName, element['group'], element['isName'], element['origin'], elementPep])
		pep4msa.append(familyPep)
	return pep4msa

def getLongPep(families):
	pep4msa = []
	for family in families:
		family_pep = []
		for element in family:
			element_pep = long_pep_v2(element[6:-2])
			if element_pep != '':
				if element[0][2][3] == '' or element[0][2][3] == '-':
					groupname = '-'
				else:
					groupname = element[0][2][3].strip()
				familyname = element[0][2][2].strip()
				elementname = element[0][0][0].strip()
				originname = element[1][1][0].strip()
				family_pep.append([familyname, groupname, elementname, originname, element_pep])
		pep4msa.append(family_pep)
	return pep4msa

def readClustersByCDhit(clusterFile):
	clusters = []
	fp = open(clusterFile, 'r')
	i = 0
	cluster = []
	for line in fp:
		if line[0] == '>':
			if len(cluster) > 0:
				clusters.append((i, cluster))
				cluster = []
				i += 1
			continue

		index2 = line.rfind('|')
		left = line[:index2]
		index1 = left.rfind('|')
		isName = left[index1+1:]
		cluster.append(isName)
	clusters.append((i, cluster))
	fp.close()
	return clusters

def writeFamilyPepFile(family):
	familyName = family[0]
	familyNameOnDisk = familyName.replace('/', '_')
	#familyPepFile = os.path.join(constants.rootpath, familyNameOnDisk, familyNameOnDisk + '.long.pep.all')
	familyPepFile = os.path.join(constants.rootpath, familyNameOnDisk, familyNameOnDisk + '.faa')

	# output multiple sequences in fasta format into a file, one file each family
	family_fasta = ''
	for element in family[1]:
		if 'tpase' in element.keys():
			pepSeq = element['tpase']
		else:
			continue
		if len(pepSeq) == 0:
			continue

		groupname =  element['group']
		elementname = element['isName']
		originname =  element['origin']
		#header = '_'.join((familyName, groupname, elementname)) + '| ' + originname
		header = '|'.join((familyName, groupname, elementname)) + '| ' + originname
		element_fasta = tools.fasta_format(header, pepSeq)

		family_fasta += element_fasta + '\n'

	#if not os.path.isfile(familyPepFile):
	fp = open(familyPepFile, 'w')
	fp.write(family_fasta)
	fp.close()
	return familyPepFile

# mfamilyFeatures: [(familyName, familyFeatures), ..., (familyName, familyFeatures)]
# familyFeatures: [isFeatures, ..., isFeatures], elements in family
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}, one element
def writeFamilySeqfile(family):
	familyName, familyfeatures = family
	familyNameOnDisk = familyName.replace('/', '_')
	familyPepFile = os.path.join(constants.rootpath, familyNameOnDisk, familyNameOnDisk + '.fna')
	filec = []
	for isfeatures in familyfeatures:
		seqid = '|'.join([familyName, isfeatures['group'], isfeatures['isName']])
		des = isfeatures['origin']
		header = ' '.join([seqid, des])
		fasta = tools.fastaFormat(header, isfeatures['isSeq'])
		filec.append(fasta)
	return '\n'.join(filec)

def outputISseq(mfamilyFeatures):
	filec4isfna = []
	for family in mfamilyFeatures:
		filec4isfna.append(writeFamilySeqfile(family))
	with open('isfinder.is.fna', 'w') as fp:
		fp.write('\n'.join(filec4isfna) + '\n')

def callCDhit4hierarch(familyPepFile, lowIdent=False):
	cdhit = cdhit
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
	cmd_line = ' '.join([clstr_rev, clstrFile90, clstrFile60])
	do_cmd = shlex.split(cmd_line)
	fp = open(clstrFile9060, 'w')
	subprocess.check_call(do_cmd, shell=False, universal_newlines=False, stdout=fp)
	fp.close()

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
	cmd_line = ' '.join([psicdhit, options])
	do_cluster = shlex.split(cmd_line)
	subprocess.check_call(do_cluster, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL, cwd=os.path.dirname(input))

	# clstr_rev.pl nr9060.clstr nr30.clstr > nr906030.clstr
	clstrFile30 = '.'.join([output, 'clstr'])
	clstrFile906030 = '.'.join([clstrFile9060, os.path.basename(clstrFile30)])
	cmd_line = ' '.join([clstr_rev, clstrFile9060, clstrFile30])
	do_cmd = shlex.split(cmd_line)
	fp = open(clstrFile906030, 'w')
	subprocess.check_call(do_cmd, shell=False, universal_newlines=False, stdout=fp)
	fp.close()

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
	familyPepFile = writeFamilyPepFile(family)

	# example: ISH3.faa.nr100.nr90.clstr.ISH3.faa.nr100.nr90.nr60.clstr.ISH3.faa.nr100.nr90.nr60.nr30.clstr
	dir, basename = os.path.split(familyPepFile)
	nr90 = '.'.join([basename, 'nr100.nr90.clstr'])
	nr9060 = '.'.join([basename, 'nr100.nr90.nr60.clstr'])
	nr906030 = '.'.join([basename, 'nr100.nr90.nr60.nr30.clstr'])

	lowIdent = lowIdent
	if lowIdent == False:
		clusterFile = '.'.join([dir, nr90, nr9060])
	else:
		clusterFile = '.'.join([dir, nr90, nr9060, nr906030])

	if os.path.isfile(clusterFile):
		print('Skip (cd-hit and psi-cd-hit) clustering {}'.format(familyPepFile))
	else:
		clusterFile = callCDhit4hierarch(familyPepFile, lowIdent)

	clustersByCDhit = readClustersByCDhit(clusterFile)
	clusters = []
	isNames = {isFeatures['isName']:isFeatures for isFeatures in family[1]}
	for g in clustersByCDhit:
		clusterName = family[0] + '_' + str(g[0])
		cluster = []
		commonNames = set(isNames.keys()) & set(g[1])
		for isName in commonNames:
			isFeatures = isNames[isName]
			isFeatures['clusterName'] = clusterName
			cluster.append(isFeatures)
		clusters.append((clusterName, cluster))
		
	return clusters

# Return all clusters of IS elements where clustering is based on transposase similarity.
# Return mCluster
# mCluster: [cluster, ..., cluster]
# cluster: (clusterName, members)
# members: [isFeatures, ..., isFeatures]
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
#
# mfamilyFeatures: [(familyName, familyFeatures), ..., (familyName, familyFeatures)]
# familyFeatures: [isFeatures, ..., isFeatures]
# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
#
def getCluster(mfamilyFeatures):
	mCluster = []
	for family in mfamilyFeatures:
		# clusters: [cluster, ..., cluster]
		# cluster: [isFeatures, ..., isFeatures]
		#mCluster.extend(doClusterByGroup(family[1]))
		#mCluster.extend(doClusterByElement(family[1]))
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

# get the lengths of IS element, IR, DR, and the activity of transposition
def getGroupIrLen(element):
	irLen = element['irLen'][-1]
	return (element['group'], irLen)

# orfs: [] or [[start, end], ..., [start, end]]
def getDistOrfs(orfs):
	dist = []
	if len(orfs) < 2:
		return dist
	if len(orfs) == 3 and (orfs[2][0] == orfs[0][0] and orfs[2][1] == orfs[1][1]):
		return dist

	# The orfs in complementary strand will be written as (start, end)cs in ISfinder,
	# Refer to ISStma11 and ISCausp2 for the different format of orfs in complementary strand. 
	#
	# For example, ISStma11 has the orf string: 
	# 446(144-1484)564(3261-1567) CS92(3550-3272) CS116(3913-3566) CS135(3988-4395)
	# Its orfs is [[144,1484],[3261,1567],[3550,3272], [3913,3566], [3988,4395]].
	# Among its five orfs, [3261,1567],[3550,3272], [3913,3566] are in complementary 
	# strand and their start and end coordinates are registered as the reverse order, 
	# namely, start > end. We reversed the start and end coordinates of such orfs in 
	# complementary strand for the convenience in the following distance caculation.
	#
	# Regarding ISCausp2, the start and end coordinates of orfs in complementary strand 
	# are not registered as the reverse order.
	orfsCS = [(min(orf), max(orf)) for orf in orfs]

	left, right = [], []
	neighbors = list(itertools.chain.from_iterable(orfsCS))[1:-1]
	for item in zip(*[iter(neighbors)]*2):
		dist.append(item[1]-item[0]-1)
	return dist

# summarize the features in each family
def featureStat(mfamilyFeatures, fp):
	#fmtStr = '{:<15} {:<13} {:<13} {:>7} {:>7} {:>5} {:>5} {:>5} {:>3} {:>10} {:>10} {:>10} {:>10} {:>4} {:>4} {:>9} {:>6}'
	#fmtStr = '{:<15} {:<13} {:>7} {:>7} {:>5} {:>5} {:>5} {:>3} {:>10} {:>10} {:>10} {:>10} {:>4} {:>4} {:>9} {:>6}'
	fmtStr = '{:<15} {:<13} {:>7} {:>7} {:>7} {:>5} {:>5} {:>5} {:>3} {:>10} {:>10} {:>10} {:>10} {:>4} {:>4} {:>9} {:>6}'
	print('', file=fp)
	print('# summary of the isfinder features', file=fp)
	#print(fmtStr.format('family', 'IS', 'group', 'seqLen', 'isLen', 'irId', 'irLen', 'drLen', 'act', 
	#print(fmtStr.format('family', 'IS', 'seqLen', 'isLen', 'irId', 'irLen', 'drLen', 'act', 
	print(fmtStr.format('cluster', 'IS', 'pepLen', 'seqLen', 'isLen', 'irId', 'irLen', 'drLen', 'act', 
				'lOrf2Ter', 'rOrf2Ter', 'lLinkerLen', 'rLinkerLen', 'nOrf', 'orfs', 'nOrfByPep', 'nIsPep'), file=fp)
	print('-' * 67, file=fp)
	dist = {}
	for family in mfamilyFeatures:
		#family[1].sort(key = getGroupIrLen)
		#family[1].sort(key = operator.itemgetter('group', 'drLen'))
		family[1].sort(key = operator.itemgetter('group', 'act'))
		for element in family[1]:
			familyName, isName, groupName = family[0], element['isName'], element['group']
			#familyName, isName, groupName = family[0], element['isName'], element['group']
			#familyName, clusterName = familyName.rsplit('_', 1)
			isLen, drLen = element['isLen'], element['drLen'] 
			if len(element['irLen']) > 1:
				irId, irLen = element['irLen'][0], element['irLen'][1]
			else:
				irId = irLen = element['irLen'][0]
			act = element['act']
			orfs = element['orf']

			# distance between neighboring orfs in one IS element
			dist[(familyName, isName)] = getDistOrfs(orfs)

			# amino acid sequence length of Tpase in IS element, based on IS_PEP
			pepLen = len(element['tpase'])

			# nucleic acid sequence length of IS element, based on IS_SEQ
			seqLen = len(element['isSeq'])

			# Estimated the length of linker between IR and ORF when assuming IRs alwasy start from termini 
			# of IS element.
			# Note: IR does not always starts from terminus though it is true for most IS elements. 
			if len(orfs) == 0:
				lOrfMin = rOrfMax = -1
			else:
				lOrfMin = min([min(item) for item in orfs])
				rOrfMax = max([max(item) for item in orfs])
			elementLen = max(isLen, seqLen)
			#if lOrfMin > elementLen or rOrfMax > elementLen:
			#	continue

			lOrf2Ter = lOrfMin - 1
			lLinkerLen = lOrf2Ter - irLen

			rOrf2Ter = elementLen - rOrfMax
			rLinkerLen = rOrf2Ter - irLen

			# number of ORFs
			nOrf = max(len(orfs), element['nOrfByPep'], element['nOrfByComment'])

			#print(fmtStr.format(familyName, clusterName, isName, #groupName, 
			print(fmtStr.format(familyName, isName, #groupName, 
						pepLen, seqLen, isLen, irId, irLen, drLen, act, 
						lOrf2Ter, rOrf2Ter, lLinkerLen, rLinkerLen, 
						nOrf, len(orfs), element['nOrfByPep'], element['nIsPep']), 
						file=fp)
	fmtStrOrfDist = '{:<15} {:<13} {:>7}'
	print('Distance between neighboring orfs', file=fp)
	print(fmtStrOrfDist.format('familyName', 'isName', 'distance'), file=fp)
	print('-' * 67, file=fp)
	for accid in dist:
		for item in dist[accid]:
			print(fmtStrOrfDist.format(accid[0], accid[1], item), file=fp)
	
def buildTpaseHMM(args):
	tpaseSeqFile = args['tpaseSeqFile']

	logFile = '.'.join(args['program'], 'log'])
	fp_log = open(logFile, 'w')
	cmdline = args['cmdline']
	print(cmdline + '\n' + datetime.datetime.now().ctime() + '\n')
	fp_log.write(cmdline + '\n' + datetime.datetime.now().ctime() + '\n\n')

	# get tranposase sequences from fasta file
	mfamilyFeatures = getTpase(tpaseSeqFile)

	
	
	# Based on Tpase peptide sequence, do cluster analysis on IS elements within same family
	# in order to divide the diverse family into multiple subfamilies or groups. 

	# mfamilyFeatures: [family, ..., family]
	# family: (familyName, familyFeatures)
	# familyFeatures: [isFeatures, ..., isFeatures]
	# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
	#
	# mCluste: [cluster, ..., cluster]
	# cluster: (clusterName, members)
	# members: [isFeatures, ..., isFeatures]
	# isFeatures: {'isName': isName, ..., 'nOrfByComment': nOrfByComment}
	#
	#mfamily_seq = getLongPep(families)
	#mfamily_seq = getLongPep1(mfamilyFeatures)
	if lowIdent == False:
		print('Clustering each family at sequence identity 60%')
	elif lowIdent == True:
		print('Clustering each family at sequence identity 30%')
	mCluster = getCluster(mfamilyFeatures)
	mfamilyFeatures = mCluster

	outputClusters(mCluster, fpIsfinderTable)

	# output summary of isfinder features
	featureStat(mfamilyFeatures, fpIsfinderTable)
	fpIsfinderTable.close()
	print('Write ISfinder feature summary into', isfinderFeaturesFile)

	print('\nFinish retrieving and summarizing features of IS elements from isfinder database')
	print(datetime.datetime.now().ctime() + '\n')

	#pep_tag = 'long.pep'
	#seq_profile_tag = pep_tag + '.' + profile_tag
	args2concurrent4seq = []

	DOphmmer = False
	DOhmmsearch = False
	fasta4singleMemberClusters = []
	for cluster in mCluster:
		# the first member of cluster
		familyname = cluster[1][0]['familyName']

		familyNameOnDisk = familyname.replace('/', '_')
		family_path = os.path.join(rootpath, familyNameOnDisk)

		clusterName = cluster[0]
		clusterNameOnDisk = clusterName.replace('/', '_')

		#seq_filename = clusterNameOnDisk + '.' + seq_profile_tag
		#msa_filename = clusterNameOnDisk + '.' + seq_profile_tag + '.msa'
		seq_filename = '.'.join([clusterNameOnDisk, 'faa'])
		msa_filename = '.'.join([clusterNameOnDisk, 'faa.msa'])
		seq_pathfilename = os.path.join(family_path, seq_filename)
		msa_pathfilename = os.path.join(family_path, msa_filename)

		# output multiple sequences in fasta format into a file, one file each family or cluster
		cluster_fasta = []
		nseq = 0
		for element in cluster[1]:
			if 'tpase' in element.keys():
				pepSeq = element['tpase']
			else:
				continue
			if len(pepSeq) == 0:
				continue

			familyname = element['familyName']
			groupname =  element['group']
			elementname = element['isName']
			originname =  element['origin']
			#header = '_'.join((clusterName, familyname, groupname, elementname)) + '| ' + originname
			header = '|'.join((clusterName, familyname, groupname, elementname)) + '| ' + originname
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
		fp_seq_pathfilename = open(seq_pathfilename, 'w')
		fp_seq_pathfilename.write('\n'.join(cluster_fasta))
		fp_seq_pathfilename.close()
		print('Write peptide sequences from cluster {} of family {} into {}'.format(
				clusterName, familyname, seq_pathfilename))

	# sequence file holding multiple peptide sequences used as input for phmmer
	#clustersSeqFile4phmmer = os.path.join(rootpath, '.'.join(['clusters', seq_profile_tag, 'faa']))
	clustersSeqFile4phmmer = os.path.join(rootpath, '.'.join(['clusters', 'single', 'faa']))
	fp_clustersSeqFile4phmmer = open(clustersSeqFile4phmmer, 'w')
	fp_clustersSeqFile4phmmer.write('\n'.join(fasta4singleMemberClusters))
	fp_clustersSeqFile4phmmer.close()
	print('Write peptide sequences from single-member clusters into', clustersSeqFile4phmmer)

	# Build multiple alignment and HMM models
	#hmms_file = os.path.join(rootpath, '.'.join(['clusters', pep_tag, profile_tag, 'hmm']))
	hmms_file = os.path.join(rootpath, '.'.join(['clusters', 'faa', 'hmm']))
	print('Building multiple alignment and HMM models.')
	buildMSA(args2concurrent4seq, constants.nproc, fp_log)
	buildHMM(args2concurrent4seq, hmms_file, fp_log)

	fp_log.close()
	print('Log into', logFile)
	print('\n' + datetime.datetime.now().ctime() + '\n')


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

	buildTpaseHMM(args)

