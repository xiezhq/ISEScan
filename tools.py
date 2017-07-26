import io
import csv
import glob
import os.path
import constants
import re
import sys
import itertools
import operator # for getbds4opt4start()
import subprocess, shlex
import errno # for makedir()


# check character string
def is_None_empty_whitespace(mystr):
	if mystr and mystr.strip():
		# mystr is not None AND mystr is not empty or whitespaces
		return False
	# mystr is None OR mystr is empty or whitespaces
	return True

# This code makes the assumption that anything that can be iterated over will contain other elements, and should not be 
# considered a leaf in the "tree". If an attempt to iterate over an object fails, then it is not a sequence, and hence 
# certainly not an empty sequence (thus False is returned). Finally, this code makes use of the fact that all returns 
# True if its argument is an empty sequence.
# Note: seq should not contain any character string except ''.
def isEmpty(seq):
	try:
		return all(map(isEmpty, seq))
	except TypeError:
		return False


def read_file(file):
	fp = open(file, 'r')
	content = fp.read()
	fp.close()

	return content

# Write content into file filePath
def write2file(filePath = None, content = None):
	directory = os.path.dirname(filePath)
	if not os.path.exists(directory):
		try:
			os.makedirs(directory)
		except OSError as error:
			if error.errno != errno.EEXIST:
				raise
	with open(filePath, 'w') as fp:
		fp.write(content)

# Split an multiple-sequence fasta file containing multiple fasta sequences into multiple individual fasta files
# with the file names with sequence id included.
def split_tandem_fasta(huge_fasta_file, output_path):
	mfastaFileName = os.path.basename(huge_fasta_file)
	fp_fasta = open("/dev/null", "r")
	with open(huge_fasta_file, "r") as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			if line[0] == '>':
				fp_fasta.close()
				seqid = line[1:].split(maxsplit=1)[0]
				fasta_file_name = '.'.join([mfastaFileName, seqid])
				fasta_file = os.path.join(output_path, fasta_file_name)
				fp_fasta= open(fasta_file, "w")
			fp_fasta.write(line+'\n')
	fp_fasta.close()


# This function returns a generator, using a generator comprehension. The generator returns the string sliced, 
# from 0 + a multiple of the length of the chunks, to the length of the chunks + a multiple of the length of the chunks.
# example:
# s = 'CATTCGTCT', tuple(chunkstring(s, 3)) = ('CAT', 'TCG', 'TCT')
# s = 'CATTCGTC', tuple(chunkstring(s, 3)) = ('CAT', 'TCG', 'TC')
# 
def chunkstring(string, length):
	return (string[i : i + length] for i in range(0, len(string), length))


# Translate gene sequence into peptide sequence.
# Return: 'MKLQFP.....RPQ', peptide sequence represented by single letter
#
# gene2pepTable: (table1, ..., table11, ...)
# table11: {'starts': (initialcodon, ..., initialcondon), codon: aa, ..., codon: aa}
# codon: 'TTT' or 'TTC' or 'TTG' or ... or 'TGA' or 'GGG'
# aa: 'F' or 'L' or ... or '*' or 'G', amino acid when codon is not the initial codon
#
# geneSeq: gene sequence eg. 'CATTCGTCT', it must be the length of 3*integer 
#	or the last one (or two) letter will be removed.
# 
def gene2pep(table='11', geneSeq=''):
	if len(geneSeq) < 9:
		print('Warning: mini ORF with length < 9 will not be translated', geneSeq)
		return ''

	pep = []
	first3 = geneSeq[:3]
	geneCodons = chunkstring(geneSeq[3:], 3)

	# translate the first codon in gene sequence
	if first3 in constants.gene2pepTable[table]['starts']:
		# it is a start codon
		pep.append('M')
	else:
		# it is not a start condon
		aa = constants.gene2pepTable[table][first3]
		if aa == '*':
			print('hello stop codon', first3)
		else:
			pep.append(aa)
	
	# stop translation if the first codon is the stop codon
	#if '*' in pep:
	#	return '*'

	# translate the left codon but stop translation if stop codon is found.
	for codon in geneCodons:
		# the last codon might be a single-letter string or double-letter string.
		if len(codon) != 3:
			continue

		aa = constants.gene2pepTable[table][codon]
		if aa == '*':
			print('hello stop codon', codon)
		#	break
		else:
			pep.append(aa)
	return ''.join(pep)
	

# Format sequence into FASTA format and return it with header
# header: messeage for header line in fasta format
# seq: sequence to format into FASTA with lineWidth single-letter codes per line
def fasta_format(header, seq):
	header = '>' + header
	fasta_seq = ''
	lineWidth = constants.fastaLineWidth
	#for row in chunkstring(seq, lineWidth):
	#	fasta_seq += row + '\n'

	fasta_seq = '\n'.join(chunkstring(seq, lineWidth))
	fasta = '{}\n{}\n'.format(header, fasta_seq)
	#fasta = header + '\n' + fasta_seq
	return fasta

# header, seq: character string
def fastaFormat(header, seq):
	header = '>' + header
	fasta_seq = ''
	fasta_seq = '\n'.join(chunkstring(seq, constants.fastaLineWidth))
	fasta = '\n'.join([header, fasta_seq])
	return fasta

# get both compound identifier and content of fasta sequence
# seqs: [(id, seq), ..., (id, seq)]
# seq: character string
def getFasta(fastaFile):
	fp = open(fastaFile, 'r')
	seqs = []
	seq = []
	id = ''
	for line in fp:
		if line.replace('\b','').strip() == '':
			continue # remove blank line
		if line[0] == '>':
			if len(seq) > 0:
				seqs.append((id, ''.join(seq)))
			id = line[1:].split(maxsplit=1)[0]
			seq = []
			continue
		seq.append(line.strip())
	if len(seq) > 0:
		seqs.append((id, ''.join(seq)))
	else:
		print('No sequence in', fastaFile)
	fp.close()
	return seqs

# seqs: [(id, seq), ..., (id, seq)]
# seq: character string
def getFasta_idseq(fastaFile):
	fp = open(fastaFile, 'r')
	seqs = []
	seq = []
	id = ''
	for line in fp:
		if line.strip() == '':
			continue # remove blank line
		if line[0] == '>':
			if len(seq) > 0:
				seqs.append((id, ''.join(seq)))
				seq = []
			id = line[1:].split(maxsplit=1)[0]
			continue
		seq.append(line.strip())
	if len(id) > 0:
		seqs.append((id, ''.join(seq)))
	else:
		print('No sequence in', fastaFile)
	fp.close()
	return seqs

# get both header and content of fasta sequence
# seqs: [(header, seq), ..., (header, seq)]
# seq: character string
def getFastaFull(file):
	seqs = []
	with open(file, 'r') as fp:
		seq = []
		for line in fp:
			#line = line.replace('\b', '').strip()
			line = line.strip()
			if line == '':
				continue # remove blank line
			if line[0] == '>':
				if len(seq) > 0:
					seqs.append((header, ''.join(seq)))
					seq = []
				header = line[1:] # header
				continue
			seq.append(line)
		seqs.append((header, ''.join(seq)))
	return seqs

# write IS elements in a genome into a csv file
# genome_is: list, [] or [['ISDge6', 'IS5', 'ISL2', '-', '-', '-', '-', '36785', '36563', '802', '1024', '802', 'Partial'],
#		['ISDge2', 'IS1', '.', '+', '37705', '38452', '748', '-', '-', '-', '-', '', '-'], ...]
# write a empty file if genome_is is empty
def output_csv(csvfile, genome_is):
	with open(csvfile, 'w', newline='') as f:
		writer = csv.writer(f, lineterminator = '\n', quoting = csv.QUOTE_MINIMAL)
		writer.writerows(genome_is)


# get isfinder genome annotation from csv files
def isfinder_IS_in_genome(isfinder_genome_file):
	isfinder_genome = []

	if not os.path.isfile(isfinder_genome_file):
		return isfinder_genome

	fp = open(isfinder_genome_file, newline='')
	reader = csv.reader(fp, delimiter=',')
	for row in reader:
		isfinder_genome.append(row)
	fp.close()
	return isfinder_genome

# Group elements in a sequence by keys indexed by (i, ..., j).
# sequence: (seq1, ..., seq, ..., seqn)
# seqn: [s1, ..., sn]
# index: (i, ..., j), both i and j are integer
# keys: (key1, ..., key, ..., keyn)
# key: (seq[index_i], ..., seq[index_j])
# item_id: dictionary, {key1: [seq1, ..., seqm], ..., keyn: [seq3, ..., seqn]}
def group_by_key(sequence, index):
	item_id = {}
	for item in sequence:
		key = []
		for i in index:
			key.append(item[i])
		key = tuple(key)
		if key not in item_id:
			item_id[key] = [item]
		else:
			item_id[key].append(item)
	return item_id


# Linearly rescale data values having observed oldMin and oldMax into a new arbitrary range newMin to newMax
# newValue = (newMax - newMin)/(oldMax - oldMin) * (value - oldMin) + newMin
# or
# newValue = a * value + b, when a = (newMax - newMin)/(oldMax - oldMin) and b = newMin - a * oldMin	
def rescale(data, newMin, newMax):
	oldMin, oldMax = min(data), max(data)
	n = len(data)
	# all values are equal
	if not (oldMin < oldMax):
		return (newMin,) * n
	a = (newMax - newMin)/(oldMax - oldMin)
	b = newMin - a * oldMin
	newData = []
	for value in data:
		newData.append(a * value + b)	
	return newData


# return True if there is any digit in string
def hasNumbers(string):
	return any(char.isdigit() for char in string)

# return True if there is a pair of round brackets in string
# Note: return False if ') precede ('.
def hasBrackets(s):
	start = s.find('(')
	if start == -1:
		return False
	start += 1
	return ')' in s[start:]

# return string within '(' and ')'
def extract(s):
	start = s.find('(')
	if start == -1:
		return ''
	start += 1
	end = s.find(')', start)
	if end == -1:
		return ''
	else:
		return s[start:end]

# Parse alignment file output by water in EMBOSS and return the alignment
# waterFile: character string, output returned by water
#
# collumns: seq1, seq2:
#	1->6 	sequence number of the first residue in alignment
#	7	space character
#	8->-1	aligned sequence with gaps
# collumns: marker
#	1->6	space characters
#	7	space character
#	8->-1	character string composed of '|' (match), '.' (mismatch) and ' ' (gap)
def getAlignByWater(waterFile):
	align = []
	# if waterFile is file, then open it to readin
	# elif waterFile is not a file, then just process it as regular character string.
	if os.path.isfile(waterFile):
		fp = open(waterFile)
		waterFile = fp.read()
		fp.close()
	else:
		waterFile = waterFile.decode()
	scoreMarker = '# Score:'
	indexScore = waterFile.find(scoreMarker) + len(scoreMarker)
	score = float(waterFile[indexScore: waterFile.find('\n', indexScore)])
	#print(waterFile)
	#print('hello', waterFile[indexScore: waterFile.find('\n', indexScore)])
	lines = [line for line in waterFile.split('\n') if line.strip() != '' and line[0] != '#']
	if len(lines) > 0:
		start1, end1 = int(lines[0][14:20]), int(lines[-3][-6:])
		end2, start2 = int(lines[2][14:20]), int(lines[-1][-6:])
		#print('hello1', lines[0][14:20], lines[-3][-6:], lines[2][14:20], lines[-1][-6:])
		seq1 = marker = seq2 = ''
		for tlines in zip(*[iter(lines)]*3):
			seq1 += tlines[0][21:-7]
			marker += tlines[1][21:]
			seq2 += tlines[2][21:-7]
		align = [score, (seq1, seq2, marker), (start1, end1, start2, end2)]

	return align


# Given a DNA sequence strand, Return its complementary sequence strand
# mode: 
# '1', seq is composed of only uppercase letters; 
# '2', seq is composed of only lowercase letters; 
# '3', seq can be the mixture of lower and upper case letters.
def complementDNA(seq, mode):
	if mode == '1':
		return seq.translate(str.maketrans(constants.na1u, constants.na2u))
	elif mode == '2':
		return seq.translate(str.maketrans(constants.na1l, constants.na2l))
	elif mode == '3':
		return seq.translate(str.maketrans(constants.na1ul, constants.na2ul))
	else:
		print('Error: incorrect mode in complementDNA', file=sys.stderr)
		exit(0)
	
# Return a cleaned DNA sequence copy where non-standard bases are replaced by 'N'.
def cleanDNA(seq):
	bases = []
	#stdBases = 'ATCGRYN'
	stdBases = 'ATCG'
	for base in seq:
		if base.upper() in stdBases:
			bases.append(base)
		else:
			bases.append('N')
	return ''.join(bases)


# convert cigar sting returned by SSW python interface into the pairs with one number and one character
# for example:
# Return [(4, 'M'), (2, 'I'), (8, 'M'), (1, 'D'), (10, 'M'), (6, 'S')] if cigarString == '4M2I8M1D10M6S'
def parseCigarString(cigarString):
	return [(int(pair[:-1]), pair[-1]) for pair in re.findall(r'\d+[MIDNSHP=X]', cigarString)]


# print a alignment into a string
# seq1, seq2: sequences of two aligned DNA strand segment
# cigarStr: '4M2I8M1D10M6S'
# cigarPair: [(4, 'M'), (2, 'I'), (8, 'M'), (1, 'D'), (10, 'M'), (6, 'S')]
# 	Note: for details of cigar sting, Please check the document "The SAM Format Specification", 
# 	http://samtools.github.io/hts-specs/SAMv1.pdf, particularly the "An example" section and 
#	"CIGAR: CIGAR string" section.
# 
def buildAlignment(sequence1, sequence2, align, cigarStr):
	begin1, end1, begin2, end2 = align.ref_begin+1, align.ref_end+1, align.query_begin+1, align.query_end+1
	header = {	'conflict': False,
			'score': align.score,
			'begin1': begin1,
			'end1': end1,
			'begin2': begin2,
			'end2': end2}
	seq1name = 'seq1'
	seq2name = 'seq2'
	line1 = '{:<10} '.format(seq1name)
	line2 = '{:<10} '.format(' ')
	line3 = '{:<10} '.format(seq2name)
	index1, index2 = 0, 0
	line1 += '{:>8} '.format(begin1 + index1)
	line3 += '{:>8} '.format(begin2 + index2)
	line2 += '{:>8} '.format(' ')
	
	# build alignment from the range defined by 
	# sequence1[align.ref_begin: align.ref_end+1] and sequence2[align.query_begin: align.query_end+1]
	seq1 = sequence1[align.ref_begin: align.ref_end+1]
	seq2 = sequence2[align.query_begin: align.query_end+1]
	cigarPair = parseCigarString(cigarStr)
	for pair in cigarPair:
		if pair[1] == 'I':
			line1 += '-' * pair[0]
			line3 += seq2[index2: index2+pair[0]]
			line2 += ' ' * pair[0]
			index2 += pair[0]
		elif pair[1] == 'D':
			line1 += seq1[index1: index1+pair[0]]
			line3 += '-' * pair[0]
			line2 += ' ' * pair[0]
			index1 += pair[0]
		elif pair[1] == 'M':
			s1 = seq1[index1: index1+pair[0]]
			s2 = seq2[index2: index2+pair[0]]
			line1 += s1
			line3 += s2
			for c1, c2 in zip(s1, s2):
				if c1 == c2:
					line2 += '|'
				else:
					line2 += '*'
			index1 += pair[0]
			index2 += pair[0]
		elif pair[1] != 'S':
			e = cigarStr + ' submittedSeq1:' + sequence1 + ' submittedSeq2:' + sequence2
			raise RuntimeError(e)

		#print(pair)
		#print(line1)
		#print(line2)
		#print(line3)
	line1 += ' {:<8}'.format(align.ref_begin + index1)
	line3 += ' {:<8}'.format(align.query_begin + index2)
	if align.ref_begin + index1 != align.ref_end+1 or align.query_begin + index2 != align.query_end+1:
		header['conflict'] = True
		header['end1'] = align.ref_begin + index1
		header['end2'] = align.query_begin + index2
		'''
		fmt = 'Warning: alignment path conflicts with aligned range reported by SSW!\n'
		fmt += ' seq1Begin:{} seq1Alinged:{} seq1End:{} seq2Begin:{} seq2Alinged:{} seq2End:{}\n'
		fmt += ' cigarStr:{}\n alignedSeq1:{}\n alignedSeq2:{}\n submittedSeq1:{}\n submittedSeq2:{}'
		w = fmt.format(
				align.ref_begin+1, index1, align.ref_end+1, 
				align.query_begin+1, index2, align.query_end+1,
				cigarStr, 
				seq1, seq2, 
				sequence1, sequence2)
		print(w)
		'''
	return (header, line1, line2, line3)

# Build the match line between two aligned sequences.
# Return matchLine
# matchLine: character string, for example, '|*|*  **|'
# seq1: character string, for example, 'TAGG--AGC'
# seq2: character string, for example, 'TGGACGGCC'
def buildMatchLine(seq1, seq2):
	matchLine = ''
	for c1, c2 in zip(seq1, seq2):
		if c1 == c2:
			matchLine += '|'
		elif c1 == '-' or c2 == '-':
			matchLine += ' '
		else:
			matchLine += '*'
	return matchLine

	
# Shorten ir to the reasonable length
# Rules to shorten:
# Discard the aligned blocks when irId/irLen < 0.7 with irLen < constants.stringenShortestIR 
#							or irLen > constants.stringentLongestIR
# Discard the aligned blocks when irId/irLen < 0.6 with constants.stringentShortestIR <= irLen <= constants.stringentLongestIR
#
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
#
def shortenIR(ir):
	if len(ir) == 0:
		return ir
	line = buildMatchLine(ir[8], ir[9])
	# for example line = '||||||*| ||*||| || ||| |   | |||| *| |*||*|||||||||||||||'
	# blocks = [list(g) for k,g in itertools.groupby(line)]
	# blocks = [	['|', '|', '|', '|', '|', '|'], 
	#		['*'], 
	#		['|'], 
	#		[' '], 
	#		['|', '|'], 
	#		['*'], 
	#		['|', '|', '|'], 
	#		[' '], 
	#		['|', '|'], 
	#		[' '], 
	#		['|', '|', '|'], 
	#		[' '], 
	#		['|'], 
	#		[' ', ' ', ' '], 
	#		['|'], 
	#		[' '], 
	#		['|', '|', '|', '|'], 
	#		[' '], 
	#		['*'], 
	#		['|'], 
	#		[' '], 
	#		['|'], 
	#		['*'], 
	#		['|', '|'], 
	#		['*'], 
	#		['|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|', '|']]
	# len(blocks) = 26
	#ids = []
	irId = 0
	irLen = 0
	for k,g in itertools.groupby(line):
		g = list(g)
		blockLen = len(g)
		irLen += blockLen
		if g[0] == '|':
			irId += blockLen
			#ids.append((irId, irLen))
			#
			# break loop, namely, stop growing alignment
			if irId/irLen < constants.optIrIdentity:
				break
			elif (irLen < constants.stringentShortestIR or irLen > constants.stringentLongestIR) and irId/irLen < constants.stringentIrIdentity:
				break
	# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
	# Build shorten ir
	score = irId + irLen
	nGaps = line[:irLen].count(' ')
	seq1 = ir[8][:irLen]
	seq2 = ir[9][:irLen]
	gap1 = seq1.count('-')
	gap2 = seq2.count('-')
	end1 = ir[4] + irLen - gap1 - 1
	end2 = ir[6] + irLen - gap2 - 1
	return [score, irId, irLen, nGaps, ir[4], end1, ir[6], end2, seq1, seq2]

# Filter ir by applying more stringent criteria:
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
# optIRsim, stringentIRsim: float, constants.optIrIdentity and constants.stringentIrIdentity
#
# 1) irId/irLen must be greater than stringentIrIdentity (0.7) if irLen < 5 or irLen > 55, else irId/irLen > optIrIdentity (0.6)
#
def filterIRbyCutoff(ir, optIRsim, stringentIRsim):
	if len(ir) == 0:
		return [] 

	sim = ir[1]/ir[2]
	if sim < optIRsim:
		ir = []
	elif (ir[2] < 5 or ir[2] > 55) and sim < stringentIRsim:
			ir = []
	return ir



# Find the number of matches in core regions which is composed of only consecutive matches.
# seq1: character string, for example, 'TAGGG--AGGC'
# seq2: character string, for example, 'TGGGACGGCGC'
# irIdCore: integer, number of matches in core regions
def getIrIdCore(seq1, seq2):
	matchLine = buildMatchLine(seq1, seq2)
	# search consecutive matches in matchLine
	irIdCore = 0

	# Option1:
	#	1) core region composed of at least two consecutive matched bases 
	#for region in re.findall(r'\|{2,}', matchLine):
	#
	# Option2:
	#	1) core region, composed of at least three consecutive matched bases 
	for region in re.findall(r'\|{3,}', matchLine):
	#
		irIdCore += len(region)
	
	return irIdCore

# Return an empirical score for an ir. The greater the score is, the better the ir is.
# ir: [] or [score, irId, irLen, nGaps, start1, end1, start2, end2, seq1, seq2]
#
def irScore(ir):
	if len(ir) == 0:
		# set a very negative value as the score of no TIR
		score = -9999.9
	else:
		irIdCore = getIrIdCore(ir[-2], ir[-1])

		# mismatch = irLen - nGaps - irId
		# irIdNonCore = irId - irIdCore

		# Option0:
		#	score	= 3(irIdCore - nGaps) + irIdNonCore - mismatch
		#		= 2(irIdCore - nGaps) + irId - (nGaps + mismatch)
		#		= 2(irIdCore - nGaps) + 2*irId - iLen
		#		= 2(irIdCore + irId - nGaps) - irLen
		score = 2 * (irIdCore + ir[1] - ir[3]) - ir[2] 

		# score	= irIdCore - nGaps + x

		# Option1: x = irIdNonCore - mismatch
		# 	score	= irIdCore - nGaps + x
		#		= irIdCore - nGaps + (irIdNonCore - mismatch) 
		#		= irIdCore + irIdNonCore - (nGaps + mismatch) 
		#		= irId - (nGaps + mismatch)
		#		= 2*irId - irLen
		#		= irId - (nGaps + mismatch)
		# Here, score == 0 means irId/irLen == 50%
		#score = 2 * ir[1] - ir[2]

		#
		# Option2: x = (irIdNonCore - mismatch) / (irIdNonCore + mismatch + 1)
		#	Here, irIdNonCore == mismatch means x contributes zero to final score.
		#
		# 	Pay special attention to the following cases: 
		#		irIdCore == nGaps and/or irIdNonCore == mismatch == 0
		#
		#	score	= irIdCore - nGaps + x
		#	socre	= irIdCore - nGaps + 
		#			(irId - irIdCore - (irLen-nGaps-irId)) / (irLen - irIdCore - nGaps + 1)
		#		= irIdCore - nGaps + 
		#			(2*irId - irIdCore - irLen + nGaps) 
		#			/ (irLen - irIdCore - nGaps + 1)
		#		= irIdCore - nGaps + 
		#			(2*irId - irLen - (irIdCore - nGaps)) 
		#			/ (irLen -2*nGaps - (irIdCore - nGaps) + 1)
		#
		#	set y = irIdCore - nGaps, then
		#	score	= y + (2*irId - irLen - y) / (irLen - 2*nGaps - y + 1)
		#y = irIdCore - ir[3]
		#score = y + (2*ir[1] - ir[2] - y) / (ir[2] - 2*ir[3] - y + 1)

		#
		# Option3: x = - (irLen - match) /irLen 
		#	score	= irIdCore -nGaps + x
		#		= irIdCore -nGaps - (irLen - irId)/irLen
		#		= irIdCore -nGaps - (1 - irId/irLen)
		#		= irIdCore + irId/irLen - nGaps - 1
		#score = irIdCore + ir[1]/ir[2] - ir[3] - 1

		#
		# Option4: x = (irIdNonCore - mismatch)/2
		# 	score	= 2*(irIdCore - nGaps) + irIdNonCore - mismatch
		#		= 3*irIdCore - 2*nGaps - mismatch
		#		= 3*irIdCore - nGaps - (irLen - irId)
		#score = 3*irIdCore - ir[3] - ir[2] + ir[1]


	return score

# Build name of matrix file holding match and mismatch values
# Matrix file example: 
# EDNAFULL.2.6.IR.water, EDNAFULL.3.4.IR.water, EDNAFULL.3.6.IR.water, EDNAFULL.5.4.IR.water
def convert2matrixFile(match, mismatch, dir2matrix):
	return os.path.join(dir2matrix, 'EDNAFULL.{}.{}.IR.water'.format(match, - mismatch))

# File name, EDNAFULL.2.6.IR.water, means match and mismatch are 2 and -6 respectively.
# example: 
# matrixFile: EDNAFULL.2.6.IR.water
# Return: (2, -6)
def resolveMatrixFileName(matrixFile):
	fileName = os.path.basename(matrixFile)
	if fileName == '':
		print('Error: matrix file required', matrixFile, file=sys.stderr)
		exit(0)
	else:
		digits = fileName.split('.',3)[1:3]
	return (int(digits[0]), - int(digits[1]))


# Convert (gapopen, gapextend, matrixfile) to (gapopen, gapextend, match, mismatch)
# filters4water:	[filter4water, ..., filter4water]
# filter4water:		(gapopen, gapextend, matrixfile)
# Return:		filters
# filters:		[filter, ..., filter]
# filter:		(gapopen, gapextend, match, mismatch)
def convertFilters4water(filters4water):
	filters = []
	for filter in filters4water:
		match, mismatch = resolveMatrixFileName(filter[2])
		# convert mismatch from a negative value to a positive value
		filters.append((filter[0], filter[1], match, -mismatch))
	return filters

# Convert (gapopen, gapextend, match, mismatch) to (gapopen, gapextend, matrixfile)
# filters:		[filter, ..., filter]
# filter:		(gapopen, gapextend, match, mismatch), here mismatch > 0
# dir4embossdata:	directory holding the matrix files used by water program
# filters4water:	[filter4water, ..., filter4water]
# filter4water:		(gapopen, gapextend, matrixfile)
def convertFilters2water(filters, dir2matrix):
	filters4water = []
	for filter in filters:
		matrixFile = convert2matrixFile(filter[2], -filter[3], dir2matrix)
		filters4water.append((filter[0], filter[1], matrixFile))
	return filters4water

# Return filters shared by filters4ssw and filtersFromWater
def commonFilters(filters4ssw, filtersFromWater):
	filters = []
	for filter1 in filters4ssw:
		for filter2 in filtersFromWater:
			if filter1 == filter2:
				filters.append(filter1)
	return filters[-1:]

# A single measure that trades off precision versus recall is the F measure, 
# which is the weighted harmonic mean of precision and recall:
# Fbeta = (1 + beta**2) * recall * precision / (recall + beta**2 * precision)
# Fbeta attaches beta times as much importance to recall as precision.
# Reference:
# Christopher D. Manning, Prabhakar Raghavan and Hinrich Schütze, 
# Introduction to Information Retrieval, Cambridge University Press. 2008
#
# In our practice, we use beta = 1 or 2, 
# F1 = 2 * recall * precision / (recall + precision)
# F2 = 5 * recall * precision / (recall + 4*precision)
# Here, precision = 1 - fdr, recall = sensitivity, fdr(false discovery rate) = false postive hit rate = nfp/nhits
def fmeasure(recall, precision, beta):
	denominator = recall + beta**2 * precision
	if denominator > 0:
		return (1 + beta**2) * recall * precision / denominator 
	else:
		return 0.0

 
# Check if hit piece is overlapped with isfinder annotation
# 0: hit and isfinder are located on different strands
# overlap/min(b-a, d-c): overlap > 0 if it is overlap, else overlap <= 0
# Note: IS element does not depends on strand while gene depends on strand.
# 	So, we don't need to compare strand when examing the overlap between 
#	two IS elements.
#
# The is_overlap() is retained as the alternative of is_overlap_min() or is_overlap_max() 
# for the code compatibility.
def is_overlap(hit_strand, hit_begin, hit_end, strand, begin, end):
	#if hit_strand != strand:
	#	return 0.0

	# Simply ignore strand
	a, b, c, d = hit_begin, hit_end, begin, end
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	if overlap > 0:
		#return float(overlap) / (min(b-a, d-c) + 1)
		return float(overlap) / (max(b-a, d-c) + 1)
		#return float(overlap) / (b-a + 1)
	else:
		return 0.0

def is_overlap_min(hit_strand, hit_begin, hit_end, strand, begin, end):
	# Simply ignore strand
	a, b, c, d = hit_begin, hit_end, begin, end
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	if overlap > 0:
		return float(overlap) / (min(b-a, d-c) + 1)
	else:
		return 0.0

def is_overlap_max(hit_strand, hit_begin, hit_end, strand, begin, end):
	# Simply ignore strand
	a, b, c, d = hit_begin, hit_end, begin, end
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	if overlap > 0:
		return float(overlap) / (max(b-a, d-c) + 1)
	else:
		return 0.0

# Check if two ORFs (genes) are overlapped
# 0: hit and isfinder are located on different strands
# overlap/min(b-a, d-c): overlap > 0 if it is overlap, else overlap <= 0
#
def orf_overlap(orf1, orf2):
	hit_strand, hit_begin, hit_end = orf1
	strand, begin, end = orf2
	if hit_strand != strand:
		return 0.0
	a, b, c, d = hit_begin, hit_end, begin, end
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	if overlap > 0:
		return float(overlap) / (min(b-a, d-c) + 1)
	else:
		return 0.0

# Check if piece1 is overlapped with piece2.
# p1: (a, b), a <= b, float or int
# p2: (c, d), c <= d, float or int
#
def overlap(p1, p2):
	a, b = p1
	c, d = p2
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	if overlap > 0:
		#return float(overlap) / (min(b-a, d-c) + 1)
		return float(overlap) / (max(b-a, d-c) + 1)
		#return float(overlap) / (b-a + 1)
	else:
		return 0.0

def overlap_min(p1, p2):
	a, b = p1
	c, d = p2
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	if overlap > 0:
		return float(overlap) / (min(b-a, d-c) + 1)
	else:
		return 0.0

def intersection(p1, p2):
	a, b = p1
	c, d = p2
	# a <=b and c <= d are satisfied
	overlap = min(b, d) - max(a, c) + 1
	return overlap

def intergap(p1, p2):
	a, b = p1
	c, d = p2
	# a <=b and c <= d are satisfied
	gap = max(a, c) - min(b, d) - 1
	return gap

# makeblastdb -dbtype nucl -in output4FragGeneScan1.19_illumina_5/NC_002754.1.fna.ffn -out blastdb/NC_002754.1.fna.ffn
def seq2blastdb(seqFile, db):
	cmd = constants.makeblastdb
	cmdline = [cmd, '-dbtype nucl', '-in', seqFile, '-out', db]
	do_cmd = shlex.split(' '.join(cmdline))
	subprocess.check_call(do_cmd, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL)

# Search all IS elements (Tpase) ORFs against IS element (Tpase ORF) database.
# command: blastn -query /home/data/insertion_sequence/output4FragGeneScan1.19_illumina_5/NC_002754.1.fna.ffn \
#		-db /home/data/insertion_sequence/blastdb/NC_002754.1.fna.ffn -out blast.6 \
#		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
#		nident qlen slen' -perc_identity 90 -dust no
# Note1: There is 0 hit when querying 1471545_1471675_- against NC_011891.1.orf.fna because blastn by default 
#	to enable filtering query sequence with DUST to mask low complex repeat sequence in query. 
#	So we must disable it with '-dust no'.
# Note2: blastn-short is BLASTN program optimized for sequences shorter than 50 bases. We will use it in blastn search 
#	when dealing with tir sequence as tir is usually shorter than 55 bases.  
def doBlastn(query, db, out, strand='both', task='megablast', perc_ident=100):
	blast = constants.blastn
	outfmt = shlex.quote('6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen')
	#perc_identity = str(constants.sim4iso * 100)
	perc_identity = str(perc_ident)
	if task == 'blastn-short':
		wordsize = '7'
		#wordsize = '5' # default value for blastn-short is 7 but we use smaller value because the length of tir is usually among 5<= length <=55.
		#wordsize = '4' # wordsize must >= 4
	elif task == 'megablast':
		wordsize = '28' # default value for megablast
	else:
		wordsize = '11' # default value for blastn
	cmd = [blast, '-query', query, '-db', db, '-out', out, '-outfmt', outfmt, '-perc_identity', perc_identity, 
			'-strand', strand, '-dust', 'no', '-task', task, '-word_size', wordsize]
	do_cmd = shlex.split(' '.join(cmd))
	if subprocess.call(do_cmd, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL) != 0:
		e = 'Fail to run {}'.format(cmd)
		raise RuntimeError(e)

def blastnSearch(query, db, out, strand='both', task='megablast'):
	blast = constants.blastn
	outfmt = shlex.quote('6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen')
	#perc_identity = str(constants.sim4iso * 100)
	perc_identity = str(constants.SIM4ISO)
	if task == 'blastn-short':
		wordsize = '7'
		#wordsize = '5' # default value for blastn-short is 7 but we use smaller value because the length of tir is usually among 5<= length <=55.
		#wordsize = '4' # wordsize must >= 4
	elif task == 'megablast':
		wordsize = '28' # default value for megablast
	else:
		wordsize = '11' # default value for blastn
	cmd = [blast, '-query', query, '-db', db, '-out', out, '-outfmt', outfmt, '-perc_identity', perc_identity, '-strand', strand, '-dust', 'no', '-task', task, '-word_size', wordsize]
	do_cmd = shlex.split(' '.join(cmd))
	if subprocess.call(do_cmd, shell=False, universal_newlines=False, stdout=subprocess.DEVNULL) != 0:
		e = 'Fail to run {}'.format(cmd)
		raise RuntimeError(e)

# Search all IS elements (Tpase) ORFs against IS element (Tpase ORF) database.
# Return out
#
# command: blastn -db db -perc_identity 90 -strand strand -dust no -task task -word_size wordsize -num_threads nthreads \
#		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen'
# Note1: query and outfile are stdin and stdou by default, respectively
# Note2: blastn-short is BLASTN program optimized for sequences shorter than 50 bases. We will use it in blastn search 
#	when dealing with tir sequence as tir is usually shorter than 55 bases.  
#
def doBlastnOnStream(query, db, strand='both', task='megablast', perc_ident=100, nthreads=1):
	blast = constants.blastn
	outfmt = shlex.quote('6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen')
	perc_identity = str(perc_ident)
	num_threads = str(nthreads)
	if task == 'blastn-short':
		wordsize = '7'
		#wordsize = '5' # default value for blastn-short is 7 but we use smaller value because the length of tir is usually among 5<= length <=55.
		#wordsize = '4' # wordsize must >= 4
	elif task == 'megablast':
		wordsize = '28' # default value for megablast
	else:
		wordsize = '11' # default value for blastn
	cmd = [blast, 
		'-db', db, '-perc_identity', perc_identity, '-strand', strand, '-dust', 'no', 
		'-task', task, '-word_size', wordsize, '-num_threads', num_threads,
		'-outfmt', outfmt
		]
	do_cmd = shlex.split(' '.join(cmd))
	blastn = subprocess.Popen(do_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			universal_newlines=True)
	out, err = blastn.communicate(input=query)
	return (out, err)

# Search protein sequence against protein database.
# command: blastp -db db -evalue 1e-10 -task task -num_threads nthreads \
#		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen'
# Note1: query and outfile are stdin and stdou by default, respectively
def doBlastpOnStream(query, db, task='blastp', e_value=1e-10, nthreads=1):
	blast = constants.blastp
	outfmt = shlex.quote('6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen')
	evalue = str(e_value)
	num_threads = str(nthreads)
	cmd = [blast,
		'-db', db, '-evalue', evalue, '-task', task, '-num_threads', num_threads,
		'-outfmt', outfmt
		]
	do_cmd = shlex.split(' '.join(cmd))
	blastn = subprocess.Popen(do_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			universal_newlines=True)
	out, err = blastn.communicate(input=query)
	return (out, err)

# blastn -query query -subject subject ....
def doBlastn2seqOnStream(query, subject, strand='both', task='megablast', perc_ident=100):
	blast = constants.blastn
	outfmt = shlex.quote('6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen')
	perc_identity = str(perc_ident)
	if task == 'blastn-short':
		wordsize = '7'
		#wordsize = '4' # wordsize must >= 4
	elif task == 'megablast':
		wordsize = '28' # default value for megablast
	else:
		wordsize = '11' # default value for blastn
	cmd = [blast, 
		'-subject', subject, '-perc_identity', perc_identity, '-strand', strand, '-dust', 'no', 
		'-task', task, '-word_size', wordsize, '-outfmt', outfmt
		]
	do_cmd = shlex.split(' '.join(cmd))
	blastn = subprocess.Popen(do_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
			universal_newlines=True)
	out, err = blastn.communicate(input=query)
	return (out, err)

# Get IS element copy number from the file output by blast search
# The output is produed by blastn with options:
# -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send nident qlen slen' \
#		-perc_identity 90
# For example:
# command: blastn -query /home/data/insertion_sequence/output4FragGeneScan1.19_illumina_5/NC_002754.1.fna.ffn \
#		-db /home/data/insertion_sequence/blastdb/NC_002754.1.fna.ffn -out blast.6 \
#		-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
#		nident qlen slen' -perc_identity 90 -dust no
# Read blast output and return hitpairs
# hitpairs: [hitpair, ..., hitpair]
# hitpair: {'qseqid': qseqid, 'sseqid': sseqid, \
#		'pident': pident, 'length': length, 'nident': nident, 'qlen': qlen, 'slen': slen}
#
# Input: 
# file: blast output file
# min4coverage: coverage cutoff = alignedLength/min(qlen,slen)
# seqtype: 'regular' for non-short sequence, 'short' for short sequence like tir in IS element
def getBlastout(file, min4coverage):
	fp = open(file, 'r')
	hits = []
	for line in fp:
		words = line.split()
		pident = float(words[2])
		# retain only hit with exact silva full-length gene sequence
		if pident < 100:
			continue

		length = int(words[3])
		nident = int(words[12])
		qlen = int(words[13])
		slen = int(words[14])

		# retain only long gene sequence
		if length < 1200:
			continue

		# retain only hit with exact silva full-length gene sequence
		#if nident != slen:
		#	continue

		# check the coverage of aligned length vs min(querylen, subjectlen),
		#if length < int(min(qlen,slen) * min4coverage):
		#	continue

		hit = {}
		hit['qseqid'] = words[0]
		hit['sseqid'] = words[1]
		hit['pident'] = float(words[2])
		hit['length'] = length
		hit['nident'] = nident

		# coordinates of alignment in query sequence and subject sequence
		hit['qstart'] = int(words[6])
		hit['qend'] = int(words[7])
		hit['sstart'] = int(words[8])
		hit['send'] = int(words[9])

		hit['qlen'] = qlen
		hit['slen'] = slen

		hits.append(hit)

	return hits

def getBlastResult(file, min4coverage):
	fp = open(file, 'r')
	hits = []
	for line in fp:
		words = line.split()
		qlen = int(words[13])
		slen = int(words[14])
		pident = float(words[2])
		length = int(words[3])
		# check the coverage of aligned length vs min(querylen, subjectlen)
		if length < int(min(qlen,slen) * min4coverage):
			continue
		hit = {}

		# words[0] example: 'gi|15896971|ref|NC_002754.1|_1_4458_+', 'NC_002754.1_1_4458_+', 
		#	'1_4458|15_4443_+', '15_4443_+', '0', '11'
		hit['qseqid'] = words[0]

		# words[1] example: 'gi|15896971|ref|NC_002754.1|_1_4458_+' or 'NC_002754.1_1_4458_'
		#	'1_4458|15_4443_+', '15_4443_+', '0', '11'
		hit['sseqid'] = words[1]

		hit['pident'] = float(words[2])
		hit['length'] = length
		hit['nident'] = int(words[12])
		hit['qlen'] = qlen
		hit['slen'] = slen
		hits.append(hit)

	return hits

# get results produced by blastp
def getBlastpResultOnStream(filec):
	fp = io.StringIO(filec)
	hits = []
	ids4query = set()
	for line in fp:
		line = line.strip()
		words = line.split()
		# -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident qlen slen'
		qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,nident,qlen,slen = words
		# Only get the first hit for each query sequence as the hits related the same query are ordered by evalue and we just need the IS family name from the subject hit
		# with the smallest evalue.
		if qseqid in ids4query:
			continue
		else:
			ids4query.add(qseqid)
		hit = {}
		hit['qseqid'] = qseqid
		hit['sseqid'] = sseqid
		hit['length'] = length
		hit['qstart'] = qstart
		hit['qend'] = qend
		hit['sstart'] = int(sstart)
		hit['send'] = int(send)
		hit['evalue'] = float(evalue)
		hit['nident'] = int(nident)
		hit['qlen'] = qlen
		hit['slen'] = int(slen)
		hit['pident'] = float(pident)
		hits.append(hit)
	return hits

def getBlastResult4dnaOnStream(filec):
	fp = io.StringIO(filec)
	hits = []
	for line in fp:
		words = line.split()

		# check the coverage of aligned length vs min(querylen, subjectlen)
		# ignore the alignment with aligned length < isMin (usually isMin = 400, namely, shortest IS element)
		length = int(words[3])
		#if length < constants.isMin:
		#	continue

		#items = words[0].rsplit('_', maxsplit=6)
		#minLen4is = int(items[1]) # minimal length of IS in the specific family which ORF belongs to

		# word[0]: seqid_family_cluster_seqbegin_seqend_orfBegin_orfEnd_orfStrand
		items = words[0].rsplit('_', maxsplit=7)
		family = items[1]
		# cluster = items[2]
		minLen4is = constants.minMaxLen4is[family][0]
		if length < minLen4is:
			continue

		seqbegin = int(items[3])
		seqend = int(items[4])
		orfBegin = int(items[5])
		orfEnd = int(items[6])
		orfStrand = items[7]

		qlen = int(words[13])
		if seqend-seqbegin+1 != qlen:
			e = 'Error: different lengths for query sequence {}: {}({}-{}+1) {}'.format(
				words[0], seqend-seqbegin+1, seqend, seqbegin, qlen)
			RuntimeError(e)
		start = int(words[6])
		end = int(words[7])
		move = start - 1
		if orfStrand == '+':
			qstart = seqbegin + move
			qend = qstart + length - 1
		else:
			qend = seqend - move
			qstart = qend - length + 1

		orfLen = orfEnd - orfBegin + 1
		intersect = intersection((orfBegin, orfEnd), (qstart, qend))
		if intersect < 1:
			continue
		copy = False
		#minLen4orf4tpase = constants.minMax4tpase[family][0]
		#maxLen4orf4tpase = constants.minMax4tpase[family][1]
		minLen4orf4pep = constants.minMax4tpase[family][2]
		#
		# Note: it means gene prediction may be not accuracy when predicted ORF is longer 
		#	than the longest tpase ORF in IS family.
		#if intersect < orfLen*constants.minOverlap4orf2aligned:
		# 
		# alignment overlapped with whole orf or part of orf
		#if orfLen <= maxLen4orf and intersect >= orfLen*constants.minOverlap4orf2aligned:
		#if orfLen <= maxLen4orf and intersect >= orfLen*0.9:
		#
		# alignment covers the full orf, namely, length >= orfLen and intersect == orfLen
		#if orfLen <= maxLen4orf and qstart <= orfBegin < orfEnd <= qend:
		#if orfLen <= maxLen4orf and intersect == orfLen:
		#if orfLen <= maxLen4orf and intersect == orfLen:
		#	copy = True
		# alignment overlapped with part of orf
		#if orfLen > maxLen4orf and intersect >= orfLen*constants.minOverlap4orf2aligned4longORF:
		#if orfLen > maxLen4orf and intersect >= orfLen*0.5:
		#if orfLen > maxLen4orf and intersect >= minLen4orf:

		if orfLen >= minLen4orf4pep:
			if intersect >= minLen4orf4pep:
				copy = True
		else:
			if intersect >= orfLen:
				copy = True


		if copy == False:
			#print('intersection({}) between {} and {} in alignment is less than threshold: \n{}'.format(
			#	intersect, (orfBegin, orfEnd, orfStrand), (qstart, qend), line))
			continue
		#else:
			#print('hello copy', intersect, (orfBegin, orfEnd, orfStrand), (qstart, qend))

		hit = {}

		# words[0] example: 'gi|556503834|ref|NC_000913.3|_IS200/IS605_12_1502643_1505379_1503157_1504865_+'
		hit['qseqid'] = words[0]
		hit['orfBegin'] = orfBegin
		hit['orfEnd'] = orfEnd
		# words[1] example: 'gi|15896971|ref|NC_002754.1|
		hit['sseqid'] = words[1]
		hit['length'] = length
		hit['qstart'] = qstart
		hit['qend'] = qend
		hit['sstart'] = int(words[8])
		hit['send'] = int(words[9])
		hit['nident'] = int(words[12])
		hit['qlen'] = qlen
		hit['slen'] = int(words[14])
		hit['pident'] = float(words[2])

		hits.append(hit)
	return hits

# create dir if it does not exist yet
def makedir(dir):
	if not os.path.exists(dir):
		try:
			os.makedirs(dir)
		except OSError as error:
			if error.errno != errno.EEXIST:
				raise
	

# Get the genome coordinates of each gene from proteome file
# Return genes
#
# genes: [(strand, begin, end), ..., (strand, begin, end)]
# strand: character, '+' or '-'
# begin, end: int, genome coordinates
def get_proteome(proteome_file):
	genes = []
	fp = open(proteome_file, 'r')
	for line in fp:
		if line[0] != '>':
			continue
		location = line.strip().rsplit('_', 3)[-3:]
		genes.append((location[2], int(location[0]),int(location[1])))
	fp.close()
	return genes


# cdss: [(seqid, id,seq), ...]
# seqid: sequence identifier, e.g. SRS075404_LANL_scaffold_1, C3691328
# id: cds identifier, e.g. SRS075404_LANL_scaffold_1_1_414_+, C3691328_7626_8378_-
def getcds(file):
	cdss = []
	seq = []
	with open(file, 'r') as fp:
		for line in fp:
			line = line.strip()
			if line == '':
				continue # remove blank line
			if line[0] == '>':
				if len(seq) > 0:
					cdss.append((seqid, id,''.join(seq)))
				id = line[1:].split(maxsplit=1)[0]
				seqid = id.rsplit('_',maxsplit=3)[0]
				seq = []
				continue
			seq.append(line)
	cdss.append((seqid, id,''.join(seq)))
	return cdss


# read genebank .fna file and return the identifier of sequence
def rdGbFna(file):
	fp = open(file, 'r')
	for line in fp:
		if line.replace('\b','').strip() == '':
			continue # remove blank line
		if line[0] == '>':
			seqid = line.strip().split(maxsplit=1)[0][1:]
			break
	fp.close()
	return seqid

# read genebank .faa file and return the peptide sequences
# prot: [p, ..., p]
# p: {'id': id, 'pid': pid, 'seq': seq, ...}
def rdGbFaa(file):
	prot = []
	seq = []
	fp = open(file, 'r')
	for line in fp:
		if line.replace('\b','').strip() == '':
			continue # remove blank line
		if line[0] == '>':
			if len(seq) > 0:
				p['seq'] = ''.join(seq)
				prot.append(p)
			p = {}
			p['id'] = line[1:].split(maxsplit=1)[0]
			p['pid'] = p['id'].split('|',2)[1]
			seq = []
			continue
		seq.append(line.strip())
	fp.close
	p['seq'] = ''.join(seq)
	prot.append(p)
	return prot

# read genebank .ptt file and return the locations and PIDs
# geneLoc: [loc, ..., loc]
# loc: {'loc': [begin, end, strand], 'pid': pid}
# begin, end: character string
# pid: character string
def rdGbPtt(file):
	geneLoc = []
	fp = open(file, 'r')
	# remove heading lines in .ptt file
	for line in fp:
		if line[:8] == 'Location':
			break
	# get data
	for line in fp:
		loc = {}
		items = line.split(maxsplit=4)
		locations = items[0].split('..', 1)
		
		loc['loc'] = [locations[0], locations[1], items[1]]
		loc['pid'] = items[3]
		geneLoc.append(loc)

	fp.close
	return geneLoc

# Convert GeneBank protein info (NC_000913.faa and NC_000913.ptt) 
# into FragGeneScan protein file format(NC_000913.fna.faa)
def gb2fgs4protein(gbFna, gbFaa, gbPtt, fgs):
	# make directories
	makedir(os.path.dirname(fgs))

	fp4fgs = open(fgs, 'w')
	# seqid: eg. 'gi|15644634|ref|NC_000915.1|'
	seqid = rdGbFna(gbFna)
	prot = rdGbFaa(gbFaa)
	geneLoc = rdGbPtt(gbPtt)
	for p in prot:
		for loc in geneLoc:
			if p['pid'] == loc['pid']:
				locStr = '_'.join(loc['loc'])
				fasta = fastaFormat('_'.join([seqid, locStr]), p['seq'])
				print(fasta, file=fp4fgs)
				break
	fp4fgs.close()

# gbk: {'accver':accver, 'gi':gi, 'prots':prots, 'seq':seq, ...}
# accver: accession.version, eg. 'NC_000915.1'
# gi: gi, eg. '15644634', character string
# prots: [prot, ..., prot]
#   prot: {'orf':orf, 'pep': pep}
#   orf: (start, end, strand), genome coordinates of cds in .gbk file, 
#	start and end are int, strand is character
#   pep: amino acid sequence, character string
# seq: nucleic acid sequence
def rdGbk(gbkFile):
	gbk = {}
	gbk['prots'] = []
	fp4gbkFile = open(gbkFile, 'r')
	for line in fp4gbkFile:
		if 'VERSION     ' == line[:12]:
			items = line[12:].split(maxsplit=2)
			gbk['accver'] = items[0]
			if len(items) > 1:
				gbk['gi'] = items[1]
			continue
		if '     CDS             ' == line[:21]:
			if line[21] == 'c' and line[32].isdigit():
				strand = '-'
				coord = line[32:].strip().strip(')').strip().split('..', maxsplit=1)
			elif line[21].isdigit():
				strand = '+'
				coord = line[21:].strip().split('..', maxsplit=1)
			# not normal CDS, for example:
			# '     CDS             join(58474..59052,59052..59279)'
			# '     CDS             complement(join(380844..381259,382591..382872))'
			else:
				continue
			while True:
				nextLine = next(fp4gbkFile)
				if '                     /translation=' == nextLine[:34]:
					nopep = False
					break
				# special case in NC_000913.gbk: 
				# CDS(272847..273954) has no translated protein sequence.
				elif len(nextLine[:21].strip()) > 0:
					nopep = True
					break
			if nopep == True:
				continue
			prot = {}
			start, end = [int(x) for x in coord]
			prot['orf'] = (start, end, strand)
			gbk['prots'].append(prot)
			pep = []
			if nextLine[36:].strip()[-1] != '"':
				pep.append(nextLine[36:].strip())
			# short peptide sequence hold by a single line in .gbk file
			else:
				gbk['prots'][-1]['pep'] = nextLine[36:].strip()[:-1]
				continue
			nextLine = next(fp4gbkFile).strip()
			while nextLine[-1] != '"':
				pep.append(nextLine)
				nextLine = next(fp4gbkFile).strip()
			pep.append(nextLine[:-1])
			gbk['prots'][-1]['pep'] = ''.join(pep)
			continue
		if 'ORIGIN' == line.strip():
			seq = []
			nextLine = next(fp4gbkFile)
			while nextLine[:2] != '//':
				seq.append(nextLine[10:].strip().replace(' ',''))
				nextLine = next(fp4gbkFile)
			gbk['seq'] = ''.join(seq).upper()
			fp4gbkFile.close()
			return gbk


# Convert GeneBank protein info (NC_000913.gbk) into FragGeneScan protein file format (NC_000913.fna.faa)
def gbk2fgs4protein(fnaFile, gbkFile, fgsFile):
	# make directories
	makedir(os.path.dirname(fgsFile))

	fp4fgsFile = open(fgsFile, 'w')
	# seqid: eg. 'gi|15644634|ref|NC_000915.1|'
	seqid = rdGbFna(fnaFile)
	# gbk: {'accver':accver, 'gi':gi, 'prots':prots, 'seq':seq, ...}
	gbk = rdGbk(gbkFile)
	if seqid.strip('|').rsplit('|', maxsplit=1)[1] != gbk['accver']:
		print('{} and {} may not the same sequence'.format(fnaFile, gbkFile))
	for prot in gbk['prots']:
		orfStr = '_'.join([str(x) for x in prot['orf']])
		header = '_'.join([seqid, orfStr])
		fasta = fastaFormat(header, prot['pep'])
		print(fasta, file=fp4fgsFile)
	fp4fgsFile.close()


# For each orf, replace seqid with accid
def seqid2accid(mhits):
	mhitsNew = {}
	for seqid in mhits:
		# seqid: 'gi|220910783|ref|NC_011886.1|'
		# accid: 'NC_011886.1'
		accid = seqid.rstrip('|').rsplit('|',1)[-1]
		hitsNew = []
		for hit in mhits[seqid]:
			begin, end, strand = hit['orf'][1:]
			hit['orf'] = (accid, begin, end, strand)
			hitsNew.append(hit)
		mhitsNew[accid] = hitsNew
	return mhitsNew

# For each orf, replace seqid with fileid
def seqid2fileid(mhits):
	mhitsNew = {}
	for seqid in mhits:
		# seqid: 'gi|220910783|ref|NC_011886.1|'
		# accid: 'NC_011886.1'
		accid = seqid.rstrip('|').rsplit('|',1)[-1]
		fileid = accid.split('.', maxsplit=1)[0]
		hitsNew = []
		for hit in mhits[seqid]:
			begin, end, strand = hit['orf'][1:]
			hit['orf'] = (fileid, begin, end, strand)
			hitsNew.append(hit)
		mhitsNew[fileid] = hitsNew
	return mhitsNew

# replace non ATCTU letter with unknown in seq
# seq: single letter DNA/RNA sequence
def qc4fna(seq, unknown='N'):
	return re.sub('[^ATCGU]', unknown, seq.upper())

# dnaFiles: [(file, org), ..., (file, org)]
def rdDNAlist(dnaListFile):
	dnaFiles = []
	fp = open(dnaListFile, 'r')
	for line in fp:
		file = line.strip()
		if file == '' or file[0] == '#':
			continue
		dirs = os.path.dirname(file)
		org = os.path.basename(dirs)
		dnaFiles.append((file, org))
	fp.close()
	return dnaFiles

# Read and return the summarization in file created by outputIndividual() in pred.py
# return: [] or [nis, %genome, bps4is, len4DNA, familySum]
# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
#
# Note: return [1] if sumFileByOrg is the *.sum file created by outputIS4multipleSeqOneFile() in pred.py.
def getSumByOrg(sumFileByOrg, org):
	if os.path.isfile(sumFileByOrg):
		fp2sumByOrg = open(sumFileByOrg, 'r')
	else:
		print('In getSumByOrg() in tools.py: no such file to read', sumFileByOrg)
		return []


	familySum = {}
	sumByOrg = []
	flag4hmpformat = False
	for line in fp2sumByOrg:
		if 'dnaLen' in line:
			flag4hmpformat = True
			break
				
		if line[0] == '#':
			continue # remove comment line
		if line[:5] == 'total':
			sumByOrg = line.split()
			break
		if line[:6] == 'family' or line[:6] == '------':
			continue # title line
		if line.replace('\b','').strip() == '':
			continue # remove blank line
		items = line.split()
		familySum[items[0]] = [int(items[1]), float(items[2]), int(items[3])]
	fp2sumByOrg.close()

	if flag4hmpformat == True:
		# Return and notify the calling code in parent function that the *.sum file (sumFileByOrg) to read
		# is the file format which requires getSumByOrg4hmp() instead of getSumByOrg() to process.
		return [1]
	elif len(sumByOrg) > 0:
		return [int(sumByOrg[1]), float(sumByOrg[2]), int(sumByOrg[3]), int(sumByOrg[4]), familySum]
	else:
		return []

# Read and return the summarization in file written by outputIS4multipleSeqOneFile() in pred.py
# return: [] or [nis, %genome, bps4is, seqLen4bps, familySum, seqLen]
# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
def getSumByOrg4hmp(sumFileByOrg, org):
	if os.path.isfile(sumFileByOrg):
		fp2sumByOrg = open(sumFileByOrg, 'r')
	else:
		print('In getSumByOrg4hmp() in tools.py: no such file to read', sumFileByOrg)
		return []

	familySum = {}
	sumByOrg = []
	dnaLen4is = {}
	# dnaLen4is: {seqid:seqlen, ...}
	for line in fp2sumByOrg:
		line = line.strip()
		if line[0] == '#' or 'family ' in line:
			continue # remove comment line and title line
		if line.replace('\b','').strip() == '':
			continue # remove blank line
		items = line.split()
		if items[1] == 'total':
			#sumByOrg = line.split()
			sumByOrg = items
			break
		seqid = items[0]
		family = items[1]
		nis = int(items[2])
		bps4is = int(items[4])
		seqlen4is = int(items[5])
		#familySum[items[0]] = [int(items[1]), float(items[2]), int(items[3])]
		#print('hello', sumFileByOrg, items[1:5])
		if family not in familySum.keys():
			familySum[family] = [0, 0.0, 0]
		familySum[family] = [familySum[family][0]+nis, 0.0, familySum[family][2]+bps4is]
		if seqid not in dnaLen4is.keys():
			dnaLen4is[seqid] = seqlen4is
	fp2sumByOrg.close()
	if len(sumByOrg) > 0:
		seqLen4is = sum(dnaLen4is.values())
		dnaLen = int(sumByOrg[5])
		for family,value in familySum.items():
			familySum[family][1] = 100 * familySum[family][2] / dnaLen # percentage
		return [int(sumByOrg[2]), float(sumByOrg[3]), int(sumByOrg[4]), seqLen4is, familySum, dnaLen]
	else:
		return []

# file4orgs: {org: files, ..., org: files}
# files: [file4fileid, ..., file4fileid]
def sum4org4hmp(file4orgs, dir4prediction=constants.dir4prediction):
	for org in sorted(file4orgs.keys()):
		# Get summarization from .sum file of each seqid
		# sum4is: {fileid: sum4file, ..., fileid: sum4file}
		# sum4file: [] or [nis, %genome, bps4is, len4DNA, familySum]
		# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
		sum4is = {}
		sum4is4genome = {}
		sum4is4plasmid = {}
		sum4is4phage = {}
		path = os.path.join(dir4prediction, org)
		makedir(path)
		for file4fileid in file4orgs[org]:
			fileid = os.path.basename(file4fileid)
			file = os.path.join(path, '.'.join([fileid, 'sum']))
			sum4is4all = getSumByOrg4hmp(file, fileid)
			# sum4is4all: [] or [nis, %genome, bps4is, dnaLen4is, familySum, dnaLen]
			# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
			if len(sum4is4all) == 0:
				dnaLen = 0
			else:
				dnaLen = sum4is4all[5]
			dnaType = 2 # chromosome DNA

			genome = 0
			genome4is = 0
			plasmid = 0
			plasmid4is = 0
			phage = 0
			phage4is = 0
			dnaLen4phage, dnaLen4plasmid, dnaLen4genome = 0, 0, 0
			if dnaType == 0:
				phage = 1
				dnaLen4phage = dnaLen
			elif dnaType == 1:
				plasmid = 1
				dnaLen4plasmid = dnaLen
			else:
				genome = 1
				dnaLen4genome = dnaLen
			
			# set values for DNA without IS element
			if len(sum4is4all) == 0:
				marker4is = 0
				# nis, %genome, bps4is, dnaLen4is, familySum = 0, 0, 0, dnaLen4is, {}
				sum4is4all = [0, 0, 0, 0, {}]
			else:
				marker4is = 1

			sum4is4all2ext = [genome4is, genome, plasmid4is, plasmid, phage4is, phage]
			sum4is4genome[fileid] = [0, 0, 0, 0, {}, dnaLen4genome] + sum4is4all2ext
			sum4is4plasmid[fileid] = [0, 0, 0, 0, {}, dnaLen4plasmid] + sum4is4all2ext
			sum4is4phage[fileid] = [0, 0, 0, 0, {}, dnaLen4phage] + sum4is4all2ext

			if marker4is == 1:
				if dnaType == 0:
					phage4is = 1
					sum4is4phage[fileid][:5] = sum4is4all
					sum4is4phage[fileid][5+5] = phage4is
				elif dnaType == 1:
					plasmid4is = 1
					sum4is4plasmid[fileid][:5] = sum4is4all
					sum4is4plasmid[fileid][5+3] = plasmid4is
				else:
					genome4is = 1
					sum4is4genome[fileid][:5] = sum4is4all
					sum4is4genome[fileid][5+1] = genome4is
			sum4is[fileid] = sum4is4all[:5] + [dnaLen, genome4is, genome, plasmid4is, plasmid, phage4is, phage]

		# output summarization for all DNAs from organism even there is no IS element found 
		# in the organism.
		sumfile4org = os.path.join(path, 'organism.sum')
		output4sumFull(sum4is, sumfile4org)
		'''
		sumfile4org = os.path.join(path, 'organism4genome.sum')
		output4sumFull(sum4is4genome, sumfile4org)
		sumfile4org = os.path.join(path, 'organism4plasmid.sum')
		output4sumFull(sum4is4plasmid, sumfile4org)
		sumfile4org = os.path.join(path, 'organism4phage.sum')
		output4sumFull(sum4is4phage, sumfile4org)
		'''

familyNames = [
		'IS1',
		'IS110', 
		'IS1182',
		'IS1380',
		'IS1595',
		'IS1634',
		'IS200/IS605',
		'IS21',
		'IS256',
		'IS3',
		'IS30',
		'IS4',
		'IS481',
		'IS5',
		'IS6',
		'IS607',
		'IS630',
		'IS66',
		'IS701',
		'IS91',
		'IS982',
		'ISAS1',
		'ISAZO13',
		'ISH3',
		'ISKRA4',
		'ISL3',
		'ISNCY',
		'new',
		]

# output summarization for IS elements for an organism or multiple organisms
# sum: {seqid: sum4seq, ..., seqid: sum4seq}
# sum4seq: [] or [nis, %genome, bps4is, dnaLen4is, familySum, dnaLen, ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage]
# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
def output4sum(sum4is, outfile):
	fmt4title4families = ' {:>11} {:>11} {:>10}'
	fmt4families = ' {:>11} {:>11.2g} {:>10}'
	fp = open(outfile, 'w')
	fmt4title = '{:<90} {:>6} {:>7} {:>10} {:>10}'
	fmt = '{:<90} {:>6} {:>7.2g} {:>10} {:>10}' 

	# print headline in table
	fp.write(fmt4title.format('organism', 'nIS', '%genome', 'bps4IS', 'dnaLen4is'))
	familyNames.sort()
	for family in familyNames:
		fp.write(fmt4title4families.format(family, '%genome', 'bps4IS'))
	fp.write('\n')

	# print data in table
	# initialize variables
	nis2sum = 0
	bps4is2sum = 0
	dnaLen2sum = 0
	nis4family2sum = {}
	bps4is4family2sum = {}
	for family in familyNames:
		nis4family2sum[family] = 0
		bps4is4family2sum[family] = 0

	# sum: {seqid: sum4seq, ..., seqid: sum4seq}
	# sum4seq: [] or [nis, %genome, bps4is, len4DNA, familySum]
	# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
	for org in sorted(sum4is.keys()):
		sum4org = sum4is[org]
		# initialize data for each family
		nis4org, percent4org, bps4is, dnaLen = 0, 0.0, 0, 0
		familySum = {}
		for family in familyNames:
			familySum[family] = [0, 0.0, 0]
		# get available data for each family
		if len(sum4org) > 0:
			nis4org, percent4org, bps4is, dnaLen =  sum4org[:4]
			for family, value in sum4org[4].items():
				familySum[family] = value

		fp.write(fmt.format(org, nis4org, percent4org, bps4is, dnaLen))
		nis2sum += nis4org 
		bps4is2sum += bps4is
		dnaLen2sum += dnaLen

		# write data for each IS element family
		for family in sorted(familySum.keys()):
			nis, percent, bps4is = familySum[family]
			fp.write(fmt4families.format(nis, percent, bps4is))
			nis4family2sum[family] += nis
			bps4is4family2sum[family] += bps4is
		fp.write('\n')
	
	# print total summary in table
	fp.write(fmt.format('total', nis2sum, (bps4is2sum/dnaLen2sum)*100, bps4is2sum, dnaLen2sum))
	for family in sorted(familySum.keys()):
		fp.write(fmt4families.format(nis4family2sum[family], 
			(bps4is4family2sum[family]/dnaLen2sum)*100, 
			bps4is4family2sum[family]))
	fp.write('\n')

	fp.close()

# output summarization for IS elements for an organism or multiple organisms
# sum: {seqid: sum4seq, ..., seqid: sum4seq}
# sum4seq: [] 
#	or [nis, %genome, bps4is, dnaLen4is, familySum, dnaLen, ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage]
# familySum: {family: [nis, %genome, bps4is, norg4is], ..., family: [nis, %genome, bps4is, norg4is]}
def output4sumFull(sum4is, outfile):
	fmt4title4families = ' {:>11} {:>13} {:>15} {:>13}'
	fmt4families = ' {:>11} {:>13.2g} {:>15} {:>13}'
	fp = open(outfile, 'w')
	fmt4title = '{:<90} {:>6} {:>7} {:>15} {:>15} {:>15} {:>10} {:>7} {:>11} {:>8} {:>9} {:>6}'
	fmt = '{:<90} {:>6} {:>7.2g} {:>15} {:>15} {:>15} {:>10} {:>7} {:>11} {:>8} {:>9} {:>6}' 

	# print headline in table
	fp.write(fmt4title.format(
		'organism', # name of species
		'nIS', # number of ISs occuring in the specific species
		'%genome', # bps4is / dnaLen
		'bps4is', # bps covered by IS
		'dnaLen4is', # bps of DNA sequence where IS occurs
		'dnaLen', # bps of DNA sequences in the specific species
		'ngenome4is', # number of chromosome DNA sequences where IS occurs
		'ngenome', # number of chromosome DNA sequences in the specific species
		'nplasmid4is', # number of plasmid DNA sequences where IS occurs
		'nplasmid', # number of plasmid DNA sequences in the specific species
		'nphage4is', # number of phage DNA sequences where IS occurs
		'nphage' # number of phage DNA sequences in the specific species
		))
	familyNames.sort()
	for family in familyNames:
		fp.write(fmt4title4families.format(
			family, # family name
			family+'_%', # family_bps / dnaLen
			family+'_bps', # bps covered by the specific family
			family+'_s' # number of species where the specific family occurs
			))
	fp.write('\n')

	# print data in table
	# initialize variables
	nis2sum, bps4is2sum, dnaLen4is2sum = 0, 0, 0
	dnaLen2sum = 0
	ngenome4is2sum, ngenome2sum, nplasmid4is2sum, nplasmid2sum, nphage4is2sum, nphage2sum = 0, 0, 0, 0, 0, 0
	nis4family2sum = {}
	bps4is4family2sum = {}
	nstatus4is4family2sum = {}
	for family in familyNames:
		nis4family2sum[family] = 0
		bps4is4family2sum[family] = 0
		nstatus4is4family2sum[family] = 0

	# sum: {seqid: sum4seq, ..., seqid: sum4seq}
	# sum4seq: [] or [nis, %genome, bps4is, len4DNA, familySum]
	# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
	for org in sorted(sum4is.keys()):
		sum4org = sum4is[org]
		# initialize data for each family
		nis4org, percent4org, bps4is, dnaLen4is = 0, 0.0, 0, 0
		dnaLen = 0
		ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage = 0, 0, 0, 0, 0, 0
		familySum = {}
		for family in familyNames:
			familySum[family] = [0, 0.0, 0]
		# get available data for each family
		if len(sum4org) > 0:
			nis4org, percent4org, bps4is, dnaLen4is = sum4org[:4]
			dnaLen, ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage = sum4org[5:]
			for family, value in sum4org[4].items():
				familySum[family] = value

		fp.write(fmt.format(org, nis4org, percent4org, bps4is, dnaLen4is, 
			dnaLen, ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage))
		nis2sum += nis4org 
		bps4is2sum += bps4is
		dnaLen4is2sum += dnaLen4is
		dnaLen2sum += dnaLen
		ngenome4is2sum += ngenome4is
		ngenome2sum += ngenome
		nplasmid4is2sum += nplasmid4is
		nplasmid2sum += nplasmid
		nphage4is2sum += nphage4is
		nphage2sum += nphage

		# write data for each IS element family
		for family in sorted(familySum.keys()):
			nis, percent, bps4is = familySum[family]
			if nis > 0:
				status = 1
			else:
				status = 0
			fp.write(fmt4families.format(nis, 
				percent, bps4is, status
				))
			nis4family2sum[family] += nis
			bps4is4family2sum[family] += bps4is
			nstatus4is4family2sum[family] += status
		fp.write('\n')
	
	# print total summary in table
	if dnaLen2sum == 0:
		percentByBps = 0
	else:
		percentByBps = (bps4is2sum/dnaLen2sum)*100
	fp.write(fmt.format('total', nis2sum, percentByBps, bps4is2sum, dnaLen4is2sum,
		dnaLen2sum, ngenome4is2sum, ngenome2sum, nplasmid4is2sum, nplasmid2sum, nphage4is2sum, nphage2sum))
	for family in sorted(familySum.keys()):
		if dnaLen2sum == 0:
			percentByBps2sum = 0
		else:
			percentByBps2sum = (bps4is4family2sum[family]/dnaLen2sum)*100
		fp.write(fmt4families.format(nis4family2sum[family], 
			percentByBps2sum, 
			bps4is4family2sum[family],
			nstatus4is4family2sum[family]
			))
	fp.write('\n')

	fp.close()

# Read and return the summarization in file written by output4sum(sum4is, outfile)
# return: [] or [nis, %genome, bps4is, len4DNA, familySum]
# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
def getSum(sumFileByOrg, org):
	if os.path.isfile(sumFileByOrg):
		fp = open(sumFileByOrg, 'r')
	else:
		print('No valid IS element was found for', org)
		return []

	familySum = {}
	sumByAll = []
	for line in fp:
		if line[:5] == 'total':
			sumByAll = line.split()
			break
		if line[:8] == 'organism':
			familys = line.split()[5:][0::3] # get IS family names
	familySum = {}
	if len(sumByAll) > 0:
		nis = int(sumByAll[1])
		percent = float(sumByAll[2])
		bps4is = int(sumByAll[3])
		len4DNA = int(sumByAll[4])
		data4familys = sumByAll[5:]
		for i, family in enumerate(familys):
			familySum[family] = [int(data4familys[i*3]), float(data4familys[i*3+1]), int(data4familys[i*3+2])]
		return [nis, percent, bps4is, len4DNA, familySum]
	else:
		return []

# Read and return the summarization in file written by output4sum(sum4is, outfile)
# return: [] 
#	or [nis, %genome, bps4is, dnaLen4is, familySum, dnaLen, ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage]
# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
def getSumFull(sumFileByOrg, org):
	if os.path.isfile(sumFileByOrg):
		fp = open(sumFileByOrg, 'r')
	else:
		print('In getSumFull() in tools.py: no valid IS element was found for', org)
		return []

	familySum = {}
	sumByAll = []
	for line in fp:
		if line[:5] == 'total':
			sumByAll = line.split()
			break
		if line[:8] == 'organism':
			familys = line.split()[12:][0::4] # get IS family names
	familySum = {}
	if len(sumByAll) > 0:
		nis = int(sumByAll[1])
		percent = float(sumByAll[2])
		bps4is = int(sumByAll[3])
		dnaLen4is = int(sumByAll[4])
		dnaLen = int(sumByAll[5])
		ngenome4is = int(sumByAll[6])
		ngenome = int(sumByAll[7])
		nplasmid4is = int(sumByAll[8])
		nplasmid = int(sumByAll[9])
		nphage4is = int(sumByAll[10])
		nphage = int(sumByAll[11])
		data4familys = sumByAll[12:]
		for i, family in enumerate(familys):
			familySum[family] = [int(data4familys[i*4]), float(data4familys[i*4+1]), int(data4familys[i*4+2])]
		return [nis, percent, bps4is, dnaLen4is, familySum, dnaLen, ngenome4is, ngenome, nplasmid4is, nplasmid, nphage4is, nphage]
	else:
		return []

# Return metainfo in fileid
# metainfo: {'dnaType': a, 'dnaLen': dnaLen}, a can be 0,1,2 to represent phage or plasmid or genome DNA
def meta4genome(dir, org, fileid):
	#fnafile = os.path.join(dir, org, fileid+'.fna')
	fnafile = os.path.join(dir, org, fileid)
	seqs = getFastaFull(fnafile)
	# seqs: [(header, seq), ..., (header, seq)]
	# get the first sequence in the fasta file which may contains one or multiple sequences
	header, seq = seqs[0]

	metainfo = {}
	metainfo['dnaLen'] = len(seq)
	if 'phage' in header.lower():
		metainfo['dnaType'] = 0
	elif 'plasmid' in header.lower():
		metainfo['dnaType'] = 1
	else:
		metainfo['dnaType'] = 2

	return metainfo

# Get mDNA from .fna file list like bacteria.fna.list
# Return: [mDNA, fileids]
# mDNA:	{seqid: (org, fileid, sequence), ..., seqid: (org, fileid, sequence)}
# fileids: [(fileid, org), ..., (fileid, org)]
# Input: file
# file: .fna file list, for example, bacteria.fna.list
def fnaFileList2mDNA(filelist):
	fileids = []
	mDNA = {}
	dnaFiles = rdDNAlist(filelist)
	for item in dnaFiles:
		dnafile, org = item
		filename = os.path.basename(dnafile)
		#fileid = filename.rsplit('.', 1)[0]
		fileid = filename
		fileids.append((fileid, org))

		# seq: [(id, seq), ..., (id, seq)]
		seqs = getFasta(dnafile)
		# simply get the first sequence in fasta file.
		if len(seqs) > 0 and len(seqs[0]) > 0:
			mDNA[seqs[0][0]] = (org, fileid, seqs[0][1])
		else:
			print('Warning: no sequence found in', dnafile)
	dir4data = os.path.dirname(os.path.dirname(dnafile))
	return [mDNA, dir4data]


def sum4org(mDNA, dir4data, dir4prediction=constants.dir4prediction):
	orgid = {}
	for seqid in sorted(mDNA.keys()):
		org, fileid, seq = mDNA[seqid]
		if org not in orgid.keys():
			orgid[org] = []

		# Get meta info from  .fna file.
		# Here we simply look for info for genome type.
		# metainfo: {'dnaType': a}, a can be 0,1,2 to represent phage or plasmid or genome DNA
		metainfo = meta4genome(dir4data, org, fileid)

		orgid[org].append([fileid, seqid, metainfo['dnaType'], metainfo['dnaLen']])

	for org in sorted(orgid.keys()):
		# Get summarization from .sum file of each seqid
		# sum4is: {seqid: sum4seq, ..., seqid: sum4seq}
		# sum4seq: [] or [nis, %genome, bps4is, len4DNA, familySum]
		# familySum: {family: [nis, %genome, bps4is], ..., family: [nis, %genome, bps4is]}
		sum4is = {}
		sum4is4genome = {}
		sum4is4plasmid = {}
		sum4is4phage = {}
		path = os.path.join(dir4prediction, org)
		makedir(path)
		for item in orgid[org]:
			fileid, seqid, dnaType, dnaLen = item
			file = os.path.join(path, '.'.join([fileid, 'sum']))

			# process file created by outputIndividual() in pred.py
			sum4is4all = getSumByOrg(file, seqid)
			# if sum4is4all is created from the new fileid.sum format (the format requested by hmp),
			# we only use the first five items in sum4is4all. The sum4is4all created from the old 
			# seqid.sum format requested by bacteria contains only five items which are the same
			# for either fileid.sum or seqid.sum. The 4th item in the sum4is4all created from fileid.sum
			# format is the total length of all sequences with IS elments in the fasta file named by 
			# fileid, and the 6th item is the total length of all sequences in the fasta file, and the
			# fasta file contains one or more sequences. So, the 4th item is same as 6th item in 
			# sum4is4all for the fasta file with one sequence, for example, the fasta file for bactieral 
			# complete genome sequence. The fasta file for whole genome shotgun sequence usually contains
			# more than one contig sequences, for example, the fasta file for HMP metagenome data.
			if len(sum4is4all) == 1:
				# process file created by outputIS4multipleSeqOneFile() in pred.py
				print(file, 'is created by outputIS4multipleSeqOneFile and we hence use getSumByOrg4hmp() to process.')
				sum4is4all = getSumByOrg4hmp(file, seqid)
				sum4is4all = sum4is4all[:5]

			genome = 0
			genome4is = 0
			plasmid = 0
			plasmid4is = 0
			phage = 0
			phage4is = 0
			dnaLen4phage, dnaLen4plasmid, dnaLen4genome = 0, 0, 0
			if dnaType == 0:
				phage = 1
				dnaLen4phage = dnaLen
			elif dnaType == 1:
				plasmid = 1
				dnaLen4plasmid = dnaLen
			else:
				genome = 1
				dnaLen4genome = dnaLen
			
			# set values for DNA without IS element
			if len(sum4is4all) == 0:
				marker4is = 0
				# nis, %genome, bps4is, dnaLen4is, familySum = 0, 0, 0, dnaLen4is, {}
				sum4is4all = [0, 0, 0, 0, {}]
			else:
				marker4is = 1

			sum4is4all2ext = [genome4is, genome, plasmid4is, plasmid, phage4is, phage]
			sum4is4genome[seqid] = [0, 0, 0, 0, {}, dnaLen4genome] + sum4is4all2ext
			sum4is4plasmid[seqid] = [0, 0, 0, 0, {}, dnaLen4plasmid] + sum4is4all2ext
			sum4is4phage[seqid] = [0, 0, 0, 0, {}, dnaLen4phage] + sum4is4all2ext

			if marker4is == 1:
				if dnaType == 0:
					phage4is = 1
					sum4is4phage[seqid][:5] = sum4is4all
					sum4is4phage[seqid][5+5] = phage4is
				elif dnaType == 1:
					plasmid4is = 1
					sum4is4plasmid[seqid][:5] = sum4is4all
					sum4is4plasmid[seqid][5+3] = plasmid4is
				else:
					genome4is = 1
					sum4is4genome[seqid][:5] = sum4is4all
					sum4is4genome[seqid][5+1] = genome4is
			sum4is[seqid] = sum4is4all + [dnaLen, genome4is, genome, plasmid4is, plasmid, phage4is, phage]

		# output summarization for all DNAs from organism even there is no IS element found 
		# in the organism.
		sumfile4org = os.path.join(path, 'organism.sum')
		output4sumFull(sum4is, sumfile4org)
		sumfile4org = os.path.join(path, 'organism4genome.sum')
		output4sumFull(sum4is4genome, sumfile4org)
		sumfile4org = os.path.join(path, 'organism4plasmid.sum')
		output4sumFull(sum4is4plasmid, sumfile4org)
		sumfile4org = os.path.join(path, 'organism4phage.sum')
		output4sumFull(sum4is4phage, sumfile4org)
	

# compute distance between vectors u and v
def distFunction(u, v):
	start1, end1 = u
	start2, end2 = v
	intersect = min(end1, end2) - max(start1, start2) + 1
	if intersect > 0:
		d = 1/intersect
	else:
		d = 10 - intersect
	return d

def distFunctionByoverlap_min(p1, p2):
	a, b = p1
	c, d = p2
	# a <=b and c <= d are satisfied
	intersect = min(b, d) - max(a, c) + 1
	if intersect > 0:
		overlap = float(intersect) / (min(b-a,d-c)+1)
	else:
		overlap = 0.0
	return 1 - overlap

# Get the distribution of a series of integers along a series of windows defined by item and cutoff.
# algorithm: 
# For n items in a list, ilist, we create n windows, (i-cutoff, i+cutoff) where the windows satisfy i-cutoff >= 1.
# Count the number of items in each window.
# Note: some items in ilist might be equal, and the number of windows should be less than n, namely, nwindows <= n.
#
# ilist: [i, ...], i is either left or right boundary 
#
# Return n4windows
# n4windows: {k:n4window, ...}
# k: the boundary, namely, i in ilist
# n4window: number of items (varied boundaries, namely, potential copies) in the window
#
def ncopyByCutoff(ilist, cutoff=0):
	ilist.sort()
	gs = itertools.groupby(ilist)
	windows = {}
	kgs = []
	for k,g in gs:
		# requirement: k >= cutoff
		if k <= cutoff:
			start = 1
		else:
			start = k - cutoff
		end = k + cutoff
		windows[k] = (start, end)
		kgs.append([k,list(g)])
	# n4windows: {k:n4window, ...}
	n4windows = {}
	for k4win,window in windows.items():
		n4windows[k4win] = 0
		for k,g in kgs:
			if window[0] <= k <= window[1]:
				n4windows[k4win] += len(g)
	return n4windows

def getbds4opt4start(n4windows, bds):
	ks4ncopy = list(n4windows.items())
	ks4ncopy.sort(key=operator.itemgetter(1), reverse=True)
	# ks4ncopy: [(k,n4window), ...], sorted by n4window
	bds4opt4k = []
	# bds4opt4k: [bd, ...], boundaries with the most number of items in windows
	# Note: there might be more than one window (k), which have the most number of items in the window
	# bd: (start, end)
	for bd in bds:
		if bd[0] == ks4ncopy[0][0]:
			bds4opt4k.append(bd)
	return bds4opt4k

def getWindowKey4abundance(ilist):
	cutoffs = set()
	ilist.sort()
	# get all possible cutoff values (distance between two items) between any two items in ilist.
	for pair in itertools.combinations(ilist, 2):
		cutoffs.add(pair[1]-pair[0])
	# n4windows4cutoffs: {k:n4win, ...}
	# n4win: is the total number of items in the window centered at k under all possible cutoffs.
	n4windows4cutoffs = {}
	for cutoff in cutoffs:
		# n4windows: {k:n, ...}
		# k: item in ilist
		# n: number of items in the window centered at k under a cutoff
		n4windows = ncopyByCutoff(ilist, cutoff)
		for k,n in n4windows.items():
			if k not in n4windows4cutoffs.keys():
				n4windows4cutoffs[k] = 0
			n4windows4cutoffs[k] += n
	return n4windows4cutoffs
	
def consensusBoundaryByCutoffBySeparated(bds):
	starts = []
	ends = []
	for bd in bds:
		starts.append(bd[0])
		ends.append(bd[1])

	# Get the number of items in each window centered at the specific key (item in starts) of n4windows
	n4windows = getWindowKey4abundance(starts)

	# Sort start (left) boundaries to ensure the first boundary with the most number of items in window
	# is the most left one when multiple items are maximal. 
	# It ensure that the representative bd is the longest bd.
	#
	# Get the key (the specific items with same value in starts) of the window with 
	# the most amount of items within the windows.
	startboundary = max(sorted(n4windows.keys()), key = lambda x: n4windows[x])

	# Get the number of items in each window centered at the specific key (item in ends) of n4windows
	n4windows = getWindowKey4abundance(ends)

	# Sort end (right) boundaries to ensure the first boundary with the most number of items in window
	# is the most right one when multiple items are maximal.
	# It ensure that the representative bd is the longest bd.
	#
	# Get the key (the specific items with same value in ends) of the window with 
	# the most amount of items within the windows.
	endboundary = max(sorted(n4windows.keys(), reverse=True), key = lambda x: n4windows[x])

	return (startboundary, endboundary)

def consensusBoundaryByCutoffByCombined(bds, cutoff=0):
	starts = []
	ends = []
	for bd in bds:
		starts.append(bd[0])
		ends.append(bd[1])

	n4windows = ncopyByCutoff(starts, cutoff)
	# n4windows: {k:n4window, ...}
	# k: the integer, namely, element in left starts (boundaries) 
	# n4window: number of items (integers) in the window

	'''
	# Sort start (left) boundaries to ensure the first boundary with the most number of items in window
	# is the most left one when multiple items are maximal. 
	# It ensure that the representative bd is the longest bd.
	n4windowsSorted = sorted(n4windows.keys())
	# Get the key (the specific items with same value in starts) of the window with 
	# the most amount of items within the windows.
	startboundary = max(n4windowsSorted, key = lambda x: n4windows[x])
	'''
	bds4opt4start = getbds4opt4start(n4windows, bds)


	n4windows = ncopyByCutoff(ends, cutoff)

	'''
	# Sort end (right) boundaries to ensure the first boundary with the most number of items in window
	# is the most right one when multiple items are maximal.
	# It ensure that the representative bd is the longest bd.
	n4windowsSorted = sorted(n4windows.keys(), reverse=True)
	# Get the key (the specific items with same value in ends) of the window with 
	# the most amount of items within the windows.
	endboundary = max(n4windowsSorted, key = lambda x: n4windows[x])
	'''
	bds4opt4end = getbds4opt4start(n4windows, bds)

	# Get the bds with both optimal starts and ends
	set4bds4opt4start = set(bds4opt4start)
	set4bds4opt4end = set(bds4opt4end)
	commonbds = set4bds4opt4start & set4bds4opt4end
	# commonbds: {bd, ...}
	ncommonbds = len(commonbds)
	if ncommonbds > 0:
		commonbdsSortByLen = sorted(commonbds, key = lambda x: x[1]-x[0], reverse = True)
		# Get the representative bd, namely, the longest one among the bds with 
		# both optimal starts and ends (commonbds). 
		startboundary, endboundary = commonbdsSortByLen[0]
	elif cutoff == 0:
		# When ncommonbds == 0, cutoff == 0,
		# get the representative bd, namely, the longest one among the bds wit
		# both optimal starts and ends (commonbds).
		noncommonbds = bds4opt4start + bds4opt4end
		noncommonbdsSortByLen = sorted(noncommonbds, key = lambda x: x[1]-x[0], reverse = True)
		startboundary, endboundary = noncommonbdsSortByLen[0]
	else:
		# When ncommonbds == 0, cutoff > 0,
		# determine the representative bd by more strict cutoff, namely, cutoff - 1
		startboundary, endboundary = consensusBoundaryByCutoff(bds4opt4start + bds4opt4end, cutoff - 1)

	return (startboundary, endboundary)

# based on the bool value of constants.intersected2remove, choose the measure and threshold for 
# clustering and removing intersected IS elements in the same genome sequence
# bd1, bd2: [start, end], boundary of IS element
def chooseMeasure(bd1, bd2):
	intersect = intersection(bd1, bd2)
	if constants.intersected2remove == True:
		measure = intersect
		threshold = constants.min4intersect
	else:
		measure = intersect / min(bd1[1]-bd1[0]+1, bd2[1]-bd2[0]+1)
		threshold = constants.overlap2removeRedundancy
	return (measure, threshold)


def getNewick(node, newick, parentDist, leafNames):
	if node.is_leaf():
		return '{}:{:.2f}{}'.format(leafNames[node.id], parentDist - node.dist, newick)
	else:
		if len(newick) > 0:
			newick = '):{:.2f}{}'.format(parentDist - node.dist, newick)
		else:
			newick = ');'
		newick = getNewick(node.get_left(), newick, node.dist, leafNames)
		newick = getNewick(node.get_right(), ',{}'.format(newick), node.dist, leafNames)
		newick = '({}'.format(newick)
		return newick

# Convert a tree returned by scipy.cluster.hierarchy.to_tree() into a tree in newick format 
# It is a recursive implementation which resultes stack overflow on a big data set. In python3,
# the maximum depth of recursive call is set to 1000.
#
# tree: tree structure returned by scipy.cluster.hierarchy.to_tree() from scipy package
# leafNames: [name_1, ..., name_m], names of the corresponding observations in m observation vectors in n dimensions (m*n array)
def linkageTree2newick(tree, leafNames):
	return getNewick(tree, '', tree.dist, leafNames)

def linkageTree2newick_iter(tree, leafNames):
	newick = ''
	return newick
