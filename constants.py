import os.path

# Config packages
#
# Config the paths of required packages in order that ISEScan can find the right commands
# on your computer.
# 
'''
# FragGeneScan
FragGeneScan = '/u/zhiqxie/informatics/inst/FragGeneScan1.19/run_FragGeneScan.pl'
# Hmmer
phmmer = '/u/zhiqxie/informatics/inst/hmmer-3.1b2/bin/phmmer'
hmmsearch = '/u/zhiqxie/informatics/inst/hmmer-3.1b2/bin/hmmsearch'
# Blast 
blastn = '/l/ncbi-blast/bin/blastn'
makeblastdb = '/l/ncbi-blast/bin/makeblastdb'
'''
FragGeneScan = '/N/u/zhiqxie/Mason/informatics/inst/FragGeneScan1.19/run_FragGeneScan.pl'
phmmer = '/N/u/zhiqxie/Mason/informatics/inst/hmmer-3.1b2/bin/phmmer'
hmmsearch = '/N/u/zhiqxie/Mason/informatics/inst/hmmer-3.1b2/bin/hmmsearch'
blastn = 'blastn'
makeblastdb = 'makeblastdb'
#
# Config packages


# When translateGenome is True, pipepline will predict and translate genes 
# from genome sequence into protein database (file in fasta format) 
# using FragGeneScan program.
translateGenome = True
# When translateGenome is False, pipeline will use the protein database (.faa file)
# from ncbi genome database, and the corresponding .fna and .ptt files in the same folder
# are required to map proteins to genome locations.
#translateGenome = False

# remove ISs with length < isMin if True
# isMin is defined in this file, it is 400. The full-length IS element is usually longer than 400 bp.
removeShortIS = True
#removeShortIS = False

# set temporary directory used by ISEScan
#tmpdir = '/N/u/zhiqxie/Karst/is/isescan/tmpdir'
tmpdir = '/N/dc2/scratch/zhiqxie/insertion_sequence/tmpdir'

#path2results = '/N/u/zhiqxie/Karst/is/isescan/results'
path2results = '/N/dc2/scratch/zhiqxie/insertion_sequence/results4hmp'
#path2results = ''
dir4prediction = os.path.join(path2results, 'prediction')

# peptide sequences of single-member clusters, which is used by phmmer in hmmer
#file4clusterSeqFile4phmmer = 'clusters.single.faa'
file4clusterSeqFile4phmmer = '/N/u/zhiqxie/Karst/is/isescan/clusters.single.faa'
# profile HMMs of multiple-member clusters, which is used by hmmsearch in hmmer
#file4clusterHMM = 'clusters.faa.hmm'
file4clusterHMM = '/N/u/zhiqxie/Karst/is/isescan/clusters.faa.hmm'

# blast database will be put here
dir4blastout = os.path.join(path2results, 'blastout')

# Optimal values for SSW to find TIR in ISfinder database
# (gapopen, gapextend, match, mismatch)
#
# Optimal filter when aligning two sequences with length = maxLenIR
filters4ssw4isMax = [(1, 10, 4, 5)] # giving the greatest number of matched IS elements and 
				# the greatest number of matched best IS elements
filters4ssw4trial = [(2, 6, 2, 2)] # trial filter to stop alignment from creating the consecutive gaps

minMaxLen4is =	{
		'IS1': (732, 4601),
		'IS110': (969, 4105),
		'IS1182': (1330, 1980),
		'IS1380': (1474, 4160),
		'IS1595': (701, 7915), 
		# IS1016V5 (272 bp) is a deleted variant of IS1016V6: 242/711 bp. IS1016V4 (672 bp) is a deleted variant of IS1016V6: 673/711 bp.
		# Then ISMha1 (701 bp) is the shortest member in family IS1595.

		'IS1634': (1511, 2089),
		'IS200/IS605': (407, 2223),
		'IS21': (1924, 3533),
		'IS256': (1124, 1629),
		'IS3': (435, 1814),
		'IS30': (1027, 8273),
		'IS4': (521, 5396),
		'IS481': (553, 3451),
		'IS5': (789, 5396),
		'IS6': (696, 1648),
		'IS607': (1415, 2607),
		'IS630': (895, 2009),
		'IS66': (1364, 3481), 
		# IS867 has about 75 % homology with IS866. IS866 is 2716 bp.
		# Then ISMno3 (1364 bp) is the shortest member in family IS66.
		'IS701': (1016, 2207),
		'IS91': (712, 2604),
		'IS982': (845, 1282),
		'ISAS1': (1139, 3041),
		'ISAZO13': (1284, 2171),
		'ISH3': (1200, 1509),
		'ISKRA4': (1164, 3746),
		'ISL3': (536, 9109),
		'ISNCY': (786, 3989),
		}

# peptide and ORF lengths of tpases in ISfinder
# shortest tpase ORF (bp), longest tpase ORF (bp)
# shortest peptide ORF (bp) among all peptides in IS_PEP record for each IS element in ISfinder,
# To be added: shortest tpase (aa), longest tpase (aa), 
# ORF = tpase + stopcodon
minMax4tpase =	{
		'IS1': (666, 1119, 252),
		'IS110': (603, 1380, 156),
		'IS1182': (822, 1731, 570),
		'IS1380': (1158, 1554, 1158),
		'IS1595': (576, 1158, 426), 
		'IS1634': (1314, 1875, 1314),
		'IS200/IS605': (366, 1482, 147),
		'IS21': (882, 1758, 231),
		'IS256': (990, 1389, 990),
		'IS3': (441, 1581, 120),
		'IS30': (540, 1419, 189),
		'IS4': (570, 1629, 219),
		'IS481': (447, 1794, 447),
		'IS5': (360, 1908, 75),
		'IS6': (528, 1062, 246),
		'IS607': (768, 1653, 453),
		'IS630': (510, 1194, 318),
		'IS66': (354, 1695, 165), 
		'IS701': (921, 1410, 921),
		'IS91': (648, 1548, 648),
		'IS982': (627, 981, 429),
		'ISAS1': (594, 1329, 189),
		'ISAZO13': (1203, 2094, 513),
		'ISH3': (573, 1206, 549),
		'ISKRA4': (1047, 1719, 114),
		'ISL3': (414, 1716, 408),
		'ISNCY': (573, 1815, 123),
		}

# allowed minimal and maximal and optimal values of the length of TIR sequence for each family
# Here, the optimal values are the empirical parameter based on the observations on ISfinder database.
# The 4th collumn is marker indicating whether the family always has TIR or not, 1 for yes and 0 for no 
# and -1 for either (in the family, some members have tir but others have no tir).
minMax4tir = {
		'IS1': (8, 67, 14, 1),
		'IS110': (2, 31, 14, -1), 
		'IS1182': (8, 44, 10, 1),
		'IS1380': (7, 39, 10, 1),
		'IS1595': (10, 43, 15, 1),
		'IS1634': (11, 32, 12, 1),
		'IS200/IS605': (10000, 0, 10000, 0), # prevent program from finding any tir with irLen > 0
		'IS200/IS605_8': (11, 11, 11, 1), # cluster 8 (cdhit30) of IS200/IS605 has tir with irLen == 0 or irLen == 11
		#'IS200/IS605': (11, 11, 11, -1), # cluster 8 (cdhit30) of IS200/IS605 has tir with irLen == 0 or irLen == 11
		'IS21': (8, 76, 10, 1),
		'IS256': (8, 48, 15, 1),
		'IS3': (7, 54, 10, -1),
		'IS30': (11, 50, 12, 1),
		'IS4': (8, 67, 12, 1),
		'IS481': (5, 52, 10, 1),
		'IS5': (7, 45, 14, 1),
		'IS6': (12, 36, 14, 1),
		'IS607': (12, 46, 12, -1),
		'IS630': (3, 92, 11, 1),
		'IS66': (11, 144, 11, 1),
		'IS701': (12, 38, 12, 1),
		'IS91': (11, 21, 11, -1),
		'IS982': (11, 35, 11, 1),
		'ISAS1': (12, 34, 12, 1),
		'ISAZO13': (18, 48, 18, 1),
		'ISH3': (11, 31, 15, 1),
		'ISKRA4': (15, 40, 18, 1),
		'ISL3': (6, 50, 11, 1),
		'ISNCY': (4, 52, 13, -1),
	}
# ssw will use minMax4tir[2] as minimal length of the alignement of two tir sequences 
# if useOPTtir == True else minMax[0] as minimal length of the alignment of two tir sequences.
#useOPTtir = True
useOPTtir = False

# the minimum of rations of irId/irLen
minIrIdentity = 0.4
# optimal ration of irId/irLen
optIrIdentity = 0.6
# stringent irId/irLen, which is usually required when irLen < 5(stringentShortestIR) or irLen > 55(stringentLongestIR)
stringentIrIdentity = 0.7

# maximum distance (bp) between two neighboring orfs (including +/- strand) within one IS element
# Statistics from isfinder:
# 764 IS elements with multiple ORFs with clear coordinates in ORF records,
# 405 with distBetweenORFs >=0, 
# 1/405 with dist >= 1000, 6/405(1%) with dist >= 500, 14/405(3%) with dist >= 400, 
# 22/405(5%) with dist >= 300, 31/405(8%) with dist >= 250, 44/405(11%) with dist >= 200, 
# 90/405(22%) with dist >= 100, 202/405(50%) with dist >= 55, 214/405(53%) with dist >= 50
#
# not to merge
#maxDistBetweenOrfs = -1
# merge ORFs with gap = 0, ('NC_000913.3', 4518418, 4519014, '+') and ('NC_000913.3', 4519015, 4519224, '+')
#maxDistBetweenOrfs = 0
# merge ORFs with gap <= 100 bps
maxDistBetweenOrfs = 100

# In isfinder, 3891 IS elements with both lORF2TER and rORF2TER >= 0, 
# 36/3891(1%) with lORF2TER >= 500, 177/3891(5%) with lORF2TER >= 250,
# 51/3891 with rORF2TER >= 500, 232/3891 with rORF2TER >= 250
# ~99% IS elements in ISfinder has lORF2TER/rORF2TER less than 500 bps
# ~95% IS elements in ISfinder has lORF2TER/rORF2TER less than 250 bps
#
# switch maxDist4ter2orf between 500 and 250
#maxDist4ter2orf = 250
maxDist4ter2orf = 500
outerDist4ter2tpase = (150,500)

# Minimum distance (bp) from near ends of IS element to the nearest ORF, namely, 
# the length of the shortest linker between TIR and the nearest ORF.
minDist4ter2orf = -150
#minDist4ter2orf = -50
#
# There is no linkder (space) between TIR and the nearest ORF.
#minDist4ter2orf = 1

# The strand does not matter when extracting two terminal sequences to align, namely,
# which sequence is the first sequence in pairwise alignement does not make sense.
#splitAlign2orf = True
splitAlign2orf = False

# IS elements with identicalBases/lengthOfAlignment > sim4iso are regarded as the same IS element (isoform)
# ISO records in ISfinder: 
# Isoforms have been defined as elements which share in the first instance more than 95% identity 
# at the level of their transposase protein sequence or otherwise 90% at the DNA level.
#sim4iso = 0.85
sim4iso = 0.9
#
# SIM4ISO = sim4iso * 100, used by blastn search to get copy number of hit
#SIM4ISO = 85
SIM4ISO = 90
#
# similarity cutoff for protein sequence
aaSim4iso = 0.95
aaSIM4ISO = 95

# Two neighboring sequences with overlap >= min4overlap are deemed overlapped.
min4overlap = 0.5 # 50%

# Two sequences with intersect >= min4intersect are deemed intersect.
#min4intersect = 100 # 100 bp, namely, 33 aa or so.
min4intersect = 1 # 1 bp.

# two neighoring segments with overlap > overlap2removeRedundancy are considered overlapped (redundant)
overlap2removeRedundancy = 0.5 # 50%
#overlap2removeRedundancy = 0.99999999999 # 100%

# use min4intersect if True else overlap2removeRedundancy as the threshold to 
# turn on clustering and remove intersected ISs/hits except the representative in a cluster.
#intersected2remove = True
intersected2remove = False

# hits with evalue <= min4evalue are defined as the final hits.
min4evalue = 1e-10
#min4evalue = 1e-5

# more strict evalue and tir are required for single copy hits
evalue4singleCopy = 1e-50
irSim4singleCopy = 0.85 # irId/irLen

# E-value cutoff for filtering hits returned by HMM search
evalue2filterHMMhits = min4evalue
#evalue2filterHMMhits = 10 # do not filter out any hits returned by HMM search

# width of line in fasta file created by us
fastaLineWidth = 60

# complementary table for DNA
#------------------------------------------
# Code	Represents		Complement
# A 	Adenine			T
# G 	Guanine			C
# C 	Cytosine 		G
# T 	Thymine 		A
# Y 	Pyrimidine (C or T) 	R
# R 	Purine (A or G) 	Y
# W 	weak (A or T) 		W
# S 	strong (G or C) 	S
# K 	keto (T or G) 		M
# M 	amino (C or A) 		K
# D 	A, G, T (not C) 	H
# V 	A, C, G (not T) 	B
# H 	A, C, T (not G) 	D
# B 	C, G, T (not A) 	V
# X/N 	any base 		X/N
# - 	Gap 			-
#------------------------------------------
#na1u = 'ATCGN'
na1u = 'ATCGNRYWSKMDVHBX'
#na2u = 'TAGCN'
na2u = 'TAGCNYRWSMKHBDVX'
#na1l = 'atcgn'
na1l = 'atcgry'
#na2l = 'tagcn'
na2l = 'tagcyr'
#na1ul = 'ATCGNatcgn'
na1ul = 'ATCGRYatcgry'
#na2ul = 'TAGCNtagcn'
na2ul = 'TAGCYRtagcyr'

# The Genetic Codes
# Refer to http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
# The Bacterial, Archaeal and Plant Plastid Code (transl_table=11).
table11 = {	
		'starts': ('TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'),

		'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
		'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
		'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
		'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',

		'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 
		'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
		'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
		'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',

		'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
		'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
		'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
		'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',

		'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 
		'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
		'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
		'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
		}
gene2pepTable = {'11': table11}

# default number of processes to use in calculation if it is not given
nproc = 16
#nproc = 1
# default number of threads to use in calculation if it is not given
#nthread = 2
nthread = 16
