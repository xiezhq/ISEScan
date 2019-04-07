import argparse

import constants

#dir = '/nfs/omics/home/data/insertion_sequence/results4MGEScan-IS/prediction'
#filename = 'benchmark_20160613_113508' # E. coli K_12 MG1655 data set
#filename = 'benchmark_20160613_095054' # ISBrowser data set
#ifile = os.path.join(dir, filename)

def evaluateFDR(args4evaluateFDR):
	ifile = args4evaluateFDR['benchfile']
	with open(ifile, 'r') as fp:
		reachfdr = False
		lines = []
		for line in fp:
			if reachfdr == False:
				if '# Full list of false discovery hits' in line:
					reachfdr = True
				else:
					continue
			if '# Full list of hits' in line:
				break
			line = line.strip()
			if len(line) < 1 or line[0] == '#':
				continue
			if line[:3] == 'NC_' or line[:5] == 'irSim':
				lines.append(line)
	# process the lines related to false discovery hits
	mmge = []
	# mmge: [mge, ...]
	# mge: {'seqid':seqid, 'start':start, 'end':end, 'family':family, 'evalue':evalue, 'ncopy4is':ncopy4is, 
	#	'orfStart':orfStart, 'orfEnd':orfEnd, 'tir':[]}
	# tir: [irId, irLen, gap, start1, end1, start2, end2]
	for linepair in zip(*[iter(lines)]*2):
		line1, line2 = linepair
		mge = {}
		items = line1[:109].strip().split()
		mge['seqid'], mge['strand'], mge['start'], mge['end'], mge['family'] = items
		mge['start'], mge['end'] = int(mge['start']), int(mge['end'])
		mge['family'] = mge['family'].split('_',maxsplit=1)[0]
		mge['evalue'] = float(line1[197:207])
		mge['ncopy4is'] = int(line1[301:309])
		mge['orfStart'] = int(line1[372:384])
		mge['orfEnd'] = int(line1[385:397])
		index1 = line2.find('tir')
		items4tir = line2[index1+2:].split(',')
		if len(items4tir) > 1:
			tir = [int(item) for item in items4tir[1:8]]

		else:
			tir = []
		mge['tir'] = tir
		mmge.append(mge)
	# evaluate false discovery hits and classify IS elements into full-length element with high confidence 
	# or unknown case, based on three parameters: 
	# isLen < maxLen4is4family, end1 < orfstar and orfend < start2, 10 <= tir
	mges = []
	for mge in mmge:
		family = mge['family'].split('_', maxsplit=1)[0]
		isLen = mge['end'] - mge['start'] + 1
		tir = mge['tir']
		# good IS elements in fdr list: 
		# ncopy4is >= 2, minLen4family <= isLen <= maxLen4family, start <= orfStart and orfEnd >= end
		if mge['ncopy4is'] < 2 or isLen > constants.minMaxLen4is[family][1] or mge['orfStart'] < mge['start'] or mge['orfEnd'] > mge['end']:
		# good IS elements in fdr list: 
		# ncopy4is >= 2, minLen4family <= isLen <= maxLen4family
		#if mge['ncopy4is'] < 2 or isLen > constants.minMaxLen4is[family][1]:
			mge['isType'] = 'n'
		else:
			mge['isType'] = 'y'

		mges.append(mge)
	# output false discovery hits with isType labels
	print('{:<10} {:<12} {:>6} {:>6} {:>6} {:>4} {:>5} {:>8} {:>9} {:>6} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}'.format(
		'seqid', 'family', 'minLen', 'maxLen', 'isLen', 'irId', 'irLen', 'ncopy4is', 'evalue', 'isType', 'isStart', 'start1', 'end1', 'orfStart', 'orfEnd', 'start2', 'end2', 'isEnd'))
	for mge in mges:
		family = mge['family']
		minLen, maxLen = constants.minMaxLen4is[family]
		tir = mge['tir']
		if len(tir) > 0:
			irId, irLen, gap, start1, end1, start2, end2 = tir
		else:
			irId, irLen, gap, start1, end1, start2, end2 = 0, 0, 0, 0, 0, 0, 0 
		print('{:<10} {:<12} {:>6} {:>6} {:>6} {:>4} {:>5} {:>8} {:>9.2g} {:>6} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}'.format(
			mge['seqid'], mge['family'], minLen, maxLen, mge['end']-mge['start']+1, irId, irLen, # 'seqid', 'family', 'isLen', 'irId', 'irLen', 
			mge['ncopy4is'], mge['evalue'], mge['isType'], # 'evalue', 'isType', 
			mge['start'], start1, end1, mge['orfStart'], # 'isStart', 'start1', 'end1', 'orfStart', 
			mge['orfEnd'], start2, end2, mge['end'], # 'orfEnd', 'start2', 'end2', 'isEnd'
			))


if __name__ == "__main__":
	descriptionStr = 'Summarize the false discovery results in ISEScan bechmark'
	parser = argparse.ArgumentParser(description = descriptionStr)
	helpStr = 'benchmark result file'
	parser.add_argument('benchfile', help = helpStr)
	args = parser.parse_args()
	args4evaluateFDR = {
				'benchfile':		args.benchfile,
				}
	evaluateFDR(args4evaluateFDR)
