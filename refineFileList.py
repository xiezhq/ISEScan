import os.path

#dnalistfile = '../hmp.ref.final.assembly.file.list.omics'
dnalistfile = '../bacteria.fna.list.omics'

outfile = 'fna.out.list.omics'

fp4outfile = open(outfile, 'r')

fnanames = []
for line in fp4outfile:
	line = line.strip()
	filename = os.path.basename(line)
	fnanames.append(filename[:-4]) # remove '.out'
fp4outfile.close()

dnalistfileNew = dnalistfile + '.new'
fp4dnalistfileNew = open(dnalistfileNew, 'w')
fp4dnalistfile = open(dnalistfile, 'r')
for line in fp4dnalistfile:
	skipflag = False
	for fnaname in fnanames:
		if fnaname in line:
			skipflag = True
			break
	if skipflag == True:
		continue
	fp4dnalistfileNew.write(line)
fp4dnalistfile.close()
fp4dnalistfileNew.close()

