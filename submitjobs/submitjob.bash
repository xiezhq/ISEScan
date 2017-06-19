#!/bin/bash

#dir4job='/N/dc2/scratch/zhiqxie/insertion_sequence/job4bacteria'
dir4job='/N/dc2/scratch/zhiqxie/insertion_sequence/job4hmpHMASM'
path2proteome='/N/dc2/scratch/zhiqxie/insertion_sequence/output4FragGeneScan1.19_illumina_5'
path2hmm='/N/dc2/scratch/zhiqxie/insertion_sequence/output4hmmsearch_illumina_5_cdhit30'
#jobscript='/N/u/zhiqxie/Karst/submitjobs/hmp.hmasm.pbs'

hmplist="$1"
jobscript="$2"

while read seqfile; do
	seqfilename=$(basename "$seqfile")
	echo qsub -d "$dir4job" -N "$seqfilename" -o "$dir4job/$seqfilename".log -v seqfile="$seqfile",path2proteome="$path2proteome",path2hmm="$path2hmm" "$jobscript"
	qsub -d "$dir4job" -N "$seqfilename" -o "$dir4job/$seqfilename".log -v seqfile="$seqfile",path2proteome="$path2proteome",path2hmm="$path2hmm" "$jobscript"
done < "$hmplist"
