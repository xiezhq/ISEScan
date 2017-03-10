# ISEScan
> ISEScan is a python pipeline to identify IS (Insertion Sequence) elements in genome.

ISEScan was developed using Python3. It 1) scanes genome (or metagenome) in fasta format; 2) predicts/translates (using FragGeneScan) genome into proteome; 3) searches the pre-built pHMMs (profile Hidden Markov Models) of transposases (two files shipped with ISEScan; clusters.faa.hmm and clusters.single.faa) against the proteome and identifies the transposase gene in genome; 4) then extends the identified transposase gene into the complete IS (Insertion Sequence) elements based on the common characteristics shared by the known IS elements reported by literatures and database; 5) finally reports the identified IS elements in a few result files (a list of IS elements, sequences of IS elements in fasta format, annotation file in GFF3 format).

## Installation

### Linux:

1. Download the latest ISEScan from https://github.com/xiezhq/ISEScan. The downloaded package is automatically saved as master.zip.

2. Use unzip command to uncompress the zip file:	
`unzip master.zip`

## Pre-required packages and libraries

* Python 3.3.3 or later
* numpy-1.8.0 or later
* scipy-0.13.1 or later
* fastcluster, latest version recommended, https://pypi.python.org/pypi/fastcluster
* FragGeneScan 1.19 or later, https://github.com/COL-IU/FragGeneScan
* HMMER-3.1b2 or later, http://hmmer.org/download.html
* BLAST 2.2.31 or later
* SSW Library, the latest version is not tested with ISEScan and the tested version of SSW library is shipped with ISEScan, please find it at ssw201507 subdirectory.
  * To use the shipped SSW library in ISEScan, please go to ssw201507 and then compile the codes by gcc:
  ```sh
  cc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h
  ```
  * And then copy sswlib.so to the directory of ISEScan and set the search path as:
  ```sh
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:libssw.so
  ```
  * The latest SSW library can be found at https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library.
* biopython 1.62 or later (required by SSW library)

## Configure ISEScan

1. Open constants.py, and find two lines marked with 'Config packages'
2. Modify the paths to FragGeneScan (and phmmer, hmmsearch, blastn, makeblastdb)
3. Save and close constants.py

## Usage example

### Let's try an example, NC_012624.fna.

* The command below scans NC_012624.fna (the genome sequence from Sulfolobus_islandicus_Y_N_15_51), and outputs all results in prediction directory:
```sh
python3 isescan.py NC_012624.fna proteome hmm
```

* Wait for its finishing. It may take a while as ISEScan uses the HMMER to scan the genome sequences and it will use 496 profile HMM models to scan each protein sequence (predicted by FragGeneScan) in the genome sequence. HMMER searching is usually more sensitive but slower than the regular BLAST searching for remote homologs.

* After ISEScan finish running, you can find three important files in prediction directory, NC_012624.fna.sum, NC_012624.fna.gff, NC_012624.fna.is.fna. The summarization of IS copies for each IS family is in NC_012624.fna.sum, NC_012624.fna.gff list each IS element copy and its TIR. NC_012624.fna.is.fna holds the nucleic acid sequence of each IS element copy.

### Tips:
* ISEScan will run much faster if you run it on the same genome sequence more than once (e.g., trying different optimal parameters of near and far regions (see our paper [...] for the definitions of near and far regions)) to search for IS elements in your genome). The reason is that it skips either FragGeneScan or both FragGeneScan and phmer/hmmsearch steps which are most time-consuming steps in ISEScan pipeline.
* If you prefer ISEScan recalculating the the results, you can simply remove the proteome file and HMMER search results which are related to your genome sequence file name. For example, you can delete NC_012624.fna.faa in proteome directory and clusters.faa.hmm.NC_012624.fna.faa and clusters.single.faa.NC_012624.fna.faa in hmm directory, and then rerun it:   
`python3 isescan.py NC_012624.fna proteome hmm`

## Release History 
* 1.2
  * CHANGE: pHMMs `clusters.faa.hmm` and `clusters.single.faa`, both files are now built upon the curated ACLAME dataset (ACLAME is a mobile genetic element database.)
* 1.1.1
  * Add option in `constants.py` to report either complete IS elements or both complete and partial IS elements
* 1.0
  * The first proper release

## Meta
Zhiqun Xie – @col-iu](https://col-iu.wikispaces.com/Zhiqun+Xie) – zhiqxie@iu.edu

Distributed under the GNU General Public License.

[https://github.com/xiezhq](https://github.com/xiezhq)
