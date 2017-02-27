import argparse
import tools

parser = argparse.ArgumentParser(description='Split a fasta file containing multiple sequences into multiple individual fasta files')
parser.add_argument('input_fasta')
parser.add_argument('output_directory')
args = parser.parse_args()
tools.split_tandem_fasta(args.input_fasta, args.output_directory)
