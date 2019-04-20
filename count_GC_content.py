"""
count_GC_content

input:
1. Input FASTA file
2. frame size

output:
distribution of GC-content
"""

from optparse import OptionParser
from Bio import SeqIO

def count_GC(seq):
	print seq.count('G') + seq.count('C')


def main(infile, framesize):
	with open(infile) as inhandle:
		for record in SeqIO.parse(inhandle, "fasta"):
			count = 0
			seq = ""
			for c in record.seq:
				seq += c
				count += 1
				if count%framesize == 0:
					count_GC(seq)
					count = 0
					seq = ""


parser = OptionParser()
parser.add_option("-i", "--infile", help="Input FASTA file")
parser.add_option("-s", "--framesize", help="frame size")
opt, args = parser.parse_args()

main(opt.infile, eval(opt.framesize))