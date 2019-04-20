"""
infer_orf_crds

Infers ORF coordinates from Liscovitch-Brauer transcriptome data

Input:
File with transcriptome (both '+' and '-' strands)

Output:
File with coordinates on '+' strands
"""

#from sys import stderr
from Bio import SeqIO
from optparse import OptionParser


def main(transcr_file, spec_id):
	print "#seq_id\torfstart\torfend\tstrand"
	with open(transcr_file) as tfile:
		for record in SeqIO.parse(tfile, "fasta"):
			descr = record.description
			descr = descr.split()
			orfstart = eval(descr[2])
			orfend = eval(descr[4])
			if (orfend - orfstart + 1) % 3 != 0:
				print descr
			strand = descr[6]
			if strand == '+':
				orfstart -= 1
				print "orfs_" + spec_id + '_' + descr[0] + '\t' + str(orfstart) + '\t' + str(orfend) + '\t' + strand
			elif strand == '-':
				l = len(record.seq)
				orfstart_n = l - orfend
				orfend_n = l - orfstart + 1
				print "orfs_" + spec_id + '_' + descr[0] + '\t' + str(orfstart_n) + '\t' + str(orfend_n) + '\t' + strand


parser = OptionParser()
parser.add_option("-i", "--transcr_file", help="File with transcriptome (both '+' and '-' strands)")
parser.add_option("-s", "--spec_id", help="Species ID")
opt, args = parser.parse_args()


main(opt.transcr_file, opt.spec_id)