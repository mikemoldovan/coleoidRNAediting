"""
Separate assembly file for blast
Input:
1. Trinity assembly file
2. number of desired threads N 

Output
1. N fasta files with transcripts
2. N scripts for the qsub
"""

from optparse import OptionParser
from os import system

def countlines(fafile):
	count = 0
	for s in open(fafile):
		if s[0] == '>':
			count += 1
	return count


def construct_datafiles(fafile, transcript_num, thread_num, header):
	single_file_num = transcript_num/thread_num
	count = 0
	filecount = 0
	outfasta = open("temp{}.fa".format(filecount),'w')
	for s in open(fafile):
		if s[0] == '>':
			count += 1
		if count == single_file_num:
			outfasta.close()
			filecount += 1
			outfasta = open("temp{}.fa".format(filecount),'w')
			count = 0
		outfasta.write(s)
	outfasta.close()

	for i in range(filecount+1):
		blastsh = open("blastx_{}.sh".format(i), 'w')
		blastsh.write(header)
		blastsh.write("export PATH=$PATH:~/tools/ncbi-blast-2.7.1+/bin/\n")
		blastsh.write("blastx -query temp{}.fa -db ~/databases/uniprot/uniprot_sprot.fasta  -evalue 1e-6 -out blastout{}.txt -outfmt 6 -num_alignments 1\n".format(i,i))
		blastsh.close()


parser = OptionParser()
parser.add_option("-i", "--fafile", help="File in FASTA format, TRINITY output")
parser.add_option("-n", "--thread_num", help="Number of partitions (threads)")
opt, args = parser.parse_args()


header = """#! /bin/bash
#$ -l h_rt=96:00:00
#$ -cwd
#$ -l mem_free=100M
#$ -p -1000
"""

transcript_num = countlines(opt.fafile)
construct_datafiles(opt.fafile, transcript_num, eval(opt.thread_num), header)
