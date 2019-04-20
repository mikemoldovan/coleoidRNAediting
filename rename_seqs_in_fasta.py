"""
rename_seqs
"""

from Bio import SeqIO
from optparse import OptionParser
from os import listdir


def rename_seqs(fafile_name):
	prefix = '.'.join(fafile_name.split('/')[-1].split('.')[:-1])
	count = 0
	with open(fafile_name) as inhandle, open(prefix+"_cus.fa", 'w') as outhandle, open(prefix+"_genetable.txt", 'w') as genetable:
		for record in SeqIO.parse(inhandle, "fasta"):
			outhandle.write(">"+prefix+'_'+str(count)+'\n')
			s = ""
			for c in record.seq:
				s += c
			outhandle.write(s + '\n')
			genetable.write(prefix+'_'+str(count)+'\t'+record.id+'\n')
			count += 1


parser = OptionParser()
parser.add_option("-i", "--fadir", help="directory with fasta files")
opt, args = parser.parse_args()

for f in listdir(opt.fadir):
	if f.split('.')[-1] == "fa" and "_cus" not in f:
		rename_seqs(opt.fadir + f)
