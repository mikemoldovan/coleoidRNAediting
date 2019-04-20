"""
count_mutations
Builds alignments and counts mutations in alignment blocks,
separately counting mutations in editing sites

Input:
1. Editing sites table for organism '1'
2. Editing sites table for organism '2'
3. Transcriptome fasta for organism '1'
4. Transcriptome fasta for organism '2'
5. Proteinortho file
6. Job prefix

Output:
1. Alignment blocks (nucleotide) for transcriptomes of two
organisms
2. Three relative mutation scores (A->G, A->C, A->T, both-sided)
"""

from os import system
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq


class AlignBlock():
	def __init__(self, gene_id1=None, gene_id2=None, q_start=None, q_seq=None, q_end=None, q_sites=None, s_start=None, s_seq=None, s_end=None, s_sites=None):
		self.gene_id1 = gene_id1
		self.gene_id2 = gene_id2
		self.q_start = q_start
		self.q_seq = q_seq
		self.q_end = q_end
		self.q_sites = q_sites
		self.s_start = s_start
		self.s_seq = s_seq
		self.s_end = s_end
		self.s_sites = s_sites

	def print_str(self, filevar):
		same_directions = True
		same_directions_flag = 'y'
		if (self.q_start > self.q_end and self.s_start < self.s_end) or (self.q_start < self.q_end and self.s_start > self.s_end):
			same_directions = False
			same_directions_flag = 'n'
		self.q_seq = self.q_seq[min(self.q_start, self.q_end) - 1: max(self.q_start, self.q_end)]
		if same_directions:
			self.s_seq = self.s_seq[min(self.s_start, self.s_end) - 1: max(self.s_start, self.s_end)]
			align_sites_2 = ','.join(map(str,[i - min(self.s_end, self.s_start) for i in self.s_sites]))
		else:
			self.s_seq = self.s_seq.reverse_complement()[-max(self.s_start, self.s_end): -min(self.s_start, self.s_end) + 1]
			align_sites_2 = [i - min(self.s_end, self.s_start) for i in self.s_sites]
			l = len(self.s_seq)
			align_sites_2 = [l - i for i in align_sites_2]
		align_sites_1 = ','.join(map(str,[i - min(self.q_start, self.q_end) for i in self.q_sites]))
		if not align_sites_1:
			align_sites_1 = '*'
#		align_sites_2 = ','.join(map(str,[i - min(self.s_end, self.s_start) for i in self.s_sites]))
		if not align_sites_2:
			align_sites_2 = '*'
		outstr = [self.gene_id1, self.gene_id2, self.q_start, self.q_seq, self.q_end, self.s_start, self.s_seq, self.s_end, align_sites_1, align_sites_2, same_directions_flag]
		outstr = '\t'.join(map(str, outstr))
		filevar.write(outstr + '\n')



def genetable_dict_build(fafile_name):
	genetable_dict = dict()
	with open('.'.join(fafile_name.split('.fa')[:-1]) + "_genetable.txt") as gentab:
		for s in gentab:
			s = s.strip().split()
			genetable_dict[s[0]] = s[1] # new name : old name
	return genetable_dict


def editing_site_dict_build(ed_site_table):
	editing_site_dict = dict()
	with open(ed_site_table) as edsites:
		for s in edsites:
			if s[0] == '#':
				continue
			s = s.strip().split()
			try:
				editing_site_dict[s[0]].append(eval(s[3]))
			except:
				editing_site_dict[s[0]] = [eval(s[3])]
	return editing_site_dict

def proteinortho_file_parse(proteinortho_file, fafile_name_1, fafile_name_2):
	homol_dict = dict()
	genetable_dict_1 = genetable_dict_build(fafile_name_1)
	genetable_dict_2 = genetable_dict_build(fafile_name_2)
	with open(proteinortho_file) as po_file:
		for s in po_file:
			s = s.strip().split('\t')
			if s[0][0] == '#':
				name1 = '.'.join(fafile_name_1.split('/')[-1].split('.')[:-1]) + '_cus.fa'
				name2 = '.'.join(fafile_name_2.split('/')[-1].split('.')[:-1]) + '_cus.fa'
				index_1 = s.index(name1)
				index_2 = s.index(name2)
				continue
			if (s[index_1] != '*') and (s[index_2] != '*'):
				genename_1 = genetable_dict_1[s[index_1]]
				genename_2 = genetable_dict_2[s[index_2]]
				homol_dict[genename_1] = genename_2
	return homol_dict

# Leave in editing_site_dict only keys corresponding to sequences that have
# orthologs 

def editing_site_dict_clean(editing_site_dict, homol_dict):
	editing_site_dict_2 = dict()
	for k in editing_site_dict.keys():
		if homol_dict.get(k):
			editing_site_dict_2[k] = editing_site_dict[k]
	return editing_site_dict_2


def build_blast_inp_fasta(fafile_name_1, editing_site_dict, prefix):
	indict = False
	with open(fafile_name_1) as fainp, open("blastinp_{}.fa".format(prefix),'w') as faout:
		for s in fainp:
			if s[0] == '>':
				if editing_site_dict.get(s.strip().split()[0][1:]):
					indict = True
				else:
					indict = False
			if indict:
				faout.write(s)


def tblastx_run(fafile_name_2, prefix):
	system("export PATH=$PATH:~/tools/ncbi-blast-2.7.1+/bin/")
	system("tblastx -query blastinp_{}.fa -db {}  -max_target_seqs 1 -evalue 1e-6 -outfmt 7 -out blastout_{}.txt".format(prefix, fafile_name_2, prefix))


def between(value, a, b):
	return (value >= min(a,b) and value <= max(a,b))

def overlap(alblock, alblock_list):
	for al in alblock_list:
		if between(alblock.q_start, al.q_start, al.q_end):
			return True
		if between(alblock.q_end, al.q_start, al.q_end):
			return True
		if between(al.q_start, alblock.q_start, alblock.q_end):
			return True
		if between(al.q_end, alblock.q_start, alblock.q_end):
			return True
		if between(alblock.s_start, al.s_start, al.s_end):
			return True
		if between(alblock.s_end, al.s_start, al.s_end):
			return True
		if between(al.s_start, alblock.s_start, alblock.s_end):
			return True
		if between(al.s_end, alblock.s_start, alblock.s_end):
			return True
	return False

def blastout_parse(prefix, editing_site_dict, editing_site_dict_2):
	blout_name = "blastout_{}.txt".format(prefix)
	gapped_hits = "gapped_hits_{}.txt".format(prefix)
	align_coord_dict = dict()
	editing_positions = dict()
	with open(blout_name) as blastinp, open(gapped_hits, 'w') as gapped:
		for s in blastinp:
			if s[0] == '#':
				continue
			s_arr = s.strip().split()
			if eval(s_arr[5]) != 0:
				gapped.write(s)
				continue
			edited_positions_1 = []
			for position in editing_site_dict.get(s_arr[0]):
				if position >= eval(s_arr[6]) and position <= eval(s_arr[7]):
					edited_positions_1.append(position)
			edited_positions_2 = []
			if editing_site_dict_2.get(s_arr[1]):
				for position in editing_site_dict_2.get(s_arr[1]):
					if position >= eval(s_arr[8]) and position <= eval(s_arr[9]):
						edited_positions_2.append(position)
			alblock = AlignBlock(
				gene_id1=s_arr[0],
				gene_id2=s_arr[1],
				q_start = eval(s_arr[6]),
				q_end = eval(s_arr[7]),
				s_start = eval(s_arr[8]),
				s_end = eval(s_arr[9]),
				q_sites = edited_positions_1,
				s_sites = edited_positions_2)
			if not align_coord_dict.get((s_arr[0], s_arr[1])):
				align_coord_dict[(s_arr[0], s_arr[1])] = [alblock]
			else:
				if not overlap(alblock, align_coord_dict[(s_arr[0], s_arr[1])]):
					align_coord_dict[(s_arr[0], s_arr[1])].append(alblock)
	return align_coord_dict


def make_fadict(fafile_name):
	with open(fafile_name) as fasta_handle:
		fadict = dict()
		for record in SeqIO.parse(fasta_handle, "fasta"):
			fadict[record.id] = record.seq
	return fadict

def print_aligns(align_coord_dict, fafile_name_1, fafile_name_2, prefix):
	fadict1 = make_fadict(fafile_name_1)
	fadict2 = make_fadict(fafile_name_2)
	with open("align_blocks_{}.txt".format(prefix), 'w') as outfile:
		for (gene_id1, gene_id2) in align_coord_dict.keys():
			for align_block in align_coord_dict[(gene_id1, gene_id2)]:
				align_block.q_seq = fadict1[gene_id1]
				align_block.s_seq = fadict2[gene_id2]
				align_block.print_str(filevar = outfile)

#def main()


parser = OptionParser()
parser.add_option("-a", "--edsites_table_1", help="Editing sites table for organism 1")
parser.add_option("-b", "--edsites_table_2", help="Editing sites table for organism 2")
parser.add_option("-1", "--fasta_path_1", help="FASTA with transcripts for organism 1")
parser.add_option("-2", "--fasta_path_2", help="FASTA with transcripts for organism 2")
parser.add_option("-p", "--proteinortho_file", help="file with proteinortho results")
parser.add_option("-r", "--prefix", help="job prefix")
opt, args = parser.parse_args()

editing_site_dict_raw_1 = editing_site_dict_build(opt.edsites_table_1)
editing_site_dict_raw_2 = editing_site_dict_build(opt.edsites_table_2)
homol_dict_1 = proteinortho_file_parse(opt.proteinortho_file, opt.fasta_path_1, opt.fasta_path_2)
homol_dict_2 = proteinortho_file_parse(opt.proteinortho_file, opt.fasta_path_2, opt.fasta_path_1)
editing_site_dict_1 = editing_site_dict_clean(editing_site_dict_raw_1, homol_dict_1)
editing_site_dict_2 = editing_site_dict_clean(editing_site_dict_raw_2, homol_dict_2)
build_blast_inp_fasta(opt.fasta_path_1, editing_site_dict_1, opt.prefix)
tblastx_run(opt.fasta_path_2, opt.prefix)


align_coord_dict = blastout_parse(opt.prefix, editing_site_dict_1, editing_site_dict_2)
print_aligns(align_coord_dict, opt.fasta_path_1, opt.fasta_path_2, opt.prefix)









