"""
count_dn_ds

Count (dN(ed)/dS(ed))/(dN/dS) for a pair of aligned transcriptomes from a multiple alignment file
for different editing levels

Input:
1. Input multal file
2. Editing level step (in percent)
3. Confidence interval
4. Pair of spec ids (comma separated)
5. Nucleotides for which dNdS to be counted
6. Table with ORF coordinates (infer_orf_crds.py output)

"""

#import scipy.stats as st
from optparse import OptionParser
from Bio.Seq import Seq
from multal_lib import align_parse
from rand_arr_generate import *

class Outfile_line():
	def __init__(self,
				 edlevel_threshold,
				 dN_e = 0,
				 dN = 0,
				 dS_e = 0,
				 dS = 0):
		self.dN_e = dN_e
		self.dN = dN
		self.dS_e = dS_e
		self.dS = dS
		self.edlevel_threshold = edlevel_threshold

	def p_and_confint_print(self, confint_p, spec_ids, permut_n):
		edsite_mutations = self.dN_e + self.dS_e
		non_edsite_mutations = self.dN + self.dS
		p_nsyn_mut_e = float(self.dN_e)/edsite_mutations
		p_nsyn_mut = float(self.dN)/non_edsite_mutations

		rand_arr_p_nsyn_mut_e = rand_arr_generate(permut_n, edsite_mutations, p_nsyn_mut_e)
		rand_arr_p_syn_mut_e = rand_arr_generate(permut_n, edsite_mutations, 1 - p_nsyn_mut_e)
		rand_arr_p_nsyn_mut = rand_arr_generate(permut_n, edsite_mutations, p_nsyn_mut)
		rand_arr_p_syn_mut = rand_arr_generate(permut_n, edsite_mutations, 1 - p_nsyn_mut)

		dN_dS_n_arr = []
		for i in range(permut_n):
			dnds_e = float(rand_arr_p_nsyn_mut_e[i])/rand_arr_p_syn_mut_e[i]
			dnds = float(rand_arr_p_nsyn_mut[i])/rand_arr_p_syn_mut[i]
			dN_dS_n_arr.append(dnds_e/dnds)
		dN_dS_n_arr = sorted(dN_dS_n_arr)

		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(spec_ids,
													self.edlevel_threshold,
													self.dN_e,
													self.dN,
													self.dS_e,
													self.dS,
													float(self.dN_e)/self.dS_e,
													float(self.dN)/self.dS,
													(float(self.dN_e)/self.dS_e)/(float(self.dN)/self.dS),
													dN_dS_n_arr[int(permut_n*(1 - confint_p)/2)],
													dN_dS_n_arr[int(permut_n*(1 - (1 - confint_p)/2))])

#"spec_id\tedlevel_threshold\tdN_e\tdN\tdS_e\tdS\tdN_e_dS_e\tdN_dS\tdN_dS_n"

def edlevel_dict_init(edlevel_step):
	edlevel_dict = dict()
	for i in range(0, 100, edlevel_step):
		edlevel_dict[i] = Outfile_line(edlevel_threshold = i)
	return edlevel_dict


def make_orf_crd_dict(orf_crd_table):
	orf_crd_dict = dict()
	with open(orf_crd_table) as orf_cd:
		for s in orf_cd:
			if s.startswith('#'):
				continue
			s = s.strip().split()
			orf_crd_dict[s[0]] = [eval(s[1]), eval(s[2])]
	return orf_crd_dict


def main(alignment_file, edlevel_step, confint_p, spec_ids_str, permut_n, orf_crd_table, nucl_list):
	nucl_list = list(nucl_list)
	spec_ids_list = spec_ids_str.split(',')
	edlevel_dict = edlevel_dict_init(edlevel_step)
	orf_crd_dict = make_orf_crd_dict(orf_crd_table)

	for align_obj in align_parse(alignment_file):
		if (spec_ids_list[0] not in align_obj.species_list) or (spec_ids_list[1] not in align_obj.species_list):
			continue
		align_length = len(align_obj.align_dict[spec_ids_list[0]])
		seq_id = align_obj.seqinfo_dict[spec_ids_list[0]].keys()[0]
		orf_crds = orf_crd_dict[seq_id]
#		print orf_crds, len(align_obj.align_dict[spec_ids_list[0]])
		for i in range(orf_crds[0], orf_crds[1]):
			if align_obj.align_dict[spec_ids_list[0]][i] != 'A':
				if not align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
					continue
			if align_obj.align_dict[spec_ids_list[1]][i] == "-":
				continue

			if align_obj.edinfo_dict[spec_ids_list[1]].get(i+1):
				continue

			if align_obj.align_dict[spec_ids_list[1]][i] in nucl_list:
				codon_start = orf_crds[0] + ((i - orf_crds[0])/3)*3
				codon_end = codon_start + 3
				shift = (i - orf_crds[0]) % 3
				codon = align_obj.align_dict[spec_ids_list[0]][codon_start:codon_end]
				codon = list(codon)
				
				for j in range(3):
					if align_obj.edinfo_dict[spec_ids_list[0]].get(codon_start + j + 1):
						codon[j] = 'A'
				
				codon_subst = codon[:]
				codon_subst[shift] = align_obj.align_dict[spec_ids_list[1]][i]
				codon = ''.join(codon)
				codon_subst = ''.join(codon_subst)
				aacid = str(Seq(codon).translate())
				aacid_subst = str(Seq(codon_subst).translate())
				
				if align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
					edlevel = align_obj.edinfo_dict[spec_ids_list[0]][i+1].edlevel
					l = len(edlevel_dict.keys())
					edlevels = sorted(edlevel_dict.keys())
					for j in range(1,l):
						if (edlevel < edlevels[j]) and (edlevel >= edlevels[j - 1]):
							if aacid == aacid_subst:
								edlevel_dict[edlevels[j - 1]].dS_e += 1
							else:
								edlevel_dict[edlevels[j - 1]].dN_e += 1
				else:
					for j in edlevel_dict.keys():
						if aacid == aacid_subst:
							edlevel_dict[j].dS += 1
						else:
							edlevel_dict[j].dN += 1

#	print "spec_id\tedlevel_threshold\tdN_e\tdN\tdS_e\tdS\tdN_e_dS_e\dN_dS\tdN_dS_n\tdN_dS_n_low\tdN_dS_n_high"

	for i in sorted(edlevel_dict.keys()):
		try:
			edlevel_dict[i].p_and_confint_print(confint_p, spec_ids_str, permut_n)
		except ZeroDivisionError:
			pass


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-e", "--edlevel_step", help="Editing level step (in percent)")
parser.add_option("-c", "--confint_p", help="Confidence interval (def. 0.95)", default="0.95")
parser.add_option("-s", "--spec_ids_str", help="Pair of spec ids (comma ceparated)")
parser.add_option("-n", "--permut_n", help="number of permutations (def. 100000)", default="100000")
parser.add_option("-o", "--orf_crd_table", help="Table with ORF coordinates (infer_orf_crds.py output)")
parser.add_option("-d", "--nucl_list", help="Nucleotides mutations to ehich to be counted")
opt, args = parser.parse_args()


main(opt.alignment_file, eval(opt.edlevel_step), eval(opt.confint_p), opt.spec_ids_str, eval(opt.permut_n), opt.orf_crd_table, opt.nucl_list)

