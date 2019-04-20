"""
double_statcount_r0

Obtain stats for a pair of aligned transcriptomes from a multiple alignment file
for different editing levels

Input:
1. Input multal file
2. Editing level step (in percent)
3. Confidence interval
4. Pair of spec ids (comma ceparated)

"""

#import scipy.stats as st
from optparse import OptionParser
from multal_lib import align_parse
from rand_arr_generate import *

class Outfile_line():
	def __init__(self,
				 edlevel_threshold,
				 A_ed2G=0,
				 A2G=0,
				 A_ed2CT=0,
				 A2CT=0,
				 n_A=0,
				 n_A_ed=0):
		self.A_ed2G = A_ed2G
		self.A2G = A2G
		self.A_ed2CT = A_ed2CT
		self.A2CT = A2CT
		self.n_A = n_A
		self.n_A_ed = n_A_ed
		self.edlevel_threshold = edlevel_threshold

	def p_and_confint_print(self, confint_p, spec_ids, permut_n):
		p_A2G = float(self.A2G)/self.n_A
		p_Aed2G = float(self.A_ed2G)/self.n_A_ed
		p_A2CT = float(self.A2CT)/self.n_A
		p_A_ed2CT = float(self.A_ed2CT)/self.n_A_ed

		rand_arr_p_A2G = rand_arr_generate(permut_n, self.n_A, p_A2G)
#		print p_A2G, rand_arr_p_A2G
		rand_arr_p_A2G = [float(rand_arr_p_A2G[i])/self.n_A for i in range(permut_n)]
#		print rand_arr_p_A2G
		rand_arr_p_Aed2G = rand_arr_generate(permut_n, self.n_A_ed, p_Aed2G)
#		print rand_arr_p_Aed2G
		rand_arr_p_Aed2G = [float(rand_arr_p_Aed2G[i])/self.n_A_ed for i in range(permut_n)]
#		print rand_arr_p_Aed2G
		rf_G_arr = sorted([float(rand_arr_p_Aed2G[i])/rand_arr_p_A2G[i] for i in range(permut_n)])
#		print rf_G_arr

		rand_arr_p_A2CT = rand_arr_generate(permut_n, self.n_A, p_A2CT)
		rand_arr_p_A2CT = [float(rand_arr_p_A2CT[i])/self.n_A for i in range(permut_n)]
		rand_arr_p_Aed2CT = rand_arr_generate(permut_n, self.n_A_ed, p_A_ed2CT)
		rand_arr_p_Aed2CT = [float(rand_arr_p_Aed2CT[i])/self.n_A_ed for i in range(permut_n)]
		rf_CT_arr = sorted([float(rand_arr_p_Aed2CT[i])/rand_arr_p_A2CT[i] for i in range(permut_n)])

		extr_num = int(permut_n*((1-confint_p)/2))

		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(spec_ids,
													self.edlevel_threshold,
													p_Aed2G,
													p_A2G,
													p_A_ed2CT,
													p_A2CT,
													p_Aed2G/p_A2G,
													rf_G_arr[extr_num],
													rf_G_arr[-extr_num],
													p_A_ed2CT/p_A2CT,
													rf_CT_arr[extr_num],
													rf_CT_arr[-extr_num])


def edlevel_dict_init(edlevel_step):
	edlevel_dict = dict()
	for i in range(0, 100, edlevel_step):
		edlevel_dict[i] = Outfile_line(edlevel_threshold = i)
	return edlevel_dict


def main(alignment_file, edlevel_step, confint_p, spec_ids_str, permut_n, select_sites):
	spec_ids_list = spec_ids_str.split(',')
	edlevel_dict = edlevel_dict_init(edlevel_step)
	for align_obj in align_parse(alignment_file):
		if (spec_ids_list[0] not in align_obj.species_list) or (spec_ids_list[1] not in align_obj.species_list):
			continue
		align_length = len(align_obj.align_dict[spec_ids_list[0]])
		for i in range(align_length):
			if align_obj.align_dict[spec_ids_list[0]][i] != 'A':
				continue
			if align_obj.align_dict[spec_ids_list[1]][i] == "-":
				continue

			for j in edlevel_dict.keys():
				edlevel_dict[j].n_A += 1

			if align_obj.align_dict[spec_ids_list[1]][i] == 'G':
				for j in edlevel_dict.keys():
					edlevel_dict[j].A2G += 1
			elif align_obj.align_dict[spec_ids_list[1]][i] in ('C','T'):
				for j in edlevel_dict.keys():
					edlevel_dict[j].A2CT += 1

			if not align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
				continue
			if select_sites == 'syn' and align_obj.edinfo_dict[spec_ids_list[0]][i+1].codon_change != 'syn':
				continue
			if select_sites == 'nsyn' and align_obj.edinfo_dict[spec_ids_list[0]][i+1].codon_change == 'syn':
				continue

			edlevel = align_obj.edinfo_dict[spec_ids_list[0]][i+1].edlevel
			for j in edlevel_dict.keys():
				if j < edlevel:
					edlevel_dict[j].n_A_ed += 1

			if align_obj.edinfo_dict[spec_ids_list[1]].get(i+1):
				continue

			if align_obj.align_dict[spec_ids_list[1]][i] == 'G':
				for j in edlevel_dict.keys():
					if j < edlevel:
						edlevel_dict[j].A_ed2G += 1
			elif align_obj.align_dict[spec_ids_list[1]][i] in ('C','T'):
				for j in edlevel_dict.keys():
					if j < edlevel:
						edlevel_dict[j].A_ed2CT += 1

	print "spec_id\tedlevel_thr\tp_Aed2G\tp_A2G\tp_A_ed2CT\tp_A2CT\trf_G\trf_G_low\trf_G_high\trf_CT\trf_CT_low\trf_CT_high"

	for i in sorted(edlevel_dict.keys()):
		try:
			edlevel_dict[i].p_and_confint_print(confint_p, spec_ids_str, permut_n)
		except ZeroDivisionError:
			pass


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-e", "--edlevel_step", help="Editing level step (in percent)")
parser.add_option("-s", "--spec_ids_str", help="Pair of spec ids (comma ceparated)")
parser.add_option("-c", "--confint_p", help="Confidence interval (def. 0.95)", default="0.95")
parser.add_option("-n", "--permut_n", help="number of permutations (def. 100000)", default="100000")
parser.add_option("-l", "--select_sites", help="Which sites to be selected? 'syn' for synonimous, 'nsyn' for nonsynonimous or 'all' (def.'all')", default="all")
opt, args = parser.parse_args()


main(opt.alignment_file, eval(opt.edlevel_step), eval(opt.confint_p), opt.spec_ids_str, eval(opt.permut_n), opt.select_sites)

