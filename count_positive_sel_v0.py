"""
count_positive_sel_v0

Obtain values of the function dN(E_nsyn -> G)/dS(A_syn -> G) with
	combinatorical normalizing

Input:
1. Input multal file
2. Editing level step (in percent)
3. Confidence interval
4. Pair of spec ids (comma separated)

"""

#import scipy.stats as st
from optparse import OptionParser
from Bio.Seq import Seq
from multal_lib import align_parse
from rand_arr_generate import *

class Outfile_line():
	def __init__(self,
				 edlevel_threshold,
				 E_nsyn2G = 0,
				 A_syn2G = 0,
				 n_E_nsyn = 0,
				 n_A_syn = 0,
				 alpha_nsyn = 1,
				 alpha_syn = 1):
		self.edlevel_threshold = edlevel_threshold
		self.E_nsyn2G = E_nsyn2G
		self.A_syn2G = A_syn2G
		self.n_E_nsyn = n_E_nsyn
		self.n_A_syn = n_A_syn
		self.alpha_nsyn = alpha_nsyn
		self.alpha_syn = alpha_syn

	def p_and_confint_print(self, confint_p, spec_ids, permut_n):
		p_E_nsyn2G = float(self.E_nsyn2G)/self.n_E_nsyn
		p_A_syn2G = float(self.A_syn2G)/self.n_A_syn

		midvalue_dnds = (p_E_nsyn2G/p_A_syn2G)*(float(self.alpha_nsyn)/self.alpha_syn)

		rand_arr_p_E_nsyn2G = rand_arr_generate(permut_n, self.n_E_nsyn, p_E_nsyn2G)
		rand_arr_p_A_syn2G = rand_arr_generate(permut_n, self.n_A_syn, p_A_syn2G)

		dnds_arr = []

		for i in range(permut_n):
			p_E_nsyn2G_t = float(rand_arr_p_E_nsyn2G[i])/self.n_E_nsyn
			p_A_syn2G_t = float(rand_arr_p_A_syn2G[i])/self.n_A_syn
			dnds_arr.append(p_E_nsyn2G_t/p_A_syn2G_t)

		dnds_arr = sorted(dnds_arr)
		comb_coeff = float(self.alpha_nsyn)/self.alpha_syn
		dnds_arr = [dnds_arr[i]*comb_coeff for i in range(permut_n)]

		extr_num = int(permut_n*((1-confint_p)/2))

		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(spec_ids,
													self.edlevel_threshold,
													self.E_nsyn2G,
													self.n_E_nsyn,
													self.A_syn2G,
													self.n_A_syn,
													self.alpha_nsyn,
													self.alpha_syn,
													comb_coeff,
													midvalue_dnds,
													dnds_arr[extr_num],
													dnds_arr[-extr_num])


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


def get_codon(align_obj, species, orf_start_crd, crd):
	codon_start = orf_start_crd + ((crd - orf_start_crd)/3)*3
	codon_end = codon_start + 3
	shift = (crd - orf_start_crd) % 3
	codon = align_obj.align_dict[species][codon_start:codon_end]
	codon = list(codon)
	for i in range(3):
		if align_obj.edinfo_dict[species].get(codon_start + i + 1):
			codon[i] = 'A'
	aacid = str(Seq(''.join(codon)).translate())
	return codon, aacid, shift


def syn(codon):
	val = 0
	codon_list = list(codon)
	aacid = str(Seq(codon).translate())
	for i in range(3):
		codon_new = codon_list[:]
		if codon_new[i] == 'A':
			codon_new[i] = 'G'
			aacid_new = str(Seq(''.join(codon_new)).translate())
			if aacid == aacid_new:
				val += 1
	return val


def nsyn(codon):
	val = 0
	codon_list = list(codon)
	aacid = str(Seq(codon).translate())
	for i in range(3):
		codon_new = codon_list[:]
		if codon_new[i] == 'A':
			codon_new[i] = 'G'
			aacid_new = str(Seq(''.join(codon_new)).translate())
			if aacid != aacid_new:
				val += 1
	return val


def alpha_count(codon_usage_dict):
	codon_num = 0
	sum_syn = 0
	sum_nsyn = 0
	for k in codon_usage_dict.keys():
		codon_num += codon_usage_dict[k]
	for k in codon_usage_dict.keys():
		if 'A' in k:
			sum_syn += syn(k)*float(codon_usage_dict[k])/codon_num
			sum_nsyn += nsyn(k)*float(codon_usage_dict[k])/codon_num
	alpha_nsyn = 1/sum_nsyn
	alpha_syn = 1/sum_syn
	return alpha_nsyn, alpha_syn


def main(alignment_file, 
		 edlevel_step, 
		 confint_p, 
		 spec_ids_str, 
		 permut_n,  
		 orf_crd_table):

	spec_ids_list = spec_ids_str.split(',')
	edlevel_dict = edlevel_dict_init(edlevel_step)
	orf_crd_dict = make_orf_crd_dict(orf_crd_table)
	codon_usage_dict = dict()

	for align_obj in align_parse(alignment_file):
		if (spec_ids_list[0] not in align_obj.species_list) or (spec_ids_list[1] not in align_obj.species_list):
			continue

		align_length = len(align_obj.align_dict[spec_ids_list[0]])
		seq_id = align_obj.seqinfo_dict[spec_ids_list[0]].keys()[0]
		orf_crds = orf_crd_dict[seq_id]

		for i in range(orf_crds[0], orf_crds[1]):
			edsite = False

			let_1 = align_obj.align_dict[spec_ids_list[0]][i]
			let_2 = align_obj.align_dict[spec_ids_list[1]][i]

			if let_2 == "-":
				continue
			if let_1 != "A":
				if not align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
					continue
				else:
					let_1 = 'A'

			if align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
				let_1 = 'A'
				edsite = True

			if align_obj.edinfo_dict[spec_ids_list[1]].get(i+1):
				let_2 = 'A'

			if edsite:
				edlevel = align_obj.edinfo_dict[spec_ids_list[0]][i+1].edlevel

				if align_obj.edinfo_dict[spec_ids_list[0]][i+1].codon_change == 'syn':
					continue

				for j in edlevel_dict.keys():
					if j <= edlevel:
						edlevel_dict[j].n_E_nsyn += 1

				if let_2 == 'G':
					for j in edlevel_dict.keys():
						if j <= edlevel:
							edlevel_dict[j].E_nsyn2G += 1

			else:
				if let_2 != 'G' and let_2 != 'A':
					continue

				codon, aacid, shift = get_codon(align_obj, spec_ids_list[0], orf_crds[0], i)
				str_codon = ''.join(codon)
				if codon_usage_dict.get(str_codon):
					codon_usage_dict[str_codon] += 1
				else:
					codon_usage_dict[str_codon] = 1
				codon_new = codon[:]
				codon_new[shift] = let_2
				aacid_new = str(Seq(''.join(codon_new)).translate())

				codon_new_2_G = codon[:]
				codon_new_2_G[shift] = 'G'
				aacid_new_2_G = str(Seq(''.join(codon_new_2_G)).translate())

				if aacid != aacid_new_2_G:
					continue

				for j in edlevel_dict.keys():
					edlevel_dict[j].n_A_syn += 1

				if let_2 == 'G':
					for j in edlevel_dict.keys():
						edlevel_dict[j].A_syn2G += 1

	alpha_nsyn, alpha_syn = alpha_count(codon_usage_dict)

	for i in sorted(edlevel_dict.keys()):
		edlevel_dict[i].alpha_nsyn = alpha_nsyn
		edlevel_dict[i].alpha_syn = alpha_syn

	print "spec_id\tedlevel_thr\tE_nsyn2G\tn_E_nsyn\tA_syn2G\tn_A_syn\talpha_nsyn\talpha_syn\tcomb_coeff\tmidvalue_dnds\tdnds_low\tdnds_high"

	for i in sorted(edlevel_dict.keys()):
		try:
			edlevel_dict[i].p_and_confint_print(confint_p, spec_ids_str, permut_n)
		except ZeroDivisionError:
			pass


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-e", "--edlevel_step", help="Editing level step (in percent, def. 5)", default = '5')
parser.add_option("-s", "--spec_ids_str", help="Pair of spec ids (comma separated)")
parser.add_option("-c", "--confint_p", help="Confidence interval (def. 0.95)", default="0.95")
parser.add_option("-n", "--permut_n", help="number of permutations (def. 100000)", default="100000")
parser.add_option("-o", "--orf_crd_table", help="Table with ORF coordinates (infer_orf_crds.py output)")
opt, args = parser.parse_args()


main(opt.alignment_file, 
	 eval(opt.edlevel_step), 
	 eval(opt.confint_p), 
	 opt.spec_ids_str, 
	 eval(opt.permut_n), 
	 opt.orf_crd_table)

#main(opt.alignment_file, eval(opt.edlevel_step), eval(opt.confint_p), opt.spec_ids_str, eval(opt.permut_n), opt.select_sites)
"""
Apolipoprotein B mRNA-editing enzyme, catalytic polypeptide-like 3
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,oct -o bim_orf_crds.txt > dnds_2_bim_oct.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,bim -o oct_orf_crds.txt > dnds_2_oct_bim.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,squ -o sep_orf_crds.txt > dnds_2_sep_squ.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,sep -o squ_orf_crds.txt > dnds_2_squ_sep.txt

Darwin, C. R. 1872. The origin of species by means of natural selection, or the preservation of favoured races in the struggle for life. London: John Murray. 6th edition; with additions and corrections. Eleventh thousand.

Hedges, S.B., Dudley, J., and Kumar, S. (2006). TimeTree: a public knowledgebase of divergence times among organisms. Bioinformatics 22, 2971–2972.
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,sep -o bim_orf_crds.txt > dnds_2_bim_sep.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,squ -o bim_orf_crds.txt > dnds_2_bim_squ.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,sep -o oct_orf_crds.txt > dnds_2_oct_sep.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,squ -o oct_orf_crds.txt > dnds_2_oct_squ.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,bim -o sep_orf_crds.txt > dnds_2_sep_bim.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,bim -o squ_orf_crds.txt > dnds_2_squ_bim.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,oct -o sep_orf_crds.txt > dnds_2_sep_oct.txt
python ../scripts/count_positive_sel_v0.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,oct -o squ_orf_crds.txt > dnds_2_squ_oct.txt

python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,oct -l nsyn -m nsyn -o bim_orf_crds.txt > bim_oct_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,bim -l nsyn -m nsyn -o oct_orf_crds.txt > oct_bim_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,squ -l nsyn -m nsyn -o sep_orf_crds.txt > sep_squ_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,sep -l nsyn -m nsyn -o squ_orf_crds.txt > squ_sep_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,sep -l nsyn -m nsyn -o bim_orf_crds.txt > bim_sep_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,squ -l nsyn -m nsyn -o bim_orf_crds.txt > bim_squ_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,sep -l nsyn -m nsyn -o oct_orf_crds.txt > oct_sep_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,squ -l nsyn -m nsyn -o oct_orf_crds.txt > oct_squ_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,bim -l nsyn -m nsyn -o sep_orf_crds.txt > sep_bim_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,bim -l nsyn -m nsyn -o squ_orf_crds.txt > squ_bim_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,oct -l nsyn -m nsyn -o sep_orf_crds.txt > sep_oct_nsyn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,oct -l nsyn -m nsyn -o squ_orf_crds.txt > squ_oct_nsyn.txt

python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,oct -l all  -m all  -o bim_orf_crds.txt > bim_oct_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,bim -l all  -m all  -o oct_orf_crds.txt > oct_bim_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,squ -l all  -m all  -o sep_orf_crds.txt > sep_squ_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,sep -l all  -m all  -o squ_orf_crds.txt > squ_sep_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,sep -l all  -m all  -o bim_orf_crds.txt > bim_sep_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,squ -l all  -m all  -o bim_orf_crds.txt > bim_squ_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,sep -l all  -m all  -o oct_orf_crds.txt > oct_sep_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,squ -l all  -m all  -o oct_orf_crds.txt > oct_squ_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,bim -l all  -m all  -o sep_orf_crds.txt > sep_bim_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,bim -l all  -m all  -o squ_orf_crds.txt > squ_bim_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,oct -l all  -m all  -o sep_orf_crds.txt > sep_oct_all.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,oct -l all  -m all  -o squ_orf_crds.txt > squ_oct_all.txt
"""
