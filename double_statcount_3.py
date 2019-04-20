"""
double_statcount_r0

Obtain stats for a pair of aligned transcriptomes from a multiple alignment file
for different editing levels for different normalizing coefficients

Input:
1. Input multal file
2. Editing level step (in percent)
3. Confidence interval
4. Pair of spec ids (comma ceparated)

"""

#import scipy.stats as st
from optparse import OptionParser
from Bio.Seq import Seq
from multal_lib import align_parse
from rand_arr_generate import *

class Outfile_line():
	def __init__(self,
				 edlevel_threshold,
				 E2G=0,
				 A2G=0,
				 E2Y=0,
				 A2Y=0,
				 n_AG=0,
				 n_AY=0,
				 n_E=0):
		self.E2G = E2G
		self.A2G = A2G
		self.E2Y = E2Y
		self.A2Y = A2Y
		self.n_AG = n_AG
		self.n_AY = n_AY
		self.n_E = n_E
		self.edlevel_threshold = edlevel_threshold

	def p_and_confint_print(self, confint_p, spec_ids, permut_n, syn_nsyn_asites):
#		print self.A2G, self.n_A
		p_E2G = float(self.E2G)/self.n_E
		p_E2Y = float(self.E2Y)/self.n_E
		if syn_nsyn_asites == "all":
			p_A2G = float(self.A2G)/self.n_AG
			p_A2Y = float(self.A2Y)/self.n_AG
		else:
			p_A2G = float(self.A2G)/self.n_AG
			if syn_nsyn_asites == "syn":
				p_A2Y = float(self.A2Y)/self.n_AY
			elif syn_nsyn_asites == "nsyn":
				p_A2Y = float(self.A2Y)/self.n_AY

		rand_arr_p_A2G = rand_arr_generate(permut_n, self.n_AG, p_A2G)
		rand_arr_p_A2G = [float(rand_arr_p_A2G[i])/self.n_AG for i in range(permut_n)]
		rand_arr_p_E2G = rand_arr_generate(permut_n, self.n_E, p_E2G)
		rand_arr_p_E2G = [float(rand_arr_p_E2G[i])/self.n_E for i in range(permut_n)]
		rf_G_arr = sorted([float(rand_arr_p_E2G[i])/rand_arr_p_A2G[i] for i in range(permut_n)])

		if syn_nsyn_asites == "all":
			rand_arr_p_A2Y = rand_arr_generate(permut_n, self.n_AG, p_A2Y)
			rand_arr_p_A2Y = [float(rand_arr_p_A2Y[i])/self.n_AG for i in range(permut_n)]
		else:
			rand_arr_p_A2Y = rand_arr_generate(permut_n, self.n_AY, p_A2Y)
			if syn_nsyn_asites == "syn":
				rand_arr_p_A2Y = [float(rand_arr_p_A2Y[i])/self.n_AY for i in range(permut_n)]
			elif syn_nsyn_asites == "nsyn":
				rand_arr_p_A2Y = [float(rand_arr_p_A2Y[i])/self.n_AY for i in range(permut_n)]
		rand_arr_p_E2Y = rand_arr_generate(permut_n, self.n_E, p_E2Y)
		rand_arr_p_E2Y = [float(rand_arr_p_E2Y[i])/self.n_E for i in range(permut_n)]
		rf_Y_arr = sorted([float(rand_arr_p_E2Y[i])/rand_arr_p_A2Y[i] for i in range(permut_n)])

		extr_num = int(permut_n*((1-confint_p)/2))

		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(spec_ids,
													self.edlevel_threshold,
													p_E2G,
													p_A2G,
													p_E2Y,
													p_A2Y,
													p_E2G/p_A2G,
													rf_G_arr[extr_num],
													rf_G_arr[-extr_num],
													p_E2Y/p_A2Y,
													rf_Y_arr[extr_num],
													rf_Y_arr[-extr_num])


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


def main(alignment_file, 
		 edlevel_step, 
		 confint_p, 
		 spec_ids_str, 
		 permut_n, 
		 syn_nsyn_esites, 
		 syn_nsyn_asites, 
		 orf_crd_table):

	spec_ids_list = spec_ids_str.split(',')
	edlevel_dict = edlevel_dict_init(edlevel_step)
	orf_crd_dict = make_orf_crd_dict(orf_crd_table)

	for align_obj in align_parse(alignment_file):
		if (spec_ids_list[0] not in align_obj.species_list) or (spec_ids_list[1] not in align_obj.species_list):
			continue

#		print spec_ids_list

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

				if syn_nsyn_esites == 'syn' and align_obj.edinfo_dict[spec_ids_list[0]][i+1].codon_change != 'syn':
					continue
				if syn_nsyn_esites == 'nsyn' and align_obj.edinfo_dict[spec_ids_list[0]][i+1].codon_change == 'syn':
					continue

				for j in edlevel_dict.keys():
					if j <= edlevel:
						edlevel_dict[j].n_E += 1

				if let_2 == 'G':
					for j in edlevel_dict.keys():
						if j <= edlevel:
							edlevel_dict[j].E2G += 1
				elif let_2 in ('C','T'):
					for j in edlevel_dict.keys():
						if j <= edlevel:
							edlevel_dict[j].E2Y += 1

			else:
				codon, aacid, shift = get_codon(align_obj, spec_ids_list[0], orf_crds[0], i)
				codon_new = codon[:]
				codon_new[shift] = let_2
				aacid_new = str(Seq(''.join(codon_new)).translate())

				codon_new_2_G = codon[:]
				codon_new_2_G[shift] = 'G'
				aacid_new_2_G = str(Seq(''.join(codon_new_2_G)).translate())
#				codon_new_2_C = codon[:]
#				codon_new_2_C[shift] = 'C'
#				aacid_new_2_C = str(Seq(''.join(codon_new_2_C)).translate())
#				codon_new_2_T = codon[:]
#				codon_new_2_T[shift] = 'T'
#				aacid_new_2_T = str(Seq(''.join(codon_new_2_T)).translate())
				if aacid == aacid_new_2_G:
					syn = True
				else:
					syn = False

				if syn_nsyn_asites == 'syn' and not syn:
					continue
				elif syn_nsyn_asites == 'nsyn' and syn:
					continue
				elif syn and syn_nsyn_asites == 'syn':
					for j in edlevel_dict.keys():
						edlevel_dict[j].n_AG += 1
						edlevel_dict[j].n_AY += 1
				elif not syn and syn_nsyn_asites == 'nsyn':
					for j in edlevel_dict.keys():
						edlevel_dict[j].n_AG += 1
						edlevel_dict[j].n_AY += 1
				elif syn_nsyn_asites == 'all':
					for j in edlevel_dict.keys():
						edlevel_dict[j].n_AG += 1

#				if aacid == aacid_new_2_C and syn_nsyn_asites == 'syn':
#					if aacid == aacid_new_2_G:
#						for j in edlevel_dict.keys():
#							edlevel_dict[j].n_AY += 1
#				elif aacid != aacid_new_2_C and syn_nsyn_asites == 'nsyn':
#					if aacid != aacid_new_2_G:
#						for j in edlevel_dict.keys():
#							edlevel_dict[j].n_AY += 1
#
#				if aacid == aacid_new_2_T and syn_nsyn_asites == 'syn':
#					if aacid == aacid_new_2_G:
#						for j in edlevel_dict.keys():
#							edlevel_dict[j].n_AY += 1
#				elif aacid != aacid_new_2_T and syn_nsyn_asites == 'nsyn':
#					if aacid != aacid_new_2_G:
#						for j in edlevel_dict.keys():
#							edlevel_dict[j].n_AY += 1

				if let_2 == 'G':
					for j in edlevel_dict.keys():
						edlevel_dict[j].A2G += 1
				elif let_2 in ('C','T'):
					for j in edlevel_dict.keys():
						edlevel_dict[j].A2Y += 1

	print "spec_id\tedlevel_thr\tp_E2G\tp_A2G\tp_E2Y\tp_A2Y\trf_G\trf_G_low\trf_G_high\trf_Y\trf_Y_low\trf_Y_high"

	for i in sorted(edlevel_dict.keys()):
		try:
			edlevel_dict[i].p_and_confint_print(confint_p, spec_ids_str, permut_n, syn_nsyn_asites)
		except ZeroDivisionError:
			pass


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-e", "--edlevel_step", help="Editing level step (in percent, def. 5)", default = '5')
parser.add_option("-s", "--spec_ids_str", help="Pair of spec ids (comma separated)")
parser.add_option("-c", "--confint_p", help="Confidence interval (def. 0.95)", default="0.95")
parser.add_option("-n", "--permut_n", help="number of permutations (def. 100000)", default="100000")
parser.add_option("-l", "--syn_nsyn_esites", help="""Which editing sites to be selected? 
	'syn' for synonimous, 'nsyn' for nonsynonimous or 'all' (def.'all')""", default="all")
parser.add_option("-m", "--syn_nsyn_asites", help="""Which kind of substitutions to be selected? 
	'syn' for synonimous, 'nsyn' for nonsynonimous or 'all' (def.'all')""", default="all")
parser.add_option("-o", "--orf_crd_table", help="Table with ORF coordinates (infer_orf_crds.py output)")
opt, args = parser.parse_args()


main(opt.alignment_file, 
	 eval(opt.edlevel_step), 
	 eval(opt.confint_p), 
	 opt.spec_ids_str, 
	 eval(opt.permut_n), 
	 opt.syn_nsyn_esites, 
	 opt.syn_nsyn_asites, 
	 opt.orf_crd_table)

#main(opt.alignment_file, eval(opt.edlevel_step), eval(opt.confint_p), opt.spec_ids_str, eval(opt.permut_n), opt.select_sites)
"""
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,oct -l syn  -m syn  -o bim_orf_crds.txt > bim_oct_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,bim -l syn  -m syn  -o oct_orf_crds.txt > oct_bim_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,squ -l syn  -m syn  -o sep_orf_crds.txt > sep_squ_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,sep -l syn  -m syn  -o squ_orf_crds.txt > squ_sep_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,sep -l syn  -m syn  -o bim_orf_crds.txt > bim_sep_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/bim_versus_all_multal.txt -s bim,squ -l syn  -m syn  -o bim_orf_crds.txt > bim_squ_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,sep -l syn  -m syn  -o oct_orf_crds.txt > oct_sep_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/oct_versus_all_multal.txt -s oct,squ -l syn  -m syn  -o oct_orf_crds.txt > oct_squ_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,bim -l syn  -m syn  -o sep_orf_crds.txt > sep_bim_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,bim -l syn  -m syn  -o squ_orf_crds.txt > squ_bim_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/sep_versus_all_multal.txt -s sep,oct -l syn  -m syn  -o sep_orf_crds.txt > sep_oct_syn.txt
python ../scripts/double_statcount_3.py -a ../mutanalysis/squ_versus_all_multal.txt -s squ,oct -l syn  -m syn  -o squ_orf_crds.txt > squ_oct_syn.txt

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
