"""
substitution_matrix_6

Count substitution matrix for two sister species and an outgroup using
only protein-coding sequences.
Substitution matrix can be counted for three groups of sites depending
on the --flag variable:
--flag=all: substitution matrix is counted for all sites
--flag=syn: only synonimous substitutions are counted
--flag=nsyn: only non-synonimous substitutions are counted

The matrix is counted only for protein-coding regions

Input:
1. Input multal file
2. Pair of species considered outgroup
3. Pair of species spec1 and spec2, which ancestor is considered.
	Matrices are counted for spec1
4. Confidence intervals for matrix values

Output:
substitution matrix 
"""

from optparse import OptionParser
from Bio.Seq import Seq
from multal_lib import align_parse
from subst_mat_lib import *

#subst_mat {let_ancestral : let_descendal : number of substitutions}


def anc_letter(let1, let2):
	if let1 == '-' and let2 == '-':
		return '-'
	if let1 == '-':
		return let2
	if let2 == '-':
		return let1
	if let1 == let2:
		return let1
	else:
		return '-'


def codon_change_A_or_E(select_sites, codon_change):
	if select_sites == "syn" and codon_change == "syn":
		let = 'E'
	elif select_sites == "nsyn" and codon_change != "syn":
		let = 'E'
	elif select_sites == "all":
		let = 'E'
	else:
		let = 'A'
	return let


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


def build_subst_mat(alignment_file,
					outgr_str, 
					spec_str,
					flag,
					orf_crd_table):
	outgroup_ref = outgr_str.split(',')
	specs = spec_str.split(',')
	orf_crd_dict = make_orf_crd_dict(orf_crd_table)

	subst_mat = subst_mat_init()
	n_E = 0
	n_A = 0
	n_G = 0
	n_Y = 0
	n_E_anc = 0
	n_A_anc = 0
	n_G_anc = 0
	n_Y_anc = 0
	n_E2G = 0
	n_A2G = 0
	n_E2Y = 0
	n_A2Y = 0
	n_G2E = 0
	n_G2A = 0
	n_Y2E = 0
	n_Y2A = 0

	for align_obj in align_parse(alignment_file):
		if specs[0] not in align_obj.species_list:
			continue
		if specs[1] not in align_obj.species_list:
			continue
		if (outgroup_ref[0] not in align_obj.species_list) and (outgroup_ref[1] not in align_obj.species_list):
			continue
		outgroup = [spec_id for spec_id in outgroup_ref if spec_id in align_obj.species_list]
		seq_id = align_obj.seqinfo_dict[specs[0]].keys()[0]
		orf_crds = orf_crd_dict[seq_id]
		align_length = len(align_obj.align_dict[specs[0]])

		for i in range(orf_crds[0], orf_crds[1]):
			if len(outgroup) == 1:
				outgroup_let = align_obj.align_dict[outgroup[0]][i]
				if align_obj.edinfo_dict[outgroup[0]].get(i+1):
					outgroup_let = codon_change_A_or_E(flag, align_obj.edinfo_dict[outgroup[0]][i+1].codon_change)
			elif len(outgroup) == 2:
				let1 = align_obj.align_dict[outgroup[0]][i]
				if align_obj.edinfo_dict[outgroup[0]].get(i+1):
					let1 = codon_change_A_or_E(flag, align_obj.edinfo_dict[outgroup[0]][i+1].codon_change)
				let2 = align_obj.align_dict[outgroup[1]][i]
				if align_obj.edinfo_dict[outgroup[1]].get(i+1):
					let2 = codon_change_A_or_E(flag, align_obj.edinfo_dict[outgroup[1]][i+1].codon_change)
				outgroup_let = anc_letter(let1, let2)
			if outgroup_let == '-':
				continue

			let1 = align_obj.align_dict[specs[0]][i]
			dna_let1 = let1
			if align_obj.edinfo_dict[specs[0]].get(i+1):
				let1 = codon_change_A_or_E(flag, align_obj.edinfo_dict[specs[0]][i+1].codon_change)
				dna_let1 = 'A'
				
			let2 = align_obj.align_dict[specs[1]][i]
			if align_obj.edinfo_dict[specs[1]].get(i+1):
				let2 = codon_change_A_or_E(flag, align_obj.edinfo_dict[specs[1]][i+1].codon_change)

			if (let1 == '-') or (let2 == '-'):
				continue

			if outgroup_let != let2:
				continue

			subst_mat[outgroup_let][let1] += 1

			if flag == "all":
				if let1 == 'A':
					n_A += 1
					if outgroup_let == 'G':
						n_G2A += 1
					elif outgroup_let in ('C', 'T'):
						n_Y2A += 1
				elif let1 == 'E':
					n_E += 1
					if outgroup_let == 'G':
						n_G2E += 1
					elif outgroup_let in ('C', 'T'):
						n_Y2E += 1
				elif let1 == 'G':
					n_G += 1
				elif let1 in ('C', 'T'):
					n_Y += 1
				if outgroup_let == 'A':
					n_A_anc += 1
					if let1 == 'G':
						n_A2G += 1
					elif let1 in ('C', 'T'):
						n_A2Y += 1
				elif outgroup_let == 'E':
					n_E_anc += 1
					if let1 == 'G':
						n_E2G += 1
					elif let1 in ('C', 'T'):
						n_E2Y += 1
				elif outgroup_let == 'G':
					n_G_anc += 1
				elif outgroup_let in ('C', 'T'):
					n_Y_anc += 1

			elif flag == "syn":
				codon, aacid, shift = get_codon(align_obj, specs[0], orf_crds[0], i)
				if let1 == 'A' or let1 == 'E':
					codon_new = codon[:]
					codon_new[shift] = 'G'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					if aacid != aacid_new:
						continue
					if let1 == 'A':
						n_A += 1
						if outgroup_let == 'G':
							n_G2A += 1
						elif outgroup_let in ('C', 'T'):
							n_Y2A += 1
					elif let1 == 'E':
						n_E += 1
						if outgroup_let == 'G':
							n_G2E += 1
						elif outgroup_let in ('C', 'T'):
							n_Y2E += 1

				elif let1 in ('C', 'T'):
					codon[shift] = 'A'
					aacid = str(Seq(''.join(codon)).translate())
					codon_new = codon[:]
					codon_new[shift] = 'G'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					if aacid == aacid_new:
						n_Y += 1
						if outgroup_let == 'E':
							n_E2Y += 1
						elif outgroup_let == 'A':
							n_A2Y += 1

				elif let1 == 'G':
					codon_new = codon[:]
					codon_new[shift] = 'A'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					if aacid == aacid_new:
						n_G += 1
						if outgroup_let == 'E':
							n_E2G += 1
						elif outgroup_let == 'A':
							n_A2G += 1

				if outgroup_let == 'E':
					n_E_anc += 1
				elif outgroup_let == 'A':
					codon_new = codon[:]
					codon_new[shift] = 'A'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					codon_new_G = codon[:]
					codon_new_G[shift] = 'G'
					aacid_new_G = str(Seq(''.join(codon_new)).translate())
					if aacid_new == aacid_new_G:
						n_A_anc += 1
				elif outgroup_let == 'G':
					codon_new = codon[:]
					codon_new[shift] = 'A'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					codon_new_G = codon[:]
					codon_new_G[shift] = 'G'
					aacid_new_G = str(Seq(''.join(codon_new)).translate())
					if aacid_new == aacid_new_G:
						n_G_anc += 1
				elif outgroup_let in ('C','T'):
					codon_new = codon[:]
					codon_new[shift] = 'A'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					codon_new_G = codon[:]
					codon_new_G[shift] = 'G'
					aacid_new_G = str(Seq(''.join(codon_new)).translate())
					if aacid_new == aacid_new_G:
						n_Y += 1


			elif flag == "nsyn":
				codon, aacid, shift = get_codon(align_obj, specs[0], orf_crds[0], i)
				if let1 == 'A' or let1 == 'E':
					codon_new = codon[:]
					codon_new[shift] = 'G'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					if aacid == aacid_new:
						continue
					if let1 == 'A':
						n_A += 1
						if outgroup_let == 'G':
							n_G2A += 1
						elif outgroup_let in ('C', 'T'):
							n_Y2A += 1
					elif let1 == 'E':
						n_E += 1
						if outgroup_let == 'G':
							n_G2E += 1
						elif outgroup_let in ('C', 'T'):
							n_Y2E += 1

				elif let1 in ('C', 'T'):
					codon[shift] = 'A'
					aacid = str(Seq(''.join(codon)).translate())
					codon_new = codon[:]
					codon_new[shift] = 'G'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					if aacid != aacid_new:
						n_Y += 1
						if outgroup_let == 'E':
							n_E2Y += 1
						elif outgroup_let == 'A':
							n_A2Y += 1

				elif let1 == 'G':
					codon_new = codon[:]
					codon_new[shift] = 'A'
					aacid_new = str(Seq(''.join(codon_new)).translate())
					if aacid != aacid_new:
						n_G += 1
						if outgroup_let == 'E':
							n_E2G += 1
						elif outgroup_let == 'A':
							n_A2G += 1

	print "@\tn_E\tn_A\tn_G\tn_Y\tn_E2G\tn_A2G\tn_E2Y\tn_A2Y\tn_G2E\tn_G2A\tn_Y2E\tn_Y2A"
	print '*', n_E, n_A, n_G, n_Y, n_E2G, n_A2G, n_E2Y, n_A2Y, n_G2E, n_G2A, n_Y2E, n_Y2A
	subst_mat_print(subst_mat)
#	print control_str
	return subst_mat


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-u", "--outgr_str", help="Pair of species ids considered outgroup (comma separated)")
parser.add_option("-s", "--spec_str", help="Pair of species spec1 and spec2, which ancestor is considered (comma separated)")
parser.add_option("-f", "--flag", help="Which sites to be selected? 'syn' for synonimous, 'nsyn' for nonsynonimous or 'all' (def.'all')", default="all")
parser.add_option("-o", "--orf_crd_table", help="Table with ORF coordinates (infer_orf_crds.py output)")
opt, args = parser.parse_args()

subst_mat = build_subst_mat(opt.alignment_file,
							opt.outgr_str, 
							opt.spec_str,
							opt.flag,
							opt.orf_crd_table)
#subst_mat_n = subst_mat_normalize(subst_mat)
#print "---------------"
#subst_mat_print(subst_mat_n)

#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/sep_versus_all_multal.txt -u oct,bim -s sep,squ -f nsyn -o sep_orf_crds.txt > ../substmats_new/sep_nsyn.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/squ_versus_all_multal.txt -u oct,bim -s squ,sep -f nsyn -o squ_orf_crds.txt > ../substmats_new/squ_nsyn.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/oct_versus_all_multal.txt -u sep,squ -s oct,bim -f nsyn -o oct_orf_crds.txt > ../substmats_new/oct_nsyn.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/bim_versus_all_multal.txt -u sep,squ -s bim,oct -f nsyn -o bim_orf_crds.txt > ../substmats_new/bim_nsyn.txt
#
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/sep_versus_all_multal.txt -u oct,bim -s sep,squ -f syn -o sep_orf_crds.txt > ../substmats_new/sep_syn.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/squ_versus_all_multal.txt -u oct,bim -s squ,sep -f syn -o squ_orf_crds.txt > ../substmats_new/squ_syn.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/oct_versus_all_multal.txt -u sep,squ -s oct,bim -f syn -o oct_orf_crds.txt > ../substmats_new/oct_syn.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/bim_versus_all_multal.txt -u sep,squ -s bim,oct -f syn -o bim_orf_crds.txt > ../substmats_new/bim_syn.txt
#
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/sep_versus_all_multal.txt -u oct,bim -s sep,squ -f all -o sep_orf_crds.txt > ../substmats_new/sep_all.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/squ_versus_all_multal.txt -u oct,bim -s squ,sep -f all -o squ_orf_crds.txt > ../substmats_new/squ_all.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/oct_versus_all_multal.txt -u sep,squ -s oct,bim -f all -o oct_orf_crds.txt > ../substmats_new/oct_all.txt
#python ../scripts/substitution_matrix_6.py -a ../mutanalysis/bim_versus_all_multal.txt -u sep,squ -s bim,oct -f all -o bim_orf_crds.txt > ../substmats_new/bim_all.txt
