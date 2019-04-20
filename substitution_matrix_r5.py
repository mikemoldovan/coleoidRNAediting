"""
substitution_matrix_r4

Input:
1. Input multal file
2. Pair of species considered outgroup
3. Pair of species spec1 and spec2, which ancestor is considered.
	Matrices are counted for spec1
4. Confidence intervals for matrix values

Output:
1. subst. matrix mean
2. subst. matrix CI low values
3. subst. matrix CI high values
4. normalized subst. matrix mean
5. normalized subst. matrix CI low values
6. normalized subst. matrix CI high values
7. entrenchment coefficient
"""

from optparse import OptionParser
from Bio.Seq import Seq
from multal_lib import align_parse
from subst_mat_lib import *

#subst_mat {let_ancestral : let_descendal : number of cases}

def make_orf_crd_dict(orf_crd_table):
	orf_crd_dict = dict()
	with open(orf_crd_table) as orf_cd:
		for s in orf_cd:
			if s.startswith('#'):
				continue
			s = s.strip().split()
			orf_crd_dict[s[0]] = [eval(s[1]), eval(s[2])]
	return orf_crd_dict


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

def build_subst_mat(alignment_file, outgr_str, spec_str, select_sites, orf_crd_table):
	outgroup_ref = outgr_str.split(',')
	specs = spec_str.split(',')
	subst_mat = subst_mat_init()
	orf_crd_dict = make_orf_crd_dict(orf_crd_table)
#	control_str = ""
#	count = 0
	for align_obj in align_parse(alignment_file):
		if specs[0] not in align_obj.species_list:
			continue
		if specs[1] not in align_obj.species_list:
			continue
		if (outgroup_ref[0] not in align_obj.species_list) and (outgroup_ref[1] not in align_obj.species_list):
			continue
#		print align_obj.seqinfo_dict 
		seq_id = align_obj.seqinfo_dict[specs[0]].keys()[0]
		orf_crds = orf_crd_dict[seq_id]
		outgroup = [spec_id for spec_id in outgroup_ref if spec_id in align_obj.species_list]
		align_length = len(align_obj.align_dict[specs[0]])

		for i in range(align_length):
			if len(outgroup) == 1:
				outgroup_let = align_obj.align_dict[outgroup[0]][i]
				if align_obj.edinfo_dict[outgroup[0]].get(i+1):
					outgroup_let = codon_change_A_or_E(select_sites, align_obj.edinfo_dict[outgroup[0]][i+1].codon_change)
			elif len(outgroup) == 2:
				let1 = align_obj.align_dict[outgroup[0]][i]
				if align_obj.edinfo_dict[outgroup[0]].get(i+1):
					let1 = codon_change_A_or_E(select_sites, align_obj.edinfo_dict[outgroup[0]][i+1].codon_change)
				let2 = align_obj.align_dict[outgroup[1]][i]
				if align_obj.edinfo_dict[outgroup[1]].get(i+1):
					let2 = codon_change_A_or_E(select_sites, align_obj.edinfo_dict[outgroup[1]][i+1].codon_change)
				outgroup_let = anc_letter(let1, let2)
			if outgroup_let == '-':
				continue

			let1 = align_obj.align_dict[specs[0]][i]
			if align_obj.edinfo_dict[specs[0]].get(i+1):
				let1 = codon_change_A_or_E(select_sites, align_obj.edinfo_dict[specs[0]][i+1].codon_change)
				
			let2 = align_obj.align_dict[specs[1]][i]
			if align_obj.edinfo_dict[specs[1]].get(i+1):
				let2 = codon_change_A_or_E(select_sites, align_obj.edinfo_dict[specs[1]][i+1].codon_change)

			if (let1 == '-') or (let2 == '-'):
				continue

			if outgroup_let == let2:
				if select_sites == "all":
					subst_mat[outgroup_let][let1] += 1
					continue
				if i < orf_crds[0] or i >= orf_crds[1]:
					inorf = False
				else:
					inorf = True

#				if not inorf and select_sites == "syn":
#					subst_mat[outgroup_let][let1] += 1
#					continue
				if not inorf:
					continue

				codon_start = orf_crds[0] + ((i - orf_crds[0])/3)*3
				codon_end = codon_start + 3
				codon = align_obj.align_dict[specs[0]][codon_start:codon_end]
				codon_pos = (i - orf_crds[0]) % 3
				codon = list(codon)
				for j in range(3):
					if align_obj.edinfo_dict[specs[0]].get(codon_start + j + 1):
						codon[j] = 'A'

				codon_anc = codon[:]

				if outgroup_let == 'E':
					codon_anc[codon_pos] = 'A'
				else:
					codon_anc[codon_pos] = outgroup_let

				codon_anc = ''.join(codon_anc)
				codon = ''.join(codon)

				aacid = str(Seq(codon).translate())
				aacid_anc = str(Seq(codon_anc).translate())
				if aacid == aacid_anc and select_sites == "syn":
					subst_mat[outgroup_let][let1] += 1
				if aacid != aacid_anc and select_sites == "nsyn":
					subst_mat[outgroup_let][let1] += 1
#					if let1 == 'E' or outgroup_let == 'E':
#						print codon, codon_anc, aacid, aacid_anc

#	print count
	subst_mat_print(subst_mat)
#	print control_str
	return subst_mat


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-o", "--outgr_str", help="Pair of species ids considered outgroup (comma separated)")
parser.add_option("-s", "--spec_str", help="Pair of species spec1 and spec2, which ancestor is considered (comma separated)")
parser.add_option("-l", "--select_sites", help="Which sites to be selected? 'syn' for synonimous, 'nsyn' for nonsynonimous or 'all' (def.'all')", default="all")
parser.add_option("-c", "--orf_crd_table", help="Table with ORF coordinates (infer_orf_crds.py output)")

opt, args = parser.parse_args()

subst_mat = build_subst_mat(opt.alignment_file, opt.outgr_str, opt.spec_str, opt.select_sites, opt.orf_crd_table)
#subst_mat_n = subst_mat_normalize(subst_mat)
#print "---------------"
#subst_mat_print(subst_mat_n)

python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/sep_versus_all_multal.txt -o oct,bim -s sep,squ -l nsyn -c sep_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/sep_nsyn.txt
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/sep_versus_all_multal.txt -o oct,bim -s sep,squ -l syn  -c sep_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/sep_syn.txt 
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/sep_versus_all_multal.txt -o oct,bim -s sep,squ -l all  -c sep_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/sep_all.txt 

python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/squ_versus_all_multal.txt -o oct,bim -s squ,sep -l nsyn -c squ_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/squ_nsyn.txt
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/squ_versus_all_multal.txt -o oct,bim -s squ,sep -l syn  -c squ_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/squ_syn.txt 
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/squ_versus_all_multal.txt -o oct,bim -s squ,sep -l all  -c squ_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/squ_all.txt 

python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/oct_versus_all_multal.txt -o squ,sep -s oct,bim -l nsyn -c oct_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/oct_nsyn.txt
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/oct_versus_all_multal.txt -o squ,sep -s oct,bim -l syn  -c oct_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/oct_syn.txt 
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/oct_versus_all_multal.txt -o squ,sep -s oct,bim -l all  -c oct_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/oct_all.txt 

python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/bim_versus_all_multal.txt -o squ,sep -s bim,oct -l nsyn -c bim_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/bim_nsyn.txt
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/bim_versus_all_multal.txt -o squ,sep -s bim,oct -l syn  -c bim_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/bim_syn.txt 
python ../scripts/substitution_matrix_r5.py -a ../mutanalysis/bim_versus_all_multal.txt -o squ,sep -s bim,oct -l all  -c bim_orf_crds.txt > ../mutanalysis/substitution_matrices_r5/bim_all.txt 

