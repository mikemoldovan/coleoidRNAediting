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
from multal_lib import align_parse
from subst_mat_lib import *

#subst_mat {let_ancestral : let_descendal : number of cases}


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

def build_subst_mat(alignment_file, confint_p, outgr_str, spec_str, select_sites):
	outgroup_ref = outgr_str.split(',')
	specs = spec_str.split(',')
	subst_mat = subst_mat_init()
#	control_str = ""
#	count = 0
	for align_obj in align_parse(alignment_file):
		if specs[0] not in align_obj.species_list:
			continue
		if specs[1] not in align_obj.species_list:
			continue
		if (outgroup_ref[0] not in align_obj.species_list) and (outgroup_ref[1] not in align_obj.species_list):
			continue
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
				subst_mat[outgroup_let][let1] += 1

#	print count
	subst_mat_print(subst_mat)
#	print control_str
	return subst_mat


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-o", "--outgr_str", help="Pair of species ids considered outgroup (comma separated)")
parser.add_option("-s", "--spec_str", help="Pair of species spec1 and spec2, which ancestor is considered (comma separated)")
parser.add_option("-c", "--confint_p", help="Confidence interval (def. 0.95)", default="0.95")
parser.add_option("-l", "--select_sites", help="Which sites to be selected? 'syn' for synonimous, 'nsyn' for nonsynonimous or 'all' (def.'all')", default="all")

opt, args = parser.parse_args()

subst_mat = build_subst_mat(opt.alignment_file, eval(opt.confint_p), opt.outgr_str, opt.spec_str, opt.select_sites)
#subst_mat_n = subst_mat_normalize(subst_mat)
#print "---------------"
#subst_mat_print(subst_mat_n)
