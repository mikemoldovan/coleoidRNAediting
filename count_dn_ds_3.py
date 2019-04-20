"""
count_dn_ds_3

Count (dN(ed)/dS(ed))/(dN/dS) for a pair of aligned transcriptomes from a multiple alignment file
for different editing levels

* dn/ds is counted separately for mutations to G and mutations to Y
* synonimous substitutions are considered for four and six times redundant codons
* non-synonimous substitutions are considered for all codons

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
from scipy.stats import binom
from Bio.Seq import Seq
from multal_lib import align_parse
from rand_arr_generate import *
from sys import stdout

aminoacids = "VAGPTLSR"

class Outfile_line():
	def __init__(self,
				 edlevel_threshold,
				 total_N_A = 0,
				 total_S_A = 0,
				 total_N_E = 0,
				 total_S_E = 0,
				 dN_e = 0,
				 dN = 0,
				 dS_e = 0,
				 dS = 0):
		self.dN_e = dN_e
		self.dN = dN
		self.dS_e = dS_e
		self.dS = dS
		self.edlevel_threshold = edlevel_threshold
		self.total_N_A = total_N_A
		self.total_S_A = total_S_A
		self.total_N_E = total_N_E
		self.total_S_E = total_S_E

	def p_and_confint_print(self, confint_p, spec_ids, permut_n, dnds_a):
#		n = self.dN_e + self.dS_e
#		if n != 0:
#			p = float(self.dN_e)/n
#		else:
#			p = 'NA'
#			dN_e_low, dN_e_high = binom.interval(confint_p, n, p)
#			dS_e_low = n - dN_e_low
#			dS_e_high = n - dN_e_high

		if self.dS_e > 0:
			dnds_e = (float(self.dN_e)/self.total_N_E)/(float(self.dS_e)/self.total_S_E)
		else:
			dnds_e = 'NA'
#		if dS_e_low > 0:
#			dnds_e_low = (float(dN_e_low)/self.total_N_E)/(float(dS_e_low)/self.total_S_E)
#		else:
#			dnds_e_low = 'NA'
#		if dS_e_high > 0:
#			dnds_e_high = (float(dN_e_high)/self.total_N_E)/(float(dS_e_high)/self.total_S_E)
#		else:
#			dnds_e_high = 'NA'
#		print self.dN, self.total_N_A, self.dS, self.total_S_A
#		dnds_a = (float(self.dN)/self.total_N_A)/(float(self.dS)/self.total_S_A)
		if dnds_e != 'NA':
			dnds_ratio = dnds_e/dnds_a
		else:
			dnds_ratio = 'NA'
#		if dnds_e_low != 'NA':
#			dnds_ratio_low = dnds_e_low/dnds_a
#		else:
#			dnds_ratio_low = 'NA'
#		if dnds_e_high != 'NA':
#			dnds_ratio_high = dnds_e_high/dnds_a
#		else:
#			dnds_ratio_high = 'NA'

		print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(spec_ids,
													self.edlevel_threshold,
													self.dN_e,
													self.dN,
													self.dS_e,
													self.dS,
													self.total_N_E,
													self.total_N_A,
													self.total_S_E,
													self.total_S_A,
													dnds_e,
													dnds_a,
													dnds_ratio) #,
#													dnds_ratio_low,
#													dnds_ratio_high)


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


def main(alignment_file, edlevel_step, confint_p, spec_ids_str, permut_n, orf_crd_table, nucl_list):

	global aminoacids

	nucl_list = list(nucl_list)
	aminoacids = list(aminoacids)
	spec_ids_list = spec_ids_str.split(',')
	edlevel_dict = edlevel_dict_init(edlevel_step)
	edlevels = sorted(edlevel_dict.keys())
	l = len(edlevels)
	orf_crd_dict = make_orf_crd_dict(orf_crd_table)

	for align_obj in align_parse(alignment_file):
		if (spec_ids_list[0] not in align_obj.species_list) or (spec_ids_list[1] not in align_obj.species_list):
			continue
		align_length = len(align_obj.align_dict[spec_ids_list[0]])
		seq_id = align_obj.seqinfo_dict[spec_ids_list[0]].keys()[0]
		orf_crds = orf_crd_dict[seq_id]

		for i in range(orf_crds[0], orf_crds[1]):
			if align_obj.align_dict[spec_ids_list[0]][i] != 'A':
				if not align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
					continue
			if align_obj.align_dict[spec_ids_list[1]][i] == "-":
				continue

			if align_obj.edinfo_dict[spec_ids_list[1]].get(i+1):
				let2 = 'A'
			else:
				let2 = align_obj.align_dict[spec_ids_list[1]][i]

			let1 = 'A'

			codon, aacid, shift = get_codon(align_obj, spec_ids_list[0], orf_crds[0], i)
			
			syn = True

			for nucl in nucl_list:
				codon_new = codon[:]
				codon_new[shift] = nucl
				if aacid != str(Seq(''.join(codon_new)).translate()):
					syn = False

			if syn and aacid not in aminoacids:
				continue

			if align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
				edlevel = align_obj.edinfo_dict[spec_ids_list[0]][i+1].edlevel
				if syn:
					for j in range(1, l):
						if (edlevel < edlevels[j]) and (edlevel >= edlevels[j - 1]):
							edlevel_dict[edlevels[j - 1]].total_S_E += 1
				else:
					for j in range(1, l):
						if (edlevel < edlevels[j]) and (edlevel >= edlevels[j - 1]):
							edlevel_dict[edlevels[j - 1]].total_N_E += 1
			else:
				if syn:
					for j in range(1, l):
						edlevel_dict[edlevels[j - 1]].total_S_A += 1
				else:
					for j in range(1, l):
						edlevel_dict[edlevels[j - 1]].total_N_A += 1

			if let1 == let2:
				continue

			if let2 in nucl_list:
				codon_new = codon[:]
				codon_new[shift] = let2
				aacid_new = str(Seq(''.join(codon_new)).translate())
				if aacid == aacid_new:
					syn = True
				else:
					syn = False

				if align_obj.edinfo_dict[spec_ids_list[0]].get(i+1):
					edlevel = align_obj.edinfo_dict[spec_ids_list[0]][i+1].edlevel
					if syn:
						for j in range(1, l):
							if (edlevel < edlevels[j]) and (edlevel >= edlevels[j - 1]):
								edlevel_dict[edlevels[j - 1]].dS_e += 1
					else:
						for j in range(1, l):
							if (edlevel < edlevels[j]) and (edlevel >= edlevels[j - 1]):
								edlevel_dict[edlevels[j - 1]].dN_e += 1
				else:
					if syn:
						for j in range(1, l):
							edlevel_dict[edlevels[j - 1]].dS += 1
					else:
						for j in range(1, l):
							edlevel_dict[edlevels[j - 1]].dN += 1

	s = edlevel_dict[sorted(edlevel_dict.keys())[0]]
	dnds_a = (float(s.dN)/s.total_N_A)/(float(s.dS)/s.total_S_A)
	for i in sorted(edlevel_dict.keys()):
		edlevel_dict[i].p_and_confint_print(confint_p, spec_ids_str, permut_n, dnds_a)


parser = OptionParser()
parser.add_option("-a", "--alignment_file", help="Input multal file")
parser.add_option("-e", "--edlevel_step", help="Editing level step (in percent)")
parser.add_option("-c", "--confint_p", help="Confidence interval (def. 0.95)", default="0.95")
parser.add_option("-s", "--spec_ids_str", help="Pair of spec ids (comma separated)")
parser.add_option("-n", "--permut_n", help="number of permutations (def. 100000)", default="100000")
parser.add_option("-o", "--orf_crd_table", help="Table with ORF coordinates (infer_orf_crds.py output)")
parser.add_option("-d", "--nucl_list", help="Nucleotides mutations to which to be counted")
opt, args = parser.parse_args()


main(opt.alignment_file, eval(opt.edlevel_step), eval(opt.confint_p), opt.spec_ids_str, eval(opt.permut_n), opt.orf_crd_table, opt.nucl_list)
