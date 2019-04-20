"""
editing_site_info
"""

from optparse import OptionParser


def mutations_count(positions, seq1, seq2, edited=True):
	A2G = 0
	transversions = 0
	for pos in positions:
		if seq1[pos] == 'A' and seq2[pos] == 'G':
			A2G += 1
		elif edited and (seq1[pos] == 'G' and seq2[pos] == 'G'):
			A2G += 1
		elif edited and (seq1[pos] == 'T' and seq2[pos] == 'C'):
			A2G += 1
		elif edited and (seq1[pos] == 'C' and seq2[pos] == 'C'):
			A2G += 1
		elif edited and ((seq1[pos] in ('A', 'G') and seq2[pos] in ('T', 'C')) or (seq2[pos] in ('A', 'G') and seq1[pos] in ('T', 'C'))):
			transversions += 1
		elif seq1[pos] == 'A' and seq2[pos] in ('T','C'):
			transversions += 1
	return A2G, transversions



parser = OptionParser()
parser.add_option("-i", "--infile", help="count_mutations output")
opt, args = parser.parse_args()

transitions_edited = 0
transversions_edited = 0
transitions_total = 0
transversions_total = 0
total_edited = 0
total_sites = 0
new_edsites = 0

edlet_dict = dict()

with open(opt.infile) as inhandle:
	for alignment in inhandle:
		alignment = alignment.strip().split()
		if alignment[-1] != 'y':
			continue
		if alignment[8] != '*':
			edsites_1 = map(eval, alignment[8].split(','))
		else:
			edsites_1 = []
		if alignment[9] != '*':
			edsites_2 = map(eval, alignment[9].split(','))
		else:
			edsites_2 = []
		specific_edsites = [i for i in edsites_1 if i not in edsites_2]
		new_edsites += len(specific_edsites)
		A2G_all, transversions_all = mutations_count(range(len(alignment[3])), alignment[3], alignment[6], edited=False)
		A2G_ed, transversions_ed = mutations_count(specific_edsites, alignment[3], alignment[6], edited=False)

		transitions_edited += A2G_ed
		transversions_edited += transversions_ed
		transitions_total += A2G_all
		transversions_total += transversions_all
		total_edited += len(edsites_1)
		total_sites += alignment[3].count('A')


print "transitions_edited", transitions_edited
print "transversions_edited", transversions_edited
print "transitions_total", transitions_total
print "transversions_total", transversions_total
print "total_edited", total_edited
print "total_sites", total_sites
print "new_edsites", new_edsites


print "A2G", (float(transitions_edited)/total_edited)/(float(transitions_total)/total_sites)
print "A2OTHER", (float(transversions_edited)/total_edited)/(float(transversions_total)/total_sites)
print "A2G rate", float(transitions_total)/total_sites
print "A2OTHER rate", float(transversions_total)/total_sites


