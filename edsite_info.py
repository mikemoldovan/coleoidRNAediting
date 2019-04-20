"""
editing_site_info
"""

from optparse import OptionParser


def new_edsite_num(edsite_string_1, edsite_string_2):
	if edsite_string_1 == '*':
		return 0
	edsite_string_1 = edsite_string_1.split(',')
	if edsite_string_2 == '*':
		return len(edsite)
	edsite_string_2 = edsite_string_2.split(',')
	count = 0
	for i in edsite_string_1:
		if i not in edsite_string_2:
			count += 1
	return count

parser = OptionParser()
parser.add_option("-i", "--infile", help="count_mutations output")
opt, args = parser.parse_args()


transitions_edited = 0
transversions_edited = 0
transitions_total = 0
transversions_total = 0
total_edited = 0
total_sites = 0
became_unedited = 0

edlet_dict = dict()

with open(opt.infile) as infile:
	for s in infile:
		s = s.strip().split()
		seq1 = s[3]
		seq2 = s[6]
		if eval(s[2]) > eval(s[4]):
			continue

		if len(s) < 10:
			continue

#		if s[10] == 'n':
#			continue
		try:
			if s[8] == '*':
				edsites = []
			else:
				edsites = map(eval, s[8].split(','))
		except:
			continue
		try:
			if s[9] == '*':
				edsites_2 = []
			else:
				edsites_2 = map(eval, s[9].split(','))
		except:
				continue
		total_edited += len(edsites)
		if len(seq1) != len(seq2):
			continue
		for i in range(len(seq1)):
			if seq1[i] == 'A':
				total_sites += 1
			if (seq1[i] == 'A' and seq2[i] == 'G'):
				transitions_total += 1
			elif (seq1[i] == 'A' and seq2[i] in ('T','C')):
				transversions_total += 1
			if i in edsites:
				if i in edsites_2:
					continue
				try:
					edlet_dict[seq1[i]] += 1
				except:
					edlet_dict[seq1[i]] = 1
				if (seq1[i] == 'A' and seq2[i] == 'G'):
					transitions_edited+= 1
				elif (seq1[i] == 'T' and seq2[i] == 'C'):
					transitions_edited+= 1
				elif (seq1[i] == 'A' and seq2[i] in ('T','C')):
					transversions_edited += 1
				elif (seq1[i] == 'T' and seq2[i] in ('A','G')):
					transversions_edited += 1
		try:
			became_unedited += new_edsite_num(s[8], s[9])
		except:
			pass

for k in edlet_dict.keys():
	print k, edlet_dict[k]

print "||A[edited]->A[non-edited]||", became_unedited
print "||A[edited]->G||", transitions_edited
print "||A[edited]->{T,C}||", transversions_edited
print "||A->G||", transitions_total
print "||A->{T,C}||", transversions_total
print "||A[edited]||", total_edited

print (float(transitions_edited)/total_edited)/(float(transitions_total)/total_sites)
print ((float(transversions_edited)/total_edited)/(float(transversions_total)/total_sites))
print float(transitions_total)/total_sites
print float(transversions_total)/total_sites
