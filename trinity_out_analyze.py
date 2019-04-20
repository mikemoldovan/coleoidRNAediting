"""
Get assembly info
Input:
Assembly in FASTA format
"""

from optparse import OptionParser

def length_dist_obtain(fafile):
	length_dist = []
	with open(fafile) as f:
		for s in f:
			if s[0] == '>':
				s = s.split()
				length_dist.append(eval(s[1][4:]))
	return sorted(length_dist)[::-1]


def length_dist_analysis(length_dist, step):
	l = len(length_dist)
	s = sum(length_dist)
	print "total transcripts:", l
	print "mean length:", float(s)/l 
	print "max length:", length_dist[0]
	print "min length:", length_dist[-1]
	curr = 0.1
	curr_sum = 0
	for i in length_dist:
		curr_sum += i
		if float(curr_sum)/s > curr:
			print "N{}: {}".format(curr*100, i)
			curr += 0.1
	for i in range(0, l, step):
		print "K{}: {}".format(i, length_dist[i])



parser = OptionParser()
parser.add_option("-i", "--fafile", help="File in FASTA format, TRINITY output")
parser.add_option("-s", "--step", help="Step in read quality analysis (recommended 10000)")
opt, args = parser.parse_args()


length_dist = length_dist_obtain(opt.fafile)
length_dist_analysis(length_dist, eval(opt.step))