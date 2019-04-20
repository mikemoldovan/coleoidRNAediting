"""
Determine read quality in a given directory
Input:
-wd
"""

from optparse import OptionParser
from os import system, listdir

def qual_analysis(working_dir):
	for file_ in listdir(working_dir):
		if file_.split('.')[-1] == "fastq" or file_.split('.')[-1] == "fq":
			dirname = '.'.join(file_.split('.')[:-1])
			try:
				system("mkdir " + working_dir + dirname)
			except:
				pass
#			print "mkdir " + working_dir + dirname
			system("/home/mmoldovan/tools/FastQC/fastqc " + working_dir + file_ + ' -o ' + working_dir + dirname)
#			print "/home/mmoldovan/tools/FastQC/fastqc " + working_dir + filename + ' -o ' + working_dir + dirname
			system("mv " + working_dir + dirname + '/*.html ' + working_dir)
#			print "mv " + working_dir + dirname + '/*.html ' + working_dir


parser = OptionParser()
parser.add_option("-w", "--wd", help="working directory")
opt, args = parser.parse_args()

qual_analysis(opt.wd)

