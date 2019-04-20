"""
NOT USED

sra_processing

Input: 
1. File in format: run_id <f/fr>
f is 'forward', 'fr' is forward-reverse
2. working directory

Output:
List of HTML files with quality scores

workflow:
1. fastq download
2. fastq analysis

/home/mmoldovan/tools/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump --split-3 -I DRR001452.fastq 
/home/mmoldovan/tools/FastQC/fastqc DRR001452.fastq -o DRR001452.fastq_qual
"""

from optparse import OptionParser
from os import system, listdir

def infile_parse(infile_name):
	sra_id_dict = dict()
	with open(infile_name) as infile:
		for s in infile:
			s = s.strip().split()
			sra_id_dict[s[0]] = s[1]
	return sra_id_dict


def fastq_download(sra_id_dict, working_dir):
	system("cd "+working_dir)
	for k in sra_id_dict.keys():
		if sra_id_dict[k] == 'f':
			system("/home/mmoldovan/tools/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  "+k)
#			print "/home/mmoldovan/tools/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump  "+k
#			print "/home/mmoldovan/tools/FastQC/fastqc "+k
		else:
			system("/home/mmoldovan/tools/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump "+"--split-3 "+k)
#			print "/home/mmoldovan/tools/sratoolkit.2.8.2-1-centos_linux64/bin/fastq-dump "+"--split-3 "+k
#			print "/home/mmoldovan/tools/FastQC/fastqc "+"--split-3 "+k


def qual_analysis(working_dir):
	for file_ in listdir(working_dir):
		if file_.split('.')[-1] == "fastq":
			dirname = '.'.join(file_.split('.')[:-1])
			system("mkdir " + working_dir + dirname)
#			print "mkdir " + working_dir + dirname
			system("/home/mmoldovan/tools/FastQC/fastqc " + working_dir + filename + ' -o ' + working_dir + dirname)
#			print "/home/mmoldovan/tools/FastQC/fastqc " + working_dir + filename + ' -o ' + working_dir + dirname
			system("mv " + working_dir + dirname + '/*.html ' + working_dir)
#			print "mv " + working_dir + dirname + '/*.html ' + working_dir


parser = OptionParser()
parser.add_option("-i", "--infile", help="File in format: run_id <f/fr> f is 'forward', 'fr' is forward-reverse")
parser.add_option("-w", "--wd", help="working directory")
opt, args = parser.parse_args()

sra_id_dict = infile_parse(opt.infile)
fastq_download(sra_id_dict, opt.wd)
qual_analysis(opt.wd)

