#!/usr/bin/python
import glob,os,re,shutil,argparse

# Input arguments

parser = argparse.ArgumentParser(description='options for slurm file generation')
parser.add_argument('-module', type=str, action="store", dest="module", default="Anaconda3")
parser.add_argument('-partition', type=str, action="store", dest="partition", default="express")
parser.add_argument('-statsfile', type=str, action="store", dest="statsfile", default="linearmodel") # prefix for this file
parser.add_argument('-gfile', type=str, action="store", dest="gfile", default="genotyping")
parser.add_argument('-efile', type=str, action="store", dest="efile", default="/data/project/thymelab/behaviordata/socialbehavior/socialscripts/socialtestrun")
parser.add_argument('-prefix', type=str, action="store", dest="prefix", default="../testlog") # slow-speed data prefix
parser.add_argument('-hprefix', type=str, action="store", dest="hprefix", default="../hsmovie") # high-speed data prefix
parser.add_argument('-pfile', type=str, action="store", dest="pfile", default="/data/project/thymelab/behaviordata/socialbehavior/socialscripts/PlotParameters")
parser.add_argument('-esfile', type=str, action="store", dest="esfile", default="/data/project/thymelab/behaviordata/socialbehavior/socialscripts/sectionsfile")
parser.add_argument('-scriptpath', type=str, action="store", dest="scriptpath", default="/data/project/thymelab/behaviordata/socialbehavior/socialscripts/")#/data/project/thymelab/scripts/behavior_published2020/ProcessMotion/") # path for the analysis suite
parser.add_argument('-other', type=str, action="store", dest="other", default="") # args to add at the end. YOU CANNOT INCLUDE THE INITIAL DASH FOR A NEW ARG WITH ARGPARSER, WILL BE ADDED BY THIS SCRIPT"

args = parser.parse_args()
module = args.module
partition = args.partition
statsfile = args.statsfile
gfile = args.gfile
efile = args.efile
prefix = args.prefix
hprefix = args.hprefix
pfile = args.pfile
esfile = args.esfile
scriptpath = args.scriptpath
other = args.other


def make_fastqc_file(fd,worc,dirlist):
	ffile = open('newmergesubmission_'+ worc +"_" + fd + '.slurm', 'w')
	ffile.write("#!/bin/bash\n")
	ffile.write("#SBATCH -p express # Partition to submit to\n")
	ffile.write("#SBATCH -n 1 # Number of cores requested\n")
	ffile.write("#SBATCH -N 1 # Ensure that all cores are on one machine\n")
	ffile.write("#SBATCH -t 20 # Runtime in minutes\n")
	ffile.write("#SBATCH --mem=32000 # Memory per cpu in MB (see also --mem-per-cpu)\n")
	ffile.write("#SBATCH -o " + statsfile + "_" + worc + "_" + fd + "_%A_%a.out # Standard out goes to this file\n")
	ffile.write("#SBATCH -e hostname_%A_%a.err # Standard err goes to this filehostname\n")
	ffile.write("module load " + module + "\n")
	ffile.write("cd outputnewmerge_" + worc + "_" + fd + '\n')
	if os.path.exists("outputnewmerge_" + worc + "_" + fd):
		print("Output directories already generated.")
	else:
		os.mkdir("outputnewmerge_" + worc+ "_"+ fd)
	ffile.write("python " + scriptpath + "processmotiondata.py -t ")
	ffile.write("\"" + prefix + ".timestamp1.Fri, Apr 8, 2022")
	ffile.write("\" -e \"" + efile + "\" -c \"" + prefix + ".centroid1.Fri, Apr 8, 2022")
	ffile.write("\" -d \"" + prefix + ".motion1.Fri, Apr 8, 2022")
	ffile.write("\" -m \"" + hprefix + "Fri, Apr 8, 2022")
	ffile.write("_\" -g \"../arid1b_hetandwt_vs_arid1b_hom_scripted_inputgenotypeids") # THESE ARE JUST PLACEHOLDER FILES FOR THE CODE, NOT USED WHEN DOING MERGING
	ffile.write("\" -s \"" + esfile)
	ffile.write("\" -j \"" + pfile + "\"")
	ffile.write(" -n 20 -r \"../rois_string\" -graphmulti \"" + dirlist + "\"")
	if other != "":
		ffile.write(" -" + other)
	return 'newmergesubmission_'+ worc +"_" + fd +'.slurm'

current_directory = os.getcwd()
geno_dirs = [os.path.basename(d)
	for root, dirs, files in os.walk(current_directory)
	for d in dirs]
geno_dirs = list(set(geno_dirs))
geno_dirs = [item for item in geno_dirs if item.startswith("outputnew")]
print(geno_dirs)

plate_dirs = [d for d in os.listdir(current_directory)
            if os.path.isdir(os.path.join(current_directory, d))]
print(plate_dirs)
for geno in geno_dirs:
	genosuffix = "_".join(geno.split("_")[1:])
	print(genosuffix)
	cdirlist = []
	wdirlist = []
	for dirs in plate_dirs:
		if "clear" in dirs:
			dirname = "../" + dirs + "/" + geno
            		if os.path.isdir("./"+dirs+"/"+geno):
				cdirlist.append(dirname)
			#print(dirname)
	for dirs in plate_dirs:
		if "white" in dirs:
			dirname = "../" + dirs + "/" + geno
            		if os.path.isdir("./"+dirs+"/"+geno):
				wdirlist.append(dirname)
	cdirlist.sort()
	wdirlist.sort()
	make_fastqc_file(genosuffix,"white",",".join(wdirlist))
	make_fastqc_file(genosuffix,"clear",",".join(cdirlist))


sfile = open("jobsubmission.sh", 'w')
sfile.write("#!/bin/bash\n")
all_files = os.listdir(current_directory)
filtered_files = [file for file in all_files if file.startswith("newmergesubmission_") and file.endswith("slurm")]
for filename in filtered_files:
	sfile.write("sbatch ")
	sfile.write(filename)
	sfile.write("\nsleep 1\n")
os.system("chmod +x jobsubmission.sh")
