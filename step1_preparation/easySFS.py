# -*- coding: utf-8 -*-
#!/usr/bin/env python

'''
easysfs.py from https://raw.githubusercontent.com/isaacovercast/easySFS/master/easySFS.py
modified by annabel beichman -- november 2018
so that sites that are 0/1 across a population get excluded
and so that there isn't an interactive portion of the script that stops it running remotely
this script only retains bi-allelic SNPs.
~~~easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py~~~
modified by meixi lin -- june 2020
1. check input function not allowing vcf file has samples not recorded in pop
modified by meixi lin -- Sat Jan  2 00:10:22 2021
1. remove functions and modules not used in the analyses
2. change the fsc 2d projection to match the dadi specifications
#### note: you must have your projection values be in the same order as the populations are in your popMap file ########
'''

from __future__ import print_function
import matplotlib
matplotlib.use('PDF')
# from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from collections import OrderedDict
from itertools import combinations
import pandas as pd
import numpy as np
import argparse
import shutil
import gzip
import dadi
import sys
import os
import datetime

def dadi_preview_projections(dd, pops, ploidy, unfold):
	msg = """
	Running preview mode. We will print out the results for # of segregating sites
	for multiple values of projecting down for each population. The dadi
	manual recommends maximizing the # of seg sites for projections, but also
	a balance must be struck between # of seg sites and sample size.

	For each population you should choose the value of the projection that looks
	best and then rerun easySFS with the `--proj` flag.
	"""
	print(msg)
	for pop in pops:
		print(pop)
		seg_sites = {}
		## Calculate projections for up to 2 x the number of samples,
		## so this is assuming populations are diploid.
		## The +1 makes it possible to see preview including all samples per pop
		## ML 20220109: remember `max(range(0,11)) == 10`
		for x in range(ploidy, ploidy*len(pops[pop])+1):
			fs =  dadi.Spectrum.from_data_dict(dd, [pop], [x], polarized=unfold)
			s = fs.S()
			seg_sites[x] = round(s) # the value was rounded ... not necessarily the exact same value
			print("({}, {})".format(x, round(s)), end="\t")
		print("")
		print("\n")
	return None

def dadi_oneD_sfs_per_pop(dd, pops, proj, unfold, outdir, dtype):
	dadi_dir = os.path.join(outdir, "dadi")
	fsc_dir = os.path.join(outdir, "fastsimcoal2")
	plot_dir = os.path.join(outdir, "ProjectionPlots")
	M_or_D = "D" if unfold else "M"
	for i, pop in enumerate(pops):
		print("Doing 1D sfs - {}".format(pop))
		dadi_sfs_file = os.path.join(dadi_dir, pop+"-"+str(proj[i])+".sfs")

		fs = dadi.Spectrum.from_data_dict(dd, [pop], [proj[i]], mask_corners=True, polarized=unfold)

		## Do int bins rather than float
		if dtype == "int":
			dat = np.rint(np.array(fs.data))
			fs = dadi.Spectrum(dat, data_folded=fs.folded, mask=fs.mask, fill_value=0, dtype=int)

		fs.to_file(dadi_sfs_file)

		## Convert each 1D sfs to fsc format
		fsc_oneD_filename = os.path.join(fsc_dir, pop+"-"+str(proj[i])+"_{}AFpop0.obs".format(M_or_D))
		with open(fsc_oneD_filename, 'w') as outfile:
			outfile.write("1 observation\n")
			outfile.write("\t".join(["d0_"+str(x) for x in xrange(proj[i]+1)]) + "\n")
			## Grab the fs data from the dadi sfs
			with open(dadi_sfs_file) as infile:
				outfile.write(infile.readlines()[1])
				outfile.write("\n")

		## Plot the 1D SFS
		todaysdate=datetime.datetime.today().strftime('%Y%m%d')
		oneD_plot=plt.figure()
		oneD_plotname=os.path.join(plot_dir, pop+"-"+str(proj[i])+".sfs."+todaysdate+".pdf")
		dadi.Plotting.plot_1d_fs(fs, show=False)
		oneD_plot.savefig(oneD_plotname)
		plt.close(oneD_plot)
	return None


def dadi_twoD_sfs_combinations(dd, pops, proj, unfold, outdir, prefix, dtype, verbose):
	dadi_dir = os.path.join(outdir, "dadi")
	fsc_dir = os.path.join(outdir, "fastsimcoal2")
	plot_dir = os.path.join(outdir, "ProjectionPlots")
	M_or_D = "D" if unfold else "M"
	## All combinations of pairs of populations
	popPairs = list(combinations(pops, 2))
	## All combinations of corresponding projection values
	## This is hackish, it really relies on the fact that `combinations`
	## is generating combinations in the same order for pops and projs
	projPairs = list(combinations(proj, 2))
	## Make a dict for pop indices. this is a mapping of population labels
	## to values (ie. {'pop1':1, 'pop2',2}) for labeling the fsc file names
	popidx = {}
	for i, pop in enumerate(pops):
		popidx[pop] = i
	if verbose: print("Population pairs - {}".format(popPairs))
	if verbose: print("Projections for each pop pair - {}".format(projPairs))
	for i, pair in enumerate(popPairs):
		outpair = "-".join([(str(pair[0]) + str(projPairs[i][0])), (str(pair[1]) + str(projPairs[i][1]))])
		print("Doing 2D sfs - {}".format(outpair))
		dadi_joint_filename = os.path.join(dadi_dir, outpair +".sfs")
		fs = dadi.Spectrum.from_data_dict(dd, list(pair), list(projPairs[i]), polarized=unfold)

		## Do int bins rather than float
		if dtype == "int":
			dat = np.rint(np.array(fs.data))
			fs = dadi.Spectrum(dat, data_folded=fs.folded, mask=fs.mask, fill_value=0, dtype=int)

		fs.to_file(dadi_joint_filename)

		## Convert each 2D sfs to fsc format
		## Specify the filenames more clearly, dimensions should match dadi now
		## NB: FSC joint format file names look like this: (NOTE the flip in order to comply to fsc formatting)
		## <prefix>.<pop_first>-<pop_second>_jointMAFpop<second>_<first>.obs
		## Where the first pop specified is listed in the columns and first in filename (pop0; ENP48) and
		## the second pop specified is listed in the rows and second in filename (pop1; GOC30).
		fsc_twoD_filename = os.path.join(fsc_dir, prefix+".{}_joint{}AFpop{}_{}.obs".format(outpair,M_or_D, popidx[pair[1]], popidx[pair[0]]))
		with open(fsc_twoD_filename, 'w') as outfile:
			outfile.write("1 observation\n")
			## Format column headers (i.e. d0_0 d0_1 d0_2 .. d0_n for deme 0 up to sample size of n)
			outfile.write("\t" + "\t".join(["d{}_".format(popidx[pair[0]]) + str(x) for x in xrange(projPairs[i][0]+1)]) + "\n")

			## Format row headers
			row_headers = ["d{}_".format(popidx[pair[1]]) + str(x) for x in xrange(projPairs[i][1]+1)]
			## Read in the joint fs from dadi and format it nice for fsc (Read in to circumvent the float precision variations)
			## dadi's column and fsc's column were flipped (transpose)
			with open(dadi_joint_filename) as infile:
				## Get the second line of the dadi-style sfs which contains the data
				fsc_data = infile.readlines()[1].split()
				## Convert to an ndarray in numpy by column, shape should be the transpose
				fsc_data = np.asarray(fsc_data, dtype = 'str')
				fsc_data = np.reshape(fsc_data, newshape = fs.data.shape).transpose()

			if not len(row_headers) == fsc_data.shape[0]:
				print("FSC Joint SFS failed for - {}".format(pair))
				print("Row headers - {}".format(row_headers))
				print("Row data - {}".format(fsc_data))
				print("Len header-data\n{}-{}".format(len(row_headers), len(fsc_data)))
				sys.exit(1)
			else:
				pass
			## Write out each row to the file
			for rowid, row_head in enumerate(row_headers):
				outrow = fsc_data[rowid, :].tolist()
				outfile.write(row_head + "\t" + " ".join(outrow) + "\n")

		## Plot the 2D SFS
		todaysdate=datetime.datetime.today().strftime('%Y%m%d')
		twoD_plot=plt.figure()
		twoD_plotname=os.path.join(plot_dir, outpair+".sfs."+todaysdate+".pdf")
		dadi.Plotting.plot_single_2d_sfs(fs, vmin = 0.1)
		twoD_plot.savefig(twoD_plotname)
		plt.close(twoD_plot)
	return None

def make_datadict(genotypes, pops, maxHetFilter,verbose=False,ploidy=2):
	dd = {}
	hetFailSiteCounter = dict.fromkeys(pops.keys(),0)
	## Get genotype counts for each population
	for row in genotypes.iterrows():
		## iterrows() returns a tuple for some reason
		row = row[1]

		calls = {}
		for pop in pops.keys():
			## If there is a bunch of info associated w/ each snp then
			## just carve it off for now.
			pop_genotypes = [row[x].split(":")[0] for x in pops[pop]]
			ref_count = sum([x == "0" or x == "0/0" or x == "0|0" for x in pop_genotypes]) * ploidy
			alt_count = sum([x == "1" or x == "1/1" or x == "1|1" for x in pop_genotypes]) * ploidy
			## Haploids shouldn't have hets in the vcf
			het_count = sum([x == "1/0" or x == "0/1" or x == "1|0" or x == "0|1" for x in pop_genotypes])

			ref_count += het_count
			alt_count += het_count
			# 20181121: AB adding a condition to account for sites that are 0/1 across all individuals in a population. don't want to include those sites as called, since they indicate some sort of problem. 20181122, AB updated her modification to account for lines with ./., which would mean that 0/1 might be all the calls, but might not equal the length of the genotypes. Instead adding a hom ref and hom alt count and saying if they are both zero (no 0/0 or 0/1 sites, but homAlt isn't zero it means that all available calls for that pop are 0/1)
			# AB adding a hom count
			# homRef_count=sum([x == "0" or x == "0/0" or x == "0|0" for x in pop_genotypes])
			# homAlt_count=sum([x == "1" or x == "1/1" or x == "1|1" for x in pop_genotypes])
			called_gts = sum([x!="./." for x in pop_genotypes])
			#noCall_count=sum([x == "./." for x in pop_genotypes])
			#if homRef_count==0 and homAlt_count==0 and het_count!=0:
				#print("found an all 0/1 site for "+str(pop)+str(pop_genotypes))
				#calls[pop] =(0,0) # set it as though it's no-call for that population
			if het_count !=0 and het_count >= called_gts*float(maxHetFilter):
				# if verbose:
				print(row["#CHROM"]+"-"+row["POS"] + "/" + pop +": found a site with >="+str(float(maxHetFilter)*100)+"% of all calls hets. het count = "+str(het_count)+" genotypes: "+str(pop_genotypes) +"\n  dadi call would be: " +str(ref_count)+","+str(alt_count)+" converting to: 0,0")
				hetFailSiteCounter[pop] += 1
				calls[pop] = (0,0) # set it as though it's no-call for that population

			else:
				calls[pop] = (ref_count, alt_count)

		dd[row["#CHROM"]+"-"+row["POS"]] =\
			{"segregating":[row["REF"], row["ALT"]],\
			"calls":calls,\
			"outgroup_allele":row["REF"]}
	print("Total sites failed maxHetFilter - {}".format(str(hetFailSiteCounter)))
	return dd


def read_input(vcf_name, all_snps=False, verbose=False):

	## use gzip?
	if vcf_name.endswith(".gz"):
		ofunc = gzip.open
	else:
		ofunc = open
	# read all the lines in a gzipped file
	infile = ofunc(vcf_name, 'rt')
	lines = infile.readlines()
	infile.close()

	for line in lines:
		if line.startswith("#CHROM"):
			header = line

	## Just get the data lines, not the comments
	lines = [x for x in lines if not x.startswith('#')]
	if verbose:
		print("  Number of snps in input file: {}".format(len(lines)))

	## Randomly select one snp per locus
	if not all_snps:
		print("  Sampling one snp per locus")
		loci_nums = set([x.split()[0] for x in lines])
		loc_dict = {}
		for loc in loci_nums:
			loc_dict[loc] = []

		## populate the loci dict
		for line in lines:
			loc_dict[line.split()[0]].append(line)

		lines = []
		for loc in loc_dict.values():
			line = np.random.choice(loc, 1)[0]
			lines.append(line)

		## Sanity check.
		## Some snp calling pipelines use the vcf Chrom/pos information to
		## convey read/snp info per locus (ipyrad), some just assign
		## all snps to one chrom and use pos/ID (tassel).
		## If the user chooses to randomly sample one snp per block and the
		## VCF doesn't use Chrom to indicate RAD loci then it'll just
		## sample one snp for the whole dataset.
		if len(loc_dict) == 1:
			msg = """
	VCF file uses non-standard Chrom/pos information.
	We assume that Chrom indicates RAD loci and pos indicates snps within each locus
	The VCF file passed does not have rad locus info in the Chrom field.

	You can re-run the easySFS conversion with the `-a` flag to use all snps in the conversion."""
			sys.exit(msg)

		if verbose:
			print("  Using n independent snps: {}".format(len(lines)))


	## lines now here has a list of either all snps in the input
	## or a subset that includes only one snp per locus
	genotypes = pd.DataFrame([x.split() for x in lines], columns=header.split())
	return genotypes


def get_inds_from_input(vcf_name, verbose):
	# Read in the vcf file and grab the line with the individual names
	# Add the 'U' to handle opening files in universal mode, squashes the
	# windows/mac/linux newline issue.
	## use gzip?
	indnames = []
	if vcf_name.endswith(".gz"):
		ofunc = gzip.open
	else:
		ofunc = open
	try:
		with ofunc(vcf_name, 'rt') as infile:
			for line in infile:
				if line.startswith('#'):
					if line.startswith('#CHROM'):
						row = line.strip().split()
						# VCF file format spec says that if the file contains genotype
						# data then "FORMAT" will be the last column header before
						# sample IDs start
						startcol = row.index('FORMAT')
						indnames = [x for x in row[startcol+1:]]
					else:
						pass
				else:
					break
	except Exception as inst:
		msg = """
	Problem reading individuals from input VCF file."""
		print(msg)
		print("Error - {}".format(inst))
		raise

	if not indnames:
		raise Exception("No sample names found in the input vcf. Check vcf file formatting.")
	return indnames


def check_inputs(ind2pop, indnames, pops):
	## Make sure all samples are present in both pops file and VCF, give the user the option
	## to bail out if something is goofy
	pop_set = set(ind2pop.keys())
	vcf_set = set(indnames)
	## Return error if vcf_set is not a subset of pop_set (ie. samples without pop designation)
	if not vcf_set.issubset(pop_set):
		print("Samples in VCF not present in pops file: {}\n".format(", ".join(vcf_set.difference(pop_set))))
		sys.exit('Samples in VCF without pop designation')
	## remove pop_set values
	if not pop_set == vcf_set:
		print("\nSamples in pops file not present in VCF: {}\n".format(", ".join(pop_set.difference(vcf_set))))
		## Remove the offending samples from ind2pop
		# in this case "pop" is popping off the samples that are in the diff bet popset and vcfset
		map(ind2pop.pop, pop_set.difference(vcf_set))
		## Remove the offending samples from the pops dict
		# here k is population name (e.g. ENP; GOC) v is the list of individuals belonging to that population
		for k,v in pops.items():
			for ind in pop_set.difference(vcf_set):
				# try a safer change:
				if ind in v:
					v.remove(ind) # doesn't return anything, just removes the entry from v; ab modified to have this happen in two steps (correct)
					pops[k] = v
		for k,v in pops.items():
			if not v:
				print("Empty population, removing - {}\n".format(k))
				pops.pop(k)
		for k,v in pops.items():
			print(str(len(v)) + ' Surviving individuals for {0}: {1}\n'.format(k,v))

	return ind2pop, indnames, pops


def get_populations(pops_file, verbose=False):
	# Here we need to read in the individual population
	# assignments file and do this:
	# - populate the locs dictionary with each incoming population name
	# - populate another dictionary with individuals assigned to populations
	# Add the 'U' to handle opening files in universal mode, squashes the
	# windows/mac/linux newline issue.

	try:
		with open(pops_file, 'rU') as popsfile:
			ind2pop = {}
			pops = OrderedDict()

			lines = popsfile.readlines()
			## Get all the populations
			for line in lines:
				pops.setdefault(line.split()[1], [])

			for line in lines:
				ind = line.split()[0]
				pop = line.split()[1]
				ind2pop[ind] = pop
				pops[pop].append(ind)

		print("Processing {} populations - {}".format(len(pops), pops.keys()))
		if(verbose):
			for pop,ls in pops.items():
				print(pop, ls)

	except Exception as inst:
		msg = """
	Problem reading populations file. The file should be plain text with one
	individual name and one population name per line, separated by any amount of
	white space. There should be no header line in this file.
	An example looks like this:

		ind1    pop1
		ind2    pop1
		ind3    pop2
		ind4    pop2"""
		print(msg)
		print("    File you specified is: {} ".format(pops_file))
		print("    Error - {}".format(inst))
		raise

	return ind2pop, pops


def parse_command_line():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog="""\n""")

	parser.add_argument("-a", dest="all_snps", action='store_true',
		help="Keep all snps within each RAD locus (ie. do _not_ randomly sample 1 snp per locus).")

	parser.add_argument("-i", dest="vcf_name", required=True,
		help="name of the VCF input file being converted")

	parser.add_argument("-p", dest="populations", required=True,
		help="Input file containing population assignments per individual")

	parser.add_argument("--proj", dest="projections",
		help="List of values for projecting populations down to different sample sizes")

	parser.add_argument("--preview", dest="preview", action='store_true',
		help="Preview the number of segragating sites per population for different projection values.")

	parser.add_argument("-o", dest="outdir", default='output',
		help="Directory to write output SFS to")

	parser.add_argument("--ploidy", dest="ploidy", type=int, default=2,
		help="Specify ploidy. Default is 2. Only other option is 1 for haploid.")

	parser.add_argument("--prefix", dest="prefix",
		help="Prefix for all output SFS files names.")

	parser.add_argument("--unfolded", dest="unfolded", action='store_true',
		help="Generate unfolded SFS. This assumes that your vcf file is accurately polarized.")

	parser.add_argument("--dtype", dest="dtype", default="float",
		help="Data type for use in output sfs. Options are `int` and `float`. Default is `float`.")

	parser.add_argument("--GQ", dest="GQual",
		help="minimum genotype quality tolerated", default=20)

	parser.add_argument("-f", dest="force", action='store_true',
		help="Force overwriting directories and existing files.")

	parser.add_argument("-v", dest="verbose", action='store_true',
		help="Set verbosity. Dump tons of info to the screen")

	parser.add_argument("-maxHetFilter", dest="maxHetFilter", default=1.0,
		help="Fraction of called genotypes per population that are heterozygous (0/1). e.g. -maxHetFilter 0.8 would exclude any site that has >=80 percent of called genotypes 0/1 within a population. Default is 1.0 which only removes sites that are all 0/1 within the population. If your SFS is U-shaped after projection, I recommend lowering the max threshold to .7-.9.")

	## if no args then return help message
	if len(sys.argv) == 1:
		parser.print_help()
		sys.exit(1)

	## parse args
	args = parser.parse_args()
	return args

def init(args):
	## Set up output directory and output prefix (overwrite if force)
	outdir = args.outdir
	if os.path.exists(outdir) and args.force == False:
		print("\nOutput directory exists. Use -f to override.\n")
		sys.exit(1)
	if os.path.exists(outdir):
		shutil.rmtree(outdir)
	os.mkdir(outdir)
	os.mkdir(os.path.join(outdir, "dadi"))
	os.mkdir(os.path.join(outdir, "fastsimcoal2"))
	os.mkdir(os.path.join(outdir, "ProjectionPlots"))

	if not args.prefix:
		prefix = args.vcf_name.split('/')[-1].split('.')[0]
	else:
		prefix = args.prefix
	if args.verbose:
		print("Prefix - {}".format(prefix))

	return outdir, prefix


def main():
	args = parse_command_line()

	if args.verbose:
		print("Input Arguments:\n\t{}".format(args))

	## Set up output directory and output prefix
	if args.preview:
		if args.verbose: print("Doing preview so skipping directory initialization")
	else:
		outdir, prefix = init(args)

	## Get populations and populations assignments for individuals
	## ind2pop - a dictionary mapping individuals to populations
	## pops - a dictionary of populations and all inds w/in them
	ind2pop, pops = get_populations(args.populations, args.verbose)

	## Read in the names of individuals present in the vcf file
	indnames = get_inds_from_input(args.vcf_name, args.verbose)

	## Check whether inds exist in the population mapping and input vcf
	## files. Give the user an opportunity to bail if there is a mismatch.
	ind2pop, indnames, pops = check_inputs(ind2pop, indnames, pops)

	## Reads the vcf and returns a pandas dataframe
	genotypes = read_input(args.vcf_name, all_snps=args.all_snps, verbose=args.verbose)
	## Convert dataframe to dadi-style datadict
	dd = make_datadict(genotypes, pops=pops, ploidy=args.ploidy, verbose=args.verbose, maxHetFilter=args.maxHetFilter)
	## Don't write the datadict to the file for preview mode
	if not args.preview:
		with open(os.path.join(args.outdir, "datadict.txt"), 'w') as outfile:
			for x,y in dd.items():
				outfile.write(x+str(y)+"\n")

	## Do preview of various projections to determine good values
	if args.preview:
		dadi_preview_projections(dd, pops, ploidy=args.ploidy, unfold=args.unfolded)
		sys.exit()

	elif args.projections:
		## Validate values passed in for projecting
		proj = [int(x) for x in args.projections.split(",")]
		if not len(pops) == len(proj):
			msg = "You must pass in the same number of values for projection as you have populations specified"
			msg += "\n\nN pops = {}\nN projections = {}\nProjections = {}".format(len(pops), len(proj), proj)
			sys.exit(msg)

		## Create 1D sfs for each population
		dadi_oneD_sfs_per_pop(dd, pops, proj=proj, unfold=args.unfolded, outdir=outdir, dtype=args.dtype)

		## Create pairwise 2D sfs for each population
		dadi_twoD_sfs_combinations(dd, pops, proj=proj, unfold=args.unfolded, outdir=outdir, prefix=prefix, dtype=args.dtype, verbose=args.verbose)
	else:
		print("Either --preview or --proj must be specified.")

if __name__ == "__main__":
	sys.exit(main())
