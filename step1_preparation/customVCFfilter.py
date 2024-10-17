# -*- coding: utf-8 -*-
'''
Custom script for additional site- and genotype-level filtering of a VCF file
Author: Meixi Lin
2022-12-29 13:42:52

Input_1 = VCF file (unphased and gzipped)
Input_2 = Output VCF file path
Input_3 = CSV Table with maximum depth for each individuals

Output_1 = filtered VCF file. Written to the ovcf_file provided.
- If GATK standard filter has been applied, no additional analyses are performed.
- Additional Filtered sites are marked as FAIL_[?] in the 7th (FILTER) column
- Sites that pass go on to genotype filtering
- Genotypes that fail filters are added a FT field in the format together with the VCf file
- No recalculation on the INFO field was performed

NOTE on the script's activity:
1. get list of samples
2. add new headers for filter names in the VCF file
3. loop through each line of record and look for:
	1. sitefilters already applied in GATKfilter --> append sitefilters and continue;
	2. no valid INFO --> FAIL_noINFO --> quit program
	3. fail if not monomorphic or simple SNP or no VariantType annotation --> FAIL_mutType --> quit program
	---new filters---
	4. reference N --> FAIL_refN
	5. reference not A/T/C/G --> FAIL_badRef
	6. alternate not A/T/C/G/. --> FAIL_badAlt
	7. INFO is . --> FAIL_noINFO
	8. In the GT field:
		2. No genotype_DP information --> FAIL_noDPi
		3. No genotype_GQ/RGQ information --> FAIL_noGQi
	If any sitefilter exists --> continue;
	--new genotype filters---
	9. Perform individual sample GT filter (using function GTfilter):
		Passing genotypes have to suffice ALL of these conditions
		1. Biallelic allele: ('0/0','0/1','1/1', '0|0', '0|1', '1|0', '1|1') --> GTf
		2. GQ/RGQ is not '.' and GQ >=6 --> Q6 (By definition of -10logP --> P = 0.748).
		3. DP is not '.' and 4 <= DP <= maxD --> d4 Dm
		Otherwise, Add genotype filters --> next sample record
	10. Finish filtering, further site filter based on the modified GT:
		1. If sites become missing --> FAIL_noGT
		Notes: did not perform WARN_missing and WARN_excessHet since there might be individuals filtered out.


The GT filters will not be analyzed if any of the FAIL_[?] sitefilters are applied earlier). If there is no continue, the script will keep analyzing (e.g. site filters will all add up before GTfilter analyses)
'''

###########################################################
# import packages and input arguments
import os
import sys
import gzip
import csv
import argparse
import datetime

###########################################################
## define functions

############################################################
# input functions
# Read csv based filter files
def readDPfile(DP_file):
	with open(DP_file) as f:
		ff=filter(None,csv.reader(f))
		dd={str(rows[0]):int(rows[1]) for rows in ff}
	sys.stdout.write('INFO: genotype maxD from {0}.\n'.format(DP_file))
	sys.stdout.write('\n'.join(['\t{0} = {1}'.format(*tup) for tup in dd.items()])+'\n')
	return dd

def ExistingFile(fname):
	"""Return *fname* if existing file, otherwise raise ValueError."""
	if os.path.isfile(fname):
		return fname
	else:
		raise ValueError("%s must specify a valid file name" % fname)

def parse_customVCFfilter_input():
	'''
	Parse command-line arguments
	'''
	parser = argparse.ArgumentParser(description="This python script performs custom vcf filtering and outputs an uncompressed vcf file.")

	parser.add_argument(
		"vcf_file", type=ExistingFile,
		help="REQUIRED. Path to the input vcf file, need to be gzipped.")
	parser.add_argument(
		"ovcf_file", type=str,
		help="REQUIRED. Path to the output vcf file, will not be gzipped.")
	parser.add_argument(
		"maxD_file", type=ExistingFile,
		help="REQUIRED. Path to the input maximum depth file, should be csv, e.g.'sample,40'.")

	args = parser.parse_args()
	vcf_file = args.vcf_file
	ovcf_file = args.ovcf_file
	maxD_file = args.maxD_file

	# remove output if exists
	if os.path.isfile(ovcf_file):
		os.remove(ovcf_file)
		sys.stderr.write('WARNING: Removing existing output {0}.\n'.format(ovcf_file))

	return vcf_file, ovcf_file, maxD_file

############################################################
# VCFHEADER function
def format_filterheader(filtertypes, filternames):
	filteriter = zip(filtertypes, filternames)
	filterheader_list = []
	for filtertype, filtername in filteriter:
		filterheader = '##FILTER=<ID={0},Description="{1}">\n'.format(filtertype, filtername)
		filterheader_list.append(filterheader)
	sys.stdout.write('INFO: new filter headers.\n')
	for filterheader in filterheader_list:
		sys.stdout.write('\t'+filterheader)
	return filterheader_list

def format_commandheader(vcf_file, ovcf_file, maxD_file):
	currentfile = sys.argv[0]
	currenttime = datetime.datetime.today().ctime()
	commandheader = "##CustomFilterCmd=<Command=\"python2.7.15 {0} {1} {2} {3}\",Time=\"{4}\">\n".format(currentfile, vcf_file, ovcf_file, maxD_file, currenttime)
	return commandheader

def make_newvcfheader(vcf_file, filterheader_list, commandheader):
	# initiate list to store the header
	newvcfheader = []
	prevline0 = ''
	insertfilter = False
	insertformat = False
	insertcommand = False
	# the genotype filter FORMAT
	formatheader = '##FORMAT=<ID=FT,Number=.,Type=String,Description="Genotype-level filter">\n'
	# load data in text mode
	with gzip.open(vcf_file, 'rt') as inVCF:
	# Add new header lines for filters being added - for GATK compatibility
		for line0 in inVCF:
			if line0.startswith('#'):
				if line0.startswith('##FORMAT') and not insertfilter:
					# insert the new filter headers before FORMAT field
					newvcfheader.extend(filterheader_list)
					insertfilter = True
				if line0.startswith('##FORMAT=<ID=GQ') and prevline0.startswith('##FORMAT=<ID=DP') and not insertformat:
					newvcfheader.append(formatheader)
					insertformat = True
				if line0.startswith("##INFO") and not insertcommand:
					# insert the commandheader
					newvcfheader.append(commandheader)
					insertcommand = True
				# append the header
				newvcfheader.append(line0)
				prevline0 = line0
			else:
				# only read the header portion
				break
	# check that both items have been inserted
	if insertfilter and insertformat and insertcommand:
		return newvcfheader
	else:
		raise ValueError('Unable to find ##FORMAT or ##contig in vcf header.')

############################################################
# VCF LINES check
def check_siteline(line0,nfields):
	# split lines by tab delimiter
	line=line0.strip().split('\t')
	# check that the fields are right
	if len(line) != nfields:
		raise ValueError('Wrong num of fields in VCF line {0}.\n'.format(line0))
	else:
		return None

############################################################
# SAMPLE functions
# get sample names from vcf file
def get_samples(vcf_file):
	samples=[]
	with gzip.open(vcf_file, 'rt') as inVCF:
		for line0 in inVCF:
			if line0.startswith('#CHROM'):
				for ii in line0.strip().split('\t')[9:]:
					samples.append(ii)
				break
	return samples

# check the DP file has the same samples as the samples in the VCF
def check_samples(samples, maxDPDict):
	dictsamples=maxDPDict.keys()
	dictsamples.sort()
	if samples==dictsamples:
		return None
	else:
		raise ValueError('inVCF sample names not matching maxDPDict')

############################################################
# INFO functions
# split INFO field with more tolerance to INFO fields with Type=Flag
# @line7: should be called from line[7], a str in the INFO field. CANNOT BE A str with '.'
def split_info(line7):
	if line7!='.':
		infodict={}
		for x in line7.split(';'):
			xx=x.split('=')
			if len(xx)==1:
				yy=[str(xx[0]), True]
			elif len(xx)==2:
				yy=[str(xx[0]), str(xx[1])]
			else:
				raise ValueError('Wrong info field in {0}.\n'.format(x))
			infodict[yy[0]]=yy[1]
		return infodict
	else:
		raise ValueError('Wrong info format in {0}.\n'.format(line7))

# combine the infodict modified during the genotype filtering
# @infodict: should be derived from split_info(line[7]) function
def combine_info(infodict):
	infolist=[]
	for key,val in sorted(infodict.items()):
		if type(val) is bool:
			infolist.append(str(key))
		else:
			infolist.append('{0}={1}'.format(key, val))
	return infolist

############################################################
# SITE FILTERS
def skip_bygatkfilter(line0):
	line=line0.strip().split('\t')
	if line[6] not in ('.', 'PASS'):
		skipthis = True
	else:
		### Check if sites that passed previous GATK filters contain any malformed sites:
		quitfilter = []
		### Must have VariantType in INFO
		### Only accept sites that are monomorphic or simple SNPs
		if line[7]!='.' and line[3]!='N':
			infodict=split_info(line[7])
			if 'VariantType' not in infodict or \
				infodict['VariantType'] not in ('NO_VARIATION', 'SNP'):
				quitfilter.append('FAIL_mutType')
		if quitfilter != []:
			raise ValueError('Bad GATK filters {0} in VCF line {1}.\n'.format(';'.join(quitfilter), line0))
		else:
			skipthis = False
	return skipthis

def custom_sitefilter(line0):
	line=line0.strip().split('\t')
	### Site filtering:
	# Should not have any original that have already been applied in GATK
	if line[6] not in ('.', 'PASS'):
		raise ValueError('Already had GATK filters {0} in VCF line {1}.\n'.format(line[6], line0))
	# Add new filters
	sitefilter=[]
	### Reference must not be N
	if line[3]=='N':
		sitefilter.append('FAIL_refN')
	### Check REF allele
	if line[3] not in ('A','C','G','T'):
		sitefilter.append('FAIL_badRef')
	### ALT allele must be one of ATCG.
	if line[4] not in ('A','C','G','T','.'):
		sitefilter.append('FAIL_badAlt')
	### Must have INFO
	if line[7]=='.':
		sitefilter.append('FAIL_noINFO')
	### Check format field
	formatfields=line[8].split(':')
	### Check DP field (genotype format)
	if 'DP' not in formatfields:
		sitefilter.append('FAIL_noDPi')
	### Check GQ or RGQ field
	if not ('GQ' in formatfields or 'RGQ' in formatfields):
		sitefilter.append('FAIL_noGQi')
	return sitefilter

def writeline_sitefilter(line0,sitefilter):
	line=line0.strip().split('\t')
	line[6] = ';'.join(sitefilter)
	outline0 = '{0}\n'.format('\t'.join(line))
	return outline0

############################################################
# GENOTYPE filters

# integer values encoded in genotypes
def gtint2int(gtint):
	if gtint == '.':
		return 0
	else:
		return int(gtint)

def format_GT(gtdict):
	gtdict['GT'] = str(gtdict['GT'])
	gtdict['DP'] = gtint2int(gtdict['DP'])
	if 'GQ' in gtdict.keys():
		gtdict['GQ'] = gtint2int(gtdict['GQ'])
	if 'RGQ' in gtdict.keys():
		gtdict['RGQ'] = gtint2int(gtdict['RGQ'])
	return gtdict

def split_GT(line0, gtentry, formatfields):
	# type check
	gt = gtentry.split(':')
	if len(formatfields) > len(gt):
		# some does not have the full number of fields
		gt.extend(['.'] * (len(formatfields) - len(gt)))
	elif len(formatfields) == len(gt):
		pass
	else:
		raise ValueError('Wrong GT format {0} in VCF line {1}'.format(gtentry, line0))
	# format the fields
	gtdict = dict(zip(formatfields, gt))
	gtdict = format_GT(gtdict)
	return gtdict

### gtdict is the genotype dictionary for that individual
### minD is the minimum genotype depth in integer
### maxD is the maximum genotype depth for the certain individual
def GTfilter(gtdict, minD, maxD):
	gtfilter=[]
	if gtdict['GT']in ('./.', '.|.'):
		return gtfilter
	else:
		# check for malformed genotypes
		if gtdict['GT'] not in ('0/0','0/1','1/1', '0|0', '0|1', '1|0', '1|1'):
			gtfilter.append('GTf')
		# check for depth
		if gtdict['DP'] < int(minD):
			gtfilter.append('d4')
		# check for depth
		if gtdict['DP'] > int(maxD):
			gtfilter.append('Dm')
		# check for GQ or RGQ (should have at least one)
		if 'GQ' in gtdict.keys():
			if gtdict['GQ'] < 6:
				gtfilter.append('Q6')
		else:
			if gtdict['RGQ'] < 6:
				gtfilter.append('Q6')
		return gtfilter

# format out gt entry
def writefield_gtfilter(gtentry,gtft,outformatfields):
	outgtentry=gtentry.split(':')
	FTpos = outformatfields.index('FT')
	# Check if the gtentry does not have AD or DP fields
	if len(outgtentry) < FTpos:
		outgtentry.extend(['.']*(FTpos-len(outgtentry)))
		gtft = 'GTl' # change the genotype filter to GTl, genotype had not sufficient length
		print(outgtentry)
	outgtentry.insert(FTpos,gtft)
	outgttext=':'.join(outgtentry)
	return outgttext

def postgt_sitefilter(excesshet, missing, ncalled, nsamples):
	# double check that nsamples = missing + ncalled
	if nsamples != (missing + ncalled):
		raise ValueError('nsamples != (missing + ncalled)')
	sitefilter=[]
	if int(missing) == int(nsamples):
		sitefilter=['FAIL_noGT']
	# else:
	# 	if int(missing) > 0.5*(nsamples):
	# 		sitefilter.append('WARN_miss50')
	# 	if int(excesshet) > 0.75*(ncalled):
	# 		sitefilter.append('WARN_het75')
	return sitefilter

def custom_gtfilter(line0, samples, minD, maxDPDict):
	line=line0.strip().split('\t')
	outformatfields = line[8].split(':')
	outformatfields.insert(outformatfields.index('DP')+1,'FT')
	excesshet = 0
	missing = 0
	ncalled = 0
	for ii in range(0,len(samples)):
		gtentry=line[ii+9]
		gtdict = split_GT(line0=line0, gtentry=gtentry, formatfields = line[8].split(':'))
		# perform genotype filters
		gtfilter = GTfilter(gtdict,minD=minD,maxD=maxDPDict[samples[ii]])
		if gtfilter == []:
			gtft = '.'
		else:
			gtft = ';'.join(gtfilter)
		# output filtered genotype
		# convert gt filter for malformatted GT (`0/0\t`) to `GTl`
		outgtentry=writefield_gtfilter(gtentry,gtft,outformatfields)
		line[ii+9]=outgtentry
		# # count for excessive het
		# if gtdict['GT'] == '0/1' and gtft == '.':
		# 	excesshet +=1
		# count for missing genotypes
		if gtdict['GT'] in ('./.', '.|.') or gtft != '.':
			missing +=1
		# count for called genotypes:
		if gtdict['GT'] not in ('./.', '.|.') and gtft == '.':
			ncalled +=1
	# change format fields
	line[8] = ':'.join(outformatfields)
	# post gtfilter site filters
	sitefilter = postgt_sitefilter(excesshet, missing, ncalled, len(samples))
	# Should not have any original that have already been applied in GATK
	if line[6] not in ('.', 'PASS'):
		raise ValueError('Already had GATK filters {0} in VCF line {1}.\n'.format(line[6], line0))
	if sitefilter != []:
		line[6] = ';'.join(sitefilter)

	outline0 = '{0}\n'.format('\t'.join(line))
	return outline0

############################################################
## main
def main():
	############################################################
	# read inputs
	vcf_file, ovcf_file, maxD_file = parse_customVCFfilter_input()

	############################################################
	# define constants
	# Note that these depths includes possibly bad reads
	# Set minimum genotype depth at 4 (lowest tolerable read numbers)
	minD=4
	sys.stdout.write('INFO: genotype minD = {0}.\n'.format(minD))
	# Set maximum genotype depth for each sample (2.5x mean coverage)
	maxDPDict=readDPfile(maxD_file) # output a dictionary with integer values

	# Get list of samples from given vcf files
	samples=get_samples(vcf_file)
	# check if the maxDPDict has the same values as the samples
	check_samples(samples, maxDPDict)

	# Set the available filters that will be appended
	filter_types = ['FAIL_refN',
					'FAIL_badRef',
					'FAIL_badAlt',
					'FAIL_noINFO',
					'FAIL_noDPi',
					'FAIL_noGQi',
					'FAIL_noGT',
					'GTl',
					'GTf',
					'd4',
					'Dm',
					'Q6']
	filter_names = ['Reference N',
					'Bad reference alleles Not A or T or C or G',
					'Bad alternative alleles Not A or T or C or G or .',
					'No valid INFO field',
					'No DP in FORMAT field',
					'No GQ or RGQ in FORMAT field',
					'No genotype after custom filters',
					'gt did not have at least GT:DP or GT:AD:DP fields',
					'gtGT not in 0/0 or 0/1 or 1/1 or 0|0 or 0|1 or 1|0 or 1|1',
					'gtDP < 4',
					'gtDP > maxDsample',
					'gtGQ or gtRGQ < 6']

	filterheader_list = format_filterheader(filter_types, filter_names)
	# Set the custom filter command line
	commandheader = format_commandheader(vcf_file, ovcf_file, maxD_file)

	# Set new vcf headers
	newvcfheader = make_newvcfheader(vcf_file, filterheader_list, commandheader)

	############################################################
	# MAIN perform filtering line by line and output new vcf file
	with open(ovcf_file, 'wt') as outVCF:
		# first write headers
		for hh in newvcfheader:
			outVCF.write(hh)
		# now parse input VCF lines
		with gzip.open(vcf_file, 'rt') as inVCF:
			for line0 in inVCF:
				if line0.startswith('#'):
					continue
				# perform site level filtering
				check_siteline(line0,nfields=(9+len(samples)))
				# if the line has previous gatk filter and is properly filtered
				if skip_bygatkfilter(line0):
					outVCF.write(line0)
					continue
				# otherwise perform custom site filters
				else:
					sitefilter = custom_sitefilter(line0)
					### If any custom filters failed, write out line and continue
					if sitefilter!=[]:
						outline0 = writeline_sitefilter(line0,sitefilter)
						outVCF.write(outline0)
						continue
					else:
					### Genotype filtering:
						outline0 = custom_gtfilter(line0,samples,minD,maxDPDict)
						outVCF.write(outline0)

###########################################################
# run main
if __name__ == "__main__":
	sys.exit(main())
