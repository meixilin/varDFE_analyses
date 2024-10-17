# -*- coding: utf-8 -*-
'''
Title: Direct estimation of SFS from a VCF file
Author: Meixi Lin (meixilin@ucla.edu)
Date: 2022-01-09 17:08:29
'''

###########################################################
## import packages

import sys, os, argparse, datetime

# import matplotlib before dadi
import matplotlib
matplotlib.use('PDF')
# from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

import dadi

# suppress numpy future warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

###########################################################
## def functions
def fmt_rightnow():
    rightnow = '[' + datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ']'
    return rightnow

def print_args(parsed_args):
    '''
    Print parsed arguments to stdout
    '''
    sys.stdout.write('{0} Beginning execution of {1} in directory {2}\n'.format(
            fmt_rightnow(), sys.argv[0], os.getcwd()))
    sys.stdout.write(str(parsed_args)+'\n')

def parse_directSFS():
    '''
    Parse command-line arguments
    '''
    parser = argparse.ArgumentParser(description='Directly estimate SFS from vcf. Based on dadi v2.')

    parser.add_argument(
        "vcf", type=str,
        help="REQUIRED. Path to the VCF file.")

    parser.add_argument(
        "popfile", type=str,
        help="REQUIRED. Path to the population maps.")

    parser.add_argument(
        "sfsfile", type=str,
        help="REQUIRED. Path to the output SFS file.")

    parser.add_argument(
        "--proj", type=int, required=True,
        help="REQUIRED. Projection value (in diploids) to project the SFS.")

    parser.add_argument(
        "--pop", type=str, required=True,
        help="REQUIRED. The population to generate SFS.")

    parser.add_argument(
        "--unfold", default=False, action='store_true', required=False,
        help="Unfold the spectrum (default is to fold ie. unfold=False).")

    args = parser.parse_args()
    print_args(args)
    return args

def plot_1DSFS(fs,args):
    ## Plot the 1D SFS in the same directory as sfsfile
    oneD_plot=plt.figure()
    oneD_plotname=args.sfsfile.replace('.sfs','.pdf')
    dadi.Plotting.plot_1d_fs(fs, show=False)
    oneD_plot.savefig(oneD_plotname)
    plt.close(oneD_plot)
    return None

###########################################################
## def variables

###########################################################
## main
def main():
    args = parse_directSFS()

    # make datadict from vcf file
    dd = dadi.make_data_dict_vcf(vcf_filename=args.vcf, popinfo_filename=args.popfile)

    # get fs from datadict with given projection
    fs =  dadi.Spectrum.from_data_dict(data_dict=dd, pop_ids=[args.pop], projections=[args.proj], polarized=args.unfold)
    # output to file
    fs.to_file(args.sfsfile)
    # plot the 1d sfs
    plot_1DSFS(fs=fs,args=args)
    sys.stdout.write('{0} Success!\n'.format(fmt_rightnow()))

if __name__ == "__main__":
    sys.exit(main())
