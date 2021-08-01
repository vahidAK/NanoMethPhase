from scipy import stats
import argparse
import os
import gzip
import bz2
from collections import defaultdict
import numpy as np
np.random.seed(123)
parser = argparse.ArgumentParser(description='resampling and binomial test for count data ASE')
parser.add_argument("--input", "-i", action="store", type=str,required=True,help="The path to the input file."
                    "Column 1 must be DMR or interval, 2 read counts for reference allele or alternative"
                    " allele (in case for trio, maternal or paternal count) for a SNV at DMR and 3 total read counts.")
parser.add_argument("--p_binom", "-p", action="store", type=float,default= 0.5
                    ,required=False,help="P val of binom. Default is 0.5")
parser.add_argument("--type_binom", "-t", action="store", type=str,default= "two-sided"
                    ,required=False,help="{two-sided, greater, less}. Default is two-sided")
parser.add_argument("--combine_pval", "-c", action="store", type=str,default= "fisher"
                    ,required=False,help="{fisher, pearson, tippett, stouffer, mudholkar_george}. Default is fisher")


args = parser.parse_args()
input_data= os.path.abspath(args.input)
p_binom= args.p_binom
type_binom= args.type_binom
combine_pval= args.combine_pval


def openfile(file):
    '''
    Opens a file
    '''
    if file.endswith('.gz'):
        opened_file = gzip.open(file,'rt')
    elif file.endswith('bz') or file.endswith('bz2'):
        opened_file = bz2.open(file,'rt')
    else:
        opened_file = open(file,'rt')
    return opened_file

def binomial(ref,total,p_binom,type_binom):
    pvals= list()
    for i,j in zip(ref,total):
        pvals.append(stats.binom_test(i,j, p=p_binom,alternative= type_binom))
    return pvals

def combine_pvals(pvals,combine_pval):
    return stats.combine_pvalues(pvals, method= combine_pval)

with openfile(input_data) as file:
    next(file)
    SNV_dict_ref=defaultdict(list)
    SNV_dict_total=defaultdict(list)
    for line in file:
        line=line.rstrip().split('\t')
        SNV_dict_ref[line[0]].append(int(line[1]))
        SNV_dict_total[line[0]].append(int(line[2]))
print("DMR\tNum_SNV_At_DMR\tTest_statistics\tAdjustedPval")
for key in SNV_dict_ref.keys():
    num_sample= len(SNV_dict_ref[key])
    pvals= binomial(SNV_dict_ref[key],SNV_dict_total[key],p_binom,type_binom)
    print(key,str(num_sample),'\t'.join(map(str,combine_pvals(pvals,combine_pval))),sep='\t')
