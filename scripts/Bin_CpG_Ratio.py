import argparse
import gzip
import bz2 
import os
import pysam
parser = argparse.ArgumentParser(description="Number of CpG and Ratio At bins")
parser.add_argument("--reference", "-r", action="store", required=True,type=str,help="")
parser.add_argument("--bin_size", "-bs", action="store",required=False,default= 3000,type=int,help="Number CpG in interval")

args = parser.parse_args()
reference= os.path.abspath(args.reference)
bin_size= args.bin_size
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

ref_chroms= set()
with openfile(reference) as ref:
    ref_chroms= set()
    for line in ref:
        line=line.rstrip()
        if line.startswith('>'):
            line=line.split()[0]
            line=line.split('>')[1]
            ref_chroms.add(line)
    ref_chroms= sorted(ref_chroms)

fasta= pysam.FastaFile(reference)

print("Chrom\tStar\tEnd\tNumber_of_CpG\tCpG_precentage")
for chrom in ref_chroms:
    sequence= fasta.fetch(reference=chrom)
    sequence= sequence.upper()
    count_dict=dict()
    for i in range(0,len(sequence),bin_size):
        count_dict[(chrom,i+1,i+bin_size)]= sequence[i:i+bin_size+1].count('CG')
    for key,val in count_dict.items():
        print('\t'.join(map(str,key))+'\t'+
              str(val)+'\t'+str(round((val/bin_size)*100,2)))
    
    
