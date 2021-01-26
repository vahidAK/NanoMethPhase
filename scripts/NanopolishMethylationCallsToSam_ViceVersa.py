import argparse
import pysam
import os
import gzip
import bz2
import re
parser = argparse.ArgumentParser(description='Convert Nanopolish Methylation Call to sam and vice versa')
parser.add_argument('--input','-i', action='store',type=str,required=True, help='Input methylation call tsv or '
                    'sam/bam converted mathylation call file')
parser.add_argument('--tsvTosam','-ts', action='store_true',required=False
                    , help='Convert methylation call tsv to sam')
parser.add_argument('--samTotsv','-st', action='store_true',required=False
                    , help='Convert sam converted methylation call to tsv')
args = parser.parse_args()

input_file= os.path.abspath(args.input)

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

if args.tsvTosam and args.samTotsv:
    raise Exception("Pleas choos only one of the tsvTosam or samTotsv flags based"
                    " on your input")
if not args.tsvTosam and not args.samTotsv:
    raise Exception("Pleas choos one of the tsvTosam or samTotsv flags based"
                    " on your input")

check_dup_reads= set()
check_dup_reads_increment= 0
if args.tsvTosam:
    with openfile(input_file) as file:
        next(file) # ignor header
        for line in file:
            line=line.rstrip().split('\t')
            check_dup_reads.add(line[4])
            if line[1] == "+":
                strand= "0"
            else:
                strand= "16"
            if line[4] not in check_dup_reads:
                check_dup_reads_increment= 0
                outwrite= [line[4]+"_"+str(check_dup_reads_increment),strand,line[0],str(int(line[2])+1),"60",
                           str(len(line[10]))+"M","*","0","0",line[10],
                           "?"*len(line[10]),"ll:f:{}".format(line[5]),
                           "lm:f:{}".format(line[6]),"lu:f:{}".format(line[7]),
                           "nc:i:{}".format(line[8]),"nm:i:{}".format(line[9])]
                
                print('\t'.join(outwrite))
            else:
                check_dup_reads_increment += 1
                outwrite= [line[4]+"_"+str(check_dup_reads_increment),strand,line[0],str(int(line[2])+1),"60",
                           str(len(line[10]))+"M","*","0","0",line[10],
                           "?"*len(line[10]),"ll:f:{}".format(line[5]),
                           "lm:f:{}".format(line[6]),"lu:f:{}".format(line[7]),
                           "nc:i:{}".format(line[8]),"nm:i:{}".format(line[9])]
                
                print('\t'.join(outwrite))

if args.samTotsv:
    print("chromosome\tstrand\tstart\tend\tread_name\tlog_lik_ratio\tlog_lik_methylated"
      "\tlog_lik_unmethylated\tnum_calling_strands\tnum_motifs\tsequence")
    bam=pysam.AlignmentFile(input_file, 'rb')
    bamiter= bam.fetch(until_eof=True)
    for read in bamiter:
        if read.is_reverse:
            strand="-"
            
        else:
            strand="+"
  
        read_id=read.query_name.split("_")[:-1][0]
        ref_name= str(read.reference_name)
        read_seq= read.query_sequence.upper()
        cg_indexes= [(j.start()) for j in re.finditer("CG", read_seq)]
        start=read.reference_start
        end= start + (cg_indexes[-1] - cg_indexes[0])
        llr= round(read.get_tag("ll"),2)
        llm= round(read.get_tag("lm"),2)
        llu= round(read.get_tag("lu"),2)
        ncs= read.get_tag("nc")
        nm= read.get_tag("nm")
        out_write= [ref_name, strand, start, end, read_id, llr, llm, 
                    llu, ncs, nm, read_seq]
        print('\t'.join(map(str,out_write)))
