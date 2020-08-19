import multiprocessing as mp
import argparse
import os
import pysam
import sys
from collections import defaultdict
from tqdm import tqdm
from itertools import repeat
import statistics 

parser = argparse.ArgumentParser(description='Calculate reads un assigned, filtered, HP1, HP2.')
parser.add_argument("--hp1_bam", "-hp1", action="store", type=str,required=False,help="The path to the HP1 bam file.") 
parser.add_argument("--hp2_bam", "-hp2", action="store", type=str,required=False,help="The path to the HP2 bam file.")
parser.add_argument("--per_read", "-pr", action="store", type=str,required=False,help="The path to per read infor file.")
parser.add_argument("--bam", "-b", action="store", type=str,required=True,help="The path to the alignment bam file") 
parser.add_argument("--min_base_quality", "-mbq", action="store", type=int,required=False,default= 7,help="minimum base quality. Default is 7")
parser.add_argument("--mappingQuality", "-mq", action="store", type=int,default= 20,required=False,help="Cutt off for filtering out low quality mapped reads from bam. Default is 20")
parser.add_argument("--out_read", "-or", action="store", type=str,required=False,help="If you wish to also output read IDs and category as well give the path to the output "
                    "fo this file. Note: main stats will be output to the stdout")
args = parser.parse_args()

def main():
    bam_file= os.path.abspath(args.bam)
    HP1= os.path.abspath(args.hp1_bam)
    HP2= os.path.abspath(args.hp2_bam)
    per_read= os.path.abspath(args.per_read)
    mq= args.mappingQuality
    mbq= args.min_base_quality
    if args.out_read is not None:
        out_reads= open(os.path.abspath(args.out_read),'w')
        out_reads.write("ReadID\tChromosome\tStrand\tReadLength\tCategory\n")
    all_reads= defaultdict(list)
    all_reads_info= dict()
    HP1_reads= defaultdict(list)
    HP2_reads= defaultdict(list)
    assigned_reads= defaultdict(list)
    bam = pysam.AlignmentFile(bam_file, 'rb')
    bamiter = bam.fetch(until_eof=True)
    for read in bamiter:
        mp_quality = read.mapping_quality
        if (read.is_unmapped or mp_quality < mq or 
            read.is_secondary or read.is_supplementary or 
            read.is_qcfail or read.is_duplicate):
            continue
        ref_name = str(read.reference_name)
        if 'chr' not in ref_name.lower():
            ref_name = 'chr'+ref_name
        read_id= read.query_name
        all_reads[ref_name].append(read_id)
        if args.out_read is not None:
            read_len = read.query_alignment_length
            if read.is_reverse:
                strand = "-"
            else:
                strand = "+"
            if read_id not in all_reads_info:
                all_reads_info[read_id]= [ref_name,strand,read_len]
            else:
                if all_reads_info[read_id][2] < read_len:
                    all_reads_info[read_id]= [ref_name,strand,read_len]
    bam = pysam.AlignmentFile(HP1, 'rb')
    bamiter = bam.fetch(until_eof=True)
    for read in bamiter:
        ref_name = str(read.reference_name)
        if 'chr' not in ref_name.lower():
            ref_name = 'chr'+ref_name
 
        HP1_reads[ref_name].append(read.query_name)
    bam = pysam.AlignmentFile(HP2, 'rb')
    bamiter = bam.fetch(until_eof=True)
    for read in bamiter:
        ref_name = str(read.reference_name)
        if 'chr' not in ref_name.lower():
            ref_name = 'chr'+ref_name
 
        HP2_reads[ref_name].append(read.query_name)
    with open(per_read) as info:
        for line in info:
            if line.startswith("#"):
                continue
            line= line.rstrip().split("\t")
            read_id= line[3]
            chrom= line[0]
            info_line= line[8].split(',')
            for pos_phred in info_line:
                pos,phred= pos_phred.split(':')[0:2]
                if int(phred) >= mbq:
                    assigned_reads[chrom].append(read_id)
                    break
    final_output= dict()
    for chrom, reads in all_reads.items():
        all_reads= set(reads)
        unassigned_reads= all_reads - set(assigned_reads[chrom])
        hp1_reads= set(HP1_reads[chrom])
        hp2_reads= set(HP2_reads[chrom])
        filtered_reads= set(assigned_reads[chrom]) - set(HP1_reads[chrom] + HP2_reads[chrom])
        final_output[chrom]= [len(all_reads),
                    len(unassigned_reads), 
                    len(hp1_reads),
                    len(hp2_reads),
                    len(filtered_reads)]
        if args.out_read is not None:
            for unassign in unassigned_reads:
                out_reads.write(unassign+'\t'+
                                '\t'.join(map(str,all_reads_info[unassign]
                                +['Unassigned']))+'\n')
            for hp1 in hp1_reads:
                out_reads.write(hp1+'\t'+
                                '\t'.join(map(str,all_reads_info[hp1]
                                +['HP1']))+'\n')
            for hp2 in hp2_reads:
                out_reads.write(hp2+'\t'+
                                '\t'.join(map(str,all_reads_info[hp2]
                                +['HP2']))+'\n')
            for filtered in filtered_reads:
                out_reads.write(filtered+'\t'+
                                '\t'.join(map(str,all_reads_info[filtered]
                                +['Filtered']))+'\n')
    if args.out_read is not None:
        out_reads.close()
    print("Chromosome\tAll_passed_reads\tUnassignedReads\tHP1Reads\tHP2Reads\tFilteredReads")
    for chrom,info in final_output.items():
        print(chrom+'\t'+'\t'.join(map(str,info)))
        
        
if __name__ == '__main__':
    main()
