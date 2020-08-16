import multiprocessing as mp
import argparse
import os
import pysam
import sys
from collections import defaultdict
from tqdm import tqdm
from itertools import repeat
import statistics 

parser = argparse.ArgumentParser(description='Calculate Mean of Qualities accross SNVs.')
parser.add_argument("--input", "-i", action="store", type=str,required=False,help="The path to the bed or vcf tru positive or false positive SNVs.") 
parser.add_argument("--bam", "-b", action="store", type=str,required=True,help="The path to the alignment bam file") 
parser.add_argument("--reference", "-r", action="store", type=str,required=True, default=None,help="The path to the reference file. File must be indexed by samtools faidx"
                          "give the path to reference file. Fasta file must be already indexed using samtools faidx.")
parser.add_argument("--mappingQuality", "-mq", action="store", type=int,default= 20,required=False,help="Cutt off for filtering out low quality mapped reads from bam. Default is 20")
parser.add_argument("--threads", "-t", action="store", type=int,required=False, default=4,help="Number of threads")
parser.add_argument("--chunk_size", "-rc", action="store", type=int,required=False, default=50,help="Number of reads send to each processes for parrallel processing. Default is 10.\n"
                           "If you do not have enough memory use low numbers.")

args = parser.parse_args()
if not args.input and not args.all:
    raise Exception("You did not give an input and also did not specify --all option. One of them must be given")

def window_mutation(window_list, bam_file,mq,reference):
    bam= pysam.AlignmentFile(bam_file, 'rb')
    out_dict= defaultdict(list)
    for window in window_list:
        chrom,windowstart,windowend,base_pos = window# thsese are 9-mers
        try:
            fasta= pysam.FastaFile(reference)
        except:
            raise Exception("Cannot load reference file.")
        try:
            try:
                window_seq= fasta.fetch(reference=chrom, start=windowstart, end=windowend)
            except:
                window_seq= fasta.fetch(reference=chrom[3:], start=windowstart, end=windowend)
        except:
#            try: Do no use this part because may causes lines with less dimention and error
#                window_seq= fasta.fetch(reference=chrom, start=windowstart+4, end=windowend)
#            except:
                continue
                
        window_seq=str(window_seq.upper())
        try:
            sam_pileup= bam.pileup(chrom, windowstart, windowend,
                                      truncate= True,min_base_quality = 0,
                                      min_mapping_quality = mq)
        except:
            sam_pileup= bam.pileup(chrom[3:], windowstart, windowend,
                                      truncate= True,min_base_quality = 0,
                                      min_mapping_quality = mq)

        for pileupcolumn in sam_pileup:
            pileupcolumn.set_min_base_quality(0)
            for pileupread in pileupcolumn.pileups:
                if not pileupread.alignment.is_supplementary:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    phred= pileupread.alignment.query_qualities[pileupread.query_position]
                    read_base= pileupread.alignment.query_sequence[pileupread.query_position]
                    read_base= read_base.upper()
                    if read_base != window_seq:
                        out_dict[(chrom,base_pos,"NotRef")].append(phred)
                    elif read_base == window_seq:
                        out_dict[(chrom,base_pos,"Ref")].append(phred)
                        
    return out_dict

                                
def main():
    bam_file= os.path.abspath(args.bam)
    mq= args.mappingQuality
    threads= args.threads
    chunk= args.chunk_size
    reference= os.path.abspath(args.reference)
    window_list=list()
    feed_list= list()
    bed_file= os.path.abspath(args.input)
    if not bed_file.endswith('.bed') and not bed_file.endswith('.vcf'):
        raise TypeError("Input file is not a vcf or bed file.")
    bed= open(bed_file)
    if bed_file.endswith('.bed'):
        next(bed)
    all_lines= []
    for line in bed:
        if line.startswith('#'):
            continue
        all_lines.append(1)
    all_lines = [all_lines[x:x+chunk]
                         for x in range(0, len(all_lines), chunk)]
    all_lines = [all_lines[x:x+threads]
                         for x in range(0, len(all_lines), threads)]
    bed.close()
    bed= open(bed_file)
    if bed_file.endswith('.bed'):
        next(bed)
    with tqdm(total=len(all_lines),desc="Processing: ",
                      bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                      ) as pbar:
        print("Chromosom\tZeroBasedPosition\tRef/Alt\tAverageQuality\tNumOfBases")
        for line in bed:
            if line.startswith('#'):
                continue
            line=line.rstrip().split('\t')
            chrom= line[0]
            if 'chr' not in chrom.lower():
                chrom= 'chr'+chrom
            if bed_file.endswith('.bed'):
                windowinfo= (chrom, int(line[1]),int(line[1])+1,int(line[1]))
            else:
                windowinfo= (chrom, int(line[1])-1,int(line[1]),int(line[1])-1)
            if windowinfo not in window_list:
                window_list.append(windowinfo)
                if len(window_list) > chunk:
                    feed_list.append(window_list)
                    window_list= list()
                if len(feed_list)== threads:
                    p = mp.Pool(threads)
                    results = p.starmap(window_mutation,
                            list(zip(feed_list,
                                     repeat(bam_file), repeat(mq),repeat(reference))))
                    p.close()
                    p.join()
                    for out_dict in results:
                        if out_dict is not None:
                            for key,val in out_dict.items():
                                if len(val) > 1:
                                    mean_qual= statistics.mean(val)
                                else:
                                    mean_qual= val[0]
                                print('\t'.join(map(str,key))+'\t'+str(mean_qual)+
                                  '\t'+str(len(val)))
                    feed_list = []            
                    pbar.update(1)                    
                    
        else:
            if feed_list or window_list:
                feed_list.append(window_list)
                p = mp.Pool(len(feed_list))
                results = p.starmap(window_mutation,
                            list(zip(feed_list,
                                     repeat(bam_file), repeat(mq),repeat(reference))))
                p.close()
                p.join()
                for out_dict in results:
                    if out_dict is not None:
                        for key,val in out_dict.items():
                            if len(val) > 1:
                                mean_qual= statistics.mean(val)
                            else:
                                mean_qual= val[0]
                            print('\t'.join(map(str,key))+'\t'+str(mean_qual)+
                                  '\t'+str(len(val)))
                pbar.update(1)                    
                feed_list = []
    bed.close()
    sys.stderr.write("Job Finished\n")

if __name__ == '__main__':
    main()
