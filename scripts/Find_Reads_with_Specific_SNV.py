import multiprocessing as mp
import argparse
import os
import pysam
from collections import defaultdict
from tqdm import tqdm
from itertools import repeat

parser = argparse.ArgumentParser(description='Extracting reads with the alt base at a SNV and print 9bp reference window.')
parser.add_argument("--input", "-i", action="store", type=str,required=True,help="The path to the input bed like file or input vcf (File must end with .vcf). " 
                    "NOTE: Input must contain only SNVs." 
                    "If bed like is given, format must be chromosome\tZeroBased_mutation_position\tAltBase.") 
parser.add_argument("--bam", "-b", action="store", type=str,required=True,help="The path to the alignment bam file") 
parser.add_argument("--reference", "-r", action="store", type=str,required=True, default=None,help="The path to the reference file. File must be indexed by samtools faidx")
parser.add_argument("--mappingQuality", "-mq", action="store", type=int,default= 20,required=False,help="Cutt off for filtering out low quality mapped reads from bam. Default is 20")
parser.add_argument("--threads", "-t", action="store", type=int,required=False, default=4,help="Number of threads")
parser.add_argument("--chunk_size", "-rc", action="store", type=int,required=False, default=50,help="Number of reads send to each processes for parrallel processing. Default is 50.\n"
                           "If you do not have enough memory use low numbers.")

args = parser.parse_args()
if not args.input and not args.all:
    raise Exception("You did not give an input and also did not specify --all option. One of them must be given")

def window_mutation(window_list, bam_file,mq,reference):
    bam= pysam.AlignmentFile(bam_file, 'rb')
    read_dict= defaultdict(list)
    for window in window_list:
        chrom,windowstart,windowend,base_pos,mutation = window# thsese are 9-mers
        try:
            fasta= pysam.FastaFile(reference)
        except:
            raise Exception("Cannot load reference file.")
        try:
            try:
                window_seq= fasta.fetch(reference=chrom, start=windowstart, end=windowend)
                ref_base= fasta.fetch(reference=chrom, start=base_pos, end=base_pos+1)
            except:
                window_seq= fasta.fetch(reference=chrom[3:], start=windowstart, end=windowend)
                ref_base= fasta.fetch(reference=chrom, start=base_pos, end=base_pos+1)
        except:
#            try: Do no use this part because may causes lines with less dimention and error
#                window_seq= fasta.fetch(reference=chrom, start=windowstart+4, end=windowend)
#            except:
                continue
                
        window_seq=str(window_seq.upper())
        try:
            pileupcolumns= bam.pileup(chrom, windowstart, windowend,
                                      truncate= True,min_base_quality = 0,
                                      min_mapping_quality = mq)
        except:
            pileupcolumns= bam.pileup(chrom[3:], windowstart, windowend,
                                      truncate= True,min_base_quality = 0,
                                      min_mapping_quality = mq)

        for pileupcolumn in pileupcolumns:
            pileupcolumn.set_min_base_quality(0)
            if pileupcolumn.pos == base_pos:
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    if pileupread.alignment.is_reverse:
                        strand = "-"
                    else:
                        strand = "+"
                    read= pileupread.alignment.query_name
                    read_len= pileupread.alignment.query_alignment_length
                    phred= pileupread.alignment.query_qualities[pileupread.query_position]
                    read_base= pileupread.alignment.query_sequence[pileupread.query_position]
                    read_base= read_base.upper()
                    if read_base == mutation:
                        read_dict[(*window,window_seq)].append((ref_base.upper(),read,read_len,
                              phred,strand))

    return read_dict

                                
def main():
    bam_file= os.path.abspath(args.bam)
    mq= args.mappingQuality
    threads= args.threads
    chunk= args.chunk_size
    reference= os.path.abspath(args.reference)
    window_list=list()
    feed_list= list()
    bed_file= os.path.abspath(args.input)
    if not bed_file.endswith(".vcf"):
        bed= open(bed_file)
        next(bed)
        all_lines= []
        for line in bed:
            all_lines.append(1)
        all_lines = [all_lines[x:x+chunk]
                             for x in range(0, len(all_lines), chunk)]
        all_lines = [all_lines[x:x+threads]
                             for x in range(0, len(all_lines), threads)]
        bed.close()
    else:
        bed= open(bed_file)
        all_lines= []
        for line in bed:
            if line.startswith("#"):
                continue
            all_lines.append(1)
        all_lines = [all_lines[x:x+chunk]
                             for x in range(0, len(all_lines), chunk)]
        all_lines = [all_lines[x:x+threads]
                             for x in range(0, len(all_lines), threads)]
        bed.close()
    if not bed_file.endswith(".vcf"):
        bed= open(bed_file)
        next(bed)
    else:
        bed= open(bed_file)
    print("Chromosome\tWindowStart_ZeroBased\tWindowEnd_ZeroBased\t"
          "Mutation_Pos_ZeroBased\tMutationBase\tWindow\tRefBase\tRead\t"
          "Read_Length\tBasePhred\tStrand")
    snv_count= defaultdict(int)
    with tqdm(total=len(all_lines),desc="Processing: ",
                      bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                      ) as pbar:
        for line in bed:
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            chrom= line[0]
            if not bed_file.endswith(".vcf"):
                windowinfo= (chrom, int(line[1])-4,int(line[1])+5,int(line[1]),line[2])
            else:
                windowinfo= (chrom, int(line[1])-4-1,int(line[1])+5-1,int(line[1])-1,line[4])
            if windowinfo not in window_list:
                window_list.append(windowinfo)
                if len(window_list) >= chunk:
                    feed_list.append(window_list)
                    window_list= list()
                if len(feed_list)== threads:
                    p = mp.Pool(threads)
                    results = p.starmap(window_mutation,
                            list(zip(feed_list,
                                     repeat(bam_file), repeat(mq),repeat(reference))))
                    p.close()
                    p.join()
                    for read_dict in results:
                        if read_dict is not None:
                            for key,val in read_dict.items():
                                for read in val:
                                    print('\t'.join(map(str,key))+'\t'+'\t'.join(map(str,read)))
                    
                    pbar.update(1)                    
                    feed_list = []
        else:
            feed_list.append(window_list)
            p = mp.Pool(len(feed_list))
            results = p.starmap(window_mutation,
                    list(zip(feed_list,
                             repeat(bam_file), repeat(mq),repeat(reference))))
            p.close()
            p.join()
            for read_dict in results:
                        if read_dict is not None:
                            for key,val in read_dict.items():
                                for read in val:
                                    print('\t'.join(map(str,key))+'\t'+'\t'.join(map(str,read)))
            pbar.update(1)                    
            feed_list = []
    bed.close()

if __name__ == '__main__':
    main()
