import multiprocessing as mp
import argparse
import os
import pysam
from collections import defaultdict
from tqdm import tqdm
from itertools import repeat
import warnings
parser = argparse.ArgumentParser(description='Finding reads mapped to ref/alt at each SNV')
parser.add_argument("--input", "-i", action="store", type=str,required=True,help="The path to the input SNV vcf file. File must contain only SNVs.") 
parser.add_argument("--bam", "-b", action="store", type=str,required=True,help="The path to the alignment bam file") 
parser.add_argument("--mappingQuality", "-mq", action="store", type=int,default= 0
                    ,required=False,help="Cutt off for filtering out low quality mapped reads from bam. Default is 0")
parser.add_argument("--baseQuality", "-bq", action="store", type=int,default= 0
                    ,required=False,help="Cutt off to skip low quality bases in read. Default is 0")
parser.add_argument("--window", "-w",
                          action="store",
                          type=str,
                          required=False,
                          help=("if you want to only fase read for a region "
                                "or chromosome. You must insert region like "
                                "this chr1 or chr1:1000-100000."))
parser.add_argument("--threads", "-t", action="store", type=int,required=False, default=4,help="Number of threads")
parser.add_argument("--chunk_size", "-rc", action="store", type=int,required=False, 
                    default=50,help="Number of reads send to each processes for parrallel processing. Default is 10.\n"
                           "If you do not have enough memory use low numbers.")

args = parser.parse_args()

def window_mutation(window_list, bam_file,mq,baseQuality):
    bam= pysam.AlignmentFile(bam_file, 'rb')
    read_dict_alt= defaultdict(list)
    read_dict_ref= defaultdict(list)
    for window in window_list:
        chrom,base_pos,ref_base,mutation = window
        try:
            pileupcolumns= bam.pileup(chrom, base_pos-1, base_pos+1,
                                  truncate= True,min_base_quality = baseQuality,
                                  min_mapping_quality = mq)
        except:
            warnings.warn("{}:{} does not exist in the alignment file.".format(chrom,base_pos+1))
            continue
        for pileupcolumn in pileupcolumns:
            pileupcolumn.set_min_base_quality(baseQuality)
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
                        read_dict_alt[window].append((read,read_len,
                              phred,strand,"alt"))
                    if read_base == ref_base:
                        read_dict_ref[window].append((read,read_len,
                              phred,strand,"ref"))

    return (read_dict_alt, read_dict_ref)

                                
def main():
    bam_file= os.path.abspath(args.bam)
    mq= args.mappingQuality
    threads= args.threads
    chunk= args.chunk_size
    baseQuality= args.baseQuality
    window_list=list()
    feed_list= list()
    vcf_file= os.path.abspath(args.input)
    vcf= open(vcf_file)
    all_lines= []
    for line in vcf:
        if line.startswith("#"):
            continue
        line=line.rstrip().split('\t')
        if (line[9].startswith('0|1') or line[9].startswith('1|0') or
            line[9].startswith('1/0') or line[9].startswith('0/1')):
            all_lines.append(1)
    all_lines = [all_lines[x:x+chunk]
                         for x in range(0, len(all_lines), chunk)]
    all_lines = [all_lines[x:x+threads]
                         for x in range(0, len(all_lines), threads)]
    vcf.close()
    vcf= open(vcf_file)
    print("Chromosome\tSNV_Pos_ZeroBased\tRef_Base\tMutationBase\t"
          "RefCount\tAltCount")
    allelic_count_ref= defaultdict(list)
    allelic_count_alt= defaultdict(list)
    if args.window is not None:
        window= args.window.replace(',','')
        window_chrom= window.split(':')[0]
        if ":" not in window and "-" not in window:
            window_start= 0
            window_end= 9000000000
        elif ":" in window and "-" in window:
            window_start = int(window.split(':')[1].split('-')[0])
            window_end = int(window.split(':')[1].split('-')[1])
        else:
            raise Exception("Enter the window correctly.")
    with tqdm(total=len(all_lines),desc="Processing: ",
                      bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                      ) as pbar:
        for line in vcf:
            windowinfo= tuple()
            if line.startswith("#"):
                continue
            line=line.rstrip().split('\t')
            chrom= line[0]
            if args.window is not None:
                if (window_chrom == line[0] and 
                    int(line[1]) >= window_start and int(line[1]) <= window_end 
                        and (line[9].startswith('0|1') or line[9].startswith('1|0') or
                        line[9].startswith('1/0') or line[9].startswith('0/1'))):
                    windowinfo= (chrom,int(line[1])-1,line[3].upper(),line[4].upper())
            else:
                if (line[9].startswith('0|1') or line[9].startswith('1|0') or
                         line[9].startswith('1/0') or line[9].startswith('0/1')):
                    windowinfo= (chrom,int(line[1])-1,line[3].upper(),line[4].upper())
            if windowinfo and windowinfo not in window_list:
                window_list.append(windowinfo)
                if len(window_list) >= chunk:
                    feed_list.append(window_list)
                    window_list= list()
                if len(feed_list)== threads:
                    p = mp.Pool(threads)
                    results = p.starmap(window_mutation,
                            list(zip(feed_list,
                                     repeat(bam_file), repeat(mq),repeat(baseQuality))))
                    p.close()
                    p.join()
                    for read_dict_alt,read_dict_ref in results:
                        if read_dict_alt is not None:
                            for key,val in read_dict_alt.items():
                                for read in val:
                                    allelic_count_alt[key].append(read[0])
                        if read_dict_ref is not None:
                            for key,val in read_dict_ref.items():
                                for read in val:
                                    allelic_count_ref[key].append(read[0])
                    all_key= set(list(allelic_count_ref.keys())+list(allelic_count_alt.keys()))
                    for key in all_key:
                        ref_count= len(set(allelic_count_ref[key]) - set(allelic_count_alt[key]))
                        alt_count= len(set(allelic_count_alt[key]) - set(allelic_count_ref[key]))
                        print('\t'.join(map(str,key))+'\t'+str(ref_count)+'\t'+str(alt_count))
                    allelic_count_ref= defaultdict(list)
                    allelic_count_alt= defaultdict(list)
                    pbar.update(1)                    
                    feed_list = []
        else:
            feed_list.append(window_list)
            p = mp.Pool(len(feed_list))
            results = p.starmap(window_mutation,
                    list(zip(feed_list,
                             repeat(bam_file), repeat(mq),repeat(baseQuality))))
            p.close()
            p.join()
            for read_dict_alt,read_dict_ref in results:
                if read_dict_alt is not None:
                    for key,val in read_dict_alt.items():
                        for read in val:
                            allelic_count_alt[key].append(read[0])
                if read_dict_ref is not None:
                    for key,val in read_dict_ref.items():
                        for read in val:
                            allelic_count_ref[key].append(read[0])
            all_key= set(list(allelic_count_ref.keys())+list(allelic_count_alt.keys()))
            for key in all_key:
                ref_count= len(set(allelic_count_ref[key]) - set(allelic_count_alt[key]))
                alt_count= len(set(allelic_count_alt[key]) - set(allelic_count_ref[key]))
                print('\t'.join(map(str,key))+'\t'+str(ref_count)+'\t'+str(alt_count))
            allelic_count_ref= defaultdict(list)
            allelic_count_alt= defaultdict(list)
            pbar.update(1)                    
            feed_list = []
    vcf.close()

if __name__ == '__main__':
    main()
