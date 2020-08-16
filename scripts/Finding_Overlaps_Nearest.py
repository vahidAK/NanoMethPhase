import tabix
import argparse
import os
import warnings
import gzip
import bz2 
from collections import defaultdict
parser = argparse.ArgumentParser(description="Mapping methylation results (or any genomic ranges) to different genomic regions in bed format.")
parser.add_argument("--notstrandlevel", "-nsl", action="store_true", required=False,help="Specify this option if you want to mapp not at strand level.")
parser.add_argument("--strandlevel_tss", "-sltss", action="store_true", required=False,help="Specify this option if you want to mapp at strand level and your annotation file has TSS"
                    " in  the 4th column and strand in fifth column. Usualy data created by web tool BioMart.")
parser.add_argument("--strandlevel", "-sl", action="store_true", required=False,help="Specify this option if you want to mapp at strand level and your annotation file has strand"
                    " in  the 4th column. This will consider the end as TSS when strand is negative and start as TSS when strand is positive")
parser.add_argument("--input", "-i", action="store", type=str,required=True,help="The path to the tab delimited methylation frequency file or any file you want to find its overlap to annotation file. file must be indexed by tabix "
                    "NOTE: Strand level information in this file is not considered for mapping.")
parser.add_argument("--annot", "-a", action="store", type=str,required=True,help="The path to the tab delimited bed file containing genomic regions."
                    "NOTE: The default is strand level annotation which also works for Upstream and downstream TSS. The format of annot file must be chrom\tstart\tend\tTSS\tstrand\tetc"
                    "If your file is not at strand level select --notstrandlevel option. you also can select ustream and downstram. The format of annot file must be chrom\tstart\tend\tetc")
#parser.add_argument("--not_aggregate", "-us", action="store_true", required=False,help="Do not aggregate multiple sites that mappe to the same annotation in a line, rather print in separate lines.")
parser.add_argument("--upstream", "-us", action="store", type=int,required=False, default=0,help="Include int bp upstreab of TSS. If you specify this it will include int upstream up to end of TSS/downstream option. Default is 0")
parser.add_argument("--downstream_tss", "-dstss", action="store", type=int,required=False,help="Include downstreab of TSS.if you include this if only consider from TSS/upstream option up to this point."
                    "Default is until the end of transcript. Note that this option overwrites the --downstream_transcript.")
parser.add_argument("--downstream", "-ds", action="store", type=int,default= 0,required=False,help="Include downstreab or after the end of gene/transcript/etc.if you include this if consider from TSS/upstream option up to this point after the end of gene/transcript. Default is zero")
parser.add_argument("--print_info", "-pi", action="store_true",help="Print overlapping information. This includes how many bases have overlap and overlap direction and the distance to original intersect (e.g if you select 1mb upstream this distance would be the distance to the original intersect)")
parser.add_argument("--distance_nearest", "-dn", action="store_true",help="If you wish to find nearest interval in annotation file to intervalls in the input select this option. Note: in this case annotation file must be cordinate sorted (sort -k1,1 -k2,2n) and indexed using tabix. In this case strand level information will not be used.")

args = parser.parse_args()

meth= os.path.abspath(args.input)
annot= os.path.abspath(args.annot)
upstream= args.upstream
downstream= args.downstream

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


if not args.distance_nearest:
    try:
        tb = tabix.open(meth)
    except:
        raise Exception("Methylation file is not indexed. please first idex the file bt tabix")
    with openfile(meth) as a:
        header_input= next(a).rstrip()

    with openfile(annot) as ann:
        head_ann= next(ann).rstrip()
        if not args.print_info:
            print(head_ann+'\t'+header_input)
        else:
            print(head_ann+'\tOvelapWith\tOverlap(bp)\tOverlap(precent)\tDistance\t'+header_input)
        for line in ann:
            line=line.rstrip()
            if args.strandlevel_tss:
                chrom,start_trans,end_trans, TSS, strand= line.split('\t')[0:5]
                if strand == "-" or strand == "-1":
                    original_start= int(start_trans)
                    original_end= int(TSS)
                    end= int(TSS) + upstream
                    start= int(start_trans) - downstream
                    if args.downstream_tss:
                        start= int(TSS) - args.downstream_tss # it will overwrite int(start_trans) - downstream
                elif strand == "+" or strand == "1":
                    original_start= int(TSS)
                    original_end= int(end_trans)
                    start= int(TSS) - upstream
                    end= int(end_trans) + downstream
                    if args.downstream_tss:
                        end= int(TSS) + args.downstream_tss
                else:
                    raise Exception("Strand info not in 5th column.")
                    
            elif args.strandlevel:
                chrom,start_trans,end_trans, strand= line.split('\t')[0:4]
                if strand == "-" or strand == "-1":
                    TSS= end_trans
                    original_start= int(start_trans)
                    original_end= int(TSS)
                    end= int(TSS) + upstream
                    start= int(start_trans) - downstream
                    if args.downstream_tss:
                        start= int(TSS) - args.downstream_tss # it will overwrite int(start_trans) - downstream
                elif strand == "+" or strand == "1":
                    TSS= start_trans
                    original_start= int(TSS)
                    original_end= int(end_trans)
                    start= int(TSS) - upstream
                    end= int(end_trans) + downstream
                    if args.downstream_tss:
                        end= int(TSS) + args.downstream_tss
                else:
                    raise Exception("Strand info not in 4th column.")
                    
            elif args.notstrandlevel:
                chrom,start_gene,end_gene= line.split('\t')[0:3]
                original_start= int(start_gene)
                original_end= int(end_gene)
                start= int(start_gene)
                end= int(end_gene)
                start= start - upstream
                end= end + downstream
            else:
                raise Exception("You did not specify any strand level mapping options.")
    
            if 'chr' not in chrom.lower():
                chrom= 'chr'+chrom
            try:
                try:
                    records = tb.query(chrom, start, end)
                except:
                    records = tb.query(chrom[3:], start, end)
            except:
                warnings.warn("{}:{}-{} does not exist in the Methylatin file. Skipping the line".format(line.split('\t')[0],line.split('\t')[1],line.split('\t')[2]))
                continue
            for record in records:
                if not args.print_info:
                    print(line+'\t'+'\t'.join(record))
                else:
                    start_file,end_file= record[1:3]
                    start_file,end_file= int(start_file),int(end_file) 
                    if ((start_file <= original_start and end_file >= original_start) or
                        (start_file >= original_start and start_file <= original_end)):
                        distance= 0
                    elif end_file < original_start:
                        distance= original_start - end_file
                    elif start_file > original_end:
                        distance= start_file - original_end
                    if start_file >= start and end_file <= end:
                        overlapbp= end_file - start_file
                        overlappr= 100
                        overlapsite= "all"
                    elif start_file <= start and end_file >= end:
                        overlapbp= end - start
                        overlappr= 100
                        overlapsite= "all"
                    elif start_file <= start and end_file <= end:
                        overlapbp= end_file - start
                        overlappr= round((overlapbp/(end_file - start_file))*100,1)
                        overlapsite= "start"
                    elif start_file >= start and end_file >= end:
                        overlapbp= end - start_file
                        overlappr= round((overlapbp/(end_file - start_file))*100,1)
                        overlapsite= "end"
                    print(line+'\t'+'\t'.join(map(str,[overlapsite,
                                                  overlapbp,
                                                  overlappr,
                                                  distance]))+'\t'+
                            '\t'.join(record))
        
else:
    with openfile(annot) as b:
        header_annot= next(b).rstrip()
    try:
        tb = tabix.open(annot)
    except:
        raise Exception("Methylation file is not indexed. please first idex the file bt tabix")
    with openfile(meth) as a:
        header_input= next(a).rstrip()
        print(header_input+'\tdistance\t'+header_annot)
        for line in a:
            chrom,start,end= line.split('\t')[0:3]
            start,end= int(start), int(end)
            if 'chr' not in chrom.lower():
                chrom= 'chr'+chrom
            try:
                try:
                    records = tb.query(chrom, start, start+1000000000)
                except:
                    records = tb.query(chrom[3:], start, start+1000000000)
            except:
                warnings.warn("{}:{} does not exist in the Methylatin file. Skipping the line".format(line.split('\t')[0],line.split('\t')[1]))
                continue
            nearest_interval= ""
            distance= None
            for record in records:
                nearest_interval= record
                nearest_start, nearest_end= int(record[1]),int(record[2])
                distance= nearest_start - end
                break
            if distance is not None and distance <= 0:
                try:
                    records = tb.query(chrom, start, end)
                except:
                    records = tb.query(chrom[3:], start, end)
                distance= 0
                for record in records:
                    print(line.rstrip(),str(distance),'\t'.join(nearest_interval), sep='\t')
            else: #Check for upstream
                if distance is None:
                    try:
                        upstream_records = tb.query(chrom, 0, start)
                    except:
                        upstream_records = tb.query(chrom[3:], 0, start)
                    for upstream_record in upstream_records:
                        nearest_upstream_interval= upstream_record
                        upstream_distance= start - int(upstream_record[2])
                        if upstream_distance < 0:
                            upstream_distance = 0
                    print(line.rstrip(),str(upstream_distance),'\t'.join(nearest_upstream_interval), sep='\t')
                else:
                    check_start= start - distance
                    if check_start < 0:
                        check_start= 0
                    try:
                        upstream_records = tb.query(chrom, check_start, start)
                    except:
                        upstream_records = tb.query(chrom[3:], check_start, start)
                    upstream_distance = distance+1
                    for upstream_record in upstream_records:
                        nearest_upstream_interval= upstream_record
                        upstream_distance= start - int(upstream_record[2])
                        if upstream_distance < 0:
                            upstream_distance = 0
                    if upstream_distance < distance:
                        print(line.rstrip(),str(upstream_distance),'\t'.join(nearest_upstream_interval), sep='\t')
                    else:
                        print(line.rstrip(),str(distance),'\t'.join(nearest_interval), sep='\t')

