import argparse
import gzip
import bz2 
import os
import tabix
import warnings
import numpy as np
from statistics import mean
parser = argparse.ArgumentParser(description="")
parser.add_argument("--input", "-i", action="store", required=True,type=str,help="Input file"
                    ". The first thre columns are cordinate and forth column is methylation freq. If multiple sample"
                    " merge files to have frequency on the forth column and after that for each sample.")
parser.add_argument("--ignore_intervals", "-ii", action="store", required=True,type=str,help="Input file"
                    " with intervals you want to ignore in bed format and indexed by tabix")
parser.add_argument("--dmr_file", "-d", action="store", required=True,type=str,help="")
parser.add_argument("--out_prefix", "-o", action="store", required=True,type=str,help="")
parser.add_argument("--method", "-m", action="store", required=False,default="freq",type=str,help="freq\num. "
                    "Either report methyllation frequency of CpGs (freq) or report number of CpG more than and less that a ratio"
                    ". Default is freq")
parser.add_argument("--thresholds", "-t", action="store", required=False,default="0.3,0.7",type=str,help="If you have selected"
                    "num method give the lower and upper bound thresholds (comma separated)."
                    ". Default is 0.3,0.7")
parser.add_argument("--distances", "-ds", action="store", required=False,default="5000,500000",type=str,help=""
                    "Distance from interval to grab adjacent sites. Default is 5000,500000."
                    "Which means grab CpGs from 5000 to 500000bp distance from interval.")


args = parser.parse_args()
input_file= os.path.abspath(args.input)
ignore_intervals= os.path.abspath(args.ignore_intervals)
dmr_file= os.path.abspath(args.dmr_file)
out_prefix= os.path.abspath(args.out_prefix)
method= args.method
if method != "freq" or method != "num":
    Exception("Enter method correctly freq/num")
thresholds= sorted([float(args.thresholds.split(',')[0]),
                    float(args.thresholds.split(',')[1])])

distances= sorted([int(args.distances.split(',')[0]),
                    int(args.distances.split(',')[1])])

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
tb = tabix.open(input_file)
tb_ignore= tabix.open(ignore_intervals)
with openfile(input_file) as a:
    samples= next(a).rstrip().split('\t')[3:]

with openfile(dmr_file) as dmr:
    next(dmr)
    header=("chr_start_end\t"+'\t'.join(i+"_DMR" for i in samples)+'\t'+
                                        '\t'.join(i+"_Before" for i in samples)+'\t'+
                                            '\t'.join(i+"_After" for i in samples))
    if method == "num":
        out1= open(out_prefix+"_BeforeAndAfter.tsv",'w')
        out2= open(out_prefix+"_BeforeAndAfter_Mean.tsv",'w')
        out1.write(header+'\n')
        out2.write("chr_start_end\tMeanAtDMR\tMeanAtAdjacent\n")
    if method == "freq":
        out1= open(out_prefix+"_BeforeAndAfter.tsv",'w')
        out1.write(header+'\n')
    for line in dmr:
        line= line.rstrip().split("\t")
        chrom= line[0]
        start= int(line[1])
        end= int(line[2])
        before_start= start-distances[1]
        before_end= start-distances[0]
        if before_start < 0:
            before_start= 0
        if before_end < 0:
            before_end= 0
        after_start= start+distances[0]
        after_end= start+distances[1]
        if before_end >= 0:
            try:
                while len(list(tb_ignore.query(chrom, before_start, before_end))) > 0:
                    before_start -= distances[0]
                    before_end -= distances[0]
                    if before_start < 0:
                        before_start= 0
                    if before_end < 0:
                        before_end= 0
                    pass
            except:
                pass
        try:
            while len(list(tb_ignore.query(chrom, after_start, after_end))) > 0:
                after_start += distances[0]
                after_end += distances[0]
        except:
            pass
#        if method == "freq":
#            out1= open(out_prefix+"_BeforeAndAfter_"+chrom+"_"+str(start)+"-"+str(end)+".tsv",'w')
#            out1.write(header+'\n')
        try:
            records_dmr= list(tb.query(chrom, start, end))
            records_before= list(tb.query(chrom, before_start, before_end))
            records_before= records_before[-len(records_dmr):]
            records_after= list(tb.query(chrom, after_start, after_end))
            records_after= records_after[0:len(records_dmr)]
        except:
            warnings.warn("No CpG at {} interval or adjacent site. Skipping it".format('_'.join(line[0:3])))
            continue
        if len(records_dmr) == 0:
            warnings.warn("No CpG at {}. Skipping it".format('_'.join(line[0:3])))
            continue

        elif len(records_dmr) == len(records_before) == len(records_after):
            adjacent= np.append(np.asarray(records_before),np.asarray(records_after)[:,3:], axis=1)
        elif len(records_dmr) == len(records_before) != len(records_after):
            records_before= list(tb.query(chrom, before_start, before_end))
            records_before= records_before[-((2*(len(records_dmr)))-3):]
            if (2*(len(records_dmr)))-3 == len(records_before):
                adjacent= np.asarray(records_before)
            else:
                warnings.warn("Not enough CpG at adjacent for {}. Skipping it".format('_'.join(line[0:3])))
                
        elif len(records_dmr) == len(records_after) != len(records_before):
            records_after= list(tb.query(chrom, after_start, after_end))
            records_after= records_after[0:(2*(len(records_dmr)))-3]
            if (2*(len(records_dmr)))-3 == len(records_after):
                adjacent= np.asarray(records_after)
            else:
                warnings.warn("Not enough CpG at adjacent for {}. Skipping it".format('_'.join(line[0:3])))
        else:
            warnings.warn("Not enough CpG at adjacent for {}. Skipping it".format('_'.join(line[0:3])))
#        if len(records_dmr) != len(records_before) or len(records_dmr) != len(records_after):
#            warnings.warn("Afte or before number of CpGs are not equl to DMRs CpGs at {}".format(chrom+":"+str(start)+"-"+str(end)))
        records_dmr= np.asarray(records_dmr)          
        if method == "freq":
#            for i,j in zip(np.asarray(records_dmr),adjacent):
            mean_dmr= np.nanmean(np.asarray(records_dmr[:,3:],dtype=np.float64),axis = 0)
            mean_adjacent= np.nanmean(np.asarray(adjacent[:,3:],dtype=np.float64),axis = 0)
            out1.write(':'.join(line[0:3])+'\t'+'\t'.join(map(str,mean_dmr))+'\t'+'\t'.join(map(str,mean_adjacent))+'\n')
#                out1.write('_'.join(i[0:3])+'\t'+'\t'.join(i[3:])+'\t'+'\t'.join(j[3:])+'\n')
#            out1.close()
        else:
            all_cpgs_DMR= np.count_nonzero(np.asarray(records_dmr[:,3:],dtype=np.float64) >= 0,axis = 0)
            samples_num_DMR= ((all_cpgs_DMR - (
                              np.count_nonzero(np.asarray(records_dmr[:,3:],dtype=np.float64) < thresholds[0],axis = 0) +
                              np.count_nonzero(np.asarray(records_dmr[:,3:],dtype=np.float64) > thresholds[1],axis = 0))) / len(records_dmr))
                              
            all_cpgs_Adjacent= np.count_nonzero(np.asarray(adjacent[:,3:],dtype=np.float64) >= 0,axis = 0)
            samples_num_Adjacent= ((all_cpgs_Adjacent - (
                          np.count_nonzero(np.asarray(adjacent[:,3:],dtype=np.float64) < thresholds[0],axis = 0) +
                          np.count_nonzero(np.asarray(adjacent[:,3:],dtype=np.float64) > thresholds[1],axis = 0))) / len(records_dmr))
            out1.write(':'.join(line[0:3])+'\t'+'\t'.join(map(str,samples_num_DMR))+'\t'+'\t'.join(map(str,samples_num_Adjacent))+'\n')
            out2.write(':'.join(line[0:3])+'\t'+ str(mean(list(samples_num_DMR)))
            +'\t'+ str(mean(list(samples_num_Adjacent)))+'\n')
                
    out1.close()
    if method=="num":
        out2.close()

