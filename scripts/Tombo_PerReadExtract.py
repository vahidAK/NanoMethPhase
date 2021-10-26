from tombo import tombo_helper, tombo_stats, resquiggle
import argparse
import os

parser = argparse.ArgumentParser(description="Tombo per read stat extraction")
parser.add_argument("--input", "-i", action="store", type=str, required=True, 
                    help="Path to tombo per-read stat hdf5 file")
parser.add_argument("--region_file", "-r", action="store", type=str, required=False,
                    default=None, help="Path to a bed like file that includes region you wish to "
                    " extract per read stat. The three first column must be "
                    "chromosome\tstart\tend and file must have NO HEADER."
                    " By defaut this script will extract "
                    "per read for chr1 to chr22 plus chrX and chrY.")
args = parser.parse_args()

file= os.path.abspath(args.input)
regions= list()
if args.region_file is None:
    regions= [('chr'+str(x),0,500000000) 
                  for x in range(1,23)] + [('chrX',0,500000000),('chrY',0,500000000)] 
else:
    region_file= os.path.abspath(args.region_file)
    with open(region_file) as bed:
        for line in bed:
            line=line.rstrip().split('\t')
            regions.append((line[0],int(line[1]),int(line[2])))


sample_per_read_stats = tombo_stats.PerReadStats(file)
for region in regions:
    chrom,start,end= region
    reg_data = tombo_helper.intervalData(chrm=chrom, start=start, end=end, strand='+')
    reg_per_read_stats = sample_per_read_stats.get_region_per_read_stats(reg_data)
    if reg_per_read_stats is not None:
        for i in reg_per_read_stats:
            print(chrom+'\t'+'+'+'\t'+'\t'.join(map(str,i)))
    reg_data = tombo_helper.intervalData(chrm=chrom, start=start, end=end, strand='-')
    reg_per_read_stats = sample_per_read_stats.get_region_per_read_stats(reg_data)
    if reg_per_read_stats is not None:
        for i in reg_per_read_stats:
            print(chrom+'\t'+'-'+'\t'+'\t'.join(map(str,i)))

sample_per_read_stats.close()
