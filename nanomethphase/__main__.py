#! /usr/bin/env python3
# coding=utf-8

# Copyright (C) 2020  Vahid Akbari

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
NanoMethPhase: Phasing DNA methylation data using Oxford Nanopore
Technology sequencing.
"""

__author__ = "Vahid Akbari"
__email__ = "vakbari@bcgsc.ca"
__copyright__ = "Copyright (C) 2020, " + __author__
__license__ = "GPLv3"
__collaborator__ = "Jean-Michel Garant"

import os
import re
import glob
import gzip
import bz2
import argparse
import warnings
import textwrap
import subprocess
import multiprocessing as mp
import sys
from collections import defaultdict
from itertools import repeat
import statistics
import pysam
import tabix
import math
from tqdm import tqdm

def getChromsFromBAM(filename):
    chroms = []
    stats = pysam.idxstats(filename)
    for row in stats.split("\n"):
        fields = row.split("\t")
        if fields[0] != '*' and fields[0] != '':
            chroms.append(fields[0])
    return chroms


def get_base_info(feed_list,
                  alignment_file,
                  chrom):
    read_HP_list= list()
    samfile = pysam.AlignmentFile(alignment_file, 'rb')
    for info in feed_list:
        HP,position,ref,alt= info
        try:
            sam_pileup= samfile.pileup(chrom,
                                       position,
                                       position+1,
                                       truncate=True)
        except:#The cordiniate is not found so ignore it
            warnings.warn("Variant {} {} did not find or do not have any map "
                          "reads in the alignment file. Check if correct"
                          "bam file is given or bam is indexed and not corrupted."
                          " Skipping it.".format(chrom,
                                                   position+1))
            continue
        for pileupcolumn in sam_pileup:
            pileupcolumn.set_min_base_quality(0)
            if pileupcolumn.pos == position:
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del or pileupread.is_refskip:
                        continue
                    read_id= pileupread.alignment.query_name
                    read_mq= pileupread.alignment.mapping_quality
                    read_id= pileupread.alignment.query_name
                    read_len= pileupread.alignment.query_alignment_length
                    read_start= pileupread.alignment.reference_start
                    read_end= pileupread.alignment.reference_end
                    flag= pileupread.alignment.flag
                    phred= pileupread.alignment.query_qualities[
                                                pileupread.query_position]
                    
                    if pileupread.alignment.is_supplementary:
                        suppl_flag= str(flag)+":yes"
                    else:
                        suppl_flag= str(flag)+":no"
                    if pileupread.alignment.is_reverse:
                        strand = "-"
                    else:
                        strand = "+"
                    key_per_read = (chrom,read_start,read_end,
                                    read_id,strand,suppl_flag,
                                    str(read_len)+":"+str(read_mq))
                    val= [position,phred]
                    read_base= None
                    if HP == '1|0' or HP == '0|1':
                        if len(ref) == 1 and len(alt)==1: # dealing with mismatches
                            read_base= pileupread.alignment.query_sequence[
                                    pileupread.query_position]
                        elif len(ref) > 1 and len(alt) == 1: #dealing with deletions
                            if pileupread.indel < 0 and abs(pileupread.indel) == len(ref) -1:
                                read_base= pileupread.alignment.query_sequence[
                                        pileupread.query_position]
                            elif pileupread.indel == 0:
                                read_base = pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(ref)]
                            else:
                                continue
                        elif len(ref) == 1 and len(alt) > 1: #dealing with insertions
                            if pileupread.indel > 0 and pileupread.indel == len(alt) -1:
                                read_base= pileupread.alignment.query_sequence[
                                                pileupread.query_position:pileupread.query_position+len(alt)]
                            elif pileupread.indel == 0:
                                read_base= pileupread.alignment.query_sequence[
                                        pileupread.query_position]
                            else:
                                continue
                                
                    if read_base is not None:
                        read_base= read_base.upper()
                        if HP == '1|0' and read_base == alt:
                            read_HP_list.append([(*key_per_read,"HP1"),
                                        ':'.join(map(str,val + [alt]))])
                        elif HP == '1|0' and read_base == ref:
                            read_HP_list.append([(*key_per_read,"HP2"),
                                        ':'.join(map(str,val + [ref]))])
                        elif HP == '0|1' and read_base == alt:
                            read_HP_list.append([(*key_per_read,"HP2"),
                                        ':'.join(map(str,val + [alt]))])
                        elif HP == '0|1' and read_base == ref:
                            read_HP_list.append([(*key_per_read,"HP1"),
                                        ':'.join(map(str,val + [ref]))])
    samfile.close()
    return read_HP_list



def bam_info_extractor(read,
                       reference,
                       fasta):
    if read.is_reverse:
        strand = "-"
    else:
        strand = "+"
    read_id = read.query_name
    start = read.reference_start
    end = read.reference_end
    true_ref_name = read.reference_name
    rnext= read.next_reference_name
    pnext= read.next_reference_start
    tlen= read.template_length
    cigar = read.cigartuples
    base_qualities = read.query_qualities
    flag = read.flag
    if read.query_sequence:
        read_seq = read.query_sequence
    ref_seq = ""
    ref_len = ""
    if reference is not None:
        try:
            ref_seq = fasta.fetch(reference=true_ref_name,
                                  start=start,
                                  end=end)
        except:
                warnings.warn("Reference genome sequence was not found "
                              "for this read: {} at this cordinates {}:{}-{}. "
                              "Skipping the read".format(read_id,
                                                         true_ref_name,
                                                         start,
                                                         end))
    if ((read_seq and cigar and base_qualities) and
    (cigar != "*" or cigar is not None) and
                                      base_qualities is not None):
        read_seq = read_seq.upper()
        read_len = read.query_alignment_length
        ref_seq= ref_seq.upper()
        ref_len= len(ref_seq)
        all_tags= read.get_tags(with_value_type=True)
        return (true_ref_name, strand, flag, read_id, read_seq ,
                read_len, cigar, rnext, pnext, tlen, base_qualities ,
                start, end, ref_seq, ref_len, all_tags)
    else:
        warnings.warn("{} does not have a read sequence,CIGAR"
                      ", or base quality information. "
                      "Skipping the read".format(read_id))


def get_indel(vcf):
    indels= set()
    with openfile(vcf) as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue
            line= line.rstrip().split('\t')
            if (line[9].startswith("0|1") or
                line[9].startswith("1|0") or
                line[9].startswith("0/1") or
                line[9].startswith("1/0")):
                if len(line[3]) > 1 or len(line[4]) > 1:
                    indels.add((line[0],str(int(line[1])-1)))
    return indels


def outexist_phase(overwrite,
                   out1,
                   out2):
    '''
    Checks if selected output for phase module already exist
    and return error if exist and args.overwrite did not specify
    '''
    if (not overwrite
            and (os.path.isfile(out1+".bam") or os.path.isfile(out2+".bam"))):
        raise FileExistsError("The selected output bam files [{}, {}] already "
                              "exist. Select --overwrite option if you want "
                              "to overwrite them".format(out1+".bam",
                                                         out2+".bam"))
    if (not overwrite
            and (os.path.isfile(out1+"_Converted2Bisulfite.bam")
                 or os.path.isfile(out2+"_Converted2Bisulfite.bam"))):
        raise FileExistsError("The selected output bam2bis files [{}, {}] "
                              "already exist. Select --overwrite option if "
                              "you want to overwrite "
                              "them".format(out1+"_Converted2Bisulfite.bam",
                                            out2+"_Converted2Bisulfite.bam"))
    if (not overwrite
            and (os.path.isfile(out1+"_MethylCall.tsv")
                 or os.path.isfile(out2+"_MethylCall.tsv"))):
        raise FileExistsError("The selected output MethylCall files [{}, {}] "
                              "already exist. Select --overwrite option if "
                              "you want to overwrite "
                              "them".format(out1+"_MethylCall.tsv",
                                            out2+"_MethylCall.tsv"))
    if (not overwrite
            and (os.path.isfile(out1+"_MethylFrequency.tsv")
                 or os.path.isfile(out2+"_MethylFrequency.tsv"))):
        raise FileExistsError("The selected output Methylation Frequency "
                              "files [{}, {}] already exist. Select "
                              "--overwrite option if you want to overwrite "
                              "them".format(out1+"_MethylFrequency.tsv",
                                            out2+"_MethylFrequency.tsv"))
    if (not overwrite
            and (os.path.isfile(out1+"_HP2_PerReadInfo.tsv"))):
        raise FileExistsError("The selected output per read info file "
                              "{} already exists. Select "
                              "--overwrite option if you want to overwrite "
                              "it".format(out1+"_HP2_MethylFrequency.tsv"))


def outformats_phase(outformat,reference,MethylCallfile):
    '''
    Make a list of selected outputs by user and show a massage to tell which
    output format will be created and also reads appropreate files.
    Also returns error if the output format
    specified wrongly.
    '''
    selected_outformat= []
    fasta=tb= ""

    if 'bam2bis' in outformat and (reference is None or MethylCallfile is None):
        raise TypeError("You have selected bam2bis output format but you did "
                        "not specify reference file and/or processed"
                        " MethylCallfile. Reference file must be also "
                        "indexed by samtools faidx.")
    elif ('bam2bis' in outformat and
          reference is not None and
          MethylCallfile is not None):
        if not os.path.isfile(os.path.abspath(MethylCallfile)+".tbi"):
            raise Exception("Could not find index file for methylation call file.")
        tb = tabix.open(os.path.abspath(MethylCallfile))
    if reference is not None:
        try:
            reference = os.path.abspath(reference)
            fasta = pysam.FastaFile(reference)
        except:
            raise Exception("Cannot load reference file.")

    if 'methylcall' in outformat and MethylCallfile is None:
        raise TypeError("You have selected methylcall output format but you "
                        "did not specify methylation call file from "
                        "methyl_call_processor module.")
    elif 'methylcall' in outformat and MethylCallfile is not None:
        MethylCallfile = os.path.abspath(MethylCallfile)
        if not os.path.isfile(os.path.abspath(MethylCallfile)+".tbi"):
            raise Exception("Could not find index file for methylation call file.")
        tb = tabix.open(os.path.abspath(MethylCallfile))
    if 'bam' in outformat:
        selected_outformat.append('bam')
    if 'bam2bis' in outformat:
        selected_outformat.append("bam2bis")
    if 'methylcall' in outformat:
        selected_outformat.append("methylcall")
    if selected_outformat:
        sys.stderr.write("NanoMethPhase selected output format(s):"
                         " {}\n".format(','.join(selected_outformat)))
    else:
        raise TypeError("Please select appropriate output format: "
                        "bam,bam2bis,methylcall. Multiple formats must be "
                        "seperated by comma.")
    return fasta, tb


def outputs(outformat,
            bam,
            out1,
            out2):
    '''
    Creates output file based on selected output formats by user
    '''
    outHP1Sam = outHP2Sam = outHP12BisSam = ""
    outHP22BisSam = outCall1 = outCall2 = outFreq1 = outFreq2 = ""
    if 'bam' in outformat:
        outHP1Sam = pysam.AlignmentFile(out1+".bam", "wb", template=bam)
        outHP2Sam = pysam.AlignmentFile(out2+".bam", "wb", template=bam)
    if 'bam2bis' in outformat:
        outHP12BisSam = pysam.AlignmentFile(out1+"_Converted2Bisulfite.bam",
                                            "wb", template=bam)
        outHP22BisSam = pysam.AlignmentFile(out2+"_Converted2Bisulfite.bam",
                                            "wb", template=bam)
    if 'methylcall' in outformat:
        outCall1 = open(out1+"_MethylCall.tsv", 'w')
        outCall2 = open(out2+"_MethylCall.tsv", 'w')
        outFreq1 = open(out1+"_MethylFrequency.tsv", 'w')
        outFreq2 = open(out2+"_MethylFrequency.tsv", 'w')
        headerCall = "\t".join(["chromosome",
                                "start",
                                "end",
                                "strand",
                                "read_name",
                                "llr_Or_DeltaProb\n"])
        outCall1.write(headerCall)
        outCall2.write(headerCall)
        headerFreq = "\t".join(["chromosome",
                                "start",
                                "end",
                                "strand",
                                "NumOfAllCalls",
                                "NumOfModCalls",
                                "MethylFreq\n"])
        outFreq1.write(headerFreq)
        outFreq2.write(headerFreq)

    return (outHP1Sam, outHP2Sam,outHP12BisSam,outHP22BisSam,
            outCall1,outCall2,outFreq1,outFreq2)


def openalignment(alignment_file,
                  window):
    '''
    Opens a bam/sam file and creates bam iterator
    '''
    bam = pysam.AlignmentFile(alignment_file, 'rb')
    if window is not None:
        window_chrom = window.split(':')[0]
        if len(window.split(':')) == 2:
            window_margin= window.split(':')[1].split('-')
            if len(window_margin) == 2:
                window_start = int(window_margin[0])
                window_end = int(window_margin[1])
                bamiter = bam.fetch(window_chrom, window_start, window_end)
                count= bam.count(window_chrom, window_start, window_end)
            else:
                window_start = int(window_margin[0])
                bamiter = bam.fetch(window_chrom, window_start)
                count= bam.count(window_chrom, window_start)
        else:
            try:
                bamiter = bam.fetch(window_chrom)
                count= bam.count(window_chrom)
            except:
                warnings.warn("{} is not in the alignment "
                              "file Or alignment is truncated or corrupt indexed."
                              " Skipping it".format(window_chrom))
                count= 0
                bamiter= ""
    else:
        bamiter = bam.fetch(until_eof=True)
        count = None
    return bamiter, bam, count


def alignmentwriter(result,
                    output):
    '''
    Writes the results of converting reads to bisulfite format
    to a bam file
    '''
    (HP, read_id, flag, true_ref_name ,
     start, mp_quality, cigar, RNEXT ,
     PNEXT, TLEN, ref_seq, QUAL, all_tags) = result
    out_samRead = pysam.AlignedSegment(output.header)
    out_samRead.query_name = read_id
    out_samRead.cigarstring = str(len(ref_seq))+'M'
    out_samRead.query_sequence = ref_seq
    out_samRead.flag = flag
    out_samRead.reference_name = true_ref_name
    out_samRead.reference_start = start
    out_samRead.mapping_quality = mp_quality
    if HP != 'NON':
        all_tags= [(HP[0:2], int(HP[-1]),"i")]+all_tags
    if len(all_tags) >= 1:
        out_samRead.set_tags(all_tags)
    output.write(out_samRead)


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


def FrequencyCalculator(file_path):
    """
    Calculating methylation frequency for each phased methylation call
    file. The output methylation frequency file include fractional
    methylation which can be used for differential methylation analysis
    and detection of differentially methylated regions (DMR)
    """
    dict_mod = defaultdict(int)
    dict_all = defaultdict(int)
    with open(file_path) as input_file:
        next(input_file)  # To exclude header
        for line in input_file:
            line = line.rstrip().split('\t')
            if line[3] == '+':
                key = tuple(line[0:4])
            else:
                key = tuple([line[0],
                                str(int(line[1])+1),
                                str(int(line[2])+1),
                                line[3]])
            dict_all[key] += 1
            if float(line[5]) > 0:
                dict_mod[key] += 1
    return dict_mod, dict_all


def methcall2bed(readlist,
                 motif,
                 callthresh,
                 tool):
    """
    This function converts methylation call file to a bed
    format and also splits multi-group CpG sites, in case of nanopolish, to single group.
    """
    read_list= list()
    if tool == "nanopolish":
        for read in readlist:
            methylated_sites = []
            unmethylated_sites= []
            llr_methylated= []
            llr_unmethylated= []
            for line in read:
                line = line.split('\t')
                num_sites = int(line[9])
                logratio = float(line[5])/int(line[9])
                sequence= line[10].upper()
                strand= line[1]
                cpg_pos= int(line[2])
                read_id= line[4]
                # Skipping ambiguous call in methyl call file
                if abs(logratio) < callthresh:
                    continue
                chrom = line[0]
                if int(line[9]) > 1:  # Check if the line includes multi-group CpG
                    splited_groupIndexes = [(j.start())
                                            for j in re.finditer(motif, sequence)]
                    if logratio > 0:
                        methylated_sites.append(cpg_pos)
                        llr_methylated.append(str(logratio))
                    else:
                        unmethylated_sites.append(cpg_pos)
                        llr_unmethylated.append(str(logratio))
                    for splited_groupIndex in splited_groupIndexes[1:]:
                        position = cpg_pos + (splited_groupIndex
                                                 - splited_groupIndexes[0])
                        if logratio > 0:
                            methylated_sites.append(position)
                            llr_methylated.append(str(logratio))
                        else:
                            unmethylated_sites.append(position)
                            llr_unmethylated.append(str(logratio))
                else:
                    if logratio > 0:
                        methylated_sites.append(cpg_pos)
                        llr_methylated.append(str(logratio))
                    else:
                        unmethylated_sites.append(cpg_pos)
                        llr_unmethylated.append(str(logratio))
            all_positions= sorted(methylated_sites + unmethylated_sites)
            if all_positions:
                if not methylated_sites:
                    methylated_sites.append('NA')
                    llr_methylated.append('NA')
                if not unmethylated_sites:
                    unmethylated_sites.append('NA')
                    llr_unmethylated.append('NA')
                append_info= (chrom,
                              str(all_positions[0]),
                              str(all_positions[-1]+1),
                                strand, read_id,
                                ','.join(llr_methylated),
                                ','.join(llr_unmethylated),
                                ','.join(map(str,methylated_sites)),
                                ','.join(map(str,unmethylated_sites)))
                read_list.append(append_info)

    elif tool == "deepsignal":
        for read in readlist:
            methylated_sites = []
            unmethylated_sites= []
            prob_methylated= []
            prob_unmethylated= []
            for line in read:
                line = line.split('\t')
                strand= line[2]
                cpg_pos= int(line[1])
                if strand == "-":
                    cpg_pos= cpg_pos - 1 # to make it like nanopolish
                read_id= line[4]
                chrom = line[0]
                deltaprob = float(line[7]) - float(line[6])
                # Skipping ambiguous call in methyl call file
                if abs(deltaprob) < callthresh:
                    continue
                if deltaprob > 0:
                    methylated_sites.append(cpg_pos)
                    prob_methylated.append(str(deltaprob))
                else:
                    unmethylated_sites.append(cpg_pos)
                    prob_unmethylated.append(str(deltaprob))
            all_positions= sorted(methylated_sites + unmethylated_sites)
            if all_positions:
                if not methylated_sites:
                    methylated_sites.append('NA')
                    prob_methylated.append('NA')
                if not unmethylated_sites:
                    unmethylated_sites.append('NA')
                    prob_unmethylated.append('NA')
                append_info= (chrom,
                              str(all_positions[0]),
                              str(all_positions[-1]+1),
                              strand, read_id,
                              ','.join(prob_methylated),
                              ','.join(prob_unmethylated),
                              ','.join(map(str,methylated_sites)),
                              ','.join(map(str,unmethylated_sites)))
                read_list.append(append_info)
        
    elif tool == "megalodon":
        for read in readlist:
            methylated_sites = []
            unmethylated_sites= []
            prob_methylated= []
            prob_unmethylated= []
            for line in read:
                line = line.split('\t')
                strand= line[2]
                if strand == "1":
                    strand = "+"
                elif strand == "-1":
                    strand = "-"
                cpg_pos= int(line[3])
                if strand == "-":
                    cpg_pos= cpg_pos - 1 # to make it like nanopolish
                read_id= line[0]
                chrom = line[1]
                deltaprob= math.exp(float(line[4])) - (1 - math.exp(float(line[4])))
                # Skipping ambiguous call in methyl call file
                if abs(deltaprob) < callthresh:
                    continue
                if deltaprob > 0:
                    methylated_sites.append(cpg_pos)
                    prob_methylated.append(str(deltaprob))
                else:
                    unmethylated_sites.append(cpg_pos)
                    prob_unmethylated.append(str(deltaprob))
            all_positions= sorted(methylated_sites + unmethylated_sites)
            if all_positions:
                if not methylated_sites:
                    methylated_sites.append('NA')
                    prob_methylated.append('NA')
                if not unmethylated_sites:
                    unmethylated_sites.append('NA')
                    prob_unmethylated.append('NA')
                append_info= (chrom,
                              str(all_positions[0]),
                              str(all_positions[-1]+1),
                              strand, read_id,
                              ','.join(prob_methylated),
                              ','.join(prob_unmethylated),
                              ','.join(map(str,methylated_sites)),
                              ','.join(map(str,unmethylated_sites)))
                read_list.append(append_info)
                
    return read_list


def read2bis(read_sam_list):
    """
    This function converts a read based on information in processed
    MethylCallFile to a bisulfite format read for nice visualization by
    IGV.
    """
    motif = 'CG'
    sam_list = read_sam_list[2:]
    ref_seq = sam_list[-5]
    strand = read_sam_list[1]
    HP = read_sam_list[0]
    all_methylated = sam_list[-2]
    all_unmethylated = sam_list[-1]
    all_tags= sam_list[-3]
    all_sites = [(j.start()) for j in re.finditer(motif, ref_seq)]
    ref_seq = list(ref_seq)
    if strand == '-' and motif == 'CG':
        offset = 1
        ambigbase = 'N'
        unmodified = 'A'
    else:
        offset = 0
        ambigbase = 'N'
        unmodified = 'T'
    for site in all_sites:
        if site not in all_methylated:
            if site+offset < len(ref_seq):
                ref_seq[site+offset] = ambigbase

    for site in all_unmethylated:
        if site+offset < len(ref_seq):
            ref_seq[site+offset] = unmodified
    ref_seq = ''.join(ref_seq)
    return [HP]+sam_list[0:-5]+[ref_seq]+[sam_list[-4]]+[all_tags]


def vcf2dict_phase(vcf, 
                   window):
    """
    This function converts the input vcf file to haplotype1 and
    haplotype2 dictionaries to be used for read phasing.
    """
    vcf_file = openfile(vcf)
    vcf_dict= defaultdict(list)
    for line in vcf_file:
        if line.startswith("#"):
            continue
        line_list = line.rstrip().split('\t')
        chrom = line_list[0]
        pos = int(line_list[1])-1#VCF file is 1-based
        if line_list[9].startswith('1|0') or line_list[9].startswith('0|1'):
            if window is None:
                vcf_dict[chrom].append((line_list[9].split(':')[0],
                                        pos,line_list[3].upper(),
                                        line_list[4].upper()))
            elif len(window.split(':')) == 1:
                if (chrom == window.split(':')[0]):
                    vcf_dict[chrom].append((line_list[9].split(':')[0],
                                            pos,line_list[3].upper(),
                                            line_list[4].upper()))
            else:
                if (chrom == window.split(':')[0]):
                    if len(window.split(':')[1].split('-')) == 1:
                        window_start= int(window.split(':')[1])
                        if pos > window_start - 1000000:
                            vcf_dict[chrom].append((line_list[9].split(':')[0],
                                                    pos,line_list[3].upper(),
                                                    line_list[4].upper()))
                    elif len(window.split(':')[1].split('-')) == 2:
                        window_start= int(window.split(':')[1].split('-')[0])
                        window_end= int(window.split(':')[1].split('-')[1])
                        if (pos > window_start - 1000000 and
                            pos < window_end + 1000000):
                            vcf_dict[chrom].append((line_list[9].split(':')[0],
                                                    pos,line_list[3].upper(),
                                                    line_list[4].upper()))
                    else:
                        raise Exception("Given window {} is not valid"
                                        ". Please give a valid window."
                                        "".format(window))
    vcf_file.close()
    return vcf_dict


def per_read_variant(vcf_dict,
                     bam_file,
                     chunk,
                     threads,
                     perReadinfo,
                     window):
    chrom_list = sorted(list(vcf_dict.keys()))
    for chrom in chrom_list:
        per_read_hp = defaultdict(lambda: defaultdict(list))
        bamiter, bam, count = openalignment(bam_file, chrom)
        if count > 0:
            feed_list= list(vcf_dict[chrom])
            feed_list = [feed_list[x:x+chunk]
                                 for x in range(0, len(feed_list),
                                                chunk)]
            feed_list = [feed_list[x:x+threads]
                                 for x in range(0, len(feed_list),
                                                threads)]
            description= "Tagging variants to reads from {}: ".format(chrom)
            with tqdm(total=len(feed_list),
                desc=description,
                bar_format="{l_bar}{bar} [ Estimated time left: {remaining} ]"
                                  ) as pbar:
                for vcf_info_list in feed_list:
                    p= mp.Pool(len(vcf_info_list))
                    results= p.starmap(get_base_info,
                                       list(zip(vcf_info_list,
                                                repeat(bam_file),
                                                repeat(chrom))))
                    p.close()
                    p.join()
                    for result in results:
                        if result is not None:
                            for read_info in result:
                                key,val= read_info
                                per_read_hp[key[0:-1]][key[-1]].append(val)
                    pbar.update(1)
            for key in list(per_read_hp.keys()):
                out_to_write= list()
                try:
                    hp1_variants= per_read_hp[key]['HP1']
                    hp1_count= str(len(hp1_variants))
                except:
                    hp1_count= '0'
                if hp1_count == '0':
                    hp1_variants= ['NA']
                    
                try:
                    hp2_variants= per_read_hp[key]['HP2']
                    hp2_count= str(len(hp2_variants))
                except:
                    hp2_count= '0'   
                if hp2_count == '0':
                    hp2_variants= ['NA']
                    
                if window is None or len(window.split(':')) == 1:
                    out_to_write= list(map(str,key)) + [','.join(sorted(hp1_variants)), 
                                                        ','.join(sorted(hp2_variants))]
                elif len(window.split(':')[1].split('-')) == 1:
                    window_start= int(window.split(':')[1])
                    if key[2] >= window_start:
                        out_to_write= list(map(str,key)) + [','.join(sorted(hp1_variants)), 
                                                            ','.join(sorted(hp2_variants))]
                else:
                    window_start= int(window.split(':')[1].split('-')[0])
                    window_end= int(window.split(':')[1].split('-')[1])
                    if key[2] >= window_start and key[1] <= window_end:
                        out_to_write= list(map(str,key)) + [','.join(sorted(hp1_variants)), 
                                                            ','.join(sorted(hp2_variants))]
                if out_to_write:
                    perReadinfo.write('\t'.join(out_to_write)+'\n')
        else:
            warnings.warn("{} in vcf file does not have any mapped reads in alignment "
                          "file Or alignment is truncated or corrupt indexed."
                          " Skipping it".format(chrom))
        per_read_hp.clear()
        bam.close()
    perReadinfo.close()
    

def main_methyl_call_processor(args):
    """
    This is the methyl_call_processor module which converts input
    methylation call to a bed format
    file for downstream steps.
    """
    MethylCallfile = os.path.abspath(args.MethylCallfile)
    tool,callthresh= args.tool_and_callthresh.split(":")
    tool= tool.lower()
    threads = args.threads
    chunk = args.chunk_size
    if args.motif.lower() == "cpg":
        args.motif = 'CG'
    else:
        raise TypeError("Please select motif type correctly (cpg)")
    if tool=="nanopolish":
        readID_index= 4
        start_index= 2
        strand_index= 1
        modprob_index= 5
    elif tool=="deepsignal":
        readID_index= 4
        start_index= 1
        strand_index= 2
        modprob_index= 7
    elif tool=="megalodon":
        readID_index= 0
        start_index= 3
        strand_index= 2
        modprob_index= 4
    else:
        raise TypeError("Please select tool correctly")
    meth = openfile(MethylCallfile)
    all_lines = 0
    for line in meth:
        all_lines += 1
    meth.close()
    feedlist = []
    chunklist = []
    readlist= []
    tqdm_add= 0
    meth = openfile(MethylCallfile)
    prev_info= next(meth).rstrip().split('\t')
    try: #check if first line is header
        prev_readID= prev_info[readID_index]
        prev_start= int(prev_info[start_index])
        prev_strand= prev_info[strand_index]
        prev_modprob= float(prev_info[modprob_index])
        readlist.append('\t'.join(prev_info))
        all_lines= all_lines - 1
    except: #skip first line
        prev_info= next(meth).rstrip().split('\t')
        prev_readID= prev_info[readID_index]
        prev_start= int(prev_info[start_index])
        prev_strand= prev_info[strand_index]
        readlist.append('\t'.join(prev_info))
        all_lines= all_lines - 2
    with tqdm(total=all_lines,
              desc="MethylCallProcessor: ", bar_format=
              "{l_bar}{bar} [ Estimated time left: {remaining} ]") as pbar:
        for line in meth:
            tqdm_add += 1
            line = line.rstrip()
            line_info= line.split('\t')
            start= int(line_info[start_index])
            if (line_info[readID_index] == prev_readID and
                line_info[strand_index] == prev_strand and
                abs(start -  prev_start) < 100000):
                prev_readID = line_info[readID_index]
                prev_strand = line_info[strand_index]
                prev_start= start
                readlist.append(line)
            else:
                chunklist.append(readlist)
                readlist = []
                readlist.append(line)
                prev_readID = line_info[readID_index]
                prev_strand = line_info[strand_index]
                prev_start= start
            if len(chunklist) == chunk:
                feedlist.append(chunklist)
                chunklist = []
            if len(feedlist) == threads:
                p = mp.Pool(threads)
                results = p.starmap(methcall2bed,
                                    list(zip(feedlist,
                                             repeat(args.motif),
                                             repeat(float(callthresh)),
                                             repeat(tool))))
                p.close()
                p.join()
                for result in results:
                    if result is not None:
                        for processed_line in result:
                            sys.stdout.write('\t'.join(processed_line)+'\n')
                feedlist = []
                pbar.update(tqdm_add)
                tqdm_add= 0
        else:
            chunklist.append(readlist)
            feedlist.append(chunklist)
            p = mp.Pool(len(feedlist))
            results = p.starmap(methcall2bed,
                                list(zip(feedlist,
                                         repeat(args.motif),
                                         repeat(float(callthresh)),
                                         repeat(tool))))
            p.close()
            p.join()
            for result in results:
                if result is not None:
                    for processed_line in result:
                        sys.stdout.write('\t'.join(processed_line)+'\n')
            feedlist = []
            pbar.update(tqdm_add)
    meth.close()
    

def main_phase(args):
    """
    This is the phase module which phase the nanopore reads and
    methylation data to corresponding haplotype using vcf file and
    processed methylation call file.
    """
    motif = args.motif
    if motif.lower() == "cpg":
        motif = 'CG'
    else:
        raise Exception("Courrently only CpG motif is supported")
    hapRatio = args.hapratio
    minvariant= args.min_variant
    vcf = os.path.abspath(args.vcf)
    if args.per_read is not None:
        per_read_file= os.path.abspath(args.per_read)
    bam_file = os.path.abspath(args.bam)
    threads = args.threads
    MinBaseQuality = args.min_base_quality
    AverageBaseQuality = args.average_base_quality
    chunk = args.chunk_size
    out = os.path.abspath(args.output)
    reference = args.reference
    out1 = out + '_NanoMethPhase_HP1'
    out2 = out + '_NanoMethPhase_HP2'
    MethylCallfile = args.methylcallfile
    MappingQuality = args.mapping_quality
    outformat = args.outformat.lower().split(',')
    fasta, tb= outformats_phase(outformat,reference,MethylCallfile)
    outexist_phase(args.overwrite,out1,out2)# Check if output files exist
    if args.per_read is None:
        perReadinfo= open(out1+"_HP2_PerReadInfo.tsv", 'w')
        perReadinfo.write("#Chromosome\tReadRefStart\tReadRefEnd\tReadID\t"
                          "Strand\tReadFlag:Is_Supplementary\tReadLength:MappingQuality\t"
                          "Position:BaseQuality:HP1_Variants\t"
                          "Position:BaseQuality:HP2_Variants\n")
        vcf_dict = vcf2dict_phase(vcf,
                                  args.window)
        per_read_variant(vcf_dict,
                            bam_file,
                            chunk,
                            threads,
                            perReadinfo,
                            args.window)
        per_read= openfile(out1+"_HP2_PerReadInfo.tsv")
    else:
        per_read= openfile(per_read_file)
    chrom_list= set()
    read_dictHP1 = defaultdict(list)
    read_dictHP2 = defaultdict(list)
    ignore_phased_indels= set()
    if not args.include_indels:
        ignore_phased_indels= get_indel(vcf)
    for line in per_read:
        if line.startswith("#"):
            continue
        line= line.rstrip().split('\t')
        if int(line[6].split(':')[1]) < MappingQuality:
            continue
        if not args.include_supplementary and line[5].split(':')[1]=="yes":
            continue
        key= (line[3],line[4],
                   line[5].split(':')[0],
                   line[6].split(':')[0])
        chrom_list.add(line[0])
        for haplotype in line[7].split(','):
            if haplotype != 'NA':
                position,phred,base= haplotype.split(':')
                if ((line[0],position) not in ignore_phased_indels and
                    int(phred) >= MinBaseQuality):
                    read_dictHP1[key].append(int(phred))
        for haplotype in line[8].split(','):
            if haplotype != 'NA':
                position,phred,base= haplotype.split(':')
                if ((line[0],position) not in ignore_phased_indels and
                    int(phred) >= MinBaseQuality):
                    read_dictHP2[key].append(int(phred))
    per_read.close()
    if read_dictHP1 or read_dictHP2:
        all_read = 0
        high_qual_reads = 0
        h1 = 0
        h2 = 0
        variant_tagged_reads= 0
        h1_bam2bis= 0
        h2_bam2bis= 0
        bamiter, bam_for_write, count = openalignment(bam_file,
                                                        args.window)
        (outHP1Sam, outHP2Sam, outHP12BisSam, outHP22BisSam,
                                   outCall1,outCall2,outFreq1,
                                   outFreq2)= outputs(outformat,
                                                      bam_for_write,out1,out2)
        bam_for_write.close()
        sys.stderr.write("Read Seperation Process Started\n")
        if args.window is not None:
            chrom_list= [args.window]
        else:
            chrom_list= sorted(chrom_list)
        for chrom in chrom_list:
            read_sam_list= list()
            bamiter, bam, count= openalignment(bam_file, chrom)
            description= "Processing reads from {}: ".format(chrom)
            with tqdm(total=count,
                      desc= description) as pbar:
                for read in bamiter:
                    pbar.update(1)
                    methylated_sites = list()
                    unmethylated_sites = list()
                    llr_methylated = list()
                    llr_unmethylated = list()
                    all_read += 1
                    mp_quality = read.mapping_quality
                    if not args.include_supplementary:
                        if (read.is_unmapped or mp_quality < MappingQuality or
                                read.is_secondary or read.is_supplementary or
                                read.is_qcfail or read.is_duplicate):
                            continue
                    else:
                        if (read.is_unmapped or mp_quality < MappingQuality or
                                read.is_secondary or
                                read.is_qcfail or read.is_duplicate):
                            continue
                    high_qual_reads += 1
                    (true_ref_name, strand, flag, read_id ,
                     read_seq, read_len, cigar, rnext, pnext, tlen ,
                     base_qualities, start, end, ref_seq, ref_len,
                                      all_tags) = bam_info_extractor(read,
                                                                     reference,
                                                                     fasta)
                    key = (read_id,strand,str(flag),str(read_len))
                    if (key not in read_dictHP1 and
                        key not in read_dictHP2):
                        continue
                    if key in read_dictHP1 or key in read_dictHP2:
                        variant_tagged_reads += 1
                    hp1_quals = read_dictHP1[key]
                    hp1_count = len(hp1_quals)
                    if hp1_count > 1:
                        mean_hp1_quals= statistics.mean(hp1_quals)
                    elif hp1_count == 1:
                        mean_hp1_quals= hp1_quals[0]
                    else:
                        mean_hp1_quals= 0
                    hp2_quals= read_dictHP2[key]
                    hp2_count = len(hp2_quals)
                    if hp2_count > 1:
                        mean_hp2_quals= statistics.mean(hp2_quals)
                    elif hp2_count == 1:
                        mean_hp2_quals= hp2_quals[0]
                    else:
                        mean_hp2_quals= 0
                    hp1_mean_count= len([x for x in hp1_quals if
                                         x >= AverageBaseQuality])
                    hp2_mean_count= len([x for x in hp2_quals if
                                         x >= AverageBaseQuality])
                    if (hp1_count > hp2_count
                            and hp2_count/hp1_count <= hapRatio
                            and hp1_count >= minvariant
                            and (mean_hp1_quals >= mean_hp2_quals
                                  or hp1_mean_count >= hp2_mean_count)):
                        h1 += 1
                        if 'bam' in outformat:
                            outHP1Sam.write(read)
                        if 'methylcall' in outformat or 'bam2bis' in outformat:
                            try:
                                records = tb.query(true_ref_name, start, end)
                            except:
                                warnings.warn("{}:{}-{} does not exist in the "
                                              "MethylCallFile."
                                              "Skipping it".format(read_id,
                                                                       start,
                                                                       end))
                                continue
                            for record in records:
                                if read_id == record[4] and strand == record[3]:
                                    if record[7] != 'NA':
                                        methylated_sites += map(int,
                                                            record[7].split(','))
                                        llr_methylated += map(float,
                                                            record[5].split(','))
                                    if record[8] != 'NA':
                                        unmethylated_sites += map(int,
                                                            record[8].split(','))
                                        llr_unmethylated += map(float,
                                                            record[6].split(','))
                            methylcall_dict= dict()
                            for i,j in zip(methylated_sites +
                                           unmethylated_sites,
                                           llr_methylated +
                                           llr_unmethylated):
                                if (i >= start and i <= end ):
                                    if i not in methylcall_dict:
                                        methylcall_dict[i]= [record[0],
                                                            i,
                                                            i+1,
                                                            strand,
                                                            read_id, j]
                                    elif abs(j) > abs(methylcall_dict[i][-1]):
                                        methylcall_dict[i]= [record[0],
                                                            i,
                                                            i+1,
                                                            strand,
                                                            read_id, j]
                            if 'methylcall' in outformat:
                                for key,val in methylcall_dict.items():
                                    outCall1.write('\t'.join(map(str,val))+'\n')
                            if 'bam2bis' in outformat:
                                h1_bam2bis += 1
                                read_sam_list.append(['HP1', strand, read_id,
                                    flag, true_ref_name, start,
                                    mp_quality, ref_len, rnext ,
                                    pnext, tlen, ref_seq, '*',
                                    all_tags,
                                    [i - start for i in methylcall_dict.keys()
                                        if methylcall_dict[i][-1] > 0],
                                    [i - start for i in methylcall_dict.keys()
                                        if methylcall_dict[i][-1] <= 0]])
                    elif (hp2_count > hp1_count
                          and hp1_count/hp2_count <= hapRatio
                          and hp2_count >= minvariant
                          and (mean_hp2_quals >= mean_hp1_quals
                                  or hp2_mean_count >= hp1_mean_count)):
                        h2 += 1
                        if 'bam' in outformat:
                            outHP2Sam.write(read)
                        if 'methylcall' in outformat or 'bam2bis' in outformat:
                            try:
                                records = tb.query(true_ref_name, start, end)
                            except:
                                warnings.warn("{}:{}-{} does not exist in the "
                                              "MethylCallFile. This read might not present any detected methylation."
                                              "Skipping it".format(read_id,
                                                                   start,
                                                                   end))
                                continue
                            for record in records:
                                if read_id == record[4] and strand == record[3]:
                                    if record[7] != 'NA':
                                        methylated_sites += map(int,
                                                            record[7].split(','))
                                        llr_methylated += map(float,
                                                            record[5].split(','))
                                    if record[8] != 'NA':
                                        unmethylated_sites += map(int,
                                                            record[8].split(','))
                                        llr_unmethylated += map(float,
                                                            record[6].split(','))
                            methylcall_dict= dict()
                            for i,j in zip(methylated_sites +
                                           unmethylated_sites,
                                           llr_methylated +
                                           llr_unmethylated):
                                if (i >= start and i <= end ):
                                    if i not in methylcall_dict:
                                        methylcall_dict[i]= [record[0],
                                                            i,
                                                            i+1,
                                                            strand,
                                                            read_id, j]
                                    elif abs(j) > abs(methylcall_dict[i][-1]):
                                        methylcall_dict[i]= [record[0],
                                                            i,
                                                            i+1,
                                                            strand,
                                                            read_id, j]
                            if 'methylcall' in outformat:
                                for key,val in methylcall_dict.items():
                                    outCall2.write('\t'.join(map(str,val))+'\n')
                            if 'bam2bis' in outformat:
                                h2_bam2bis += 1
                                read_sam_list.append(['HP2', strand, read_id,
                                    flag, true_ref_name, start,
                                    mp_quality, ref_len, rnext ,
                                    pnext, tlen, ref_seq, '*',
                                    all_tags,
                                    [i - start for i in methylcall_dict.keys()
                                        if methylcall_dict[i][-1] > 0],
                                    [i - start for i in methylcall_dict.keys()
                                        if methylcall_dict[i][-1] <= 0]])
                    if len(read_sam_list) == (threads * chunk):
                        p = mp.Pool(threads)
                        results = p.map(read2bis, read_sam_list)
                        p.close()
                        p.join()
                        for result in results:
                            if result is not None:
                                if result[0] == "HP1":
                                    alignmentwriter(result, outHP12BisSam)
                                else:
                                    alignmentwriter(result, outHP22BisSam)
                        read_sam_list= list()
                else:
                    if read_sam_list:
                        p = mp.Pool(threads)
                        results = p.map(read2bis, read_sam_list)
                        p.close()
                        p.join()
                        for result in results:
                            if result is not None:
                                if result[0] == "HP1":
                                    alignmentwriter(result, outHP12BisSam)
                                else:
                                    alignmentwriter(result, outHP22BisSam)
            bam.close()
        if 'bam' in outformat:
            outHP1Sam.close()
            outHP2Sam.close()
            sys.stderr.write("Phased Bam files are ready\n")
        if 'methylcall' in outformat:
            outCall1.close()
            outCall2.close()
            dict_mod, dict_all = FrequencyCalculator(out1+"_MethylCall.tsv")
            for key, val in dict_all.items():
                modCall = dict_mod[key]
                freq = modCall/val
                outwrite = '\t'.join(list(key) + [str(val),
                                     str(modCall), str(freq)+'\n'])
                outFreq1.write(outwrite)
            outFreq1.close()
            dict_mod, dict_all = FrequencyCalculator(out2+"_MethylCall.tsv")
            for key, val in dict_all.items():
                modCall = dict_mod[key]
                freq = modCall/val
                outwrite = '\t'.join(list(key) + [str(val),
                                     str(modCall), str(freq)+'\n'])
                outFreq2.write(outwrite)
            outFreq2.close()
            sys.stderr.write("Phased Methylation Call and "
                                     "frequency files are ready\n")
        if 'bam2bis' in outformat:
            outHP12BisSam.close()
            outHP22BisSam.close()
            sys.stderr.write("bam2bis output files are ready\n")
        if not args.include_supplementary:
            sys.stderr.write("Job Finished.\n"
                                         "Number of all reads at processed chroms: {}\n"
                                         "Number of nonsuplementary/"
                                         "nonsecondary/notPCRdup/qcPassed"
                                         " mapped reads with quality more "
                                         "than {} at processed chroms: {}.\n"
                                         "Number of reads with at least one"
                                         " tagged phased variant: {}\n"
                                         "Number of HP1 reads: {}\n"
                                         "Number of HP2 reads: {}\n"
                                         "".format(all_read,
                                                   MappingQuality,
                                                   high_qual_reads,
                                                   variant_tagged_reads,
                                                   h1,
                                                   h2))
        else:
            sys.stderr.write("Job Finished.\n"
                                         "Number of all reads at processed chroms: {}\n"
                                         "Number of "
                                         "nonsecondary/notPCRdup/qcPassed"
                                         " mapped reads with quality more "
                                         "than {} at processed chroms: {}.\n"
                                         "Number of reads with at least one"
                                         " tagged phased variant: {}\n"
                                         "Number of HP1 reads: {}\n"
                                         "Number of HP2 reads: {}\n"
                                         "".format(all_read,
                                                   MappingQuality,
                                                   high_qual_reads,
                                                   variant_tagged_reads,
                                                   h1,
                                                   h2))
    else:
        sys.stderr.write("There is no phased variant in your vcf file or "
                                     "Noe reads could be tagged.\n")
    if reference is not None:
        fasta.close()

def main_bam2bis(args):
    motif = args.motif
    if motif.lower() == "cpg":
        motif = 'CG'
    else:
        raise Exception("Courrently only CpG motif is supported")
    bam_file = os.path.abspath(args.bam)
    threads = args.threads
    chunk = args.chunk_size
    out = os.path.abspath(args.output)
    reference = os.path.abspath(args.reference)
    MethylCallfile = os.path.abspath(args.methylcallfile)
    MappingQuality = args.mapping_quality
    fasta = pysam.FastaFile(reference)
    if not os.path.isfile(MethylCallfile+".tbi"):
        raise Exception("Could not find index file for methylation call file.")
    tb = tabix.open(MethylCallfile)
    if args.methylation:
        if not args.overwrite and (os.path.isfile(out+"MethylationCall.tsv") or
                                              os.path.isfile(out+
                                              "MethylationFrequency.tsv" )):
            raise FileExistsError("The selected output methylation call and"
                                              " frequency files {} already "
                                              "exist. Select --overwrite option"
                                              " if you want to overwrite them"
                                              "".format(out))
        outCall= open(out+"MethylationCall.tsv", 'w')
        outCall.write("\t".join(["chromosome",
                                "start",
                                "end",
                                "strand",
                                "read_name",
                                "llr_Or_DeltaProb\n"]))
        outFreq= open(out+"MethylationFrequency.tsv", 'w')
        outFreq.write("\t".join(["chromosome",
                                "start",
                                "end",
                                "strand",
                                "NumOfAllCalls",
                                "NumOfModCalls",
                                "MethylFreq\n"]))
    if not args.overwrite and os.path.isfile(out+"_Converted2Bisulfite.bam"):
        raise FileExistsError("The selected output bam file {} already "
                                              "exists. Select --overwrite "
                                              "option if you want to overwrite"
                                              " it".format(out))
    high_quality_reads= 0
    converted_reads= 0
    all_read= 0
    bamiter, bam, counts= openalignment(bam_file, args.window)
    outBisbam= pysam.AlignmentFile(out+"_Converted2Bisulfite.bam",
                                            "wb", template=bam)
    bam.close()
    if args.window is None:
        chroms= sorted(getChromsFromBAM(bam_file))
    else:
        chroms= [args.window]
    for chrom in chroms:
        read_sam_list= list()
        bamiter, bam, counts= openalignment(bam_file, chrom)
        description= "Converting reads from {}: ".format(chrom)
        if counts > 0:
            with tqdm(total=counts,
                      desc= description) as pbar:
                for read in bamiter:
                    mp_quality = read.mapping_quality
                    all_read += 1
                    methylated_sites = list()
                    unmethylated_sites = list()
                    llr_methylated = list()
                    llr_unmethylated = list()
                    pbar.update(1)
                    if not args.include_supplementary:
                        if (read.is_unmapped or mp_quality < MappingQuality or
                                read.is_secondary or read.is_supplementary or
                                read.is_qcfail or read.is_duplicate):
                            continue
                    else:
                        if (read.is_unmapped or mp_quality < MappingQuality or
                                read.is_secondary or
                                read.is_qcfail or read.is_duplicate):
                            continue
                    high_quality_reads += 1
                    (true_ref_name, strand, flag, read_id ,
                     read_seq, read_len, cigar, rnext, pnext, tlen ,
                     base_qualities, start, end, ref_seq, ref_len,
                     all_tags) = bam_info_extractor(read,
                                                    reference,
                                                    fasta)
                    try:
                        records = tb.query(true_ref_name, start, end)
                    except:
#                                warnings.warn("{}:{}-{} does not exist in the "
#                                              "MethylCallFile.This read might not present any detected methylation."
#                                              "Skipping it".format(read_id,
#                                                                       start,
#                                                                       end))
                        continue
                    for record in records:
                        if read_id == record[4] and strand == record[3]:
                            if record[7] != 'NA':
                                methylated_sites += map(int,
                                                    record[7].split(','))
                                llr_methylated += map(float,
                                                    record[5].split(','))
                            if record[8] != 'NA':
                                unmethylated_sites += map(int,
                                                    record[8].split(','))
                                llr_unmethylated += map(float,
                                                    record[6].split(','))
                    methylcall_dict= dict()
                    for i,j in zip(methylated_sites +
                                   unmethylated_sites,
                                   llr_methylated +
                                   llr_unmethylated):
                        if (i >= start and i <= end ):
                            if i not in methylcall_dict:
                                methylcall_dict[i]= [record[0],
                                                    i,
                                                    i+1,
                                                    strand,
                                                    read_id, j]
                            elif abs(j) > abs(methylcall_dict[i][-1]):
                                methylcall_dict[i]= [record[0],
                                                    i,
                                                    i+1,
                                                    strand,
                                                    read_id, j]
                    if args.methylation:
                        for key,val in methylcall_dict.items():
                            outCall.write('\t'.join(map(str,val))+'\n')
                    read_sam_list.append(['NON', strand, read_id,
                        flag, true_ref_name, start,
                        mp_quality, ref_len, rnext ,
                        pnext, tlen, ref_seq, '*',
                        all_tags,
                        [i - start for i in methylcall_dict.keys()
                            if methylcall_dict[i][-1] > 0],
                        [i - start for i in methylcall_dict.keys()
                            if methylcall_dict[i][-1] <= 0]])
                    if len(read_sam_list) == (threads * chunk):
                        converted_reads += len(read_sam_list)
                        p = mp.Pool(threads)
                        results = p.map(read2bis, read_sam_list)
                        p.close()
                        p.join()
                        for read in results:
                            if read is not None:
                                alignmentwriter(read, outBisbam)
                        read_sam_list= list()
                else:
                    if read_sam_list:
                        converted_reads += len(read_sam_list)
                        p = mp.Pool(threads)
                        results = p.map(read2bis,read_sam_list)
                        p.close()
                        p.join()
                        for read in results:
                            if read is not None:
                                alignmentwriter(read, outBisbam)
        bam.close()
    outBisbam.close()
    if args.methylation:
        outCall.close()
        dict_mod, dict_all = FrequencyCalculator(out+"MethylationCall.tsv")
        for key, val in dict_all.items():
            modCall = dict_mod[key]
            freq = modCall/val
            outwrite = '\t'.join(list(key)
                                 + [str(val), str(modCall), str(freq)+'\n'])
            outFreq.write(outwrite)
        outFreq.close()
    if not args.include_supplementary:
        sys.stderr.write("Job Finished.\n"
                         "Number of all reads: {}\n"
                         "Number of nonsuplementary/"
                         "nonsecondary/notPCRdup/qcPassed"
                         " mapped reads with quality more "
                         "than {}: {}.\n"
                         "Number of converted reads: {}\n"
                         "".format(all_read, MappingQuality,
                                   high_quality_reads, converted_reads))
    else:
        sys.stderr.write("Job Finished.\n"
                         "Number of all reads: {}\n"
                         "Number of "
                         "nonsecondary/notPCRdup/qcPassed"
                         " mapped reads with quality more "
                         "than {}: {}.\n"
                         "Number of converted reads: {}\n"
                         "".format(all_read, MappingQuality,
                                   high_quality_reads, converted_reads))
    fasta.close()

def main_dma(args):
    """
    This is the DMA module which does differential methylation analysis
    using DSS R package to find differentially methylated regions.
    """
    if os.path.isdir(os.path.abspath(args.case)):
        cases= []
        for (dirpath, dirnames, filenames) in os.walk(os.path.abspath(args.case)):
            for filename in filenames:
                cases.append(dirpath+'/'+filename)
    else:    
        cases= [os.path.abspath(args.case)]
    
    if os.path.isdir(os.path.abspath(args.control)):
        controls= []
        for (dirpath, dirnames, filenames) in os.walk(os.path.abspath(args.control)):
            for filename in filenames:
                controls.append(dirpath+'/'+filename)
    else:    
        controls= [os.path.abspath(args.control)]
        
    out_dir = os.path.abspath(args.out_dir)
    out_prefix = out_dir+'/'+(args.out_prefix)
    coverage = args.coverage
    columns = args.columns
    Rscript = args.Rscript #  os.path.abspath(args.Rscript)
    script = os.path.abspath(args.script_file)
    dis_merge = args.dis_merge
    minlen = args.minlen
    minCG = args.minCG
    smoothing_span = args.smoothing_span
    smoothing_flag = args.smoothing_flag.upper()
    equal_disp = args.equal_disp.upper()
    pval_cutoff = args.pval_cutoff
    delta_cutoff = args.delta_cutoff
    pct_sig = args.pct_sig
    # check if outputs exist
    if args.columns:
        check_outs = [x for x in glob.glob("{}*{}".format(out_prefix+"_ReadyForDSS_case",".tsv"))]
        if check_outs and not args.overwrite:
            raise FileExistsError("The selected output files {} already "
                                  "exist. Select --overwrite option if you "
                                  "want to overwrite them or use a different "
                                  "prefix".format(check_outs))
    else:
        check_outs = [x for x in glob.glob("{}*DM*.txt".format(out_prefix))]
        if check_outs and not args.overwrite:
            raise FileExistsError("The selected output files {} already "
                                  "exist. Select --overwrite option if you "
                                  "want to overwrite them or use a different "
                                  "prefix".format(check_outs))
    if args.columns:
        ready_cases = []
        ready_controls = []
        columns = args.columns
        columns = list(map(int, columns.split(',')))
        if len(columns) == 4:
            col1, col2, col3, col4 = columns
            out_putNumber = 0
            for case in cases:
                out_putNumber += 1
                case_out = open("{}_ReadyForDSS_case{}.tsv".format(
                    out_prefix, out_putNumber), "w")
                with openfile(case) as case_file:
                    next(case_file)  # Exclude header
                    case_out.write("Chromosome\tPosition\tAllCalls\tModCalls\n")
                    for line in case_file:
                        line = line.rstrip().split('\t')
                        cov = int(line[col3-1])
                        if cov >= coverage:
                            chrom = line[col1-1]
                            start = line[col2-1]
                            mod_sites = round(cov * float(line[col4-1]))
                            case_out.write('\t'.join([chrom,
                                                      start,
                                                      str(cov),
                                                      str(mod_sites)+'\n']))
                case_out.close()
                ready_cases.append("{}_ReadyForDSS_case{}.tsv".format(
                    out_prefix, out_putNumber))
            out_putNumber = 0
            for control in controls:
                out_putNumber += 1
                control_out = open("{}_ReadyForDSS_control{}.tsv".format(
                    out_prefix, out_putNumber), "w")
                with openfile(control) as control_file:
                    next(control_file)
                    control_out.write("Chromosome\tPosition\tAllCalls\tModCalls\n")
                    for line in control_file:
                        line = line.rstrip().split('\t')
                        cov = int(line[col3-1])
                        if cov >= coverage:
                            chrom = line[col1-1]
                            start = line[col2-1]
                            mod_sites = round(cov * float(line[col4-1]))
                            control_out.write('\t'.join([chrom,
                                                         start,
                                                         str(cov),
                                                         str(mod_sites)+'\n']))
                control_out.close()
                ready_controls.append("{}_ReadyForDSS_control{}.tsv".format(
                    out_prefix, out_putNumber))
        elif len(columns) == 5:
            col1, col2, col3, col4, col5 = columns
            out_putNumber = 0
            for case in cases:
                out_putNumber += 1
                case_out = open("{}_ReadyForDSS_case{}.tsv".format(
                    out_prefix, out_putNumber), "w")
                cov_dict = defaultdict(int)
                mod_sites_dict = defaultdict(int)
                with openfile(case) as case_file:
                    next(case_file)
                    case_out.write("Chromosome\tPosition\tAllCalls\tModCalls\n")
                    for line in case_file:
                        line = line.rstrip().split('\t')
                        cov = int(line[col4-1])
                        if cov >= coverage:
                            strand = line[col3-1]
                            start = int(line[col2-1])
                            chrom = line[col1-1]
                            if strand == "-":
                                start = start - 1
                            key = (chrom,start)
                            cov_dict[key] += cov
                            mod_sites_dict[key]+= round(cov*float(line[col5-1]))
                    for key in sorted(mod_sites_dict.keys()):
                        cov = cov_dict[key]
                        chrom = key[0]
                        start = key[1]
                        mod_sites = mod_sites_dict[key]
                        case_out.write('\t'.join([chrom,
                                                  str(start),
                                                  str(cov),
                                                  str(mod_sites)+'\n']))
                    case_out.close()
                ready_cases.append("{}_ReadyForDSS_case{}.tsv".format(
                    out_prefix, out_putNumber))
            out_putNumber = 0
            for control in controls:
                out_putNumber += 1
                control_out = open("{}_ReadyForDSS_control{}.tsv".format(
                    out_prefix, out_putNumber), "w")
                cov_dict = defaultdict(int)
                mod_sites_dict = defaultdict(int)
                with openfile(control) as control_file:
                    next(control_file)
                    control_out.write("Chromosome\tPosition\tAllCalls\tModCalls\n")
                    for line in control_file:
                        line = line.rstrip().split('\t')
                        cov = int(line[col4-1])
                        if cov >= coverage:
                            strand = line[col3-1]
                            start = int(line[col2-1])
                            chrom = line[col1-1]
                            if strand == "-":
                                start = start - 1
                            key = (chrom,start)
                            cov_dict[key] += cov
                            mod_sites_dict[key] += round(cov*float(line[col5-1]))
                    for key in sorted(mod_sites_dict.keys()):
                        cov = cov_dict[key]
                        chrom = key[0]
                        start = key[1]
                        mod_sites = mod_sites_dict[key]
                        control_out.write('\t'.join([chrom,
                                                     str(start),
                                                     str(cov),
                                                     str(mod_sites)+'\n']))
                    control_out.close()
                ready_controls.append("{}_ReadyForDSS_control{}.tsv".format(
                    out_prefix,
                    out_putNumber))
        else:
            raise TypeError("Please enter columns correctly")
        ready_cases = ','.join(ready_cases)
        ready_controls = ','.join(ready_controls)
    else:
        ready_cases = ','.join(cases)
        ready_controls = ','.join(controls)
    subprocess.call(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(Rscript,
                                                           script,
                                                           ready_cases,
                                                           ready_controls,
                                                           out_prefix,
                                                           dis_merge,
                                                           minlen,
                                                           minCG,
                                                           smoothing_span,
                                                           smoothing_flag,
                                                           pval_cutoff,
                                                           delta_cutoff,
                                                           pct_sig,
                                                           equal_disp),
        shell=True)


def phase_parser(subparsers):
    """
    Specific argument parser for phase command.
    """
    sub_phase = subparsers.add_parser("phase",
                                      add_help=False,
                                      help="Phasing reads and Methylation.",
                                      description="Phasing reads and "
                                      "Methylation")
    sp_input = sub_phase.add_argument_group("required arguments")
    sp_input.add_argument("--bam", "-b",
                          action="store",
                          type=str,
                          required=True,
                          help="The path to the cordinate sorted bam file.")
    sp_input.add_argument("--vcf", "-v",
                          action="store",
                          type=str,
                          required=True,
                          help="The path to the phased vcf file.")
    sp_input.add_argument("--output", "-o",
                          action="store",
                          type=str,
                          required=True,
                          help=("The path to directory and prefix to save "
                                "files. e.g path/to/directory/prefix"))
    sp_input = sub_phase.add_argument_group("conditional required arguments based"
                                            " on selected output format(s)")
    sp_input.add_argument("--reference", "-r",
                          action="store",
                          type=str,
                          required=False,
                          default=None,
                          help=("The path to the reference file in case you "
                                "selected bam2bis output format. Fasta file "
                                "must be already indexed using samtools "
                                "faidx."))
    sp_input.add_argument("--methylcallfile", "-mc",
                          action="store",
                          type=str,
                          required=False,
                          default=None,
                          help=("If you have selected methylcall or bam2bis "
                                "output format to phase methylation or make mock bisulfite bams, "
                                "give the path to the bgziped and indexed methylation "
                                "call file from methyl_call_processor Module."))
    sp_input = sub_phase.add_argument_group("General optional arguments")
    sp_input.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    sp_input.add_argument("--outformat", "-of",
                          action="store",
                          type=str,
                          required=False,
                          default="bam2bis,methylcall",
                          help=("What type of output you want (bam,bam2bis,"
                                "methylcall). Default is bam2bis,methylcall."
                                "bam: outputs phased reads to seperate bam "
                                "files."
                                "bam2bis: outputs phased reads to seperate "
                                "bam files converted to bisulfite bam format "
                                "for visualisation in IGV."
                                "methylcall: outputs phased methylcall and "
                                "methylation frequency files for seperate "
                                "haplotypes. You can select any format and "
                                "multiple or all of them seperated by comma."))
    sp_input.add_argument("--per_read", "-pr",
                          action="store",
                          type=str,
                          required=False,
                          default=None,
                          help="If it is your second try and you have per read info"
                               " file from the first try you can specify the per-read"
                               " file to make the process faster. This also enables you "
                               "to try different threshols of options (-mv, -mbq, mq, -hr, -abq),"
                               " include/exclude indels, include/exclude supp reads.")
    sp_input.add_argument("--min_variant", "-mv",
                          action="store",
                          type=int,
                          required=False,
                          default=1,
                          help=("minimum number of phased variants must a read "
                                "have to be phased. Default is 1"))
    sp_input.add_argument("--hapratio", "-hr",
                          action="store",
                          type=float,
                          required=False,
                          default=0.75,
                          help=("0-1 .The maximum ratio (# of SNVs from one halotype over the other)"
                                " between haplotypes to tag as H1 or H2 (H2/H1 when #SNVs at H1 > H2"
                                " and H1/H2 when #SNVs at H2 > H1). Default is 0.75"))
    sp_input.add_argument("--mapping_quality", "-mq",
                          action="store",
                          type=int,
                          required=False,
                          default=20,
                          help=("An integer value to specify minimum"
                                " read mapping quality. "
                                "Default is 20"))
    sp_input.add_argument("--min_base_quality", "-mbq",
                          action="store",
                          type=int,
                          required=False,
                          default=7,
                          help=("Only include bases with this minimum base quality"
                                ". Default is 7."))
    sp_input.add_argument("--average_base_quality", "-abq",
                          action="store",
                          type=int,
                          required=False,
                          default=20,
                          help=("Minimum quality that variants tagged to a haplotype "
                                "should have compare to the other haplotype when average "
                                "of qualities is not informative. "
                                "This will be used "
                                "when the average base quality of variants mapped"
                                " to two haplotypes for one read is not informative and"
                                " decision cannot be made "
                                "(e.g. when 10 variants of HP1 mapped to a read "
                                "with average quality of 30, but only one variant "
                                "from HP2 mapped to the same read with bq=35). "
                                "Then, instead of considering average of qualities, "
                                "the tool will count number of variants in both haplotypes "
                                " that meet the given average_base_quality and uses"
                                " the counts to make decision."
                                " Default is 20."))
    sp_input.add_argument("--include_indels", "-ind",
                          action="store_true",
                          required=False,
                          help="Also include indels for read phasing to haplotypes.")
    sp_input.add_argument("--include_supplementary", "-is",
                          action="store_true",
                          required=False,
                          help="Also include supplementary reads")
    sp_input.add_argument("--motif", "-mt",
                          action="store",
                          type=str,
                          required=False,
                          default="cpg",
                          help=("The motif you called methylation for (cpg), "
                                "Currently just cpg."))
    sp_input.add_argument("--window", "-w",
                          action="store",
                          type=str,
                          required=False,
                          help=("if you want to only phase read for a region "
                                "or chromosome. You must insert region like "
                                "this chr1 or chr1:1000-100000."))
    sp_input.add_argument("--threads", "-t",
                          action="store",
                          type=int,
                          required=False,
                          default=4,
                          help="Number of parallel processes. Default is 4")
    sp_input.add_argument("--chunk_size", "-cs",
                          action="store",
                          type=int,
                          required=False,
                          default=100,
                          help=("Number of reads send to each proccessor. "
                                "Default is 100"))
    sp_input.add_argument("--overwrite", "-ow",
                          action="store_true",
                          required=False,
                          help="If output files exist overwrite them")
    sub_phase.set_defaults(func=main_phase)


def methyl_call_processor_parser(subparsers):
    """
    Specific argument parser for methyl calling command.
    """
    sub_methyl_call_processor = subparsers.add_parser(
        "methyl_call_processor",
        add_help=False,
        help=("Preparing methylation call file for methylation phasing"
              " or mock bisulfite conversion of bam file."),
        description=("Preparing methylation call file for methylation phasing. "
                     "Extended usage: nanomethphase methyl_call_processor -mc "
                     "[FILE] | sort -k1,1 -k2,2n -k3,3n | bgzip > "
                     "[FILE].bed.gz && tabix -p bed [FILE].bed.gz"))
    smp_input = sub_methyl_call_processor.add_argument_group("required arguments")
    smp_input.add_argument("--MethylCallfile", "-mc",
                           action="store",
                           type=str,
                           required=True,
                           default=None,
                           help=("The path to the per-read methylation call "
                                 "file."))
    smp_input = sub_methyl_call_processor.add_argument_group("optional arguments")
    smp_input.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    smp_input.add_argument("--tool_and_callthresh", "-tc",
                           action="store",
                           type=str,
                           required=False,
                           default="nanopolish:2",
                           help=("Software you have used for methylation calling "
                                 "(nanoplish, megalodon, deepsignal):"
                                 "methylation call threshold for considering a site as "
                                 "methylated, unmethylated or ambiguous in methylation call file. "
                                 "Default is nanopolish:2 which is when methylation"
                                 " calling performed by nanopolish and a CpG with llr >= 2 will be considered "
                                 "as methylated and llr <= -2 as unmethylated, anything "
                                 "in between will be considered as ambiguous call."
                                 "For megalodon, call thresold will be delta probability "
                                 "(prob_methylated - prob_unmethylated) "
                                 "which is e^mod_log_prob - (1 - e^mod_log_prob). "
                                 "For example, with a call threshold of 0.6 (0.8-0.2) CpGs "
                                 "between 0.8 and 0.2 probability"
                                 " will be considered as ambiguous and >=0.8 as methylated "
                                 "and <=0.2 as unmethylated.\n"
                                 "For deepsignal, as for megalodon, this call threshold is delta probability"
                                 " (prob_methylated - prob_unmethylated)."
                                 "NOTE: Megalodon per-read text file must be for only 5mC CpGs "
                                 ". Do not use per-read text file where there are predictions"
                                 " for 2 or more modifications (e.g. 5mC and 5hmC)"))
    smp_input.add_argument("--motif", "-mf",
                           action="store",
                           type=str,
                           required=False,
                           default="cpg",
                           help=("The motif you called methylation for (cpg), "
                                 "Currently just cpg."))
    smp_input.add_argument("--threads", "-t",
                           action="store",
                           type=int,
                           required=False,
                           default=4,
                           help="Number of parallel processes. Default is 4")
    smp_input.add_argument("--chunk_size", "-cs",
                           action="store",
                           type=int,
                           required=False,
                           default=100,
                           help=("Number of reads send to each proccessor. "
                                 "Default is 100"))
    sub_methyl_call_processor.set_defaults(func=main_methyl_call_processor)


def bam2bis_parser(subparsers):
    """
    Specific argument parser for phase command.
    """
    sub_bam2bis = subparsers.add_parser(
        "bam2bis",
        add_help=False,
        help=("Convert a bam file to a bisulfite format for visualization "
              "in IGV"),
        description=("Convert a bam file to a bisulfite format for nice "
                    "visualization in IGV"))
    sbb_input = sub_bam2bis.add_argument_group("required arguments")
    sbb_input.add_argument("--bam", "-b",
                          action="store",
                          type=str,
                          required=True,
                          help="The path to the cordinate sorted bam file.")
    sbb_input.add_argument("--reference", "-r",
                          action="store",
                          type=str,
                          required=True,
                          help=("The path to the reference file. Fasta file "
                                "must be already indexed using samtools "
                                "faidx."))
    sbb_input.add_argument("--methylcallfile", "-mc",
                          action="store",
                          type=str,
                          required=True,
                          help=("The path to the the bgziped and indexed "
                                "methylation call file from "
                                "methyl_call_processor Module."))
    sbb_input.add_argument("--output", "-o",
                          action="store",
                          type=str,
                          required=True,
                          help=("The path to the output directory and desired"
                                " prefix."))
    sbb_input = sub_bam2bis.add_argument_group("optional arguments")
    sbb_input.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    sbb_input.add_argument("--mapping_quality", "-mq",
                          action="store",
                          type=int,
                          required=False,
                          default=20,
                          help=("An integer value to specify minimum"
                                " mapping quality of the read. "
                                "Default is 20"))
    sbb_input.add_argument("--motif", "-mt",
                          action="store",
                          type=str,
                          required=False,
                          default="cpg",
                          help=("The motif you called methylation for (cpg), "
                                "Currently just cpg."))
    sbb_input.add_argument("--methylation", "-met",
                          action="store_true",
                          required=False,
                          help="Output methylation call and frequency for "
                          "converted reads.")
    sbb_input.add_argument("--include_supplementary", "-is",
                          action="store_true",
                          required=False,
                          help="Also include supplementary reads")
    sbb_input.add_argument("--window", "-w",
                          action="store",
                          type=str,
                          required=False,
                          help=("if you want to only convert reads for a "
                                "region or chromosome. You must insert region "
                                "like this chr1 or chr1:1000-100000."))
    sbb_input.add_argument("--threads", "-t",
                          action="store",
                          type=int,
                          required=False,
                          default=4,
                          help="Number of parallel processes. Default is 4")
    sbb_input.add_argument("--chunk_size", "-cs",
                          action="store",
                          type=int,
                          required=False,
                          default=100,
                          help=("Number of reads send to each proccessor. "
                                "Default is 100"))
    sbb_input.add_argument("--overwrite", "-ow",
                          action="store_true",
                          required=False,
                          help="If output files exist overwrite it")
    sub_bam2bis.set_defaults(func=main_bam2bis)


def dma_parser(subparsers):
    """
    Specific argument parser for analysis of differential methylation
    command.
    """
    sub_dma = subparsers.add_parser(
        "dma",
        add_help=False,
        help=("Differential Methylation analysis for two group only (to find  "
              "DMRs using phased frequency results) using DSS R package."),
        description=("Differential Methylation analysis for two group only "
                     "(to find DMRs using phased frequency results) using DSS "
                     "R package.\n"))
    sdma_input = sub_dma.add_argument_group("required arguments")
    sdma_input.add_argument("--case", "-ca",
                            action="store",
                            type=str,
                            required=True,
                            help=("The path to the tab delimited input "
                                  "methylation frequency or ready input case "
                                  "file(s) (First rwo is header which will be ignored). "
                                  "If multiple files, files must be "
                                  "in the same directory and give the path to the directory."))
    sdma_input.add_argument("--control", "-co",
                            action="store",
                            type=str,
                            required=True,
                            help=("The path to the tab delimited input "
                                  "methylation frequency or ready input "
                                  "control file(s) (First rwo is header which will be ignored). " 
                                  "If multiple files, files must be "
                                  "in the same directory and give the path to the directory."))
    sdma_input.add_argument("--out_dir", "-o",
                            action="store",
                            type=str,
                            required=True,
                            help="The path to the output directory")
    sdma_input.add_argument("--out_prefix", "-op",
                            action="store",
                            type=str,
                            required=True,
                            help="The prefix for the output files")
    sdma_input = sub_dma.add_argument_group("General optional arguments.")
    sdma_input.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    sdma_input.add_argument("--columns", "-c",
                            action="store",
                            type=str,
                            required=False,
                            help=("Comma seperated Columns in the methylation "
                                  "frequency files that include the following "
                                  "information, respectively:\n"
                                  "chromosome\tstart\tstrand\tcoverage\t"
                                  "methylation_frequency.\n"
                                  "If the methylation frequency file does not "
                                  "have strand level information then just "
                                  "enter columns number for\n"
                                  "chromosome\tstart\tcoverage\t"
                                  "methylation_frequency.\n"
                                  "Default is that your input files are "
                                  "already in a format required by DSS so you "
                                  "do not need to select any column.\n"
                                  "If you giving as input NanoMethPhase "
                                  "frequency files select this:"
                                  "--columns 1,2,4,5,7\n"))
    sdma_input.add_argument("--coverage", "-cov",
                            action="store",
                            type=int,
                            default=1,
                            required=False,
                            help=("Minimum coverage cutoff. Default is 1. It is "
                                  "recommended that do not filter for "
                                  "coverage as DSS R package will take care "
                                  "of it."))
    sdma_input.add_argument("--Rscript", "-rs",
                            action="store",
                            type=str,
                            required=False,
                            default="Rscript",
                            help="The path to a particular instance of "
                                 "Rscript to use.")
    sdma_input.add_argument("--script_file", "-sf",
                            action="store",
                            type=str,
                            required=False,
                            default=os.path.join(os.path.dirname(
                                                    os.path.realpath(__file__)
                                                        ),
                                                 "DSS_DMA.R"),
                            help="The path to the DSS_DMA.R script file."
                            " By default the script that was shipped during nanomethphase"
                            "installation will be used.")
    sdma_input.add_argument("--overwrite", "-ow",
                            action="store_true",
                            required=False,
                            help="If output files exist overwrite them.")
    sdma_input = sub_dma.add_argument_group("optional arguments that will be used in "
                                            "DSS DMLtest function.")
    sdma_input.add_argument("--smoothing_span", "-sms",
                            action="store",
                            type=int,
                            default=500,
                            required=False,
                            help=("The size of smoothing window, in "
                                  "basepairs. Default is 500."))
    sdma_input.add_argument("--smoothing_flag", "-smf",
                            action="store",
                            type=str,
                            default="TRUE",
                            required=False,
                            help=("TRUE/FALSE. A flag to indicate whether to apply smoothing"
                                  " in estimating mean methylation levels. It is "
                                  "recommended to use smoothing TRUE for "
                                  "whole-genome BS-seq data, and "
                                  "smoothing FALSE for sparser data such "
                                  "like from RRBS or hydroxyl-methylation "
                                  "data (TAB-seq). see -ed option and DSS R package details "
                                  " for more information. Default is TRUE."))
    sdma_input.add_argument("--equal_disp", "-ed",
                            action="store",
                            type=str,
                            default="FALSE",
                            required=False,
                            help=("TRUE/FALSE. A flag to indicate whether the "
                                  "dispersion in two groups are deemed equal or not. "
                                  "Default is FALSE "
                                  "and the dispersion shrinkages are performed "
                                  "on two conditions independently. "
                                  "More info on -ed and -smf: When there is no biological "
                                  "replicate in one or both treatment groups, "
                                  "users can either (1) specify "
                                  "equal.disp TRUE, which assumes both groups "
                                  "have the same dispersion, then the data "
                                  "from two groups are combined and used as "
                                  "replicates to estimate dispersion; or (2) "
                                  "specify smoothing TRUE, which uses the "
                                  "smoothed means (methylation levels) to "
                                  "estimate dispersions via a shrinkage "
                                  "estimator. This smoothing procedure uses "
                                  "data from neighboring CpG sites as "
                                  "\"pseudo-replicate\" for estimating "
                                  "biological variance."))
    sdma_input = sub_dma.add_argument_group("optional arguments that will be used in "
                                            "DSS callDML and callDMR functions.")
    sdma_input.add_argument("--pval_cutoff", "-pvc",
                            action="store",
                            type=float,
                            default=0.001,
                            required=False,
                            help=("A threshold of p-values for calling DMLs and DMRs. "
                                  "When delta is not specified, Loci with p-values less than this "
                                  "threshold will be picked as DML and also joint to form "
                                  "the DMRs. See DSS R package 'details' for more "
                                  "information for this regarding DMLs and DMRs. "
                                  "Default is 0.001."))
    sdma_input.add_argument("--delta_cutoff", "-dc",
                            action="store",
                            type=float,
                            default=0,
                            required=False,
                            help=("A threshold for defining DMLs and DMRs. In DML "
                                  "detection procedure, a hypothesis test "
                                  "that the two groups means are equal is "
                                  "conducted at each CpG site. Here if "
                                  "'delta' is specified, the function will "
                                  "compute the posterior probability that the "
                                  "difference of the means are greater than "
                                  "delta, and then call DML and construct DMR based on "
                                  "that. This only works when the test "
                                  "results are from 'DMLtest', which is for "
                                  "two-group comparison. See DSS R package for "
                                  "more information. Default is 0."))
    sdma_input = sub_dma.add_argument_group("optional arguments that will be used in "
                                            "DSS callDMR function.")
    sdma_input.add_argument("--minlen", "-ml",
                            action="store",
                            type=int,
                            default=100,
                            required=False,
                            help=("Minimum length (in basepairs) required for "
                                  "DMR. Default is 100 bps."))
    sdma_input.add_argument("--minCG", "-mcg",
                            action="store",
                            type=int,
                            default=15,
                            required=False,
                            help=("Minimum number of CpG sites required for "
                                  "DMR. Default is 15."))
    
    sdma_input.add_argument("--pct_sig", "-pct",
                            action="store",
                            type=float,
                            default=0.5,
                            required=False,
                            help=("In all DMRs, the percentage of CG sites "
                                  "with significant p-values (less than "
                                  "p.threshold) must be greater than this "
                                  "threshold. Default is 0.5. This is mainly "
                                  "used for correcting the effects of merging "
                                  "of nearby DMRs."))
    sdma_input.add_argument("--dis_merge", "-dm",
                            action="store",
                            type=int,
                            default=100,
                            required=False,
                            help=("When two DMRs are very close to each other "
                                  "and the distance (in bps) is less than "
                                  "this number, they will be merged into one. "
                                  "Default is 100 bps."
                                  " See dma section notes for more details."))
    sub_dma.set_defaults(func=main_dma)


def main():
    """
    Docstring placeholder.
    """
    parser = argparse.ArgumentParser(
        prog="nanomethphase",
        description="NanoMethPhase: For phasing Nanopore Reads and "
                    "Methylation")
    parser.add_argument('--version', action='version', version='%(prog)s 1.1.0_dev')
    subparsers = parser.add_subparsers(title="Modules")
    methyl_call_processor_parser(subparsers)
    phase_parser(subparsers)
    dma_parser(subparsers)
    bam2bis_parser(subparsers)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
