NanoMethPhase
=============

Phase long reads and CpG methylations from Oxford Nanopore Technologies.

## Installation

Using pypi repository

```
pip install nanomethphase
```

Using conda

```
TBD
```

## Creating a dedicated conda environment

Environment file available in the [git repository](https://svn.bcgsc.ca/bitbucket/users/vakbari/repos/nanomethphase/browse)

```
git clone https://jmgarant@svn.bcgsc.ca/bitbucket/scm/~vakbari/nanomethphase.git
cd nanomethphase
conda env create -f ens/environment.yaml
```  
  
# Quickstart
If you have your methylation call data and phased vcf file you can get the haplotype methylome via:  
  
1- Processing and indexing methylation call file  
  
`nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > MethylationCall.bed.gz && tabix -p bed MethylationCall.bed.gz`
  
2- Getting haplotype methylome:  
  
 `nanomethphase  phase -mc MethylationCall.bed.gz -o Test_methylome -of bam,methylcall,bam2bis -b sorted.bam -r hg38.fa -v Phased.vcf -t 64`  
  
You can select 3 output options:  

bam: output phased bam files  
methylcall: output phased methylation call and frequency files  
bam2bis: output mock whole-genome bisulfite converted bam files  

# Full Tutorial

In order to get the phased methylome you also need the following third-party software:  
[Nanopolish](https://github.com/jts/nanopolish) : To call CpG methylation.  
[Clair](https://github.com/HKU-BAL/Clair) or other variant callers: To call variants for your sample. Alternatively, you might already have variant calling data for example from Illumina sequencing.  
[WhatsHap](https://github.com/whatshap/whatshap): To phase single nucleotide variants.  

## Methylation Calling
1- indexing fastq file and fast5 files:  
NOTE: Fastqs must be merged to a single file  
`nanopolish index -d /path/to/fast5s_directory/.fastq`  

2- Methylation calling for CpG from each read:  
`nanopolish call-methylation -t <number of threads> -q cpg -r /path/to/fastq_fromstep-1/fastq.fastq -b /path/to/sorted_and_indexed/bam.bam -g /path/to/reference.fa > /path/to/MethylationCall.tsv`  

## Phasing

