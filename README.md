NanoMethPhase
=============

Phase long reads and CpG methylations from Oxford Nanopore Technologies.

## Installation
Using [pypi repository](https://pypi.org/project/nanomethphase/)

```
pip install nanomethphase
```
**NOTE:** NanoMethPhase needs python 3.7  
Using [Docker image](https://hub.docker.com/r/jmgarant/nanomethphase)

It ships with complementary softwares SNVoter, Nanopolish, Clair and WhatsHap.

```bash
docker pull jmgarant/nanomethphase

# usage example:
docker run -t jmgarant/nanomethphase nanomethphase
docker run -t jmgarant/nanomethphase snvoter
docker run -t jmgarant/nanomethphase nanopolish
docker run -t jmgarant/nanomethphase clair
docker run -t jmgarant/nanomethphase whatsapp
```

From [source](https://github.com/vahidAK/NanoMethPhase.git)

```
git clone https://github.com/vahidAK/NanoMethPhase.git
cd NanoMethPhase
./nanomethphase.py
```

## Creating a dedicated conda environment

Environment file available in the
[git repository](https://github.com/vahidAK/NanoMethPhase)

```
git clone https://github.com/vahidAK/NanoMethPhase.git
conda env create -f NanoMethPhase/envs/environment.yaml
```

# Quickstart

If you have your methylation call data and phased vcf file you can get the
haplotype methylome via:

1- Processing and indexing methylation call file

`nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > MethylationCall.bed.gz && tabix -p bed MethylationCall.bed.gz`

2- Getting haplotype methylome:

`nanomethphase  phase -mc MethylationCall.bed.gz -o Test_methylome -of bam,methylcall,bam2bis -b sorted.bam -r hg38.fa -v Phased.vcf -t 64`

You can select 3 output options:

bam: output phased bam files

methylcall: output phased methylation call and frequency files

bam2bis: output mock whole-genome bisulfite converted bam files

# Full Tutorial

In order to get the phased methylome you also need the following third-party
software:

[Nanopolish](https://github.com/jts/nanopolish) : To call CpG methylation.

[Clair](https://github.com/HKU-BAL/Clair) or other variant callers: To call
variants for your sample. Alternatively, you might already have variant calling
data for example from Illumina sequencing.

[WhatsHap](https://github.com/whatshap/whatshap): To phase single nucleotide
variants.

## Methylation Calling

1- indexing fastq file and fast5 files:

NOTE: Fastqs must be merged to a single file

`nanopolish index -d /path/to/fast5s_directory/.fastq`

2- Methylation calling for CpG from each read:

`nanopolish call-methylation -t <number of threads> -q cpg -r /path/to/fastq_fromstep-1/fastq.fastq -b /path/to/sorted_and_indexed/bam.bam -g /path/to/reference.fa > /path/to/MethylationCall.tsv`

For the full tutorial please refer to
[Nanopolish](https://github.com/jts/nanopolish) page on GitHub.

## Variant Calling

We have used Clair to call variants. However, you may call variants with other
tools or your variant data may come from Illumina or other methods.

You can call variants for each chromosome using the following command and the
concatenate all files:

`for i in chr{1..22} chrX chrY; do callVarBam --chkpnt_fn <path to model file> --ref_fn <reference_genome.fa> --bam_fn <sorted_indexed.bam> --ctgName $i --sampleName <your sample name> --call_fn $i".vcf" --threshold 0.2 --samtools <path to executable samtools software> --pypy <path to executable pypy > --threads <number of threads>`

For the full tutorial please refer to [Clair](https://github.com/HKU-BAL/Clair)
page on GitHub.

After variant calling, you can select only SNVs which will be used for phasing:
`awk '$4 != "." && $5 != "." && length($4) == 1 && length($5) == 1' && $6 > <the variant calling quality threshold> variants.vcf > HighQualitySNVs.vcf`

If you are calling variants from low coverage nanopore data (<30x) using Clair, you can also use our other tool [SNVoter](https://github.com/vahidAK/SNVoter) to improve SNV detection.

## Phasing of detected SNVs

If you have your SNVs data available you need to phase them using
[WhatsHap](https://github.com/whatshap/whatshap).

`whatshap phase --ignore-read-groups --reference reference.fa -o HighQualitySNVs_whatshap_phased.vcf HighQualitySNVs.vcf sorted_indexed.bam`

For the full tutorial please refer to
[WhatsHap](https://github.com/whatshap/whatshap) page on GitHub.

If you have Trio data (Father, Mother, Child) you can use the script
[Trio_To_PhaseVCF_4FemaleChild.sh](https://github.com/vahidAK/NanoMethPhase/tree/master/scripts)
or
[Trio_To_PhaseVCF_4MaleChild.sh](https://github.com/vahidAK/NanoMethPhase/tree/master/scripts)
script to make a mock phased vcf file and use it as input for NanoMethPhase.

## Detecting Haplotype Methylome

1- First you need to phase process methylation call file from Nanopolish.

`nanomethphase methyl_call_processor -mc MethylationCall.tsv -t 20 | sort -k1,1 -k2,2n -k3,3n | bgzip > MethylationCall.bed.gz && tabix -p bed MethylationCall.bed.gz`

2- Getting haplotype methylome:

`nanomethphase  phase -mc MethylationCall.bed.gz -o Test_methylome -of bam,methylcall,bam2bis -b sorted.bam -r hg38.fa -v Phased.vcf -t 64`

If your are not using called SNVs from nanopore data, and they come from, for
example, short-read sequencing, we recommend using -mbq 0 in the above code. 

You can select 3 output options:

bam: output phased bam files

methylcall: output phased methylation call and frequency files

bam2bis: output mock whole-genome bisulfite converted bam files

3- Differential Methylation Analysis:

`nanomethphase dma -c 1,2,4,5,7 -ca <path to output methylation frequency for haplotype1> -co <path to output methylation frequency for haplotype2> -rs <path to executable Rscript> -sf /projects/vakbari_prj/scratch/BitBucket/NanoMethPhase/DSS_DMA.R -o . -op Aicardi_NanoMethPhase_DMA`

To findout about different modules run:

`nanomethphase -h`

For a full list of options and help for each module run:

`nanomethphase <module name> -h`
  
We have included an example data in the Example_Data folder which you can use for a quick detection of haplotype methylome on 1Mb of chr21.
