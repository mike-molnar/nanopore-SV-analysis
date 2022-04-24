import os
import re

configfile: "config.yaml"

# Get the path information for the programs and scripts
root_dir = os.getcwd()
conda_dir = config["conda_dir"]
scripts_dir = config["workflow_dir"] + "/scripts"

# Declare the variables and files needed for the reference genome
reference = config["workflow_dir"] + "/reference/GRCh38.fa"
chromosome_sizes = config["workflow_dir"] + "/reference/GRCh38_chromosome_sizes.tsv"
high_frequency_kmers = config["workflow_dir"] + "/reference/GRCh38_high_frequency_kmers.txt"
genome_gaps = config["workflow_dir"] + "/reference/GRCh38_gaps.bed"
genes = config["workflow_dir"] + "/reference/GRCh38_genes.bed"
genome_length = 3100000000
chromosomes = [chr1 ,chr2 ,chr3 ,chr4 ,chr5 ,chr6 ,chr7 ,chr8 ,chr9 ,chr10 ,chr11 ,chr12 ,chr13,
               chr14 ,chr15 ,chr16 ,chr17 ,chr18 ,chr19 ,chr20 ,chr21 ,chr22 ,chrX ,chrY]

# Get the mean coverage for the dataset
def get_coverage(sample_name):
    file_name = str(sample_name) + "/analysis/nanoplot/" + str(sample_name) + "NanoStats.txt"
    if os.path.exists(file_name):
        f = open(file_name)
        pattern = "Total bases"
        for line in f:
            if re.search(pattern, line):
                total_bases = line.split(":")[1].strip().replace(',', '')
                return float(total_bases)/3100000000

# Get the list of regions for the workflow
def get_regions():
    f = open(config["workflow_dir"] + "/reference/GRCh38_regions_list.txt")
    regions = list()
    for line in f:
        regions.append(line.rstrip())
    return regions
    
regions = get_regions()

# Include the rules for the workflow
include: config["workflow_dir"] + "/rules/analysis.smk"
include: config["workflow_dir"] + "/rules/mapping.smk"
include: config["workflow_dir"] + "/rules/structural_variants.smk"

# Default rule to make everything
rule all:
    input:
        expand("{sample}/analysis/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.translocations.bedpe", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.CNVs.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.1000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.small_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.large_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosome}.pdf", normal=config["normals"], sample=config["tumors"], chromosome=chromosomes)

rule filter_fastq:
    input:
        expand("{sample}/fastq/{sample}.fastq.gz", sample=config["samples"])

rule mapping:
    input:
        expand("{sample}/mapped/{sample}.b_allele_frequency.bed", sample=config["samples"])
        
rule run_SV_analysis:
    input:
        expand("{sample}/analysis/structural_variants/{sample}.insertions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.deletions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.duplications.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.inversions.bed", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.translocations.bedpe", sample=config["samples"]),
        expand("{sample}/analysis/structural_variants/{sample}.CNVs.bed", sample=config["samples"])

rule run_coverage_analysis:
    input:
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.1000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed", sample=config["samples"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.small_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.large_chr.pdf", normal=config["normals"], sample=config["tumors"]),
        expand("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosome}.pdf", normal=config["normals"], sample=config["tumors"], chromosome=chromosomes)
