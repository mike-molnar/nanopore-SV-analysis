# Get FASTQ stats from NanoPlot
rule nanoplot:
    input:
        "{sample}/fastq/{sample}.fastq.gz"
    output:
        protected("{sample}/analysis/nanoplot/{sample}LengthvsQualityScatterPlot_kde.png"),
        protected("{sample}/analysis/nanoplot/{sample}NanoStats.txt")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        scale="--loglength"
    log:
        "{sample}/analysis/logs/analysis/{sample}.nanoplot.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/{sample}.nanoplot.txt"
    threads: 4
    shell:
        "{conda_dir}/NanoPlot {params.scale} -t {threads} -p {wildcards.sample} \
        --title {wildcards.sample} --fastq {input} -o {wildcards.sample}/analysis/nanoplot &> {log}"
        

# Find depth of coverage at each genomic position in the reference       
rule samtools_depth:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        temp("{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/temp_files/samtools_depth/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/samtools_depth/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools depth -r {wildcards.chromosomes} {input.bam} > {output} 2> {log}"

# Calculate average depth of coverage in 500,000bp windows
rule coverage_by_window_500000:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_500000/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="500000",
        coverage = lambda wildcards: get_coverage(wildcards.sample)
    log:
        "{sample}/analysis/logs/temp_files/samtools_depth/window_500000/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/samtools_depth/window_500000/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -c {params.coverage} &> {log}"
        
rule cat_coverage_by_window_500000:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_500000/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/analysis/coverage/samtools_depth/{sample}.depth.500000_window.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.500000_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

# Calculate average depth of coverage in 100,000bp windows
rule coverage_by_window_100000:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_100000/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="100000",
        coverage = lambda wildcards: get_coverage(wildcards.sample)
    log:
        "{sample}/analysis/logs/temp_files/samtools_depth/window_100000/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/samtools_depth/window_100000/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -c {params.coverage} &> {log}"
        
rule cat_coverage_by_window_100000:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_100000/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100000_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/analysis/coverage/samtools_depth/{sample}.depth.100000_window.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.100000_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

# Calculate average depth of coverage in 10,000bp windows
rule coverage_by_window_10000:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_10000/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="10000",
        coverage = lambda wildcards: get_coverage(wildcards.sample)
    log:
        "{sample}/analysis/logs/temp_files/samtools_depth/window_10000/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/samtools_depth/window_10000/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -c {params.coverage} &> {log}"
        
rule cat_coverage_by_window_10000:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_10000/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/analysis/coverage/samtools_depth/{sample}.depth.10000_window.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.10000_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

# Calculate average depth of coverage in 1,000bp windows
rule coverage_by_window_1000:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_1000/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="1000",
        coverage = lambda wildcards: get_coverage(wildcards.sample)
    log:
        "{sample}/analysis/logs/temp_files/samtools_depth/window_1000/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/samtools_depth/window_1000/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -c {params.coverage} &> {log}"
        
rule cat_coverage_by_window_1000:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_1000/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.1000_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/analysis/coverage/samtools_depth/{sample}.depth.1000_window.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.1000_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

# Calculate average depth of coverage in 100bp windows
rule coverage_by_window_100:
    input:
        depth = "{sample}/temp_files/samtools_depth/{sample}.{chromosomes}.tsv",
        stats = "{sample}/analysis/nanoplot/{sample}NanoStats.txt"
    output:
        temp("{sample}/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        window="100",
        coverage = lambda wildcards: get_coverage(wildcards.sample)
    log:
        "{sample}/analysis/logs/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/samtools_depth/window_100/{sample}.{chromosomes}.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/coverage_by_window.py -i {input.depth} \
        -o {output} -w {params.window} -c {params.coverage} &> {log}"

rule cat_coverage_by_window_100:
    input:
        expand("{{sample}}/temp_files/samtools_depth/window_100/{{sample}}.{chromosomes}.bed", chromosomes=chromosomes)
    output:
        protected("{sample}/analysis/coverage/samtools_depth/{sample}.depth.100_window.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/analysis/coverage/samtools_depth/{sample}.depth.100_window.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/samtools_depth/{sample}.depth.100_window.txt"
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"

rule plot_genome_coverage:
    input:
        normal = "{normal}/analysis/coverage/samtools_depth/{normal}.depth.10000_window.bed",
        sample = "{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed"
    output:
        protected("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.pdf")
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        genome="hg38"
    log:
        "{sample}/analysis/logs/analysis/coverage/plots/{normal}/{sample}.depth.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/plots/{normal}/{sample}.depth.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploter_coverage.R {input.normal} {input.sample} \
        {output} {params.genome} {wildcards.normal} {wildcards.sample} &> {log}"

rule plot_genome_coverage_small_chr:
    input:
        normal = "{normal}/analysis/coverage/samtools_depth/{normal}.depth.10000_window.bed",
        sample = "{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed"
    output:
        protected("{sample}/analysis/coverage/plots/{normal}/{sample}.small_chr.pdf")
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        genome="hg38"
    log:
        "{sample}/analysis/logs/analysis/coverage/plots/{normal}/{sample}.small_chr.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/plots/{normal}/{sample}.small_chr.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploter_coverage_small_chr.R {input.normal} {input.sample} \
        {output} {params.genome} {wildcards.normal} {wildcards.sample} &> {log}"
        
rule plot_genome_coverage_large_chr:
    input:
        normal = "{normal}/analysis/coverage/samtools_depth/{normal}.depth.10000_window.bed",
        sample = "{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed"
    output:
        protected("{sample}/analysis/coverage/plots/{normal}/{sample}.large_chr.pdf")
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        genome="hg38"
    log:
        "{sample}/analysis/logs/analysis/coverage/plots/{normal}/{sample}.large_chr.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/plots/{normal}/{sample}.large_chr.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploter_coverage_large_chr.R {input.normal} {input.sample} \
        {output} {params.genome} {wildcards.normal} {wildcards.sample} &> {log}"
        
rule plot_chr_coverage:
    input:
        normal_cov = "{normal}/analysis/coverage/samtools_depth/{normal}.depth.10000_window.bed",
        sample_cov = "{sample}/analysis/coverage/samtools_depth/{sample}.depth.10000_window.bed",
        normal_baf = "{normal}/mapped/{normal}.b_allele_frequency.bed",
        sample_baf = "{sample}/mapped/{sample}.b_allele_frequency.bed"
    output:
        protected("{sample}/analysis/coverage/plots/{normal}/{sample}.depth.{chromosomes}.pdf")
    params:
        memory_per_thread="32G",
        run_time="0:1:0:0",
        genome="hg38"
    log:
        "{sample}/analysis/logs/analysis/coverage/plots/{normal}/{sample}.{chromosomes}.plot_chr_coverage.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/coverage/plots/{normal}/{sample}.{chromosomes}.plot_chr_coverage.txt"
    threads: 1
    shell:
        "{conda_dir}/Rscript {scripts_dir}/karyoploter_coverage_chr.R {input.normal_cov} {input.sample_cov} \
        {input.normal_baf} {input.sample_baf} {wildcards.normal} {wildcards.sample} {params.genome} \
        {wildcards.chromosomes} {output} &> {log}"
        
# Filter structural variants
rule filter_structural_variants:
    input:
        INS = "{sample}/structural_variants/{sample}.insertions.bed",
        DEL = "{sample}/structural_variants/{sample}.deletions.bed",
        DUP = "{sample}/structural_variants/{sample}.duplications.bed",
        INV = "{sample}/structural_variants/{sample}.inversions.bed",
        trans = "{sample}/structural_variants/{sample}.translocations.bedpe",
        depth = "{sample}/analysis/coverage/samtools_depth/{sample}.depth.500000_window.bed",
        split_alignments = "{sample}/structural_variants/{sample}.split_alignments.bed",
        low_map = "{sample}/mapped/{sample}.low_mapping_regions.bed"
    output:
        INS = "{sample}/analysis/structural_variants/{sample}.insertions.bed",
        DEL = "{sample}/analysis/structural_variants/{sample}.deletions.bed",
        DUP = "{sample}/analysis/structural_variants/{sample}.duplications.bed",
        INV = "{sample}/analysis/structural_variants/{sample}.inversions.bed",
        trans = "{sample}/analysis/structural_variants/{sample}.translocations.bedpe",
        CNVs = "{sample}/analysis/structural_variants/{sample}.CNVs.bed"
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        coverage = lambda wildcards: get_coverage(wildcards.sample)
    log:
        "{sample}/analysis/logs/analysis/structural_variants/filter_structural_variants.log"
    benchmark:
        "{sample}/analysis/benchmarks/analysis/structural_variants/filter_structural_variants.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/filter_SVs.py -ins {input.INS} -del {input.DEL} -dup {input.DUP} \
        -inv {input.INV} -trans {input.trans} -depth {input.depth} -split {input.split_alignments} \
        -low_map {input.low_map} -gaps {genome_gaps} -ins_out {output.INS} -del_out \
        {output.DEL} -dup_out {output.DUP} -inv_out {output.INV} -trans_out {output.trans} -cnv_out \
        {output.CNVs} -cov {params.coverage} &> {log}"
        
