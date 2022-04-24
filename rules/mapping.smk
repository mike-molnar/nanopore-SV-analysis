# Filter out short or low quality reads
rule nanofilt:
    input:
        "{sample}/fastq/{sample}.fastq"
    output:
        protected("{sample}/fastq/{sample}.fastq.gz")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        quality="7",
        min_length="500"
    log:
        "{sample}/analysis/logs/fastq/{sample}.nanofilt.log"
    benchmark:
        "{sample}/analysis/benchmarks/fastq/{sample}.nanofilt.txt"
    threads: 1
    shell:
        "{conda_dir}/NanoFilt --logfile {log} -q {params.quality} -l {params.min_length} \
        {input} | {conda_dir}/bgzip -c > {output} 2>> {log}"

# Index the zipped fastq file
rule index_fastq:
    input:
        "{sample}/fastq/{sample}.fastq.gz"
    output:
        "{sample}/fastq/{sample}.fastq.gz.gzi"
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/fastq/{sample}.bgzip_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/fastq/{sample}.bgzip_index.txt"
    threads: 1
    shell:
        "{conda_dir}/bgzip -r {input} 2> {log}"

# Map the reads to the reference genome
rule winnowmap:
    input:
        fastq = "{sample}/fastq/{sample}.fastq.gz",
        index = "{sample}/fastq/{sample}.fastq.gz.gzi"
    output:
        temp("{sample}/mapped/{sample}.bam")
    params:
        memory_per_thread="8G",
        run_time="7:0:0:0",
        preset_options="-ax map-ont",
        include_MD_tag="--MD"
    log:
        "{sample}/analysis/logs/mapped/{sample}.winnowmap.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.winnowmap.txt"
    threads: 16
    shell:
        "{conda_dir}/winnowmap {params.preset_options} {params.include_MD_tag} \
        -t {threads} -W {high_frequency_kmers} {reference} {input.fastq} 2> {log} | \
        {conda_dir}/samtools sort -o {output} &>> {log}"

# Index the mapped reads
rule winnowmap_index:
    input:
        "{sample}/mapped/{sample}.bam"
    output:
        temp("{sample}/mapped/{sample}.bam.bai")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.winnowmap_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.winnowmap_index.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools index {input} &> {log}"
        
# Find SNVs from the mapped reads
rule longshot:
    input:
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        temp("{sample}/temp_files/{sample}.{regions}.vcf")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        thresholds="-c 4 -C 80 -e 1 -I 7000"
    log:
        "{sample}/analysis/logs/temp_files/{sample}.{regions}.longshot.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/{sample}.{regions}.longshot.txt"
    threads: 1
    shell:
        "{conda_dir}/longshot {params.thresholds} -r {wildcards.regions} --bam {input.bam} \
        --ref {reference} --out {output} &> {log}"

# Zip the VCFs from longshot
rule zip_vcf_region:
    input:
        "{sample}/temp_files/{sample}.{regions}.vcf"
    output:
        temp("{sample}/temp_files/{sample}.{regions}.vcf.gz")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/temp_files/{sample}.{regions}.zip_vcf_region.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/{sample}.{regions}.zip_vcf_region.txt"
    threads: 1
    shell:
        "{conda_dir}/bgzip --threads {threads} -c {input} > {output} 2> {log}"
        
# Index the zipped VCFs
rule index_vcf_region:
    input:
        "{sample}/temp_files/{sample}.{regions}.vcf.gz"
    output:
        temp("{sample}/temp_files/{sample}.{regions}.vcf.gz.tbi")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        file_type="vcf"
    log:
        "{sample}/analysis/logs/temp_files/{sample}.{regions}.index_vcf_region.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/{sample}.{regions}.index_vcf_region.txt"
    threads: 1
    shell:
        "{conda_dir}/tabix -p {params.file_type} {input} &> {log}"
                
# Merge the vcf files
rule merge_longshot_vcf:
    input:
        zip = expand("{{sample}}/temp_files/{{sample}}.{regions}.vcf.gz", regions=regions),
        index = expand("{{sample}}/temp_files/{{sample}}.{regions}.vcf.gz.tbi", regions=regions)
    output:
        protected("{sample}/mapped/{sample}.longshot.vcf.gz")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        concat="-a -O v",
        sort="-O z"
    log:
        "{sample}/analysis/logs/mapped/{sample}.merge_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.merge_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/bcftools concat {params.concat} {input.zip} 2> {log} | \
        {conda_dir}/bcftools sort {params.sort} -o {output} &>> {log}"
        
# Merge the vcf files
rule index_longshot_vcf:
    input:
        "{sample}/mapped/{sample}.longshot.vcf.gz"
    output:
        protected("{sample}/mapped/{sample}.longshot.vcf.gz.tbi")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        file_type="vcf"
    log:
        "{sample}/analysis/logs/mapped/{sample}.index_longshot_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.index_longshot_vcf.txt"
    threads: 1
    shell:
        "{conda_dir}/tabix -p {params.file_type} {input} &> {log}"

# Tag the bam file with haplotype information
rule haplotag_bam:
    input:
        vcf = "{sample}/mapped/{sample}.longshot.vcf.gz",
        tbi = "{sample}/mapped/{sample}.longshot.vcf.gz.tbi",
        bam = "{sample}/mapped/{sample}.bam",
        bai = "{sample}/mapped/{sample}.bam.bai"
    output:
        temp("{sample}/temp_files/whatshap/{regions}/{sample}.{regions}.phased.bam")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0",
        RG="--ignore-read-groups"
    log:
        "{sample}/analysis/logs/temp_files/whatshap/{regions}/{sample}.{regions}.phased.log"
    benchmark:
        "{sample}/analysis/benchmarks/temp_files/whatshap/{regions}/{sample}.{regions}.phased.txt"
    threads: 1
    shell:
        "{conda_dir}/whatshap haplotag {params.RG} --regions {wildcards.regions} --reference {reference} \
        {input.vcf} {input.bam} 2> {log} | {conda_dir}/samtools sort -o {output} &>> {log}"

# Merge the phased bam files
rule merge_bam:
    input:
        expand("{{sample}}/temp_files/whatshap/{regions}/{{sample}}.{regions}.phased.bam", regions=regions)
    output:
        protected("{sample}/mapped/{sample}.phased.bam")
    params:
        memory_per_thread="16G",
        run_time="1:0:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.merge_bam.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.merge_bam.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools merge {output} {input} &> {log}"
        
# Index the phased bam file
rule haplotag_bam_index:
    input:
        "{sample}/mapped/{sample}.phased.bam"
    output:
        protected("{sample}/mapped/{sample}.phased.bam.bai")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.haplotag_index.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.haplotag_index.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools index {input} &> {log}"
                 
# Extract the allele counts from the vcf
rule extract_longshot_vcf:
    input:
        zip = "{sample}/mapped/{sample}.longshot.vcf.gz",
        index = "{sample}/mapped/{sample}.longshot.vcf.gz.tbi"
    output:
        temp("{sample}/mapped/{sample}.longshot.tsv")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        view_param="'PASS'",
        query_param="'%CHROM\t%POS\t%AC\n'"
    log:
        "{sample}/analysis/logs/mapped/{sample}.extract_longshot_vcf.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.extract_longshot_vcf.txt"
    threads: 1
    shell:
        "zcat {input.zip} 2> {log} | {conda_dir}/bcftools view -f {params.view_param} \
        2>> {log} | {conda_dir}/bcftools query -f {params.query_param} > {output} 2>> {log}"
        
# Calcualte the allele counts from the vcf
rule calculate_b_allele_frequency:
    input:
        "{sample}/mapped/{sample}.longshot.tsv"
    output:
        protected("{sample}/mapped/{sample}.b_allele_frequency.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.calculate_b_allele_frequency.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.calculate_b_allele_frequency.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/calculate_b_allele_frequency.py -i {input} -o {output} 2> {log}"
       
# Find reads with low mapping quality
rule find_low_mapped_reads:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.sam")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        quality="30"
    log:
        "{sample}/analysis/logs/mapped/{sample}.find_low_mapped_reads.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.find_low_mapped_reads.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools view -h {input.bam} 2> {log} | \
        awk '$5<{params.quality} {{print $0}}' > {output} 2>> {log}"

# Find the coverage of regions with low mapping quality reads
rule find_low_mapping_coverage:
    input:
        "{sample}/mapped/{sample}.low_mapped_reads.sam"
    output:
        temp("{sample}/mapped/{sample}.low_mapped_reads.cov"),
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/mapped/{sample}.find_low_mapping_coverage.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.find_low_mapping_coverage.txt"
    threads: 1
    shell:
        "{conda_dir}/samtools view -S -b -h {input} 2> {log} | \
        {conda_dir}/samtools depth - > {output} 2>> {log}"
      
# Make a bed file of the regions with low mapping quality      
rule make_low_mapped_bed:
    input:
        "{sample}/mapped/{sample}.low_mapped_reads.cov"
    output:
        protected("{sample}/mapped/{sample}.low_mapping_regions.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        slop="500",
        merge_length="1000",
        distance="100",
        min_coverage="2"
    log:
        "{sample}/analysis/logs/mapped/{sample}.make_low_mapped_bed.log"
    benchmark:
        "{sample}/analysis/benchmarks/mapped/{sample}.make_low_mapped_bed.txt"
    threads: 1
    shell:
        "{conda_dir}/SURVIVOR bincov {input} {params.distance} {params.min_coverage} 2> {log} | \
        {conda_dir}/bedtools slop -i stdin -g {chromosome_sizes} -b {params.slop} 2>> {log} | \
        {conda_dir}/bedtools merge -i stdin -d {params.merge_length} > {output} 2>> {log}"
