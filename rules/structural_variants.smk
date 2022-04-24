# Find structural variants with cuteSV
rule cuteSV:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        protected("{sample}/structural_variants/cuteSV/{sample}.cuteSV.vcf")
    params:
        memory_per_thread="4G",
        run_time="1:0:0:0",
        working_dir="{sample}/structural_variants/cuteSV/",
        min_support="3",
        max_splits="2",
        max_cluster_bias_DEL="100",
        diff_ratio_merging_DEL="0.3"
    log:
        "{sample}/analysis/logs/structural_variants/cuteSV/{sample}.cuteSV.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/cuteSV/{sample}.cuteSV.txt"
    threads: 8
    shell:
        "{conda_dir}/cuteSV --diff_ratio_merging_DEL {params.diff_ratio_merging_DEL} \
        --min_support {params.min_support} --max_split_parts {params.max_splits} \
        --threads {threads} --max_cluster_bias_DEL {params.max_cluster_bias_DEL} {input.bam} \
        {reference} {output} {params.working_dir} &> {log}"

# Convert the cuteSV VCF file to a standardized BEDPE
rule vcf2bedpe_cuteSV:
    input:
        "{sample}/structural_variants/cuteSV/{sample}.cuteSV.vcf"
    output:
        protected("{sample}/structural_variants/cuteSV/{sample}.cuteSV.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/structural_variants/cuteSV/{sample}.vcf2bedpe_cuteSV.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/cuteSV/{sample}.vcf2bedpe_cuteSV.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/vcf2bedpe_cuteSV.py {input} {output} &> {log}"

# Find structural variants with sniffles
rule sniffles:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        snf = protected("{sample}/structural_variants/sniffles/{sample}.sniffles.snf"),
        vcf = protected("{sample}/structural_variants/sniffles/{sample}.sniffles.vcf")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0",
        min_reads = "2",
        max_splits = "1"
    log:
        "{sample}/analysis/logs/structural_variants/sniffles/{sample}.sniffles.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/sniffles/{sample}.sniffles.txt"
    threads: 8
    shell:
        "{conda_dir}/sniffles --input {input.bam} --snf {output.snf} --vcf {output.vcf} --non-germline \
         --minsupport {params.min_reads} --max-splits-base {params.max_splits} -t {threads} &> {log}"
                
# Standardize the sniffles BEDPE
rule vcf2bedpe_sniffles:
    input:
        "{sample}/structural_variants/sniffles/{sample}.sniffles.vcf"
    output:
        protected("{sample}/structural_variants/sniffles/{sample}.sniffles.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/structural_variants/svim/{sample}.vcf2bedpe_svim.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/svim/{sample}.vcf2bedpe_svim.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/vcf2bedpe.py {input} {output} &> {log}"
        
# Find structural variants with SVIM
rule svim:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        protected("{sample}/structural_variants/svim/variants.vcf")
    params:
        memory_per_thread="32G",
        run_time="1:0:0:0",
        working_dir="{sample}/structural_variants/svim"
    log:
        "{sample}/analysis/logs/structural_variants/svim/{sample}.svim.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/svim/{sample}.svim.txt"
    threads: 1
    shell:
        "{conda_dir}/svim alignment {params.working_dir} {input.bam} {reference} &> {log}"

# Convert the SVIM VCF file to a standardized BEDPE     
rule vcf2bedpe_svim:
    input:
        "{sample}/structural_variants/svim/variants.vcf"
    output:
        protected("{sample}/structural_variants/svim/{sample}.svim_all_calls.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0"
    log:
        "{sample}/analysis/logs/structural_variants/svim/{sample}.vcf2bedpe_svim.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/svim/{sample}.vcf2bedpe_svim.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/vcf2bedpe.py {input} {output} &> {log}"

# Normalize the SVIM calls to match sniffles and cuteSV 
rule normalize_svim_calls:
    input:
        "{sample}/structural_variants/svim/{sample}.svim_all_calls.bedpe"
    output:
        protected("{sample}/structural_variants/svim/{sample}.svim.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk="'NR==1 || $9>1'",
        fix_columns="'{{if($1 > $4) {{print $4,$5,$6,$1,$2,$3,$7,$8,$9}} else {{print $0}}}}' OFS='\t'",
        remove_duplicates="'!seen[$0]++'"
    log:
        "{sample}/analysis/logs/structural_variants/svim/{sample}.normalize_svim_calls.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/svim/{sample}.normalize_svim_calls.txt"
    threads: 1
    shell:
        "awk {params.awk} {input} 2> {log} | awk {params.fix_columns} 2>> {log} | \
        awk {params.remove_duplicates} > {output} 2>> {log}"

#=================================================================
# Merge SVs from the three callers into files for each type of SV
#=================================================================

# Merge the insertions into a single file
rule make_INS_bed_file:
    input:
        cuteSV = "{sample}/structural_variants/cuteSV/{sample}.cuteSV.bedpe",
        sniffles = "{sample}/structural_variants/sniffles/{sample}.sniffles.bedpe",
        svim = "{sample}/structural_variants/svim/{sample}.svim.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.insertions.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk_INS="'$7==\"INS\" {{print $1,$2,$3,$7,$8,$9}}' OFS='\t'",
        sort="-k1,1 -k2,2n -k3,3n",
        merge_length="100",
        columns="-c 4,5,6,6",
        output_type="-o distinct,mean,mean,count"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.make_INS_bed_file.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.make_INS_bed_file.txt"
    threads: 1
    shell:
        "awk {params.awk_INS} {input.cuteSV} {input.sniffles} {input.svim} 2> {log} | \
        sort {params.sort} 2>> {log} | {conda_dir}/bedtools merge -i stdin -d {params.merge_length} \
        {params.columns} {params.output_type} > {output} 2>> {log}"

# Merge the deletions into a single file
rule make_DEL_bed_file:
    input:
        cuteSV = "{sample}/structural_variants/cuteSV/{sample}.cuteSV.bedpe",
        sniffles = "{sample}/structural_variants/sniffles/{sample}.sniffles.bedpe",
        svim = "{sample}/structural_variants/svim/{sample}.svim.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.deletions.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk_DEL="'$7==\"DEL\" {{print $1,$2,$6,$7,$8,$9}}' OFS='\t'",
        sort="-k1,1 -k2,2n -k3,3n",
        merge_length="100",
        columns="-c 4,5,6,6",
        output_type="-o distinct,mean,mean,count"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.make_DEL_bed_file.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.make_DEL_bed_file.txt"
    threads: 1
    shell:
        "awk {params.awk_DEL} {input.cuteSV} {input.sniffles} {input.svim} 2> {log} | \
        sort {params.sort} 2>> {log} | {conda_dir}/bedtools merge -i stdin -d {params.merge_length} \
        {params.columns} {params.output_type} > {output} 2>> {log}"

# Merge the duplications into a single file
rule make_duplications_bed_file:
    input:
        cuteSV = "{sample}/structural_variants/cuteSV/{sample}.cuteSV.bedpe",
        sniffles = "{sample}/structural_variants/sniffles/{sample}.sniffles.bedpe",
        svim = "{sample}/structural_variants/svim/{sample}.svim.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.duplications.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        grep_DUP="-h 'DUP'",
        select_columns="-f1,2,6,7-9",
        sort="-k1,1 -k2,2n -k3,3n",
        merge_length="100",
        columns="-c 4,5,6,6",
        output_type="-o distinct,mean,mean,count"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.make_duplications_bed_file.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.make_duplications_bed_file.txt"
    threads: 1
    shell:
        "grep {params.grep_DUP} {input.cuteSV} {input.sniffles} {input.svim} 2> {log} | \
        cut {params.select_columns} 2>> {log} | sort {params.sort} 2>> {log} | \
        {conda_dir}/bedtools merge -i stdin -d {params.merge_length} {params.columns} \
        {params.output_type} > {output} 2>> {log}"

# Merge the inversions into a single file
rule make_INV_bed_file:
    input:
        cuteSV = "{sample}/structural_variants/cuteSV/{sample}.cuteSV.bedpe",
        sniffles = "{sample}/structural_variants/sniffles/{sample}.sniffles.bedpe",
        svim = "{sample}/structural_variants/svim/{sample}.svim.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.inversions.bed")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        grep_INV="-h 'INV'",
        select_columns="-f1,2,6,7-9",
        sort="-k1,1 -k2,2n -k3,3n",
        merge_length="100",
        columns="-c 4,5,6,6",
        output_type="-o distinct,mean,mean,count"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.make_inversions_bed_file.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.make_inversions_bed_file.txt"
    threads: 1
    shell:
        "grep {params.grep_INV} {input.cuteSV} {input.sniffles} {input.svim} 2> {log} | \
        cut {params.select_columns} 2>> {log} | sort {params.sort} 2>> {log} | \
        {conda_dir}/bedtools merge -i stdin -d {params.merge_length} {params.columns} \
        {params.output_type} > {output} 2>> {log}"

# Merge the translocations into a single file        
rule make_translocations_bedpe_file:
    input:
        cuteSV = "{sample}/structural_variants/cuteSV/{sample}.cuteSV.bedpe",
        sniffles = "{sample}/structural_variants/sniffles/{sample}.sniffles.bedpe",
        svim = "{sample}/structural_variants/svim/{sample}.svim.bedpe"
    output:
        protected("{sample}/structural_variants/{sample}.translocations.bedpe")
    params:
        memory_per_thread="8G",
        run_time="0:1:0:0",
        awk_BND_TRA="'$7==\"BND\" || $7==\"TRA\"'",
        fix_columns="'{{if($1 > $4) {{print $4,$5,$6,$1,$2,$3,$7,$8,$9}} else {{print $0}}}}' OFS='\t'",
        select_columns="-f1-7,9",
        sort="-k1,1 -k2,2n -k3,3n"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.make_translocations_bed_file.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.make_translocations_bed_file.txt"
    threads: 1
    shell:
        "awk {params.awk_BND_TRA} {input.cuteSV} {input.sniffles} {input.svim} 2> {log} | \
        awk {params.fix_columns} 2>> {log} | cut {params.select_columns} 2>> {log} | \
        sort {params.sort} > {output} 2>> {log}"

# Merge split alignments
rule merge_split_alignments:
    input:
        bam = "{sample}/mapped/{sample}.phased.bam",
        bai = "{sample}/mapped/{sample}.phased.bam.bai"
    output:
        protected("{sample}/structural_variants/{sample}.merge_split_alignments.bed")
    params:
        memory_per_thread="8G",
        run_time="1:0:0:0"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.merge_split_alignments.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.merge_split_alignments.txt"
    threads: 1
    shell:
        "{conda_dir}/python {scripts_dir}/merge_split_alignments.py \
        --read-to-reference-bam {input.bam} --output-bedpe {output} &> {log}"

# Make a bed file of regions from the split alignments
rule split_alignments:
    input:
        "{sample}/structural_variants/{sample}.merge_split_alignments.bed"
    output:
        protected("{sample}/structural_variants/{sample}.split_alignments.bed")
    params:
        memory_per_thread="4G",
        run_time="0:1:0:0",
        sort="-k1,1 -k2,2n -k3,3n",
        merge_length="100",
        columns="5,5",
        output_type="mean,count"
    log:
        "{sample}/analysis/logs/structural_variants/{sample}.split_alignments.log"
    benchmark:
        "{sample}/analysis/benchmarks/structural_variants/{sample}.split_alignments.txt"
    threads: 1
    shell:
        "sort {params.sort} {input} 2> {log} | {conda_dir}/bedtools merge -i stdin \
        -d {params.merge_length} -c {params.columns} -o {params.output_type} > {output} 2>> {log}"
