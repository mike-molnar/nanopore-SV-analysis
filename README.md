# Nanopore-SV-Analysis
Structural variant filtering and analysis of Nanopore human WGS data.

![nanopore-SV-analysis](https://user-images.githubusercontent.com/39533525/165219632-9dd89a98-53dd-4abd-9701-df7fe990a21a.png)

## Installation instructions

Download the latest code from GitHub:

```
git clone https://github.com/mike-molnar/nanopore-SV-analysis.git
```

Before running the workflow, you will need to download the reference genome. I have not included the download as part of the workflow because it is designed to run on a cluster that may not have internet access.  You can use a local copy of GRCh38 if you have one, but the chromosomes must be named `chr1, chr2, ...` , and the reference can only contain the autosomes and sex chromosomes. To download the reference genome and index it, change to the reference directory of the workflow and run the `download_reference.sh` script:

```
cd /path/to/nanopore-SV-analysis/reference
chmod u+x download_reference.sh
./download_reference.sh
```

## Dependencies

There are many dependencies so it is best to create a new Conda environment using the provided YAML file:

```
conda env create -n nanopore-SV-analysis -f nanopore-SV-analysis.yml
conda activate nanopore-SV-analysis
```

Below is a list of the Conda dependencies:
- BCFtools
- bedtools
- karyoploteR
- cuteSV
- Longshot
- NanoFilt v2.8.0
- NanoPlot v1.20.0
- pybedtools
- PySAM
- PyVCF
- seaborn v0.10.0
- Snakemake
- Sniffles2
- SURVIVOR
- SVIM
- WhatsHap
- Winnowmap2

## To run on a grid engine

Copy the `Snakefile` and `config.yaml` files to the directory that you want to run the workflow.  Modify the information in `config.yaml` for your sample names and FASTQ locations. There are a few different grid engines, so the exact format to run the workflow may be different for your particular grid engine:

```
snakemake --jobs 500 --rerun-incomplete --keep-going --latency-wait 60 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -q queue_name -P project_name -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_rt={params.run_time} -b y"
```

You will have to replace `queue_name` and `project_name` with the necessary values to run on your cluster.
