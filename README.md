# nanopore-SV-analysis
Structural variant filtering and anlaysis of Nanopore human WGS data.

## Installation instructions

Download the latest code from GitHub:

```
git clone https://github.com/mike-molnar/nanopore-SV-analysis.git
```

Before running the workflow, you will need to download the reference genome. I have not included the download as part of the workflow because it is designed to run on a cluster that may not have internet access.  You can use a local copy of GRCh38 if you have one, but the chromosomes must be named `chr1, chr2, ...` , and the reference can only contain the autosomes and sex chromosomes. To download the reference genome and index it, change to the reference directory of the workflow and run the script:

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

## To run on a grid engine

Copy the `Snakefile` and `config.yaml` files to the directory that you want to run the workflow. Modify the input in `config.yaml` for the FASTQs and sample names for the samples. There are a few different grid engines, so the exact format to run the workflow may be different for your particular grid engine:

```
snakemake --jobs 500 --rerun-incomplete --keep-going --latency-wait 60 --cluster "qsub -cwd -V -o snakemake.output.log -e snakemake.error.log -q queue_name -P project_name -pe smp {threads} -l h_vmem={params.memory_per_thread} -l h_rt={params.run_time} -b y" all_but_assembly
```

You will have to replace `queue_name` and `project_name` with the necessary values to run on your grid.
