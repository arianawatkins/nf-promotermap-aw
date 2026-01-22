# scbirlab/nf-promotermap

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-promotermap/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-promotermap** maps Illumina sequences to bacterial genomes and calls peaks.

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

The pipeline carries out the following steps, given a [sample sheet (see below)](#sample-sheet):

1. Downloads reference genome and annotations from NCBI
2. Trims adapter sequences from Illumina reads using `cutadapt`
3. Aligns to reference genome using either `bowtie2` or `minimap2`
4. Plot gene start coverage with `deeptools`.
5. Call peaks across all samples with `MACS3`.
6. Annotate peaks with nearest genes.
7. Generate FASTA of peak sequences.
8. Calculate coverage of each peak for each bin.
9. Calculate coverage variance across bins.
10. Calculate per-base coverage within each peak for each bins and mean and variance across bins.

### Work in progress

- [ ] Identify elements associated with strength and variance.
- [ ] Identify common sequence motifs in those elements. 

### Other steps

1. Get FASTQ quality metrics with `fastqc`.
2. Calculate coverage and other with `samtools`.  
3. Compile the logs of processing steps into an HTML report with `multiqc`.

## Requirements

### Software

You need to have Nextflow and either Anaconda, Singularity, or Docker installed on your system.

#### Crick users (and other HPC users)

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow Singularity
```

#### Everyone else: installing Nextflow 

Otherwise, if it's your first time using Nextflow on your system and you have Conda installed, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet) and, optionally, 
a [`nextflow.config` file](#inputs) in the directory where you want the 
pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-promotermap -latest
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.3`, 
you can do so using

```bash 
nextflow run scbirlab/nf-promotermap -r v0.0.3
```

For help, use `nextflow run scbirlab/nf-promotermap --help`.

The first time you run the pipeline, the software dependencies 
in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed
- `fastq_dir`: path to where FASTQ files are stored
- `control_label`: the bin ID (from [sample sheet](#sample-sheet)) of background controls

The following parameters have default values which can be overridden if necessary.

- `inputs = "inputs"` : The folder containing your inputs.
- `outputs = "outputs"` : The folder to containing the pipeline outputs.
- `trim_qual = 5`: Minimum base-call quality for trimming.
- `min_length = 9`: Discard reads shorter than this number of bases after trimming.
- `mapper = "bowtie2"`: Alignment tool.

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "/path/to/sample-sheet.csv"
    inputs = "/path/to/inputs"
    fastq_dir = "/path/to/fastq"
    control_label = "U" // bin_id of your background control

    mapper = "minimap2"
}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-promotermap \
    --sample_sheet /path/to/sample-sheet.csv \
    --inputs /path/to/inputs \
    --fastq_dir /path/to/fastq \
    --control_label U \
    --mapper minimap2
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which FASTQ files belong to which sample.

The file must have a header with the column names below (in any order), and one line per sample to be processed. 
You can have additional columns eith extra information if you like.

- `expt_id`: Unique name of a peak-calling experiment. Peaks will be called across all samples with the same experiment ID.
- `sample_id`: Unique name of the sample within an experiment. FASTQ files under the same sample ID will be combined.
- `bin_id`:  Unique name of a bin within an experiment. Sample IDs under the same bin will be pooled before coverage analysis.
- `fastq_pattern`: Partial filename that matches at least both R1 and R2 FASTQ files for a sample in the `fastq_dir` ([defined above](#inputs)).
- `genome_accession`: The [NCBI assembly accession](https://www.ncbi.nlm.nih.gov/datasets/genome/) number for the genome for alignment and annotation. This number starts with "GCF_" or "GCA_".
- `adapter_read1_3prime`: [the 3' adapter on the forward read to trim](#cutadapt-format). The adapter itself and sequences _downstream_ will be removed.
- `adapter_read2_3prime`: [the 3' adapter on the reverse read to trim](#cutadapt-format). The adapter itself and sequences _downstream_ will be removed.
- `adapter_read1_5prime`: [the 5' adapter on the forward read to trim](#cutadapt-format). The adapter itself and sequences _upstream_ will be removed.
- `adapter_read2_5prime`: [the 5' adapter on the reverse read to trim](#cutadapt-format). The adapter itself and sequences _upstream_ will be removed.

Here is an example of the sample sheet:

| expt_id | sample_id   | bin_id | fastq_pattern | genome_accession | adapter_read1_3prime  | adapter_read2_3prime | adapter_read1_5prime | adapter_read2_5prime | 
| ------- | ----------- | ------ | ------------- | ---------------- | --------------------- | -------------------- | -------------------- | -------------------- |
| expt-01 | 01-Unsorted | U      | G5512A22_R    | GCF_904425475.1  | ATTAACCTCCTAATCGTGCGT | CTACCGCCTTGCTGCTGCGT | ACGCAGCAGCAAGGCGG    | ACGCACGATTAGGA       |
| expt-01 | 01-Red1     | Red1   | G5512A23_R    | GCF_904425475.1  | ATTAACCTCCTAATCGTGCGT | CTACCGCCTTGCTGCTGCGT | ACGCAGCAGCAAGGCGG    | ACGCACGATTAGGA       |


#### Cutadapt format

Read more [here](https://cutadapt.readthedocs.io/en/stable/guide.html#specifying-adapter-sequences).

### Example inputs

You can find some examples in the `test` directory of this repository.

## Outputs

Outputs are saved in the directory specified by `--outputs` (`outputs` by default). 
They are organised into these directories:

- `bigwig`: Coverage bigwig files
- `coverage`: Coverage of peaks per bin
- `genome`: Reference genomes and annotations
- `mapped`: BAM files of mapped Illumina reads
- `multiqc`: HTML reports from the outputs of intermediate steps
- `peaks`: peak calls
- `samtools`: Coverage and other metrics
- `trimmed`: Trimming logs and FASTQ files.

## Issues, problems, suggestions

If you run into problems not covered here, add to the 
[issue tracker](https://www.github.com/scbirlab/nf-promotermap/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
- [deeptools](https://deeptools.readthedocs.io/en/develop/index.html)
- [fastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MACS3](https://macs3-project.github.io/MACS/index.html)
- [minimap2](https://lh3.github.io/minimap2/minimap2.html)
- [multiQC](https://multiqc.info/)
- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [samtools](http://www.htslib.org/doc/samtools.html)
