# MethylSeq pipeline of TOBI lab

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.4-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)

## Introduction

> **Note:** This pipeline has been modified using the original nf-core/methylseq pipeline (version 2.6.0) to include the optical removal duplicates and some other modules in the Bismark subworkflow. You can find the original pipeline at the following [nfcore repository](https://nf-co.re/methylseq/2.6.0).

**MethylSeq** is a bioinformatics analysis pipeline used for Methylation (Bisulfite) sequencing data. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker / Singularity containers making installation trivial and results highly reproducible.

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/methylseq/results).

## Pipeline Summary

The pipeline allows you to choose between running either [Bismark](https://github.com/FelixKrueger/Bismark) or [bwa-meth](https://github.com/brentp/bwa-meth) / [MethylDackel](https://github.com/dpryan79/methyldackel).
Choose between workflows by using `--aligner bismark` (default, uses bowtie2 for alignment), `--aligner bismark_hisat` or `--aligner bwameth`.

| Step                                         | Bismark workflow      | bwa-meth workflow     |
| -------------------------------------------- | --------------------- | --------------------- |
| Generate Reference Genome Index _(optional)_ | Bismark               | bwa-meth              |
| Merge re-sequenced FastQ files               | cat                   | cat                   |
| Raw data QC                                  | FastQC                | FastQC                |
| Adapter sequence trimming                    | Trim Galore!          | Trim Galore!          |
| Align Reads                                  | Bismark               | bwa-meth              |
| Filter Non Conversion                        | Bismark               | -                     |
| Deduplicate Alignments                       | Bismark               | Picard MarkDuplicates |
| Removal of optical duplicates                | Picard MarkDuplicates | -                     |
| Extract methylation calls                    | Bismark               | MethylDackel          |
| Sample report                                | Bismark               | -                     |
| Summary Report                               | Bismark               | -                     |
| Alignment QC                                 | Qualimap              | Qualimap              |
| Sample complexity _(optional)_               | Preseq                | Preseq                |
| CpGs-level saturation _(optional)_           | methurator            | methurator            |
| Project Report                               | MultiQC               | MultiQC               |

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
SRR389222_sub1,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub1.fastq.gz
SRR389222_sub2,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub2.fastq.gz
SRR389222_sub3,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/SRR389222_sub3.fastq.gz
Ecoli_10K_methylated,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R1.fastq.gz,https://github.com/nf-core/test-datasets/raw/methylseq/testdata/Ecoli_10K_methylated_R2.fastq.gz
```

Each row represents a fastq file (single-end) or a pair of fastq files (paired end).

Now, you can run the pipeline using:

```bash
nextflow run main.nf --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/methylseq/2.6.0/docs/usage), the [parameter documentation](https://nf-co.re/methylseq/2.6.0/parameters) and the section below for the parameters that have been added to the modifed version of the pipeline.

## Parameters included in the modified version

The original [nf-core/methylseq pipeline](https://nf-co.re/methylseq/2.6.0/) has been modified to include modules that are important when processing RRBS/cfRRBS samples. The pipeline has been modified by @edogiuili and is currently maintained by Edoardo Giuili ([@edogiuili](https://github.com/edogiuili)) and Sofie Van de Velde ([@sofvdvel](https://github.com/sofvdvel)).

### Removal of optical duplicates

#### `--remove_optic_duplicates`

If specified, it removes the optical duplicates. It can be used together with `--sequencer` (see below).

#### `--sequencer`

If the flag `--remove_optic_duplicates` has been specified, `--sequencer` will be by default set to **NovaSeq**. Alternatively, you can specify **HiSeq** or **NextSeq**. This will change the OpticalDupsPixelDistance within PicardMarkDuplicates step (**NovaSeq**: 12000, **HiSeq**: 2500, **NextSeq**: 100).

### Filter non conversion (Bismark) module

#### `--filter_non_conversion`

If specified, it filters out all those reads that have a methylation value >= than a preset threshold in non-CG context where you expect a very low methylation level (<5% usually). For more information please consult the [Bismark usage](https://felixkrueger.github.io/Bismark/bismark/filter_nonconverted_reads/).

#### `--minimum_count`

Minimum number of methylation sites for a read to be filtered out (def. 3).

#### `--percentage_cutoff`

Minimum methylation percentage for a read to be filtered out (def. 90%).

### CpG-level sequencing saturation using methurator

#### `--run_methurator`

If specified (by default: true), it will run [methurator](https://vibtobilab.github.io/methurator/latest/) to compute a CpG-level sequencing saturation analysis.

#### `--t_max`

Maximum extrapolation factor (t) value (def. 10).

#### `--minimum_coverage`

Minimum number of counts to define a new CpG (def. 1). More than 1 value can be specified (e.g. '1,3,5').

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/methylseq/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/methylseq/output).

## Credits

These scripts were originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.

- Main authors of original nf-core/methylseq pipeline:
  - Phil Ewels ([@ewels](https://github.com/ewels/))
- Maintainers:
  - Felix Krueger ([@FelixKrueger](https://github.com/FelixKrueger))
  - Sateesh Peri ([@Sateesh_Peri](https://github.com/sateeshperi))
  - Edmund Miller ([@EMiller88](https://github.com/emiller88))
- Contributors:
  - Rickard Hammarén ([@Hammarn](https://github.com/Hammarn/))
  - Alexander Peltzer ([@apeltzer](https://github.com/apeltzer/))
  - Patrick Hüther ([@phue](https://github.com/phue/))

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
