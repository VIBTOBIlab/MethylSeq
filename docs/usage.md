# MethylSeq: Usage

## Introduction

MethylSeq is two pipelines in one, chosen with `--aligner`:

- `bismark` (default): [Bismark](https://github.com/FelixKrueger/Bismark) with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as the aligner.
- `bismark_hisat`: Bismark with [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) as the aligner. HISAT2 also supports splice-aware alignment, e.g. for [SLAM-seq](https://science.sciencemag.org/content/360/6390/800) experiments, with a file of known splice sites given with `--known_splices`.
- `bwameth`: [bwa-meth](https://github.com/brentp/bwa-meth) for the alignment and [MethylDackel](https://github.com/dpryan79/methyldackel) for the methylation calls.

The modules added to the original [nf-core/methylseq](https://nf-co.re/methylseq/2.6.0) pipeline are described in [Modules added to nf-core/methylseq](#modules-added-to-nf-coremethylseq). All the parameters are listed on the [parameters](parameters.md) page.

## Samplesheet

The samples are given in a comma-separated samplesheet with a header row and 3 columns:

```bash
--input '[path to samplesheet file]'
```

| Column | Description |
| --- | --- |
| `sample` | Sample name, without spaces. |
| `fastq_1` | Path to the FastQ file of read 1. The file has to be gzipped and have the extension `.fastq.gz` or `.fq.gz`. |
| `fastq_2` | Path to the FastQ file of read 2, for paired-end data. Leave it empty for single-end data. |

The pipeline detects whether a sample is single-end or paired-end from the `fastq_2` column, so a samplesheet can mix both:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_S1,CONTROL_S1_R1.fastq.gz,CONTROL_S1_R2.fastq.gz
TREATMENT_S2,TREATMENT_S2_R1.fastq.gz,
```

An [example samplesheet](https://github.com/VIBTOBIlab/MethylSeq/blob/main/assets/samplesheet.csv) is provided with the pipeline.

### Multiple runs of the same sample

The FastQ files of a sample sequenced more than once (e.g. on several lanes) are concatenated before the analysis. The pipeline groups the rows by the sample name **without its last `_` suffix**, so in the samplesheet below the three rows are merged into one sample called `CONTROL`:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_L002,CONTROL_L002_R1.fastq.gz,CONTROL_L002_R2.fastq.gz
CONTROL_L003,CONTROL_L003_R1.fastq.gz,CONTROL_L003_R2.fastq.gz
CONTROL_L004,CONTROL_L004_R1.fastq.gz,CONTROL_L004_R2.fastq.gz
```

!!! warning "Sample names are shortened"

    The last `_` suffix is removed from every sample name, also when a sample has a single row: `CONTROL_S1` becomes `CONTROL`. Samples whose names only differ after the last `_` (e.g. `PATIENT_A` and `PATIENT_B`) are merged into one. Give each sample a unique name before its last `_`, e.g. `PATIENTA_S1` and `PATIENTB_S2`.

## Reference genome

The reference genome is given by its [iGenomes](https://nf-co.re/docs/usage/reference_genomes) name (`--genome`, e.g. `GRCh38`) or as a FASTA file (`--fasta`). The pipeline builds the Bismark or bwa-meth index and the FASTA index when they are not given. To reuse them in later runs, save them with `--save_reference` and give them with `--bismark_index` (or `--bwa_meth_index`) and `--fasta_index`.

## Library types

The library preparation determines how the reads are trimmed and aligned. The presets below set the right options of Trim Galore! and Bismark:

| Parameter | Library |
| --- | --- |
| `--rrbs` | Reduced representation bisulfite sequencing (MspI digested material), e.g. RRBS and cfRRBS. |
| `--pbat` | Post-bisulfite adapter tagging (PBAT) libraries. |
| `--em_seq` | Enzymatic methyl-seq (EM-seq) libraries. |
| `--single_cell` | Single-cell bisulfite libraries. |
| `--accel`, `--cegx`, `--epignome`, `--zymo` | Trimming presets for the Accel, CEGX, Epignome and Zymo kits. |

!!! note "RRBS"

    With `--rrbs`, Trim Galore! removes the filled-in cytosines at the MspI sites, and the deduplication and Preseq are skipped, since the reads of an RRBS library start at the same positions. The optical duplicates are still removed (see [`--remove_optic_duplicates`](#--remove_optic_duplicates)).

## Modules added to nf-core/methylseq

The original nf-core/methylseq pipeline has been extended with modules that are important to process RRBS and cfRRBS samples.

### Removal of optical duplicates

Optical duplicates are copies of the same cluster on the flow cell, which are not removed by the deduplication of RRBS data. In the Bismark workflow, they are marked and removed with Picard MarkDuplicates before the deduplication.

#### `--remove_optic_duplicates`

Remove the optical duplicates (default: `true`). Use `--remove_optic_duplicates false` to keep them.

#### `--sequencer`

Sequencer used, which sets the maximum distance between optical duplicates (`--OPTICAL_DUPLICATE_PIXEL_DISTANCE`) of Picard MarkDuplicates:

| `--sequencer` | Pixel distance |
| --- | --- |
| `NovaSeq` (default) | 12000 |
| `HiSeq` | 2500 |
| `NextSeq` | 100 |

### Filtering of reads with an incomplete conversion

Reads with a high methylation in non-CG context (CHG and CHH), where methylation is expected to be low (usually < 5%), most likely come from an incomplete bisulfite conversion. In the Bismark workflow, they can be removed with [`filter_non_conversion`](https://felixkrueger.github.io/Bismark/bismark/filter_nonconverted_reads/) after the alignment.

#### `--filter_non_conversion`

Filter out the reads with an incomplete conversion (default: `false`).

#### `--minimum_count`

Minimum number of cytosines in non-CG context in a read before the `--percentage_cutoff` filter is applied (default: `3`).

#### `--percentage_cutoff`

Minimum methylation percentage in non-CG context for a read to be filtered out (default: `90`).

### CpG-level sequencing saturation

[methurator](https://vibtobilab.github.io/methurator/latest/) estimates how many new CpGs would be covered with more sequencing, from the deduplicated alignments of each sample. With `--rrbs`, it runs in RRBS mode.

#### `--run_methurator`

Run the sequencing saturation analysis (default: `true`). Use `--run_methurator false` to skip it.

#### `--t_max`

Maximum extrapolation factor `t`, i.e. up to how many times the current sequencing depth the saturation is estimated (default: `10`).

#### `--minimum_coverage`

Minimum number of reads covering a CpG to count it as covered (default: `1`). Several values can be given, e.g. `'1,3,5'`, to get one curve per value.

## Running the pipeline

The typical command to run the pipeline is:

```bash
nextflow run VIBTOBIlab/MethylSeq --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile docker
```

This runs the pipeline with the `docker` profile (see [`-profile`](#-profile)) and creates the following files in the working directory:

```bash
work                # Directory containing the Nextflow working files
<OUTDIR>            # Finished results in the specified location (defined with --outdir)
.nextflow.log       # Log file from Nextflow
# Other Nextflow hidden files, e.g. history of pipeline runs and old logs.
```

To reuse the same parameters in several runs, write them in a YAML or JSON file and give it with `-params-file`:

```bash
nextflow run VIBTOBIlab/MethylSeq -profile docker -params-file params.yaml
```

```yaml title="params.yaml"
input: ./samplesheet.csv
outdir: ./results
genome: GRCh38
rrbs: true
```

!!! warning

    Do not use `-c <file>` to give parameters. Configuration files given with `-c` must only be used to [tune the resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) of the processes, other infrastructural settings or the arguments of the modules.

### Updating the pipeline

Nextflow caches the pipeline code the first time it runs it, and keeps using the cached version afterwards. Update it to the latest version with:

```bash
nextflow pull VIBTOBIlab/MethylSeq
```

### Reproducibility

Give a version of the pipeline with `-r` (e.g. `-r v2.6.4`), so that the same code and software versions are used every time. The versions are listed on the [releases page](https://github.com/VIBTOBIlab/MethylSeq/releases), and the version used is reported at the bottom of the MultiQC report.

## Core Nextflow arguments

!!! note

    These options are part of Nextflow and use a _single_ hyphen, while the pipeline parameters use a double hyphen.

### `-profile`

Configuration profile to use. Profiles set the software packaging (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer or Conda) and can be combined, e.g. `-profile test,docker`. They are loaded in order, so later profiles overwrite earlier ones. Without `-profile`, all software has to be installed and available on the `PATH`, which is _not_ recommended.

The pipeline also loads the institutional profiles of [nf-core/configs](https://github.com/nf-core/configs), e.g. for HPC clusters.

- `test`: a small test dataset with all the parameters needed, to check the setup.
- `docker`, `singularity`, `podman`, `shifter`, `charliecloud`, `apptainer`: run the processes in containers with the corresponding tool.
- `conda`: run the processes in Conda environments. Use it only when containers cannot be used.

### `-resume`

Restart a run from where it stopped, reusing the results of the steps whose inputs did not change. Give a run name (e.g. `-resume [run-name]`) to resume a specific run; `nextflow log` lists the previous runs. See [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html) for how Nextflow decides what to reuse.

### `-c`

Path to a custom configuration file, e.g. to change the resources of a process. See the [nf-core documentation](https://nf-co.re/docs/usage/configuration) for more information.

### `-bg`

Run Nextflow in the background, so that the pipeline keeps running when you log out. Alternatively, use `screen` or `tmux`, or submit Nextflow itself as a job to your cluster scheduler.

## Custom configuration

### Resource requests

Each process has default CPUs, memory and time requests. Most processes that fail because they ran out of resources are retried with twice, then three times the original requests. To change the requests of a process, write them in a configuration file and give it with `-c`:

```groovy title="custom.config"
process {
    withName: 'NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN' {
        cpus   = 16
        memory = 64.GB
    }
}
```

See [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) on the nf-core website for more details.

### Custom tool arguments

The arguments of a tool can be changed with the `ext.args` of its process, e.g. to add arguments that are not exposed as pipeline parameters. See [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) on the nf-core website.

### Custom containers

To use another version of a tool, override the `container` (or `conda`) of its process in a configuration file. See [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) on the nf-core website.

## Nextflow memory requirements

The Java virtual machine of Nextflow can request a large amount of memory. Limit it by adding the following line to your environment, e.g. in `~/.bashrc`:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
