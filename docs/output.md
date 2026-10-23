# MethylSeq: Output

## Introduction

This page describes the files produced by the pipeline. The directories below are created in the output directory (`--outdir`) and all paths are relative to it. Most plots shown here are taken from the MultiQC report, which summarises the results of all samples at the end of the pipeline.

The files depend on the workflow that is run: `<aligner>` in the paths below is the value of `--aligner` (`bismark`, `bismark_hisat` or `bwameth`).

## Pipeline overview

- [Merged FastQ files](#merged-fastq-files): concatenation of the re-sequenced samples
- [FastQC](#fastqc): raw read QC
- [Trim Galore!](#trim-galore): adapter and quality trimming
- [Reference genome](#reference-genome): genome indices
- [Alignment](#alignment): alignment of the reads to the reference genome
- [Filtering of reads with an incomplete conversion](#filtering-of-reads-with-an-incomplete-conversion): Bismark workflow only
- [Removal of optical duplicates](#removal-of-optical-duplicates): Bismark workflow only
- [Deduplication](#deduplication): removal or marking of the duplicated reads
- [Sorted and indexed alignments](#sorted-and-indexed-alignments): final BAM files
- [Methylation extraction](#methylation-extraction): cytosine methylation calls
- [Bismark reports](#bismark-reports): sample and summary reports
- [Qualimap](#qualimap): alignment QC
- [Preseq](#preseq): library complexity
- [methurator](#methurator): CpG-level sequencing saturation
- [MultiQC](#multiqc): report of the QC of all samples
- [Pipeline information](#pipeline-information): reports of the execution of the pipeline

## Merged FastQ files

<details markdown="1">
<summary>Output files</summary>

- `cat/`
  - `<sample>.merged.fastq.gz`: concatenated FastQ files of a sample sequenced more than once (see [multiple runs of the same sample](usage.md#multiple-runs-of-the-same-sample)).

</details>

## FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report with the quality metrics of the raw reads.
  - `zips/*_fastqc.zip`: archive with the FastQC report, a tab-delimited data file and the plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics of the sequenced reads: the distribution of the quality scores, the per base sequence content (%A/T/G/C), the adapter contamination and the overrepresented sequences. See the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/) for more information.

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

!!! note

    These plots are generated from the raw reads, which may contain adapter sequences and low-quality regions. The FastQC reports of the trimmed reads are in `trimgalore/fastqc/`.

## Trim Galore!

<details markdown="1">
<summary>Output files</summary>

- `trimgalore/`
  - `*_val_1.fq.gz`, `*_val_2.fq.gz` (paired-end) or `*_trimmed.fq.gz` (single-end): trimmed reads. Only saved with `--save_trimmed`.
  - `logs/*_trimming_report.txt`: trimming report, with the parameters used.
  - `fastqc/*_fastqc.html`, `fastqc/zips/*_fastqc.zip`: FastQC reports of the trimmed reads.

</details>

[Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) removes the adapter sequences and the low-quality regions of the reads with [Cutadapt](https://github.com/marcelm/cutadapt), then runs FastQC on the trimmed reads. MultiQC reports the percentage of bases removed in the _General Statistics_ table, together with a plot of where the reads were trimmed.

## Reference genome

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/reference_genome/`
  - Bismark or bwa-meth index and FASTA index (`.fai`) built by the pipeline. Only saved with `--save_reference`.

</details>

## Alignment

Bismark and bwa-meth convert _in silico_ all the cytosines of the reads to thymines and align them to a three-letter reference genome, which avoids a methylation-specific alignment bias.

<details markdown="1">
<summary>Output files: Bismark workflow</summary>

- `<aligner>/alignments/`
  - `*.bam`: aligned reads. Only saved with `--save_align_intermeds`, or with `--rrbs` (or `--skip_deduplication`) when the optical duplicates are not removed.
  - `logs/*_report.txt`: summary statistics of the alignment.
  - `unmapped/*.fq.gz`: unmapped reads. Only saved with `--unmapped`.

</details>

<details markdown="1">
<summary>Output files: bwa-meth workflow</summary>

- `<aligner>/alignments/`
  - `*.bam`, `*.sorted.bam`, `*.sorted.bam.bai`: aligned, sorted and indexed reads. Only saved with `--save_align_intermeds`; the sorted reads are also saved with `--rrbs`.
  - `samtools_stats/*.flagstat`: number of reads aligned in each way.
  - `samtools_stats/*.stats`: statistics of the aligned reads.

</details>

## Filtering of reads with an incomplete conversion

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/filter_non_converted/`
  - `*.non-conversion_filtering.txt`: number of reads filtered out.
  - `*.nonCG_filtered.bam`, `*.nonCG_removed_seqs.bam`: alignments without the filtered reads, and the filtered reads. Only saved with `--save_align_intermeds`.

</details>

Only run in the Bismark workflow with `--filter_non_conversion` (see [filtering of reads with an incomplete conversion](usage.md#filtering-of-reads-with-an-incomplete-conversion)).

## Removal of optical duplicates

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/picard_markoptical_deduplicated/`
  - `picard_metrics/*.MarkDuplicates.metrics.txt`: number of optical duplicates found by Picard MarkDuplicates.
  - `*.bam`: alignments without the optical duplicates. Only saved with `--save_align_intermeds`.
- `<aligner>/sorted_queryname_reads/`
  - `*.bam`: alignments sorted by read name. Only saved with `--save_align_intermeds`.

</details>

Only run in the Bismark workflow, unless `--remove_optic_duplicates false` is given (see [removal of optical duplicates](usage.md#removal-of-optical-duplicates)).

## Deduplication

Duplicated reads, with the same alignment positions, are technical copies of the same DNA fragment. This step is skipped with `--rrbs` or `--skip_deduplication`.

<details markdown="1">
<summary>Output files: Bismark workflow</summary>

- `<aligner>/deduplicated/`
  - `logs/*.deduplication_report.txt`: summary statistics of the deduplication.
  - `*.deduplicated.bam`: alignments without the duplicated reads. Only saved with `--save_align_intermeds`.

</details>

<details markdown="1">
<summary>Output files: bwa-meth workflow</summary>

- `<aligner>/deduplicated/`
  - `*.markdup.sorted.bam`, `*.markdup.sorted.bam.bai`: alignments with the duplicated reads marked, and their index.
  - `picard_metrics/*.metrics.txt`: summary statistics of the deduplication.

</details>

!!! note

    In the bwa-meth workflow, Picard MarkDuplicates only marks the duplicated reads, without removing them from the BAM file.

## Sorted and indexed alignments

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/sorted_indexed_reads/`
  - `*.sorted.bam`, `*.sorted.bam.bai`: final alignments of the Bismark workflow, sorted by coordinate and indexed. These are the alignments used for the methylation extraction, Qualimap and methurator.

</details>

## Methylation extraction

The methylation extractor takes the aligned reads and generates files with the methylation calls of the cytosines. The files may differ slightly with `--comprehensive` or `--non_directional`.

The file names contain the strand the reads aligned to:

- `OT`: original top strand
- `OB`: original bottom strand
- `CTOT`: complementary to the original top strand
- `CTOB`: complementary to the original bottom strand

<details markdown="1">
<summary>Output files: Bismark workflow</summary>

- `<aligner>/methylation_calls/`
  - `methylation_calls/*_context_*.txt.gz`: individual methylation calls, in one file per cytosine context and strand.
  - `methylation_coverage/*.bismark.cov.gz`: methylation level and coverage of each cytosine.
  - `bedGraph/*.bedGraph.gz`: methylation levels in [bedGraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) format, with 0-based start and 1-based end coordinates.
  - `mbias/*.M-bias.txt`: methylation bias along the reads. See the [Bismark documentation](https://felixkrueger.github.io/Bismark/bismark/methylation_extraction/#m-bias-plot) for more information.
  - `splitting_report/*_splitting_report.txt`: summary statistics of the methylation extraction.
- `bismark/coverage2cytosine/`: genome-wide cytosine reports. Only created with `--cytosine_report` or `--nomeseq`.
  - `reports/*_report.txt.gz`: methylation of every cytosine of the genome.
  - `summaries/*cytosine_context_summary.txt`: summary of the cytosine report per context.
  - `coverage/*.cov.gz`: coverage files (NOMe-seq only).

</details>

!!! note

    `CTOT` and `CTOB` are only aligned with `--non_directional`.

<details markdown="1">
<summary>Output files: bwa-meth workflow</summary>

- `methyldackel/`
  - `*.bedGraph`: methylation levels of the CpGs in [bedGraph](http://genome.ucsc.edu/goldenPath/help/bedgraph.html) format.
  - `*.methylKit`: methylation levels in [methylKit](https://bioconductor.org/packages/release/bioc/html/methylKit.html) format. Only created with `--methyl_kit`.
  - `mbias/*.mbias.txt`: methylation bias along the reads.

</details>

## Bismark reports

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/reports/`
  - `*_report.html`: report of a sample, with the alignment, deduplication and methylation extraction statistics.
- `<aligner>/summary/`
  - `bismark_summary_report.html`, `bismark_summary_report.txt`: summary report of all samples.

</details>

## Qualimap

<details markdown="1">
<summary>Output files</summary>

- `qualimap/<sample>/`
  - `qualimapReport.html`: Qualimap HTML report.
  - `genome_results.txt`, `raw_data_qualimapReport/*.txt`: statistics in text format.

</details>

[Qualimap BamQC](http://qualimap.conesalab.org/doc_html/analysis.html#bam-qc) computes statistics of the aligned reads, e.g. the insert sizes and the coverage. Its statistics are not specific to bisulfite data. A file of target regions can be given with `--bamqc_regions_file`.

## Preseq

<details markdown="1">
<summary>Output files</summary>

- `<aligner>/preseq/`
  - `*.lc_extrap.txt`: values of the complexity curve, plotted in the MultiQC report.
  - `log/*.log`: log of Preseq.

</details>

[Preseq](http://smithlabresearch.org/software/preseq/) estimates the complexity of a library, i.e. how many new unique reads would be obtained with more sequencing. A flat curve means that the library is saturated and that more sequencing would add few unique reads. The dashed line shows a perfectly complex library, where all reads are unique.

These are predictions only. The MultiQC plot can show extreme sequencing depths on the x axis: click and drag from the left of the plot to zoom in on realistic values.

Preseq is skipped with `--rrbs` or `--skip_preseq`.

## methurator

<details markdown="1">
<summary>Output files</summary>

- `methurator/`
  - `summaries/<sample>/methurator_summary.yml`: metadata, model parameters and extrapolation results of a sample.
  - `plots/*.html`: interactive saturation plot of each sample.

</details>

[methurator](https://vibtobilab.github.io/methurator/latest/) estimates the CpG-level sequencing saturation of each sample, i.e. how many new CpGs would be covered with more sequencing. It is skipped with `--run_methurator false`.

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/<aligner>/`
  - `multiqc_report.html`: standalone HTML report, to open in a web browser.
  - `multiqc_data/`: statistics parsed from the tools of the pipeline.
  - `multiqc_plots/`: static images of the plots of the report.

</details>

[MultiQC](http://multiqc.info) summarises the QC of all samples in a single HTML report, together with the versions of the tools used. Further statistics are in `multiqc_data/`. See <http://multiqc.info> for how to use the MultiQC reports.

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.html`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `methylseq_software_mqc_versions.yml`. The `pipeline_report*` files are only created with `--email` or `--email_on_fail`.
  - Parameters used by the run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) reports the execution of the pipeline, e.g. the commands, run times and resource usage of the processes, which helps troubleshoot failed runs.
