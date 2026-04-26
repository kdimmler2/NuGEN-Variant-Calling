# NuGEN Variant Calling Prototype

This repository contains an independently written prototype Snakemake workflow for raw sequencing data processing and variant calling. The workflow was developed while I was building and testing approaches for processing sequencing data (custom targeted sequencing data, in particular) in my thesis project.

## Overview

The main workflow is defined in:

```text
variant_calling.smk
```

The workflow includes steps for:

- detecting paired-end and single-end FASTQ inputs
- converting FASTQ files to unmapped BAMs
- marking adapters
- aligning reads to the reference genome
- sorting and merging aligned BAM files
- recalibrating base quality scores
- generating per-sample GVCFs
- joint genotyping
- variant filtering
- variant annotation with SnpEff
- generating QC outputs

## Notes on workflow context

This repository represents an initial prototype workflow that I wrote myself. It was not the final production workflow used for my thesis analysis.

For the final analysis, I later modified a more robust lab pipeline to better fit the targeted sequencing assay used in the project. Those later modifications included:

- handling both paired-end and single-end data
- adding amplicon clipping for targeted sequencing data
- removing the duplicate-marking step, as appropriate for targeted sequencing data
- Changing the behavior of final variant filtration, as appropriate for targeted sequencing data

This prototype is still useful as a code example because it shows how I structured a raw sequencing data processing workflow in Snakemake and handled mixed input types.

## Input data

The workflow expects FASTQ files organized into separate paired-end and single-end input directories defined in the configuration file.

Input paths and reference files are provided through:

```text
config.yaml
```

Expected input categories include:

- paired-end FASTQ files
- single-end FASTQ files
- reference genome files
- known sites for BQSR
- A list of target loci for targeted sequencing data
- Defined regions around each target locus for qualimap
- SnpEff database files

## Single-end and paired-end handling

The workflow identifies whether each sample is paired-end or single-end based on the available FASTQ files. This allows the same workflow to process mixed sequencing layouts.

Paired-end samples are identified from R2 FASTQ files so that each sample/lane is only added once. Single-end samples are identified from the single-end FASTQ directory, with duplicate files skipped.

## Main workflow stages

1. **Input discovery**  
   FASTQ files are scanned and sample layout information is stored for downstream rules.

2. **uBAM creation and adapter marking**  
   FASTQ files are converted to unmapped BAM format and adapters are marked.

3. **Alignment and BAM processing**  
   Reads are aligned to the reference genome, sorted, merged, and prepared for variant calling.

4. **Base quality score recalibration**  
   BQSR is performed using known variant sites.

5. **GVCF generation**  
   Per-sample GVCFs are generated with GATK HaplotypeCaller.

6. **Joint genotyping and filtering**  
   GVCFs are combined for cohort-level genotyping, followed by variant filtering.

7. **Annotation and QC**  
   Filtered variants are annotated with SnpEff, and QC outputs are generated throughout the workflow.

## Software requirements

This prototype workflow depends on common bioinformatics tools, including:

- Snakemake
- GATK
- BWA
- Samtools
- Qualimap
- SnpEff

A formal Conda environment file has not yet been added to this prototype repository.

## How to run

Example command:

```bash
snakemake -s nugen_pipeline.smk \
    --profile slurm.nugen \
    --keep-going \
    --rerun-incomplete \
    --rerun-triggers mtime \
```

Alternatively, on an HPC system, this workflow can be run from an interactive session or submitted through a scheduler script, depending on cluster configuration.

    sbatch nugen_pipeline.slurm

## Limitations

This repository is a prototype and is not intended to represent the final production pipeline used for the thesis analysis. In particular, this workflow does not include some targeted-sequencing-specific steps used later in the project, such as amplicon clipping, altered filtering methods, and removal of duplicate-marking.

Future improvements would include:

- adding a formal environment file
- adding a small test dataset
- documenting expected file naming conventions more fully
- adapting the workflow for targeted sequencing-specific processing
- Adding additional QC steps, such as fastqc
- Adding the ability to do variant filtration through VQSR or hard filtering
- Adding intervals for faster parallel processing of sequence data
