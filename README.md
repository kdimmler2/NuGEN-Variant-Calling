# Variant Calling in Targeted Sequencing Data

This is a Snakemake workflow based on the GATK best practices, adapted for targeted sequencing data. After creating GVCFs, it subsets each GVCF to contain only target regions (100 bp +/- each target locus).

This repository is under construction. A step for Ampliconclip needs to be added, and the variant filtering parameters need to be edited based on new findings in the literature for targeted sequencing data.
