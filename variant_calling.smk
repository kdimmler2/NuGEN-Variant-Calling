import os
import gzip

localrules: gvcf_list,
            vcf_list

#Create a dictionary to identify FASTQs that are paired-end and single-end
layouts = {}

# Identify paired-end FASTQ files.
# Only R2 files are used here so each paired-end sample/lane is added once.
for root, dirs, files in os.walk(config['paired_fastqs']):
    for file in files:
        if file.endswith('.fastq.gz') and 'R2' in file:
            Mnum, sample, lane, read, end = file.split('_')
            test_list.append(f'{Mnum}')
            layouts[f'{Mnum}_{sample}_{lane}'] = 'paired'

# Identify single-end FASTQ files.
# Duplicate files are skipped, and each remaining sample is labeled as single-end.
for root, dirs, files in os.walk(config['single_fastqs']):
    for file in files:
        if file.endswith('.fastq.gz'):
            if 'Dup' in file:
                pass
            else:
                Mnum, sample, read, end = file.split('_')
                layouts[f'{Mnum}_{sample}'] = 'single'

include: 'src/utils.py'

PE_FASTQS,READ, = glob_wildcards(os.path.join('paired_fastqs','{pe_fastq}_{read}_001.fastq.gz'))
SE_FASTQS, = glob_wildcards(os.path.join('single_fastqs','{se_fastq}.fastq.gz'))

# Target outputs for the full workflow.
# Some intermediate outputs are included intentionally because the workflow
# branches across QC, per-sample processing, joint genotyping, filtering,
# and annotation steps.
rule all:
    input:
        expand('Data_Preprocessing/uBAMs/{sample}/{sample}.fastqtosam.bam', sample = layouts.keys()),
        expand('Data_Preprocessing/Marked_uBAMs/{sample}/{sample}.markilluminaadapters.bam', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/{sample}.aligned.bam', sample = layouts.keys()),
        expand('Data_Preprocessing/Marked_BAMs/{sample}/{sample}.MarkDuplicates.bam', sample = layouts.keys()),
        expand('Data_Preprocessing/BQSR/BaseRecalibrator/{sample}/{sample}.table', sample = layouts.keys()),
        expand('Data_Preprocessing/BQSR/ApplyBQSR/{sample}/{sample}.recalibrated.bam', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/qualimap/report.pdf', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/qualimap/report.html', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/qualimap/genom_results.txt', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/qualimap/raw_data_qualimapReport/', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/qualimap/images_qualimapReport/', sample = layouts.keys()),
        expand('Data_Preprocessing/BAMs/{sample}/qualimap/css/', sample = layouts.keys()),
        'Data_Preprocessing/qualimap.list',
        'Data_Preprocessing/qualimap/',
        expand('Variant_Calling/HaplotypeCaller/{sample}/{sample}.g.vcf.gz', sample = layouts.keys()),
        expand('Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf', sample = layouts.keys()),
        expand('Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf.gz', sample = layouts.keys()),
        expand('Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf.gz.tbi', sample = layouts.keys()),
        expand('Variant_Calling/GenotypeGVCFs/equine_{interval}/output.vcf.gz', interval=[str(i).zfill(3) for i in range(0,228+1)]),
        'Variant_Calling/GatherVcfsCloud/output.vcf',
        'Variant_Calling/SelectVariants/SNPs/snps_only.vcf.gz',
        'Variant_Calling/SelectVariants/Indels/indels_only.vcf.gz',
        'Variant_Calling/VariantFiltration/indels_filtered.vcf.gz',
        'Variant_Calling/ExcessHet/excess_het.vcf.gz',
        'Variant_Calling/MakeSitesOnly/sitesonly.vcf.gz',
        'Variant_Calling/VariantRecalibrator/output.recal',
        'Variant_Calling/VariantRecalibrator/output.tranches',
        'Variant_Calling/ApplyVQSR/vqsr_recal.vcf.gz',
        'Variant_Calling/GatherVcfsCloud/output.vcf',
        'Variant_Calling/SnpEff/indels.snpeff.vcf',
        'Variant_Calling/SnpEff/indels.snpeff.csv',
        'Variant_Calling/SnpEff/indels.snpeff.html',
        'Variant_Calling/SnpEff/indels.snpeff.genes.txt',
        'Variant_Calling/SnpEff/snps.snpeff.vcf',
        'Variant_Calling/SnpEff/snps.snpeff.csv',
        'Variant_Calling/SnpEff/snps.snpeff.html',
        'Variant_Calling/SnpEff/snps.snpeff.genes.txt',
        'Variant_Calling/SnpEff/snps.snpeff.vcf.gz',
        'Variant_Calling/SnpEff/indels.snpeff.vcf.gz',
        'Variant_Calling/SnpEff/snps.snpeff.vcf.gz.tbi',
        'Variant_Calling/SnpEff/indels.snpeff.vcf.gz.tbi',

rule fastq_to_sam:
    input:
        unpack(get_fastqs),
    output:
        ubam = 'Data_Preprocessing/uBAMs/{sample}/{sample}.fastqtosam.bam',
    resources:
        time    = 120,
        mem_mb  = 60000,
        cpus    = 4,
    run:
        infile = gzip.open(str(input.r1), 'rt')
        head = infile.readline().strip()

        d2 = {}

        tmp = head.split(':')

        d2['sample_name'] = wildcards.sample
        d2['read_group'] = f'{wildcards.sample}_A'
        d2['plat_unit'] = f'{tmp[2]}.{tmp[3]}.{tmp[9]}'
        d2['library'] = f'Illumina.{wildcards.sample}'

        # the layouts dictionary was created at the top of this script
        # to identify single-end and paired-end samples for appropriate processing here
        if layouts[wildcards.sample] == 'single':
            shell(f'''
                java -jar src/picard.jar FastqToSam \
                    --FASTQ {{input.r1}} \
                    --OUTPUT {{output.ubam}} \
                    --READ_GROUP_NAME {d2['read_group']}\
                    --SAMPLE_NAME {{wildcards.sample}} \
                    --LIBRARY_NAME {d2['library']} \
                    --PLATFORM ILLUMINA \
                    --PLATFORM_UNIT {d2['plat_unit']}
                    ''')
        elif layouts[wildcards.sample] == 'paired':
            shell(f'''
                java -jar src/picard.jar FastqToSam \
                    --FASTQ {input.r1} \
                    --FASTQ2 {input.r2} \
                    --OUTPUT {{output.ubam}} \
                    --READ_GROUP_NAME {d2['read_group']}\
                    --SAMPLE_NAME {{wildcards.sample}} \
                    --LIBRARY_NAME {d2['library']} \
                    --PLATFORM ILLUMINA \
                    --PLATFORM_UNIT {d2['plat_unit']} 
                    ''')

rule mark_adapters:
    input:
        ubam = {rules.fastq_to_sam.output.ubam},
    output:
        marked = 'Data_Preprocessing/Marked_uBAMs/{sample}/{sample}.markilluminaadapters.bam',
        metrics = 'Data_Preprocessing/Marked_uBAMs/{sample}/{sample}.txt/',
    resources:
        time    = 120,
        mem_mb    = 60000,
        cpus    = 4,
    shell:
        '''
            java -jar src/picard.jar MarkIlluminaAdapters \
            I={input.ubam} \
            M={output.metrics} \
            O={output.marked}
        '''

rule align:
    input:
        ubam = {rules.mark_adapters.output.marked},
    output:
        bam = 'Data_Preprocessing/BAMs/{sample}/{sample}.aligned.bam',
    params:
        ref = config[ref],
    threads: 16
    resources:
        time    = 1440,
        mem_mb    = 60000,
        cpus    = 4,
    shell:
        '''
            java -Xmx8G -jar src/picard.jar SamToFastq \
            I={input.ubam} \
            FASTQ=/dev/stdout \
            CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true | \
 
            bwa mem -M -t 7 -p {params.ref} /dev/stdin | \

            java -Xmx16G -jar src/picard.jar MergeBamAlignment \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM={input.ubam} \
            OUTPUT={output.bam} \
            R={params.ref} CREATE_INDEX=true ADD_MATE_CIGAR=true \
            CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
            INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS
        '''

rule mark_duplicates:
    input:
        bam = {rules.align.output.bam},
    output:
        dup_marked_bam = 'Data_Preprocessing/Marked_BAMs/{sample}/{sample}.MarkDuplicates.bam',
        metrics = 'Data_Preprocessing/Marked_BAMs/{sample}/{sample}.txt',
    params:
        tmp_dir = '/dev/shm/{sample}.mark_dup.tmp',
    threads: 4
    resources:
        time    = 360,
        mem_mb    = 60000,
        cpus    = 4,
    shell:
        '''
            mkdir -p {params.tmp_dir}
            gatk MarkDuplicates \
                -I {input.bam} \
                -O {output.dup_marked_bam} \
                -M {output.metrics} \
                --TMP_DIR {params.tmp_dir}
        '''

rule bqsr:
    input:
        dup_marked_bam = {rules.mark_duplicates.output.dup_marked_bam},
    output:
        table = 'Data_Preprocessing/BQSR/BaseRecalibrator/{sample}/{sample}.table',
    params:
        ref = config[ref],
        known_sites = config[known_sites],
    resources:
        time    = 360,
        mem_mb    = 60000,
        cpus    = 4,
    shell:
        '''
            gatk BaseRecalibrator \
                -I {input.dup_marked_bam} \
                -R {params.ref} \
                --known-sites {params.known_sites} \
                -O {output.table}
        '''

rule apply_bqsr:
    input:
        dup_marked_bam = {rules.mark_duplicates.output.dup_marked_bam},
        table = {rules.bqsr.output.table},
    output:
        recal_bam = 'Data_Preprocessing/BQSR/ApplyBQSR/{sample}/{sample}.recalibrated.bam',
    params:
        ref = config[ref],
    threads: 4
    resources:
        time    = 360,
        mem_mb    = 60000,
        cpus    = 4,
    shell:
        '''
            gatk ApplyBQSR \
                -R {params.ref} \
                -I {input.dup_marked_bam} \
                --bqsr-recal-file {input.table} \
                -O {output.recal_bam}
        '''

rule qualimap:
    input:
        final_bam = {rules.apply_bqsr.output.recal_bam},
    output:
        report_pdf = 'Data_Preprocessing/BAMs/{sample}/qualimap/report.pdf',
        report_html = 'Data_Preprocessing/BAMs/{sample}/qualimap/report.html',
        report_txt = 'Data_Preprocessing/BAMs/{sample}/qualimap/genom_results.txt',
        raw_dat = 'Data_Preprocessing/BAMs/{sample}/qualimap/raw_data_qualimapReport/',
        figures = 'Data_Preprocessing/BAMs/{sample}/qualimap/images_qualimapReport/',
        css = 'Data_Preprocessing/BAMs/{sample}/qualimap/css/',
    params:
        out_dir = 'Data_Preprocessing/BAMs/{sample}/qualimap/',
        regions = config[qualimap_regions],
    threads: 12
    resources:
        time    = 120,
        mem_mb  = 24000,
    shell:
        '''
            unset DISPLAY

            qualimap bamqc \
                -gff {params.regions} \
                -bam {input.final_bam} \
                -outdir {params.out_dir} \
                -outformat PDF:HTML \
                -nt 12 \
                --java-mem-size=24G
        '''

rule qualimap_list:
    input:
        stats_dir = expand('Data_Preprocessing/BAMs/{sample}/qualimap/', sample = layouts.keys()),
    output:
        stats_list = 'Data_Preprocessing/qualimap.list',
    threads: 6
    resources:
        time    = 120,
        mem_mb  = 24000,
    run:
        outfile = open('qualimap.list', 'wt')
        samples = []

        for root,dirs,files in os.walk('Data_Preprocessing/BAMs/'):
            for dir in dirs:
                path = os.path.join(root, dir)
                split = path.split('/')
                if split[2] not in samples:
                    samples.append(split[2])
                    print(split[2] + '\t' +
                    split[0] + '/' +
                    split[1] + '/' +
                    split[2] + '/' +
                    'qualimap', file=outfile)

        oufile.close()

rule qualimap_multi:
    input:
        stats_list = {rules.qualimap_list.output.stats_list},
    output:
        multi_report = 'Data_Preprocessing/qualimap/',
    threads: 12
    resources:
        time    = 120,
        mem_mb  = 60000,
    shell:
        '''
            qualimap multi-bamqc \
                -d {input.stats_list} \
                -outdir {output.multi_report} \
                --java-mem-size=60G
        '''


rule hap_calling:
    input:
        #bam = '/scratch.global/dimml002/wgs_BAMs/{sample}.goldenPath.bam',
        bam = 'Data_Preprocessing/BQSR/ApplyBQSR/{sample}/{sample}.recalibrated.bam', 
    output:
        gvcf = 'Variant_Calling/HaplotypeCaller/{sample}/{sample}.g.vcf.gz',
    params:
        ref = config[ref],
    threads: 12
    resources:
        time    = 2400,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk --java-options '-Xmx4g' HaplotypeCaller  \
                -R {params.ref} \
                -I {input.bam} \
                -O {output.gvcf} \
                --disable-read-filter NotDuplicateReadFilter \
                --emit-ref-confidence GVCF
        '''

rule subset_targets_gvcfs:
    input:
        gvcf = {rules.hap_calling.output.gvcf},
    output:
        subset_gvcf = 'Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf',
    params:
        target_file = config[target_file],
    shell:
        '''
            bcftools view \
                {input.gvcf} \
                -R {params.target_file} \
                 > {output.subset_gvcf}

        '''

rule zip_and_index_gvcfs:
    input:
        gvcf = {rules.subset_targets_gvcfs.output.subset_gvcf},
    output:
        zipped = 'Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf.gz',
        indexed = 'Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf.gz.tbi',
    shell:
        '''
            bgzip {input.gvcf}
            gatk IndexFeatureFile -I {output.zipped}
        '''

# Create a list of GVCF paths for GenomicsDBImport.
rule gvcf_list:
    input:
    #    gvcfs = expand('Variant_Calling/HaplotypeCaller/{sample}/{sample}.g.vcf.gz', sample = layouts.keys()),
        gvcfs = expand('Variant_Calling/targets_only_gvcfs/{sample}/{sample}.g.vcf.gz', sample = layouts.keys()),
    output:
        gvcf_list_output = 'Variant_Calling/gvcf_list.txt',
    run:
        with open(output.gvcf_list_output, 'w') as f:
            for i in input.gvcfs:
                horse = os.path.basename(i).split('.')[0]
                # gvcf = i.replace('done','g.vcf.gz')
                f.write(horse + '\t' + i + '\n') 
        

rule genomics_db_import:
    input:
        gvcf_list = {rules.gvcf_list.output.gvcf_list_output},
        intervals = '/panfs/jay/groups/27/mccuem/shared/Projects/RER/GenIntervals/Intervals_Renamed/equine_{interval}.interval_list',
    output:
        db = directory('Variant_Calling/GenomicsDBImport/equine_{interval}.interval_list/'),
    threads: 16
    resources:
        time    = 1440,
        mem_mb     = 120000,
    shell:
        '''
            gatk --java-options '-Xmx50g' \
            GenomicsDBImport \
                --genomicsdb-workspace-path {output.db} \
                --batch-size 50 \
                --intervals {input.intervals} \
                --sample-name-map {input.gvcf_list} \
                --reader-threads 5 \
                --genomicsdb-shared-posixfs-optimizations true
        '''

rule genotype_gvcfs:
    input:
        db = {rules.genomics_db_import.output.db},
        interval = '/panfs/jay/groups/27/mccuem/shared/Projects/RER/GenIntervals/Intervals_Renamed/equine_{interval}.interval_list/',
    output:
        vcf = 'Variant_Calling/GenotypeGVCFs/equine_{interval}/output.vcf.gz',
    params:
        ref = config[ref],
    threads: 16
    resources:
        time    = 1440,
        mem_mb    = 60000,
    shell:
        '''
            gatk --java-options '-Xms50g' \
            GenotypeGVCFs \
                -R {params.ref} \
                -O {output.vcf} \
                -G StandardAnnotation \
                --only-output-calls-starting-in-intervals \
                -V gendb://{input.db} \
                -L {input.interval}
        '''

# Create a list of per-interval VCFs for GatherVcfsCloud.
rule vcf_list:
    input:
        vcfs = expand('Variant_Calling/GenotypeGVCFs/equine_{interval}/output.vcf.gz', interval=[str(i).zfill(3) for i in range(0,228+1)]),
    output:
        vcf_list = 'Variant_Calling/vcf.list',
    run:
        outfile = open(output.vcf_list, 'wt')
        interval=[str(i).zfill(3) for i in range(0,228+1)]
        for num in interval:
            print('Variant_Calling/GenotypeGVCFs/equine_' + num + '/output.vcf.gz', file=outfile)
        outfile.close()

rule gather_vcfs:
    input:
        vcf_list = {rules.vcf_list.output.vcf_list}, 
    output:
        gathered_vcf = 'Variant_Calling/GatherVcfsCloud/output.vcf',
    threads: 4
    resources:
        time    = 600,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk GatherVcfsCloud \
                -I {input.vcf_list} \
                -O {output.gathered_vcf}
        '''

rule select_snps:
    input:
        vcf = {rules.gather_vcfs.output.gathered_vcf},
    output:
        snp_vcf = 'Variant_Calling/SelectVariants/SNPs/snps_only.vcf.gz',
    threads: 4
    resources:
        time    = 600,
        mem_mb  = 24000,
        cpus    = 4,
    shell:
        '''
            gatk SelectVariants \
                -V {input.vcf} \
                -select-type SNP \
                -O {output.snp_vcf}
        '''

rule select_indels:
    input:
        vcf = {rules.gather_vcfs.output.gathered_vcf},
    output:
        indel_vcf = 'Variant_Calling/SelectVariants/Indels/indels_only.vcf.gz',
    threads: 4
    resources:
        time    = 600,
        mem_mb  = 24000,
        cpus    = 4,
    shell:
        '''
            gatk SelectVariants \
                -V {input.vcf} \
                -select-type INDEL \
                -select-type MIXED \
                -O {output.indel_vcf}
'''

rule hard_filter:
    input:
        indel_vcf = {rules.select_indels.output.indel_vcf},
    output:
        filtered_vcf = 'Variant_Calling/VariantFiltration/indels_filtered.vcf.gz',
    threads: 4
    resources:
        time    = 600,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk VariantFiltration \
                -V {input.indel_vcf} \
                -filter 'QD < 2.0' --filter-name 'QD2' \
                -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
                #-filter 'SOR > 3.0' --filter-name 'SOR3' \
                -filter 'FS > 200.0' --filter-name 'FS200' \
                #-filter 'MQ < 40.0' --filter-name 'MQ40' \
                #-filter 'MQRankSum < 12.5' --filter-name 'MQRankSum-12.5' \
                -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \
                -O {output.filtered_vcf}
        '''

rule annotate_indels:
    input:
        vcf = {rules.hard_filter.output.filtered_vcf},
    output:
        ann_vcf = 'Variant_Calling/SnpEff/indels.snpeff.vcf',
        csv = 'Variant_Calling/SnpEff/indels.snpeff.csv',
        stats = 'Variant_Calling/SnpEff/indels.snpeff.html',
        genes = 'Variant_Calling/SnpEff/indels.snpeff.genes.txt',
    params:
        config = config[SnpEff_config],
    threads: 4
    resources:
        time    = 600,
        mem_mb  = 24000,
        cpus    = 4,
    shell:
        '''
            java -Xmx4g -Djava.io.tmpdir=/home/mccuem/shared/Projects/HorseGenomeProject/Data/temp_files \
                -jar src/snpEff/snpEff.jar \
                -C {params.config} \
                -csvstats {output.csv} \
                -stats {output.stats} \
                EquCab3.0.99 {input.vcf} > {output.ann_vcf}

        '''

rule excess_het:
    input:
        vcf = {rules.select_snps.output.snp_vcf},
    output:
        het_vcf = 'Variant_Calling/ExcessHet/excess_het.vcf.gz',
    threads: 4
    resources:
        time    = 600,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk --java-options '-Xmx3g -Xms3g' VariantFiltration \
                -V {input.vcf} \
                --filter-expression 'ExcessHet > 54.69' \
                --filter-name ExcessHet \
                -O {output.het_vcf}
        '''

rule sites_only:
    input:
        vcf = {rules.excess_het.output.het_vcf},
    output:
        sitesonly = 'Variant_Calling/MakeSitesOnly/sitesonly.vcf.gz',
    threads: 4
    resources:
        time    = 600,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk MakeSitesOnlyVcf \
                -I {input.vcf} \
                -O {output.sitesonly}
        '''

rule variant_recalibrator:
    input:
        snp_vcf = {rules.sites_only.output.sitesonly},
    output:
        recal = 'Variant_Calling/VariantRecalibrator/output.recal',
        tranches = 'Variant_Calling/VariantRecalibrator/output.tranches',
    params:
        ref = config[ref],
        dbsnp = config[dbsnp],
        known_sites = config[known_sites],
    threads: 4
    resources:
        time    = 600,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk VariantRecalibrator \
                -R {params.ref} \
                -V {input.snp_vcf} \
                --trust-all-polymorphic \
                -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 \
                 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
                -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
                -mode SNP \
                --max-gaussians 4 \
                --resource:snpchip,known=false,training=true,truth=true,prior=12 {params.known_sites} \
                --resource:dbsnp,known=false,training=true,truth=false,prior=10 {params.dbsnp} \
                -O {output.recal} \
                --tranches-file {output.tranches}
        '''

rule apply_vqsr:
    input:
        vcf = {rules.excess_het.output.het_vcf},
        recal = {rules.variant_recalibrator.output.recal},
        tranches = {rules.variant_recalibrator.output.tranches},
    output:
        recal_vcf = 'Variant_Calling/ApplyVQSR/vqsr_recal.vcf.gz',
    threads: 4
    resources:
        time    = 600,
        mem_mb    = 24000,
        cpus    = 4,
    shell:
        '''
            gatk --java-options '-Xmx5g -Xms5g' ApplyVQSR \
                -V {input.vcf} \
                --recal-file {input.recal} \
                --tranches-file {input.tranches} \
                --truth-sensitivity-filter-level 99.7 \
                --create-output-variant-index true \
                -mode SNP \
                -O {output.recal_vcf}
        '''

rule annotate_snps:
    input:
        vcf = {rules.apply_vqsr.output.recal_vcf},
    output:
        ann_vcf = 'Variant_Calling/SnpEff/snps.snpeff.vcf',
        csv = 'Variant_Calling/SnpEff/snps.snpeff.csv',
        stats = 'Variant_Calling/SnpEff/snps.snpeff.html',
        genes = 'Variant_Calling/SnpEff/snps.snpeff.genes.txt',
    params:
        config = config[SnpEff_config],
    threads: 4
    resources:
        time    = 1440,
        mem_mb  = 24000,
        cpus    = 4,
    shell:
        '''
            java -Xmx4g -Djava.io.tmpdir=/home/mccuem/shared/Projects/HorseGenomeProject/Data/temp_files \
                -jar src/snpEff/snpEff.jar \
                -C {params.config} \
                -csvstats {output.csv} \
                -stats {output.stats} \
                EquCab3.0.99 {input.vcf} > {output.ann_vcf}

        '''

rule combine_vcfs:
    input:
        indels = {rules.hard_filter.output.filtered_vcf},
        snps = {rules.apply_vqsr.output.recal_vcf},
    output:
        combined_vcf = 'Variant_Calling/GatherVcfsCloud/Final/output.vcf',
    threads: 4
    resources:
        time    = 600,
        mem_mb  = 24000,
        cpus    = 4,
    shell:
        '''
            java -jar src/picard.jar MergeVcfs \
                I= {input.indels} \
                I= {input.snps} \
                O= {output.combined_vcf}
        '''

