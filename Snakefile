'''Snakefile for GWAS Variant and Sample QC Version 0.1'''

configfile: "config.yaml"
BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
SAMPLE = config['sample']

rule all:
    input:
        expand("temp/{sample}_filtered_PCA.{ext}", ext = ['eigenval', 'eigenvec'], sample=SAMPLE),
        expand("temp/{sample}_exclude.samples", sample=SAMPLE),
        expand("stats/{sample}_GWAS_QC.html", sample=SAMPLE)

## ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input:
        expand("Data/{sample}.{ext}", ext = BPLINK, sample=SAMPLE)
    output:
        temp(expand("temp/{{sample}}_SnpQc.{ext}", ext = BPLINK))
    params:
        indat = 'Data/{sample}',
        out = 'temp/{sample}_SnpQc'
    shell:
        'plink --bfile {params.indat} --geno 0.05 --maf 0.01 --hardy --make-bed --out {params.out}'

## ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input:
        rules.snp_qc.output
    output:
        temp(expand("temp/{{sample}}_callRate.{ext}", ext = BPLINK)),
        "temp/{sample}_callRate.irem"
    params:
        indat = rules.snp_qc.params.out,
        out = 'temp/{sample}_callRate'
    shell:
        'plink --bfile {params.indat} --mind 0.05 --make-bed --out {params.out}'

## ---- Exclude Samples with discordant sex ----
##  Use ADNI hg18 data, as the liftover removed the x chromsome data
rule sexcheck_QC:
    input:
        expand("Data/{{sample}}_xy.{ext}", ext = BPLINK)
    output:
        "temp/{sample}_SexQC.sexcheck"
    params:
        indat = 'Data/{sample}_xy',
        out = "temp/{sample}_SexQC"
    shell:
        'plink --bfile {params.indat} --check-sex --out {params.out}'

rule sex_sample_fail:
    input:
        "temp/{sample}_SexQC.sexcheck"
    output:
        "temp/{sample}_exclude.sexcheck"
    params:
        indat = "temp/{sample}_SexQC.sexcheck",
        out = "temp/{sample}_exclude.sexcheck"
    shell:
        'Rscript scripts/sexcheck_QC.R {params.indat} {params.out}'

rule sex__exclude_failed:
    input:
        "temp/{sample}_exclude.sexcheck",
        expand("temp/{{sample}}_callRate.{ext}", ext = BPLINK)
    output:
        temp(expand("temp/{{sample}}_SexExclude.{ext}", ext = BPLINK))
    params:
        indat_exclude = "temp/{sample}_exclude.sexcheck",
        indat_plink = 'temp/{sample}_callRate',
        out = 'temp/{sample}_SexExclude'
    shell:
        'plink --bfile {params.indat_plink} --remove {params.indat_exclude} --make-bed --out {params.out}'

## ---- Prune SNPs, autosome only ----
##  Pruned SNP list is used for IBD, PCA and heterozigosity calculations
rule Prune_snps:
    input:
        expand("temp/{{sample}}_SexExclude.{ext}", ext = BPLINK)
    output:
        temp(expand("temp/{{sample}}_thinned.{ext}", ext = ['prune.in', 'prune.out']))
    params:
        indat = 'temp/{sample}_SexExclude',
        out = "temp/{sample}_thinned"
    shell:
        'plink --bfile {params.indat} --autosome --indep 50 5 1.5 --out {params.out}'

## ---- Exclude Samples with interealtedness ----
rule relatedness_QC:
    input:
        expand("temp/{{sample}}_SexExclude.{ext}", ext = BPLINK),
        "temp/{sample}_thinned.prune.in"
    output:
        "temp/{sample}_IBDQC.genome"
    params:
        indat_plink = 'temp/{sample}_SexExclude',
        indat_prune_in = "temp/{sample}_thinned.prune.in",
        out = 'temp/{sample}_IBDQC'
    shell:
        'plink --bfile {params.indat_plink} --extract {params.indat_prune_in} --genome --min 0.05 --out {params.out}'

rule relatedness_sample_fail:
    input:
        "temp/{sample}_IBDQC.genome",
        "temp/{sample}_SnpQc.fam"
    output:
        "temp/{sample}_exclude.relatedness"
    params:
        indat_genome = "temp/{sample}_IBDQC.genome",
        indat_fam = "temp/{sample}_SnpQc.fam",
        out = "temp/{sample}_exclude.relatedness"
    shell:
        'Rscript scripts/relatedness_QC.R {params.indat_genome} {params.indat_fam} {params.out}'

rule relatedness_exclude_failed:
    input:
        "temp/{sample}_exclude.relatedness",
        expand("temp/{{sample}}_SexExclude.{ext}", ext = BPLINK)
    output:
        temp(expand("temp/{{sample}}_RelatednessExclude.{ext}", ext = BPLINK))
    params:
        indat_exclude = "temp/{sample}_exclude.relatedness",
        indat_plink = 'temp/{sample}_SexExclude',
        out = 'temp/{sample}_RelatednessExclude'
    shell:
        'plink --bfile {params.indat_plink} --remove {params.indat_exclude} --make-bed --out {params.out}'

## ---- Exclude Samples with outlying heterozigosity ----
rule heterozigosity_QC:
    input:
        expand("temp/{{sample}}_RelatednessExclude.{ext}", ext = BPLINK),
        "temp/{sample}_thinned.prune.in"
    output:
        "temp/{sample}_HetQC.het"
    params:
        indat_plink = 'temp/{sample}_RelatednessExclude',
        indat_prune_in = "temp/{sample}_thinned.prune.in",
        out = 'temp/{sample}_HetQC'
    shell:
        'plink --bfile {params.indat_plink} --extract {params.indat_prune_in} --het --out {params.out}'

rule heterozigosity_sample_fail:
    input:
        "temp/{sample}_HetQC.het",
    output:
        "temp/{sample}_exclude.heterozigosity"
    params:
        indat_het = "temp/{sample}_HetQC.het",
        out = "temp/{sample}_exclude.heterozigosity"
    shell:
        'Rscript scripts/heterozygosity_QC.R {params.indat_het} {params.out}'

rule heterozigosity_exclude_failed:
    input:
        "temp/{sample}_exclude.heterozigosity",
        expand("temp/{{sample}}_RelatednessExclude.{ext}", ext = BPLINK)
    output:
        temp(expand("temp/{{sample}}_HetExclude.{ext}", ext = BPLINK))
    params:
        indat_exclude = "temp/{sample}_exclude.heterozigosity",
        indat_plink = 'temp/{sample}_RelatednessExclude',
        out = 'temp/{sample}_HetExclude'
    shell:
        'plink --bfile {params.indat_plink} --remove {params.indat_exclude} --make-bed --out {params.out}'

## Prune sample dataset
rule sample_prune:
    input:
        expand("temp/{{sample}}_HetExclude.{ext}", ext = BPLINK),
        "temp/{sample}_thinned.prune.in"
    output:
        temp(expand("temp/{{sample}}_samp_thinned.{ext}", ext = BPLINK))
    params:
        indat_plink = 'temp/{sample}_HetExclude',
        indat_prune_in = "temp/{sample}_thinned.prune.in",
        out = 'temp/{sample}_samp_thinned'
    shell:
        'plink --bfile {params.indat_plink} --extract {params.indat_prune_in} --make-bed --out {params.out}'

## align 1000 genomes to fasta refrence
rule Sample_Flip:
    input:
        bim = "temp/{sample}_samp_thinned.bim",
        bed = "temp/{sample}_samp_thinned.bed",
        fam = "temp/{sample}_samp_thinned.fam",
        fasta = "Data/hg19.fa"
    output:
        temp(expand("temp/{{sample}}_samp_thinned_flipped.{ext}", ext = BPLINK))
    shell:
        "/Users/sheaandrews/Programs/flippyr/flippyr.py -p {input.fasta} {input.bim}"

## Recode sample plink file to vcf
rule Sample_Plink2Bcf:
    input:
        expand("temp/{{sample}}_samp_thinned_flipped.{ext}", ext = BPLINK)
    output:
        temp("temp/{sample}_samp_thinned_flipped.vcf.gz")
    params:
        indat = "temp/{sample}_samp_thinned_flipped",
        out = "temp/{sample}_samp_thinned_flipped"
    shell:
        'plink --bfile {params.indat} --recode vcf bgz --keep-allele-order --real-ref-alleles --out {params.out}'

## Index bcf
rule Sample_IndexBcf:
    input:
        bcf = "temp/{sample}_samp_thinned_flipped.vcf.gz"
    output:
        temp("temp/{sample}_samp_thinned_flipped.vcf.gz.csi")
    shell:
        'bcftools index -f {input.bcf}'

## ---- Principal Compoent analysis ----
##  Project ADNI onto a PCA using the 1000 Genomes dataset to identify population outliers

##  Extract a pruned dataset from 1000 genomes using the same pruning SNPs from Sample
rule Reference_prune:
    input:
        expand("Data/1000genomes_allChr.{ext}", ext = BPLINK),
        "temp/{sample}_thinned.prune.in"
    output:
        temp(expand("temp/{{sample}}_1kg_thinned.{ext}", ext = BPLINK))
    params:
        indat_plink = 'Data/1000genomes_allChr',
        indat_prune_in = "temp/{sample}_thinned.prune.in",
        out = 'temp/{sample}_1kg_thinned'
    shell:
        'plink --bfile {params.indat_plink} --filter-founders --extract {params.indat_prune_in} --make-bed --out {params.out}'


## align 1000 genomes to fasta refrence
rule Reference_flip:
    input:
        bim = "temp/{sample}_1kg_thinned.bim",
        bed = "temp/{sample}_1kg_thinned.bed",
        fam = "temp/{sample}_1kg_thinned.fam",
        fasta = "Data/hg19.fa"
    output:
        temp(expand("temp/{{sample}}_1kg_thinned_flipped.{ext}", ext = BPLINK))
    shell:
        "/Users/sheaandrews/Programs/flippyr/flippyr.py -p {input.fasta} {input.bim}"

## Recode 1kg to vcf
rule Reference_Plink2Bcf:
    input:
        expand("temp/{{sample}}_1kg_thinned_flipped.{ext}", ext = BPLINK)
    output:
        temp(expand("temp/{{sample}}_1kg_thinned_flipped.{ext}", ext = 'vcf.gz'))
    params:
        indat = "temp/{sample}_1kg_thinned_flipped",
        out = "temp/{sample}_1kg_thinned_flipped"
    shell:
        'plink --bfile {params.indat} --recode vcf bgz --keep-allele-order --real-ref-alleles --out {params.out}'

## Index bcf
rule Reference_IndexBcf:
    input:
        bcf = "temp/{sample}_1kg_thinned_flipped.vcf.gz"
    output:
        temp("temp/{sample}_1kg_thinned_flipped.vcf.gz.csi")
    shell:
        'bcftools index -f {input.bcf}'

## Merge ref and sample
rule Merge_RefenceSample:
    input:
        bcf_1kg = "temp/{sample}_1kg_thinned_flipped.vcf.gz",
        csi_1kg = "temp/{sample}_1kg_thinned_flipped.vcf.gz.csi",
        bcf_samp = "temp/{sample}_samp_thinned_flipped.vcf.gz",
        csi_samp = "temp/{sample}_samp_thinned_flipped.vcf.gz.csi"
    output:
        out = temp("temp/{sample}_1kg_merged.vcf")
    shell:
        'bcftools merge -m none {input.bcf_1kg} {input.bcf_samp} -Ov -o {output.out}'

## recode merged sample to plink
rule Plink_RefenceSample:
    input:
        vcf = "temp/{sample}_1kg_merged.vcf",
    output:
        temp(temp(expand("temp/{{sample}}_1kg_merged.{ext}", ext = BPLINK)))
    params:
        out = "temp/{sample}_1kg_merged"
    shell:
        'plink --vcf {input.vcf} --const-fid --make-bed --out {params.out}'

rule fix_fam:
    input:
        fam = "temp/{sample}_1kg_merged.fam",
    output:
        out = temp("temp/{sample}_1kg_merged_fixed.fam")
    shell:
        'scripts/fix_fam.py {input.fam} {output.out}'

## PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        expand("temp/{{sample}}_1kg_merged.{ext}", ext = BPLINK),
        "temp/{sample}_1kg_merged_fixed.fam",
        "Data/1000genomes_pops.txt",
        "Data/pops.txt"
    output:
        expand("temp/{{sample}}_1kg_merged.{ext}", ext = ['eigenval', 'eigenvec'])
    params:
        indat_plink = "temp/{sample}_1kg_merged",
        indat_fam = "temp/{sample}_1kg_merged_fixed.fam",
        indat_pop = "Data/1000genomes_pops.txt",
        indat_clust = "Data/pops.txt",
        out = "temp/{sample}_1kg_merged"
    shell:
        'plink --bfile {params.indat_plink} --fam {params.indat_fam} --pca 10 --within {params.indat_pop} --pca-clusters {params.indat_clust} --out {params.out}'

## Rscript to identify EUR population outliers
rule ExcludePopulationOutliers:
    input:
        indat_eigenval = "temp/{sample}_1kg_merged.eigenval",
        indat_eigenvec = "temp/{sample}_1kg_merged.eigenvec",
        indat_fam = "temp/{sample}_samp_thinned.fam",
        indat_1kgped = "Data/20130606_g1k.ped"
    output:
        out = "temp/{sample}_exclude.pca"
    shell:
        'Rscript scripts/PCA_QC.R {input.indat_eigenvec} {input.indat_1kgped} {input.indat_fam} {input.indat_eigenval} {output.out} '

## Run PCA to for population stratification
rule PopulationStratification:
    input:
        expand("temp/{{sample}}_samp_thinned.{ext}", ext = BPLINK),
        "temp/{sample}_exclude.pca"
    output:
        expand("temp/{{sample}}_filtered_PCA.{ext}", ext = ['eigenval', 'eigenvec'])
    params:
        indat = "temp/{sample}_samp_thinned",
        exclude = "temp/{sample}_exclude.pca",
        out = "temp/{sample}_filtered_PCA"
    shell:
        'plink --bfile {params.indat}  --remove {params.exclude} --pca 10 --out {params.out}'

rule SampleExclusion:
    input:
        SampCallRate = "temp/{sample}_callRate.irem",
        het = "temp/{sample}_exclude.heterozigosity",
        sex = "temp/{sample}_exclude.sexcheck",
        pca = "temp/{sample}_exclude.pca",
        relat = "temp/{sample}_exclude.relatedness"
    output:
        out = "temp/{sample}_exclude.samples"
    shell:
        'Rscript scripts/sample_QC.R {input.SampCallRate} {input.het} {input.sex} {input.pca} {input.relat} {output.out}'



rule GWAS_QC_Report:
    input:
        "scripts/GWAS_QC.Rmd",
        "temp/{sample}_SnpQc.hwe",
        "temp/{sample}_IBDQC.genome",
        "Data/{sample}.fam",
        "temp/{sample}_1kg_merged.eigenval",
        "temp/{sample}_1kg_merged.eigenvec",
        "temp/{sample}_samp_thinned.fam",
        "Data/20130606_g1k.ped",
        "temp/{sample}_filtered_PCA.eigenval",
        "temp/{sample}_filtered_PCA.eigenvec"
    output:
        "stats/{sample}_GWAS_QC.html"
    params:
        rwd = RWD,
        Sample = "{sample}",
        output_dir = "stats",
        Path_SexFile = "temp/{sample}_SexQC.sexcheck",
        Path_hwe = "temp/{sample}_SnpQc.hwe",
        Path_HetFile = "temp/{sample}_HetQC.het",
        Path_GenomeFile = "temp/{sample}_IBDQC.genome",
        Path_FamFile = "Data/{sample}.fam",
        Path_eigenval = "temp/{sample}_1kg_merged.eigenval",
        Path_eigenvec = "temp/{sample}_1kg_merged.eigenvec",
        Path_TargetPops = "temp/{sample}_samp_thinned.fam",
        PATH_BasePops = "Data/20130606_g1k.ped",
        Path_PopStrat_eigenval = "temp/{sample}_filtered_PCA.eigenval",
        Path_PopStrat_eigenvec = "temp/{sample}_filtered_PCA.eigenvec",
        out = "{sample}_GWAS_QC.html",
    shell:
        "R -e 'rmarkdown::render("
        """"{input[0]}", output_file = "{params.out}", output_dir = "{params.output_dir}", \
        params = list(rwd = "{params.rwd}", Sample = "{params.Sample}", \
        Path_SexFile = "{params.Path_SexFile}", \
        Path_hwe = "{params.Path_hwe}", \
        Path_HetFile = "{params.Path_HetFile}", \
        Path_GenomeFile = "{params.Path_GenomeFile}", Path_FamFile = "{params.Path_FamFile}", \
        Path_eigenval = "{params.Path_eigenval}", Path_eigenvec = "{params.Path_eigenvec}", Path_TargetPops = "{params.Path_TargetPops}", PATH_BasePops = "{params.PATH_BasePops}", \
        Path_PopStrat_eigenval = "{params.Path_PopStrat_eigenval}", Path_PopStrat_eigenvec = "{params.Path_PopStrat_eigenvec}"))' --slave
        """
