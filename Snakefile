'''Snakefile for GWAS Variant and Sample QC Version 0.1'''

configfile: "config.yaml"
BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
SAMPLE = config['sample']
FAMILY = config['family']
if config['SampleSex']:
    SAMPLESEX = config['SampleSex']
else:
    SAMPLESEX = SAMPLE
DATAIN = config['DataIn']
DATAOUT = config['DataOut']

rule all:
    input:
        expand("{DataOut}/{sample}_filtered_PCA.{ext}", ext = ['eigenval', 'eigenvec'], sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/{sample}_exclude.samples", sample=SAMPLE, DataOut=DATAOUT),
        expand("stats/{sample}_GWAS_QC.html", sample=SAMPLE)

## ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input:
        expand("{DataIn}/{{sample}}.{ext}", ext = BPLINK, DataIn=DATAIN)
    output:
        temp(expand("{DataOut}/{{sample}}_SnpQc.{ext}", ext = BPLINK, DataOut=DATAOUT)),
        expand("{DataOut}/{{sample}}_SnpQc.hwe", sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_SnpQc.frq", sample=SAMPLE, DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_SnpQc.frqx", sample=SAMPLE, DataOut=DATAOUT)
    params:
        indat = '{DataIn}/{{sample}}'.format(DataIn=DATAIN),
        out = "{DataOut}/{{sample}}_SnpQc".format(DataOut=DATAOUT),
    shell:
        'plink --bfile {params.indat} --freq --out {params.out}; \
        plink --bfile {params.indat} --freqx --out {params.out}; \
        plink --bfile {params.indat} --geno 0.05 --maf 0.01 --hardy --make-bed --out {params.out}'

## ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input:
        rules.snp_qc.output
    output:
        temp(expand("{DataOut}/{{sample}}_callRate.{ext}", ext = BPLINK, DataOut=DATAOUT)),
        expand("{DataOut}/{{sample}}_callRate.imiss", sample=SAMPLE, DataOut=DATAOUT),
        touch("{DataOut}/{{sample}}_callRate.irem".format(DataOut=DATAOUT))
    params:
        indat = rules.snp_qc.params.out,
        out = "{DataOut}/{{sample}}_callRate".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat} --mind 0.05 --missing --make-bed --out {params.out}'

## ---- Exclude Samples with discordant sex ----
##  Use ADNI hg18 data, as the liftover removed the x chromsome data
rule sexcheck_QC:
    input:
        expand("{DataIn}/{{SampleSex}}.{ext}", ext = BPLINK, DataIn=DATAIN)
    output:
        "{DataOut}/{{SampleSex}}_SexQC.sexcheck".format(DataOut=DATAOUT)
    params:
        indat = '{DataIn}/{{SampleSex}}'.format(DataIn=DATAIN),
        out = "{DataOut}/{{SampleSex}}_SexQC".format(DataOut=DATAOUT),
    shell:
        'plink --bfile {params.indat} --check-sex --out {params.out}'

rule sex_sample_fail:
    input:
        "{DataOut}/{SampleSex}_SexQC.sexcheck".format(DataOut=DATAOUT, SampleSex=SAMPLESEX)
    output:
        "{DataOut}/{{sample}}_exclude.sexcheck".format(DataOut=DATAOUT),
    shell:
        'Rscript scripts/sexcheck_QC.R {input} {output}'

rule sex__exclude_failed:
    input:
        expand("{DataOut}/{{sample}}_callRate.{ext}", ext=BPLINK, DataOut=DATAOUT),
        indat_exclude = "{DataOut}/{{sample}}_exclude.sexcheck".format(DataOut=DATAOUT)
    output:
        temp(expand("{DataOut}/{{sample}}_SexExclude.{ext}", ext = BPLINK, DataOut=DATAOUT))
    params:
        indat_plink = "{DataOut}/{{sample}}_callRate".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_SexExclude".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --remove {input.indat_exclude} --make-bed --out {params.out}'

## ---- Prune SNPs, autosome only ----
##  Pruned SNP list is used for IBD, PCA and heterozigosity calculations
rule PruneDupvar_snps:
    input:
        expand("{DataOut}/{{sample}}_SexExclude.{ext}", ext = BPLINK, DataOut=DATAOUT),

    output:
        expand("{DataOut}/{{sample}}_thinned.{ext}", ext = ['prune.in', 'prune.out'], DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_thinned.dupvar.delete", DataOut=DATAOUT)
    params:
        indat = "{DataOut}/{{sample}}_SexExclude".format(DataOut=DATAOUT),
        dupvar = "{DataOut}/{{sample}}_thinned.dupvar".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_thinned".format(DataOut=DATAOUT)
    shell:
        "plink --bfile {params.indat} --autosome --indep 50 5 1.5 --list-duplicate-vars --out {params.out}; "
        "Rscript scripts/DuplicateVars.R {params.dupvar}"

## Prune sample dataset
rule sample_prune:
    input:
        expand("{DataOut}/{{sample}}_SexExclude.{ext}", ext = BPLINK, DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_thinned.prune.in", DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_thinned.dupvar.delete", DataOut=DATAOUT),
    output:
        temp(expand("{DataOut}/{{sample}}_samp_thinned.{ext}", ext = BPLINK, DataOut=DATAOUT))
    params:
        indat_plink = "{DataOut}/{{sample}}_SexExclude".format(DataOut=DATAOUT),
        indat_prune_in = "{DataOut}/{{sample}}_thinned.prune.in".format(DataOut=DATAOUT),
        indat_dupvar = "{DataOut}/{{sample}}_thinned.dupvar.delete".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_samp_thinned".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --extract {params.indat_prune_in} --exclude {params.indat_dupvar} --make-bed --out {params.out}'

## ---- Exclude Samples with interealtedness ----
rule relatedness_QC:
    input:
        expand("{DataOut}/{{sample}}_samp_thinned.{ext}", ext = BPLINK, DataOut=DATAOUT),
    output:
        "{DataOut}/{{sample}}_IBDQC.genome".format(DataOut=DATAOUT)
    params:
        indat_plink = "{DataOut}/{{sample}}_samp_thinned".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_IBDQC".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --genome --min 0.05 --out {params.out}'

rule relatedness_sample_fail:
    input:
        indat_genome = expand("{DataOut}/{{sample}}_IBDQC.genome", DataOut=DATAOUT),
        indat_fam = "{DataOut}/{{sample}}_SnpQc.fam".format(DataOut=DATAOUT)
    params:
        Family = FAMILY
    output:
        out = "{DataOut}/{{sample}}_exclude.relatedness".format(DataOut=DATAOUT)
    shell:
        'Rscript scripts/relatedness_QC.R {input.indat_genome} {input.indat_fam} {params.Family} {output.out}'

rule relatedness_exclude_failed:
    input:
        expand("{DataOut}/{{sample}}_exclude.relatedness", DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_samp_thinned.{ext}", ext = BPLINK, DataOut=DATAOUT)
    output:
        temp(expand("{DataOut}/{{sample}}_RelatednessExclude.{ext}", ext = BPLINK, DataOut=DATAOUT))
    params:
        indat_exclude = "{DataOut}/{{sample}}_exclude.relatedness".format(DataOut=DATAOUT),
        indat_plink = "{DataOut}/{{sample}}_samp_thinned".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_RelatednessExclude".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --remove {params.indat_exclude} --make-bed --out {params.out}'

## ---- Exclude Samples with outlying heterozigosity ----
rule heterozigosity_QC:
    input:
        expand("{DataOut}/{{sample}}_RelatednessExclude.{ext}", ext = BPLINK, DataOut=DATAOUT),
    output:
        expand("{DataOut}/{{sample}}_HetQC.het", DataOut=DATAOUT)
    params:
        indat_plink = "{DataOut}/{{sample}}_RelatednessExclude".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_HetQC".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --het --out {params.out}'

rule heterozigosity_sample_fail:
    input:
        expand("{DataOut}/{{sample}}_HetQC.het", DataOut=DATAOUT),
    output:
        "{DataOut}/{{sample}}_exclude.heterozigosity".format(DataOut=DATAOUT)
    shell:
        'Rscript scripts/heterozygosity_QC.R {input} {output}'

rule heterozigosity_exclude_failed:
    input:
        expand("{DataOut}/{{sample}}_exclude.heterozigosity", DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_RelatednessExclude.{ext}", ext = BPLINK, DataOut=DATAOUT)
    output:
        temp(expand("{DataOut}/{{sample}}_HetExclude.{ext}", ext = BPLINK, DataOut=DATAOUT))
    params:
        indat_exclude = "{DataOut}/{{sample}}_exclude.heterozigosity".format(DataOut=DATAOUT),
        indat_plink = "{DataOut}/{{sample}}_RelatednessExclude".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_HetExclude".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --remove {params.indat_exclude} --make-bed --out {params.out}'

## align sample to fasta refrence
rule Sample_Flip:
    input:
        bim = expand("{DataOut}/{{sample}}_HetExclude.bim", DataOut=DATAOUT),
        bed = expand("{DataOut}/{{sample}}_HetExclude.bed", DataOut=DATAOUT),
        fam = expand("{DataOut}/{{sample}}_HetExclude.fam", DataOut=DATAOUT),
        fasta = "data/hg19.fa"
    output:
        temp(expand("{DataOut}/{{sample}}_HetExclude_flipped.{ext}", ext = BPLINK, DataOut=DATAOUT))
    shell:
        "/Users/sheaandrews/Programs/flippyr/flippyr.py -p {input.fasta} {input.bim}"

rule Sample_ChromPosRefAlt:
    input:
        bim = expand("{DataOut}/{{sample}}_HetExclude_flipped.bim", DataOut=DATAOUT),
    output:
        bim = temp("{DataOut}/{{sample}}_HetExclude_flipped_ChromPos.bim".format(DataOut=DATAOUT)),
        snplist = temp("{DataOut}/{{sample}}_thinned_snplist".format(DataOut=DATAOUT))
    shell:
        "Rscript scripts/bim_ChromPosRefAlt.R {input.bim} {output.bim} {output.snplist}"

## Recode sample plink file to vcf
rule Sample_Plink2Bcf:
    input:
        expand("{DataOut}/{{sample}}_HetExclude_flipped.{ext}", ext = BPLINK, DataOut=DATAOUT),
        expand("{DataOut}/{{sample}}_HetExclude_flipped_ChromPos.bim", DataOut=DATAOUT)
    output:
        expand("{DataOut}/{{sample}}_samp_thinned_flipped.vcf.gz", DataOut=DATAOUT)
    params:
        indat = "{DataOut}/{{sample}}_HetExclude_flipped".format(DataOut=DATAOUT),
        bim = "{DataOut}/{{sample}}_HetExclude_flipped_ChromPos.bim".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_samp_thinned_flipped".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat} --bim {params.bim} --recode vcf bgz --keep-allele-order --real-ref-alleles --out {params.out}'

## Index bcf
rule Sample_IndexBcf:
    input:
        bcf = expand("{DataOut}/{{sample}}_samp_thinned_flipped.vcf.gz", DataOut=DATAOUT)
    output:
        "{DataOut}/{{sample}}_samp_thinned_flipped.vcf.gz.csi".format(DataOut=DATAOUT)
    shell:
        'bcftools index -f {input.bcf}'

## ---- Principal Compoent analysis ----
##  Project ADNI onto a PCA using the 1000 Genomes dataset to identify population outliers

##  Extract a pruned dataset from 1000 genomes using the same pruning SNPs from Sample
## align 1000 genomes to fasta refrence
rule Reference_flip:
    input:
        bim = "data/1000genomes_allChr.bim",
        bed = "data/1000genomes_allChr.bed",
        fam = "data/1000genomes_allChr.fam",
        fasta = "data/hg19.fa"
    output:
        temp(expand("data/1000genomes_allChr_flipped.{ext}", ext=BPLINK))
    shell:
        "/Users/sheaandrews/Programs/flippyr/flippyr.py -p {input.fasta} {input.bim}"

rule Reference_ChromPosRefAlt:
    input:
        bim = "data/1000genomes_allChr_flipped.bim".format(DataOut=DATAOUT),
    output:
        bim = temp("{DataOut}/1000genomes_allChr_flipped.bim".format(DataOut=DATAOUT)),
        snplist = temp("{DataOut}/Reference_snplist".format(DataOut=DATAOUT))
    shell:
        "Rscript scripts/bim_ChromPosRefAlt.R {input.bim} {output.bim} {output.snplist}"

rule Reference_prune:
    input:
        expand("data/1000genomes_allChr_flipped.{ext}", ext=BPLINK, DataOut=DATAOUT),
        "{DataOut}/1000genomes_allChr_flipped.bim".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_thinned_snplist".format(DataOut=DATAOUT),
    output:
        expand("{DataOut}/{{sample}}_1kg_thinned_flipped.{ext}", ext=BPLINK, DataOut=DATAOUT)
    params:
        indat_plink = "data/1000genomes_allChr_flipped".format(DataOut=DATAOUT),
        bim = "{DataOut}/1000genomes_allChr_flipped.bim".format(DataOut=DATAOUT),
        indat_prune_in = "{DataOut}/{{sample}}_thinned_snplist".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_1kg_thinned_flipped".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat_plink} --bim {params.bim} --filter-founders --extract {params.indat_prune_in} --keep-allele-order --make-bed --out {params.out}'

## Recode 1kg to vcf
rule Reference_Plink2Bcf:
    input:
        expand("{DataOut}/{{sample}}_1kg_thinned_flipped.{ext}", ext=BPLINK, DataOut=DATAOUT)
    output:
        expand("{DataOut}/{{sample}}_1kg_thinned_flipped.vcf.gz", DataOut=DATAOUT)
    params:
        indat = "{DataOut}/{{sample}}_1kg_thinned_flipped".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_1kg_thinned_flipped".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat} --recode vcf bgz --keep-allele-order --real-ref-alleles --out {params.out}'

## Index bcf
rule Reference_IndexBcf:
    input:
        bcf = expand("{DataOut}/{{sample}}_1kg_thinned_flipped.vcf.gz", DataOut=DATAOUT)
    output:
        "{DataOut}/{{sample}}_1kg_thinned_flipped.vcf.gz.csi".format(DataOut=DATAOUT)
    shell:
        'bcftools index -f {input.bcf}'

## Merge ref and sample
rule Merge_RefenceSample:
    input:
        bcf_1kg = expand("{DataOut}/{{sample}}_1kg_thinned_flipped.vcf.gz", DataOut=DATAOUT),
        csi_1kg = expand("{DataOut}/{{sample}}_1kg_thinned_flipped.vcf.gz.csi", DataOut=DATAOUT),
        bcf_samp = expand("{DataOut}/{{sample}}_samp_thinned_flipped.vcf.gz", DataOut=DATAOUT),
        csi_samp = expand("{DataOut}/{{sample}}_samp_thinned_flipped.vcf.gz.csi", DataOut=DATAOUT)
    output:
        out = "{DataOut}/{{sample}}_1kg_merged.vcf".format(DataOut=DATAOUT)
    shell:
        'bcftools merge -m none {input.bcf_1kg} {input.bcf_samp} | \
        bcftools view -g \^miss \
        -Ov -o {output.out}'

## recode merged sample to plink
rule Plink_RefenceSample:
    input:
        vcf = "{DataOut}/{{sample}}_1kg_merged.vcf".format(DataOut=DATAOUT),
    output:
        expand("{DataOut}/{{sample}}_1kg_merged.{ext}", ext = BPLINK, DataOut=DATAOUT)
    params:
        out = "{DataOut}/{{sample}}_1kg_merged".format(DataOut=DATAOUT)
    shell:
        'plink --vcf {input.vcf} --const-fid --make-bed --out {params.out}'

rule fix_fam:
    input:
        fam = "{DataOut}/{{sample}}_1kg_merged.fam".format(DataOut=DATAOUT)
    output:
        out = "{DataOut}/{{sample}}_1kg_merged_fixed.fam".format(DataOut=DATAOUT)
    shell:
        'scripts/fix_fam.py {input.fam} {output.out}'

## PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        expand("{DataOut}/{{sample}}_1kg_merged.{ext}", ext = BPLINK, DataOut=DATAOUT),
        "{DataOut}/{{sample}}_1kg_merged_fixed.fam".format(DataOut=DATAOUT),
        "data/1000genomes_pops.txt",
        "data/pops.txt"
    output:
        expand("{DataOut}/{{sample}}_1kg_merged.{ext}", ext = ['eigenval', 'eigenvec'], DataOut=DATAOUT)
    params:
        indat_plink = "{DataOut}/{{sample}}_1kg_merged".format(DataOut=DATAOUT),
        indat_fam = "{DataOut}/{{sample}}_1kg_merged_fixed.fam".format(DataOut=DATAOUT),
        indat_pop = "data/1000genomes_pops.txt",
        indat_clust = "data/pops.txt",
        out = "{DataOut}/{{sample}}_1kg_merged".format(DataOut=DATAOUT),
    shell:
        'plink --bfile {params.indat_plink} --fam {params.indat_fam} --pca 10 --within {params.indat_pop} --pca-clusters {params.indat_clust} --out {params.out}'

## Rscript to identify EUR population outliers
rule ExcludePopulationOutliers:
    input:
        indat_eigenval = "{DataOut}/{{sample}}_1kg_merged.eigenval".format(DataOut=DATAOUT),
        indat_eigenvec = "{DataOut}/{{sample}}_1kg_merged.eigenvec".format(DataOut=DATAOUT),
        indat_fam = "{DataOut}/{{sample}}_samp_thinned.fam".format(DataOut=DATAOUT),
        indat_1kgped = "data/20130606_g1k.ped"
    output:
        out = "{DataOut}/{{sample}}_exclude.pca".format(DataOut=DATAOUT)
    shell:
        'Rscript scripts/PCA_QC.R {input.indat_eigenvec} {input.indat_1kgped} {input.indat_fam} {input.indat_eigenval} {output.out} '

## Run PCA to for population stratification
rule PopulationStratification:
    input:
        expand("{DataOut}/{{sample}}_samp_thinned.{ext}", ext = BPLINK, DataOut=DATAOUT),
        "{DataOut}/{{sample}}_exclude.pca".format(DataOut=DATAOUT)
    output:
        expand("{DataOut}/{{sample}}_filtered_PCA.{ext}", ext = ['eigenval', 'eigenvec'], DataOut=DATAOUT),
    params:
        indat = "{DataOut}/{{sample}}_samp_thinned".format(DataOut=DATAOUT),
        exclude = "{DataOut}/{{sample}}_exclude.pca".format(DataOut=DATAOUT),
        out = "{DataOut}/{{sample}}_filtered_PCA".format(DataOut=DATAOUT)
    shell:
        'plink --bfile {params.indat}  --remove {params.exclude} --pca 10 --out {params.out}'

rule SampleExclusion:
    input:
        SampCallRate = "{DataOut}/{{sample}}_callRate.irem".format(DataOut=DATAOUT),
        het = "{DataOut}/{{sample}}_exclude.heterozigosity".format(DataOut=DATAOUT),
        sex = "{DataOut}/{{sample}}_exclude.sexcheck".format(DataOut=DATAOUT),
        pca = "{DataOut}/{{sample}}_exclude.pca".format(DataOut=DATAOUT),
        relat = "{DataOut}/{{sample}}_exclude.relatedness".format(DataOut=DATAOUT)
    output:
        out = "{DataOut}/{{sample}}_exclude.samples".format(DataOut=DATAOUT)
    shell:
        'Rscript scripts/sample_QC.R {input.SampCallRate} {input.het} {input.sex} {input.pca} {input.relat} {output.out}'

print(FAMILY
)

rule GWAS_QC_Report:
    input:
        "scripts/GWAS_QC.Rmd",
        "{DataOut}/{{sample}}_SnpQc.hwe".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_SnpQc.frq".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_SnpQc.frqx".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_callRate.imiss".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_IBDQC.genome".format(DataOut=DATAOUT),
        "{DataIn}/{{sample}}.fam".format(DataIn=DATAIN),
        "{DataOut}/{{sample}}_1kg_merged.eigenval".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_1kg_merged.eigenvec".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_samp_thinned.fam".format(DataOut=DATAOUT),
        "data/20130606_g1k.ped",
        "{DataOut}/{{sample}}_filtered_PCA.eigenval".format(DataOut=DATAOUT),
        "{DataOut}/{{sample}}_filtered_PCA.eigenvec".format(DataOut=DATAOUT)
    output:
        "stats/{sample}_GWAS_QC.html"
    params:
        rwd = RWD,
        Sample = "{sample}",
        Family = FAMILY,
        output_dir = "stats",
        Path_SexFile = expand("{DataOut}/{SampleSex}_SexQC.sexcheck", SampleSex=SAMPLESEX, DataOut=DATAOUT),
        Path_hwe = "{DataOut}/{{sample}}_SnpQc.hwe".format(DataOut=DATAOUT),
        Path_frq = "{DataOut}/{{sample}}_SnpQc.frq".format(DataOut=DATAOUT),
        Path_frqx = "{DataOut}/{{sample}}_SnpQc.frqx".format(DataOut=DATAOUT),
        Path_imiss = "{DataOut}/{{sample}}_callRate.imiss".format(DataOut=DATAOUT),
        Path_HetFile = "{DataOut}/{{sample}}_HetQC.het".format(DataOut=DATAOUT),
        Path_GenomeFile = "{DataOut}/{{sample}}_IBDQC.genome".format(DataOut=DATAOUT),
        Path_eigenval = "{DataOut}/{{sample}}_1kg_merged.eigenval".format(DataOut=DATAOUT),
        Path_eigenvec = "{DataOut}/{{sample}}_1kg_merged.eigenvec".format(DataOut=DATAOUT),
        Path_TargetPops = "{DataOut}/{{sample}}_samp_thinned.fam".format(DataOut=DATAOUT),
        PATH_BasePops = "data/20130606_g1k.ped",
        Path_PopStrat_eigenval = "{DataOut}/{{sample}}_filtered_PCA.eigenval".format(DataOut=DATAOUT),
        Path_PopStrat_eigenvec = "{DataOut}/{{sample}}_filtered_PCA.eigenvec".format(DataOut=DATAOUT),
        out = "{sample}_GWAS_QC.html",
    shell:
        "R -e 'rmarkdown::render("
        """"{input[0]}", output_file = "{params.out}", output_dir = "{params.output_dir}", \
        params = list(rwd = "{params.rwd}", Sample = "{params.Sample}", \
        Path_SexFile = "{params.Path_SexFile}", \
        Path_hwe = "{params.Path_hwe}", Path_frq = "{params.Path_frq}", Path_frqx = "{params.Path_frqx}", Path_imiss = "{params.Path_imiss}", \
        Path_HetFile = "{params.Path_HetFile}", \
        Path_GenomeFile = "{params.Path_GenomeFile}", Family = {params.Family}, \
        Path_eigenval = "{params.Path_eigenval}", Path_eigenvec = "{params.Path_eigenvec}", Path_TargetPops = "{params.Path_TargetPops}", PATH_BasePops = "{params.PATH_BasePops}", \
        Path_PopStrat_eigenval = "{params.Path_PopStrat_eigenval}", Path_PopStrat_eigenvec = "{params.Path_PopStrat_eigenvec}"))' --slave
        """
