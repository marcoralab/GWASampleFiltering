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

# QC Steps:
QC_snp = True
QC_callRate = True
QC_sex = True


def decorate(text):
    return expand("{DataOut}/{sample}_" + text,
                  sample=SAMPLE, DataOut=DATAOUT)


rule all:
    input:
        expand("stats/{sample}_GWAS_QC.html", sample=SAMPLE)
#        expand("{DataOut}/{sample}_filtered_PCA.{ext}",
#               ext=['eigenval', 'eigenvec'],
#               sample=SAMPLE, DataOut=DATAOUT),
#        decorate("exclude.samples"),
#        # snp_qc:
#        decorate("SnpQc.hwe"),  # hwe file
#        decorate("SnpQc.frq"),  # plink frequency file
#        decorate("SnpQc.frqx"),  #
#        decorate("exclude.sexcheck")  # exclusion file for sex discordance

start = expand("{DataIn}/{{sample}}.{ext}", ext=BPLINK, DataIn=DATAIN)

# ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input: start
    output:
        temp(expand("{{DataOut}}/{{sample}}_SnpQc.{ext}", ext=BPLINK)),
        "{DataOut}/{sample}_SnpQc.hwe",
        "{DataOut}/{sample}_SnpQc.frq",
        "{DataOut}/{sample}_SnpQc.frqx",
    params:
        indat = '{DataIn}/{sample}',
        out = "{DataOut}/{sample}_SnpQc",
    shell:
        """
plink --bfile {params.indat} --freq --out {params.out}
plink --bfile {params.indat} --freqx --out {params.out}
plink --bfile {params.indat} --geno 0.05 --maf 0.01 \
--hardy --make-bed --out {params.out}"""

# ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input: rules.snp_qc.output if QC_snp else start
    output:
        temp(expand("{{DataOut}}/{{sample}}_callRate.{ext}", ext=BPLINK)),
        "{DataOut}/{sample}_callRate.imiss",
        touch("{DataOut}/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out,
        out = "{DataOut}/{sample}_callRate"
    shell:
        """
plink --bfile {params.indat} --mind 0.05 \
--missing --make-bed --out {params.out}"""

# ---- Exclude Samples with discordant sex ----
#  Use ADNI hg18 data, as the liftover removed the x chromsome data


rule sexcheck_QC:
    input:
        expand("{DataIn}/{SampleSex}.{ext}", ext=BPLINK)
    output:
        "{DataOut}/{Sample}_SexQC.sexcheck"
    params:
        indat = '{DataIn}/{{SampleSex}}'.format(DataIn=DATAIN),
        out = "{DataOut}/{Sample}_SexQC",
    shell:
        'plink --bfile {params.indat} --check-sex --out {params.out}'

rule sex_sample_fail:
    input:
        rules.sexcheck_QC.output
    output:
        "{DataOut}/{sample}_exclude.sexcheck",
    shell:
        'Rscript scripts/sexcheck_QC.R {input} {output}'

if QC_callRate:
    sexcheck_in_plink = rules.sample_callRate.output[0]
elif QC_snp:
    sexcheck_in_plink = rules.snp_qc.output
else:
    sexcheck_in_plink = start

rule sex_exclude_failed:
    input:
        plink = rules.sample_callRate.output[0],
        indat_exclude = rules.sex_sample_fail.output
    output:
        temp(expand("{{DataOut}}/{{sample}}_SexExclude.{ext}", ext=BPLINK))
    params:
        indat_plink = "{DataOut}/{sample}_callRate",
        out = "{DataOut}/{sample}_SexExclude"
    shell:
        """
plink --bfile {params.indat_plink} \
--remove {input.indat_exclude} \
--make-bed --out {params.out}"""

# ---- Prune SNPs, autosome only ----
#  Pruned SNP list is used for IBD, PCA and heterozigosity calculations

if QC_sex:
    QCd_plink = rules.sex_exclude_failed.output
elif QC_callRate:
    QCd_plink = rules.sample_callRate.output[0]
elif QC_snp:
    QCd_plink = rules.snp_qc.output
else:
    QCd_plink = start

rule PruneDupvar_snps:
    input: QCd_plink
    output:
        expand("{{DataOut}}/{{sample}}_thinned.{ext}",
               ext=['prune.in', 'prune.out']),
        "{DataOut}/{sample}_thinned.dupvar.delete"
    params:
        indat = "{DataOut}/{sample}_SexExclude",
        dupvar = "{DataOut}/{sample}_thinned.dupvar",
        out = "{DataOut}/{sample}_thinned"
    shell:
        """
plink --bfile {params.indat} --autosome --indep 50 5 1.5 \
--list-duplicate-vars --out {params.out}
Rscript scripts/DuplicateVars.R {params.dupvar}"""

# Prune sample dataset
rule sample_prune:
    input:
        expand("{{DataOut}}/{{sample}}_SexExclude.{ext}", ext=BPLINK),
        prune = "{DataOut}/{sample}_thinned.prune.in",
        dupvar = "{DataOut}/{sample}_thinned.dupvar.delete"
    output:
        temp(expand("{{DataOut}}/{{sample}}_samp_thinned.{ext}", ext=BPLINK))
    params:
        indat_plink = "{DataOut}/{sample}_SexExclude",
        out = "{DataOut}/{sample}_samp_thinned"
    shell:
        """
plink --bfile {params.indat_plink} --extract {input.prune} \
--exclude {input.dupvar} --make-bed --out {params.out}"""

# ---- Exclude Samples with interealtedness ----
rule relatedness_QC:
    input:
        expand("{{DataOut}}/{{sample}}_samp_thinned.{ext}", ext=BPLINK),
    output:
        "{DataOut}/{sample}_IBDQC.genome"
    params:
        indat_plink = "{DataOut}/{sample}_samp_thinned",
        out = "{DataOut}/{sample}_IBDQC"
    shell:
        """
plink --bfile {params.indat_plink} --genome --min 0.05 --out {params.out}"""

rule relatedness_sample_fail:
    input:
        genome = "{DataOut}/{sample}_IBDQC.genome",
        fam = "{DataOut}/{sample}_SnpQc.fam"
    params:
        Family = FAMILY
    output:
        out = "{DataOut}/{sample}_exclude.relatedness"
    shell:
        """
Rscript scripts/relatedness_QC.R {input.genome} {input.fam} \
{params.Family} {output.out}"""

rule relatedness_exclude_failed:
    input:
        exclude = "{DataOut}/{sample}_exclude.relatedness",
        plink = expand("{{DataOut}}/{{sample}}_samp_thinned.{ext}", ext=BPLINK)
    output:
        temp(expand("{{DataOut}}/{{sample}}_RelatednessExclude.{ext}",
                    ext=BPLINK))
    params:
        indat_plink = "{DataOut}/{sample}_samp_thinned",
        out = "{DataOut}/{sample}_RelatednessExclude"
    shell:
        """
plink --bfile {params.indat_plink} --remove {input.exclude} \
--make-bed --out {params.out}"""

# ---- Exclude Samples with outlying heterozigosity ----
rule heterozigosity_QC:
    input:
        expand("{{DataOut}}/{{sample}}_RelatednessExclude.{ext}", ext=BPLINK)
    output: "{DataOut}/{sample}_HetQC.het"
    params:
        indat_plink = "{DataOut}/{sample}_RelatednessExclude",
        out = "{DataOut}/{sample}_HetQC"
    shell:
        'plink --bfile {params.indat_plink} --het --out {params.out}'

rule heterozigosity_sample_fail:
    input: rules.heterozigosity_QC.output
    output: "{DataOut}/{sample}_exclude.heterozigosity"
    shell: 'Rscript scripts/heterozygosity_QC.R {input} {output}'

rule heterozigosity_exclude_failed:
    input:
        exclude = "{DataOut}/{sample}_exclude.heterozigosity",
        plink = expand("{{DataOut}}/{{sample}}_RelatednessExclude.{ext}",
                       ext=BPLINK)
    output:
        temp(expand("{{DataOut}}/{{sample}}_HetExclude.{ext}", ext=BPLINK))
    params:
        indat_plink = "{DataOut}/{sample}_RelatednessExclude",
        out = "{DataOut}/{sample}_HetExclude"
    shell:
        """
plink --bfile {params.indat_plink} --remove {input.exclude} \
--make-bed --out {params.out}"""

# align sample to fasta refrence
rule Sample_Flip:
    input:
        bim = "{DataOut}/{sample}_HetExclude.bim",
        bed = "{DataOut}/{sample}_HetExclude.bed",
        fam = "{DataOut}/{sample}_HetExclude.fam",
        fasta = "data/hg19fa"
    output:
        temp(expand("{{DataOut}}/{{sample}}_HetExclude_flipped.{ext}",
                    ext=BPLINK))
    shell:
        """
/Users/sheaandrews/Programs/flippyr/flippyr.py -p {input.fasta} {input.bim}"""

rule Sample_ChromPosRefAlt:
    input: "{DataOut}/{sample}_HetExclude_flipped.bim"
    output:
        bim = temp("{DataOut}/{sample}_HetExclude_flipped_ChromPos.bim"),
        snplist = temp("{DataOut}/{sample}_thinned_snplist")
    shell:
        """
Rscript scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"""

# Recode sample plink file to vcf
rule Sample_Plink2Bcf:
    input:
        expand("{{DataOut}}/{{sample}}_HetExclude_flipped.{ext}", ext=BPLINK),
        "{DataOut}/{sample}_HetExclude_flipped_ChromPos.bim"
    output:
        "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz"
    params:
        indat = "{DataOut}/{sample}_HetExclude_flipped",
        bim = "{DataOut}/{sample}_HetExclude_flipped_ChromPos.bim",
        out = "{DataOut}/{sample}_samp_thinned_flipped"
    shell:
        """
plink --bfile {params.indat} --bim {params.bim} --recode vcf bgz \
--keep-allele-order --real-ref-alleles --out {params.out}"""

# Index bcf
rule Sample_IndexBcf:
    input:
        bcf = "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz"
    output:
        "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz.csi"
    shell:
        'bcftools index -f {input.bcf}'

# ---- Principal Compoent analysis ----
#  Project ADNI onto a PCA using the 1000 Genomes dataset to identify
#    population outliers

#  Extract a pruned dataset from 1000 genomes using the same pruning SNPs
#    from Sample
# align 1000 genomes to fasta refrence
rule Reference_flip:
    input:
        bim = "data/1000genomes_allChr.bim",
        bed = "data/1000genomes_allChr.bed",
        fam = "data/1000genomes_allChr.fam",
        fasta = "data/hg19.fa"
    output:
        temp(expand("data/1000genomes_allChr_flipped.{ext}", ext=BPLINK))
    shell:
        "flippyr -p {input.fasta} {input.bim}"

rule Reference_ChromPosRefAlt:
    input: "data/1000genomes_allChr_flipped.bim"
    output:
        bim = temp("{DataOut}/1000genomes_allChr_flipped.bim"),
        snplist = temp("{DataOut}/Reference_snplist")
    shell:
        """
Rscript scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}
"""

rule Reference_prune:
    input:
        plink = expand("data/1000genomes_allChr_flipped.{ext}", ext=BPLINK),
        bim = "{DataOut}/1000genomes_allChr_flipped.bim",
        prune = "{DataOut}/{sample}_thinned_snplist"
    output:
        expand("{{DataOut}}/{{sample}}_1kg_thinned_flipped.{ext}", ext=BPLINK)
    params:
        indat_plink = "data/1000genomes_allChr_flipped",
        out = "{DataOut}/{sample}_1kg_thinned_flipped"
    shell:
        """
plink --bfile {params.indat_plink} --bim {input.bim} --filter-founders \
--extract {input.prune} --keep-allele-order --make-bed --out {params.out}"""


# Recode 1kg to vcf
rule Reference_Plink2Bcf:
    input:
        expand("{{DataOut}}/{{sample}}_1kg_thinned_flipped.{ext}", ext=BPLINK)
    output: "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz"
    params:
        indat = "{DataOut}/{sample}_1kg_thinned_flipped",
        out = "{DataOut}/{sample}_1kg_thinned_flipped"
    shell:
        """
plink --bfile {params.indat} --recode vcf bgz \
--real-ref-alleles --out {params.out}"""

# Index bcf
rule Reference_IndexBcf:
    input:
        bcf = "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz"
    output:
        "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz.csi"
    shell:
        'bcftools index -f {input.bcf}'

# Merge ref and sample
rule Merge_RefenceSample:
    input:
        bcf_1kg = "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz",
        csi_1kg = "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz.csi",
        bcf_samp = "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz",
        csi_samp = "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz.csi",
    output:
        out = "{DataOut}/{sample}_1kg_merged.vcf"
    shell:
        r"""
bcftools merge -m none {input.bcf_1kg} {input.bcf_samp} | \
bcftools view -g ^miss -Ov -o {output.out}"""

# recode merged sample to plink
rule Plink_RefenceSample:
    input:
        vcf = "{DataOut}/{sample}_1kg_merged.vcf"
    output:
        expand("{{DataOut}}/{{sample}}_1kg_merged.{ext}", ext=BPLINK)
    params:
        out = "{DataOut}/{sample}_1kg_merged"
    shell:
        'plink --vcf {input.vcf} --const-fid --make-bed --out {params.out}'

rule fix_fam:
    input:
        fam = "{DataOut}/{sample}_1kg_merged.fam"
    output:
        out = "{DataOut}/{sample}_1kg_merged_fixed.fam"
    shell:
        'scripts/fix_fam.py {input.fam} {output.out}'

# PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        plink = expand("{{DataOut}}/{{sample}}_1kg_merged.{ext}", ext=BPLINK),
        fam = "{DataOut}/{sample}_1kg_merged_fixed.fam",
        pop = "data/1000genomes_pops.txt",
        clust = "data/pops.txt"
    output:
        expand("{{DataOut}}/{{sample}}_1kg_merged.{ext}",
               ext=['eigenval', 'eigenvec'])
    params:
        indat_plink = "{DataOut}/{sample}_1kg_merged",
        out = "{DataOut}/{sample}_1kg_merged"
    shell:
        """
plink --bfile {params.indat_plink} --fam {input.fam} --pca 10 \
--within {input.pop} --pca-clusters {input.clust} --out {params.out}
"""

# Rscript to identify EUR population outliers
rule ExcludePopulationOutliers:
    input:
        indat_eigenval = "{DataOut}/{sample}_1kg_merged.eigenval",
        indat_eigenvec = "{DataOut}/{sample}_1kg_merged.eigenvec",
        indat_fam = "{DataOut}/{sample}_samp_thinned.fam",
        indat_1kgped = "data/20130606_g1k.ped"
    output:
        out = "{DataOut}/{sample}_exclude.pca"
    shell:
        """
Rscript scripts/PCA_QC.R {input.indat_eigenvec} {input.indat_1kgped} \
{input.indat_fam} {input.indat_eigenval} {output.out}
"""

# Run PCA to for population stratification
rule PopulationStratification:
    input:
        plink = expand("{{DataOut}}/{{sample}}_samp_thinned.{ext}",
                       ext=BPLINK),
        exclude = "{DataOut}/{sample}_exclude.pca"
    output:
        expand("{{DataOut}}/{{sample}}_filtered_PCA.{ext}",
               ext=['eigenval', 'eigenvec'])
    params:
        indat = "{DataOut}/{sample}_samp_thinned",
        out = "{DataOut}/{sample}_filtered_PCA"
    shell:
        """
plink --bfile {params.indat} --remove {input.exclude} --pca 10 \
--out {params.out}
        """

rule SampleExclusion:
    input:
        SampCallRate = "{DataOut}/{sample}_callRate.irem",
        het = "{DataOut}/{sample}_exclude.heterozigosity",
        sex = "{DataOut}/{sample}_exclude.sexcheck",
        pca = "{DataOut}/{sample}_exclude.pca",
        relat = "{DataOut}/{sample}_exclude.relatedness"
    output:
        out = "{DataOut}/{sample}_exclude.samples"
    shell:
        """
Rscript scripts/sample_QC.R {input.SampCallRate} {input.het} \
{input.sex} {input.pca} {input.relat} {output.out}'"""

print(FAMILY)

rule GWAS_QC_Report:
    input:
        script = "scripts/GWAS_QC.Rmd",
        SexFile = decorate("SexQC.sexcheck"),
        hwe = decorate("SnpQc.hwe"),
        frq = decorate("SnpQc.frq"),
        frqx = decorate("SnpQc.frqx"),
        imiss = decorate("callRate.imiss"),
        HetFile = decorate("HetQC.het"),
        GenomeFile = decorate("IBDQC.genome"),
        eigenval = decorate("1kg_merged.eigenval"),
        eigenvec = decorate("1kg_merged.eigenvec"),
        TargetPops = decorate("samp_thinned.fam"),
        BasePops = "data/20130606_g1k.ped",
        PopStrat_eigenval = decorate("filtered_PCA.eigenval"),
        PopStrat_eigenvec = decorate("filtered_PCA.eigenvec")
    output:
        "stats/{sample}_GWAS_QC.html"
    params:
        rwd = RWD,
        Family = FAMILY,
        output_dir = "stats"
    shell:
        """
R -e 'rmarkdown::render("{input.script}", \
output_file = "{output}", output_dir = "{params.output_dir}", \
params = list(rwd = "{params.rwd}", Sample = "{wildcards.sample}", \
Path_SexFile = "{input.SexFile}", Path_hwe = "{input.hwe}", \
Path_frq = "{input.frq}", Path_frqx = "{input.frqx}", \
Path_imiss = "{input.imiss}", Path_HetFile = "{input.HetFile}", \
Path_GenomeFile = "{input.GenomeFile}", Family = {params.Family}, \
Path_eigenval = "{input.eigenval}", \
Path_eigenvec = "{input.eigenvec}", \
Path_TargetPops = "{input.TargetPops}", \
PATH_BasePops = "{input.BasePops}", \
Path_PopStrat_eigenval = "{input.PopStrat_eigenval}", \
Path_PopStrat_eigenvec = "{input.PopStrat_eigenvec}"))' --slave
"""
