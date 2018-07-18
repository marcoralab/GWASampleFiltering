'''Snakefile for GWAS Variant and Sample QC Version 0.2'''

from scripts.parse_config import parser
shell.executable("/bin/bash")
configfile: "config.yaml"
BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
start, FAMILY, SAMPLE, DATAOUT = parser(config)

# QC Steps:
QC_snp = True
QC_callRate = True
QC_sex = False

# com = {'flippyr': '/Users/sheaandrews/Programs/flippyr/flippyr.py',
#        'plink': 'plink --keep-allele-order', 'plink2': 'plink',
#        'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R'}
# loads = {'flippyr': '', 'plink': '', 'bcftools': '',  'R': ''}

com = {'flippyr': 'flippyr', 'plink': 'plink --keep-allele-order',
       'plink2': 'plink', 'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R'}
loads = {'flippyr': '', 'plink': 'module load plink/1.90',
         'bcftools': 'module load bcftools/1.7',
         'R': ('module load R/3.4.3 pandoc/2.1.3; ',
               'RSTUDIO_PANDOC=$(which pandoc)')}


def decorate(text):
    return expand("{DataOut}/{sample}_" + text,
                  sample=SAMPLE, DataOut=DATAOUT)


if QC_sex:
    QCd_plink = expand("{DataOut}/{sample}_SexExclude.{ext}",
                       ext=BPLINK, sample=SAMPLE, DataOut=DATAOUT)
elif QC_callRate:
    QCd_plink = expand("{DataOut}/{sample}_callRate.{ext}",
                       ext=BPLINK, sample=SAMPLE, DataOut=DATAOUT)
elif QC_snp:
    QCd_plink = expand("{DataOut}/{sample}_SnpQc.{ext}",
                       ext=BPLINK, sample=SAMPLE, DataOut=DATAOUT)
else:
    QCd_plink = start['files']

rule all:
    input:
        expand("{DataOut}/stats/{sample}_GWAS_QC.html",
               sample=SAMPLE, DataOut=DATAOUT),
        QCd_plink

# ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input: start['files']
    output:
        temp(expand("{{DataOut}}/{{sample}}_SnpQc.{ext}", ext=BPLINK)),
        "{DataOut}/{sample}_SnpQc.hwe",
        "{DataOut}/{sample}_SnpQc.frq",
        "{DataOut}/{sample}_SnpQc.frqx",
    params:
        stem = start['stem'],
        out = "{DataOut}/{sample}_SnpQc",
        miss = config['QC']['GenoMiss'],
        MAF = config['QC']['MAF']
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.stem} --freq --out {params.out}
{com[plink]} --bfile {params.stem} --freqx --out {params.out}
{com[plink]} --bfile {params.stem} --geno {params.miss} \
--maf {params.MAF} --hardy --make-bed --out {params.out}"""

# ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input: rules.snp_qc.output if QC_snp else start['files']
    output:
        expand("{{DataOut}}/{{sample}}_callRate.{ext}", ext=BPLINK),
        "{DataOut}/{sample}_callRate.imiss",
        touch("{DataOut}/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out if QC_snp else start['stem'],
        out = "{DataOut}/{sample}_callRate"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat} --mind 0.05 \
--missing --make-bed --out {params.out}"""

# ---- Exclude Samples with discordant sex ----
#  Use ADNI hg18 data, as the liftover removed the x chromsome data


rule sexcheck_QC:
    input: start['sex']
    output:
        "{DataOut}/{sample}_SexQC.sexcheck"
    params:
        indat = start['sex_stem'],
        out = "{DataOut}/{sample}_SexQC",
    shell:
        '''
{loads[plink]}
{com[plink]} --bfile {params.indat} --check-sex --out {params.out}
'''

rule sex_sample_fail:
    input:
        rules.sexcheck_QC.output
    output:
        "{DataOut}/{sample}_exclude.sexcheck",
    shell:
        '{loads[R]}; {com[R]} scripts/sexcheck_QC.R {input} {output}'

if QC_callRate:
    sexcheck_in_plink = rules.sample_callRate.output[0]
    sexcheck_in_plink_stem = rules.sample_callRate.params.out
elif QC_snp:
    sexcheck_in_plink = rules.snp_qc.output
    sexcheck_in_plink_stem = rules.snp_qc.params.out
else:
    sexcheck_in_plink = start['files']
    sexcheck_in_plink_stem = start['stem']

rule sex_exclude_failed:
    input:
        plink = sexcheck_in_plink,
        indat_exclude = rules.sex_sample_fail.output
    output:
        temp(expand("{{DataOut}}/{{sample}}_SexExclude.{ext}", ext=BPLINK))
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = "{DataOut}/{sample}_SexExclude"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} \
--remove {input.indat_exclude} \
--make-bed --out {params.out}"""

# ---- Prune SNPs, autosome only ----
#  Pruned SNP list is used for IBD, PCA and heterozigosity calculations

rule PruneDupvar_snps:
    input: sexcheck_in_plink
    output:
        expand("{{DataOut}}/{{sample}}_thinned.{ext}",
               ext=['prune.in', 'prune.out']),
        "{DataOut}/{sample}_thinned.dupvar.delete"
    params:
        indat = sexcheck_in_plink_stem,
        dupvar = "{DataOut}/{sample}_thinned.dupvar",
        out = "{DataOut}/{sample}_thinned"
    shell:
        """
{loads[plink]}
{loads[R]}
{com[plink]} --bfile {params.indat} --autosome --indep 50 5 1.5 \
--list-duplicate-vars --out {params.out}
{com[R]}  scripts/DuplicateVars.R {params.dupvar}"""

# Prune sample dataset
rule sample_prune:
    input:
        sexcheck_in_plink,
        prune = "{DataOut}/{sample}_thinned.prune.in",
        dupvar = "{DataOut}/{sample}_thinned.dupvar.delete"
    output:
        temp(expand("{{DataOut}}/{{sample}}_samp_thinned.{ext}", ext=BPLINK))
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = "{DataOut}/{sample}_samp_thinned"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --extract {input.prune} \
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
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --genome --min 0.05 \
--out {params.out}"""

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
{loads[R]}; {com[R]}  scripts/relatedness_QC.R {input.genome} {input.fam} \
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
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --remove {input.exclude} \
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
        '''
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --het --out {params.out}'''

rule heterozigosity_sample_fail:
    input: rules.heterozigosity_QC.output
    output: "{DataOut}/{sample}_exclude.heterozigosity"
    shell: '{loads[R]}; {com[R]}  scripts/heterozygosity_QC.R {input} {output}'

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
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --remove {input.exclude} \
--make-bed --out {params.out}"""

# align sample to fasta refrence
rule Sample_Flip:
    input:
        bim = "{DataOut}/{sample}_HetExclude.bim",
        bed = "{DataOut}/{sample}_HetExclude.bed",
        fam = "{DataOut}/{sample}_HetExclude.fam",
        fasta = "data/hg19.fa"
    output:
        temp(expand("{{DataOut}}/{{sample}}_HetExclude_flipped.{ext}",
                    ext=BPLINK))
    shell:
        """
{loads[flippyr]}
{com[flippyr]} -p {input.fasta} {input.bim}"""

rule Sample_ChromPosRefAlt:
    input: "{DataOut}/{sample}_HetExclude_flipped.bim"
    output:
        bim = temp("{DataOut}/{sample}_HetExclude_flipped_ChromPos.bim"),
        snplist = temp("{DataOut}/{sample}_thinned_snplist")
    shell:
        """
{loads[R]}
{com[R]} scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"""

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
{loads[plink]}
{com[plink2]} --bfile {params.indat} --bim {params.bim} --recode vcf bgz \
--real-ref-alleles --out {params.out}"""

# Index bcf
rule Sample_IndexBcf:
    input:
        bcf = "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz"
    output:
        "{DataOut}/{sample}_samp_thinned_flipped.vcf.gz.csi"
    shell:
        '{loads[bcftools]}; {com[bcftools]} index -f {input.bcf}'

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
        """
{loads[flippyr]}
{com[flippyr]} {input.fasta} {input.bim}
sed 's/--memory 256/--memory 2048/' data/1000genomes_allChr.runPlink > \
data/1000genomes_allChr.moremem.runPlink
bash data/1000genomes_allChr.moremem.runPlink
"""

rule Reference_ChromPosRefAlt:
    input: "data/1000genomes_allChr_flipped.bim"
    output:
        bim = temp("{DataOut}/1000genomes_allChr_flipped.bim"),
        snplist = temp("{DataOut}/Reference_snplist")
    shell:
        """
{loads[R]}
{com[R]} scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}
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
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --bim {input.bim} --filter-founders \
--extract {input.prune} --make-bed --out {params.out}"""


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
{loads[plink]}
{com[plink2]} --bfile {params.indat} --recode vcf bgz \
--real-ref-alleles --out {params.out}"""

# Index bcf
rule Reference_IndexBcf:
    input:
        bcf = "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz"
    output:
        "{DataOut}/{sample}_1kg_thinned_flipped.vcf.gz.csi"
    shell:
        '{loads[bcftools]}; {com[bcftools]} index -f {input.bcf}'

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
{loads[bcftools]}
{com[bcftools]} merge -m none {input.bcf_1kg} {input.bcf_samp} | \
{com[bcftools]} view -g ^miss -Ov -o {output.out}"""

# recode merged sample to plink
rule Plink_RefenceSample:
    input:
        vcf = "{DataOut}/{sample}_1kg_merged.vcf"
    output:
        expand("{{DataOut}}/{{sample}}_1kg_merged.{ext}", ext=BPLINK)
    params:
        out = "{DataOut}/{sample}_1kg_merged"
    shell:
        '''
{loads[plink]}
{com[plink]} --vcf {input.vcf} --const-fid --make-bed --out {params.out}'''

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
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --fam {input.fam} --pca 10 \
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
{loads[R]}
{com[R]} scripts/PCA_QC.R {input.indat_eigenvec} {input.indat_1kgped} \
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
{loads[plink]}
{com[plink]} --bfile {params.indat} --remove {input.exclude} --pca 10 \
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
{loads[R]}
{com[R]}  scripts/sample_QC.R {input.SampCallRate} {input.het} \
{input.sex} {input.pca} {input.relat} {output.out}"""


def decorate2(text):
    return expand("{DataOut}/{{sample}}_" + text, DataOut=DATAOUT)


rule GWAS_QC_Report:
    input:
        script = "scripts/GWAS_QC.Rmd",
        SexFile = decorate2("SexQC.sexcheck"),
        hwe = decorate2("SnpQc.hwe"),
        frq = decorate2("SnpQc.frq"),
        frqx = decorate2("SnpQc.frqx"),
        imiss = decorate2("callRate.imiss"),
        HetFile = decorate2("HetQC.het"),
        GenomeFile = decorate2("IBDQC.genome"),
        eigenval = decorate2("1kg_merged.eigenval"),
        eigenvec = decorate2("1kg_merged.eigenvec"),
        TargetPops = decorate2("samp_thinned.fam"),
        BasePops = "data/20130606_g1k.ped",
        PopStrat_eigenval = decorate2("filtered_PCA.eigenval"),
        PopStrat_eigenvec = decorate2("filtered_PCA.eigenvec")
    output:
        "{DataOut}/stats/{sample}_GWAS_QC.html"
    params:
        rwd = RWD,
        Family = FAMILY,
        output_dir = "{DataOut}/stats"
    shell:
        """
{loads[R]}
{com[R2]} -e 'rmarkdown::render("{input.script}", \
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
