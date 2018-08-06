'''Snakefile for GWAS Variant and Sample QC Version 0.2'''

from scripts.parse_config import parser

configfile: "config.yaml"

shell.executable("/bin/bash")
shell.prefix("PATH=" + config["anaconda"] + ":$PATH; ")

BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
start, FAMILY, SAMPLE, DATAOUT = parser(config)

# QC Steps:
QC_snp = True
QC_callRate = True

# com = {'flippyr': '/Users/sheaandrews/Programs/flippyr/flippyr.py',
#        'plink': 'plink --keep-allele-order', 'plink2': 'plink',
#        'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R'}
#loads = {'flippyr': '', 'plink': '', 'bcftools': '',  'R': ''}

com = {'flippyr': 'flippyr', 'plink': 'plink --keep-allele-order',
       'plink2': 'plink', 'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R'}
loads = {'flippyr': '', 'plink': 'module load plink/1.90',
         'bcftools': 'module load bcftools/1.7',
         'R': ('module load R/3.4.3 pandoc/2.1.3 udunits/2.2.26; ',
               'RSTUDIO_PANDOC=$(which pandoc)')}


def decorate(text):
    return expand(DATAOUT + "/{sample}_" + text,
                  sample=SAMPLE)


rule all:
    input:
        expand(DATAOUT + "/stats/{sample}_GWAS_QC.html",
               sample=SAMPLE),
        # expand(DATAOUT + "/{sample}_exclude.samples",
        #        sample=SAMPLE),
        expand(DATAOUT + "/{sample}_Excluded.{ext}",
               sample=SAMPLE, ext=BPLINK)


# ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input: start['files']
    output:
        temp(expand(DATAOUT + "/{{sample}}_SnpQc.{ext}", ext=BPLINK)),
        DATAOUT + "/{sample}_SnpQc.hwe",
        DATAOUT + "/{sample}_SnpQc.frq",
        DATAOUT + "/{sample}_SnpQc.frqx",
    params:
        stem = start['stem'],
        out = DATAOUT + "/{sample}_SnpQc",
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
        expand(DATAOUT + "/{{sample}}_callRate.{ext}", ext=BPLINK),
        DATAOUT + "/{sample}_callRate.imiss",
        touch(DATAOUT + "/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out if QC_snp else start['stem'],
        out = DATAOUT + "/{sample}_callRate"
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
        DATAOUT + "/{sample}_SexQC.sexcheck"
    params:
        indat = start['sex_stem'],
        out = DATAOUT + "/{sample}_SexQC",
    shell:
        '''
{loads[plink]}
{com[plink]} --bfile {params.indat} --check-sex --out {params.out}
'''

rule sex_sample_fail:
    input:
        rules.sexcheck_QC.output
    output:
        DATAOUT + "/{sample}_exclude.sexcheck",
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

# ---- Prune SNPs, autosome only ----
#  Pruned SNP list is used for IBD, PCA and heterozigosity calculations

rule PruneDupvar_snps:
    input: sexcheck_in_plink
    output:
        expand(DATAOUT + "/{{sample}}_nodup.{ext}",
               ext=['prune.in', 'prune.out']),
        DATAOUT + "/{sample}_nodup.dupvar.delete"
    params:
        indat = sexcheck_in_plink_stem,
        dupvar = DATAOUT + "/{sample}_nodup.dupvar",
        out = DATAOUT + "/{sample}_nodup"
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
        prune = DATAOUT + "/{sample}_nodup.prune.in",
        dupvar = DATAOUT + "/{sample}_nodup.dupvar.delete"
    output:
        temp(expand(DATAOUT + "/{{sample}}_pruned.{ext}", ext=BPLINK))
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = DATAOUT + "/{sample}_pruned"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --extract {input.prune} \
--exclude {input.dupvar} --make-bed --out {params.out}"""

# ---- Exclude Samples with interealtedness ----
rule relatedness_QC:
    input: rules.sample_prune.output
    output:
        DATAOUT + "/{sample}_IBDQC.genome"
    params:
        indat_plink = DATAOUT + "/{sample}_pruned",
        out = DATAOUT + "/{sample}_IBDQC"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --genome --min 0.05 \
--out {params.out}"""

rule relatedness_sample_fail:
    input:
        genome = DATAOUT + "/{sample}_IBDQC.genome",
        fam = sexcheck_in_plink_stem + ".fam"
    params:
        Family = FAMILY
    output:
        out = DATAOUT + "/{sample}_exclude.relatedness"
    shell:
        """
{loads[R]}; {com[R]}  scripts/relatedness_QC.R {input.genome} {input.fam} \
{params.Family} {output.out}"""

# ---- Exclude Samples with outlying heterozigosity ----
rule heterozygosity_QC:
    input: rules.sample_prune.output
    output: DATAOUT + "/{sample}_HetQC.het"
    params:
        indat_plink = rules.sample_prune.params.out,
        out = DATAOUT + "/{sample}_HetQC"
    shell:
        '''
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --het --out {params.out}'''

rule heterozygosity_sample_fail:
    input: rules.heterozygosity_QC.output
    output: DATAOUT + "/{sample}_exclude.heterozigosity"
    shell: '{loads[R]}; {com[R]}  scripts/heterozygosity_QC.R {input} {output}'

# align sample to fasta refrence
rule Sample_Flip:
    input:
        bim = DATAOUT + "/{sample}_pruned.bim",
        bed = DATAOUT + "/{sample}_pruned.bed",
        fam = DATAOUT + "/{sample}_pruned.fam",
        fasta = "data/hg19.fa"
    output:
        temp(expand(DATAOUT + "/{{sample}}_pruned_flipped.{ext}",
                    ext=BPLINK))
    shell:
        """
{loads[flippyr]}
{com[flippyr]} -p {input.fasta} {input.bim}"""

rule Sample_ChromPosRefAlt:
    input: DATAOUT + "/{sample}_pruned_flipped.bim"
    output:
        bim = temp(DATAOUT + "/{sample}_flipped_ChromPos.bim"),
        snplist = temp(DATAOUT + "/{sample}_pruned_snplist")
    shell:
        """
{loads[R]}
{com[R]} scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"""

# Recode sample plink file to vcf
rule Sample_Plink2Bcf:
    input:
        expand(DATAOUT + "/{{sample}}_pruned_flipped.{ext}", ext=BPLINK),
        DATAOUT + "/{sample}_flipped_ChromPos.bim"
    output:
        DATAOUT + "/{sample}_pruned_flipped.vcf.gz"
    params:
        indat = DATAOUT + "/{sample}_pruned_flipped",
        bim = DATAOUT + "/{sample}_flipped_ChromPos.bim",
        out = DATAOUT + "/{sample}_pruned_flipped"
    shell:
        """
{loads[plink]}
{com[plink2]} --bfile {params.indat} --bim {params.bim} --recode vcf bgz \
--real-ref-alleles --out {params.out}"""

# Index bcf
rule Sample_IndexBcf:
    input:
        bcf = DATAOUT + "/{sample}_pruned_flipped.vcf.gz"
    output:
        DATAOUT + "/{sample}_pruned_flipped.vcf.gz.csi"
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
        expand("data/1000genomes_allChr_flipped.{ext}", ext=BPLINK)
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
        bim = "data/1000genomes_allChr_flipped_CPRA.bim",
        snplist = "data/Reference_snplist"
    shell:
        """
{loads[R]}
{com[R]} scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}
"""

rule Reference_prune:
    input:
        plink = expand("data/1000genomes_allChr_flipped.{ext}", ext=BPLINK),
        bim = "data/1000genomes_allChr_flipped_CPRA.bim",
        prune = DATAOUT + "/{sample}_pruned_snplist"
    output:
        expand(DATAOUT + "/{{sample}}_1kgpruned_flipped.{ext}", ext=BPLINK)
    params:
        indat_plink = "data/1000genomes_allChr_flipped",
        out = DATAOUT + "/{sample}_1kgpruned_flipped"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --bim {input.bim} --filter-founders \
--extract {input.prune} --make-bed --out {params.out}"""


# Recode 1kg to vcf
rule Reference_Plink2Bcf:
    input:
        expand(DATAOUT + "/{{sample}}_1kgpruned_flipped.{ext}", ext=BPLINK)
    output: DATAOUT + "/{sample}_1kgpruned_flipped.vcf.gz"
    params:
        indat = DATAOUT + "/{sample}_1kgpruned_flipped",
        out = DATAOUT + "/{sample}_1kgpruned_flipped"
    shell:
        """
{loads[plink]}
{com[plink2]} --bfile {params.indat} --recode vcf bgz \
--real-ref-alleles --out {params.out}"""

# Index bcf
rule Reference_IndexBcf:
    input:
        bcf = DATAOUT + "/{sample}_1kgpruned_flipped.vcf.gz"
    output:
        DATAOUT + "/{sample}_1kgpruned_flipped.vcf.gz.csi"
    shell:
        '{loads[bcftools]}; {com[bcftools]} index -f {input.bcf}'

# Merge ref and sample
rule Merge_RefenceSample:
    input:
        bcf_1kg = DATAOUT + "/{sample}_1kgpruned_flipped.vcf.gz",
        csi_1kg = DATAOUT + "/{sample}_1kgpruned_flipped.vcf.gz.csi",
        bcf_samp = DATAOUT + "/{sample}_pruned_flipped.vcf.gz",
        csi_samp = DATAOUT + "/{sample}_pruned_flipped.vcf.gz.csi",
    output:
        out = DATAOUT + "/{sample}_1kg_merged.vcf"
    shell:
        r"""
{loads[bcftools]}
{com[bcftools]} merge -m none {input.bcf_1kg} {input.bcf_samp} | \
{com[bcftools]} view -g ^miss -Ov -o {output.out}"""

# recode merged sample to plink
rule Plink_RefenceSample:
    input:
        vcf = DATAOUT + "/{sample}_1kg_merged.vcf"
    output:
        expand(DATAOUT + "/{{sample}}_1kg_merged.{ext}", ext=BPLINK)
    params:
        out = DATAOUT + "/{sample}_1kg_merged"
    shell:
        '''
{loads[plink]}
{com[plink]} --vcf {input.vcf} --const-fid --make-bed --out {params.out}'''

rule fix_fam:
    input:
        fam = DATAOUT + "/{sample}_1kg_merged.fam"
    output:
        out = DATAOUT + "/{sample}_1kg_merged_fixed.fam"
    shell:
        'scripts/fix_fam.py {input.fam} {output.out}'

# PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        plink = expand(DATAOUT + "/{{sample}}_1kg_merged.{ext}", ext=BPLINK),
        fam = DATAOUT + "/{sample}_1kg_merged_fixed.fam",
        pop = "data/1000genomes_pops.txt",
        clust = "data/pops.txt"
    output:
        expand(DATAOUT + "/{{sample}}_1kg_merged.{ext}",
               ext=['eigenval', 'eigenvec'])
    params:
        indat_plink = DATAOUT + "/{sample}_1kg_merged",
        out = DATAOUT + "/{sample}_1kg_merged"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --fam {input.fam} --pca 10 \
--within {input.pop} --pca-clusters {input.clust} --out {params.out}
"""

# Rscript to identify EUR population outliers
rule ExcludePopulationOutliers:
    input:
        indat_eigenval = DATAOUT + "/{sample}_1kg_merged.eigenval",
        indat_eigenvec = DATAOUT + "/{sample}_1kg_merged.eigenvec",
        indat_fam = DATAOUT + "/{sample}_pruned.fam",
        indat_1kgped = "data/20130606_g1k.ped"
    output:
        out = DATAOUT + "/{sample}_exclude.pca"
    shell:
        """
{loads[R]}
{com[R]} scripts/PCA_QC.R {input.indat_eigenvec} {input.indat_1kgped} \
{input.indat_fam} {input.indat_eigenval} {output.out}
"""

# Run PCA to for population stratification
rule PopulationStratification:
    input:
        plink = expand(DATAOUT + "/{{sample}}_pruned.{ext}",
                       ext=BPLINK),
        exclude = DATAOUT + "/{sample}_exclude.pca"
    output:
        expand(DATAOUT + "/{{sample}}_filtered_PCA.{ext}",
               ext=['eigenval', 'eigenvec'])
    params:
        indat = DATAOUT + "/{sample}_pruned",
        out = DATAOUT + "/{sample}_filtered_PCA"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat} --remove {input.exclude} --pca 10 \
--out {params.out}
        """

rule SampleExclusion:
    input:
        SampCallRate = DATAOUT + "/{sample}_callRate.irem",
        het = DATAOUT + "/{sample}_exclude.heterozigosity",
        sex = DATAOUT + "/{sample}_exclude.sexcheck",
        pca = DATAOUT + "/{sample}_exclude.pca",
        relat = DATAOUT + "/{sample}_exclude.relatedness"
    output:
        out = DATAOUT + "/{sample}_exclude.samples",
        out_distinct = DATAOUT + "/{sample}_exclude.distinct_samples"
    shell:
        """
{loads[R]}
{com[R]} scripts/sample_QC.R {input.SampCallRate} {input.het} \
{input.sex} {input.pca} {input.relat} {output.out} {output.out_distinct}
"""

rule Exclude_failed:
    input:
        plink = sexcheck_in_plink,
        indat_exclude = rules.SampleExclusion.output.out_distinct
    output:
        temp(expand(DATAOUT + "/{{sample}}_Excluded.{ext}", ext=BPLINK)),
        excl = temp(DATAOUT + '/{sample}_exclude.plink')
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = DATAOUT + "/{sample}_Excluded"
    shell:
        """
cat {input.indat_exclude} | sed '1d' | cut -d' ' -f1,2 > {output.excl}
{loads[plink]}
{com[plink]} --bfile {params.indat_plink} --remove {output.excl} \
--make-bed --out {params.out}"""


def decorate2(text):
    return DATAOUT + "/{sample}_" + text


rule GWAS_QC_Report:
    input:
        script = "scripts/GWAS_QC.Rmd",
        SexFile = decorate2("SexQC.sexcheck"),
        hwe = decorate2("SnpQc.hwe"),
        frq = decorate2("SnpQc.frq"),
        frqx = decorate2("SnpQc.frqx"),
        imiss = decorate2("callRate.imiss"),
        HetFile = decorate2("HetQC.het"),
        #GenomeFile = decorate2("IBDQC.genome"),
        eigenval = decorate2("1kg_merged.eigenval"),
        eigenvec = decorate2("1kg_merged.eigenvec"),
        TargetPops = decorate2("pruned.fam"),
        BasePops = "data/20130606_g1k.ped",
        PopStrat_eigenval = decorate2("filtered_PCA.eigenval"),
        PopStrat_eigenvec = decorate2("filtered_PCA.eigenvec")
    output:
        DATAOUT + "/stats/{sample}_GWAS_QC.html"
    params:
        rwd = RWD,
        Family = FAMILY,
        output_dir = DATAOUT + "/stats"
    shell:
        """
{loads[R]}
{com[R2]} -e 'rmarkdown::render("{input.script}", \
output_file = "{output}", output_dir = "{params.output_dir}", \
params = list(rwd = "{params.rwd}", Sample = "{wildcards.sample}", \
Path_SexFile = "{input.SexFile}", Path_hwe = "{input.hwe}", \
Path_frq = "{input.frq}", Path_frqx = "{input.frqx}", \
Path_imiss = "{input.imiss}", Path_HetFile = "{input.HetFile}", \
Family = {params.Family}, \
Path_eigenval = "{input.eigenval}", \
Path_eigenvec = "{input.eigenvec}", \
Path_TargetPops = "{input.TargetPops}", \
PATH_BasePops = "{input.BasePops}", \
Path_PopStrat_eigenval = "{input.PopStrat_eigenval}", \
Path_PopStrat_eigenvec = "{input.PopStrat_eigenvec}"))' --slave
"""

#Path_GenomeFile = "{input.GenomeFile}", 
