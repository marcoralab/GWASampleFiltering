'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

from scripts.parse_config_GWASampleFiltering import parser

configfile: "config/config.yaml"
do_sexqc = config['do_sexqc']

shell.executable("/bin/bash")

BPLINK = ["bed", "bim", "fam"]

start, SAMPLE, DATAOUT = parser(config)

istg = config['isTG'] if 'isTG' in config else False

# QC Steps:
QC_snp = True
QC_callRate = True


def flatten(nested):
    flat = []
    for el in nested:
        if not isinstance(el, list):
            flat.append(el)
        else:
            flat += flatten(el)
    return flat


default_ref = (("custom_ref" not in config)
               or (config['custom_ref'] is False)
               or (config['custom_ref']['file'] is False)
               or (config['custom_ref']['name'] is False))

if default_ref:
    REF = '1kG'
else:
    REF = config['custom_ref']['name']
    creftype = detect_ref_type(config['custom_ref']['file'])

extraref = (("extra_ref" in config)
            and (config['extra_ref'] is not False)
            and (config['extra_ref']['file'] is not False)
            and (config['extra_ref']['name'] is not False))

if "pca_sd" in config:
    pca_sd = int(config["pca_sd"])
else:
    pca_sd = 6

outs = {
    "exclude": expand("{dataout}/{sample}_exclude.pca",
                      sample=SAMPLE, dataout=DATAOUT),
    "repdata": expand("{dataout}/{sample}_pca.Rdata",
                      sample=SAMPLE, dataout=DATAOUT)
    }

outputs = [outs[x] for x in ['exclude', 'repdata']]
outputs = flatten(outputs)

rule all:
    input:
        outputs

# ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input: start['files']
    output:
        temp(expand("{{dataout}}/{{sample}}_SnpQc.{ext}", ext=BPLINK)),
        "{dataout}/{sample}_SnpQc.hwe",
        "{dataout}/{sample}_SnpQc.frq",
        "{dataout}/{sample}_SnpQc.frqx",
    params:
        stem = start['stem'],
        out = "{dataout}/{sample}_SnpQc",
        miss = config['QC']['GenoMiss'],
        MAF = config['QC']['MAF'],
        HWE = config['QC']['HWE']
    conda: "envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.stem} --freq --out {params.out}
plink --keep-allele-order --bfile {params.stem} --freqx --out {params.out}
plink --keep-allele-order --bfile {params.stem} --geno {params.miss} \
--maf {params.MAF} --hardy --hwe {params.HWE} --make-bed --out {params.out}
'''

# ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input: expand("{{dataout}}/{{sample}}_SnpQc.{ext}", ext=BPLINK) if QC_snp else start['files']
    output:
        expand("{{dataout}}/{{sample}}_callRate.{ext}", ext=BPLINK, dataout = DATAOUT),
        "{dataout}/{sample}_callRate.imiss",
        touch("{dataout}/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out if QC_snp else start['stem'],
        miss = config['QC']['SampMiss'],
        out = "{dataout}/{sample}_callRate"
    conda: "envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} --mind {params.miss} \
  --missing --make-bed --out {params.out}
'''

# ---- Exclude Samples with discordant sex ----
#  Use ADNI hg18 data, as the liftover removed the x chromsome data


rule sexcheck_QC:
    input: start['sex']
    output:
        "{dataout}/{sample}_SexQC.sexcheck"
    params:
        indat = start['sex_stem'],
        out = "{dataout}/{sample}_SexQC",
    conda: "envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} \
  --check-sex --aec --out {params.out}
'''

rule sex_sample_fail:
    input: rules.sexcheck_QC.output
    output: "{dataout}/{sample}_exclude.sexcheck"
    conda: 'envs/r.yaml'
    script: 'scripts/sexcheck_QC.R'

if QC_callRate:
    sexcheck_in_plink = rules.sample_callRate.output[0]
    sexcheck_in_plink_stem = rules.sample_callRate.params.out
elif QC_snp:
    sexcheck_in_plink = expand("{{dataout}}/{{sample}}_SnpQc.{ext}", ext=BPLINK)
    sexcheck_in_plink_stem = rules.snp_qc.params.out
else:
    sexcheck_in_plink = start['files']
    sexcheck_in_plink_stem = start['stem']

# ---- Prune SNPs, autosome only ----
#  Pruned SNP list is used for IBD, PCA and heterozigosity calculations

# align sample to fasta refrence
rule Sample_Flip:
    input:
        bim = sexcheck_in_plink_stem + '.bim',
        bed = sexcheck_in_plink_stem + '.bed',
        fam = sexcheck_in_plink_stem + '.fam',
        fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
    output:
        multiext("{dataout}/{sample}_flipped", ".bim", ".bed", ".fam")
    params:
        dataout = DATAOUT
    conda: "envs/flippyr.yaml"
    shell: "flippyr -p {input.fasta} -o {params.dataout}/{wildcards.sample} {input.bim}"

rule Sample_ChromPosRefAlt:
    input:
        flipped = "{dataout}/{sample}_flipped.bim"
    output:
        bim = temp("{dataout}/{sample}_flipped_ChromPos.bim"),
        snplist = temp("{dataout}/{sample}_flipped_snplist")
    conda: "envs/r.yaml"
    script: "scripts/bim_ChromPosRefAlt.R"

p_intersect = ((overlap_panel in config)
               and (config['overlap_panel'] == 'intersection'))

panel_variants = '{dataout}/panelvars_all.snp' if p_intersect else '/dev/null'
extract_sample = '--extract {} '.format(panel_variants) if p_intersect else ''

rule PruneDupvar_snps:
    input:
        fileset = rules.Sample_Flip.output,
        bim = rules.Sample_ChromPosRefAlt.output.bim,
        pvars = panel_variants
    output:
        "{dataout}/{sample}_nodup.dupvar",
        expand("{{dataout}}/{{sample}}_nodup.{ext}",
               ext=['prune.in', 'prune.out'], dataout=DATAOUT),
    params:
        indat = "{dataout}/{sample}_flipped",
        dupvar = "{dataout}/{sample}_nodup.dupvar",
        out = "{dataout}/{sample}_nodup",
        extract = extract_sample
    conda: "envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} -bim {input.bim} \
  {params.extract}--autosome --indep 50 5 1.5 \
  --list-duplicate-vars --out {params.out}
'''

rule SelectDupvar_snps:
    input: rules.PruneDupvar_snps.output[0]
    output: "{dataout}/{sample}_nodup.dupvar.delete"
    conda: "envs/r.yaml"
    script: "scripts/DuplicateVars.R"

# Prune sample dataset
rule sample_prune:
    input:
        fileset = rules.Sample_Flip.output,
        bim = rules.Sample_ChromPosRefAlt.output.bim,
        prune = "{dataout}/{sample}_nodup.prune.in",
        dupvar = rules.SelectDupvar_snps.output
    output:
        temp(expand("{{dataout}}/{{sample}}_pruned.{ext}", ext=BPLINK))
    params:
        indat = "{dataout}/{sample}_flipped",
        out = "{dataout}/{sample}_pruned"
    conda: "envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} --bim {input.bim} \
  --extract {input.prune} --exclude {input.dupvar} \
  --make-bed --out {params.out}
'''

rule sample_make_prunelist:
  input: "{dataout}/{sample}_pruned.bim"
  output: "{dataout}/{sample}_pruned.snplist"
  shell: "cut -f2 {input} > {output}"

rule Reference_prune:
    input:
        vcf = expand("reference/{{refname}}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
                     gbuild=BUILD, miss=config['QC']['GenoMiss']),
        prune = "{dataout}/{sample}_pruned.snplist",
        founders = "reference/20130606_g1k.founders" if REF == '1kG' else '/dev/urandom'
    output:
        vcf = temp("{dataout}/{sample}_{refname}pruned.vcf.gz"),
        tbi = temp("{dataout}/{sample}_{refname}pruned.vcf.gz.tbi")
    params:
        founders = "-S reference/20130606_g1k.founders " if REF == '1kG' else ''
    conda: "envs/bcftools.yaml"
    shell:
        '''
bcftools view -i 'ID=@{input.prune}' {params.founders}\
  -Oz -o {output.vcf} --force-samples {input.vcf} --threads 4
bcftools index -ft {output.vcf}
'''

rule Reference_prune_extra:
    input:
        vcf = expand("{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
                     gbuild=BUILD, miss=config['QC']['GenoMiss'], dataout = DATAOUT),
        prune = "{dataout}/{sample}_pruned.snplist"
    output:
        vcf = temp("{dataout}/eref.{sample}pruned.vcf.gz"),
        tbi = temp("{dataout}/eref.{sample}pruned.vcf.gz.tbi")
    conda: "envs/bcftools.yaml"
    shell:
        '''
bcftools view -i 'ID=@{input.prune}' \
  -Oz -o {output.vcf} --force-samples {input.vcf} --threads 4
bcftools index -ft {output.vcf}
'''

# allow for tg sample:
rule tgfam:
    input: "{dataout}/{sample}_pruned.fam"
    output: "{dataout}/{sample}_pruned_tg.fam"
    shell: '''awk '$1 = "1000g___"$1 {{print}}' {input} > {output}'''

# Recode sample plink file to vcf
rule Sample_Plink2Bcf:
    input:
        bed = "{dataout}/{sample}_pruned.bed",
        bim = "{dataout}/{sample}_pruned.bim",
        fam = rules.tgfam.output if istg else rules.tgfam.input
    output: "{dataout}/{sample}_pruned.vcf.gz"
    params:
        out = "{dataout}/{sample}_pruned"
    conda: "envs/plink.yaml"
    shell:
        '''
plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

# Index bcf
rule Sample_IndexBcf:
    input: "{dataout}/{sample}_pruned.vcf.gz"
    output: "{dataout}/{sample}_pruned.vcf.gz.csi"
    conda: "envs/bcftools.yaml"
    shell: 'bcftools index -f {input}'

# Merge ref and sample
rule Merge_RefenceSample:
    input:
        bcf_ref = "{dataout}/{sample}_{refname}pruned.vcf.gz",
        tbi_ref = "{dataout}/{sample}_{refname}pruned.vcf.gz.tbi",
        bcf_samp = "{dataout}/{sample}_pruned.vcf.gz",
        csi_samp = "{dataout}/{sample}_pruned.vcf.gz.csi",
        bcf_ext = "{dataout}/eref.{sample}pruned.vcf.gz" if extraref else '/dev/null',
        tbi_ext = "{dataout}/eref.{sample}pruned.vcf.gz.tbi" if extraref else '/dev/null'
    params:
        miss = config['QC']['GenoMiss'],
        extra = "{dataout}/eref.{sample}pruned.vcf.gz" if extraref else ''
    output:
        out = "{dataout}/{sample}_{refname}_merged.vcf"
    conda: "envs/bcftools.yaml"
    shell:
        r'''
bcftools merge -m none --threads 2 \
  {input.bcf_ref} {input.bcf_samp} {params.extra} | \
  bcftools view  -i 'F_MISSING <= {params.miss}' -Ov -o {output.out} --threads 2
'''

# recode merged sample to plink
rule Plink_RefenceSample:
    input:
        vcf = "{dataout}/{sample}_{refname}_merged.vcf"
    output:
        expand("{{dataout}}/{{sample}}_{{refname}}_merged.{ext}", ext=BPLINK)
    params:
        out = "{dataout}/{sample}_{refname}_merged"
    conda: "envs/plink.yaml"
    shell: "plink --keep-allele-order --vcf {input.vcf} --const-fid --make-bed --out {params.out}"

rule fix_fam:
    input:
        oldfam = rules.Sample_Plink2Bcf.input.fam,
        newfam = "{dataout}/{sample}_{refname}_merged.fam",
        tgped = tgped
    output: fixed = "{dataout}/{sample}_{refname}_merged_fixed.fam"
    conda: "envs/r.yaml"
    script: "scripts/fix_fam.R"

rule merge_pops:
    input:
        main = "reference/{refname}_pops.txt",
        extra = rules.Reference_prune_extra.input.vcf
    output:
        DATAOUT + '/{refname}_allpops.txt',
        DATAOUT + '/{refname}_allpops_unique.txt'
    params:
        extra_ref_code = config['extra_ref_subpop']
    conda: "envs/r.yaml"
    script: "scripts/add_extraref_pops.R"

# PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        plink = expand("{{dataout}}/{{sample}}_{{refname}}_merged.{ext}", ext=BPLINK),
        fam = rules.fix_fam.output,
        ref = DATAOUT + '/{refname}_allpops.txt' if extraref else "reference/{refname}_pops.txt",
        clust = DATAOUT + '/{refname}_allpops_unique.txt' if extraref else "reference/{refname}_pops_unique.txt"
    output:
        expand("{{dataout}}/{{sample}}_{{refname}}_merged.{ext}", ext=['eigenval', 'eigenvec'])
    params:
        indat_plink = "{dataout}/{sample}_{refname}_merged",
        out = "{dataout}/{sample}_{refname}_merged"
    conda: "envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat_plink} --fam {input.fam} \
  --pca 10 --within {input.ref} --pca-clusters {input.clust} --out {params.out}
'''

# Rscript to identify population outliers
rule ExcludePopulationOutliers:
    input:
        eigenval = expand("{{dataout}}/{{sample}}_{refname}_merged.eigenval", refname=REF),
        eigenvec = expand("{{dataout}}/{{sample}}_{refname}_merged.eigenvec", refname=REF),
        fam = rules.Sample_Plink2Bcf.input.fam,
        pops = expand(DATAOUT + '/{refname}_allpops.txt' if extraref else "reference/{refname}_pops.txt", refname=REF)
    output:
        excl = "{dataout}/{sample}_exclude.pca",
        rmd = "{dataout}/{sample}_pca.Rdata"
    params:
        superpop = config['superpop'],
        extraref = 'none' if not extraref else config['extra_ref_subpop'],
        sd = pca_sd
    conda: "envs/r.yaml"
    script: "scripts/PCA_QC.R"
