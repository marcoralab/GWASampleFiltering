'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''
# import ipdb; ipdb.set_trace()
from workflow.scripts.parse_config_GWASampleFiltering import parser
import socket
import sys
import getpass
import warnings

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from urllib.request import urlopen
from urllib.error import URLError

try:
    response = urlopen('https://www.google.com/', timeout=10)
    iconnect = True
except urllib.error.URLError as ex:
    iconnect = False


class dummyprovider:
    def remote(string_, allow_redirects = "foo"):
        return string_

if iconnect and not ('nointernet' in config and config['nointernet']):
    FTP = FTPRemoteProvider()
    HTTP = HTTPRemoteProvider()
else:
    FTP = dummyprovider
    HTTP = dummyprovider

# FTP = FTPRemoteProvider() if iconnect else dummyprovider
# HTTP = HTTPRemoteProvider() if iconnect else dummyprovider

isMinerva = "hpc.mssm.edu" in socket.getfqdn()

configfile: "config/config.yaml"
do_sexqc = config['do_sexqc']

shell.executable("/bin/bash")

BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
# import ipdb; ipdb.set_trace()
start, SAMPLE, DATAOUT = parser(config)

FAMILY = 'T' if config['family'] else 'F'
istg = config['isTG'] if 'isTG' in config else False

# QC Steps:
QC_snp = True
QC_callRate = True
union_panel_TF = True

localrules: all, download_tg_fa, download_tg_ped, download_tg_chrom


def flatten(nested):
    flat = []
    for el in nested:
        if not isinstance(el, list):
            flat.append(el)
        else:
            flat += flatten(el)
    return flat


def detect_ref_type(reffile):
    if '.vcf' in reffile or '.bcf' in reffile:
        if '{chrom}' in reffile:
            return "vcfchr"
        else:
            return "vcf"
    else:
        return "plink"


if ("custom_ref" not in config) or (config['custom_ref'] is False):
    REF = '1kG'
else:
    REF = config['custom_ref_name']
    creftype = detect_ref_type(config['custom_ref'])

extraref = False
if ("extra_ref" in config) and (config['extra_ref'] is not False):
    extraref = True
    ereftype = detect_ref_type(config['extra_ref'])
else:
    ereftype = 'none'

if "pca_sd" in config:
    pca_sd = int(config["pca_sd"])
else:
    pca_sd = 6

outs = {
    "report": expand("{dataout}/stats/{sample}_GWAS_QC.html", sample=SAMPLE, dataout = DATAOUT),
    "exclusions": expand("{dataout}/{sample}_exclude.samples", sample=SAMPLE, dataout = DATAOUT),
    "filtered": expand("{dataout}/{sample}_Excluded.{ext}",
                       sample=SAMPLE, ext=BPLINK, dataout = DATAOUT)}

outputs = [outs[x] for x in config["outputs"]]
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
    conda: "workflow/envs/plink.yaml"
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
    conda: "workflow/envs/plink.yaml"
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
    conda: "workflow/envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} \
  --check-sex --aec --out {params.out}
'''

rule sex_sample_fail:
    input: rules.sexcheck_QC.output
    output: "{dataout}/{sample}_exclude.sexcheck"
    conda: 'workflow/envs/r.yaml'
    script: 'workflow/scripts/sexcheck_QC.R'

if QC_callRate:
    sexcheck_in_plink = rules.sample_callRate.output[0]
    sexcheck_in_plink_stem = rules.sample_callRate.params.out
elif QC_snp:
    sexcheck_in_plink = expand("{{dataout}}/{{sample}}_SnpQc.{ext}", ext=BPLINK)
    sexcheck_in_plink_stem = rules.snp_qc.params.out
else:
    sexcheck_in_plink = start['files']
    sexcheck_in_plink_stem = start['stem']

# ---- Principal Compoent analysis ----
#  Project ADNI onto a PCA using the 1000 Genomes dataset to identify
#    population outliers

#  Extract a pruned dataset from 1000 genomes using the same pruning SNPs
#    from Sample
# align 1000 genomes to fasta refrence


if config['mirror'] == 'ncbi':
    tgbase = "ftps://ftp-trace.ncbi.nih.gov/1000genomes/ftp/"
elif config['mirror'] == 'ebi':
    tgbase = "ftps://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"

tgped = tgbase + "technical/working/20130606_sample_info/20130606_g1k.ped"

if config['genome_build'] in ['hg19', 'hg37', 'GRCh37', 'grch37', 'GRCH37']:
    BUILD = 'hg19'
    tgurl = (tgbase + "release/20130502/ALL.chr{chrom}." +
        "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
    tgfa = tgbase + "technical/reference/human_g1k_v37.fasta"
elif config['genome_build'] in ['hg38', 'GRCh38', 'grch38', 'GRCH38']:
    BUILD = 'GRCh38'
    tgurl = (tgbase +
        'data_collections/1000_genomes_project/release/20181203_biallelic_SNV/' +
        'ALL.chr{chrom}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz')
    tgfa = 'technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'

# import ipdb; ipdb.set_trace()
rule download_tg_chrom:
    input:
        HTTP.remote(tgurl),
        HTTP.remote(tgurl + ".tbi"),
    output:
        temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz"),
        temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi")
    shell: "cp {input[0]} {output[0]}; cp {input[1]} {output[1]}"


rule download_tg_fa:
    input:
        HTTP.remote(tgfa + ".gz") if BUILD == 'hg19' else HTTP.remote(tgfa)
    output:
        "reference/human_g1k_{gbuild}.fasta",
        "reference/human_g1k_{gbuild}.fasta.fai"
    conda: 'workflow/envs/bcftools.yaml'
    shell:
        '''
if [[ "{input[0]}" == *.gz ]]; then
  zcat {input[0]} > {output[0]} && rstatus=0 || rstatus=$?; true
  if [ $rstatus -ne 2 && $rstatus -ne 0 ]; then
    exit $rstatus
  fi
else
  cp {input[0]} {output[0]}
fi
samtools faidx {output[0]}
'''

rule download_tg_ped:
    input:
        HTTP.remote(tgped),
    output:
        "reference/20130606_g1k.ped",
    shell: "cp {input} {output}"

tgped = "reference/20130606_g1k.ped"


if REF == '1kG':
    """
    Selects everyone who is unrelated or only has third degree relatives in
    thousand genomes.
    """
    rule Reference_foundersonly:
        input: tgped
        output:
            "reference/20130606_g1k.founders"
        shell:
            r'''
awk -F "\t" '!($12 != 0 || $10 != 0 || $9 != 0 || $3 != 0 || $4 != 0) {{print $2}}' \
{input} > {output}
'''

    rule makeTGpops:
        input: tgped
        output:
            "reference/1kG_pops.txt",
            "reference/1kG_pops_unique.txt"
        shell:
            '''
awk 'BEGIN {{print "FID","IID","Population"}} NR>1 {{print $1,$2,$7}}' \
{input} > {output[0]}
cut -f7 {input} | sed 1d | sort | uniq > {output[1]}
'''
else:
    rule make_custom_pops:
        input: config['custom_pops']
        output:
            "reference/{refname}_pops.txt",
            "reference/{refname}_pops_unique.txt"
        shell:
            '''
cp {input} {output[0]}
awk 'NR > 1 {{print $3}}' {input} | sort | uniq > {output[1]}
'''


if REF == '1kG' or creftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz" if REF == '1kG' else config['custom_ref'],
            tbi = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi" if REF == '1kG' else config['custom_ref'] + '.tbi'
        output: temp("reference/{refname}.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
'''

    rule Reference_cat:
        input:
            vcfs = expand("reference/{{refname}}.{{gbuild}}.chr{chrom}.maxmiss{{miss}}.vcf.gz",
                          chrom = list(range(1, 23)))
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
bcftools concat {input.vcfs} -Oz -o {output.vcf} --threads 2
bcftools index -ft {output.vcf}
'''
elif creftype == 'vcf':
    rule Reference_prep:
        input:
            vcf = config['custom_ref'],
            tbi = config['custom_ref'] + '.tbi'
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''
else: #PLINK fileset of all chromosomes
    # align custom ref to fasta refrence
    rule Ref_Flip:
        input:
            bim = config['custom_ref'] + '.bim',
            bed = config['custom_ref'] + '.bed',
            fam = config['custom_ref'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand("reference/{{refname}}_{{gbuild}}_flipped.{ext}", ext=BPLINK))
        conda: "workflow/envs/flippyr.yaml"
        shell:
            '''
flippyr -p {input.fasta} \
  -o reference/{wildcards.refname}_{wildcards.gbuild}_flipped {input.bim}
'''

    rule Ref_ChromPosRefAlt:
        input:
            flipped = "reference/{refname}_{gbuild}_flipped.bim"
        output:
            bim = temp("reference/{refname}_{gbuild}_flipped_ChromPos.bim"),
            snplist = temp("reference/{refname}_{gbuild}_flipped_snplist")
        conda: "workflow/envs/r.yaml"
        shell: "Rscript scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"

    # Recode sample plink file to vcf
    rule Ref_Plink2Vcf:
        input:
            bim = rules.Ref_ChromPosRefAlt.output.bim,
            flipped = rules.Ref_Flip.output
        output:
            temp("reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz")
        params:
            out = "reference/{refname}_{gbuild}_unQC_maxmissUnfilt",
            inp = "reference/{refname}_{gbuild}_flipped"
        conda: "workflow/envs/plink.yaml"
        shell:
            '''
plink --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

    # Index bcf
    rule Ref_IndexVcf:
        input: "reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: "reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        shell: 'bcftools index -f {input}'

    rule Reference_prep:
        input:
            vcf = rules.Ref_Plink2Vcf.output,
            csi = rules.Ref_IndexVcf.output
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        #wildcard_constraints:
            #gbuild = "hg19|GRCh38",
            #refname = "[a-zA-Z0-9-]",
            #miss = "[0-9.]",
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''

if ereftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = config['extra_ref'],
            tbi = config['extra_ref'] + '.tbi'
        output: temp("{dataout}/extraref.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
'''

    rule Reference_cat_extra:
        input:
            vcfs = expand("{dataout}/extraref.{{gbuild}}.chr{chrom}.maxmiss{{miss}}.vcf.gz",
                          chrom = list(range(1, 23)), dataout = DATAOUT)
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        shell:
            '''
bcftools concat {input.vcfs} -Oz -o {output.vcf} --threads 2
bcftools index -ft {output.vcf}
'''
elif ereftype == 'vcf':
    rule Reference_prep_extra:
        input:
            vcf = config['extra_ref'],
            tbi = config['extra_ref'] + '.tbi'
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''
elif ereftype != 'none': #PLINK fileset of all chromosomes
    # align custom ref to fasta refrence
    rule Ref_Flip_extra:
        input:
            bim = config['extra_ref'] + '.bim',
            bed = config['extra_ref'] + '.bed',
            fam = config['extra_ref'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand("{{dataout}}/extraref_{{gbuild}}_flipped.{ext}", ext=BPLINK, dataout = DATAOUT))
        conda: "workflow/envs/flippyr.yaml"
        shell:"flippyr -p {input.fasta} -o {DATAOUT}/extraref_{wildcards.gbuild}_flipped {input.bim}"

    rule Ref_ChromPosRefAlt_extra:
        input:
            flipped = "{dataout}/extraref_{gbuild}_flipped.bim"
        output:
            bim = temp("{dataout}/extraref_{gbuild}_flipped_ChromPos.bim"),
            snplist = temp("{dataout}/extraref_{gbuild}_flipped_snplist")
        conda: "workflow/envs/r.yaml"
        shell: "R scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"

    # Recode sample plink file to vcf
    rule Ref_Plink2Vcf_extra:
        input:
            bim = rules.Ref_ChromPosRefAlt_extra.output.bim,
            flipped = rules.Ref_Flip_extra.output
        output:
            temp("{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz")
        params:
            out = "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt",
            inp = "{dataout}/extraref_{gbuild}_flipped"
        conda: "workflow/envs/plink.yaml"
        shell:
            '''
plink --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

    # Index bcf
    rule Ref_IndexVcf_extra:
        input: "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        conda: "workflow/envs/bcftools.yaml"
        shell: 'bcftools index -f {input}'

    rule Reference_prep_extra:
        input:
            vcf = rules.Ref_Plink2Vcf_extra.output,
            csi = rules.Ref_IndexVcf_extra.output
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "workflow/envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''


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
    conda: "workflow/envs/flippyr.yaml"
    shell: "flippyr -p {input.fasta} -o {params.dataout}/{wildcards.sample} {input.bim}"

rule Sample_ChromPosRefAlt:
    input:
        flipped = "{dataout}/{sample}_flipped.bim"
    output:
        bim = temp("{dataout}/{sample}_flipped_ChromPos.bim"),
        snplist = temp("{dataout}/{sample}_flipped_snplist")
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/bim_ChromPosRefAlt.R"

if config['union_panel'] and extraref:
    if extraref and ("union_extra" in config) and (config['union_extra'] != False):
        erefc = True
    PVARS = DATAOUT + '/panelvars_all.snp'
    extract_sample = "--extract {} ".format(PVARS)
else:
    PVARS = "/dev/urandom"
    extract_sample = ""


if not extraref:
    erefc = False

rule get_panelvars:
    input:
        expand("reference/{{refname}}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
               gbuild=BUILD, miss=config['QC']['GenoMiss'], refname=REF)
    output: 'reference/panelvars_{refname}.snps' if erefc else DATAOUT + '/panelvars_all.snp'
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools query -f '%ID\n' {input} > {output}"

rule merge_panelvars:
    input:
        ref = expand('reference/panelvars_{refname}.snps', refname=REF),
        eref = expand("{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
                      gbuild=BUILD, miss=config['QC']['GenoMiss'], dataout = DATAOUT)
    output: PVARS if erefc else "do.not"
    shell: 'cat {input} | sort | uniq > {output}'
#import ipdb; ipdb.set_trace()

rule PruneDupvar_snps:
    input:
        fileset = rules.Sample_Flip.output,
        bim = rules.Sample_ChromPosRefAlt.output.bim,
        pvars = PVARS
    output:
        "{dataout}/{sample}_nodup.dupvar",
        expand("{{dataout}}/{{sample}}_nodup.{ext}",
               ext=['prune.in', 'prune.out'], dataout = DATAOUT),
    params:
        indat = "{dataout}/{sample}_flipped",
        dupvar = "{dataout}/{sample}_nodup.dupvar",
        out = "{dataout}/{sample}_nodup",
        extract = extract_sample
    conda: "workflow/envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} -bim {input.bim} \
  {params.extract}--autosome --indep 50 5 1.5 \
  --list-duplicate-vars --out {params.out}
'''

rule SelectDupvar_snps:
    input: rules.PruneDupvar_snps.output[0]
    output: "{dataout}/{sample}_nodup.dupvar.delete"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/DuplicateVars.R"

# Prune sample dataset
rule sample_prune:
    input:
        fileset = rules.Sample_Flip.output,
        bim = rules.Sample_ChromPosRefAlt.output.bim,
        prune = "{dataout}/{sample}_nodup.prune.in",
        dupvar = rules.SelectDupvar_snps.output
    output:
        temp(expand("{{dataout}}/{{sample}}_pruned.{ext}",
                    ext=BPLINK, dataout = DATAOUT))
    params:
        indat = "{dataout}/{sample}_flipped",
        out = "{dataout}/{sample}_pruned"
    conda: "workflow/envs/plink.yaml"
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
    conda: "workflow/envs/bcftools.yaml"
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
    conda: "workflow/envs/bcftools.yaml"
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
    conda: "workflow/envs/plink.yaml"
    shell:
        '''
plink --bed {input.bed} --bim {input.bim} --fam {input.fam} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

# Index bcf
rule Sample_IndexBcf:
    input: "{dataout}/{sample}_pruned.vcf.gz"
    output: "{dataout}/{sample}_pruned.vcf.gz.csi"
    conda: "workflow/envs/bcftools.yaml"
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
    conda: "workflow/envs/bcftools.yaml"
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
        expand("{{dataout}}/{{sample}}_{{refname}}_merged.{ext}", ext=BPLINK, dataout = DATAOUT)
    params:
        out = "{dataout}/{sample}_{refname}_merged"
    conda: "workflow/envs/plink.yaml"
    shell: "plink --keep-allele-order --vcf {input.vcf} --const-fid --make-bed --out {params.out}"

rule fix_fam:
    input:
        oldfam = rules.Sample_Plink2Bcf.input.fam,
        newfam = "{dataout}/{sample}_{refname}_merged.fam",
        tgped = tgped
    output: fixed = "{dataout}/{sample}_{refname}_merged_fixed.fam"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/fix_fam.R"

rule merge_pops:
    input:
        main = "reference/{refname}_pops.txt",
        extra = rules.Reference_prune_extra.input.vcf
    output:
        DATAOUT + '/{refname}_allpops.txt',
        DATAOUT + '/{refname}_allpops_unique.txt'
    params:
        extra_ref_code = config['extra_ref_subpop']
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/add_extraref_pops.R"

# PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        plink = expand("{{dataout}}/{{sample}}_{{refname}}_merged.{ext}", ext=BPLINK, dataout = DATAOUT),
        fam = rules.fix_fam.output,
        ref = DATAOUT + '/{refname}_allpops.txt' if extraref else "reference/{refname}_pops.txt",
        clust = DATAOUT + '/{refname}_allpops_unique.txt' if extraref else "reference/{refname}_pops_unique.txt"
    output:
        expand("{{dataout}}/{{sample}}_{{refname}}_merged.{ext}", ext=['eigenval', 'eigenvec'], dataout = DATAOUT)
    params:
        indat_plink = "{dataout}/{sample}_{refname}_merged",
        out = "{dataout}/{sample}_{refname}_merged"
    conda: "workflow/envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat_plink} --fam {input.fam} \
  --pca 10 --within {input.ref} --pca-clusters {input.clust} --out {params.out}
'''

# Rscript to identify population outliers
rule ExcludePopulationOutliers:
    input:
        eigenval = expand("{dataout}/{{sample}}_{refname}_merged.eigenval", refname=REF, dataout = DATAOUT),
        eigenvec = expand("{dataout}/{{sample}}_{refname}_merged.eigenvec", refname=REF, dataout = DATAOUT),
        fam = rules.Sample_Plink2Bcf.input.fam,
        pops = expand(DATAOUT + '/{refname}_allpops.txt' if extraref else "reference/{refname}_pops.txt", refname=REF)
    output:
        excl = "{dataout}/{sample}_exclude.pca",
        rmd = temp("{dataout}/{sample}_pca.Rdata")
    params:
        superpop = config['superpop'],
        extraref = 'none' if not extraref else config['extra_ref_subpop'],
        sd = pca_sd
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/PCA_QC.R"
    # shell:
    #     '''
    #     Rscript scripts/PCA_QC.R -s {wildcards.sample} -p {params.superpop} \
    #     --vec {input.eigenvec} --val {input.eigenval} \
    #     -b {input.tgped} -t {input.fam} -o {output.excl} -R {output.rmd}
    #     '''

# ---- Exclude Samples with interealtedness ----
rule relatedness_sample_prep:
    input: sexcheck_in_plink
    output:
        bed = temp("{dataout}/{sample}_IBDQCfilt.bed"),
        bim = temp("{dataout}/{sample}_IBDQCfilt.bim"),
        fam = temp("{dataout}/{sample}_IBDQCfilt.fam")
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = "{dataout}/{sample}_IBDQCfilt"
    conda: "workflow/envs/plink.yaml"
    shell:
        '''
plink --bfile {params.indat_plink} \
  --geno 0.02 \
  --maf 0.02 \
  --memory 6000 \
  --make-bed --out {params.out}
'''

if config['king']:
    rule relatedness_QC:
        input:
            bed = rules.relatedness_sample_prep.output.bed,
            bim = rules.relatedness_sample_prep.output.bim,
            fam = rules.relatedness_sample_prep.output.fam
        output:
            genome = "{dataout}/{sample}_IBDQC.kingfiles"
        params:
            out = "{dataout}/{sample}_IBDQC",
            dataout = DATAOUT
        conda: "workflow/envs/king.yaml"
        shell:
            '''
king -b {input.bed} --related --degree 3 --prefix {params.out}
if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.kin*")"; then
  find {params.dataout} -name "{wildcards.sample}_IBDQC.kin*" > {output.genome}
fi
'''

    rule king_all:
        input:
            bed = rules.relatedness_sample_prep.output.bed,
            bim = rules.relatedness_sample_prep.output.bim,
            fam = rules.relatedness_sample_prep.output.fam
        output:
            genome = "{dataout}/{sample}_IBDQC.all.kingfiles",
        params:
            out = "{dataout}/{sample}_IBDQC.all",
            dataout = DATAOUT
        conda: "workflow/envs/king.yaml"
        shell:
            '''
king -b {input.bed} --kinship --ibs --prefix {params.out}
if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.all.kin*")"; then
  find {params.dataout} -name "{wildcards.sample}_IBDQC.all.kin*" > {output.genome}
fi
'''
else:
    rule relatedness_QC:
        input: rules.sample_prune.output
        output:
            genome = "{dataout}/{sample}_IBDQC.genome"
        params:
            indat_plink = "{dataout}/{sample}_pruned",
            out = "{dataout}/{sample}_IBDQC"
        conda: "workflow/envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat_plink} --genome --min 0.05 \
  --out {params.out}
'''

rule relatedness_sample_fail:
    input:
        genome = rules.relatedness_QC.output.genome,
        geno_all = rules.king_all.output if config['king'] else "/dev/null",
        fam = sexcheck_in_plink_stem + ".fam"
    params:
        Family = FAMILY,
        king = config['king'],
        threshold = 0.1875,
        geno = rules.relatedness_QC.params.out if config['king'] else rules.relatedness_QC.output.genome
    output:
        out = "{dataout}/{sample}_exclude.relatedness",
        rdat = "{dataout}/{sample}_IBDQC.Rdata"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/relatedness_QC.R"

# ---- Exclude Samples with outlying heterozigosity ----
rule heterozygosity_QC:
    input: rules.sample_prune.output
    output: "{dataout}/{sample}_HetQC.het"
    params:
        indat_plink = rules.sample_prune.params.out,
        out = "{dataout}/{sample}_HetQC"
    conda: "workflow/envs/plink.yaml"
    shell: "plink --keep-allele-order --bfile {params.indat_plink} --het --out {params.out}"

rule heterozygosity_sample_fail:
    input: rules.heterozygosity_QC.output
    output: "{dataout}/{sample}_exclude.heterozigosity"
    conda: "workflow/envs/r.yaml"
    script: 'workflow/scripts/heterozygosity_QC.R'

# Run PCA to control for population stratification
if config['pcair']:
    assert config['king'], "You must use KING for relatedness if using PCAiR!"
    rule ancestryFilt:
        input:
            plink = expand("{{dataout}}/{{sample}}_pruned.{ext}",
                           ext=BPLINK, dataout = DATAOUT),
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            temp(expand("{{dataout}}/{{sample}}_filtered_PCApre.{ext}",
                 ext=BPLINK, dataout = DATAOUT)),
        params:
            indat = "{dataout}/{sample}_pruned",
            plinkout = "{dataout}/{sample}_filtered_PCApre"
        conda: "workflow/envs/plink.yaml"
        shell:
            r'''
plink --keep-allele-order --bfile {params.indat} \
  --remove {input.exclude} \
  --make-bed --out {params.plinkout}
'''

    rule filterKING:
        input:
            king = "{dataout}/{sample}_IBDQC.all.kingfiles",
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            "{dataout}/{sample}_IBDQC.all.popfilt.kingfiles"
        params:
            indat = "{dataout}/{sample}_IBDQC.all",
        conda: "workflow/envs/r.yaml"
        script: "workflow/scripts/filterKing.R"
#         shell:
#             r'''
# Rscript workflow/scripts/filterKing.R {params.indat} {input.exclude}
# if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.all.popfilt.kin*")"; then
#   find {params.dataout} -name "{wildcards.sample}_IBDQC.all.popfilt.kin*" > {output}
# fi
# '''

    rule PCAPartitioning:
        input:
            plink = rules.ancestryFilt.output,
            king = rules.filterKING.output,
            iterative = rules.relatedness_sample_fail.output.out
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCApre.{ext}",ext=['unrel', 'partition.log'], dataout = DATAOUT)
        params:
            stem = rules.ancestryFilt.params.plinkout,
            king = rules.filterKING.params.indat + ".popfilt"
        conda: "workflow/envs/r.yaml"
        script: "workflow/scripts/PartitionPCAiR.R"

    rule stratFrq:
        input:
            plink = rules.ancestryFilt.output,
            unrel = rules.PCAPartitioning.output[0],
        output: "{dataout}/{sample}_filtered_PCAfreq.frqx"
        params:
            indat = rules.ancestryFilt.params.plinkout,
            out = "{dataout}/{sample}_filtered_PCAfreq"
        conda: "workflow/envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} --freqx \
  --within {input.unrel} --keep-cluster-names unrelated \
  --out {params.out}
'''

    rule PopulationStratification:
        input:
            plink = rules.ancestryFilt.output,
            unrel = rules.PCAPartitioning.output[0],
            frq = rules.stratFrq.output
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCA.{ext}", ext=['eigenval', 'eigenvec'], dataout = DATAOUT)
        params:
            indat = rules.ancestryFilt.params.plinkout,
            out = "{dataout}/{sample}_filtered_PCA"
        conda: "workflow/envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} --read-freq {input.frq} --pca 10 \
  --within {input.unrel} --pca-cluster-names unrelated \
  --out {params.out}
'''
else:
    rule PopulationStratification:
        input:
            plink = expand("{{dataout}}/{{sample}}_pruned.{ext}",
                           ext=BPLINK, dataout = DATAOUT),
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCA.{ext}", ext=['eigenval', 'eigenvec'], dataout = DATAOUT)
        params:
            indat = "{dataout}/{sample}_pruned",
            out = "{dataout}/{sample}_filtered_PCA"
        conda: "workflow/envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} --remove {input.exclude}
  --pca 10 --out {params.out}
'''

rule SampleExclusion:
    input:
        SampCallRate = "{dataout}/{sample}_callRate.irem",
        het = "{dataout}/{sample}_exclude.heterozigosity",
        sex = "{dataout}/{sample}_exclude.sexcheck" if do_sexqc else '/dev/null',
        pca = "{dataout}/{sample}_exclude.pca",
        relat = "{dataout}/{sample}_exclude.relatedness"
    output:
        out = "{dataout}/{sample}_exclude.samples",
        out_distinct = "{dataout}/{sample}_exclude.distinct_samples"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/sample_QC.R"

rule Exclude_failed:
    input:
        plink = sexcheck_in_plink,
        indat_exclude = rules.SampleExclusion.output.out_distinct
    output:
        temp(expand("{{dataout}}/{{sample}}_Excluded.{ext}", ext=BPLINK, dataout = DATAOUT)),
        excl = temp('{dataout}/{sample}_exclude.plink')
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = "{dataout}/{sample}_Excluded"
    conda: "workflow/envs/plink.yaml"
    shell:
        '''
cat {input.indat_exclude} | sed '1d' | cut -d' ' -f1,2 > {output.excl}
plink --keep-allele-order --bfile {params.indat_plink} --remove {output.excl} \
--make-bed --out {params.out}
'''


rule GWAS_QC_Report:
    input:
        SexFile = "{dataout}/{sample}_SexQC.sexcheck" if do_sexqc else '/dev/null',
        hwe = "{dataout}/{sample}_SnpQc.hwe",
        frq = "{dataout}/{sample}_SnpQc.frq",
        frqx = "{dataout}/{sample}_SnpQc.frqx",
        imiss = "{dataout}/{sample}_callRate.imiss",
        HetFile = "{dataout}/{sample}_HetQC.het",
        IBD_stats = "{dataout}/{sample}_IBDQC.Rdata",
        PCA_rdat = "{dataout}/{sample}_pca.Rdata",
        PopStrat_eigenval = "{dataout}/{sample}_filtered_PCA.eigenval",
        PopStrat_eigenvec = "{dataout}/{sample}_filtered_PCA.eigenvec",
        partmethod = rules.PCAPartitioning.output[1] if config["pcair"] else "/dev/null"
    output:
         "{dataout}/stats/{sample}_GWAS_QC.html"
    params:
        sample = "{sample}"
        script = "GWAS_QC.Rmd",
        rwd = RWD,
        Family = FAMILY,
        pi_threshold = 0.1875,
        output_dir = "{dataout}/stats",
        idir = "{dataout}/stats/md/{sample}",
        geno_miss = config['QC']['GenoMiss'],
        samp_miss = config['QC']['SampMiss'],
        MAF = config['QC']['MAF'],
        HWE = config['QC']['HWE'],
        superpop = config['superpop'],
        partmethod = rules.PCAPartitioning.output[1] if config["pcair"] else "none"
    conda: "workflow/envs/r.yaml"
    script: "workflow/scripts/RenderMarkdown.R"
#     shell:
#         '''
# R -e 'nm <- sample(c("Shea J. Andrews", "Brian Fulton-Howard"), \
# replace=F); nm <- paste(nm[1], "and", nm[2]); \
# rmarkdown::render("{params.script}", output_dir = "{params.output_dir}", \
# output_file = "{output}", intermediates_dir = "{params.idir}", \
# params = list(rwd = "{params.rwd}", Sample = "{wildcards.sample}", \
# auth = nm, Path_SexFile = "{input.SexFile}", Path_hwe = "{input.hwe}", \
# Path_frq = "{input.frq}", Path_frqx = "{input.frqx}", \
# Path_imiss = "{input.imiss}", Path_HetFile = "{input.HetFile}", \
# pi_threshold = {params.pi_threshold}, Family = {params.Family}, \
# Path_IBD_stats = "{input.IBD_stats}", Path_PCA_rdat = "{input.PCA_rdat}", \
# Path_PopStrat_eigenval = "{input.PopStrat_eigenval}", \
# Path_PopStrat_eigenvec = "{input.PopStrat_eigenvec}", maf = {params.MAF}, \
# hwe = {params.HWE}, missing_geno = {params.geno_miss}, \
# partmethod = "{params.partmethod}", \
# missing_sample = {params.samp_miss}, superpop = "{params.superpop}"))' --slave
# '''
