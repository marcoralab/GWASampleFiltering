'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

from scripts.parse_config_GWASampleFiltering import parser

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

configfile: "config/config.yaml"

shell.executable("/bin/bash")

BPLINK = ["bed", "bim", "fam"]

start, SAMPLE, DATAOUT = parser(config)

localrules: download_tg_fa, download_tg_ped, download_tg_chrom


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

ereftype = detect_ref_type(config['extra_ref']) if extraref else 'none'


# ---- Principal Compoent analysis ----
#  Project ADNI onto a PCA using the 1000 Genomes dataset to identify
#    population outliers

#  Extract a pruned dataset from 1000 genomes using the same pruning SNPs
#    from Sample
# align 1000 genomes to fasta refrence


if config['mirror'] == 'ncbi':
    tgbase = "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/"
elif config['mirror'] == 'ebi':
    tgbase = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"

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
    conda: 'envs/bcftools.yaml'
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
        input: config['custom_ref']['custom_pops']
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
            vcf = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz" if REF == '1kG' else config['custom_ref']['file'],
            tbi = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi" if REF == '1kG' else config['custom_ref']['file'] + '.tbi'
        output: temp("reference/{refname}.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        conda: "envs/bcftools.yaml"
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
        conda: "envs/bcftools.yaml"
        shell:
            '''
bcftools concat {input.vcfs} -Oz -o {output.vcf} --threads 2
bcftools index -ft {output.vcf}
'''
elif creftype == 'vcf':
    rule Reference_prep:
        input:
            vcf = config['custom_ref']['file'],
            tbi = config['custom_ref']['file'] + '.tbi'
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "envs/bcftools.yaml"
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
            bim = config['custom_ref']['file'] + '.bim',
            bed = config['custom_ref']['file'] + '.bed',
            fam = config['custom_ref']['file'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand("reference/{{refname}}_{{gbuild}}_flipped.{ext}", ext=BPLINK))
        conda: "envs/flippyr.yaml"
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
        conda: "envs/r.yaml"
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
        conda: "envs/plink.yaml"
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
        conda: "envs/bcftools.yaml"
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
            vcf = config['extra_ref']['file'],
            tbi = config['extra_ref']['file'] + '.tbi'
        output: temp("{dataout}/extraref.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        conda: "envs/bcftools.yaml"
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
            vcf = config['extra_ref']['file'],
            tbi = config['extra_ref']['file'] + '.tbi'
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "envs/bcftools.yaml"
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
            bim = config['extra_ref']['file'] + '.bim',
            bed = config['extra_ref']['file'] + '.bed',
            fam = config['extra_ref']['file'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand("{{dataout}}/extraref_{{gbuild}}_flipped.{ext}", ext=BPLINK, dataout = DATAOUT))
        conda: "envs/flippyr.yaml"
        shell:"flippyr -p {input.fasta} -o {DATAOUT}/extraref_{wildcards.gbuild}_flipped {input.bim}"

    rule Ref_ChromPosRefAlt_extra:
        input:
            flipped = "{dataout}/extraref_{gbuild}_flipped.bim"
        output:
            bim = temp("{dataout}/extraref_{gbuild}_flipped_ChromPos.bim"),
            snplist = temp("{dataout}/extraref_{gbuild}_flipped_snplist")
        conda: "envs/r.yaml"
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
        conda: "envs/plink.yaml"
        shell:
            '''
plink --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

    # Index bcf
    rule Ref_IndexVcf_extra:
        input: "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        conda: "envs/bcftools.yaml"
        shell: 'bcftools index -f {input}'

    rule Reference_prep_extra:
        input:
            vcf = rules.Ref_Plink2Vcf_extra.output,
            csi = rules.Ref_IndexVcf_extra.output
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        conda: "envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''

union_extraref = (extraref
                  and ("overlap_extra" in config)
                  and (config['overlap_extra'] == 'union'))



rule get_panelvars:
    input: "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz"
    output: '{dataout}/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps'
    conda: "envs/bcftools.yaml"
    shell: "bcftools query -f '%ID\n' {input} > {output}"

rule get_extravars:
    input: '{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz'
    output: '{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.snps'
    conda: "envs/bcftools.yaml"
    shell: "bcftools query -f '%ID\n' {input} > {output}"

rule union_panelvars:
    input:
        ref = expand(
            '{dataout}/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
            gbuild=BUILD, miss=config['QC']['GenoMiss'], refname=REF,
            dataout=DATAOUT),
        eref = expand("{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.snps",
                      gbuild=BUILD, miss=config['QC']['GenoMiss'],
                      dataout=DATAOUT)
    output: '{dataout}/panelvars_all.snp' if union_extraref else "do.not"
    shell: 'cat {input} | sort | uniq > {output}'

rule intersection_panelvars:
    input:
        ref = expand(
            '{dataout}/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
            gbuild=BUILD, miss=config['QC']['GenoMiss'], refname=REF,
            dataout=DATAOUT),
        eref = expand("{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.snps",
                      gbuild=BUILD, miss=config['QC']['GenoMiss'],
                      dataout=DATAOUT)
    output: '{dataout}/panelvars_all.snp' if not union_extraref else "do.not"
    conda: "envs/miller.yaml"
    shell:
        '''
mlr --tsv --implicit-csv-header --headerless-csv-output \
  join -j 1 -f {input[0]} {input[1]} > {output}
'''
