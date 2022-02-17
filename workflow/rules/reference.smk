'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

# from scripts.parse_config_GWASampleFiltering import parser

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

BPLINK = ["bed", "bim", "fam"]

localrules: download_tg_fa, download_tg_ped, download_tg_chrom


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

# ---- Principal Compoent analysis ----
#  Project ADNI onto a PCA using the 1000 Genomes dataset to identify
#    population outliers

#  Extract a pruned dataset from 1000 genomes using the same pruning SNPs
#    from Sample
# align 1000 genomes to fasta refrence


if config['mirror'] == 'ncbi':
    tgbase = "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/"
elif config['mirror'] == 'ebi':
    tgbase = "ftp.1000genomes.ebi.ac.uk/vol1/ftp/"

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
    tgfa = tgbase + 'technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'

# import ipdb; ipdb.set_trace()
rule download_tg_chrom:
    input:
        HTTP.remote(tgurl),
        HTTP.remote(tgurl + ".tbi"),
    output:
        temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz"),
        temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi")
    resources:
        mem_mb = 10000,
        time_min = 30
    shell: "cp {input[0]} {output[0]}; cp {input[1]} {output[1]}"

import time
def http_sleep(url):
    time.sleep(5)
    return HTTP.remote(url)

rule download_tg_fa:
    input:
        HTTP.remote(tgfa + ".gz") if BUILD == 'hg19' else HTTP.remote(tgfa)
    output:
        "reference/human_g1k_{gbuild}.fasta",
        "reference/human_g1k_{gbuild}.fasta.fai"
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    cache: True
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
    cache: True
    resources:
        mem_mb = 10000,
        time_min = 30
    shell: "cp {input} {output}"

tgped = "reference/20130606_g1k.ped"


if REF == '1kG':
    """
    Selects everyone who is unrelated or only has third degree relatives in
    thousand genomes.
    """

    rule Reference_find_founders:
        input: tgped
        output: "reference/20130606_g1k.founders"
        cache: True
        resources:
            mem_mb = 10000,
            time_min = 30
        shell:
            r'''
awk -F "\t" '!($12 != 0 || $10 != 0 || $9 != 0 || $3 != 0 || $4 != 0) {{print $2}}' \
{input} > {output}
'''

    rule Reference_foundersonly:
        input:
            vcf = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz",
            tbi = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi",
            founders = rules.Reference_find_founders.output
        input:
            vcf = temp("reference/1000gFounders.{gbuild}.chr{chrom}.vcf.gz"),
            tbi = temp("reference/1000gFounders.{gbuild}.chr{chrom}.vcf.gz.tbi")
        threads: 4
        resources:
            mem_mb = 4000,
            walltime = '4:00'
            conda: "../envs/bcftools.yaml"
        shell:
            r'''
bcftools view \
  -S {input.founders} \
  -s ^NA20299,NA20314,NA20274,HG01880 \
  -Oz -o {output.vcf} --force-samples --threads 4 {input.vcf}
bcftools index -ft {output.vcf}
'''

    rule makeTGpops:
        input: tgped
        output:
            "reference/1kG_pops.txt",
            "reference/1kG_pops_unique.txt"
        cache: True
        resources:
            mem_mb = 10000,
            time_min = 30
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
        cache: True
        resources:
            mem_mb = 10000,
            time_min = 30
        shell:
            '''
cp {input} {output[0]}
awk 'NR > 1 {{print $3}}' {input} | sort | uniq > {output[1]}
'''


if REF == '1kG' or creftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = "reference/1000gFounders.{gbuild}.chr{chrom}.vcf.gz" if REF == '1kG' else config['custom_ref']['file'],
            tbi = "reference/1000gFounders.{gbuild}.chr{chrom}.vcf.gz.tbi" if REF == '1kG' else config['custom_ref']['file'] + '.tbi'
        output: temp("reference/{refname}.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        conda: "../envs/bcftools.yaml"
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
        cache: True
        threads: 2
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        conda: "../envs/bcftools.yaml"
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
        cache: True
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        conda: "../envs/bcftools.yaml"
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
        resources:
            mem_mb = 10000,
            time_min = 30
        container: 'docker://befh/flippyr:0.5.3'
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
        resources:
            mem_mb = 10000,
            time_min = 30
        container: 'docker://befh/r_env_gwasamplefilt:5'
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
        resources:
            mem_mb = 10000,
            time_min = 30
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

    # Index bcf
    rule Ref_IndexVcf:
        input: "reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: "reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        cache: True
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
        cache: True
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        conda: "../envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''

rule get_panelvars:
    input: "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz"
    output: 'reference/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps'
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/bcftools.yaml"
    shell: "bcftools query -f '%ID\n' {input} > {output}"
