'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

from scripts.parse_config import parser
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

FTP = FTPRemoteProvider() if iconnect else dummyprovider
HTTP = HTTPRemoteProvider() if iconnect else dummyprovider

isMinerva = "hpc.mssm.edu" in socket.getfqdn()

configfile: "config/config.yaml"

# shell.executable("/bin/bash")
#
# if isMinerva:
#     anacondapath = sys.exec_prefix + "/bin"
#     shell.prefix(". ~/.bashrc; PATH={}:$PATH; ".format(anacondapath))
#

BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
start, FAMILY, SAMPLE, DATAOUT = parser(config)

istg = config['isTG'] if 'isTG' in config else False

# QC Steps:
QC_snp = True
QC_callRate = True
union_panel_TF = True

if isMinerva:
    com = {'flippyr': 'flippyr', 'plink': 'plink --keep-allele-order',
           'plink2': 'plink', 'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R',
           'king': 'king', 'faidx': 'samtools faidx'}
    loads = {'flippyr': 'module load plink/1.90b6.10', 'plink': 'module load plink/1.90b6.10',
             'bcftools': 'module load bcftools/1.9', 'faidx': 'module load samtools',
             'king': 'module unload gcc; module load king/2.1.6',
             'R': ('module load R/3.6.3 pandoc/2.6 udunits/2.2.26; ',
                   'RSTUDIO_PANDOC=$(which pandoc)')}
else:
    com = {'flippyr': 'flippyr',
           'plink': 'plink --keep-allele-order', 'plink2': 'plink', 'faidx': 'samtools faidx',
           'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R', 'king': 'king'}
    loads = {'flippyr': 'echo running flippyr', 'plink': 'echo running plink',
             'bcftools': 'echo running bcftools',  'R': 'echo running R',
             'king': 'echo running KING', 'faidx': 'echo running samtools faidx'}

def decorate(text):
    return expand(DATAOUT + "/{sample}_" + text,
                  sample=SAMPLE)

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

if (not "custom_ref" in config) or (config['custom_ref'] == False):
    REF = '1kG'
else:
    REF = config['custom_ref_name']
    creftype = detect_ref_type(config['custom_ref'])

extraref = False
if ("extra_ref" in config) and (config['extra_ref'] != False):
    extraref = True
    ereftype = detect_ref_type(config['extra_ref'])
else:
    ereftype = 'none'

outs = {
    "report": expand(DATAOUT + "/stats/{sample}_GWAS_QC.html", sample=SAMPLE),
    "exclusions": expand(DATAOUT + "/{sample}_exclude.samples", sample=SAMPLE),
    "filtered": expand(DATAOUT + "/{sample}_Excluded.{ext}",
        sample=SAMPLE, ext=BPLINK)}

outputs = [outs[x] for x in config["outputs"]]
outputs = flatten(outputs)

rule all:
    input:
        expand(DATAOUT + "/{sample}_pruned.snplist", sample=SAMPLE)

    #expand(DATAOUT + "/{sample}_{refname}_merged.vcf", sample=SAMPLE, refname=REF) # outputs


if config['mirror'] == 'ncbi':
    tgbase = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/"
elif config['mirror'] == 'ebi':
    tgbase = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"

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


rule download_tg_chrom:
   input:
       FTP.remote(tgurl),
       FTP.remote(tgurl + ".tbi"),
   output:
       temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz"),
       temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi")
   shell: "cp {input[0]} {output[0]}; cp {input[1]} {output[1]}"


rule download_tg_fa:
   input:
       FTP.remote(tgfa + ".gz") if BUILD == 'hg19' else FTP.remote(tgfa)
   output:
       "reference/human_g1k_{gbuild}.fasta",
       "reference/human_g1k_{gbuild}.fasta.fai"
   shell:
           """
{loads[faidx]}
if [[ "{input[0]}" == *.gz ]]; then
  zcat {input[0]} > {output[0]}
else
  cp {input[0]} {output[0]}
fi
{com[faidx]} {output[0]}"""

rule download_tg_ped:
    input:
        FTP.remote(tgped),
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
            r"""
            awk -F "\t" '!($12 != 0 || $10 != 0 || $9 != 0 || $3 != 0 || $4 != 0) {{print $2}}' \
            {input} > {output}
            """

    rule makeTGpops:
        input: tgped
        output:
            "reference/1kG_pops.txt",
            "reference/1kG_pops_unique.txt"
        shell:
            """
            awk 'BEGIN {{print "FID","IID","Population"}} NR>1 {{print $1,$2,$7}}' \
            {input} > {output[0]}
            cut -f7 {input} | sed 1d | sort | uniq > {output[1]}
            """
else:
    rule make_custom_pops:
        input: config['custom_pops']
        output:
            "reference/{refname}_pops.txt",
            "reference/{refname}_pops_unique.txt"
        shell:
            """
            cp {input} {output[0]}
            awk 'NR > 1 {{print $3}}' {input} | sort | uniq > {output[1]}
            """


if REF == '1kG' or creftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz" if REF == '1kG' else config['custom_ref'],
            tbi = "reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi" if REF == '1kG' else config['custom_ref'] + '.tbi'
        output: temp("reference/{refname}.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        shell:
            """
    {loads[bcftools]}
    {com[bcftools]} norm -m- {input.vcf} --threads 2 | \
    {com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
    {com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
    """

    rule Reference_cat:
        input:
            vcfs = expand("reference/{{refname}}.{{gbuild}}.chr{chrom}.maxmiss{{miss}}.vcf.gz",
                          chrom = list(range(1, 23)))
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        shell:
            """
    {loads[bcftools]}
    {com[bcftools]} concat {input.vcfs} -Oz -o {output.vcf} --threads 2
    {com[bcftools]} index -ft {output.vcf}
    """
elif creftype == 'vcf':
    rule Reference_prep:
        input:
            vcf = config['custom_ref'],
            tbi = config['custom_ref'] + '.tbi'
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        shell:
            """
    {loads[bcftools]}
    {com[bcftools]} norm -m- {input.vcf} --threads 2 | \
    {com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
    {com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
    {com[bcftools]} index -ft {output.vcf}
    """
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
        shell:
            """
    {loads[flippyr]}
    {com[flippyr]} -p {input.fasta} -o reference/{wildcards.refname}_{wildcards.gbuild}_flipped {input.bim}"""

    rule Ref_ChromPosRefAlt:
        input:
            flipped = "reference/{refname}_{gbuild}_flipped.bim"
        output:
            bim = temp("reference/{refname}_{gbuild}_flipped_ChromPos.bim"),
            snplist = temp("reference/{refname}_{gbuild}_flipped_snplist")
        shell:
            """
    {loads[R]}
    {com[R]} scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"""

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
        shell:
            """
    {loads[plink]}
    {com[plink2]} --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
    --real-ref-alleles --out {params.out}"""

    # Index bcf
    rule Ref_IndexVcf:
        input: "reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: "reference/{refname}_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        shell: '{loads[bcftools]}; {com[bcftools]} index -f {input}'

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
        shell:
            """
{loads[bcftools]}
{com[bcftools]} norm -m- {input.vcf} --threads 2 | \
{com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
{com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
{com[bcftools]} index -ft {output.vcf}
"""

if ereftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = config['extra_ref'],
            tbi = config['extra_ref'] + '.tbi'
        output: temp(DATAOUT + "/extraref.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        shell:
            """
    {loads[bcftools]}
    {com[bcftools]} norm -m- {input.vcf} --threads 2 | \
    {com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
    {com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
    """

    rule Reference_cat_extra:
        input:
            vcfs = expand(DATAOUT + "/extraref.{{gbuild}}.chr{chrom}.maxmiss{{miss}}.vcf.gz",
                          chrom = list(range(1, 23)))
        output:
            vcf = DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        shell:
            """
    {loads[bcftools]}
    {com[bcftools]} concat {input.vcfs} -Oz -o {output.vcf} --threads 2
    {com[bcftools]} index -ft {output.vcf}
    """
elif ereftype == 'vcf':
    rule Reference_prep_extra:
        input:
            vcf = config['extra_ref'],
            tbi = config['extra_ref'] + '.tbi'
        output:
            vcf = DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        shell:
            """
    {loads[bcftools]}
    {com[bcftools]} norm -m- {input.vcf} --threads 2 | \
    {com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
    {com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
    {com[bcftools]} index -ft {output.vcf}
    """
elif ereftype != 'none': #PLINK fileset of all chromosomes
    # align custom ref to fasta refrence
    rule Ref_Flip_extra:
        input:
            bim = config['extra_ref'] + '.bim',
            bed = config['extra_ref'] + '.bed',
            fam = config['extra_ref'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand(DATAOUT + "/extraref_{{gbuild}}_flipped.{ext}", ext=BPLINK))
        shell:
            """
    {loads[flippyr]}
    {com[flippyr]} -p {input.fasta} -o {DATAOUT}/extraref_{wildcards.gbuild}_flipped {input.bim}"""

    rule Ref_ChromPosRefAlt_extra:
        input:
            flipped = DATAOUT + "/extraref_{gbuild}_flipped.bim"
        output:
            bim = temp(DATAOUT + "/extraref_{gbuild}_flipped_ChromPos.bim"),
            snplist = temp(DATAOUT + "/extraref_{gbuild}_flipped_snplist")
        shell:
            """
    {loads[R]}
    {com[R]} scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"""

    # Recode sample plink file to vcf
    rule Ref_Plink2Vcf_extra:
        input:
            bim = rules.Ref_ChromPosRefAlt_extra.output.bim,
            flipped = rules.Ref_Flip_extra.output
        output:
            temp(DATAOUT + "/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz")
        params:
            out = DATAOUT + "/extraref_{gbuild}_unQC_maxmissUnfilt",
            inp = DATAOUT + "/extraref_{gbuild}_flipped"
        shell:
            """
    {loads[plink]}
    {com[plink2]} --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
    --real-ref-alleles --out {params.out}"""

    # Index bcf
    rule Ref_IndexVcf_extra:
        input: DATAOUT + "/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: DATAOUT + "/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        shell: '{loads[bcftools]}; {com[bcftools]} index -f {input}'

    rule Reference_prep_extra:
        input:
            vcf = rules.Ref_Plink2Vcf_extra.output,
            csi = rules.Ref_IndexVcf_extra.output
        output:
            vcf = DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        shell:
            """
{loads[bcftools]}
{com[bcftools]} norm -m- {input.vcf} --threads 2 | \
{com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
{com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
{com[bcftools]} index -ft {output.vcf}
"""

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
    shell:
        """
{loads[bcftools]}
{com[bcftools]} query -f '%ID\n' {input} > {output}
"""

rule merge_panelvars:
    input:
        ref = expand('reference/panelvars_{refname}.snps', refname=REF),
        eref = expand(DATAOUT + "/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
                      gbuild=BUILD, miss=config['QC']['GenoMiss'])
    output: PVARS if erefc else "do.not"
    shell: 'cat {input} | sort | uniq > {output}'
