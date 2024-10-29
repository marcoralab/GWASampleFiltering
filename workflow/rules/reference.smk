'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

# from scripts.parse_config_GWASampleFiltering import parser

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

from urllib.request import urlopen
from urllib.error import URLError

try:
    response = urlopen('https://www.google.com/', timeout=10)
    iconnect = True
except URLError as ex:
    iconnect = False


class dummyprovider:
    def remote(string_, allow_redirects="foo", immediate_close="bar"):
        return string_


if not 'config' in locals() or not config:
    configfile: "config/references.yaml"
    config['ref_only'] = True

if iconnect and not ('nointernet' in config and config['nointernet']):
    FTP = FTPRemoteProvider()
    HTTP = HTTPRemoteProvider()
else:
    FTP = dummyprovider
    HTTP = dummyprovider

BPLINK = ["bed", "bim", "fam"]

def map_genome_build(genome_build):
    if genome_build.lower() in ['hg19', 'hg37', 'grch37', 'b37']:
        return 'hg19'
    elif genome_build.lower() in ['hg38', 'grch38', 'b38']:
        return 'GRCh38'
    else:
        raise ValueError("Invalid genome build!")


if isinstance(config['genome_build'], str):
    BUILD = map_genome_build(config['genome_build'])
elif (isinstance(config['genome_build'], list)
      and isinstance(config['genome_build'][0], str)):
    BUILD = [map_genome_build(x) for x in config['genome_build']]
else:
    raise ValueError("Genome build must be a string or list of strings.")

localrules: download_tg_fa, download_tg_ped, download_tg_chrom, download_md5_b38, download_md5_hg19


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

tgbase = "ftp-trace.ncbi.nih.gov/1000genomes/ftp/"
tgbase_38 = "ftp.1000genomes.ebi.ac.uk/vol1/ftp/"

tgfasta = dict(
    hg19=(tgbase + "technical/reference/human_g1k_v37.fasta.gz"),
    GRCh38=(tgbase_38 + 'technical/reference/GRCh38_reference_genome/'
            + 'GRCh38_full_analysis_set_plus_decoy_hla.fa'))


if default_ref:
    REF = '1kG'
    reftype = 'vcfchr'
    ref_genotypes = "reference/1000gFounders.{gbuild}.chr{chrom}.vcf.gz"
else:
    REF = config['custom_ref']['name']
    reftype = detect_ref_type(config['custom_ref']['file'])
    ref_genotypes = config['custom_ref']['file']


if 'ref_only' in config and config['ref_only']:
    MISS = config['GenoMiss']
    rule all:
        input:
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD),
            panelvars = expand(
                'reference/panelvars_{refname}_{gbuild}'
                + '_allChr_maxmiss{miss}.snps',
                gbuild=BUILD, miss=MISS, refname=REF),
            vcf = expand(
                'reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz',
                gbuild=BUILD, miss=MISS, refname=REF),
            tgped = "reference/20130606_g1k.ped" if REF == '1kG' else [],
            pops = expand("reference/{refname}_pops.txt", refname=REF),
            pops_unique = expand("reference/{refname}_pops_unique.txt",
                                 refname=REF)


if default_ref:
    include: 'reference_1kG.smk'
else:
    rule make_custom_pops:
        input: config['custom_ref']['custom_pops']
        output:
            "reference/{refname}_pops.txt",
            "reference/{refname}_pops_unique.txt"
        #cache: True
        resources:
            mem_mb = 10000,
            time_min = 30
        shell:
            '''
cp {input} {output[0]}
awk 'NR > 1 {{print $3}}' {input} | sort | uniq > {output[1]}
'''


rule make_contig_convert_b38:
    output: 'reference/chr_name_conv.txt'
    resources:
        time_min = 4
    shell: '''
cat <(for i in {{1..22}} X Y; do echo chr$i $i; done) <(echo chrM MT) > {output}
'''


def fasta_remote(wc):
    gbuild = wc.gbuild
    if gbuild == 'GRCh38':
        return FTP.remote(tgfasta['GRCh38'], immediate_close=True)
    else:
        return HTTP.remote(tgfasta[gbuild])

def fasta_md5(wc):
    gbuild = wc.gbuild
    if gbuild == 'GRCh38':
        return '64b32de2fc934679c16e83a2bc072064'
    elif gbuild == 'hg19':
        return '45f81df94f0408d082363e34a081ed81'

rule download_tg_fa:
    input: fasta_remote
    output:
        "reference/human_g1k_{gbuild}.fasta",
        "reference/human_g1k_{gbuild}.fasta.fai"
    params:
        md5 = lambda wc: fasta_md5(wc)
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: '../envs/bcftools.yaml'
    #cache: True
    shell:
        '''
md5sum -c <(echo {params.md5} {input[0]})
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

if reftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = ref_genotypes,
            tbi = ref_genotypes + '.tbi',
            conv = 'reference/chr_name_conv.txt'
        output:
            temp('reference/{refname}.{gbuild}.chr{chrom}'
                 + '.maxmiss{miss}.vcf.gz')
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        container: "docker://befh/bcftools-htslib-samtools:1.15"
        shell:
            '''
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' \
  --threads 2 | \
bcftools annotate --rename-chrs {input.conv} --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
'''

    rule Reference_cat:
        input:
            vcfs = expand("reference/{{refname}}.{{gbuild}}"
                          + ".chr{chrom}.maxmiss{{miss}}.vcf.gz",
                          chrom=list(range(1, 23)))
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = ("reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz"
                   + ".tbi")
        #cache: True
        threads: 2
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        container: "docker://befh/bcftools-htslib-samtools:1.15"
        shell:
            '''
bcftools concat {input.vcfs} -Oz -o {output.vcf} --threads 2
bcftools index -ft {output.vcf}
'''
elif creftype == 'vcf':
    rule Reference_prep:
        input:
            vcf = config['custom_ref']['file'],
            tbi = config['custom_ref']['file'] + '.tbi',
            conv = 'reference/chr_name_conv.txt'
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = ("reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz"
                   + ".tbi")
        #cache: True
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        container: "docker://befh/bcftools-htslib-samtools:1.15"
        shell:
            '''
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' \
  --threads 2 | \
bcftools annotate --rename-chrs {input.conv} --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 \
  -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''
else:  # PLINK fileset of all chromosomes
    # align custom ref to fasta refrence
    rule Ref_Flip:
        input:
            bim = config['custom_ref']['file'] + '.bim',
            bed = config['custom_ref']['file'] + '.bed',
            fam = config['custom_ref']['file'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand("reference/{{refname}}_{{gbuild}}_flipped.{ext}",
                        ext=BPLINK))
        resources:
            mem_mb = 10000,
            time_min = 30
        container: 'docker://befh/flippyr:0.6.1'
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
        shell: '''
Rscript scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}'''

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
            csi = rules.Ref_IndexVcf.output,
            conv = 'reference/chr_name_conv.txt'
        output:
            vcf = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        #cache: True
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        container: "docker://befh/bcftools-htslib-samtools:1.15"
        shell:
            '''
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' \
  --threads 2 | \
bcftools annotate --rename-chrs {input.conv} --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 \
  -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''

rule get_panelvars:
    input: "reference/{refname}_{gbuild}_allChr_maxmiss{miss}.vcf.gz"
    output: 'reference/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps'
    resources:
        mem_mb = 10000,
        time_min = 30
    container: "docker://befh/bcftools-htslib-samtools:1.15"
    shell: "bcftools query -f '%ID\n' {input} > {output}"
