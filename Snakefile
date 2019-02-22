'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

from scripts.parse_config import parser
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
import socket
import getpass

FTP = FTPRemoteProvider()

isMinerva = "hpc.mssm.edu" in socket.getfqdn()

configfile: "./config.yaml"

shell.executable("/bin/bash")

if isMinerva:
    shell.prefix("PATH=" + config["anaconda"] + ":$PATH; ")

BPLINK = ["bed", "bim", "fam"]
RWD = os.getcwd()
start, FAMILY, SAMPLE, DATAOUT = parser(config)

# QC Steps:
QC_snp = True
QC_callRate = True

if isMinerva:
    com = {'flippyr': 'flippyr', 'plink': 'plink --keep-allele-order',
           'plink2': 'plink', 'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R',
           'king': 'king'}
    loads = {'flippyr': 'module load plink/1.90', 'plink': 'module load plink/1.90',
             'bcftools': 'module load bcftools/1.9',
             'king': 'module unload gcc; module load king/2.1.6',
             'R': ('module load R/3.4.3 pandoc/2.1.3 udunits/2.2.26; ',
                   'RSTUDIO_PANDOC=$(which pandoc)')}
else:
    com = {'flippyr': 'flippyr',
           'plink': 'plink --keep-allele-order', 'plink2': 'plink',
           'bcftools': 'bcftools', 'R': 'Rscript', 'R2': 'R', 'king': 'king'}
    loads = {'flippyr': 'echo running flippyr', 'plink': 'echo running plink',
             'bcftools': 'echo running bcftools',  'R': 'echo running R',
             'king': 'echo running KING'}

if getpass.getuser() == "sheaandrews":
    com["flippyr"] = '/Users/sheaandrews/Programs/flippyr/flippyr.py'


def decorate(text):
    return expand(DATAOUT + "/{sample}_" + text,
                  sample=SAMPLE)

localrules: all, download_tg, download_tg_chrom

def flatten(nested):
    flat = []
    for el in nested:
        if not isinstance(el, list):
            flat.append(el)
        else:
            flat += flatten(el)
return flat

outs = dict(
    report=expand(DATAOUT + "/stats/{sample}_GWAS_QC.html", sample=SAMPLE),
    exclusions=expand(DATAOUT + "/{sample}_exclude.samples", sample=SAMPLE),
    filtered=expand(DATAOUT + "/{sample}_Excluded.{ext}",
        sample=SAMPLE, ext=BPLINK)

outputs = [outs[x] for x in config["outputs"]]
outputs = flatten(outputs)

rule all:
    input: outputs


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
        MAF = config['QC']['MAF'],
        HWE = config['QC']['HWE']
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.stem} --freq --out {params.out}
{com[plink]} --bfile {params.stem} --freqx --out {params.out}
{com[plink]} --bfile {params.stem} --geno {params.miss} \
--maf {params.MAF} --hardy --hwe {params.HWE} --make-bed --out {params.out}"""

# ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input: rules.snp_qc.output if QC_snp else start['files']
    output:
        expand(DATAOUT + "/{{sample}}_callRate.{ext}", ext=BPLINK),
        DATAOUT + "/{sample}_callRate.imiss",
        touch(DATAOUT + "/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out if QC_snp else start['stem'],
        miss = config['QC']['SampMiss'],
        out = DATAOUT + "/{sample}_callRate"
    shell:
        """
{loads[plink]}
{com[plink]} --bfile {params.indat} --mind {params.miss} \
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
rule relatedness_sample_prep:
    input: sexcheck_in_plink
    output:
        bed = temp(DATAOUT + "/{sample}_IBDQCfilt.bed"),
        bim = temp(DATAOUT + "/{sample}_IBDQCfilt.bim"),
        fam = temp(DATAOUT + "/{sample}_IBDQCfilt.fam")
    params:
        indat_plink = sexcheck_in_plink_stem,
        out = DATAOUT + "/{sample}_IBDQCfilt"
    shell:
        """
{loads[plink]}
{com[plink2]} --bfile {params.indat_plink} \
  --geno 0.02 \
  --maf 0.02 \
  --memory 6000 \
  --make-bed --out {params.out}"""

if config['king']:
    rule relatedness_QC:
        input:
            bed = rules.relatedness_sample_prep.output.bed,
            bim = rules.relatedness_sample_prep.output.bim,
            fam = rules.relatedness_sample_prep.output.fam
        output:
            genome = DATAOUT + "/{sample}_IBDQC.kin0",
            kin = DATAOUT + "/{sample}_IBDQC.kin"
        params:
            out = DATAOUT + "/{sample}_IBDQC"
        shell:
            """
{loads[king]}
{com[king]} -b {input.bed} --related --degree 3 --prefix {params.out}
if [[ -f {params.out}.kin && !( -f {params.out}.kin0 ) ]]; then
  echo .kin file exists, but .kin0 does not. This is likely because there is
  echo only one FID. converting .kin to .kin0
  {loads[R]}
  {com[R]} scripts/kin2kin0.R {params.out}.kin
fi"""
else:
    rule relatedness_QC:
        input: rules.sample_prune.output
        output:
            genome = DATAOUT + "/{sample}_IBDQC.genome"
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
        genome = rules.relatedness_QC.output.genome,
        fam = sexcheck_in_plink_stem + ".fam"
    params:
        Family = FAMILY,
        king = config['king'],
        threshold = 0.1875
    output:
        out = DATAOUT + "/{sample}_exclude.relatedness",
        rdat = DATAOUT + "/{sample}_IBDQC.Rdata"
    shell:
        """
{loads[R]}; {com[R]}  scripts/relatedness_QC.R {input.genome} {params.threshold} \
{params.Family} {params.king} {output.out} {output.rdat}"""

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
        fasta = "data/human_g1k_v37.fasta"
    output:
        temp(expand(DATAOUT + "/{{sample}}_pruned_flipped.{ext}",
                    ext=BPLINK))
    shell:
        """
{loads[flippyr]}
{com[flippyr]} -p {input.fasta} {input.bim}"""

rule Sample_ChromPosRefAlt:
    input:
        flipped = DATAOUT + "/{sample}_pruned_flipped.bim"
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

tgbase = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"
tgurl = (tgbase + "release/20130502/ALL.chr{chrom}." +
    "phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
tgped = tgbase + "technical/working/20130606_sample_info/20130606_g1k.ped"
tgfa = tgbase + "technical/reference/human_g1k_v37.fasta"

predownload = True

rule download_tg_chrom:
   input:
       FTP.remote(tgurl, keep_local=True),
       FTP.remote(tgurl + ".tbi", keep_local=True),
   output:
       temp("data/1000gRaw.chr{chrom}.vcf.gz"),
       temp("data/1000gRaw.chr{chrom}.vcf.gz.tbi")
   shell: "cp {input[0]} {output[0]}; cp {input[1]} {output[1]}"

rule download_tg:
   input:
       FTP.remote(tgped, keep_local=True),
       FTP.remote(tgfa + ".gz", keep_local=True),
       FTP.remote(tgfa + ".fai", keep_local=True)
   output:
       "data/20130606_g1k.ped",
       "data/human_g1k_v37.fasta",
       "data/human_g1k_v37.fasta.fai"
   shell: "cp {input[0]} {output[0]}; zcat {input[1]} > {output[1]}; cp {input[2]} {output[2]}"

tgped = "data/20130606_g1k.ped"
if predownload:
    refraw = ["data/1000gRaw.chr{chrom}.vcf.gz",
              "data/1000gRaw.chr{chrom}.vcf.gz.tbi"]
else:
#    tgped = FTP.remote(tgped, keep_local=True)
    refraw = "/dev/urandom"

rule makeTGpops:
    input: tgped
    output:
        "data/1000genomes_pops.txt",
        "data/pops.txt"
    shell:
        """
awk 'BEGIN {{print "FID","IID","Population"}} NR>1 {{print $1,$2,$7}}' \
{input} > {output[0]}
cut -f7 {input} | sed 1d | sort | uniq > {output[1]}
"""

rule Reference_prep:
    input: refraw
    output: temp("data/1000Gchr{chrom}.maxmiss{miss}.vcf.gz")
    params:
        url = rules.download_tg_chrom.output[0] if predownload else tgurl
    shell:
        """
{loads[bcftools]}
{com[bcftools]} norm -m- {params.url} --threads 2 | \
{com[bcftools]} view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
{com[bcftools]} annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
"""

"""
Selects everyone who is unrelated or only has third degree relatives in
thousand genomes.
"""
rule Reference_foundersonly:
    input: tgped
    output:
        "data/20130606_g1k.founders"
    shell: r"""
awk -F "\t" '!($12 != 0 || $10 != 0 || $9 != 0 || $3 != 0 || $4 != 0) {{print $2}}' \
{input} > {output}
"""

rule Reference_cat:
    input:
        vcfs = expand("data/1000Gchr{chrom}.maxmiss{{miss}}.vcf.gz",
                      chrom = list(range(1, 23)))
    output:
        vcf = "data/1000genomes_allChr_maxmiss{miss}.vcf.gz",
        tbi = "data/1000genomes_allChr_maxmiss{miss}.vcf.gz.tbi"
    shell:
        """
{loads[bcftools]}
{com[bcftools]} concat {input.vcfs} -Oz -o {output.vcf} --threads 2
{com[bcftools]} index -ft {output.vcf}
"""

rule Reference_prune:
    input:
        vcf = expand("data/1000genomes_allChr_maxmiss{miss}.vcf.gz",
                     miss = config['QC']['GenoMiss']),
        prune = DATAOUT + "/{sample}_pruned_snplist",
        founders = "data/20130606_g1k.founders"
    output:
        vcf = temp(DATAOUT + "/{sample}_1kgpruned.vcf.gz"),
        tbi = temp(DATAOUT + "/{sample}_1kgpruned.vcf.gz.tbi")
    shell:
        """
{loads[bcftools]}
{com[bcftools]} view -i 'ID=@{input.prune}' -S {input.founders} \
-Oz -o {output.vcf} --force-samples {input.vcf} --threads 4
{com[bcftools]} index -ft {output.vcf}
"""

# Merge ref and sample
rule Merge_RefenceSample:
    input:
        bcf_1kg = DATAOUT + "/{sample}_1kgpruned.vcf.gz",
        tbi_1kg = DATAOUT + "/{sample}_1kgpruned.vcf.gz.tbi",
        bcf_samp = DATAOUT + "/{sample}_pruned_flipped.vcf.gz",
        csi_samp = DATAOUT + "/{sample}_pruned_flipped.vcf.gz.csi"
    params:
        miss = config['QC']['GenoMiss']
    output:
        out = DATAOUT + "/{sample}_1kg_merged.vcf"
    shell:
        r"""
{loads[bcftools]}
{com[bcftools]} merge -m none --threads 2 {input.bcf_1kg} {input.bcf_samp} | \
{com[bcftools]} view  -i 'F_MISSING <= {params.miss}' -Ov -o {output.out} --threads 2"""

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
        oldfam = DATAOUT + "/{sample}_pruned.fam",
        newfam = DATAOUT + "/{sample}_1kg_merged.fam",
        tgped = tgped
    output: DATAOUT + "/{sample}_1kg_merged_fixed.fam"
    shell:
        """
{loads[R]}
{com[R]} scripts/fix_fam.R {input.oldfam} {input.newfam} {output} {input.tgped}"""

# PCA analysis to identify population outliers
rule PcaPopulationOutliers:
    input:
        plink = expand(DATAOUT + "/{{sample}}_1kg_merged.{ext}", ext=BPLINK),
        fam = rules.fix_fam.output,
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

# Rscript to identify population outliers
rule ExcludePopulationOutliers:
    input:
        eigenval = DATAOUT + "/{sample}_1kg_merged.eigenval",
        eigenvec = DATAOUT + "/{sample}_1kg_merged.eigenvec",
        fam = DATAOUT + "/{sample}_pruned.fam",
        tgped = tgped
    output:
        excl = DATAOUT + "/{sample}_exclude.pca",
        rmd = temp(DATAOUT + "/{sample}_pca.Rdata")
    params:
        samp = "{sample}",
        superpop = config['superpop']
    shell:
        """
{loads[R]}
scripts/PCA_QC.R -s {params.samp} -p {params.superpop} \
--vec {input.eigenvec} --val {input.eigenval} \
-b {input.tgped} -t {input.fam} -o {output.excl} -R {output.rmd}
"""

# Run PCA to control for population stratification
if config['pcair']:
    assert config['king'], "You must use KING for relatedness if using PCAiR!"
    rule ancestryFilt:
        input:
            plink = expand(DATAOUT + "/{{sample}}_pruned.{ext}",
                           ext=BPLINK),
            exclude = DATAOUT + "/{sample}_exclude.pca"
        output:
            temp(expand(DATAOUT + "/{{sample}}_filtered_PCApre.{ext}",
                 ext=BPLINK)),
        params:
            indat = DATAOUT + "/{sample}_pruned",
            plinkout = DATAOUT + "/{sample}_filtered_PCApre"
        shell:
            r"""
{loads[plink]}
{com[plink]} --bfile {params.indat} \
  --remove {input.exclude} \
  --make-bed --out {params.plinkout}
"""

    rule pcair_king:
        input:
            bed = rules.relatedness_sample_prep.output.bed,
            bim = rules.relatedness_sample_prep.output.bim,
            fam = rules.relatedness_sample_prep.output.fam
        output:
            genome = DATAOUT + "/{sample}_pcairKING.kin0",
            kin = DATAOUT + "/{sample}_pcairKING.kin"
        params:
            out = DATAOUT + "/{sample}_pcairKING"
        shell:
            """
{loads[king]}
{com[king]} -b {input.bed} --kinship --prefix {params.out}
if [[ -f {params.out}.kin && !( -f {params.out}.kin0 ) ]]; then
  echo .kin file exists, but .kin0 does not. This is likely because there is
  echo only one FID. converting .kin to .kin0
  {loads[R]}
  {com[R]} scripts/kin2kin0.R {params.out}.kin
fi"""

    rule filterKING:
        input:
            king = expand(DATAOUT + "/{{sample}}_pcairKING.{ext}",
                          ext=["kin", "kin0"]),
            exclude = DATAOUT + "/{sample}_exclude.pca"
        output:
            temp(expand(DATAOUT + "/{{sample}}_pcairKING.popfilt.{ext}",
                        ext=["kin", "kin0"]))
        params:
            indat = DATAOUT + "/{sample}_pcairKING",
        shell:
            r"""
{loads[R]}
{com[R]} scripts/filterKing.R {params.indat} {input.exclude}
"""

    rule PCAPartitioning:
        input:
            plink = rules.ancestryFilt.output,
            king = rules.filterKING.output,
            iterative = rules.relatedness_sample_fail.output.out
        output:
            expand(DATAOUT + "/{{sample}}_filtered_PCApre.{ext}",
                   ext=['unrel', 'partition.log'])
        params:
            stem = rules.ancestryFilt.params.plinkout,
            king = rules.filterKING.params.indat
        shell:
            """
{loads[R]}
{com[R]} scripts/PartitionPCAiR.R {params.stem} {params.king} {input.iterative}
"""

    rule stratFrq:
        input:
            plink = rules.ancestryFilt.output,
            unrel = rules.PCAPartitioning.output[0],
        output: DATAOUT + "/{sample}_filtered_PCAfreq.frqx"
        params:
            indat = rules.ancestryFilt.params.plinkout,
            out = DATAOUT + "/{sample}_filtered_PCAfreq"
        shell:
            """
{loads[plink]}
{com[plink]} --bfile {params.indat} --freqx \
  --within {input.unrel} --keep-cluster-names unrelated \
  --out {params.out}
"""

    rule PopulationStratification:
        input:
            plink = rules.ancestryFilt.output,
            unrel = rules.PCAPartitioning.output[0],
            frq = rules.stratFrq.output
        output:
            expand(DATAOUT + "/{{sample}}_filtered_PCA.{ext}",
                   ext=['eigenval', 'eigenvec'])
        params:
            indat = rules.ancestryFilt.params.plinkout,
            out = DATAOUT + "/{sample}_filtered_PCA"
        shell:
            """
{loads[plink]}
{com[plink]} --bfile {params.indat} --read-freq {input.frq} --pca 10 \
  --within {input.unrel} --pca-cluster-names unrelated \
  --out {params.out}
"""
else:
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
        IBD_stats = decorate2("IBDQC.Rdata"),
        PCA_rdat = decorate2("pca.Rdata"),
        PopStrat_eigenval = decorate2("filtered_PCA.eigenval"),
        PopStrat_eigenvec = decorate2("filtered_PCA.eigenvec"),
        partmethod = rules.PCAPartitioning.output[1] if config["pcair"] else "/dev/null"
    output:
        DATAOUT + "/stats/{sample}_GWAS_QC.html"
    params:
        rwd = RWD,
        Family = FAMILY,
        pi_threshold = 0.1875,
        output_dir = DATAOUT + "/stats",
        idir = DATAOUT + "/stats/md/{sample}",
        geno_miss = config['QC']['GenoMiss'],
        samp_miss = config['QC']['SampMiss'],
        MAF = config['QC']['MAF'],
        HWE = config['QC']['HWE'],
        superpop = config['superpop'],
        partmethod = rules.PCAPartitioning.output[1] if config["pcair"] else "none"
    shell:
        """
{loads[R]}
{com[R2]} -e 'nm <- sample(c("Shea J. Andrews", "Brian Fulton-Howard"), \
replace=F); nm <- paste(nm[1], "and", nm[2]); \
rmarkdown::render("{input.script}", output_dir = "{params.output_dir}", \
output_file = "{output}", intermediates_dir = "{params.idir}", \
params = list(rwd = "{params.rwd}", Sample = "{wildcards.sample}", \
auth = nm, Path_SexFile = "{input.SexFile}", Path_hwe = "{input.hwe}", \
Path_frq = "{input.frq}", Path_frqx = "{input.frqx}", \
Path_imiss = "{input.imiss}", Path_HetFile = "{input.HetFile}", \
pi_threshold = {params.pi_threshold}, Family = {params.Family}, \
Path_IBD_stats = "{input.IBD_stats}", Path_PCA_rdat = "{input.PCA_rdat}", \
Path_PopStrat_eigenval = "{input.PopStrat_eigenval}", \
Path_PopStrat_eigenvec = "{input.PopStrat_eigenvec}", maf = {params.MAF}, \
hwe = {params.HWE}, missing_geno = {params.geno_miss}, \
partmethod = "{params.partmethod}", \
missing_sample = {params.samp_miss}, superpop = "{params.superpop}"))' --slave
"""
