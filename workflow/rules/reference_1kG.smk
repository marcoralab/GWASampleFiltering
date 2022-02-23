tgurls = dict(
    hg19=dict(
        vcf=(tgbase + "release/20130502/ALL.chr{chrom}."
             + "phase3_shapeit2_mvncall_integrated_v5a."
             + "20130502.genotypes.vcf.gz")),
    GRCh38=dict(
        vcf=(tgbase_38
             + 'data_collections/1000_genomes_project/'
             + 'release/20181203_biallelic_SNV/'
             + 'ALL.chr{chrom}.shapeit2_integrated_v1a.'
             + 'GRCh38.20181129.phased.vcf.gz'))
    )
tgurls = {k: {**v, 'tbi': v['vcf'] + '.tbi'} for k, v in tgurls.items()}
tgurls['ped'] = (tgbase
                 + "technical/working/20130606_sample_info/20130606_g1k.ped")


rule download_tg_chrom:
    input:
        vcf = lambda wildcards: HTTP.remote(tgurls[wildcards.gbuild]['vcf']),
        tbi = lambda wildcards: HTTP.remote(tgurls[wildcards.gbuild]['tbi'])
    output:
        vcf = temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz"),
        tbi = temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi")
    resources:
        mem_mb = 10000,
        time_min = 30
    shell: "cp {input.vcf} {output.vcf}; cp {input.tbi} {output.tbi}"

tgped = "reference/20130606_g1k.ped"

rule download_tg_ped:
    input:
        HTTP.remote(tgurls['ped']),
    output: tgped
    cache: True
    resources:
        mem_mb = 10000,
        time_min = 30
    shell: "cp {input} {output}"

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
    output:
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
