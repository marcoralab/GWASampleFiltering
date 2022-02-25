tgurls = {'hg19': {'vcf': (tgbase + "release/20130502/ALL.chr{chrom}."
                           + "phase3_shapeit2_mvncall_integrated_v5a."
                           + "20130502.genotypes.vcf.gz"),
                   'md5': tgbase + 'current.tree'},
          'GRCh38': {'vcf': (tgbase_38
                             + 'data_collections/1000_genomes_project/'
                             + 'release/20181203_biallelic_SNV/'
                             + 'ALL.chr{chrom}.shapeit2_integrated_v1a.'
                             + 'GRCh38.20181129.phased.vcf.gz'),
                     'md5': (tgbase_38
                             + 'data_collections/1000_genomes_project/'
                             + 'release/20181203_biallelic_SNV/'
                             + '20181203_biallelic_SNV_manifest.txt')}}

tgurls['ped'] = (tgbase
                 + "technical/working/20130606_sample_info/20130606_g1k.ped")

def md5_remote(wc):
    gbuild = wc.gbuild
    if gbuild == 'GRCh38':
        return
    return HTTP.remote(tgurls[wc.gbuild]['md5'])


rule download_md5_b38:
    input: FTP.remote(tgurls['GRCh38']['md5'], immediate_close=True)
    output: 'reference/md5.GRCh38.vcfs.txt'
    shell:
        r'''
awk 'match($1, "chr[0-9]+", chrom) && match($1, "[.]vcf.+", ext) \
     {{a="reference/1000gRaw.GRCh38."chrom[0]ext[0]; print $3, a}}' \
     {input} > {output}
'''

rule download_md5_hg19:
    input: HTTP.remote(tgurls['hg19']['md5'])
    output: 'reference/md5.hg19.vcfs.txt'
    shell:
        r'''
awk '{{FS="\t"}} match($1, "chr[0-9]+", chrom) && match($1, "[.]vcf.+", ext) \
     && $1 ~ "v5a[.]20130502[.]genotypes" \
     {{a="reference/1000gRaw.hg19."chrom[0]ext[0]; print $5, a}} \
     $1 ~ "20130606_g1k.ped" {{print $5, "reference/20130606_g1k.ped"}}' \
     {input} > {output}
'''


def vcf_remote(wc):
    gbuild = wc.gbuild
    vcf = tgurls[gbuild]['vcf']
    if gbuild == 'GRCh38':
        outs = {'vcf': FTP.remote(vcf, immediate_close=True),
                'tbi': FTP.remote(vcf + '.tbi', immediate_close=True)}
    else:
        outs = {'vcf': HTTP.remote(vcf), 'tbi': HTTP.remote(vcf + '.tbi')}
    return {'md5': f'reference/md5.{gbuild}.vcfs.txt', **outs}


rule download_tg_chrom:
    input: unpack(vcf_remote)
    output:
        vcf = temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz"),
        tbi = temp("reference/1000gRaw.{gbuild}.chr{chrom}.vcf.gz.tbi")
    resources:
        mem_mb = 10000,
        time_min = 30
    shell:
        '''
cp {input.vcf} {output.vcf}
cp {input.tbi} {output.tbi}
sleep 10
md5sum -c <(awk '$2 ~ "chr{wildcards.chrom}[.]"' {input.md5})
'''

tgped = "reference/20130606_g1k.ped"

rule download_tg_ped:
    input:
        ped = FTP.remote(tgurls['ped'], immediate_close=True),
        md5 = 'reference/md5.hg19.vcfs.txt'
    output: tgped
    cache: True
    resources:
        mem_mb = 10000,
        time_min = 30
    shell:
        '''
cp {input.ped} {output}
md5sum -c <(awk '$2 ~ "20130606_g1k[.]ped"' {input.md5})
'''

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
bcftools view -S {input.founders} --force-samples --threads 2 {input.vcf} | \
  bcftools view -s '^NA20299,NA20314,NA20274,HG01880' \
  -Oz -o {output.vcf} --force-samples --threads 2
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
