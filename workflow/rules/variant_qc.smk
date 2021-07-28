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
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.stem} --freq --out {params.out}
plink --keep-allele-order --bfile {params.stem} --freqx --out {params.out}
plink --keep-allele-order --bfile {params.stem} --geno {params.miss} \
--maf {params.MAF} --hardy --hwe {params.HWE} --make-bed --out {params.out}
'''

# ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input: expand("{{dataout}}/{{sample}}_SnpQc.{ext}", ext=BPLINK) if qc_type['variant'] else start['files']
    output:
        expand("{{dataout}}/{{sample}}_callRate.{ext}", ext=BPLINK, dataout = DATAOUT),
        "{dataout}/{sample}_callRate.imiss",
        touch("{dataout}/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out if qc_type['variant'] else start['stem'],
        miss = config['QC']['SampMiss'],
        out = "{dataout}/{sample}_callRate"
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} --mind {params.miss} \
  --missing --make-bed --out {params.out}
'''
