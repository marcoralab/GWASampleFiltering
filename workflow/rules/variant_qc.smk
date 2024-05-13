# ---- Exlude SNPs with a high missing rate and low MAF----
rule snp_qc:
    input: start['files']
    output:
        temp(multiext("{dataout}/{sample}_SnpQc", '.bed', '.bim', '.fam')),
        "{dataout}/{sample}_SnpQc.hwe",
        "{dataout}/{sample}_SnpQc.frq",
        "{dataout}/{sample}_SnpQc.frqx",
    params:
        stem = apply_prefix(start['stem']),
        out = apply_prefix("{dataout}/{sample}_SnpQc"),
        miss = config['QC']['GenoMiss'],
        MAF = config['QC']['MAF'],
        HWE = config['QC']['HWE']
    resources:
        mem_mb = 30000,
        time_min = 30
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --allow-no-sex --bfile {params.stem} --freq --out {params.out}
plink --keep-allele-order --allow-no-sex --bfile {params.stem} --freqx --out {params.out}
plink --keep-allele-order --allow-no-sex --bfile {params.stem} --geno {params.miss} \
--maf {params.MAF} --hardy --hwe {params.HWE} --make-bed --out {params.out}
'''

# ---- Exclude Samples with high missing rate ----
rule sample_callRate:
    input: multiext("{dataout}/{sample}_SnpQc", '.bed', '.bim', '.fam') if qc_type['variant'] else start['files']
    output:
        multiext("{dataout}/{sample}_callRate", '.bed', '.bim', '.fam'),
        touch("{dataout}/{sample}_callRate.irem")
    params:
        indat = rules.snp_qc.params.out if qc_type['variant'] else apply_prefix(start['stem']),
        miss = config['QC']['SampMiss'],
        out = apply_prefix("{dataout}/{sample}_callRate")
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --allow-no-sex --bfile {params.indat} \
  --mind {params.miss} --make-bed --out {params.out}
'''

rule sample_callRate_imiss:
    input: rules.sample_callRate.input
    output:
        "{dataout}/{sample}_imiss_callRate.imiss"
    params:
        indat = rules.sample_callRate.params.indat,
        out = apply_prefix("{dataout}/{sample}_imiss_callRate")
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --allow-no-sex --bfile {params.indat} \
  --missing --out {params.out}
'''
