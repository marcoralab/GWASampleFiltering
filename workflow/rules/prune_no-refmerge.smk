if not config['full_pipeline']:
    include: 'variant_qc.smk'

if "sampleqc_in_plink" not in locals():
    if qc_type['callrate']:
        sampleqc_in_plink = rules.sample_callRate.output[0]
        sampleqc_in_plink_stem = rules.sample_callRate.params.out
    elif qc_type['variant']:
        sampleqc_in_plink = expand("{{dataout}}/{{sample}}_SnpQc.{ext}",
                                   ext=BPLINK)
        sampleqc_in_plink_stem = rules.snp_qc.params.out
    else:
        sampleqc_in_plink = start['files']
        sampleqc_in_plink_stem = start['stem']

# ---- Prune SNPs, autosome only ----

rule PruneDupvar_snps_noancestry:
    input: sampleqc_in_plink
    output:
        "{dataout}/{sample}_nodup_noancestry.dupvar",
        expand("{{dataout}}/{{sample}}_nodup_noancestry.{ext}",
               ext=['prune.in', 'prune.out'], dataout=DATAOUT),
    params:
        indat = sampleqc_in_plink_stem,
        out = "{dataout}/{sample}_nodup_noancestry"
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} \
  --autosome --indep 50 5 1.5 \
  --list-duplicate-vars --out {params.out}
'''

rule SelectDupvar_snps_noancestry:
    input: rules.PruneDupvar_snps_noancestry.output[0]
    output: "{dataout}/{sample}_nodup_noancestry.dupvar.delete"
    container: 'docker://befh/r_env_gwasamplefilt:3'
    resources:
        mem_mb = 10000,
        time_min = 30
    script: '../scripts/DuplicateVars.R'

# Prune sample dataset
rule sample_prune_noancestry:
    input:
        fileset = sampleqc_in_plink,
        prune = "{dataout}/{sample}_nodup_noancestry.prune.in",
        dupvar = rules.SelectDupvar_snps_noancestry.output
    output:
        temp(expand("{{dataout}}/{{sample}}_pruned_noancestry.{ext}", ext=BPLINK))
    params:
        indat = sampleqc_in_plink_stem,
        out = "{dataout}/{sample}_pruned_noancestry"
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --keep-allele-order --bfile {params.indat} \
  --extract {input.prune} --exclude {input.dupvar} \
  --make-bed --out {params.out}
'''
