# ---- Exclude Samples with interealtedness ----
rule relatedness_sample_prep:
    input: sampleqc_in_plink
    output:
        bed = temp("{dataout}/{sample}_IBDQCfilt.bed"),
        bim = temp("{dataout}/{sample}_IBDQCfilt.bim"),
        fam = temp("{dataout}/{sample}_IBDQCfilt.fam")
    params:
        ins = apply_prefix(sampleqc_in_plink_stem),
        out = apply_prefix("{dataout}/{sample}_IBDQCfilt")
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --bfile {params.ins} \
  --geno 0.02 \
  --maf 0.02 \
  --memory 6000 \
  --make-bed --out {params.out}
'''

rule relatedness_sample_prep_remnoncc:
    input:
        bed = rules.relatedness_sample_prep.output.bed,
        bim = rules.relatedness_sample_prep.output.bim,
        fam = rules.relatedness_sample_prep.output.fam
    output:
        bed = temp("{dataout}/{sample}_IBDQCfilt_rm-non-cc.bed"),
        bim = temp("{dataout}/{sample}_IBDQCfilt_rm-non-cc.bim"),
        fam = temp("{dataout}/{sample}_IBDQCfilt_rm-non-cc.fam")
    resources:
        mem_mb = 10000,
        time_min = 30
    shell:
        '''
cp {input.bed} {output.bed}
cp {input.bim} {output.bim}
awk 'BEGIN {{OFS="\t"}} $6 !~ "^[12]$" {{$6=-9}}1{{print $1,$2,$3,$4,$5,$6}}' \
  {input.fam} > {output.fam}
'''

rule relatedness_QC:
    input:
        bed = rules.relatedness_sample_prep_remnoncc.output.bed,
        bim = rules.relatedness_sample_prep_remnoncc.output.bim,
        fam = rules.relatedness_sample_prep_remnoncc.output.fam
    output: "{dataout}/{sample}_IBDQC.kingfiles"
    params:
        out = apply_prefix("{dataout}/{sample}_IBDQC"),
        dataout = apply_prefix(DATAOUT)
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = '144:00'
    conda: "../envs/king.yaml"
    shell:
        '''
king -b {input.bed} --related --degree 4 --prefix {params.out} --cpus 24 > {params.out}.log
if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.kin*")"; then
  find {params.dataout} -name "{wildcards.sample}_IBDQC.kin*" > {output}
elif grep --quiet "No close relatives" {params.out}.log; then
  echo norel > {output}
fi
'''

rule king_all:
    input:
        bed = rules.relatedness_sample_prep_remnoncc.output.bed,
        bim = rules.relatedness_sample_prep_remnoncc.output.bim,
        fam = rules.relatedness_sample_prep_remnoncc.output.fam
    output: "{dataout}/{sample}_IBDQC.all.kingfiles",
    params:
        out = apply_prefix("{dataout}/{sample}_IBDQC.all"),
        dataout = apply_prefix(DATAOUT)
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = '144:00'
    conda: "../envs/king.yaml"
    shell:
        '''
king -b {input.bed} --kinship --ibs --prefix {params.out} --cpus 24 > {params.out}.log
if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.all.kin*")"; then
  find {params.dataout} -name "{wildcards.sample}_IBDQC.all.kin*" > {output}
fi
'''

rule SampleExclusion_prerelate:
    input:
        SampCallRate = "{dataout}/{sample}_callRate.irem" if qc_type['callrate'] else ancient('/dev/null'),
        het = "{dataout}/{sample}_exclude.heterozigosity" if qc_type['heterozygosity'] else ancient('/dev/null'),
        sex = "{dataout}/{sample}_exclude.sexcheck" if qc_type['sex'] else ancient('/dev/null'),
        pca = "{dataout}/{sample}_exclude.pca" if qc_type['ancestry'] else ancient('/dev/null'),
        relat = ancient('/dev/null')
    output:
        out = temp("{dataout}/{sample}_exlude.prerelate"),
        out_distinct = temp("{dataout}/{sample}_exlude.prerelate.distinct")
    resources:
        mem_mb = 10000,
        time_min = 30
    container: 'docker://befh/r_env_gwasamplefilt:5'
    script: "../scripts/sample_QC.R"

rule relatedness_sample_fail:
    input:
        genome = rules.relatedness_QC.output,
        geno_all = rules.king_all.output,
        fam = sampleqc_in_plink_stem + ".fam",
        exclude = rules.SampleExclusion_prerelate.output.out if config['relatedness_preexclude'] else ancient('/dev/null')
    params:
        Family = FAMILY,
        threshold = 0.1875,
        geno = rules.relatedness_QC.params.out
    output:
        out = "{dataout}/{sample}_exclude.relatedness",
        allrelate = "{dataout}/{sample}_exclude.relatedness_all" if FAMILY == 'T' else [],
        rdat = "{dataout}/{sample}_IBDQC.Rdata"
    threads: 8
    resources:
        mem_mb = 64000,
        walltime = '24:00'
    container: 'docker://befh/r_env_gwasamplefilt:7'
    script: '../scripts/relatedness_QC.R'
