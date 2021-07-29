# ---- Exclude Samples with interealtedness ----
rule relatedness_sample_prep:
    input: sampleqc_in_plink
    output:
        bed = temp("{dataout}/{sample}_IBDQCfilt.bed"),
        bim = temp("{dataout}/{sample}_IBDQCfilt.bim"),
        fam = temp("{dataout}/{sample}_IBDQCfilt.fam")
    params:
        indat_plink = sampleqc_in_plink_stem,
        out = "{dataout}/{sample}_IBDQCfilt"
    conda: "../envs/plink.yaml"
    shell:
        '''
plink --bfile {params.indat_plink} \
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
    params:
        indat_plink = rules.relatedness_sample_prep.params.out,
        out = "{dataout}/{sample}_IBDQCfilt_rm-non-cc"
    conda: "../envs/plink.yaml"
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
    output:
        genome = "{dataout}/{sample}_IBDQC.kingfiles"
    params:
        out = "{dataout}/{sample}_IBDQC",
        dataout = DATAOUT
    conda: "../envs/king.yaml"
    shell:
        '''
king -b {input.bed} --related --degree 4 --prefix {params.out}
if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.kin*")"; then
  find {params.dataout} -name "{wildcards.sample}_IBDQC.kin*" > {output.genome}
fi
'''

rule king_all:
    input:
        bed = rules.relatedness_sample_prep_remnoncc.output.bed,
        bim = rules.relatedness_sample_prep_remnoncc.output.bim,
        fam = rules.relatedness_sample_prep_remnoncc.output.fam
    output:
        genome = "{dataout}/{sample}_IBDQC.all.kingfiles",
    params:
        out = "{dataout}/{sample}_IBDQC.all",
        dataout = DATAOUT
    conda: "../envs/king.yaml"
    shell:
        '''
king -b {input.bed} --kinship --ibs --prefix {params.out}
if test -n "$(find {params.dataout} -name "{wildcards.sample}_IBDQC.all.kin*")"; then
  find {params.dataout} -name "{wildcards.sample}_IBDQC.all.kin*" > {output.genome}
fi
'''

rule relatedness_sample_fail:
    input:
        genome = rules.relatedness_QC.output.genome,
        geno_all = rules.king_all.output if config['king'] else "/dev/null",
        fam = sampleqc_in_plink_stem + ".fam"
    params:
        Family = FAMILY,
        king = config['king'],
        threshold = 0.1875,
        geno = rules.relatedness_QC.params.out if config['king'] else rules.relatedness_QC.output.genome
    output:
        out = "{dataout}/{sample}_exclude.relatedness",
        rdat = "{dataout}/{sample}_IBDQC.Rdata"
    conda: "../envs/r.yaml"
    script: '../scripts/relatedness_QC.R'
