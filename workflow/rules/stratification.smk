'''Snakefile for GWAS Variant and Sample QC Version 0.3.1'''

full_pipeline = config['full_pipeline'] if 'full_pipeline' in config else False

if not config['full_pipeline']:
    from scripts.parse_config_GWASampleFiltering import parser

    configfile: "config/config.yaml"

    shell.executable("/bin/bash")

    BPLINK = ["bed", "bim", "fam"]

    start, SAMPLE, DATAOUT = parser(config)

    rule all:
        input:
            expand("{dataout}/{sample}_filtered_PCA.{ext}",
                   ext=['eigenval', 'eigenvec'], sample=SAMPLE, dataout=DATAOUT)


    def enableqc(qc_type):
        if 'qc' in config:
            if (isinstance(config['qc'], list)
                    and all([isinstance(x, str) for x in config['qc']])):
                return qc_type in config['qc']
            elif isinstance(config['qc'], dict):
                if (qc_type in config['qc']
                        and isinstance(config['qc'][qc_type], bool)):
                    return config['qc'][qc_type]
                else:
                    raise Exception(
                        "Malformed QC dict: All supported QC must be present")
            else:
                raise Exception("Malformed QC list: Please provide dict OR list")
        else:
            return True


    qc_type = ['variant', 'callrate', 'popstrat', 'ancestry']
    qc_type = {x: enableqc(x) for x in qc_type}

    if config['pcair']:
        include: 'relatedness.smk'


if qc_type['ancestry']:
    if not config['full_pipeline']:
        include: 'ancestry.smk'
else:
    include: 'prune_no-refmerge.smk'

# Run PCA to control for population stratification
if config['pcair']:
    assert config['king'], "You must use KING for relatedness if using PCAiR!"
    rule ancestryFilt:
        input:
            plink = expand("{{dataout}}/{{sample}}_pruned.{ext}", ext=BPLINK),
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            temp(expand("{{dataout}}/{{sample}}_filtered_PCApre.{ext}",
                 ext=BPLINK, dataout = DATAOUT)),
        params:
            indat = "{dataout}/{sample}_pruned",
            plinkout = "{dataout}/{sample}_filtered_PCApre"
        conda: "../envs/plink.yaml"
        shell:
            r'''
plink --keep-allele-order --bfile {params.indat} \
  --remove {input.exclude} \
  --make-bed --out {params.plinkout}
'''

    rule filterKING:
        input:
            king = "{dataout}/{sample}_IBDQC.all.kingfiles",
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            "{dataout}/{sample}_IBDQC.all.popfilt.kingfiles"
        params:
            indat = "{dataout}/{sample}_IBDQC.all",
        container: 'docker://befh/r_env_gwasamplefilt:2'
        script: '../scripts/filterKing.R'

    rule PCAPartitioning:
        input:
            plink = rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output,
            king = rules.filterKING.output if qc_type['ancestry'] else "{dataout}/{sample}_IBDQC.all.kingfiles",
            iterative = rules.relatedness_sample_fail.output.out
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCApre.{ext}",ext=['unrel', 'partition.log'], dataout = DATAOUT)
        params:
            stem = rules.ancestryFilt.params.plinkout if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            king = rules.filterKING.params.indat + ".popfilt" if qc_type['ancestry'] else rules.filterKING.params.indat
        container: 'docker://befh/genesis_env_gwasamplefilt:2.1'
        script: '../scripts/PartitionPCAiR.R'

    rule stratFrq:
        input:
            plink = rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output,
            unrel = rules.PCAPartitioning.output[0],
        output: "{dataout}/{sample}_filtered_PCAfreq.frqx"
        params:
            indat = rules.ancestryFilt.params.plinkout if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            out = "{dataout}/{sample}_filtered_PCAfreq"
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} --freqx \
  --within {input.unrel} --keep-cluster-names unrelated \
  --out {params.out}
'''

    rule PopulationStratification:
        input:
            plink = rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output,
            unrel = rules.PCAPartitioning.output[0],
            frq = rules.stratFrq.output
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCA.{ext}",
                   ext=['eigenval', 'eigenvec'])
        params:
            indat = rules.ancestryFilt.params.plinkout if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            out = "{dataout}/{sample}_filtered_PCA"
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} --read-freq {input.frq} --pca 10 \
  --within {input.unrel} --pca-cluster-names unrelated \
  --out {params.out}
'''
elif qc_type['ancestry']:
    rule PopulationStratification:
        input:
            plink = expand("{{dataout}}/{{sample}}_pruned.{ext}", ext=BPLINK),
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCA.{ext}",
                   ext=['eigenval', 'eigenvec'])
        params:
            indat = "{dataout}/{sample}_pruned",
            out = "{dataout}/{sample}_filtered_PCA"
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} --remove {input.exclude} \
  --pca 10 --out {params.out}
'''
else:
    rule PopulationStratification:
        input: rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCA.{ext}",
                   ext=['eigenval', 'eigenvec'])
        params:
            indat = rules.ancestryFilt.params.plinkout if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            out = "{dataout}/{sample}_filtered_PCA"
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.indat} \
  --pca 10 --out {params.out}
'''
