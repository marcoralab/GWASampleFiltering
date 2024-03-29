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
            ins = apply_prefix("{dataout}/{sample}_pruned"),
            out = apply_prefix("{dataout}/{sample}_filtered_PCApre")
        resources:
            mem_mb = 10000,
            time_min = 60
        conda: "../envs/plink.yaml"
        shell:
            r'''
plink --keep-allele-order --bfile {params.ins} \
  --remove {input.exclude} \
  --make-bed --out {params.out}
'''

    rule filterKING:
        input:
            king = "{dataout}/{sample}_IBDQC.all.kingfiles",
            exclude = "{dataout}/{sample}_exclude.pca"
        output:
            "{dataout}/{sample}_IBDQC.all.popfilt.kingfiles"
        params:
            indat = apply_prefix("{dataout}/{sample}_IBDQC.all"),
        threads: 6
        resources:
            mem_mb = 10000,
            time_min = 60
        container: 'docker://befh/r_env_gwasamplefilt:5'
        script: '../scripts/filterKing.R'

    rule PCAPartitioning:
        input:
            plink = rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output,
            plink_unpruned = sampleqc_in_plink,
            king = rules.filterKING.output if qc_type['ancestry'] else "{dataout}/{sample}_IBDQC.all.kingfiles",
            iterative = rules.relatedness_sample_fail.output.out
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCApre.{ext}",ext=['unrel', 'partition.log'], dataout = DATAOUT)
        params:
            stem = rules.ancestryFilt.params.out if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            stem_unpruned = apply_prefix(sampleqc_in_plink_stem),
            king = rules.filterKING.params.indat + ".popfilt" if qc_type['ancestry'] else rules.filterKING.params.indat
        threads: 48
        resources:
            mem_mb = 30000,
            walltime = '100:00'
        container: 'docker://befh/genesis_env_gwasamplefilt:3.0'
        script: '../scripts/RunPCAiR.R'

    rule stratFrq:
        input:
            plink = rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output,
            unrel = rules.PCAPartitioning.output[0],
        output: "{dataout}/{sample}_filtered_PCAfreq.frqx"
        params:
            ins = rules.ancestryFilt.params.out if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            out = apply_prefix("{dataout}/{sample}_filtered_PCAfreq")
        resources:
            mem_mb = 10000,
            time_min = 60
        conda: "../envs/plink.yaml"
        shell:
            '''
if [ -f {params.ins}.pcair.eigenval ]; then
  echo "SNPRelate used instead." > {output}
else
  plink --keep-allele-order --bfile {params.ins} --freqx \
    --within {input.unrel} --keep-cluster-names unrelated \
    --out {params.out}
fi
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
            ins = rules.ancestryFilt.params.out if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            out = apply_prefix("{dataout}/{sample}_filtered_PCA")
        resources:
            mem_mb = 10000,
            time_min = 60
        conda: "../envs/plink.yaml"
        shell:
            '''
if [ -f {params.ins}.pcair.eigenval ]; then
  ln -s $(readlink -e {params.ins}.pcair.eigenval) {params.out}.eigenval
  ln -s $(readlink -e {params.ins}.pcair.eigenvec) {params.out}.eigenvec
else
  plink --keep-allele-order --bfile {params.ins} --read-freq {input.frq} --pca 10 \
    --within {input.unrel} --pca-cluster-names unrelated \
    --out {params.out}
fi
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
            ins = apply_prefix("{dataout}/{sample}_pruned"),
            out = apply_prefix("{dataout}/{sample}_filtered_PCA")
        resources:
            mem_mb = 10000,
            time_min = 30
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.ins} --remove {input.exclude} \
  --pca 10 --out {params.out}
'''
else:
    rule PopulationStratification:
        input: rules.ancestryFilt.output if qc_type['ancestry'] else rules.sample_prune_noancestry.output
        output:
            expand("{{dataout}}/{{sample}}_filtered_PCA.{ext}",
                   ext=['eigenval', 'eigenvec'])
        params:
            ins = rules.ancestryFilt.params.out if qc_type['ancestry'] else rules.sample_prune_noancestry.params.out,
            out = "{dataout}/{sample}_filtered_PCA"
        resources:
            mem_mb = 10000,
            time_min = 30
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --keep-allele-order --bfile {params.ins} \
  --pca 10 --out {params.out}
'''
