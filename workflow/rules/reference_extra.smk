BPLINK = ["bed", "bim", "fam"]

def detect_ref_type(reffile):
    if '.vcf' in reffile or '.bcf' in reffile:
        if '{chrom}' in reffile:
            return "vcfchr"
        else:
            return "vcf"
    else:
        return "plink"

ereftype = detect_ref_type(config['extra_ref'])

if ereftype == 'vcfchr':
    rule Reference_prep:
        input:
            vcf = config['extra_ref']['file'],
            tbi = config['extra_ref']['file'] + '.tbi'
        output: temp("{dataout}/extraref.{gbuild}.chr{chrom}.maxmiss{miss}.vcf.gz")
        threads: 12
        resources:
            mem_mb = 4000,
            walltime = '4:00'
        conda: "../envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output}
'''

    rule Reference_cat_extra:
        input:
            vcfs = expand("{dataout}/extraref.{{gbuild}}.chr{chrom}.maxmiss{{miss}}.vcf.gz",
                          chrom = list(range(1, 23)), dataout = DATAOUT)
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        threads: 2
        resources:
            mem_mb = 10000,
            time_min = 30
        conda: "../envs/bcftools.yaml"
        shell:
            '''
bcftools concat {input.vcfs} -Oz -o {output.vcf} --threads 2
bcftools index -ft {output.vcf}
'''

elif ereftype == 'vcf':
    rule Reference_prep_extra:
        input:
            vcf = config['extra_ref']['file'],
            tbi = config['extra_ref']['file'] + '.tbi'
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        threads: 12
        resources:
            mem_mb = 4000,
            time_min = 30
        conda: "../envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''
elif ereftype != 'none': #PLINK fileset of all chromosomes
    # align custom ref to fasta refrence
    rule Ref_Flip_extra:
        input:
            bim = config['extra_ref']['file'] + '.bim',
            bed = config['extra_ref']['file'] + '.bed',
            fam = config['extra_ref']['file'] + '.fam',
            fasta = expand("reference/human_g1k_{gbuild}.fasta", gbuild=BUILD)
        output:
            temp(expand("{{dataout}}/extraref_{{gbuild}}_flipped.{ext}", ext=BPLINK, dataout = DATAOUT))
        resources:
            mem_mb = 10000,
            time_min = 30
        container: 'docker://befh/flippyr:0.5.3'
        shell: "flippyr -p {input.fasta} -o {DATAOUT}/extraref_{wildcards.gbuild}_flipped {input.bim}"

    rule Ref_ChromPosRefAlt_extra:
        input:
            flipped = "{dataout}/extraref_{gbuild}_flipped.bim"
        output:
            bim = temp("{dataout}/extraref_{gbuild}_flipped_ChromPos.bim"),
            snplist = temp("{dataout}/extraref_{gbuild}_flipped_snplist")
        resources:
            mem_mb = 10000,
            time_min = 30
        container: 'docker://befh/r_env_gwasamplefilt:5'
        shell: "R scripts/bim_ChromPosRefAlt.R {input} {output.bim} {output.snplist}"

    # Recode sample plink file to vcf
    rule Ref_Plink2Vcf_extra:
        input:
            bim = rules.Ref_ChromPosRefAlt_extra.output.bim,
            flipped = rules.Ref_Flip_extra.output
        output:
            temp("{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz")
        params:
            out = "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt",
            inp = "{dataout}/extraref_{gbuild}_flipped"
        resources:
            mem_mb = 10000,
            time_min = 30
        conda: "../envs/plink.yaml"
        shell:
            '''
plink --bfile {params.inp} --bim {input.bim} --recode vcf bgz \
  --real-ref-alleles --out {params.out}
'''

    # Index bcf
    rule Ref_IndexVcf_extra:
        input: "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz"
        output: "{dataout}/extraref_{gbuild}_unQC_maxmissUnfilt.vcf.gz.csi"
        resources:
            mem_mb = 10000,
            time_min = 30
        conda: "../envs/bcftools.yaml"
        shell: 'bcftools index -f {input}'

    rule Reference_prep_extra:
        input:
            vcf = rules.Ref_Plink2Vcf_extra.output,
            csi = rules.Ref_IndexVcf_extra.output
        output:
            vcf = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz",
            tbi = "{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz.tbi"
        threads: 12
        resources:
            mem_mb = 4000,
            time_min = 30
        conda: "../envs/bcftools.yaml"
        shell:
            '''
for i in {1..22}; do echo "chr$i $i"; done > reference/chr_name_conv.txt;
bcftools norm -m- {input.vcf} --threads 2 | \
bcftools view -v snps --min-af 0.01:minor -i 'F_MISSING <= {wildcards.miss}' --threads 2 | \
bcftools annotate --rename-chrs reference/chr_name_conv.txt --threads 2 | \
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' --threads 6 -Oz -o {output.vcf}
bcftools index -ft {output.vcf}
'''

union_extraref = (extraref
                  and ("overlap_extra" in config)
                  and (config['overlap_extra'] == 'union'))

rule get_extravars:
    input: '{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.vcf.gz'
    output: '{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.snps'
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/bcftools.yaml"
    shell: "bcftools query -f '%ID\n' {input} > {output}"

rule union_panelvars:
    input:
        ref = expand(
            '{dataout}/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
            gbuild=BUILD, miss=config['QC']['GenoMiss'], refname=REF,
            dataout=DATAOUT),
        eref = expand("{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.snps",
                      gbuild=BUILD, miss=config['QC']['GenoMiss'],
                      dataout=DATAOUT)
    output: '{dataout}/panelvars_all.snp' if union_extraref else "do.not"
    shell: 'cat {input} | sort | uniq > {output}'

rule intersection_panelvars:
    input:
        ref = expand(
            '{dataout}/panelvars_{refname}_{gbuild}_allChr_maxmiss{miss}.snps',
            gbuild=BUILD, miss=config['QC']['GenoMiss'], refname=REF,
            dataout=DATAOUT),
        eref = expand("{dataout}/extraref_{gbuild}_allChr_maxmiss{miss}.snps",
                      gbuild=BUILD, miss=config['QC']['GenoMiss'],
                      dataout=DATAOUT)
    output: '{dataout}/panelvars_all.snp' if not union_extraref else "do.not"
    resources:
        mem_mb = 10000,
        time_min = 30
    conda: "../envs/miller.yaml"
    shell:
        '''
mlr --tsv --implicit-csv-header --headerless-csv-output \
  join -j 1 -f {input[0]} {input[1]} > {output}
'''
