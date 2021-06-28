
message("Render Final Report",
        "\n markdown: ", snakemake@input[["markdown"]],
        "\n SexFile: ", snakemake@input[["SexFile"]],
        "\n hwe: ", snakemake@input[["hwe"]],
        "\n frq: ", snakemake@input[["frq"]],
        "\n frqx: ", snakemake@input[["frqx"]],
        "\n imiss: ", snakemake@input[["imiss"]],
        "\n HetFile: ", snakemake@input[["HetFile"]],
        "\n IBD_stats: ", snakemake@input[["IBD_stats"]],
        "\n PCA_rdat: ", snakemake@input[["PCA_rdat"]],
        "\n PopStrat_eigenval: ", snakemake@input[["PopStrat_eigenval"]],
        "\n PopStrat_eigenvec: ", snakemake@input[["PopStrat_eigenvec"]],
        "\n partmethod: ", snakemake@input[["partmethod"]],
        "\n output: ", snakemake@output[[1]],
        "\n rwd: ", snakemake@params[["rwd"]],
        "\n Family: ", snakemake@params[["Family"]],
        "\n pi_threshold: ", snakemake@params[["pi_threshold"]],
        "\n output_dir: ", snakemake@params[["output_dir"]],
        "\n geno_miss: ", snakemake@params[["geno_miss"]],
        "\n MAF: ", snakemake@params[["MAF"]],
        "\n HWE: ", snakemake@params[["HWE"]],
        "\n superpop: ", snakemake@params[["superpop"]],
        "\n partmethod: ", snakemake@params[["partmethod"]]
      )

      rule GWAS_QC_Report:
          input:
              script = "scripts/GWAS_QC.Rmd",
              SexFile = "{dataout}/{sample}_SexQC.sexcheck" if do_sexqc else '/dev/null',
              hwe = "{dataout}/{sample}_SnpQc.hwe",
              frq = "{dataout}/{sample}_SnpQc.frq",
              frqx = "{dataout}/{sample}_SnpQc.frqx",
              imiss = "{dataout}/{sample}_callRate.imiss",
              HetFile = "{dataout}/{sample}_HetQC.het",
              IBD_stats = "{dataout}/{sample}_IBDQC.Rdata",
              PCA_rdat = "{dataout}/{sample}_pca.Rdata",
              PopStrat_eigenval = "{dataout}/{sample}_filtered_PCA.eigenval",
              PopStrat_eigenvec = "{dataout}/{sample}_filtered_PCA.eigenvec",
              partmethod = rules.PCAPartitioning.output[1] if config["pcair"] else "/dev/null"
          output:
               "{dataout}/stats/{sample}_GWAS_QC.html"
          params:
              rwd = RWD,
              Family = FAMILY,
              pi_threshold = 0.1875,
              output_dir = "{dataout}/stats",
              idir = "{dataout}/stats/md/{sample}",
              geno_miss = config['QC']['GenoMiss'],
              samp_miss = config['QC']['SampMiss'],
              MAF = config['QC']['MAF'],
              HWE = config['QC']['HWE'],
              superpop = config['superpop'],
              partmethod = rules.PCAPartitioning.output[1] if config["pcair"] else "none"


rmarkdown::render(
  input = snakemake@input[["markdown"]],
  clean = TRUE,
  intermediates_dir = snakemake@params[["output_dir"]],
  output_file = snakemake@output[[1]],
  output_dir = snakemake@params[["output_dir"]],
  output_format = "all",
  params = list(
    rwd = snakemake@params[['rwd']],
    path = snakemake@params[["path"]],
    chrom = snakemake@params[["chrom"]],
    cohort = snakemake@params[["cohort"]],
    maf = snakemake@params[["maf"]],
    rsq = snakemake@params[["rsq"]],
    rsq2 = snakemake@params[["rsq2"]],
    sampsize = snakemake@params[["sampsize"]],
    # out = snakemake@params[["out"]],
    outpath = snakemake@params[["output_dir"]]
  )
)
