# GWASampleFiltering Pipeline

A Snakemake pipeline to automatically perform best-practice pre-GWAS sample and variant level QC. This sample filtering pipeline is designed to identify potentially outlying, poorly genotyped or mislabeled samples using PLINK formatted, unimputed data. It even makes a pretty report for your reference!

Outputs include a list of samples to exclude (with reasons for exclusion) and an extensive report, with figures and explanations of all operations performed.

## Requirements:

This is a Snakemake pipeline with mostly R and PLINK based rules. It also uses the Flippyr Python package.

 ### The main requirements are as follows:

*   Python 3 (The future is now!)
*   R 3.4 or above
*   PLINK 1.9 (The future is later!)

The following are reccomended:

*   KING

The main package requirements are as follows:

*   Python:
    *   Flippyr (available only using pip, not Anaconda)
    *   Pandas
    *   Numpy
    *   Snakemake
    *   BioPython
*   General R reqirements:
    *   Tidyverse
        *   dplyr
        *   readr
        *   tibble
        *   tidyr
        *   stringr
        *   magrittr
        *   ggplot2
    *   rmarkdown
    *   knitr
    *   plyr


### Pipeline R reqirements:
*   Most components:
    *   kableExtra
    *   DT
*   call rates:
    *   Hmisc
    *   ggforce
    *   Ternary
    *   scales
*   PCA
    *   GENESIS (available from BioConductor)
    *   GWASTools (available from BioConductor)
    *   argparse
    *   GGally
    *   ggforce
    *   plotly
    *   gridExtra
*   Relatedness
    *   Hmisc

## Usage:

### Config files:

There are two config files with settings for imputation prep:

 * *config.yaml* contains settings for the QC and processing.
 * *cluster.yaml* contains settings for cluster execution.

Please review those files before starting.

### Running:

You can run the pipeline with `snakemake` and whatever options you wish to use.

## Pipeline protocol:

#### 1. Basic SNP Level QC:

Perform basic single nucleotide polymorphism (SNP) level QC in PLINK, based on provided cutoffs. This QC increases the power to identify true associations with disease risk by removing suboptimal markers that can increase false positives.

This step calculate allele frequencies, call-rates, and Hardy-Weinberg Equilibrium, saving the outputs.

We then filter the genotypes for other QC based on the maximum genotype missingness (typically < 0.5; ≥ 0.95 call rate), Hardy Weinberg Equilibrium (HWE; default 1e-6, but PLINK recommends 1e-50) and minimum minor allele frequency (MAF; often ≥ 0.01 or ≥ 1%) specified in the configuration file.

Violations of Hardy Weinberg Equilibrium can indicate either the presence of population substructure, or the occurence of genotyping error. It is common practice to assume that violoations are indicative of genotyping error and remove SNPs in which the HWE test statistic is highly significant.

Filtering SNPs on MAF, HWE and call-rate can be done in `PLINK 1.9` by typing the following at the shell prompt:

```bash
plink --bfile [raw-GWA-data] \
  --geno [chosen missingness threshold] \
  --maf [chosen MAF threshold] \
  --keep-allele-order \
  --autosome \
  --hardy \
  --hwe 0.000001 \
  --make-bed --out [filtered-GWA-data]
```

SNP QC is displayed in the final report. By default, this filtered output is used for all other steps, but the setting can be changed in the Snakemake rule file.

#### 2. Sample level call-rate:

Use PLINK to record the call-rate for each sample, and exclude samples with call-rates below the threshold set in the configuration file. A low genotyping call rate in a sample can be indicative of poor DNA sample quality, so samples with a call rate less than specified are excluded from further analysis.

Filtering samples on call-rate (95% here) can be done in `PLINK 1.9` by typing the following at the shell prompt:

```bash
plink --bfile raw-GWA-data  \
  --mind 0.05 \
  --make-bed --out filtered-GWA-data
```

We perform sample-level filtering after SNP level filtering, which is less stringent on the sample level, but accounts for the effect of bad probes on a given array. Not everyone performs QC this way, but we believe it to be best.

Sample call-rate QC is displayed in the final report, and failures are listed in the sample exclusion file, as well as removed in other QC tests.

Any individual failing this basic cutoff is not tested for other quality problems, partially because the tests are often sensitive to bad genotyping.

#### 3. Sex discordance

Use PLINK to find sex-discordant samples. Discordance between self-reported and genetically predicted sex likely indicates errors in sample handling, such as sample swaps. Predicted sex can be determined by calculating X chromosome heterozygosity using an F test, because biological men have one X chromosome and women have two. An F value of ~0.99 in indicates males, and an F value of ~0.03 indicates females. Furthermore, checking X chromosome heterozygosity may revel sex chromosome anomalies (~0.28 in reported females; ~0.35 in males).

Identification of individuals with discordent sex can be done in PLINK 1.9 by typing the following at the shell prompt, which will produce a list of individuals with discordent sex data.

```bash
plink --bfile raw-GWA-data  \
  --check-sex --out --out output.sexcheck
```

Since sex discordance may be due to sample swaps or to incorrect phenotyping, sex discordant samples should generally be removed unless a swap can be reliably resolved.

Sex discordance is displayed in the final report, and listed in the sample exclusion file, but does not affect other QC tests.

#### Interlude: Pruning

Many statistical tests require that all independant variables are independant of each other instead of covarying. This is a problem for genetics where variants are often passed down together during recombination (depending on their distance and other factors) and are correlated.

Pruning selects a set of mostly uncorrelated variants for statistical testing. We prune each cohort and remove duplicate variants for our relatedness, population outlier, population stratification, and heterozygosity tests.

#### 4. Relatedness
Population based cohorts are often limited to unrelated individules as associations statisitcs often assume independence across individules. Closely related samples will share more of their genome and are likely to be more phenotypically similar than than two individules chosen randomly from the population. A common measure of relatedness is identity by descent (IBD), where a kinship correlation coefficent (pi-hat) greater 0.1 suggests that samples maybe related or duplicates samples.

Relatedness QC can be performed here with either PLINK or KING. If relatedness-robust PCA is performed, KING is required.

Identifing duplicated or related samples can be done in `PLINK 1.9` by typing the following command at the shell prompt. This produce a file containing the pariwise IBS estimates for all pairs of individules with minimum IBS of 0.05. Analysis should be performed using an LD pruned snplist.

```
plink --bfile raw-GWA-data \
  --extract snplist.prune.in \
  --genome --min 0.05 --out output.ibd
```

IBD can be calculated using `KING Robust`by typing the following command at the shell prompt. Genotypes do not need to be pruned, but basic QC is recommended.

```
king -b raw-GWA-data.bed \
  --related \
  --degree 3
  --prefix output.ibd
```

To construct a sample of unrelated individules, participants with the highest number of pairwise kinship coefficents > 0.1875 are then iterativly removed.

IBD is displayed in the final report. By default, this output is only used if KING is selected for IBD and PC-AiR is selected for population stratification, but the setting can be changed in the Snakemake rule file.

#### 5. Heterozygosity

Insufficient heterozygosity can indicate inbreeding or other family substructures, while excessive heterozygosity may indicate poor sample quality.

Individuals with outlying heterozygosity rates can be identified in PLINK 1.9 by typing the following command at the shell prompt:

```
plink --bfile raw-GWA-data \
  --extract snplist.prune.in \
  --het --out output.het
```

This produces a file containing Method-of-moments F coefficient estimates, which can be used to calculate the observed heterozygosity rate in each individual. Analysis is performed using an LD pruned snplist.

We calculate a heterozygocity similarly using observed and expected counts from the PLINK output [(Observed - Expected)/N) and exclude samples that are ± 3 sd from the cohort mean.

#### 6. Population outliers

Population stratification occurs when the study population under investigation comprises several different subpopulations that differ in both genetic ancestry and in the phenotype of interest. Spurious associations can result from genetic ancestry rather than true associations of alleles with the phenotype. Principal component analysis (PCA) can be used to identify population outliers by perfoming a PCA in a reference panel such as 1000 genomes and projecting the sample of interest onto the resulting space.

To perform a PCA analysis in PLINK, we first obtain genotype data and pedigree infromation for the 1000 genomes reference panal from <http://www.internationalgenome.org/data>. Then we remove related individuals (children, siblings, second degree and other, keeping third degree relatives), and merged the reference with our sample. The alleles in the sample dataset are aligned to the same DNA strand as the reference dataset to allow the datasets to be merged correctly.

Two additional files are required, one listing the FID, IID and population for each sample in the reference dataset and the second listing the population clusters in the reference data. We generate the first automatically and provide the second.

PLINK calculates PCs for the 1000 genomes reference, then projects the samples onto those components, so related samples are not a problem.

Both datasets were LD pruned to eliminate a large degree of the redundency in the data and reduce the influence of potential chromsomal artifacts.

The following `PLINK 1.9` commands for performing the PCA analysis can be entered at the shell prompt to generate two files containing the principal component eigenvalues and Principal component eigenvectors.

```bash
plink --bfile merged-reference-sample-data \
  --pca 10 --within sample_population.txt \
  --pca-clusters population_clusters.txt \
  --out pca.output
```

We remove any samples not within 6 standard deviations of the chosen superpopulation on all 10 PCs.

Population outliers are displayed in the final report, listed in the sample exclusion file, and removed for population stratification.

#### 6. Population stratification
After excluding population outliers from the dataset, population substructure will remain due to the presence of genetic diversity within apparently homogenous populations. Within a single ethnic population, even subtle degrees of population stratification can bias results due to differences in allele frequencies between subpopulations. Principal components based on the observed genotypes in the dataset of interest can be used to capture information on substructure and be included as covariates in downstream analysis.

To obtain the principal components for the sample dataset after population outliers have been removed, type the following `PLINK 1.9` commands at the shell prompt to generate the principal component eigenvalue and eigenvector files.

```bash
plink --bfile raw-GWA-data \
  --remove fail-ancestry-QC.txt \
  --pca 10 \
  --out filter-GWA-data
```

However, ancestry can bias the dimensionality reduction, and PC's may represent relatedness as well as gross population structure. To account for this, we have an option in the configuration (on by default) to use a similar algorithm to PC-AiR, from the GENESIS R package.

First, we partition the samples into related and unrelated sets. By default, we use PC-AiR for this, but if it fails, we use the iteratively generated set from our relatedness QC.

Then we use PLINK to calculate allele frequencies for PCA weighting in the unrelated sample only, then generate PCs in plink, with that weighting, projecting the related samples onto the PCA transform for the unrelated samples. This will allow the experimenter to choose which QC filters to use and still have access to PCs for stratification.
