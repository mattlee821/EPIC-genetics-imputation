# EPIC Genetics Pipeline

This repository contains the EPIC genetics imputation pipeline. The pipeline is composed of 4 stages, each stage uses an independent Nextflow workflow:

1. `pipeline_stage1/`: study-specific raw [genotype](https://en.wikipedia.org/wiki/Genotype) preprocessing and liftover to hg38
2. `pipeline_stage2/`: [phasing](https://en.wikipedia.org/wiki/Haplotype_estimation) and [imputation](https://en.wikipedia.org/wiki/Imputation_(genetics)) against the [1000 Genomes Project](https://en.wikipedia.org/wiki/1000_Genomes_Project) high-coverage GRCh38 reference
3. `pipeline_stage3/`: post-imputation QC, [rsID](https://www.ncbi.nlm.nih.gov/snp/docs/RefSNP_about/) annotation, and conversion to PLINK2 format
4. `pipeline_stage4/`: cross-stage master report generation and final study archive

Stage-specific documentation:

- [pipeline_stage1/README.md](pipeline_stage1/README.md)
- [pipeline_stage2/README.md](pipeline_stage2/README.md)
- [pipeline_stage3/README.md](pipeline_stage3/README.md)
- [pipeline_stage4/README.md](pipeline_stage4/README.md)

Some code/scripts/functions and processes used in this pipeline are direct copies or adaptations from prior EPIC genetics pipelines by: Joshua Atkins, Marie Breeur, Aurelie Gabriel, Emilie Gerard-Marchant, Manon Knuchel...

## Stage 1: Study Harmonization and Genome-Build Standardization

- Directory: `pipeline_stage1/`
- Main purpose: Convert each raw study into a harmonized hg38 PLINK handoff
- Primary input:  raw study data plus manifests and EPIC ID linkage files
- Primary output: `analysis/<STUDY>/stage1/<STUDY>.bed/.bim/.fam`
- Submission script: `src/004_stage1.sh` 

Stage 1 converts each raw study delivery into a common hg38 PLINK handoff. Each study is processed with its own bespoke script because the original deliveries differ in array platform, allele coding, chromosome naming, manifest structure, and identifier linkage rules. All downstream matching is coordinate-based and assumes one common genome build. All downstream tools expect each study to exist as one standardized PLINK handoff.

#### What this stage does:

1. harmonizes sample identifiers to the EPIC identifier system
2. applies study-specific preprocessing such as BIM sorting, AB-allele translation, sex-chromosome handling, and build-specific logic
3. standardizes all studies onto GRCh38
4. writes one consistent study-level PLINK dataset for the next stage

#### Data used:

- raw study genotype files from source
- study-specific chip manifests and metadata files
- study-specific EPIC linkage files
- the shared EPIC ID reference
- LiftOver resources used for build conversion to GRCh38

## Stage 2: Phasing and Imputation

- Directory: `pipeline_stage2/`
- Main purpose: Export stage-1 PLINK data to per-chromosome target VCFs, phase, and impute
- Primary input:  `analysis/<STUDY>/stage1/` plus 1000G GRCh38 reference and SHAPEIT5 genetic map
- Primary output: `analysis/<STUDY>/stage2/<STUDY>_chr*_GxS.imputed.vcf.gz`
- Submission script: `src/005_stage2.sh` 

#### What this stage does:
1. performs pre-phasing sample and variant QC
2. phases autosomes and chrX against the 1000G reference
3. imputes each phased chromosome against the 1000G reference panel
4. generates per-study phasing and imputation reports

### Phasing

[Phasing](https://en.wikipedia.org/wiki/Haplotype_estimation) estimates which [alleles](https://en.wikipedia.org/wiki/Allele) sit together on the same physical [chromosome](https://en.wikipedia.org/wiki/Chromosome) copy. Standard [GWAS](https://en.wikipedia.org/wiki/Genome-wide_association_study) array data tell us which two alleles an individual has at a locus, but not which allele belongs to which parental [haplotype](https://en.wikipedia.org/wiki/Haplotype). Phasing resolves that ambiguity.

#### Why:

- imputation works best on phased haplotypes rather than unordered diploid genotypes
- phased target data align more naturally to phased reference haplotypes
- local haplotype structure carries more information than isolated variants

#### How:

- stage-1 PLINK files are first passed through a pre-phasing QC step (STAGE1_QC) that applies sex checks, heterozygosity filtering, duplicate removal, HWE, and MAF filters
- the QC-cleaned data are then exported to per-chromosome target VCFs
- SHAPEIT5 phases each study chromosome against the matching public reference chromosome
- chrX is handled separately because `PAR1`, `nonPAR`, and `PAR2` require different treatment; see [here](https://github.com/statgen/Minimac4/issues/74)

#### Data used:

- stage-1 hg38 PLINK files, exported to chromosome-specific target VCFs
- the public 1000 Genomes high-coverage GRCh38 phased reference panel
	- the 1000 Genomes panel provides the external phased haplotypes that SHAPEIT5 uses to infer study haplotypes
- the SHAPEIT5 [genetic recombination map](https://en.wikipedia.org/wiki/Genetic_map)
	- the SHAPEIT5 map supplies recombination information that helps decide where haplotypes are likely to switch

### Imputation

[Imputation](https://en.wikipedia.org/wiki/Genotype_imputation) statistically infers [genotypes](https://en.wikipedia.org/wiki/Imputation_(genetics)) at variants that were not directly typed on the original array.

#### Why:

- different studies were genotyped on different arrays, so raw variant sets are not identical
- many informative variants are absent from one or more arrays
- imputation greatly increases genomic coverage and improves cross-study comparability

#### How:

- Minimac4 compares each phased study haplotype to the corresponding reference haplotypes
- when the study haplotypes closely match the reference haplotypes, untyped nearby alleles can be inferred
- the final output is an imputed per-chromosome VCF with dosage-style genotype information and quality metrics such as `R2`

#### Data used:

- phased study VCFs from the phasing step
- the same public 1000 Genomes GRCh38 reference panel, compressed to Minimac4 `msav` format
	- the 1000 Genomes panel is the public reference source that makes imputation possible
- the GRCh38 [reference genome](https://en.wikipedia.org/wiki/Reference_genome) FASTA used during normalization and reference preparation
	- the GRCh38 FASTA ensures that REF and ALT alleles are interpreted on the correct genome sequence

## Stage 3: Sample and variant QC, variant annotation, and finalisation

- Directory: `pipeline_stage3/`
- Main purpose: Annotate rsIDs, filter imputed variants, convert to PLINK2, and perform sample QC
- Primary input: `analysis/<STUDY>/stage2/` plus [dbSNP](https://en.wikipedia.org/wiki/DbSNP) GRCh38
- Primary output: `analysis/<STUDY>/stage3/final/<STUDY>_chr*.pgen/.pvar/.psam`
- Submission script: `src/006_stage3.sh`

Stage 3 turns the chromosome-level imputed VCFs from stage 2 into analysis-ready per-chromosome PLINK2 datasets. It combines four linked operations into one workflow:

#### What this stage does:

1. variant identifier standardization
2. variant-level QC (including HWE)
3. conversion to per-chromosome PLINK2 format
4. sample-level QC (ancestry, relatedness, heterozygosity)
5. generate per-study QC reports

#### Why:

- stage-2 outputs intentionally use coordinate-style IDs, but many downstream resources still expect rsIDs where available
- imputed data need an explicit post-imputation QC layer before they are treated as analysis-ready
- chromosome-level imputed VCFs are not the most practical final format for downstream study analyses
- sample-level QC is needed to identify technical problems, unexpected relatedness, and ancestry outliers

#### How:

- each imputed chromosome is annotated against dbSNP
	- if an exact GRCh38 match exists, the variant receives its rsID
	- if no exact dbSNP match exists, the pipeline uses a fallback `CHR:POS:REF:ALT` identifier
	- a mapping file is written so the identifier transition remains auditable
- the annotated variants are then filtered on imputation quality using `INFO/R2`
- the data are also filtered on [minor allele frequency](https://en.wikipedia.org/wiki/Minor_allele_frequency)
- [Hardy-Weinberg equilibrium](https://en.wikipedia.org/wiki/Hardy%E2%80%93Weinberg_principle) filtering is applied by default on autosomes only; failing variants are written to a per-study `hwe.exclude` list for downstream use (no hard filter is applied to the final dataset)
- the filtered per-chromosome VCFs are converted into PLINK2 `pgen/pvar/psam` files; each chromosome is finalized independently — no genome-wide merge is performed
- sample QC uses a LD-pruned autosomal subset derived in parallel across all 22 autosomes; the per-chromosome pruned datasets are merged into a single small BED file for relatedness and PCA computation only
- sample QC estimates relatedness with KING, identifies [heterozygosity](https://en.wikipedia.org/wiki/Heterozygosity) outliers, and runs [PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)-based ancestry QC
- recorded sex is updated from the stage-1 `.fam` file during per-chromosome finalization
- a final sample-removal list is built from the required QC exclusions; ancestry outliers are identified but not removed by default
- Stage 3-only tables, figures, flags, and HTML reports are written under `analysis/<STUDY>/stage3/report/`

#### Data used:

- stage-2 imputed VCFs and consolidated metrics
- [dbSNP](https://en.wikipedia.org/wiki/DbSNP) GRCh38
	- dbSNP provides rsIDs where a known GRCh38 variant match exists
- the stage-1 `.fam` file as the authoritative recorded-sex source
	- the stage-1 `.fam` anchors sample sex to the original study handoff
- LD-pruned autosomal subsets derived per-chromosome in parallel and merged for relatedness and PCA-based QC
	- the merged LD-pruned autosomal subset provides a stable basis for relatedness, heterozygosity, and ancestry QC

## Stage 4: Master Report and Finalisation

- Directory: `pipeline_stage4/`
- Main purpose: Generate per-study cross-stage master HTML reports and build the final deliverable archive
- Primary input: `analysis/<STUDY>/stage3/` outputs plus stage 2 HTML reports
- Primary output: `final/<STUDY>/` with archive, master HTML, and QC review files
- Submission script: `src/007_stage4.sh`

Stage 4 assembles the outputs from all preceding stages into a single deliverable per study. It runs three processes:

#### What this stage does:

1. generate a cross-stage master HTML report integrating stage 2 and stage 3 QC summaries
2. build a per-study archive containing the finalized PLINK2 files, QC exclude lists, and all reports
3. write a combined `sample-manifest.tsv` listing every sample across all studies

#### How:

- `MASTER_REPORT` reads the stage 2 and stage 3 report trees from `analysis_root` and writes a single master HTML per study
- `FINALISE_STUDY` collects the per-chromosome PGEN files from `stage3/final/`, the QC exclude lists from `stage3/report/flags/`, and the stage 2 and stage 3 HTML reports, then packages them into a deliverable tarball under `final/<STUDY>/`
- `SAMPLE_MANIFEST` reads the finalised per-study `.psam` files, joins their sample IIDs to the `Idepic`/`Idepic_Bio` map that `003-data-epic.R` exports from `genetics_id.sas7bdat`, and writes a single `final/sample-manifest.tsv` — a shareable index of which samples exist in the genetics data (see [Sample manifest](#sample-manifest) below for the columns and usage)

#### Data used:

- per-chromosome PLINK2 files from stage 3 finalization
- HWE, relatedness, and ancestry exclude lists from stage 3 QC
- stage 2 imputation HTML report
- stage 3 sample and variant QC HTML report

## Output

All pipeline outputs are written to a date-labelled run directory on the scratch filesystem:

```
${SCRATCH}/${SCRATCH_DATE}/
├── studies/                          # per-study stage outputs
│   └── <STUDY>/
│       ├── stage1/                   # stage 1 PLINK handoffs and reports
│       ├── stage2/                   # stage 2 imputed VCFs and reports
│       └── stage3/
│           ├── final/                # per-chromosome PLINK2 pgen/pvar/psam files
│           ├── report/               # QC tables, figures, flags, HTML report, manifests
│           └── sample_review/        # het, king, and PCA files
├── stage1/
│   ├── work/                         # pipeline-level work root (stage 1)
│   ├── nxf-work/                     # Nextflow task work directory (stage 1)
│   └── .nextflow/                    # Nextflow cache and history (stage 1)
├── stage2/
│   ├── work/                         # Nextflow task work directory (stage 2)
│   ├── conda/                        # Nextflow conda environment cache (stage 2)
│   └── .nextflow/                    # Nextflow cache and history (stage 2)
├── stage3/
│   ├── work/                         # Nextflow task work directory (stage 3)
│   ├── conda/                        # Nextflow conda environment cache (stage 3)
│   └── .nextflow/                    # Nextflow cache and history (stage 3)
├── stage4/
│   ├── work/                         # Nextflow task work directory (stage 4)
│   ├── conda/                        # Nextflow conda environment cache (stage 4)
│   └── .nextflow/                    # Nextflow cache and history (stage 4)
└── final/                            # finalised outputs (007_stage4.sh)
    ├── <STUDY>/
    │   ├── <STUDY>.stage3.tar.gz     # per-chromosome PLINK2 pfiles + QC exclude lists
    │   ├── report-stage2.html
    │   ├── report-stage3.html
    │   ├── report-master.html
    │   └── review/                   # het, relatedness, and PCA files
    ├── sample-manifest.tsv           # all samples: study, IID, Idepic, Idepic_Bio, sex, pheno, overlap_keep
    └── summaries/                    # stage1/2/3 summary markdown files
```

Stage-level summary files are written directly into `studies/`:

- `studies/stage1-summary.md`
- `studies/stage2-summary.md`
- `studies/stage3-summary.md`
- `studies/<STAGE>-pipeline_info/` — Nextflow trace and execution report for each stage

### Strand and allele alignment

The finalised data carry **no residual per-variant strand ambiguity**. Strand is
resolved twice upstream: Stage 1 determines each array SNP's strand from the array
manifest and flips it to a consistent orientation before lifting over, and Stage 2
conforms every study to the 1000 Genomes GRCh38 reference panel before phasing and
imputation. As a result, **all variants and alleles are reported on the forward (+)
strand of GRCh38, with REF/ALT defined relative to the 1000 Genomes reference panel.**

When using these data in a GWAS, describe the strand accordingly — e.g.:

> Genotypes were harmonised to the GRCh38 forward strand during pre-imputation QC and
> imputed against the 1000 Genomes GRCh38 reference panel (SHAPEIT5 phasing, Minimac4
> imputation). All variants and effect estimates are reported on the forward strand of
> GRCh38, with alleles defined relative to the reference panel.

Because every variant sits on a single, unambiguous reference, palindromic (A/T, C/G)
SNPs are not a problem *within* this dataset. They only become ambiguous when aligning
to an **external** dataset of unknown strand (meta-analysis, cross-cohort merging, or
PRS with external weights). To keep that alignment safe:

- publish summary statistics with the **effect allele, other allele, and effect-allele
  frequency (EAF)** — the EAF is what lets consumers realign palindromic SNPs;
- apply the palindromic exclusion (drop A/T and C/G SNPs with MAF > ~0.4, where allele
  frequency cannot resolve a strand flip) only **at the point of external harmonisation**,
  not to the base data, since a same-panel or within-cohort analysis loses nothing by
  keeping them.

### Sample manifest

`final/sample-manifest.tsv` is a single tab-separated index of every sample in the
finalised genetics data — one row per (study, sample) — for sharing with collaborators
so they can see which of their participants have genetics data and build a keep-list for
their own analysis. It contains sample identifiers only; no genotypes.

| column | meaning |
|---|---|
| `study` | study the sample belongs to (e.g. `Brea_02`) |
| `IID` | sample ID **exactly as stored in that study's PLINK2 files** — use this to subset genotypes. It is a composite `<Idepic>_<Idepic_Bio>` (or a single `<Idepic_Bio>`) |
| `Idepic` | EPIC participant (person) ID, from `genetics_id.sas7bdat` |
| `Idepic_Bio` | EPIC biosample ID, from `genetics_id.sas7bdat`. **This is the key to match against EPIC source data**; the same participant has the same `Idepic_Bio` across studies |
| `sex` | PLINK sex code: `1` = male, `2` = female, `0`/`NA` = unknown |
| `pheno` | PLINK phenotype as stored: `1` = control, `2` = case, `NA` = missing |
| `overlap_keep` | de-duplication flag for participants genotyped in **more than one study**: `TRUE` = keep this copy (the study with the smallest total N), `FALSE` = drop this copy (duplicate held in a larger study), `NA` = participant is in this one study only. Keyed on `Idepic_Bio`; keeping `TRUE` + `NA` yields exactly one genetics sample per participant |

Both `Idepic` and `Idepic_Bio` are looked up from `genetics_id.sas7bdat` (via the map
`003-data-epic.R` writes), not parsed from the `IID`, so they are authoritative even for
samples whose `IID` stores only the `Idepic_Bio`.

#### Checking overlap with your study

You have a list of your study's participants as EPIC `Idepic_Bio` IDs (one per line in
`my_participants.txt`). Match them against the manifest:

```r
library(data.table)

manifest <- fread("sample-manifest.tsv", na.strings = "")           # keep 'NA' text literal
mine     <- fread("my_participants.txt", header = FALSE)$V1          # your Idepic_Bio IDs

overlap <- manifest[Idepic_Bio %in% mine]
overlap[, .N, by = study][order(-N)]        # how many of your participants per genetics study
uniqueN(overlap$Idepic_Bio)                 # your participants that have genetics data
```

#### Building a `samples-keep.tsv`

To use the genetics data for your participants, take the overlap and drop cross-study
duplicates (keep each participant once, in the smallest study), giving the exact samples
to extract:

```r
# keep one genetics sample per participant: TRUE (kept copy) + NA (single-study),
# dropping FALSE (the duplicate copy held in a larger study)
keep <- overlap[overlap_keep != "FALSE"]

fwrite(keep[, .(study, IID)], "samples-keep.tsv", sep = "\t")
```

(If you instead want *every* genetics sample for your participants, replicates and all,
skip the `overlap_keep` filter and use `overlap` directly.)

#### Using `samples-keep.tsv` with the genetics data

The finalised genotypes are per study (`<STUDY>/<STUDY>_chr*.pgen|pvar|psam`). Subset
each study's PLINK2 files to your kept samples with `plink2 --keep`, matching on `IID`:

```bash
# samples-keep.tsv columns: study<TAB>IID
FINAL=/data/Epic/subprojects/Depot_Genetics/sources/finalised_data   # extracted <STUDY>/ pfiles
OUT=./mystudy_genetics; mkdir -p "$OUT"

for study in $(tail -n +2 samples-keep.tsv | cut -f1 | sort -u); do
  # plink2 --keep list: one IID per line under a '#IID' header
  { echo "#IID"; awk -F'\t' -v s="$study" 'NR>1 && $1==s {print $2}' samples-keep.tsv; } > "$OUT/$study.keep"

  for pgen in "$FINAL/$study/${study}"_chr*.pgen; do
    prefix="${pgen%.pgen}"
    plink2 --pfile "$prefix" \
           --keep "$OUT/$study.keep" \
           --make-pgen \
           --out "$OUT/$(basename "$prefix").mystudy"
  done
done
```

This writes per-study, per-chromosome PLINK2 files containing only your participants,
de-duplicated across studies. From there you can merge across studies/chromosomes as your
analysis requires (e.g. `plink2 --pmerge-list`).

## How To Run 

### 0: `.env` and `tools/`

You must first ensure that you have created a `.env` file at root, see: [.env.example](.env.example)

> **Important — set `SCRATCH_DATE` before running.**
> All pipeline stages write their outputs under `${SCRATCH}/${SCRATCH_DATE}/`. You must set `SCRATCH_DATE` to a date string (e.g. `2026-05-28`) in `.env` before submitting any stage, and keep it the same value for the entire analysis run. If you start a fresh analysis, update `SCRATCH_DATE` to a new date so the new run writes to a separate directory. The `.env.example` ships with `SCRATCH_DATE="CHANGE-ME"` as a deliberate placeholder.

You must set-up a minimal `conda` environment:

```bash
bash src/000_env.sh
```

You must download and compile all tools which are not available through `conda`; we run this as a job as it takes a while:

```bash
sbatch src/000_tools.sh
```

### 1: prepare study data

We create a copy of all required EPIC genetics data files:

```bash
sbatch src/001_data-genetics.sh
```

### 2: download reference data

We download all of the required reference data:
1. 1000 Genomes NYGC 2022 high-coverage VCFs (GRCh38)
2. SHAPEIT5 GRCh38 genetic map
3. dbSNP GRCh38 VCF
4. Annovar hg38 Database

```bash
bash src/002_data-reference.sh
```

### 3: prepare EPIC data

We need to create a reference file which provides information on sex and case status for each sample as this is not provided in the raw genetics data:

```bash
Rscript src/003_data-epic.R
```

### 4: stage1, stage2, and stage3

We can only run `stage1`, `stage2`, and `stage3` sequentially as `stage2`, and `stage3` are dependent upon the prior stages handoff data. We use the same `sbatch src/00*_stage*.sh` command for each.

We can run a `stage` for all studies simultaneously:

```bash
sbatch src/004_stage1.sh
```

Before progressing to `stage2` and `stage3` and from `stage3` to finalisation, you must look at the `stage1-summary.md`, `stage2-summary.md`, `stage3-summary.md` and `stage2` and `stage3` reports to check that the studies have completed and the pre-QC, phasing and imputation, and post-QC are good.

### 5: report and finalising

With all stages finished we generate the master cross-stage HTML reports and build the final deliverable archive for each study:

```bash
bash src/007_stage4.sh
```

### testing/other

We can perform a test across a single study if needed; we use `Glbd_01` for testing as it is the smallest study:

```bash
STAGE1_SCRIPTS=process_glbd_01.py sbatch src/004_stage1.sh
sbatch src/005_stage2.sh --study Glbd_01
sbatch src/006_stage3.sh --study Glbd_01
bash src/007_stage4.sh --study Glbd_01
```

Relatedness and ancestry outlier exclusions are off by default. To enable them:

```bash
sbatch src/006_stage3.sh --exclude-related --exclude-ancestry-outliers
```

## Methods And Thresholds Summary

### Stage 1

| Method area | What is done | Active arguments / thresholds | Notes |
| --- | --- | --- | --- |
| Study-specific preprocessing | Each study runs through one bespoke script to normalize file structure, SNP coding, and study-specific quirks before ID harmonization | No single global numeric threshold; logic is controlled by per-study settings for `BUILD`, `PRE`, `PART1`, `COMPLETION`, and `PART2_BUILD35_REL` | Methods include BIM sorting, AB-allele translation, name linkage, position reset, and optional early X/Y/XY/MT exclusion |
| Manifest comparison QC | `preProcessing.py` compares the raw study PLINK data to the array manifest and writes discrepancy files | No single global cut-off; exclusion sets are driven by manifest disagreement classes and archived study logic | Drives all downstream SNP exclusion and completion steps |
| Sample ID harmonization | Sample IDs are harmonized to EPIC IDs using the shared EPIC reference plus optional study-specific link files | `--update-sex`, `--remove`, `--update-ids`, `--make-pheno` with an empty phenotype file, `--indiv-sort n` | Stage 1 currently blanks phenotype intentionally after ID harmonization |
| SNP exclusion part 1 | Initial exclusion of manifest-disagreement SNPs | No single global numeric threshold; study-specific exclusion branches decide whether to remove negative-strand, chr/pos-mismatch, allele-mismatch, or unknown-strand SNPs | This is the first major study-specific branch point |
| Completion / repair step | Remaining variants have metadata repaired before the second exclusion pass | No single global numeric threshold; study-specific completion mode decides whether to update chromosome, position, strand, alleles, or rsID-based metadata | Converts the dataset into the best possible harmonized pre-final state |
| SNP exclusion part 2 | Final cleanup on the completed dataset | Always removes misplaced mitochondrial SNPs and duplicate classes 4, 2, 3, and 1; two studies also apply an extra build-35 exclusion list | Produces the final post-QC stage-1 prefix before liftover |
| Build harmonization | Final post-QC data are lifted to GRCh38 / hg38 | Build 36 -> `hg18ToHg38`; build 37 -> `hg19ToHg38`; final PLINK step uses `--split-x b38 no-fail` | Ensures all studies enter stage 2 on a common assembly |

### Stage 2

| Method area | What is done | Active arguments / thresholds | Notes |
| --- | --- | --- | --- |
| Pre-phasing QC | Stage-1 PLINK handoffs are filtered before phasing to remove problematic samples and variants | Sex check: F-stat < 0.2 (female), > 0.8 (male); het outliers: > 3.0 SD; KING duplicate cutoff: 0.354; HWE p < 0.000001; MAF < 0.005 | Requires ≥ 100 chrX variants for sex check; sex mismatches, het outliers, and duplicates are removed; HWE and MAF filters applied on autosomes |
| Target VCF export | QC-cleaned PLINK data are exported chromosome by chromosome to target VCF | Chromosomes renamed to UCSC style; IDs rewritten to `CHROM:POS:REF:ALT`; multiallelics split | Chromosomes absent in stage 1 are skipped cleanly |
| Reference preparation | 1000 Genomes reference VCFs are normalized and converted to per-chromosome BCF / `msav` files | Reference filter includes `MAC >= 10` during prep; chrX keeps SNPs and indels only | chrX is handled as `PAR1`, `nonPAR`, and `PAR2` |
| Phasing | SHAPEIT5 phases each study chromosome against the matching reference chromosome | `phase_cpus = 4`; uses SHAPEIT5 genetic recombination map | chrX is phased block by block as `PAR1`, `nonPAR`, `PAR2`; SHAPEIT5 outputs BCF which is converted to VCF.gz |
| Imputation | Minimac4 imputes each phased chromosome against the matching `msav` reference panel | `min_r2 = 0.3`; `minimac_batch_size = 200`; `minimac_threads = 4` | Final stage-2 VCFs retain imputation INFO metrics including `R2` |
| chrX handling | chrX is processed separately from autosomes | Blocks: `PAR1`, `nonPAR`, `PAR2`; block is skipped if no overlap or too few target variants exist; `chrx_min_ratio = 0.0` | Successful chrX blocks are concatenated; empty chrX output is written if none survive |
| Task runtime / retries | Internal stage-2 tasks are submitted to Slurm through Nextflow | Partition `low_p`; retries for exit codes `137`, `140`, `143`; `maxRetries = 2`; `errorStrategy = finish` | Tolerates interruption-style failures without masking real data errors; REPORTING retried up to 2 times for transient NFS errors |

### Stage 3

| Method area | What is done | Active arguments / thresholds | Notes |
| --- | --- | --- | --- |
| rsID annotation | Each imputed chromosome is annotated against dbSNP GRCh38 | Exact dbSNP match -> keep rsID; otherwise fallback to `chr:pos:REF:ALT` | A per-chromosome mapping file is written for ID auditability |
| Variant QC: imputation quality | Variants are filtered on imputation quality | `min_r2 = 0.3` | Applied before PLINK2 conversion |
| Variant QC: allele frequency | Variants are filtered on minor allele frequency | `maf = 0.01` | Applied together with the `R2` filter |
| Variant QC: HWE | Hardy-Weinberg equilibrium filtering is applied by default on autosomes | `run_hwe = true`; `hwe_p = 0.000005`; `hwe_k = 0`; mode `midp keep-fewhet` | chrX is not HWE filtered |
| PLINK2 conversion | Filtered imputed VCFs are converted to per-chromosome `pgen/pvar/psam` files; each chromosome is finalized independently | Uses dosage import from `HDS`; chrX uses `--split-par b38` and `--lax-chrx-import` | Final output is per-chromosome; no genome-wide merge is performed |
| LD pruning for QC | Per-chromosome LD pruning is run in parallel across all 22 autosomes; pruned BED files are merged into one small genome-wide BED for sample QC only | `--indep-pairwise 1500 150 0.2`; duplicate variant IDs resolved with `--set-all-var-ids '@:#:$r:$a' --rm-dup force-first` | The merged LD-pruned BED is an intermediate used only for KING and PCA; it is not the final output |
| Sample QC: relatedness | Related or duplicate samples are identified with KING | `king_cutoff = 0.0884` | Relatedness is always computed and reported; exclusion from the final dataset only happens when `--exclude-related` is set (default: off) |
| Sample QC: heterozygosity | Heterozygosity outliers are identified from the LD-pruned autosomal dataset | `het_sd_threshold = 3.0` | Samples beyond the SD threshold are added to the removal list |
| Sample QC: ancestry identification | PCA-based ancestry outliers are identified for every study | `ancestry_pc_count = 10`; `ancestry_z_threshold = 6.0` | Identification is always performed |
| Sample QC: ancestry exclusion | Ancestry outliers are optionally excluded from the final dataset | `exclude_ancestry_outliers = false` by default; use `--exclude-ancestry-outliers` to remove them | Ancestry outlier detection always runs; only the removal step is conditional |
| Sample QC: sex update | Recorded sex is updated from the stage-1 `.fam` during per-chromosome finalization | No separate numeric threshold | Applied in FINALIZE_CHROM alongside sample removal |

### Summary Table

**N** and **Variants** are Stage 3 final counts after all QC filters. **Mean ER2** is the variant-count-weighted mean empirical dosage R² across chromosomes (leave-one-out validation of genotyped variants only; higher values indicate that imputed dosages closely match the observed genotypes at genotyped sites). **Mean R2** is the mean theoretical imputation R² (Rsq) across all imputed variants from Stage 2 (higher values indicate high-confidence imputation). **AF Pearson R** is the Pearson correlation between study and 1000 Genomes reference allele frequencies across all imputed variants (values closer to 1 indicate good allele-frequency concordance with the reference panel).

| Study | N | Variants | Mean ER2 | Mean R2 | AF Pearson R |
| --- | ---: | ---: | ---: | ---: | ---: |
| Brea_01_Erneg | 863 | 11,600,099 | 0.9306 | 0.783 | 0.8636 |
| Brea_02 | 6,990 | 11,627,978 | 0.9335 | 0.688 | 0.8978 |
| Clrt_01 | 3,611 | 11,577,614 | 0.9267 | 0.763 | 0.8650 |
| Ecvd_01 | 7,104 | 11,601,382 | 0.9014 | 0.676 | 0.8915 |
| Ecvd_02 | 5,220 | 11,603,669 | 0.8116 | 0.502 | 0.9231 |
| Ecvd_03 | 5,274 | 10,260,727 | 0.9115 | 0.463 | 0.9132 |
| Glbd_01 | 94 | 11,757,294 | 0.9439 | 0.847 | 0.8462 |
| Inte_01 | 6,855 | 11,602,895 | 0.9327 | 0.734 | 0.8970 |
| Inte_02 | 5,547 | 11,678,033 | 0.8964 | 0.670 | 0.8932 |
| Inte_03 | 5,610 | 11,231,959 | 0.9080 | 0.685 | 0.8884 |
| Kidn_01 | 301 | 11,587,919 | 0.9421 | 0.818 | 0.8612 |
| Kidn_02 | 250 | 11,661,876 | 0.9771 | 0.909 | 0.8530 |
| Lung_01 | 2,132 | 11,642,277 | 0.9290 | 0.719 | 0.8784 |
| Lymp_01 | 437 | 11,663,657 | 0.9351 | 0.813 | 0.8639 |
| Neuro_01 | 4,831 | 11,853,445 | 0.8899 | 0.626 | 0.8912 |
| Ovar_01 | 1,174 | 11,676,124 | 0.9289 | 0.734 | 0.8748 |
| Panc_01 | 577 | 11,607,669 | 0.9378 | 0.800 | 0.8672 |
| Panc_02 | 137 | 11,622,662 | 0.9489 | 0.856 | 0.8495 |
| Pros_01 | 708 | 11,570,386 | 0.9397 | 0.795 | 0.8629 |
| Pros_03 | 1,049 | 11,649,351 | 0.9291 | 0.741 | 0.8704 |
| Pros_04 | 1,422 | 11,325,549 | 0.9325 | 0.758 | 0.8872 |
| Stom_01 | 240 | 11,652,085 | 0.9400 | 0.826 | 0.8608 |
| Uadt_01 | 198 | 11,659,622 | 0.9326 | 0.800 | 0.8618 |

*N and Variants reflect the Stage 3 final dataset after R²/MAF/HWE variant filtering and sex/relatedness/heterozygosity/ancestry sample QC. Imputation metrics are from Stage 2. Per-study details are in the master reports under `report/`.*

**Total unique participants across all 23 studies: 46,350** (60,624 total sample-study pairs; 9,578 participants appear in two or more studies). Counts are by EPIC `Idepic_Bio`, so a participant genotyped in several studies — including where the sample ID is stored in different forms across studies — is counted once.

### Sample Overlap

The UpSet plot below shows the intersection sizes across studies. Each bar represents the number of participants shared by the indicated combination of studies; horizontal bars on the left show each study's total sample count. Single-study bars (one dot) represent participants unique to that study.

![Sample overlap across EPIC genetics studies](docs/img/sample_overlap_upset.png)

## Filtering Steps Reference

All numeric filters applied across the three pipeline stages are listed below. Stage 1 and Stage 2 filters are set in the respective `nextflow.config` and `params.yaml` files; Stage 3 filters are set in `src/006_stage3.sh` and passed as Nextflow parameters.

| Filter | Stage | Threshold | Description |
| --- | --- | --- | --- |
| Variant MAF (pre-phasing) | Stage 1 | MAF ≥ 0.005 | Removes rare variants from the stage-1 PLINK handoff before phasing; applied to autosomal variants to improve phasing accuracy |
| Variant HWE (pre-phasing) | Stage 1 | p ≥ 0.000001 | Removes variants with extreme Hardy-Weinberg equilibrium departure (midp) on autosomes before phasing; reduces noise from genotyping errors |
| Sex check (F-stat) | Stage 1 | F-stat < 0.2 → female; F-stat > 0.8 → male | Checks recorded sex against X chromosome F-statistic; samples with mismatching inferred sex are excluded; requires ≥ 100 chrX variants |
| Heterozygosity outliers | Stage 1 | > 3.0 SD from study mean | Removes samples with excess or deficit autosomal heterozygosity, indicating sample contamination or extreme inbreeding |
| Duplicate / MZ twin removal | Stage 1 | KING kinship > 0.354 | Removes exact duplicates and monozygotic twin pairs; the sample with more missing genotype data is removed |
| Reference panel MAC | Stage 2 | MAC ≥ 10 | Applied during 1000G reference panel preparation; removes very low-frequency variants from the reference to ensure reliably phased haplotypes |
| chrX nonPAR phasing ratio | Stage 2 | ratio ≥ 0.0 | Minimum fraction of target variants overlapping the reference in the chrX nonPAR block required for phasing to proceed; set to 0.0 so the block is always attempted when any overlap exists |
| Imputation quality (R²) | Stage 2 | R² ≥ 0.3 | Minimac4 `--min-r2` flag; imputed variants below this R² threshold are excluded from the stage-2 output VCFs |
| Empirical validation minimum N | Stage 2 | N ≥ 20 samples | Minimum study size to compute empirical dosage R² (leave-one-out validation at genotyped sites); smaller studies skip this metric |
| Dose-zero minimum N | Stage 2 | N ≥ 5 samples | Minimum study size to report dose-zero statistics in the stage-2 report |
| Summary report R² threshold | Stage 2 | R² ≥ 0.3 | Minimum R² for variants included in stage-2 summary counts and plots |
| High-quality imputation threshold | Stage 2 | R² ≥ 0.8 | R² threshold used to classify a variant as high quality in the stage-2 summary report; reported as a separate count alongside the total imputed variants |
| Imputation quality (R²) | Stage 3 | R² ≥ 0.3 | Post-imputation filter on `INFO/R2` applied before PLINK2 conversion; removes poorly imputed variants from the final dataset |
| Minor allele frequency | Stage 3 | MAF ≥ 0.01 | Applied together with the R² filter before PLINK2 conversion; removes very rare variants from the final dataset |
| Hardy-Weinberg equilibrium | Stage 3 | p ≥ 0.000005 | Applied on autosomes only (`midp keep-fewhet`); chrX is not HWE filtered; removes variants deviating from equilibrium expectation |
| Relatedness / kinship | Stage 3 | KING kinship ≥ 0.0884 | Identifies sample pairs at 2nd-degree or closer relationship; the identified count is always reported; one sample per pair is removed only when `--exclude-related` is set (default: off) |
| Heterozygosity outliers | Stage 3 | > 3.0 SD from study mean | Removes samples with excess or deficit heterozygosity from the LD-pruned autosomal dataset post-imputation |
| Ancestry outliers | Stage 3 | z-score ≥ 6.0 on 10 PCs | PCA-based ancestry outlier identification using 1000G reference principal components; the identified count is always reported; outliers removed only when `--exclude-ancestry-outliers` is set (default: off) |

## Finalised Data

Each study's deliverable is written to `final/<STUDY>/` by stage 4. This section describes what is in that directory and how to work with it.

### Directory layout

```
final/<STUDY>/
├── <STUDY>.tar.gz        # per-chromosome PLINK2 files (pgen/pvar/psam)
├── report-master.html    # cross-stage master QC report
├── report-stage2.html    # stage 2 imputation report
├── report-stage3.html    # stage 3 sample and variant QC report
└── review/
    ├── related.exclude   # sample IDs flagged as related (KING ≥ 0.0884)
    ├── ancestry.exclude  # sample IDs flagged as ancestry outliers (PCA z-score ≥ 6.0)
    ├── hwe.exclude       # variant IDs failing HWE (p < 0.000005) on autosomes
    ├── het.het           # per-sample autosomal heterozygosity statistics
    ├── king.kin0         # pairwise kinship coefficients (KING format)
    ├── pca.eigenvec      # per-sample PCA scores (10 PCs)
    └── pca.eigenval      # PCA eigenvalues
```

### Extracting the archive

Extract the tarball in the `final/<STUDY>/` directory. This produces a `<STUDY>/` subdirectory containing all per-chromosome PLINK2 files:

```bash
cd final/<STUDY>
tar -xzf <STUDY>.tar.gz
```

The extracted directory contains 23 chromosome pairs plus one sample file:

```
<STUDY>/
├── <STUDY>_chr1.pgen
├── <STUDY>_chr1.pvar
...
├── <STUDY>_chr22.pgen
├── <STUDY>_chr22.pvar
├── <STUDY>_chrX.pgen
├── <STUDY>_chrX.pvar
└── <STUDY>.psam
```

Each chromosome is a self-contained PLINK2 dataset. The `.psam` is shared across all chromosomes; all per-chromosome `pgen/pvar` files reference the same sample order.

### Loading the data in PLINK2

Reference a chromosome using `--pfile` with the chromosome-specific prefix (no extension):

```bash
plink2 --pfile <STUDY>/<STUDY>_chr1 [...]
```

All chromosomes share the same sample list, so the `.psam` from any chromosome can be used as the sample reference. PLINK2 will locate the matching `.psam` automatically from the prefix.

### Applying QC exclusions

The pipeline applies sample exclusions for heterozygosity outliers by default. Relatedness and ancestry exclusions are **off by default** and are provided as `review/` files for downstream analysts to apply as appropriate for each analysis.

The HWE exclude list is also intentionally **not hard-applied** to the final PGEN files: HWE filtering may be study-specific (e.g. case-only strata or X-linked variants), so the list is provided for analyst use rather than applied universally.

#### Sample exclusions

`review/related.exclude` and `review/ancestry.exclude` are sample-level files used with `--remove`:

```bash
# Apply relatedness exclusion to a single chromosome
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --remove review/related.exclude \
  --make-pgen \
  --out <STUDY>_chr1_unrelated
```

```bash
# Apply both relatedness and ancestry exclusions together
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --remove <(cat review/related.exclude review/ancestry.exclude | sort -u) \
  --make-pgen \
  --out <STUDY>_chr1_filtered
```

#### Variant exclusions

`review/hwe.exclude` is a variant-level file used with `--exclude`:

```bash
# Apply HWE exclusion to a single chromosome
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --exclude review/hwe.exclude \
  --make-pgen \
  --out <STUDY>_chr1_hwe
```

#### Combined sample and variant exclusions

```bash
plink2 \
  --pfile <STUDY>/<STUDY>_chr1 \
  --remove review/related.exclude \
  --exclude review/hwe.exclude \
  --make-pgen \
  --out <STUDY>_chr1_qc
```

To apply exclusions across all chromosomes in a loop:

```bash
STUDY="Glbd_01"
for chr in {1..22} X; do
  plink2 \
    --pfile ${STUDY}/${STUDY}_chr${chr} \
    --remove review/related.exclude \
    --exclude review/hwe.exclude \
    --make-pgen \
    --out ${STUDY}_qc/${STUDY}_chr${chr}
done
```

### Review files

The `review/` files are written by stage 3 and preserved in the deliverable for post-hoc analyst inspection:

| File | Contents | Use |
| --- | --- | --- |
| `related.exclude` | Sample IDs with KING kinship ≥ 0.0884 (2nd-degree or closer) | Pass to `--remove` to drop related individuals |
| `ancestry.exclude` | Sample IDs identified as ancestry outliers on 10 PCs (z-score ≥ 6.0) | Pass to `--remove` to restrict to the primary ancestry cluster |
| `hwe.exclude` | Variant IDs failing HWE (p < 0.000005) on autosomes | Pass to `--exclude` to remove HWE-failing variants |
| `het.het` | Per-sample autosomal heterozygosity (PLINK2 `.het` format) | Inspect outlier distribution; samples beyond 3 SD are already removed from the PGEN files |
| `king.kin0` | Pairwise kinship coefficients for all sample pairs (KING `.kin0` format) | Inspect relatedness graph; pairs above 0.0884 appear in `related.exclude` |
| `pca.eigenvec` | Per-sample scores on 10 PCs | Use for stratification correction in association analyses; ancestry outliers flagged in `ancestry.exclude` are visible as outliers in these scores |
| `pca.eigenval` | Variance explained by each of the 10 PCs | Useful for scree plots and deciding how many PCs to include as covariates |
