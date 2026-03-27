# Major-and-minor-changes-in-R-and-Bioconductor-since-2017

# Table of Contents
* [# R & Bioconductor: What Changed 2014–2017](https://github.com/swajid/Major-and-minor-changes-in-R-and-Bioconductor-since-2017/edit/main/README.md#r--bioconductor-what-changed-20142017)
* [# R & Bioconductor: What Changed Since 2017](https://github.com/swajid/Major-and-minor-changes-in-R-and-Bioconductor-since-2017/edit/main/README.md#r--bioconductor-what-changed-20142017) 
* [# R UI Development: Changes from 2012 to 2026](https://github.com/swajid/Major-and-minor-changes-in-R-and-Bioconductor-since-2017/edit/main/README.md#r-ui-development-2012-to-2026)

# R & Bioconductor: What Changed 2014–2017
### The Foundation Layer — Tidyverse, SummarizedExperiment, Rcpp Attributes

> **Scope:** R 3.1 → 3.4 · Bioconductor 3.0 → 3.5 · Rcpp 0.11 → 0.12 · Package Development · OOP Class Systems  
> **Context:** This era crystallized the idioms that the 2017–2026 ecosystem was built on top of.

---

## Table of Contents

1. [R Language (3.1 → 3.4)](#1-r-language-31--34)
2. [The Tidyverse Emerges](#2-the-tidyverse-emerges)
3. [Bioconductor (3.0 → 3.5)](#3-bioconductor-30--35)
4. [Rcpp (0.11 → 0.12)](#4-rcpp-011--012)
5. [Package Development](#5-package-development)
6. [OOP Class Systems](#6-oop-class-systems)
7. [Summary — What This Era Set Up](#7-summary--what-this-era-set-up)

---

## 1. R Language (3.1 → 3.4)

This period was relatively quiet at the language level — no breaking changes on the scale of R 4.0 — but a few things mattered day-to-day.

| Version | Year | Key Change |
|---------|------|-----------|
| R 3.1 | 2014 | Improved `c()` attribute preservation; fixed subtle bugs in named vector pipelines |
| R 3.2 | 2015 | `vapply()` stabilization; better `tryCatch()` behaviour |
| R 3.3 | 2016 | `lengths()` added — vectorized `length()` over lists |
| R 3.4 | 2017 | JIT bytecode compiler enabled by default; ~10–20% speedup to pure-R loops, no code changes required |

### Notable: `lengths()` for Genomics

`lengths()` sounds trivial but was genuinely useful for genomic interval list processing:

```r
# How many variants per chromosome?
chr_counts <- lengths(variants_by_chr)   # vectorized, no sapply needed
```

### Notable: The Pipe Goes Mainstream

The most significant language-adjacent event of this period was **magrittr's `%>%` pipe going mainstream** via `dplyr` adoption around 2014. By 2017 it had fundamentally changed how R code was written — this is the direct ancestor of the native `|>` added in R 4.1.

```r
# Before (~2013): nested calls
filter(select(read.csv("samples.csv"), tumor_id, purity), purity > 0.3)

# After (~2015+): pipe style
read.csv("samples.csv") %>%
  select(tumor_id, purity) %>%
  filter(purity > 0.3)
```

---

## 2. The Tidyverse Emerges

The word **"tidyverse"** was coined in 2016 by Hadley Wickham. Before that, `ggplot2`, `dplyr`, `tidyr`, and `readr` were separate packages with no unifying brand. By 2017 the `tidyverse` meta-package was on CRAN.

### Key Package Releases

| Package | What happened | Why it mattered |
|---------|--------------|-----------------|
| `dplyr` 0.4 → 0.7 | Window functions; `group_by` stabilization; `tidyeval` replacing NSE | Standard for sample metadata and variant table manipulation |
| `tidyr` | Split from `reshape2`; `gather()`/`spread()` API | Replaced `reshape2::melt()` for long/wide transformations |
| `purrr` 0.1 (2015) | Replaced `plyr`'s list operations with functional API | `map()`, `map_dfr()` replaced `lapply()` + `do.call(rbind, ...)` |
| `readr` | Replaced `read.csv`; 10x faster; no `stringsAsFactors` | Standard for loading clinical metadata and variant tables |
| `ggplot2` 2.0 | Major API stabilization | Stable enough to build on for publication figures |
| `tibble` | `data.frame` replacement with better printing | No row names by default; never silently converts strings |

### The tidyeval / NSE Shift

`dplyr` 0.7 introduced **tidyeval** (`rlang`) to replace the older non-standard evaluation (NSE) approach. This affected anyone writing functions that wrapped `dplyr` verbs:

```r
# Old NSE (dplyr < 0.7): fragile inside functions
filter_by_col <- function(df, col) {
  filter(df, col > 0)   # broken — col is evaluated in wrong environment
}

# New tidyeval (dplyr 0.7+): explicit quoting
filter_by_col <- function(df, col) {
  col <- enquo(col)
  filter(df, !!col > 0)
}
```

This was a source of significant pain for anyone who had written helper functions wrapping dplyr — essentially all such code needed rewriting.

---

## 3. Bioconductor (3.0 → 3.5)

Package count grew from ~900 to ~1,300. The infrastructure push was toward **`SummarizedExperiment` replacing `ExpressionSet`** — this transition was underway but not complete by 2017.

### 3.1 Core Infrastructure — The SE Transition Begins

```r
# ExpressionSet (old, still dominant in 2014)
eset <- ExpressionSet(
  assayData   = exprs_matrix,
  phenoData   = AnnotatedDataFrame(sample_info),
  featureData = AnnotatedDataFrame(gene_info)
)
exprs(eset)          # get expression matrix
pData(eset)          # get sample metadata

# SummarizedExperiment (new, increasingly standard by 2017)
se <- SummarizedExperiment(
  assays    = list(counts = count_matrix),
  colData   = sample_info,
  rowRanges = gene_ranges   # GRanges — genomic coordinates on rows
)
assay(se, "counts")  # get count matrix
colData(se)          # get sample metadata
```

By 2017 the community was split — Bioconductor core was pushing SE hard, but most existing packages and tutorials still used `ExpressionSet`. The full transition was completed in the 2018–2020 era.

### 3.2 GenomicRanges — Mature Interval Operations

`GRanges` and `GenomicRanges` matured significantly. By 2016 these were the standard for any BED/BAM region operation:

```r
library(GenomicRanges)

# Operations that replaced bedtools calls in R
overlaps   <- findOverlaps(tumor_regions, gene_ranges)
coverage   <- coverage(aligned_reads)
gaps       <- gaps(exon_ranges)
reduced    <- reduce(fragmented_ranges)    # merge overlapping intervals
flanked    <- flank(promoter_ranges, 500)  # extend intervals
```

`GenomicAlignments` replaced much of what people were doing with `samtools` calls piped into R — reading BAM files and computing per-base coverage became native R operations.

### 3.3 Variant Analysis Infrastructure

```r
library(VariantAnnotation)

# VCF reading stabilized as the canonical interface
vcf <- readVcf("tumor.vcf", "hg38")
geno(vcf)$GT         # genotype matrix
info(vcf)$AF         # allele frequency from INFO field
fixed(vcf)           # CHROM, POS, REF, ALT, QUAL, FILTER
```

`VariantAnnotation` stabilized as the canonical VCF parser in R, replacing ad-hoc `read.table()` approaches on VCF files.

### 3.4 Differential Expression — DESeq2 Takes Over

`DESeq2` (released late 2014) rapidly displaced `DESeq` and became the most-cited bioinformatics software package in existence within a few years:

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "tumor", "normal"))
# Note: by 2017, lfcShrink() was not yet standard — that came later
```

`limma-voom` was refined alongside this for RNA-seq, remaining the preferred method for datasets where linear modelling and array-era infrastructure mattered.

### 3.5 Copy Number — FACETS Arrives

**FACETS** (released 2015–2016) was developed at MSKCC and became the standard for allele-specific copy number estimation from tumor-normal sequencing pairs, directly feeding into IMPACT pipeline outputs. `DNAcopy` (CBS segmentation) underpins FACETS and was already a mature Bioconductor package by this point.

```r
library(facets)

# The core FACETS workflow established in this era
rcmat  <- readSnpMatrix("tumor_normal_pileup.csv")
xx     <- preProcSample(rcmat)
oo     <- procSample(xx, cval = 150)
fit    <- emcncf(oo)
# Output: dipLogR, purity, ploidy, per-segment tcn/lcn
```

### 3.6 Annotation Resources

`BSgenome` and `TxDb` annotation packages became the reference standard for genome and transcript coordinates, replacing manual GTF parsing:

```r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)

txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
exons <- exons(txdb)                        # all exons as GRanges
genes <- genes(txdb)                        # all genes as GRanges
seq   <- getSeq(BSgenome.Hsapiens.UCSC.hg38, roi_granges)  # sequence retrieval
```

`AnnotationHub` was introduced in this period — the idea of a curated, queryable hub of annotation resources was new and not yet universally adopted by 2017.

---

## 4. Rcpp (0.11 → 0.12)

Rcpp crossed **1,000 reverse dependencies on CRAN** in this period and became the de facto standard for compiled R extensions.

### Rcpp Attributes Becomes Universal

The **`// [[Rcpp::export]]`** attribute system became universal, replacing the older `.Call()` boilerplate entirely:

```cpp
// Old way (pre-Attributes): manual .Call() registration
// Required writing SEXP wrappers, registering with R_registerRoutines, etc.

// New way (Attributes): just annotate and sourceCpp() / devtools::document()
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector segment_log2r(NumericVector log2r, int min_width = 5) {
  // ... CBS-style segmentation logic
  return segments;
}
```

`sourceCpp()` for interactive prototyping became standard in this period — write C++ inline in an R session, compile and test immediately.

### RcppArmadillo and RcppEigen Mature

For matrix-heavy genomics work (PCA, NMF, regression on large variant matrices):

```cpp
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat fast_pca(arma::mat X, int n_comp) {
  arma::mat U, V;
  arma::vec s;
  arma::svd_econ(U, s, V, X);
  return U.cols(0, n_comp - 1);
}
```

### Modern C++ Plugins

Experimental support for C++11/14 via plugin syntax:

```cpp
// [[Rcpp::plugins(cpp14)]]
#include <Rcpp.h>

// [[Rcpp::export]]
auto sum_squares(std::vector<double> x) {   // auto return type (C++14)
  return std::accumulate(x.begin(), x.end(), 0.0,
    [](double acc, double v) { return acc + v * v; });
}
```

---

## 5. Package Development

`devtools` was the monolith — one package to rule them all. The split into `usethis`, `pkgbuild`, `pkgload`, etc. had not yet happened.

### The 2014–2017 Workflow

```r
# The full development cycle in this era
library(devtools)

create("mypackage")           # scaffold a new package
load_all()                    # load package for interactive testing
document()                    # run roxygen2 → generate .Rd files
test()                        # run testthat suite
check()                       # R CMD check
install_github("user/pkg")    # install from GitHub
```

### roxygen2 — No Markdown Yet

`roxygen2` was in wide use but **without Markdown support** (that arrived in 2017–2018). Documentation was written in Rd-flavoured syntax directly in comments:

```r
#' Compute GC content
#'
#' Calculates the GC fraction of a DNA sequence string.
#'
#' @param seq A \code{character} vector of DNA sequences.
#' @return A \code{numeric} vector of GC fractions between 0 and 1.
#' @examples
#' gc_content("ATGCGC")
#' @export
gc_content <- function(seq) { ... }
```

Note the `\code{}` Rd macros — in 2017+, these become backtick Markdown.

### testthat Edition 1/2

```r
# testthat in this era
test_that("gc_content returns correct value", {
  expect_equal(gc_content("ATGC"), 0.5)
  expect_error(gc_content(123))
  expect_that(gc_content("GGG"), equals(1.0))  # expect_that() — deprecated by end of era
})
```

`expect_that()` was the original generic expectation function and was deprecated by 2017 in favour of the specific `expect_*` functions. `expect_snapshot()` did not exist yet.

### CI: Travis CI

GitHub Actions did not exist (launched 2018). The standard was `.travis.yml`:

```yaml
# .travis.yml — the 2014-2017 CI standard for R packages
language: r
r:
  - release
  - devel
r_packages:
  - covr
after_success:
  - Rscript -e 'covr::codecov()'
```

### Dependency Management: packrat

`renv` did not exist. `packrat` was the dependency isolation solution — it worked by snapshotting package versions into a `packrat/` directory inside the project. It was functional but slow and not widely adopted. Most teams relied on documented `sessionInfo()` output for reproducibility instead.

---

## 6. OOP Class Systems

### The Landscape in 2014–2017

| System | Status | Used by |
|--------|--------|---------|
| **S3** | Dominant for CRAN | Base R, most user-facing packages |
| **S4** | Dominant for Bioconductor | `GRanges`, `ExpressionSet`, `SummarizedExperiment` |
| **R5 / RC** | Effectively deprecated in practice | Almost nothing new |
| **R6** | Newly released (2014), gaining traction | `shiny`, `httpuv`, infrastructure packages |
| **S7** | Does not exist | — |

### S4 in Bioconductor — The Established Pattern

S4 was already the mandated class system for Bioconductor core infrastructure. The `ExpressionSet`/`SummarizedExperiment` transition in this period demonstrated S4's strength — formal slot validation and inheritance allowed `RangedSummarizedExperiment` to cleanly extend `SummarizedExperiment`:

```r
# S4 in action — the Bioconductor idiom
setClass("CopyNumberResult",
  contains  = "SummarizedExperiment",
  slots     = list(
    purity  = "numeric",
    ploidy  = "numeric",
    dipLogR = "numeric"
  )
)

setGeneric("purity", function(x) standardGeneric("purity"))
setMethod("purity", "CopyNumberResult", function(x) x@purity)

setValidity("CopyNumberResult", function(object) {
  if (object@purity < 0 || object@purity > 1) "purity must be in [0,1]"
  else TRUE
})
```

### R6 Arrives (2014)

`R6` was released in 2014 and began replacing Reference Classes in infrastructure packages. `shiny` and `httpuv` switched from RC to R6, yielding substantial performance improvements. By 2017 it was not yet common in data analysis packages — that shift happened post-2018.

```r
library(R6)

# R6 in 2014-2017: primarily used in infrastructure, not analysis code
Counter <- R6Class("Counter",
  public = list(
    count = 0,
    increment = function() {
      self$count <- self$count + 1
      invisible(self)
    }
  )
)
```

---

## 7. Summary — What This Era Set Up

This 2014–2017 window is the **foundation layer** for everything that followed. Every major shift in the 2017–2026 ecosystem has a direct root here:

| 2014–2017 development | What it became by 2026 |
|-----------------------|------------------------|
| `%>%` pipe via magrittr | Native `\|>` pipe in R 4.1 |
| tidyverse brand coined | Universal R idiom; `tidybulk`, `plyranges`, tidy bioinformatics |
| `SummarizedExperiment` transition begins | Universal Bioconductor container; `SingleCellExperiment` extends it |
| `GRanges` matures | Foundation for all interval-based genomics in R |
| `DESeq2` displaces `DESeq` | Still the DE standard in 2026; `lfcShrink` added later |
| FACETS released (2015–2016) | Embedded in MSKCC IMPACT pipeline; basis for clinical CNV reporting |
| Rcpp Attributes universal | Direct lineage to cpp11, and contrast point for extendr/Rust |
| `devtools` as monolith | Split into `usethis`, `pkgbuild`, `pkgload`, etc. by 2019 |
| Travis CI | Replaced entirely by GitHub Actions from 2019 onward |
| `packrat` for reproducibility | Replaced by `renv` (2019) |
| R6 arrives | Becomes dominant encapsulated OOP; S7 later attempts to unify |

---

*Previous era (2017–2026): see `R_Bioconductor_Changes_Since_2017.md`*  
*Primary references: R release NEWS files · Bioconductor release notes · Wickham, "Advanced R" (1e 2014, 2e 2019) · Eddelbuettel, "Seamless R and C++ Integration with Rcpp" (2013)*



# R & Bioconductor: What Changed Since 2017
### A Reference for Computational Biologists Returning to the Ecosystem

> **Scope:** R 3.4 → 4.5+ · Bioconductor 3.5 → 3.22 · Rcpp · Rust/extendr · Package Development · OOP Class Systems  
> **Last updated:** March 2026

---

## Table of Contents

1. [R Language (Base)](#1-r-language-base)
2. [Bioconductor](#2-bioconductor)
3. [Rcpp — C++ Integration](#3-rcpp--c-integration)
4. [Rust in R — extendr / rextendr](#4-rust-in-r--extendr--rextendr)
5. [Package Development Toolchain](#5-package-development-toolchain)
6. [OOP Class Systems](#6-oop-class-systems)
7. [Quick Decision Guide](#7-quick-decision-guide)

---

## 1. R Language (Base)

### R 4.0.0 (2020) — The Big Break

This was the most significant release since R 3.0.0. Several changes silently broke old code:

| Change | Old behaviour | New behaviour |
|--------|--------------|---------------|
| `stringsAsFactors` | `TRUE` by default | `FALSE` by default |
| Raw string literals | Not available | `r"(...)"` syntax supported |
| Reference counting | NAMED mechanism | Proper reference counting; fewer unnecessary copies |
| `matrix()` | Strings stayed strings | Character columns → factors; factors → integers |

**`stringsAsFactors = FALSE` is the most practically impactful change.** Years of defensive `stringsAsFactors = FALSE` boilerplate in data loading code can now be removed. But any old code that *relied* on automatic factor conversion (e.g. genomics pipelines that expected sample metadata columns to be factors) may silently behave differently.

### R 4.1 — Pipe and Lambda

```r
# Native pipe (no magrittr required)
mtcars |> head()

# Native lambda / anonymous function shorthand
\(x) x + 1        # equivalent to function(x) x + 1
```

These feel minor but have large downstream effects — the native pipe `|>` does not support placeholder `.` by default (unlike magrittr's `%>%`), which means some pipeline patterns need rewriting. R 4.2 added `_` as a placeholder on the right-hand side.

### R 4.2 — Placeholder in Native Pipe

```r
mtcars |> lm(mpg ~ cyl, data = _)   # _ is the placeholder
```

### R 4.5 (2025)

- `grepv()` added — returns matching text (not indices), complementing `grep()`
- C++ standard bumped to C++20 by default where available (affects compiled packages)
- `grepv()` example relevant in bioinformatics: parsing variant annotation text fields

### R 4.5.3 / 4.6.0 (2026)

- R 4.5.3 ("Reassured Reassurer") released March 2026
- R 4.6.0 ("Happy Hop") scheduled April 2026
- New BLAS routines (`dgemmtr`, `zgemmtr`) introduced — **BLAS implementations are no longer freely interchangeable**. If you swap BLAS on an HPC cluster (e.g. swapping to OpenBLAS for FACETS speed gains), check that your BLAS version is ≥ OpenBLAS 0.3.29

### Other Language-Level Changes Worth Knowing

- **`|>` pipe with `_` placeholder** (4.2+): `df |> subset(select = _)` works but has restrictions vs magrittr
- **Better error messages**: R 4.3+ overhauled error and warning formatting significantly
- **`switch()` on numeric**: now warns on non-integer input — catches latent bugs in condition dispatching
- **`readline()` and UTF-8**: better Unicode handling throughout; relevant for sample name parsing

---

## 2. Bioconductor

**Package count:** ~1,300 (2017) → **2,289** (2024)  
**PubMedCentral citations:** >60,500 as of 2024

### 2.1 Core Data Infrastructure

The single biggest architectural shift: **everything is now `SummarizedExperiment`-derived.**

```r
# 2017 style: list of matrices
exprs_mat <- assays$exprs
coldata <- phenoData@data

# 2024 style: SummarizedExperiment
se <- readRDS("tumor_se.rds")
counts_mat <- assay(se, "counts")
col_metadata <- colData(se)
row_metadata <- rowData(se)
```

Key container classes and their roles:

| Class | Package | Purpose |
|-------|---------|---------|
| `SummarizedExperiment` | `SummarizedExperiment` | Rectangular assay data + row/col metadata |
| `RangedSummarizedExperiment` | same | SE + genomic ranges on rows |
| `SingleCellExperiment` | `SingleCellExperiment` | scRNA-seq, extends SE |
| `MultiAssayExperiment` | `MultiAssayExperiment` | Multi-omic patient-level data |
| `SpatialExperiment` | `SpatialExperiment` | Spatial transcriptomics |

If you wrote code in 2017 that operated on `ExpressionSet` objects (the old Bioconductor standard), those are now effectively legacy. Most new packages do not accept `ExpressionSet`.

### 2.2 Single-Cell Analysis (The Biggest Growth Area)

In 2017 at BioC2017, a `SingleCellExperiment` class was *just being proposed*. It is now the foundation of the entire single-cell Bioconductor sub-ecosystem.

The canonical workflow reference is the **OSCA book** (Orchestrating Single-Cell Analysis with Bioconductor): https://osca.bioconductor.org

Key packages:

```r
# QC and normalization
library(scater)    # QC metrics, dimension reduction, visualization
library(scran)     # normalization, HVG selection, graph-based clustering

# Cell type annotation
library(SingleR)   # automated annotation via reference datasets
library(celldex)   # curated reference datasets (ImmGen, ENCODE, etc.)

# Clustering
library(bluster)   # unified interface to clustering algorithms

# Trajectory / pseudotime
library(slingshot) # lineage inference on low-dimensional embeddings
library(TSCAN)     # tree-based pseudotime

# Differential expression (single-cell aware)
library(muscat)    # multi-sample, multi-group pseudobulk DE
```

**For cancer genomics specifically:** packages like `copykit`, `infercnv`, and `CopyNumberPlots` enable CNV inference from scRNA-seq — directly relevant if you're moving from bulk IMPACT-style pipelines toward single-cell tumor heterogeneity analysis.

### 2.3 Spatial Transcriptomics

Entirely absent in 2017. Now a full sub-ecosystem:

```r
library(SpatialExperiment)  # data container extending SCE
library(nnSVG)              # spatially variable gene detection
library(CARD)               # cell type deconvolution
library(BayesSpace)         # spatial clustering
```

### 2.4 Multi-Modal / Multi-Assay

```r
library(MultiAssayExperiment)  # container for multi-omic patient data
# Integrates WES + RNA-seq + methylation + clinical data per patient
# Directly relevant for TCGA-style and GENIE cohort analyses

mae <- MultiAssayExperiment(
  experiments = list(rnaseq = rse, cnv = cnv_se, mut = mut_se),
  colData = patient_clinical_df
)
```

`MultiAssayExperiment` has become the standard container for the kind of multi-modal data you work with in clinical cancer genomics. If you're doing integrative TCGA analyses, this is essential.

### 2.5 Differential Expression — Still DESeq2/limma, But Updated

DESeq2 and limma/voom remain the workhorses, but:

- **`DESeq2`**: now recommends `lfcShrink()` (not the old `results()`) for fold change estimates; adaptive shrinkage via `apeglm` is the default
- **`limma`**: `voomWithQualityWeights()` more widely used; better handling of low-count genes
- **`edgeR`**: `glmQLFit()` / `glmQLFTest()` now preferred over the older `glmFit()` / `glmLRT()` for small-sample datasets

### 2.6 Cloud and HPC Integration

- **AnVIL** (Terra/GCP): Bioconductor packages run natively on Terra via `AnVIL` package
- **BiocParallel**: mature parallel backend, integrates with LSF (your HPC environment), SLURM, and cloud schedulers
- **HDF5Array / DelayedArray**: on-disk array representation for datasets too large for RAM — critical for large single-cell data

```r
library(HDF5Array)
# Load a 500k-cell count matrix without loading into memory
counts <- HDF5Array("counts.h5", "matrix")
```

### 2.7 Annotation Resources — AnnotationHub and ExperimentHub

```r
library(AnnotationHub)    # access to genome annotations, GTFs, chain files
library(ExperimentHub)    # access to curated experiment datasets

ah <- AnnotationHub()
# Query for hg38 Ensembl annotation
query(ah, c("Homo sapiens", "EnsDb", "GRCh38"))
```

These replaced the pattern of manually downloading and loading annotation files. If you were doing `read.table("Homo_sapiens.GRCh38.gtf")` — this is the modern replacement.

---

## 3. Rcpp — C++ Integration

Rcpp has matured considerably since 2017. It is now at version **1.0.13+** and used by **~14,000 CRAN packages**.

### Key Changes Since 2017

**Safety improvements:**
- Rcpp now **emits a warning on out-of-bounds vector access** (will become an error in a future release). Old code that accidentally read past vector ends will now surface these bugs
- Switched from non-API `DATAPTR` to compliant API variants (`VECTOR_PTR`, `STRING_PTR`) — this matters if you're building against newer R headers
- `count` variables switched to `size_t` to avoid conversion-narrowing warnings that appear with modern compilers

**API modernisation:**
- `getRcppVersion()` function available to programmatically check API version at runtime
- C++ standard is now C++20 by default in R 4.6+ — Rcpp packages benefit from modern C++ features (concepts, ranges, `std::span`, etc.)
- The `Rcpp-modules` vignette received a major review — using `RCPP_MODULE` to expose C++ classes to R is now better documented

**Workflow:**
```cpp
// Modern Rcpp attribute style (mature since ~2013 but now universal)
// [[Rcpp::export]]
NumericVector fast_log_sum(NumericVector x) {
  return log(sum(exp(x)));  // Rcpp sugar — vectorized, readable
}
```

**Rcpp Sugar** (bringing R-like syntax into C++) is now widely used and well-tested. If you wrote raw C loops in 2017, there's almost always a sugar equivalent.

### cpp11 — The Lightweight Alternative

Introduced by Posit (~2020), `cpp11` is a header-only alternative to Rcpp that:
- Has a smaller footprint (faster compile times)
- Uses modern C++11 idioms throughout
- Is now used in `dplyr`, `vroom`, `readr`, and other tidyverse internals

```cpp
#include "cpp11.hpp"
[[cpp11::register]]
doubles log_transform(doubles x) {
  writable::doubles out(x.size());
  for (int i = 0; i < x.size(); i++) out[i] = std::log(x[i]);
  return out;
}
```

For new packages, **cpp11 is worth considering** if you want faster compilation and tighter C++11 semantics. For genomics packages with heavy linear algebra, Rcpp + RcppArmadillo/RcppEigen remains the standard.

---

## 4. Rust in R — extendr / rextendr

This is entirely new since 2017. **Rust did not have an R interface in 2017.**

### Why Rust?

Rust offers memory safety without a garbage collector — meaning you get C-level performance with compile-time guarantees against null pointer dereferences, buffer overflows, and data races. For genomics workloads (parsing large BAM/VCF files, alignment operations, string-heavy bioinformatics), Rust's `rayon` crate makes parallelism trivial and safe.

### The extendr Ecosystem

- **`extendr`** (Rust crate): the core interface between R's C API and Rust
- **`rextendr`** (R package, on CRAN since 2021): developer tooling, analogous to how `Rcpp` tools work

```r
# Setup a Rust-powered R package
usethis::create_package("mypackage")
rextendr::use_extendr()          # scaffolds src/rust/ structure
rextendr::document()             # compiles Rust, generates R wrappers
rextendr::rust_sitrep()          # diagnostic: checks Rust toolchain
```

### What Rust Code Looks Like in an R Package

```rust
// src/rust/src/lib.rs
use extendr_api::prelude::*;

/// Fast GC content calculation
/// @param seq A character vector of DNA sequences
/// @export
#[extendr]
fn gc_content(seq: Vec<String>) -> Vec<f64> {
    seq.iter().map(|s| {
        let gc = s.chars().filter(|c| *c == 'G' || *c == 'C').count();
        gc as f64 / s.len() as f64
    }).collect()
}

// Auto-generates R wrapper: gc_content <- function(seq) .Call(wrap__gc_content, seq)
extendr_module! {
    mod mypackage;
    fn gc_content;
}
```

### Bioconductor + Rust

The Bioconductor community is actively exploring Rust integration. A dedicated Bioconductor Developer Forum session on "R + Rust in Bioinformatics" was held in July 2025. Real-world bioinformatics Rust crates (e.g. from the `noodles` ecosystem for BAM/CRAM/BCF I/O) can be wrapped and called from R packages via extendr.

### Comparison: C++ (Rcpp) vs Rust (extendr)

| Dimension | Rcpp / cpp11 | extendr (Rust) |
|-----------|-------------|----------------|
| Maturity | Very mature (~14k packages) | Early adoption (~dozens of packages) |
| Memory safety | Manual; undefined behaviour possible | Compile-time guaranteed |
| Parallelism | OpenMP / TBB | `rayon` — trivially safe |
| Compilation speed | Moderate | Slow (Cargo dependency resolution) |
| Bioconductor support | Full | Emerging |
| CRAN support | Full | Supported but requires Rust toolchain on build machines |
| Learning curve | C++ knowledge needed | Rust learning curve is steep |
| Best for | Numerical/matrix work, mature ecosystem | String-heavy parsing, safe concurrency, wrapping Rust crates |

---

## 5. Package Development Toolchain

The toolchain has been completely reorganised since 2017. What was a monolithic `devtools` is now a constellation of focused packages that `devtools` wraps.

### The Core Stack

```r
install.packages(c("devtools", "usethis", "roxygen2", "testthat",
                   "pkgdown", "covr", "pak"))
```

| Package | Role | Notes since 2017 |
|---------|------|-----------------|
| `devtools` | Meta-wrapper for the whole stack | Many functions moved to `usethis`; devtools re-exports them |
| `usethis` | Project/package setup automation | Did not exist in 2017 as a standalone package |
| `roxygen2` | Documentation generation | Now supports Markdown natively; R6 class docs; inline knitr chunks |
| `testthat` | Unit testing | Edition 3 (2020) — major redesign; snapshot testing added |
| `pkgdown` | Package website generation | Mature; GitHub Actions integration standard |
| `covr` | Test coverage | Integrates with `codecov.io` |
| `pak` | Fast dependency installer | Replacement for `install.packages()` for development |

### usethis — The Biggest Workflow Change

`usethis` did not exist as a standalone package in 2017. It automates nearly everything:

```r
# Start a new package
usethis::create_package("mypkg")

# Add infrastructure
usethis::use_testthat()          # sets up testthat
usethis::use_github_actions()    # CI/CD via GitHub Actions
usethis::use_pkgdown()           # package website
usethis::use_mit_license()       # license
usethis::use_vignette("intro")   # vignette scaffold
usethis::use_cpp11()             # C++ via cpp11
usethis::use_git()               # init git repo
usethis::use_github()            # push to GitHub, set remotes

# Situational awareness
usethis::git_sitrep()            # git/GitHub status report
```

### roxygen2 — Markdown Documentation

Since ~2019, `roxygen2` supports full Markdown in doc comments. This is now the recommended approach:

```r
#' Compute allele-specific copy number
#'
#' Runs FACETS on a matched tumor-normal pair. Returns a
#' `data.frame` with columns `chrom`, `seg`, `tcn`, `lcn`.
#'
#' @param tumor_bam Path to tumor BAM file
#' @param normal_bam Path to normal BAM file
#' @param snp_vcf Path to dbSNP VCF for pileup sites
#' @returns A `data.frame` with segmentation output
#' @examples
#' \dontrun{
#'   run_facets("tumor.bam", "normal.bam", "dbsnp.vcf")
#' }
#' @export
run_facets <- function(tumor_bam, normal_bam, snp_vcf) { ... }
```

Also new: **`roxygen2` can now document R6 classes directly** — `@description`, `@field`, and `@method` tags in an R6 class definition generate proper `.Rd` documentation.

### testthat Edition 3

`testthat` underwent a major redesign with Edition 3 (declared in `DESCRIPTION` as `Config/testthat/edition: 3`):

```r
# Snapshot testing — critical for testing output format stability
test_that("facets output has expected structure", {
  result <- run_facets(...)
  expect_snapshot(result)          # creates/compares a stored snapshot
  expect_snapshot_error(bad_call()) # snapshots error messages too
})

# New expectation functions
expect_no_error(safe_function())   # asserts no error thrown
expect_no_warning(clean_fn())      # asserts no warning
```

### GitHub Actions (CI/CD)

In 2017, CI for R packages meant `.travis.yml`. That is now essentially gone. The standard is `r-lib/actions` on GitHub Actions:

```yaml
# .github/workflows/R-CMD-check.yaml
on: [push, pull_request]
jobs:
  R-CMD-check:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
    steps:
      - uses: actions/checkout@v4
      - uses: r-lib/actions/setup-r@v2
      - uses: r-lib/actions/setup-r-dependencies@v2
      - uses: r-lib/actions/check-r-package@v2
```

`usethis::use_github_actions()` adds this file automatically. **Bioconductor packages use `bioc-check` action** which runs additional `BiocCheck` linting on top of R CMD check.

### pak — Fast Package Installation

```r
# Old
install.packages("SummarizedExperiment")
BiocManager::install("SummarizedExperiment")

# Modern — pak resolves dependencies in parallel, is much faster
pak::pak("bioc::SummarizedExperiment")
pak::pak("github::user/pkg")     # GitHub packages
pak::pak(c("dplyr", "bioc::DESeq2", "github::user/dev-pkg"))
```

---

## 6. OOP Class Systems

In 2017, you had S3, S4, and (rarely used) R5/Reference Classes. The landscape has expanded.

### The Current OOP Landscape

| System | Where | Paradigm | Best For |
|--------|-------|----------|----------|
| **S3** | Base R | Functional (generics) | Simple, idiomatic R; extending base functions |
| **S4** | Base R (`methods`) | Functional (formal) | Bioconductor packages; large collaborative codebases |
| **R5 / RC** | Base R | Encapsulated | Rarely used; effectively superseded by R6 |
| **R6** | `R6` package | Encapsulated | Mutable state; API wrappers; Shiny internals |
| **S7** | `S7` package (2022+) | Functional (modern) | **New recommendation** as S3/S4 replacement |

### S4 — Bioconductor Standard (Unchanged in Role, Updated in Practice)

S4 remains the **required class system for Bioconductor core data structures**. `SummarizedExperiment`, `GRanges`, `DNAStringSet`, `SingleCellExperiment` are all S4.

```r
# Defining an S4 class
setClass("VariantCall",
  contains = "SummarizedExperiment",  # inherit from SE
  slots = list(
    pipeline_version = "character",
    tumor_purity     = "numeric"
  )
)

# Defining S4 methods (generics)
setGeneric("tumorPurity", function(x) standardGeneric("tumorPurity"))
setMethod("tumorPurity", "VariantCall", function(x) x@tumor_purity)

# Validity checking
setValidity("VariantCall", function(object) {
  if (object@tumor_purity < 0 || object@tumor_purity > 1)
    "tumor_purity must be between 0 and 1"
  else TRUE
})
```

The `methods` package interface for S4 has not changed dramatically, but tooling around it has improved: `roxygen2` can now generate S4 documentation more cleanly, and `BiocCheck` enforces S4 correctness standards.

### R6 — Encapsulated OOP (Major Rise Since 2017)

R6 went from niche to mainstream between 2017 and 2025. It's used in `shiny`, `httpuv`, `plumber`, `httr2`, and many pipeline management packages. The reference semantics (modify-in-place) make it natural for stateful objects like database connections, pipeline state, or progress trackers.

```r
library(R6)

# Example: a pipeline runner with mutable state
PipelineRunner <- R6::R6Class("PipelineRunner",
  private = list(
    .log = NULL,
    .results = NULL
  ),
  public = list(
    initialize = function() {
      private$.log <- character(0)
      private$.results <- list()
    },
    run_step = function(name, fn, ...) {
      private$.log <- c(private$.log, paste0("[", Sys.time(), "] Running: ", name))
      private$.results[[name]] <- fn(...)
      invisible(self)          # enables method chaining: runner$run_step(...)$run_step(...)
    },
    get_results = function() private$.results,
    print = function(...) {
      cat("PipelineRunner with", length(private$.log), "steps completed\n")
    }
  )
)

runner <- PipelineRunner$new()
runner$run_step("facets", run_facets, tumor_bam, normal_bam)
runner$run_step("annotate", annotate_variants, runner$get_results()$facets)
```

**Key R6 concepts:**
- `public = list(...)` — accessible externally via `$`
- `private = list(...)` — internal only; good for internal state
- `active = list(...)` — computed properties (getter/setter syntax)
- `self` — reference to current object (like Python's `self`)
- `super` — access parent class methods in inheritance
- Objects are **mutable by reference** — no copy-on-modify

### S7 — The New Unified System (2022+)

S7 is designed as the eventual successor to both S3 and S4, unifying their best features:

```r
library(S7)

# Define a class
SegmentResult <- new_class("SegmentResult",
  properties = list(
    chrom    = class_character,
    start    = class_integer,
    end      = class_integer,
    tcn      = class_double,         # total copy number
    lcn      = class_double,         # lesser copy number
    purity   = new_property(
      class = class_double,
      validator = function(self, value) {
        if (!is.na(value) && (value < 0 || value > 1))
          "purity must be between 0 and 1"
      }
    )
  )
)

# Define generics and methods
loh_status <- new_generic("loh_status", "x")
method(loh_status, SegmentResult) <- function(x) x@lcn == 0

# S7 is backward compatible with S3 generics
method(print, SegmentResult) <- function(x, ...) {
  cat(sprintf("Segment %s:%d-%d | TCN=%.1f LCN=%.1f\n",
              x@chrom, x@start, x@end, x@tcn, x@lcn))
}
```

**S7 key properties:**
- S7 objects **are** S3 objects — backward compatible
- Properties replace slots; can be dynamic (computed)
- Works with both S3 and S4 generics
- Cannot extend existing S4 classes (S7 → S4 inheritance not supported)
- **Not yet in base R** — currently a CRAN package, but intended for eventual base R inclusion

### When to Use Which

```
Writing a new Bioconductor package?
  └─ S4 (required for Bioconductor core; expected by reviewers)

Writing a new CRAN package with formal OOP?
  └─ S7 (modern, forward-looking) or R6 (if you need mutable state)

Writing an internal analysis package / pipeline framework?
  └─ R6 for stateful objects (pipeline runners, DB connections)
  └─ S3 for simple typed results with print/summary methods

Extending existing Bioconductor classes?
  └─ S4 (inherit from SummarizedExperiment, GRanges, etc.)

Need Python-like OOP for a data scientist collaborator to understand?
  └─ R6 (most familiar to Python/Java programmers)
```

---

## 7. Quick Decision Guide

### Language Feature Quick Reference

| You want to... | Old way (2017) | Modern way (2025+) |
|---------------|---------------|-------------------|
| Pipe operations | `%>%` (magrittr) | `\|>` (native) or `%>%` still fine |
| Anonymous functions | `function(x) x + 1` | `\(x) x + 1` |
| Load a string without factor coercion | `read.csv(..., stringsAsFactors=FALSE)` | Default; no flag needed |
| Write inline string with backslashes | Escape hell | `r"(C:\Users\path)"` |
| Install packages fast | `install.packages()` | `pak::pak()` |

### Bioconductor Workflow Quick Reference

| Task | Key Package(s) |
|------|---------------|
| Store RNA-seq data | `SummarizedExperiment` |
| Store scRNA-seq data | `SingleCellExperiment` + `scater` + `scran` |
| Store multi-omic patient data | `MultiAssayExperiment` |
| Store spatial transcriptomics | `SpatialExperiment` |
| Differential expression (bulk) | `DESeq2` (with `lfcShrink`), `limma::voom` |
| Differential expression (single-cell) | `muscat` (pseudobulk), `MAST` |
| CNV from single-cell | `infercnv`, `copykit` |
| Cell type annotation | `SingleR` + `celldex` |
| Genome annotations | `AnnotationHub`, `EnsDb.*` packages |
| Large matrix (out-of-memory) | `HDF5Array` + `DelayedArray` |
| Cloud-scale analysis | `AnVIL` (Terra/GCP) |

### Extension Language Quick Reference

| Use Case | Recommendation |
|----------|---------------|
| Numerical / matrix operations | Rcpp + RcppArmadillo |
| String-heavy parsing (BAM/VCF-like) | Rust (extendr) or Rcpp |
| New package, want modern C++ | cpp11 |
| Safe parallelism | Rust (rayon) |
| Wrapping existing C++ libraries | Rcpp |
| Wrapping existing Rust crates | extendr |

---

*This document synthesizes changes across R 3.4 → 4.5+, Bioconductor 3.5 → 3.22, and the surrounding ecosystem. For ongoing changes, follow:*
- *https://developer.r-project.org — R release notes*
- *https://bioconductor.org/news — Bioconductor release notes (biannual)*
- *https://extendr.rs — Rust/R ecosystem*
- *https://r-pkgs.org — R Packages (2e) by Wickham & Bryan*
- *https://adv-r.hadley.nz — Advanced R (2e) for OOP deep dives*

# R UI Development: 2012 to 2026
### From Raw `fluidPage()` to bslib Dashboards, Modules, and Serverless WebAssembly

> **Scope:** Shiny 0.x → 1.9+ · Bootstrap 2 → 5 · shinydashboard → bslib · R Markdown → Quarto · shinylive / webR  
> **Shiny has not been replaced** — it remains the dominant R UI framework, but has been substantially modernised.

---

## Table of Contents

1. [Era 1 (2012–2016): Early Shiny](#1-era-1-20122016-early-shiny)
2. [Era 2 (2016–2020): Modules, Dashboards, htmlwidgets](#2-era-2-20162020-modules-dashboards-htmlwidgets)
3. [Era 3 (2020–2023): bslib, Bootstrap 5, Quarto](#3-era-3-20202023-bslib-bootstrap-5-quarto)
4. [Era 4 (2023–2026): Serverless, ExtendedTask, AI Integration](#4-era-4-20232026-serverless-extendedtask-ai-integration)
5. [Current Best Practices (2026)](#5-current-best-practices-2026)
6. [Complete Example: DNA → Amino Acid Translator](#6-complete-example-dna--amino-acid-translator)
7. [Deployment Quick Reference](#7-deployment-quick-reference)

---

## 1. Era 1 (2012–2016): Early Shiny

Shiny launched in 2012. The programming model — reactive expressions, `ui.R` / `server.R` split, Bootstrap 2 — was revolutionary for statisticians who had never written a web app. But the aesthetics and ergonomics were rough.

### The Original App Structure

```
myapp/
├── ui.R       # UI definition
└── server.R   # server logic
```

Or the single-file `app.R` format (introduced ~2015):

```r
# app.R — the single-file format that became standard
library(shiny)

ui <- fluidPage(
  titlePanel("My App"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Sample size", 1, 100, 50)
    ),
    mainPanel(
      plotOutput("plot")
    )
  )
)

server <- function(input, output) {
  output$plot <- renderPlot({
    hist(rnorm(input$n))
  })
}

shinyApp(ui, server)
```

### What Was Painful in This Era

- **No theming beyond CSS hacks** — Bootstrap 2, grey sidebars, blue buttons. Customisation meant writing raw CSS.
- **No modules** — large apps became a single enormous `server` function with colliding input/output IDs.
- **Reactivity debugging was opaque** — no `reactlog`, errors surfaced poorly.
- **No async** — any slow computation blocked the entire session.
- **`shinydashboard` (2015)** was the first real improvement to layout — boxes, value boxes, sidebars with icons — and became immediately ubiquitous. Most bioinformatics Shiny apps built 2015–2020 use it.

```r
library(shinydashboard)

ui <- dashboardPage(
  dashboardHeader(title = "Genomics Dashboard"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("CNV", tabName = "cnv", icon = icon("dna"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "cnv",
        box(plotOutput("cnv_plot"), width = 12)
      )
    )
  )
)
```

---

## 2. Era 2 (2016–2020): Modules, Dashboards, htmlwidgets

### Shiny Modules (2016) — The Most Important Design Pattern

Modules solved the namespace collision problem in large apps. A module is a UI function + server function pair with an isolated ID namespace. This is now **mandatory practice** for any non-trivial app.

```r
# A self-contained module for sequence input
# --- module definition ---
sequenceInputUI <- function(id) {
  ns <- NS(id)   # NS() creates the namespace function
  tagList(
    textAreaInput(ns("seq"),  "Enter DNA sequence", rows = 3,
                  placeholder = "ATGCGA..."),
    actionButton(ns("submit"), "Translate", class = "btn-primary")
  )
}

sequenceInputServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # input$seq refers to "mymod-seq" in the parent app — no collision
    eventReactive(input$submit, {
      toupper(trimws(input$seq))
    })
  })
}

# --- use in app ---
ui <- fluidPage(
  sequenceInputUI("mymod")   # ID prefix is "mymod"
)

server <- function(input, output, session) {
  dna <- sequenceInputServer("mymod")
}
```

**The rule:** any UI component you expect to reuse, or any component with >~5 inputs/outputs, should be a module.

### htmlwidgets — JavaScript Visualisation in R

The `htmlwidgets` framework (2014–2015, mature by 2016) allowed R wrappers around any JavaScript visualisation library. The bioinformatics-relevant ones:

| Package | JS library | Use |
|---------|-----------|-----|
| `plotly` | Plotly.js | Interactive plots, hover data |
| `leaflet` | Leaflet.js | Genomic coordinate maps, coverage plots |
| `DT` | DataTables | Interactive sortable/filterable tables |
| `visNetwork` | vis.js | Network/pathway graphs |
| `igvShiny` | IGV.js | In-browser genome browser |

```r
library(DT)
library(plotly)

# DT — interactive variant tables
output$variant_table <- renderDT({
  datatable(variant_df,
    filter   = "top",        # per-column filters
    rownames = FALSE,
    options  = list(pageLength = 25, scrollX = TRUE)
  )
})

# plotly — interactive coverage plot
output$coverage_plot <- renderPlotly({
  p <- ggplot(coverage_df, aes(pos, depth, color = allele)) +
    geom_line()
  ggplotly(p, tooltip = c("pos", "depth"))
})
```

### `reactlog` — Reactive Debugging

`reactlog` arrived in this era, providing a visual graph of reactive dependencies. Still the primary debugging tool for reactivity problems:

```r
options(shiny.reactlog = TRUE)
shinyApp(ui, server)
# Then press Ctrl+F3 to view the reactive graph
```

---

## 3. Era 3 (2020–2023): bslib, Bootstrap 5, Quarto

### The Core Shift: `bslib` Replaces `shinydashboard`

`bslib` was released to CRAN in 2021 and is now the **recommended way to build all Shiny UIs**. It wraps Bootstrap (currently version 5) and provides a modern layout and theming system. `shinydashboard` is effectively in maintenance mode — new projects should use `bslib`.

**The key mental model change:** instead of `fluidPage()` + `sidebarLayout()`, you now use `page_*()` functions from `bslib`:

| Old function | New bslib equivalent |
|-------------|---------------------|
| `fluidPage()` | `page_fluid()` |
| `navbarPage()` | `page_navbar()` |
| `sidebarLayout()` | `page_sidebar()` |
| `dashboardPage()` (shinydashboard) | `page_sidebar()` + `layout_columns()` |
| `tabsetPanel()` | `navset_card_tab()` |

### Modern Layout with `bslib`

```r
library(shiny)
library(bslib)

ui <- page_sidebar(
  title = "Genomics Explorer",
  theme = bs_theme(
    preset    = "shiny",     # modern default look
    bootswatch = "flatly",  # or any Bootswatch theme
    base_font = font_google("Inter")
  ),

  # --- sidebar ---
  sidebar = sidebar(
    title  = "Parameters",
    width  = 280,
    fileInput("vcf",  "Upload VCF"),
    selectInput("gene", "Gene", choices = NULL),
    input_switch("loh_filter", "Show LOH only", value = FALSE),
    hr(),
    actionButton("run", "Run Analysis", class = "btn-primary w-100")
  ),

  # --- main panel: responsive columns ---
  layout_columns(
    col_widths = c(8, 4),

    card(
      full_screen = TRUE,
      card_header("CNV Profile"),
      plotOutput("cnv_plot", height = "400px")
    ),

    card(
      card_header("Summary Statistics"),
      value_box(
        title    = "Tumor Purity",
        value    = textOutput("purity"),
        showcase = bsicons::bs_icon("activity"),
        theme    = "primary"
      ),
      value_box(
        title    = "Ploidy",
        value    = textOutput("ploidy"),
        showcase = bsicons::bs_icon("layers"),
        theme    = "secondary"
      )
    )
  ),

  card(
    full_screen = TRUE,
    card_header("Variant Table"),
    DT::DTOutput("var_table")
  )
)
```

### Theming — Real-Time with `bs_themer()`

```r
server <- function(input, output, session) {
  # During development: live theme editor in a sidebar panel
  bs_themer()

  # ... rest of server
}
```

### Quarto — Replacing R Markdown for Documents and Dashboards

Quarto (2022) replaced R Markdown as the standard for reproducible documents. For Shiny specifically, Quarto enables **inline Shiny apps in documents** and the `quarto-dashboard` format:

```yaml
# _quarto.yml for a Quarto dashboard with Shiny
format:
  dashboard:
    theme: flatly
server: shiny
```

```{r}
#| context: server
output$plot <- renderPlot({ ... })
```

For static, no-server reports embedding interactive plots — `plotly` + Quarto HTML is now the standard, replacing `flexdashboard`.

---

## 4. Era 4 (2023–2026): Serverless, ExtendedTask, AI Integration

### `ExtendedTask` — Non-Blocking Async (Shiny 1.8+, 2024)

Before `ExtendedTask`, running a slow computation (e.g. a FACETS run, an alignment call) blocked the entire Shiny session — the UI froze. `promises` + `future` offered async, but were complex to wire up correctly. `ExtendedTask` makes this simple:

```r
library(shiny)
library(bslib)
library(future)
plan(multisession)   # background R processes

ui <- page_sidebar(
  sidebar = sidebar(
    actionButton("run", "Run FACETS"),
    # input_task_button gives visual feedback + prevents double-click
    bslib::input_task_button("run_task", "Run (non-blocking)")
  ),
  card(verbatimTextOutput("result"))
)

server <- function(input, output, session) {

  # ExtendedTask: runs in background, UI stays responsive
  facets_task <- ExtendedTask$new(function(tumor_bam, normal_bam) {
    future({
      run_facets_pipeline(tumor_bam, normal_bam)  # slow operation
    })
  }) |> bslib::bind_task_button("run_task")

  observeEvent(input$run_task, {
    facets_task$invoke(tumor_bam = "tumor.bam", normal_bam = "normal.bam")
  })

  output$result <- renderPrint({
    facets_task$result()   # reactive; updates when task completes
  })
}
```

### shinylive — Shiny Without a Server

`shinylive` is a new way to run Shiny entirely in the browser, without any need for a hosted server, using WebAssembly via the webR project.

This means a Shiny app can be deployed to **GitHub Pages, Netlify, or any static host** — no Shiny Server or Posit Connect required.

```r
# Export a local Shiny app as a static site
shinylive::export("myapp/", "site/")

# Deploy site/ to GitHub Pages → app runs in user's browser
# No server costs, no session limits
```

**Constraints to know:**
- No raw socket access (no `curl` — HTTP via XHR only)
- Not all packages available; must be pre-compiled WebAssembly binaries
- For a complex Shinylive app, load times decreased from over a minute to just 15 seconds with recent optimisations — startup is still slower than a traditional server app
- All data and code are visible to the client — not suitable for PHI or proprietary data

**In Quarto documents:**

````markdown
```{shinylive-r}
#| standalone: true
library(shiny)
library(bslib)

ui <- page_fluid(
  textInput("seq", "DNA Sequence"),
  verbatimTextOutput("aa")
)

server <- function(input, output, session) {
  output$aa <- renderText({ translate_dna(input$seq) })
}

shinyApp(ui, server)
```
````

### New UI Inputs (bslib 2023–2025)

| New input | Replaces | Notes |
|-----------|---------|-------|
| `input_switch()` | `checkboxInput()` | Toggle switch visual |
| `input_task_button()` | `actionButton()` | Shows spinner during task, prevents re-click |
| `tooltip()` | Custom JS | Info icons with hover text |
| `popover()` | Custom JS | Click-triggered overlay with controls |
| `accordion()` | Manual `div` folding | Collapsible input groups |
| `layout_columns()` | `column()` + `fluidRow()` | Responsive CSS grid |
| `value_box()` | `infoBox()` (shinydashboard) | Modern KPI cards |

### `updateOn = "blur"` for Text Inputs

A small but important UX change: text inputs can now defer reactivity until the user finishes typing, preventing constant re-firing on every keystroke — critical for sequence inputs:

```r
textAreaInput("seq", "DNA sequence", updateOn = "blur")
# Reactive update fires only when user clicks away or presses Enter
```

---

## 5. Current Best Practices (2026)

### App Structure — Modules from the Start

```
myapp/
├── app.R              # shinyApp(ui, server) — thin wrapper only
├── R/
│   ├── mod_input.R    # sequenceInputUI() + sequenceInputServer()
│   ├── mod_results.R  # resultsUI() + resultsServer()
│   └── utils.R        # translate_dna(), validate_sequence(), etc.
├── www/
│   └── custom.css     # minimal custom CSS (rarely needed with bslib)
└── tests/
    └── testthat/      # shinytest2 for UI testing
```

### The `app.R` Pattern

```r
# app.R — should be as thin as possible
library(shiny)
library(bslib)

# Source all module and utility files
box::use(R/mod_input, R/mod_results, R/utils)  # or source() / box

ui <- page_navbar(
  title = "My Bioinformatics App",
  theme = bs_theme(preset = "shiny"),
  nav_panel("Translate", mod_input$ui("translate")),
  nav_panel("About",     includeMarkdown("README.md"))
)

server <- function(input, output, session) {
  mod_input$server("translate")
}

shinyApp(ui, server)
```

### Testing with `shinytest2`

`shinytest2` (replaced `shinytest` ~2022) uses a headless Chromium browser to snapshot app state:

```r
# tests/testthat/test-translation.R
library(shinytest2)

test_that("DNA translation works end-to-end", {
  app <- AppDriver$new(app_dir = ".")

  app$set_inputs(seq = "ATGAAATGA")
  app$click("submit")
  app$expect_values(output = "aa_output")  # snapshot comparison
})
```

### Dependency Management with `renv`

```r
renv::init()        # initialise lockfile
renv::snapshot()    # record current package versions
renv::restore()     # reproduce environment on another machine
```

---

## 6. Complete Example: DNA → Amino Acid Translator

A full worked example demonstrating all current best practices: `bslib` layout, modules, `updateOn = "blur"`, `tooltip()`, `value_box()`, Bioconductor translation via `Biostrings`.

### File structure

```
dna_translator/
├── app.R
└── R/
    ├── mod_translator.R
    └── utils_biology.R
```

---

### `R/utils_biology.R`

```r
# Genetic code and translation utilities
# Uses base R only — no Biostrings dependency for portability.
# For production: use Biostrings::translate() on a DNAString object.

CODON_TABLE <- c(
  TTT="F", TTC="F", TTA="L", TTG="L",
  CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M",
  GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S",
  CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T",
  GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K",
  GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W",
  CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

#' Validate a DNA string
#' @param seq character — raw input from user
#' @return list(valid = logical, cleaned = character, message = character)
validate_sequence <- function(seq) {
  cleaned <- toupper(gsub("[^ATCGatcg]", "", seq))

  if (nchar(cleaned) == 0)
    return(list(valid = FALSE, cleaned = "", message = "No valid DNA bases found."))

  if (nchar(cleaned) < 3)
    return(list(valid = FALSE, cleaned = cleaned,
                message = "Sequence too short — need at least one codon (3 bp)."))

  remainder <- nchar(cleaned) %% 3
  if (remainder != 0)
    message <- sprintf("Sequence trimmed: last %d base(s) ignored (not a complete codon).", remainder)
  else
    message <- NULL

  list(valid = TRUE, cleaned = cleaned, message = message)
}

#' Translate a cleaned DNA string to amino acids
#' @param dna character — uppercase, only ATCG, length divisible by 3
#' @param stop_at_stop logical — whether to terminate at first stop codon
#' @return list(aa = character, n_codons = int, n_aa = int, has_stop = logical)
translate_dna <- function(dna, stop_at_stop = FALSE) {
  n_codons <- nchar(dna) %/% 3
  codons   <- substring(dna, seq(1, n_codons * 3 - 2, 3),
                             seq(3, n_codons * 3,     3))
  aa_vec   <- CODON_TABLE[codons]
  aa_vec[is.na(aa_vec)] <- "?"

  has_stop <- any(aa_vec == "*")

  if (stop_at_stop && has_stop)
    aa_vec <- aa_vec[seq_len(which(aa_vec == "*")[1])]

  aa_string <- paste(aa_vec, collapse = "")

  list(
    aa       = aa_string,
    codons   = codons,
    n_codons = n_codons,
    n_aa     = sum(aa_vec != "*"),
    has_stop = has_stop
  )
}
```

---

### `R/mod_translator.R`

```r
library(shiny)
library(bslib)
library(bsicons)

# ── UI ─────────────────────────────────────────────────────────────────────────

translatorUI <- function(id) {
  ns <- NS(id)

  layout_columns(
    col_widths = c(5, 7),

    # Left card: inputs
    card(
      card_header(
        "Input Sequence",
        tooltip(
          bs_icon("info-circle"),
          "Paste a raw DNA sequence. Non-ATCG characters are silently removed.
           Sequence is read 5' → 3' from the first base."
        )
      ),

      textAreaInput(
        ns("seq"),
        label     = NULL,
        rows      = 6,
        width     = "100%",
        value     = "ATGAAAGCAATTTTTCAGTTCGAGTAA",   # example
        placeholder = "Paste DNA sequence here...",
        updateOn  = "blur"   # don't fire on every keystroke
      ),

      input_switch(ns("stop_early"), "Stop at first stop codon (*)", value = TRUE),

      layout_columns(
        col_widths = c(6, 6),
        actionButton(ns("example"), "Load Example",
                     class = "btn-outline-secondary w-100"),
        input_task_button(ns("translate"), "Translate →",
                          class = "btn-primary w-100")
      ),

      uiOutput(ns("validation_msg"))
    ),

    # Right card: results
    card(
      card_header("Translation Output"),
      card_body(
        layout_columns(
          col_widths = c(4, 4, 4),
          value_box(
            title    = "Codons",
            value    = textOutput(ns("n_codons")),
            showcase = bs_icon("boxes"),
            theme    = "primary"
          ),
          value_box(
            title    = "Amino Acids",
            value    = textOutput(ns("n_aa")),
            showcase = bs_icon("link-45deg"),
            theme    = "secondary"
          ),
          value_box(
            title    = "Stop Codon",
            value    = textOutput(ns("has_stop")),
            showcase = bs_icon("octagon"),
            theme    = "warning"
          )
        ),

        hr(),

        h6("Amino Acid Sequence (single-letter code)"),
        verbatimTextOutput(ns("aa_out")),

        hr(),

        h6("Codon Table"),
        DT::DTOutput(ns("codon_table"), height = "250px")
      )
    )
  )
}

# ── Server ──────────────────────────────────────────────────────────────────────

translatorServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Load example sequence
    observeEvent(input$example, {
      updateTextAreaInput(session, "seq",
        value = "ATGAAAGCAATTTTTCAGTTCGAATTTGCCCGCGATACGCATAA"
      )
    })

    # Core reactive: only fires on button click or blur
    result <- eventReactive(input$translate, {
      req(nchar(trimws(input$seq)) > 0)

      raw   <- trimws(input$seq)
      check <- validate_sequence(raw)

      if (!check$valid) return(list(error = check$message))

      tr <- translate_dna(check$cleaned, stop_at_stop = input$stop_early)
      tr$warning <- check$message   # may be NULL
      tr$valid   <- TRUE
      tr
    }, ignoreNULL = FALSE)

    # Validation message
    output$validation_msg <- renderUI({
      r <- result()
      if (is.null(r)) return(NULL)
      if (!isTRUE(r$valid))
        return(div(class = "alert alert-danger mt-2", r$error))
      if (!is.null(r$warning))
        return(div(class = "alert alert-warning mt-2", r$warning))
      NULL
    })

    # Value boxes
    output$n_codons <- renderText({
      r <- result(); req(isTRUE(r$valid)); r$n_codons
    })
    output$n_aa <- renderText({
      r <- result(); req(isTRUE(r$valid)); r$n_aa
    })
    output$has_stop <- renderText({
      r <- result(); req(isTRUE(r$valid))
      if (r$has_stop) "Yes ✓" else "None"
    })

    # Amino acid string — wrap at 60 characters for readability
    output$aa_out <- renderText({
      r <- result(); req(isTRUE(r$valid))
      aa <- r$aa
      paste(
        substring(aa, seq(1, nchar(aa), 60),
                      pmin(seq(60, nchar(aa) + 59, 60), nchar(aa))),
        collapse = "\n"
      )
    })

    # Codon table
    output$codon_table <- DT::renderDT({
      r <- result(); req(isTRUE(r$valid))

      df <- data.frame(
        Position = seq_along(r$codons),
        Codon    = r$codons,
        AA       = CODON_TABLE[r$codons],
        Type     = ifelse(CODON_TABLE[r$codons] == "*", "Stop",
                   ifelse(CODON_TABLE[r$codons] %in% c("M"),  "Start/Met",
                                                               "Coding")),
        stringsAsFactors = FALSE
      )

      DT::datatable(df,
        rownames  = FALSE,
        selection = "none",
        options   = list(
          pageLength = 10,
          dom        = "tp",   # table + pagination only
          scrollX    = TRUE
        )
      ) |>
      DT::formatStyle("Type",
        backgroundColor = DT::styleEqual(
          c("Stop", "Start/Met", "Coding"),
          c("#ffe0e0",  "#e0f0ff",  "white")
        )
      )
    })

  })
}
```

---

### `app.R`

```r
library(shiny)
library(bslib)

# Source utilities and module
source("R/utils_biology.R")
source("R/mod_translator.R")

ui <- page_navbar(
  title = "DNA → Protein Translator",
  theme = bs_theme(
    preset     = "shiny",
    bootswatch = "flatly",
    base_font  = font_google("Inter")
  ),

  nav_panel(
    title = "Translator",
    icon  = bsicons::bs_icon("dna"),
    padding = "1.5rem",
    translatorUI("tr1")
  ),

  nav_panel(
    title = "About",
    icon  = bsicons::bs_icon("info-circle"),
    padding = "1.5rem",
    card(
      card_body(
        h4("DNA → Amino Acid Translator"),
        p("Reads sequence 5' → 3', uses standard genetic code (NCBI table 1)."),
        p("For production use, replace the codon table logic with:"),
        tags$pre(
"library(Biostrings)
dna <- DNAString('ATGAAATGA')
aa  <- translate(dna)        # returns AAString"
        )
      )
    )
  ),

  nav_spacer(),
  nav_item(
    tags$a(
      bsicons::bs_icon("github"), "Source",
      href   = "https://github.com/yourname/dna-translator",
      target = "_blank"
    )
  )
)

server <- function(input, output, session) {
  translatorServer("tr1")
}

shinyApp(ui, server)
```

### Install Dependencies

```r
install.packages(c("shiny", "bslib", "bsicons", "DT"))
# Optional: for production Biostrings-backed translation
BiocManager::install("Biostrings")
```

### Run

```r
shiny::runApp("dna_translator/")
```

### Export as Serverless (shinylive)

```r
library(shinylive)
shinylive::export("dna_translator/", "site/")
# → Deploy site/ to GitHub Pages; app runs in user's browser via WebAssembly
```

---

## 7. Deployment Quick Reference

| Option | When to use | Cost |
|--------|------------|------|
| **shinyapps.io** (Posit) | Quick sharing, small apps | Free tier available |
| **Posit Connect** | Enterprise, auth, scheduling | Licensed |
| **ShinyProxy** | Docker-based, multi-user, self-hosted | Infrastructure cost only |
| **shinylive → GitHub Pages** | Public apps, no PHI, no heavy compute | Free |
| **shinylive → Netlify** | Custom domains, static hosting | Free tier generous |
| **Docker + nginx** | Full control, any cloud | Infrastructure cost |

### Checklist Before Deploying

```
□ renv::snapshot() — lockfile committed
□ All inputs validated (req(), validate(), need())
□ No hardcoded file paths — use here::here() or app-relative paths
□ Long operations use ExtendedTask (not blocking renderPlot)
□ Modules used — no naked server-level input$ IDs at scale
□ shinytest2 snapshot tests pass
□ No PHI/credentials in source — use environment variables or config::get()
```

---

*Companion documents: `R_Bioconductor_Changes_Since_2017.md` · `R_Bioconductor_Changes_2014_2017.md`*  
*Primary references: https://shiny.posit.co · https://rstudio.github.io/bslib · https://posit-dev.github.io/r-shinylive*

