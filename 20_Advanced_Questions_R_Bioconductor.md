# 20 Advanced Questions: R, Bioconductor & Computational Biology
### Beyond the README — Applied, Architectural, and Statistical Depth

> These questions go **outside the scope of the document** and assume familiarity with its contents.  
> Topics: advanced reactivity, memory architecture, statistical edge cases, Bioconductor internals,  
> genomic analysis design, production R, and computational biology judgment calls.  
> Answers are at the bottom.

---

## Questions

**Q1.**
R uses a copy-on-modify memory model. Consider this code running in a Shiny server function:

```r
server <- function(input, output, session) {
  big_matrix <- readRDS("500mb_count_matrix.rds")

  output$plot1 <- renderPlot({ plot(big_matrix[,1]) })
  output$plot2 <- renderPlot({ plot(big_matrix[,2]) })
}
```

How many copies of `big_matrix` exist in memory when both outputs are being computed simultaneously across two browser sessions? What R mechanism governs this, and what is the correct architectural pattern to avoid it in a production multi-user Shiny app?

---

**Q2.**
You are designing an S4 class hierarchy for a clinical sequencing pipeline. You have `TumorSample` extending `SummarizedExperiment`, and `PairedSample` extending `TumorSample` (adding a `normal_bam` slot). A colleague wants to write a single `QCReport` method that dispatches differently when the input is paired vs. unpaired. They propose `setMethod("QCReport", "TumorSample", ...)` and `setMethod("QCReport", "PairedSample", ...)`. What is the S4 dispatch rule that makes this work, what is the term for it, and what happens if someone passes an object of a class that inherits from neither?

---

**Q3.**
`DESeq2` uses a negative binomial model with dispersion estimates. Explain what happens statistically when you run a DESeq2 analysis on a dataset with only 2 samples per condition (n=2 per group), why the dispersion estimates are unreliable in this scenario, and what the `DESeq2` documentation recommends as the minimum sample size for reliable results. What specific dispersion estimation strategy does DESeq2 use to partially compensate for small n?

---

**Q4.**
In the context of FACETS output you know well: what is `dipLogR` and why is it a more useful reference point than log2(tumor/normal) = 0 for interpreting copy number segments? In which specific tumor type scenario would `dipLogR` deviate most dramatically from 0, and what does a strongly negative `dipLogR` indicate about the sample?

---

**Q5.**
You have a `SingleCellExperiment` with 50,000 cells and 30,000 genes stored in a sparse matrix. You apply `logNormCounts()` from `scater`, then call `runPCA()`. A reviewer asks why you did not use raw counts for PCA directly. Give the statistical argument for log-normalisation before PCA in single-cell data, specifically addressing: (a) the mean-variance relationship in count data, (b) what happens to PCA when this relationship is not addressed, and (c) why `log1p` (log of count + 1) is used rather than `log`.

---

**Q6.**
Shiny's reactive graph has three fundamental primitives: reactive sources, reactive conductors, and reactive endpoints. A colleague writes:

```r
server <- function(input, output, session) {
  output$result <- renderText({
    Sys.sleep(3)
    paste("Result:", expensive_computation(input$x))
  })
  output$summary <- renderText({
    paste("Summary:", expensive_computation(input$x))
  })
}
```

`expensive_computation()` is called twice whenever `input$x` changes. What is the correct reactive primitive to use to cache the result, why does simply assigning `val <- expensive_computation(input$x)` inside the server function not work, and what is the invalidation rule for your chosen solution?

---

**Q7.**
`AnnotationHub` and `ExperimentHub` are described as modern replacements for manual GTF downloading. However, there is a scenario where using `AnnotationHub` in a Bioconductor package introduces a reproducibility problem for users. Describe the problem, what causes it, and the two approaches Bioconductor packages use to address it — one using a specific package type, one using caching.

---

**Q8.**
You are doing pseudobulk differential expression with `muscat` on a single-cell dataset from a clinical trial: 15 patients, 8 treated and 7 control, each with 3,000–15,000 cells. A postdoc suggests analysing at the single-cell level instead (treating each cell as an independent observation, ~150,000 observations total) for more statistical power. Give three specific reasons why this approach is statistically invalid, referencing the appropriate statistical concept for each.

---

**Q9.**
In Rcpp, what is the difference between `Rcpp::NumericVector` and `std::vector<double>`? When would you prefer one over the other inside an `[[Rcpp::export]]` function, and what is the performance cost of converting between them? What is "Rcpp Sugar" and give one example of where it avoids this conversion entirely?

---

**Q10.**
You are building a Bioconductor package and submitting for the first time. Beyond passing `R CMD check` with no errors, what are four specific additional requirements that `BiocCheck` enforces that CRAN does not? Your package wraps an existing genomics tool — what does Bioconductor require regarding the tool's availability that CRAN does not?

---

**Q11.**
In copy number analysis of tumors, the concepts of tumor purity (cellularity) and ploidy are confounded. Describe the mathematical relationship between these two variables that creates the confounding, and explain why a low-purity high-ploidy sample can produce an identical FACETS output (in terms of observed log2 ratios) to a high-purity low-ploidy sample. What is the term for these indistinguishable solutions?

---

**Q12.**
`HDF5Array` and `DelayedArray` are recommended for large single-cell matrices. However, there is a class of operations where `DelayedArray` performs dramatically worse than an in-memory matrix, even when the data fits in RAM. Describe the access pattern that causes this, why it happens at the HDF5 storage level, and what the `DelayedArray` package's `blockApply()` function is designed to do about it.

---

**Q13.**
`renv` locks package versions using a `renv.lock` file. You are using a Bioconductor package that depends on a specific version of Bioconductor itself (e.g. `BiocGenerics 0.44.0` which only exists in Bioconductor 3.16). A colleague with a different Bioconductor version runs `renv::restore()` and the install fails. What is the fundamental limitation of `renv` with respect to Bioconductor version management, and what additional tool or approach is needed to fully reproduce a Bioconductor-dependent environment?

---

**Q14.**
R's native `|>` pipe introduced in 4.1 evaluates its right-hand side differently from magrittr's `%>%` in one technically important way related to **function calls vs. expressions**. Specifically: `x %>% {. * 2}` is valid magrittr syntax. Why does this not work with `|>`, what is the underlying reason (related to how `|>` desugars), and what is the `|>`-compatible way to apply an anonymous transformation?

---

**Q15.**
You are integrating WES and RNA-seq data for 200 TCGA patients using `MultiAssayExperiment`. Not all patients have both data types — 180 have WES, 160 have RNA-seq, and 140 have both. When you call `intersectColumns()` on the MAE, what does it return and how many patients are in the result? Then describe the function that would give you all patients who have *at least one* data type, and how missing data is represented in a MAE for a patient who has WES but not RNA-seq.

---

**Q16.**
In the context of clonal hematopoiesis (CH) detection from clinical sequencing data: explain why the standard somatic variant calling pipeline (designed for tumor-normal pairs) is poorly suited for CH detection, what specific variant allele frequency (VAF) range CH variants typically fall in, and why that range creates specific statistical challenges for variant callers. What is the key statistical concept — borrowed from signal detection theory — that governs the tradeoff in tuning a CH caller?

---

**Q17.**
`Rust` via `extendr` is described as having "slow compilation (Cargo dependency resolution)." This is a real practical problem in a CRAN package context — CRAN checks packages on submission and the build must complete within time limits. Describe two strategies that Rust-based R packages use to make CRAN submission feasible, one involving pre-compiled binaries and one involving the Cargo build configuration.

---

**Q18.**
You are writing a Shiny app that queries a PostgreSQL database with patient genomic data. 100 concurrent users each trigger queries on a table with 50M rows. The naive approach uses a single database connection stored in the global environment. Describe three specific problems this causes in a multi-user Shiny context, and describe the correct architectural pattern using a specific R package designed for this — including what it provides that a single global connection does not.

---

**Q19.**
`S7` properties can be "dynamic" — computed rather than stored. This is listed as a key difference from S4 slots. Write a concrete example of a dynamic S7 property that would be useful in a genomics class — specifically one that computes a value from other stored properties — and explain why the same thing in S4 requires a different (and more verbose) pattern involving a separate generic and method.

---

**Q20.**
In a Shiny app using `bslib`, you notice that when a user resizes their browser window, your `layout_columns()` grid collapses from a 3-column layout to a single column. This is expected Bootstrap responsive behaviour. However, for a variant table that needs all columns visible at all times regardless of screen size, this is a problem. Describe: (a) why Bootstrap does this, (b) two approaches at the Shiny/bslib level to prevent column collapse for a specific component, and (c) the CSS class from Bootstrap 5 you would add manually to lock a layout to a specific number of columns on all screen sizes.

---
---
---

## Answers

**A1.**
The answer depends on R's **copy-on-modify** (also called copy-on-write) semantics. When `big_matrix` is bound to a name and passed into render functions, R does not immediately copy it — it increments a reference count. A copy is only made when the object is *modified*. Since `renderPlot` only reads the matrix (no modification), both `output$plot1` and `output$plot2` share the **same underlying memory block** within a single session. However, across **two browser sessions**, `server()` is called twice — creating two independent R environments, each loading their own copy from disk. So: **2 copies** (one per session), not 4.

The correct architectural pattern is to load the large object **once, outside `server()`**, in the global scope of `app.R` (or a `global.R` file for Shiny apps). Objects in the global environment are shared across all sessions via R's copy-on-modify semantics — sessions share the reference until they modify it:

```r
# global.R — loaded once, shared across all sessions
big_matrix <- readRDS("500mb_count_matrix.rds")

server <- function(input, output, session) {
  # big_matrix is accessed by reference — no per-session copy
  output$plot1 <- renderPlot({ plot(big_matrix[,1]) })
}
```

For truly large objects that must be writable per-session, `bigmemory` or `arrow` memory-mapped files provide shared memory across processes.

---

**A2.**
S4 uses **inheritance-based method dispatch** with a concept called **method inheritance** (formally: the class graph is searched breadth-first from the most specific class to the most general). When `QCReport` is called on a `PairedSample` object, S4 first looks for a method defined directly on `PairedSample` — finds one, and uses it. When called on a `TumorSample` object, it uses the `TumorSample` method. This is standard **single dispatch with inheritance**.

If an object of a class that inherits from neither is passed, S4 looks up the class hierarchy for that class. If no applicable method is found anywhere in the hierarchy (including the universal `"ANY"` method), S4 throws an error: `unable to find an inherited method for function 'QCReport' for signature '"UnrelatedClass"'`. The correct defensive pattern is to define a fallback:

```r
setMethod("QCReport", "ANY", function(x) {
  stop("QCReport not implemented for class: ", class(x))
})
```

The specific term for the S4 dispatch mechanism is **formal method dispatch** or **S4 dispatch**, distinguished from S3's informal dispatch by requiring explicit `setMethod()` registration and supporting multiple dispatch (dispatching on more than one argument simultaneously) via `signature()`.

---

**A3.**
With n=2 per condition, DESeq2 has only one degree of freedom for estimating within-group variance per gene — which is statistically insufficient. The **dispersion estimate** (which captures gene-wise biological and technical variability) requires replication to distinguish true signal from noise; with 2 samples, the sample variance is an extremely noisy estimator of the true dispersion parameter.

DESeq2 compensates via **empirical Bayes shrinkage of dispersion estimates** (implemented through the `estimateDispersions()` step). It fits a gene-wise dispersion to a **mean-dispersion trend** across all genes, then shrinks each gene's raw dispersion estimate toward this trend. Genes with unstable estimates (due to small n) are pulled strongly toward the global trend. This is called **dispersion shrinkage** or the "sharing of information across genes."

The DESeq2 documentation recommends **at least 3 biological replicates per condition** as a practical minimum, with 6+ per condition for reliable results in differential expression. With n=2, DESeq2 will complete but explicitly warns that results should be treated with caution and fold change estimates will be unreliable for many genes.

---

**A4.**
`dipLogR` is the **log2 copy number ratio at which the diploid (2-copy) segments cluster**. In a pure diploid tumor, this would be 0. But in a tumor that has undergone whole-genome duplication (WGD) or has high aneuploidy, the "modal" copy number state is no longer 2 — it might be 3 or 4. The FACETS algorithm estimates `dipLogR` as the center of mass of the largest cluster of segments weighted by length, which represents the most common copy number state.

`dipLogR` is more useful than log2ratio = 0 because it anchors interpretation: a segment at `dipLogR + 1` is one copy above diploid *regardless* of what the absolute log2ratio is.

A **strongly negative `dipLogR`** (e.g. -0.5 to -1.0) indicates **low tumor purity** — the bulk of the genome appears to be losing signal relative to the normal reference, meaning the tumor fraction is low and the normal cell contamination is pulling the ratio down. Alternatively, a near-universal deletion could cause this, but low purity is the more common interpretation. The scenario where `dipLogR` deviates most dramatically from 0 is a **low-purity, high-ploidy sample** (e.g. 20% purity, tetraploid tumor) — the reference diploid level is buried under a large fraction of normal diploid cells, causing the apparent "diploid" cluster to sit well below 0.

---

**A5.**
The statistical arguments:

**(a) Mean-variance relationship:** Raw RNA-seq counts follow a negative binomial distribution where variance scales with the mean (and exceeds the mean due to overdispersion). Highly expressed genes have both higher means and much higher variances than lowly expressed genes. This is called **heteroscedasticity** — the variance is not constant across the range of the data.

**(b) Effect on PCA:** PCA finds directions of maximum variance. If counts are used directly, the first principal components will be dominated by the most highly expressed genes (which have the largest absolute variance), not by the genes with the most biologically meaningful variation. The result is that PC1 often captures library size differences, not biology.

**(c) Why `log1p` rather than `log`:** Pure `log(0)` is undefined (`-Inf` in R), and single-cell data is extremely sparse — the majority of entries in a scRNA-seq count matrix are zero (a dropout). `log1p(x)` = log(x + 1) maps 0 to 0 (since log(1) = 0), preserving the sparsity structure and avoiding undefined values. The pseudocount of 1 is a pragmatic choice; alternatives exist (e.g. `log(x + 0.5)`) but `log1p` is the universal convention in `scater`/`scran`.

---

**A6.**
The correct primitive is **`reactive()`** — a reactive conductor (also called a reactive expression):

```r
server <- function(input, output, session) {
  computed <- reactive({
    expensive_computation(input$x)
  })

  output$result  <- renderText({ paste("Result:", computed()) })
  output$summary <- renderText({ paste("Summary:", computed()) })
}
```

Simply assigning `val <- expensive_computation(input$x)` inside `server()` does not work because it executes **once at session startup**, not reactively — it does not establish a reactive dependency on `input$x`, so it never re-runs when `input$x` changes.

`reactive()` works because it: (1) lazily evaluates its expression only when first called by a consumer, (2) **caches the result** and returns the cached value to all subsequent callers within the same reactive flush cycle, and (3) automatically **invalidates** when any reactive dependency inside it changes. The invalidation rule is: the reactive expression is invalidated (its cache cleared) when **any reactive input it read during its last evaluation changes**. After invalidation it does not re-execute immediately — it waits until a consumer calls it again (lazy evaluation).

---

**A7.**
The reproducibility problem: `AnnotationHub` resources are versioned but the **default behaviour is to return the latest version** of a resource. If a package calls `AnnotationHub()` and queries for an `EnsDb` annotation without pinning a specific hub ID (e.g. `ah[["AH89426"]]`), a user running the same code six months later may get a different Ensembl version, producing different gene coordinates and IDs — silently breaking reproducibility.

**Two approaches Bioconductor packages use:**

1. **AnnotationHub-backed annotation packages (ExperimentData packages):** Create a dedicated data package (submitted to Bioconductor as an `ExperimentData` package or `AnnotationData` package) that bundles a specific, versioned snapshot of the annotation. The data is fixed at package version time. Examples: `EnsDb.Hsapiens.v86`. This makes the version explicit and tied to package versioning.

2. **Pinning by Hub ID with local caching:** Pin the specific AnnotationHub record ID (e.g. `ah[["AH89426"]]`) rather than querying by metadata. AnnotationHub caches resources locally after first download; subsequent calls return the cached file, not a fresh download. Combined with a pinned ID, this ensures the same resource is returned every time:
```r
ah <- AnnotationHub()
edb <- ah[["AH89426"]]   # specific Ensembl 104 for hg38 — never changes
```

---

**A8.**
Three reasons pseudobulk is required:

1. **Non-independence (pseudo-replication):** Cells from the same patient are not independent observations — they share genetic background, treatment history, and sample processing effects. Treating 15,000 cells from patient A as 15,000 independent data points is pseudo-replication. The true unit of replication is the **patient**, not the cell. Statistical tests that assume independence will have deflated p-values and massively inflated false discovery rates.

2. **Intra-patient correlation (clustering):** Cells within a patient are correlated — their expression profiles are more similar to each other than to cells from another patient. Ignoring this within-patient correlation structure violates the assumptions of standard linear models and t-tests. The appropriate framework is a **mixed-effects model** or, equivalently, pseudobulk aggregation which collapses cells to the patient level before modelling.

3. **Confounding with cell composition:** The number of cells per patient varies (3,000–15,000 in your example). Patients with more cells contribute more "observations" in a naive single-cell analysis, weighting their contribution to the differential expression estimate disproportionately. This confounds the biological effect with sequencing depth / cell recovery efficiency. Pseudobulk aggregation (summing or averaging per patient) equalises patient contributions regardless of cell count.

The statistical concept unifying all three is **hierarchical/nested data structure** — cells are nested within patients, and ignoring the nesting violates the independent-and-identically-distributed (i.i.d.) assumption of standard tests.

---

**A9.**
`Rcpp::NumericVector` is an R object wrapper — it directly references the underlying R SEXP memory, with no copying. Accessing elements still goes through R's memory management. It supports Rcpp Sugar (vectorised operations written in R-like syntax). `std::vector<double>` is a standard C++ container with its own heap-allocated memory, completely independent of R.

**When to prefer each:**
- `Rcpp::NumericVector` at function **boundaries** (parameters and return values) — avoids copying data in/out of R's memory.
- `std::vector<double>` for **internal computation** where you need cache-friendly sequential access, STL algorithms, or to avoid R's overhead on element access.

**Conversion cost:** Converting between them involves copying all elements — O(n) time and memory. For a 10M-element vector this is non-trivial.

**Rcpp Sugar** provides vectorised R-like operations (`log()`, `sum()`, `ifelse()`, etc.) directly on `Rcpp::NumericVector` without converting to `std::vector`. Example: `return log(sum(exp(x)))` where `x` is a `NumericVector` uses Sugar throughout and never creates an intermediate `std::vector` — the entire expression is evaluated lazily using expression templates, producing a single pass over the data.

---

**A10.**
Four BiocCheck requirements that CRAN does not enforce:

1. **Bioconductor coding style:** BiocCheck enforces specific formatting — for example, lines must be ≤80 characters, functions must use `vapply()` not `sapply()` where type-safe (CRAN has no such requirement), and `1:n` in loop indices must be replaced with `seq_len(n)` or `seq_along(x)` to avoid the `1:0` bug.

2. **`NEWS` file requirement:** A `NEWS` file in the package root documenting changes is required by Bioconductor; CRAN merely recommends it.

3. **Vignette requirement:** Bioconductor requires at least one vignette (built with BiocStyle) that demonstrates a complete workflow. CRAN accepts packages with no vignettes.

4. **Namespace import discipline:** BiocCheck enforces that packages do not use `library()` or `require()` inside package functions (use `::` or `@importFrom` instead). CRAN notes this but does not enforce it as strictly.

**Regarding external tools:** Bioconductor requires that any wrapped external command-line tool (e.g. STAR, samtools, bowtie2) is available as a **Bioconductor software package** (`Rsamtools`, etc.) or via **Bioconductor's `ExternalData` or `SystemRequirements`** field. More practically: if the tool is called via `system()` or `system2()`, BiocCheck will flag it and reviewers will require that the tool either be installable via `BiocManager::install()` or that the package gracefully degrades or provides a clear `SystemRequirements` declaration with installation instructions. CRAN has no such community infrastructure requirement.

---

**A11.**
The mathematical relationship: the observed log2(tumor/normal) ratio for a segment is:

```
observed_log2R = log2( (p × CN_tumor/2 + (1-p) × 1) )
```

where `p` is purity (fraction of tumor cells) and `CN_tumor` is the true tumor copy number. The `(1-p) × 1` term represents the diploid normal cell contamination contributing 2 copies / 2 = 1 to the ratio.

The confounding arises because multiple combinations of (purity, ploidy) can produce the same set of observed log2 ratios. For example, a sample with 50% purity and true copy number 4 (tetraploid) produces the same observed ratio as a different combination at different purity. These mathematically indistinguishable solutions are called **degenerate solutions** or **purity-ploidy confounding** — FACETS and similar tools report the solution space and choose among solutions using the likelihood and the constraint that `dipLogR` must be consistent with a plausible ploidy.

A low-purity high-ploidy sample and a high-purity low-ploidy sample can produce identical log2 ratio profiles because the dilution of the tumor signal (low purity) mirrors the effect of the reference-normalised attenuation from high ploidy. This is why manual review of FACETS output — checking whether the purity and ploidy estimates are clinically plausible — remains essential.

---

**A12.**
The access pattern that causes catastrophic performance with `DelayedArray`/`HDF5Array` is **column-by-column (or row-by-row) random access** on a matrix stored with the opposite chunking orientation.

HDF5 stores data in **chunks** — contiguous blocks on disk. If a 50,000 × 30,000 matrix is stored with row-oriented chunks (each chunk contains a contiguous block of rows), then accessing column 1 requires reading partial data from thousands of separate chunks scattered across the file, each requiring a separate disk seek. This is extremely slow — O(n_chunks) disk reads for a single column.

`DelayedArray`'s **`blockApply()`** function is designed to solve this by decomposing operations into **blocks** that align with the chunk boundaries of the underlying storage:

```r
library(DelayedArray)
blockApply(hdf5_mat, FUN = colSums, grid = colGrid(hdf5_mat))
```

By iterating over blocks that match the storage layout, `blockApply()` ensures each disk chunk is read exactly once, converting O(n_chunks²) random access into O(n_chunks) sequential access. The correct use pattern is to always be aware of how the HDF5 file was written (row-chunked vs. column-chunked) and choose operations that traverse it in the same direction.

---

**A13.**
The fundamental limitation is that `renv` manages **R package versions** but does not manage the **Bioconductor release version** itself. Bioconductor packages are released in matched sets (Bioconductor 3.16, 3.17, etc.) that are tightly coupled to specific R versions, and packages from one Bioconductor release are not always installable under another release's environment. `renv::restore()` will attempt to install `BiocGenerics 0.44.0` from Bioconductor's repositories, but if the user's `BiocManager` points to a different Bioconductor release, that exact version may not be available — only the version from the current release.

The two approaches to address this:

1. **Docker/container-based reproducibility:** Package the entire R + Bioconductor environment into a Docker image using `rocker/bioconductor` as a base image. The container pins the R version, Bioconductor release, and system libraries simultaneously. `renv` inside the container then manages package-level pinning. This is the gold standard for fully reproducible Bioconductor environments.

2. **`BiocManager::install()` with `version` pinning + MRAN/CRAN snapshots:** Use `BiocManager::install(version = "3.16")` to explicitly lock the Bioconductor release, combined with `renv`'s ability to use archived CRAN snapshots (via `renv.lock`'s `Repository` field pointing to a CRAN time-snapshot like Posit Package Manager's dated snapshots). This is imperfect — it works for CRAN packages but Bioconductor archives are less consistently accessible.

---

**A14.**
The native `|>` pipe desugars directly into a **function call**: `x |> f(y)` becomes `f(x, y)`. It requires the right-hand side to be a **function call expression** — something that looks like `f(...)`. 

`x %>% {. * 2}` works in magrittr because `%>%` has special handling for `{...}` blocks, treating them as anonymous function bodies with `.` bound to the left-hand side. This is magrittr-specific syntax, not standard R.

`x |> {. * 2}` does not work because `{. * 2}` is not a function call — it is a braced expression, and `|>` does not support that syntax.

The underlying reason: `|>` is implemented as a syntactic transformation at the parser level in R, not as a function. The parser only knows how to desugar `lhs |> fun(args...)` patterns, not arbitrary expressions.

The `|>`-compatible way to apply an anonymous transformation is to use the `\(x)` lambda syntax introduced in R 4.1:

```r
x |> (\(.) . * 2)()   # wrap in lambda, immediately invoke
```

Or more readably, assign the transformation to a name:

```r
double_it <- \(x) x * 2
x |> double_it()
```

---

**A15.**
`intersectColumns()` returns an MAE subset to **only the patients (samples) that have data in ALL experiments simultaneously**. With 180 WES patients, 160 RNA-seq patients, and 140 with both: `intersectColumns()` returns **140 patients**.

The function for "at least one data type" is **`union()`** applied to the column names, but MAE does not have a single `unionColumns()` function — you would subset manually or use `colnames()` across experiments. The MAE itself already holds all 200 patients (the union) — `intersectColumns()` is a filtering operation, not a storage model.

Missing data representation: for patient 41 who has WES but not RNA-seq, the MAE stores `NA` in the **sample map** (`sampleMap`) for the RNA-seq experiment — the patient exists in `colData()` (the master patient table), but has no corresponding entry in the RNA-seq `SummarizedExperiment`. Accessing `experiments(mae)[["RNAseq"]][, "patient_41"]` would throw an error or return nothing — the patient is simply absent from that assay's column space. This is the key design principle of MAE: the master `colData` is the complete patient roster; individual assays are sparse subsets of it.

---

**A16.**
Standard somatic variant callers (Mutect2, Strelka2, VarScan in tumor-normal mode) are designed to detect variants present in a **substantial fraction of tumor cells** — typically VAF > 5–10% in a sample with known tumor purity. They apply filters aggressively tuned to eliminate germline variants and sequencing artefacts at low VAF.

CH variants typically fall in the **VAF range of 0.5%–10%** in peripheral blood (the clone has expanded but still represents a small fraction of hematopoietic cells). This range is problematic because:

1. **It overlaps with sequencing error rates** at standard sequencing depths (∼100–200× for MSK-IMPACT). At 0.5% VAF, a variant needs careful statistical modelling to distinguish from systematic sequencing errors at the same base position.

2. **Standard tumor-normal callers treat the "normal" as truly normal** — a blood sample used as the matched normal for tumor analysis is precisely where CH variants live. A true CH variant in the matched blood normal will be filtered out by the somatic caller as a "germline" or "normal" variant, making it invisible.

3. **Strand bias and position artefacts** accumulate at low VAF thresholds, requiring specialised error-suppression (duplex sequencing, UMIs) or statistical models that account for the specific error profile of the assay.

The signal detection theory concept governing the caller tuning tradeoff is the **receiver operating characteristic (ROC) curve** and specifically the **sensitivity/specificity tradeoff** — equivalently, the **false positive rate vs. false negative rate** tradeoff governed by the classification threshold (minimum VAF, minimum read support, etc.). Lowering the VAF threshold increases sensitivity for true CH variants but decreases specificity (more artefacts called). The operating point must be chosen based on the cost ratio of false positives to false negatives for the clinical application.

---

**A17.**
Two strategies:

1. **Pre-compiled binary approach (vendor/bundle):** The `rextendr` and associated infrastructure support generating **pre-compiled static library artifacts** (`libr_*.a` files) that are bundled with the package source. Rather than compiling Rust from scratch on the CRAN build server, the package ships a pre-compiled static library for each target platform (Linux x86_64, macOS ARM, macOS x86_64, Windows x86_64). The `configure` script detects the platform and links the appropriate pre-built binary. This is the approach used by packages like `hellorust` and emerging bioinformatics Rust packages — it completely eliminates Cargo resolution time at the cost of requiring platform-specific binary builds during the package release process.

2. **Cargo configuration to minimise dependency resolution:** In `Cargo.toml`, use `[profile.release]` with `opt-level = 1` or `2` (not `3`) to reduce compile time at minor performance cost. More importantly, minimise the dependency graph — prefer `no_std` crates where possible and avoid transitive dependencies with large compile footprints. The `Cargo.lock` file committed to the package source ensures Cargo does not re-resolve dependencies on each build — it simply fetches the already-resolved dependency versions from the lockfile. Combined with CRAN's build server having pre-populated Cargo caches for common crates, this can reduce build times from 5+ minutes to under 2 minutes for small packages.

---

**A18.**
Three problems with a single global database connection in multi-user Shiny:

1. **Thread safety / concurrent access:** A single PostgreSQL connection is not safe for concurrent use. If two sessions simultaneously execute queries, requests interleave, causing garbled queries, partial result sets, or connection errors. PostgreSQL connections are stateful — mid-query state (transaction, cursor position) is corrupted by interleaving.

2. **Connection drops and no reconnection logic:** A single global connection established at app startup will eventually time out (PostgreSQL's `idle_in_transaction_session_timeout`, TCP keepalive failures). With a single connection there is no pool to fall back to — the entire app loses database access until the server is restarted.

3. **No per-session query isolation:** A global connection means all sessions share the same database transaction context. One session's uncommitted transaction blocks reads in another session under certain isolation levels, causing sessions to lock each other.

The correct pattern uses **`pool`** (the `pool` R package by Posit):

```r
library(pool)
library(DBI)

# Create once at app startup (global.R)
db_pool <- dbPool(
  drv      = RPostgres::Postgres(),
  dbname   = "genomics",
  host     = Sys.getenv("DB_HOST"),
  maxSize  = 10    # max 10 simultaneous connections
)

onStop(function() poolClose(db_pool))

# In server — pool checks out a connection per query, returns it automatically
server <- function(input, output, session) {
  output$table <- renderDT({
    dbGetQuery(db_pool, "SELECT * FROM variants WHERE gene = $1",
               params = list(input$gene))
  })
}
```

`pool` provides: automatic connection checkout/return, health checking with reconnection on stale connections, configurable pool size limits, and per-query isolation. Each query atomically checks out a connection, executes, and returns it — concurrent sessions never share a connection.

---

**A19.**
A dynamic S7 property that computes from other stored properties — for example, LOH (loss of heterozygosity) status computed from total copy number and lesser copy number:

```r
library(S7)

CopyNumberSegment <- new_class("CopyNumberSegment",
  properties = list(
    chrom = class_character,
    tcn   = class_double,    # total copy number — stored
    lcn   = class_double,    # lesser copy number — stored

    # Dynamic property — computed, not stored
    loh = new_property(
      class = class_logical,
      getter = function(self) self@lcn == 0
    ),

    # Another dynamic property
    log2_ratio = new_property(
      class = class_double,
      getter = function(self) log2(self@tcn / 2)
    )
  )
)

seg <- CopyNumberSegment(chrom = "17", tcn = 3.0, lcn = 0.0)
seg@loh        # TRUE — computed, not stored
seg@log2_ratio # 0.585 — computed, not stored
```

In **S4**, achieving the same thing requires: (1) defining a generic with `setGeneric("loh", function(x) standardGeneric("loh"))`, (2) defining a method with `setMethod("loh", "CopyNumberSegment", function(x) x@lcn == 0)`, and (3) the property is only accessible via `loh(seg)` as a function call — not `seg@loh` (which in S4 would try to access a slot named "loh" and fail if no such slot exists). S7's `getter` syntax unifies property access under `@` regardless of whether the value is stored or computed, which is more ergonomic and mirrors how Python `@property` decorators work.

---

**A20.**
**(a) Why Bootstrap collapses columns:** Bootstrap uses a **12-column grid system** with responsive breakpoints (xs, sm, md, lg, xl, xxl corresponding to screen widths from <576px to ≥1400px). By default, `layout_columns()` with `col_widths = c(4, 4, 4)` translates to Bootstrap's `col-md-4` classes — meaning: 3 columns on medium screens and above, but 1 column (full width) on small/extra-small screens. This is intentional responsive design for mobile usability.

**(b) Two approaches in bslib/Shiny:**

First, use `col_widths` with the `breakpoints` argument to specify behaviour at every breakpoint explicitly, preventing collapse at small screens:

```r
layout_columns(
  col_widths = breakpoints(sm = c(4, 4, 4), xs = c(4, 4, 4)),
  card1, card2, card3
)
```

Second, wrap the specific component in `div()` with an inline style or a Bootstrap class that forces horizontal scrolling rather than column stacking — this is the pragmatic approach for a wide data table:

```r
div(style = "overflow-x: auto;",
  DT::DTOutput("variant_table")
)
```

Combined with `options = list(scrollX = TRUE)` in `DT::datatable()`, this keeps the table columns fixed while allowing horizontal scrolling within the card on narrow screens.

**(c) Bootstrap 5 CSS class to lock column width on all screen sizes:**

Use the **`col-`** prefix *without* a breakpoint infix — but Bootstrap 5 has no true "all sizes" prefix. The closest is **`col-xs-*`** (Bootstrap 3 legacy) which doesn't exist in Bootstrap 5. The correct Bootstrap 5 approach is to use `col-` (no breakpoint), which applies to all screen sizes:

However, more practically, use `col-sm-4` at the smallest breakpoint you care about — `sm` (≥576px) means the 3-column layout applies from 576px upward and never collapses below that. Alternatively, within bslib, pass `fill = FALSE` and fixed pixel widths to `layout_columns()` to opt out of Bootstrap's responsive grid entirely and use a fixed CSS layout.

---