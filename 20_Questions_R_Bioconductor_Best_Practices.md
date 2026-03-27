# 20 Questions: R, Bioconductor & Shiny — Best Practices & What Changed
### Based on: github.com/swajid/Major-and-minor-changes-in-R-and-Bioconductor-since-2017

> Questions focus on **current best practices** and **how they differ from older approaches**.  
> Answers are at the bottom. Try to answer each before looking.

---

## Questions

**Q1.**
You have a data frame loading function written in 2015 that never used `stringsAsFactors = FALSE`. You run it under R 4.0+ and your downstream code silently breaks because it expected factor columns. What changed in R 4.0.0 that caused this, and what is the current best practice?

---

**Q2.**
A colleague hands you this 2016-era code:

```r
filter_samples <- function(df, col) {
  dplyr::filter(df, col > 0.3)
}
```

It doesn't work correctly when called as `filter_samples(samples, purity)`. What is the problem and what is the modern fix?

---

**Q3.**
You need to write a fast GC-content calculator in a new R package. You are choosing between Rcpp, cpp11, and Rust via extendr. The function is string-heavy — iterating over every base in millions of sequences. Which tool does the README recommend for this use case and why? What specific Rust feature makes safe parallelism easy?

---

**Q4.**
In 2017 your Bioconductor package accepted `ExpressionSet` objects. A reviewer in 2025 says your package is using a legacy container. What should you migrate to, and what are the three accessor functions that replace `exprs()`, `pData()`, and `fData()`?

---

**Q5.**
You are building a Shiny app for a lab that will run FACETS on user-uploaded BAM files — a process that takes 2–5 minutes. In early Shiny this would freeze the entire UI. What is the current (Shiny 1.8+) best practice for handling this, and what bslib component pairs with it to give the user visual feedback?

---

**Q6.**
Describe the full modern R package development workflow from scratch. You have an idea for a package called `cnvtools`. List, in order, the `usethis` and related commands you would run to: create the package, set up testing, set up CI, set up a package website, and initialise git.

---

**Q7.**
You are writing documentation for an exported S4 method in a Bioconductor package in 2026. What documentation system should you use, and what is one specific syntax difference between the old Rd-macro style and the current recommended style for inline code formatting?

---

**Q8.**
A postdoc wants to share their Shiny variant annotation app publicly via a URL — no server budget, no Posit Connect. The data is already public (no PHI). What technology enables this, what does it use under the hood, and what is one hard constraint they need to know before committing to it?

---

**Q9.**
What is the difference between `S4` and `S7` in terms of: (a) where they live, (b) backward compatibility with S3, and (c) the one thing S7 cannot do that S4 can regarding inheritance?

---

**Q10.**
You wrote a dplyr pipeline in 2016 using `%>%` extensively. A junior collaborator asks whether they should use `%>%` or `|>` in new code. What are two concrete behavioural differences between them that matter for deciding, based on what the README documents?

---

**Q11.**
In the DNA → Amino Acid translator Shiny example, `updateOn = "blur"` is used on the sequence text input. What is the default behaviour without this option, why is it problematic for sequence inputs specifically, and what is the alternative trigger for `updateOn`?

---

**Q12.**
The README describes Shiny modules as "mandatory practice for any non-trivial app." What specific problem do modules solve, what is the key function used inside the UI half to create a namespaced ID, and what is the key function used inside the server half to wire up the module?

---

**Q13.**
You have a 500,000-cell single-cell RNA-seq count matrix that is too large to load into RAM. What Bioconductor package and class does the README recommend for this, and what is the function call pattern for loading it?

---

**Q14.**
Compare the old CI/CD standard for R packages (2014–2017 era) with the current standard. Name the old system, why it is no longer used, and what replaced it. What `usethis` command sets up the new system automatically?

---

**Q15.**
You have a DESeq2 analysis from 2017 that calls `results(dds)` to get fold changes. A reviewer flags this as outdated. What is the current recommended approach for obtaining fold change estimates from DESeq2, and what shrinkage method is currently the default?

---

**Q16.**
The README distinguishes four OOP systems in active use: S3, S4, R6, and S7. For each of the following use cases, state which system is recommended:
- (a) A new Bioconductor package extending `SummarizedExperiment`
- (b) A stateful pipeline runner object that logs steps and modifies itself
- (c) A new CRAN package needing formal OOP without Bioconductor dependency
- (d) Simple typed results with `print()` and `summary()` methods

---

**Q17.**
In the 2014–2017 era the dependency management tool was `packrat`. What replaced it and when, and what are the three core commands of the modern tool that correspond to: initialise, record current state, and reproduce on another machine?

---

**Q18.**
`bslib` introduces `page_sidebar()`, `layout_columns()`, `card()`, and `value_box()`. Map each of these to what you would have used in the `shinydashboard` / base Shiny era, and identify which `shinydashboard` function is now considered to be in maintenance mode.

---

**Q19.**
The README says Rcpp now "emits a warning on out-of-bounds vector access." What was the old behaviour, why does this matter for genomics code specifically, and what does the README say this warning will eventually become?

---

**Q20.**
You are deploying a Shiny app that processes clinical genomic data (tumor purity estimates, variant calls). A colleague suggests using shinylive/GitHub Pages for free hosting. Based on the README's deployment guidance, give two specific reasons why this is the wrong choice for this use case, and name two deployment options that would be appropriate instead.

---
---
---

## Answers

**A1.**
In R 4.0.0, the default for `stringsAsFactors` changed from `TRUE` to `FALSE`. Before 4.0, `data.frame()` and `read.table()`/`read.csv()` silently converted character columns to factors. Any code that relied on this implicit conversion — such as downstream subsetting by integer level, or passing a factor to a function expecting one — would silently behave differently or produce wrong results. The current best practice is simply to remove all defensive `stringsAsFactors = FALSE` arguments from new code (they are no longer needed), but to audit old pipelines for any code that implicitly depended on columns being factors.

---

**A2.**
The problem is non-standard evaluation (NSE). In `dplyr` 0.7+, the old NSE approach broke when `col` was passed as a bare name from an outer function — it was evaluated in the wrong environment. The modern fix uses `tidyeval`: capture the column argument with `enquo()` and inject it with `!!`:

```r
filter_samples <- function(df, col) {
  col <- enquo(col)
  dplyr::filter(df, !!col > 0.3)
}
```

This pattern was introduced with `rlang` in dplyr 0.7 (2017) and is required for any function that programmatically wraps dplyr verbs.

---

**A3.**
The README recommends **Rust via extendr** for string-heavy parsing workloads (BAM/VCF-like operations). The reasoning is that Rust provides compile-time memory safety guarantees — no buffer overflows or undefined behaviour — that C++ cannot guarantee. The specific Rust feature that makes safe parallelism trivial is the **`rayon` crate**, which provides data-parallel iterators; adding `.par_iter()` instead of `.iter()` parallelises loops with no data races because the borrow checker enforces thread safety at compile time.

---

**A4.**
Migrate to `SummarizedExperiment` (or a subclass appropriate to your data type, such as `RangedSummarizedExperiment` if you have genomic coordinates on rows). The replacement accessors are:

- `exprs(eset)` → `assay(se, "counts")` (or the named assay of your choice)
- `pData(eset)` → `colData(se)`
- `fData(eset)` → `rowData(se)`

The README notes that most new Bioconductor packages do not accept `ExpressionSet` at all as of the 2018–2020 transition period.

---

**A5.**
The current best practice is `ExtendedTask` (introduced in Shiny 1.8+, 2024). An `ExtendedTask` wraps a `future()` call and runs the slow computation in a background R process, leaving the reactive graph — and therefore the UI — fully responsive. The bslib component that pairs with it is `input_task_button()`, which replaces a plain `actionButton()` and automatically: (a) shows a spinner while the task runs, and (b) prevents the user from clicking again until the task completes, avoiding duplicate job submissions.

---

**A6.**
The modern sequence of commands:

```r
usethis::create_package("cnvtools")   # scaffold package structure
usethis::use_testthat()               # set up testthat edition 3
usethis::use_github_actions()         # add R-CMD-check.yaml for CI
usethis::use_pkgdown()                # configure package website
usethis::use_git()                    # initialise local git repo
usethis::use_github()                 # create GitHub remote and push
```

Additionally, `usethis::use_mit_license()` (or another licence) and `usethis::use_vignette("intro")` would be standard. The key point is that in 2017 all of this was done manually or via `devtools` functions; `usethis` as a standalone package did not exist.

---

**A7.**
Use **`roxygen2`** with Markdown support enabled (standard since ~2019; declared via `Roxygen: list(markdown = TRUE)` in `DESCRIPTION`). The key syntax difference: the old Rd-macro style used `\code{myfunction}` for inline code formatting; the modern Markdown style uses backticks `` `myfunction` ``. Additionally, `roxygen2` now uses `@returns` (not `@return`) as the preferred tag, and full Markdown — headings, bold, links, code blocks — is supported natively in all documentation fields.

---

**A8.**
The technology is **shinylive**, which runs Shiny entirely in the browser using **WebAssembly via the webR project** — no R server is required at all. The app is exported as static files and hosted on GitHub Pages or Netlify. The one hard constraint they must know: **all source code and data are visible to the client** — there is no server-side privacy. This makes it unsuitable for PHI, credentials, or any proprietary data. Additionally, not all R packages are available as pre-compiled WebAssembly binaries, so package availability must be verified at `repo.r-wasm.org`.

---

**A9.**
(a) **Location:** S4 lives in base R (the `methods` package); S7 is currently a CRAN package (intended for eventual base R inclusion but not there yet as of 2026).

(b) **S3 backward compatibility:** S7 objects *are* S3 objects — S7 is built on top of S3, so existing S3 code works with S7 objects. S4 objects are not S3 objects; they use a separate dispatch system.

(c) **Inheritance limitation:** S7 classes **cannot extend existing S4 classes**. This means if you want to inherit from `SummarizedExperiment`, `GRanges`, or any other Bioconductor core class (all of which are S4), you must use S4. S7 → S4 inheritance is explicitly not supported.

---

**A10.**
Two concrete differences:

1. **Placeholder syntax:** The magrittr `%>%` pipe uses `.` as a placeholder (e.g. `x %>% lm(y ~ ., data = .)`), while the native `|>` pipe uses `_` as a placeholder (R 4.2+), but with restrictions — `_` must appear on the right-hand side only in a named argument position (e.g. `mtcars |> lm(mpg ~ cyl, data = _)`).

2. **No placeholder by default:** `|>` has no placeholder at all before R 4.2, so many pipeline patterns that relied on `.` simply do not work with `|>` and need rewriting. The README notes that some pipeline patterns require rewriting when switching from `%>%` to `|>`.

---

**A11.**
The default behaviour without `updateOn = "blur"` is `updateOn = "change"`: the reactive fires immediately on every character the user types. For sequence inputs this is problematic because it triggers `validate_sequence()` and `translate_dna()` on every single keystroke — so typing a 50-base sequence would trigger 50 reactive evaluations, producing errors for invalid partial sequences and potentially making the UI sluggish. With `updateOn = "blur"`, the reactive fires only when the text input **loses focus** (user clicks away) or when the user **presses Enter** (or Cmd/Ctrl+Enter for `textAreaInput`).

---

**A12.**
Modules solve the **input/output ID namespace collision problem** in large Shiny apps. Without modules, every input and output ID in the app must be globally unique — in large apps this becomes impossible to manage, and copy-pasting UI components creates silent ID collisions. 

The key function in the UI half is **`NS(id)`**, which returns a namespace function; every ID inside the UI is wrapped with `ns("myid")` to produce namespaced IDs like `"mymod-myid"`.

The key function in the server half is **`moduleServer(id, function(input, output, session) { ... })`**, which wraps the server logic and ensures `input$`, `output$`, and `session` are scoped to the module's namespace automatically.

---

**A13.**
The README recommends the **`HDF5Array`** package, which provides the `HDF5Array` class (extending `DelayedArray`). The loading pattern is:

```r
library(HDF5Array)
counts <- HDF5Array("counts.h5", "matrix")
```

This creates an object that behaves like an in-memory matrix but reads data from disk on demand. Operations on it are lazy — they are only evaluated when results are actually needed — which allows standard Bioconductor workflows to operate on data far larger than available RAM.

---

**A14.**
The old system was **Travis CI** (`.travis.yml`). It is no longer the standard primarily because GitHub acquired its own CI/CD platform, **GitHub Actions**, which launched in 2018 and became dominant by ~2020. Travis CI also changed its pricing model, removing the free tier for open-source projects. The replacement is **`r-lib/actions`** on GitHub Actions, using a standard `R-CMD-check.yaml` workflow that tests across Ubuntu, macOS, and Windows. The `usethis` command to set it up automatically is:

```r
usethis::use_github_actions()
```

---

**A15.**
The current recommendation is to use **`lfcShrink()`** rather than calling `results()` directly for fold change estimates. The default shrinkage method is **`apeglm`** (adaptive t prior shrinkage), which was not available in the 2017-era DESeq2:

```r
res_shrunk <- lfcShrink(dds,
                        coef = "condition_tumor_vs_normal",
                        type = "apeglm")
```

The rationale is that raw log2 fold changes from `results()` are noisy for low-count genes; `lfcShrink()` pulls extreme estimates toward zero, producing more reliable rankings and better volcano plots. The older `normal` shrinkage method (the 2017 default within `results()`) is now considered inferior.

---

**A16.**
- **(a) Extending SummarizedExperiment → S4.** Required for Bioconductor core; `SummarizedExperiment` itself is S4, and inheritance must use `setClass(contains = "SummarizedExperiment")`.
- **(b) Stateful pipeline runner → R6.** Reference semantics (modify-in-place) make R6 natural for objects that accumulate state across method calls. `invisible(self)` enables method chaining.
- **(c) New CRAN package with formal OOP → S7.** The README explicitly lists S7 as the "new recommendation" for CRAN packages needing formal OOP. It is forward-looking and backward compatible with S3.
- **(d) Simple typed results with print/summary → S3.** The README recommends S3 as the default for "simple typed results with print/summary methods" — it is idiomatic R and the lowest-friction option.

---

**A17.**
`packrat` was replaced by **`renv`**, released in **2019**. The three core commands:

```r
renv::init()        # initialise — creates renv.lock and activates the project library
renv::snapshot()    # record — writes current package versions to renv.lock
renv::restore()     # reproduce — installs packages from renv.lock on another machine
```

The README notes that `packrat` was "functional but slow and not widely adopted" — most teams in the 2014–2017 era relied on `sessionInfo()` output for reproducibility instead, which is not machine-actionable.

---

**A18.**
Mappings:

| bslib (current) | Old equivalent |
|----------------|---------------|
| `page_sidebar()` | `sidebarLayout()` in base Shiny, or `dashboardPage()` in shinydashboard |
| `layout_columns()` | `column()` + `fluidRow()` in base Shiny |
| `card()` | `box()` in shinydashboard |
| `value_box()` | `infoBox()` / `valueBox()` in shinydashboard |

The function now considered to be in **maintenance mode** is `dashboardPage()` from `shinydashboard` — the README states new projects should use `bslib` instead, and that `bslib` is "a viable alternative to shinydashboard."

---

**A19.**
The **old behaviour** was silent — reading past the end of an Rcpp vector was undefined behaviour in C++, meaning it could return garbage data, corrupt memory, or crash with no informative error message. For genomics code this matters because off-by-one errors are common when iterating over genomic windows, codon positions, or index arrays, and silent corruption of numerical results is worse than a visible error. The README states that the warning **"will become an error in a future Rcpp release"** — meaning old code that accidentally relied on out-of-bounds access will eventually fail loudly rather than silently, which is the safer and intended behaviour.

---

**A20.**
Two reasons shinylive/GitHub Pages is wrong for clinical genomic data:

1. **No data privacy:** The README explicitly states "All data and code are visible to the client — not suitable for PHI or proprietary data." Tumor purity estimates, variant calls, and patient-linked genomic data are PHI and cannot be sent to the client browser.

2. **All computation runs client-side:** Shinylive shifts compute to the user's browser. For genomic analyses, this means sensitive data would need to leave the server and be transmitted to — and processed on — the user's machine, which violates the architecture required for PHI handling.

Appropriate alternatives from the README's deployment table:
- **Posit Connect** — enterprise-grade, supports authentication, scheduling, and access control.
- **ShinyProxy** — Docker-based, self-hosted, multi-user; data stays on institutional infrastructure and never leaves the server.

