# current_exercises — manifest

What has actually been created, with provenance. One row per file.
The full build list is [`EXERCISES.md`](EXERCISES.md); the rules are
[`../agentGuide/EXERCISE_RULES.md`](../agentGuide/EXERCISE_RULES.md).

| # | File | Source | Source date | Modified? | Supersedes |
|---|---|---|---|---|---|
| 2 | `linux/intro_bash_linux.md` | `kenya2026/exercises/Day1/IntroToBash.ipynb` | 2026-08-18 | yes — converted ipynb -> md | `kenya2026/exercises/post_course/day1_morning_bash_linux.ipynb` (2026-08-25) |
| 1 | `linux/intro_linux.md` | `summer2025/BriefIntro2Linux.md` | 2025-07-23 | yes — course branding removed | `summer2024/BriefIntro2Linux.md` (2024-07-04) |

| 9 | `shiny/needleman_wunsch_dna.R` | `BSA/NW_DNA.R` | 2025-09-02 | no — copied verbatim | none |
| 10 | `shiny/needleman_wunsch_blosum50.R` | `BSA/needleman_wunsch_shiny_app_blosum_50.r` | 2025-08-30 | yes — self-`source()` URL repointed | none |
| 11 | `shiny/dotplot.R` | `BSA/dotplotShiny.R` | 2025-08-30 | yes — self-`source()` URL repointed | none |
| 3 | `shiny/stats_binomial.R` | `stat_molbio/binom.R` | 2026-01-16 | yes — `library(shiny)` + `shinyApp()` added | none |
| 4 | `shiny/stats_normal.R` | `stat_molbio/normal.R` | 2026-01-16 | yes — `library(shiny)` + `shinyApp()` added | none |
| 5 | `em_algorithms/em_algorithm.ipynb` | `advBinf/exercises/advBinf_EM_algorithm.ipynb` | 2026-09-09 | yes — 3 bugs fixed, figure added, headings restructured, typos, citation | none |
| 6+7 | `em_algorithms/haplotype_frequencies.ipynb` | `advBinf/exercises/solution_haplotype_frequencies.ipynb` | 2025-09-12 | yes — rebuilt as one scaffolded notebook; print bug fixed | `advBinf/exercises/haplotype_frequencies.ipynb` (2025-09-12) |
| 12 | `ngs/ngs_intro_human.ipynb` | `chinacourse2026/Day2_Morning_NGSintro_human.ipynb` | 2026-09-14 | yes — quizzes, questions, figure, typos, paths | 5 older copies (see EXERCISES.md) |
| 13 | `ngs/ngs_intro_animal.ipynb` | `kenya2026/exercises/Day1/Kenya2026_NGSintro.ipynb` | 2026-08-17 | yes — quizzes, questions, typos, paths | 6 older copies (see EXERCISES.md) |
| — | `ngs/quiz/*.json` (10 files) | new + `kenya2024/.../quiz{1..4}.json` | 2026-09-16 | new quiz bank | kenya2024 quiz1-4 |
| 14 | `ngs/ngs_inference_human.ipynb` | `summer2025/exercises/Day2_NGS_Inference.ipynb` | 2025-08-04 | yes — split, quizzes, questions, bug fixes, paths | `summer2024/exercises/NGS_inference.ipynb`, `summer2023/NGSinference/` |
| 19 | `em_algorithms/pca_em_human.ipynb` | `advBinf/exercises/advBinf_PCA_EM.ipynb` | 2026-09-16 | yes — 3 quizzes, 14 question blocks | none |
| 24 | `em_algorithms/admixture_em_human.ipynb` | `advBinf/exercises/advBinf_admixture_EM.ipynb` | 2026-09-14 | yes — 2 quizzes | none |
| 37 | `em_algorithms/sfs_model.ipynb` | `advBinf/exercises/advBinf_SFSmodel.ipynb` | 2026-09-16 | yes — quiz bank moved in, typo | `advBinf/exercises/SFS.md` (2024-09-20) |

## Notes

- **#2 `intro_linux.md`** — checked before copying: contains no quizzes, no data
  paths, and no references under `/course`, `/davidData` or a home directory.
  The only external links are the MobaXterm and XQuartz download pages, which
  were kept.

  Course-specific references removed (R11):
  - dropped the `Summer Course in Population Genetics` / `University of Copenhagen`
    title block
  - `For this course you are not required to know linux` -> `You are not required
    to know linux`
  - `enough to get by for the entire week` -> `enough to get by for the exercises`

  Kept deliberately: the server name `emily.popgen.dk` (shared infrastructure,
  not course-specific) and the phrase "ask when you are here".

- **#3 `stats_binomial.R`, #4 `stats_normal.R`** — these are **Shiny apps**, so
  they live in `shiny/`, not with the statistics material. My first detection
  pass missed them: it grepped for `shinyApp|library(shiny)|runApp`, and these
  files contain none of those — they define `ui <- fluidPage(...)` and
  `server <- function(input, output)` and then simply stop.

  That is also a bug: without `library(shiny)` and a `shinyApp()` call, sourcing
  either file defines two objects and launches nothing. Both were added, plus
  the `source(...)` launch comment the other Shiny apps carry. All five apps in
  `shiny/` now have both.

- **#9, #10, #11 are all Shiny apps** and live in `shiny/`, not `shared/` (R13).
  All three were found by `grep -rliE 'shinyApp|library\(shiny\)|runApp'`; note
  that `BSA/NW_DNA.R` is a Shiny app too despite its plain name.

  Renamed for clarity now that the folder supplies the "shiny" part, and because
  the two Needleman-Wunsch apps needed telling apart:
  - `BSA/NW_DNA.R` -> `shiny/needleman_wunsch_dna.R` (DNA, match/mismatch/gap)
  - `BSA/needleman_wunsch_shiny_app_blosum_50.r` -> `shiny/needleman_wunsch_blosum50.R`
    (protein, BLOSUM50)
  - `BSA/dotplotShiny.R` -> `shiny/dotplot.R`

  None has a fixed data path — the dotplot app takes its FASTA through a
  `fileInput` upload. #10 and #11 each open with a commented `source()` line
  giving their own raw-GitHub URL so students can launch them in one line; those
  URLs were repointed at `current_exercises/shiny/`. Nothing else changed.

  Note: those URLs only resolve once `current_exercises/` is pushed to
  `popgenDK/courses` on `main`.

- **#12 `motif_discovery.R` — DEFERRED by request (2026-09-16), not copied.** The data has now been
  **found**: `/davidData/data/BSA` (1.2 G) is the `COURSES/BIO/BSA` workgroup
  share that gets mounted into home directories as `~/work/COURSES/BIO/BSA`, so
  the script's paths were never a laptop path after all. Both files are present:
  `motif_discovery/PUM2.top500.fa` and `pHMM/globins4.fasta`.

  The data is now **copied** to `data/BSA/` (2026-09-16, 84/84 files verified),
  so the three path references can be rewritten to
  `data/BSA/motif_discovery/PUM2.top500.fa` and `data/BSA/pHMM/globins4.fasta`.

  One blocker remains: `BSA/motif_discovery_ex.R` is an instructor
  working/solution draft rather than a student-facing exercise — it has the
  answers inline as comments and lines 100-102 are not valid R.

- **Helper libraries live in the data folder, not here (2026-09-16).**
  `admixFun.R`, `newPlotPlink.R` and `online.R` are function libraries that
  exercises `source()`, not exercises, so all three sit in `data/scripts/`:

  | File | Source | Sourced by |
  |---|---|---|
  | `data/scripts/admixFun.R` | `chinaCourse2025/assets/admixFun.R` | #25, #26 |
  | `data/scripts/newPlotPlink.R` | `summer2023/InfererPopStructure/newPlotPlink.R` | #27, #46, #47, #53 |
  | `data/scripts/online.R` | `kenya2024/online.R` | #18, #53, #54 |

  All copied verbatim. `current_exercises/` holds only exercises.

- **`data/geneticMap/`** (480 M, 85/85 files) copied from
  `/course/scripts/geneticMap` — a dependency of `newPlotPlink.R` that the
  original data survey missed, because it is referenced from a helper library
  rather than from any notebook.

- **Date correction** — `BSA/needleman_wunsch_shiny_app_blosum_50.r` is
  2025-08-30, not 2025-09-02 as first recorded. The original survey globbed
  `*.R` and missed the lowercase `.r` extension.

- **#2 `intro_bash_linux.md`** — converted from `IntroToBash.ipynb` to terminal
  markdown on 2026-09-16. 67 notebook cells -> 547 lines of markdown with 29
  fenced `bash` blocks.

  The notebook form fought the material: `less`, `nano`, `top` and `man` are all
  interactive and cannot run in a Jupyter cell, so the original worked around
  them. Those workarounds are removed, because in a terminal the commands work:
  - dropped "To run the selected code cell, hit `Shift + Enter`"
  - section 2.1 — dropped "Because we are not in the terminal, I will copy the
    file content here" and its canned text; the reader now runs `less` / `cat`
  - section 2.2 — dropped "that is not possible in a notebook" and the pasted
    nano shortcut list; the reader now opens `nano` and reads its bottom bar
  - section 6.5 — dropped the pasted `top` output; the reader now runs `top`
  - removed two `attachment:` images (`tux.png`, `solution.png`) that were
    embedded in the `.ipynb` and cannot survive as markdown

  Added a short note at the top saying it is run in a terminal. Data path
  repointed from `/course/popgenmsc26/exercises/linux/Exercises.zip` to
  `/course/data/popgenmsc26_exercises/linux/Exercises.zip` (verified present).
  The CC-BY-SA license and Kristian Rother attribution are preserved verbatim.

- **#12 / #13 the two NGS intro notebooks** — built 2026-09-16.

  | | human | animal |
  |---|---|---|
  | cells | 68 -> 91 | 74 -> 101 |
  | quizzes | 4 -> 8 | 4 -> 8 |
  | code cells followed by a question or quiz | 28/31 | 29/31 |

  The three or four code cells without a question are environment setup
  (`ROOT_PATH=`, `setwd`, `os.chdir`) where a question would be noise.

  **Figures.** The pipeline-position figures `ngs_files1-4.png` (FASTQ -> QC ->
  mapping/BAM -> VCF, hosted on popgen.dk) are now in both notebooks. The human
  one was missing `ngs_files2` at the FastQC step; added.

  **Quiz bank** in `ngs/quiz/`, loaded by raw-GitHub URL (the pattern the
  kenya2026 notebook already used). Six are shared between the two notebooks
  (`fastq_format`, `sam_format`, `sam_to_bam`, `depth`, `vcf`, `flags`) and four
  are track-specific because the answers are dataset-dependent
  (`fastq_counts_*`, `fastqc_*`).

  All numeric answers were computed from the copied data, not guessed:
  animal 334024 lines / 83506 reads / 150 bp (matches the old kenya2024 quiz,
  which validates it); human 687520 lines / 171880 reads / 100 bp, median base
  quality 40, with *Per sequence GC content* the only FastQC warning (read out
  of the FastQC zip).

  **Typos fixed** (about 40 occurrences, both notebooks): defermines->determines,
  swich->switch, jypiter->Jupyter, chromome->chromosome, refence/referecne->reference,
  concensus->consensus, extact->extract, likelely->likely,
  heterzygoes/heterozygoes->heterozygous, particually->particularly,
  usefull->useful, easiy->easily, truely->truly, faciliate->facilitate,
  seperator->separator, contrains->contains, postion->position,
  "tell tell you"->"tell you", "are are called"->"are called",
  "Is is possible"->"Is it possible", "How many lines to you have"->"do you have",
  "In this exercise will cover"->"we will cover", a stray leading comma, and
  "where you only the reads where both ends maps"->"where you keep only the reads
  where both ends map".

  **Paths** repointed to the data store: `/davidData/data/course/kenya2026/anders`
  -> `/course/data/kenya2026_anders`, `/course/chinacourse2026/shared` ->
  `/course/data/chinacourse2026_shared`. Both verified present. The human setup
  cell no longer copies `quiz*.json` from the old shared folder.

  **Caveat:** the quiz URLs only resolve once `current_exercises/` is pushed to
  `popgenDK/courses` on `main` (same caveat as the shiny scripts, R12).

- **#6+7 `haplotype_frequencies.ipynb` — two notebooks became one (2026-09-16).**
  The old pair was `haplotype_frequencies.ipynb` (the exercise) and
  `solution_haplotype_frequencies.ipynb` (the answers). The exercise half was five
  empty code cells under the instruction "make your own EM algorithm", which was too
  hard to work through, so the new notebook is **rebuilt from the solution** and there
  is no separate solution notebook.

  It keeps the solution's structure, variable names, comments and all of its LaTeX
  derivations, and blanks only the lines that carry the idea, marked `????`:
  17 blanks across 9 numbered tasks, 1-2 lines each, with the surrounding loops and
  bookkeeping given. `countHaplotypes()` and `calculate9Q()` are given in full — they
  are index bookkeeping, not the learning point. The tasks build up in the solution's
  own order: allele frequencies -> haplotype frequencies with the haplotypes visible
  -> `likeG` -> log likelihood -> fast log likelihood -> E step -> M step + EM run ->
  fast EM.

  Made easier than the old exercise by:
  - splitting one "write an EM algorithm" instruction into 9 tasks
  - giving every function signature, loop and comment, so only the key expression is missing
  - stating above each task which formula in the text it implements
  - adding self-checks the student can use without an answer sheet: the posterior must
    sum to 1, the fast likelihood must equal the slow one, the fast EM must iterate
    identically to the slow one, and the EM must converge to the Task 3 estimate
  - printing the true haplotype frequencies at the end of both EM runs to compare against

  Changed from the source solution:
  - **bug fixed:** the "haplotype frequencies" cell computed `estFreqHap` and then
    printed `round(freqHap,4)`, the *true* frequencies, so the estimate it had just
    calculated was never shown and appeared to be exact. It now prints both, labelled.
  - dropped the Colab link to the solution that opened the old exercise (R11)
  - `#Estimating` -> `# Estimating` in the title, which was not rendering as a heading
  - dropped a trailing empty code cell
  - `emStepFast(theta,data)` -> `emStepFast(theta,dataTab)`, since the body reads
    `dataTab`; the old signature only worked because `dataTab` was a global
  - removed the unused `N <- nrow(dataTab)` from the fast likelihood (`nrow` of a
    3x3 table is 3, not the sample size, so it read as a bug)

  No data paths and no quizzes, so R4 and R10 need nothing. Kernel left as `ir`
  (a plain R notebook, not SoS). Committed **without outputs**, per the advBinf
  convention in `RUN_DATA_NOTEBOOKS.md`.

  Verified by filling in all 17 blanks in a scratchpad copy and running it through
  `Rscript`: it reproduces the source solution's numbers exactly — log likelihood
  -1670.391 at convergence and theta 0.0300/0.2010/0.5655/0.2035 against a true
  0.03/0.20/0.55/0.22. No `.html` render is shipped; say the word if you want one
  rendered from a filled-in copy as the published solution, the way
  `advBinf/exercises/*.html` works.

- **NGS intro repointed at a clean data folder (2026-09-16).** Exercises should
  not read from a folder named after the course they came from, so the data was
  curated into `data/NGSintro/{animal,human,software}` and both notebooks now
  point there. The working directories lost their course names too:

  | | before | after |
  |---|---|---|
  | animal data | `data/kenya2026_anders/NGSintro_day1` | `data/NGSintro/animal` |
  | animal work dir | `~/kenya2026/NGSintro` | `~/ngs_intro_animal` |
  | human data | `data/chinacourse2026_shared/data/NGSIntro` | `data/NGSintro/human` |
  | human work dir | `~/sysu2026_day2_ngsintro` | `~/ngs_intro_human` |

  The human setup cell's four-way `ROOT_PATH`/`TOOL_PATH`/`SHARED_PATH`/
  `INPUT_PATH` split collapsed to `DATA` + `SOFTWARE`. `picard.jar` was
  byte-identical in both course folders, so one copy is shared. `FASTQC=fastqc`
  is now a variable rather than a hardcoded tool path.

  Two FastQC screenshots still loaded from `kenya2024/` and `chinaCourse2025/`.
  They were copied to `ngs/fastqc_report_animal.png` and
  `ngs/fastqc_report_human.png` so the exercises are self-contained.

  **Fixed along the way:** `chr21.fa.gz` had no `.fai`/`.gzi` index in the
  original course folder, although `samtools tview` and `bcftools mpileup -f`
  both need one. Generated in the clean folder; all 7 index files are now
  present for both references.

- **#5 `em_algorithm.ipynb`** — copied from the 2026-09-09 advBinf notebook and
  corrected. Unlike #6+7 this one was already a single fully worked notebook with
  no blanks and no separate solution, so the structure and the pedagogy are
  unchanged: 23 cells, same order, same code. Kernel left as `ir`, committed
  without outputs. No data paths and no quizzes, so R4 and R10 need nothing.

  **Three real bugs fixed:**
  - the two-coin EM initialised `theta <- c(0.5,0.6)` while the text above it says
    to start at $\theta_A=0.6, \theta_B=0.5$. Coin A therefore converged to 0.52 and
    coin B to 0.80 — the labels swapped relative to the text, relative to Figure 1
    of the source article the exercise tells students to compare against, and
    relative to the notebook's own Bonus 2 cell, which does start `(0.6,0.5,0.9)`.
    Now `c(0.6,0.5)`, giving $\theta_A=0.797$, $\theta_B=0.520$.
  - the GATK genotype-likelihood formula wrote $g_{TT}$ with a plain `e` instead of
    `\epsilon`, inconsistent with the $g_{CC}$ and $g_{CT}$ lines above it. The R code
    was always right; only the formula was wrong.
  - the E-step Bayes formula used curly typographic apostrophes in $z'$, which
    MathJax renders literally instead of as a prime.

  **Rendering fixes** (per the math pitfalls in `RUN_DATA_NOTEBOOKS.md`):
  - three multi-line formulas were inline `$...$`, one of them containing a `\\`
    line break, which MathJax cannot do inline. Made them display `$$...$$`.
  - a literal Unicode `∝` inside math replaced with `\propto`.

  **Structure:** the notebook was titled "Two-Coin EM Example" but contains three
  parts, with a second `#` heading at "EM algorithm examples" and a third at
  "Coin toss - how to find the maximum likelihood for a binomial". Retitled
  `# The EM algorithm` with a 3-item contents list, and the two later `#` headings
  demoted to `# 1.`/`# 2.`/`# 3.` sections. **Cell order is unchanged** — note that
  part 3 is the gentlest material (plain binomial ML, no latent variable, no EM) and
  still sits last, so the contents list says to read it first if likelihoods are new.
  Say the word if you want it moved to the front instead.

  **Other changes:**
  - the coin example is now cited — Do CB & Batzoglou S (2008), *What is the
    expectation maximization algorithm?*, Nature Biotechnology 26:897-899
    (**confirmed correct, 2026-09-16**). The exercise asked students to "compare with
    the review" and to read "the figure in the article" without ever naming either.
  - **a figure was added** (2026-09-16), two cells after the two-coin EM run: our own
    version of that paper's Figure 1 walk-through, generated in R from the notebook's
    own numbers rather than copied from the paper. Left panel, the E-step posterior
    q(Z_i) per sequence at the first iteration; right panel, the two estimates
    converging over 10 iterations.

    It is drawn rather than embedded for two reasons: the published figure is
    copyright Nature Biotechnology and `current_exercises/` is headed for a public
    GitHub repo, and a generated figure cannot fall out of sync with the code above
    it. The repo's other notebooks link images from `popgen.dk/albrecht/open/`, so if
    you would rather show the real figure, put a copy there and the markdown cell
    becomes a one-line image link.

    The E-step weights it draws — 0.45, 0.80, 0.73, 0.35, 0.65 — match the first
    iteration of the paper's Figure 1b exactly, which independently confirms both the
    reference and the corrected `c(0.6,0.5)` initialisation above. Colours are the
    blue/orange categorical pair `#2a78d6`/`#eb6834`, checked colourblind-safe
    (worst all-pairs CVD ΔE 24.7).
  - $\pi$ was introduced as a parameter and fixed to 0.5 in the same sentence; it now
    says it is assumed here and estimated in Bonus 2.
  - the conditional likelihoods $P(X_i|Z_i)$ were labelled "(complete-data) likelihood",
    which would be $P(X_i,Z_i)$.
  - one cell contained nothing but the fragment `$p(X|\theta)$`, evidently a stub for
    the Bonus 1 answer. It is now an answer prompt.
  - `set.seed(1)` added to the simulated coin tosses in part 3. Without it `k` changes
    on every run, so a published html render would not match a student's own numbers.
  - dead commented-out `#abline(v=theta_save)` removed (no such variable).
  - ~15 typos, including 4 in printed output (`Genotypes likehooods:`,
    `log-likelhoods`) and 3 garbled sentences.

  Verified by extracting all 10 code cells and running them through `Rscript`: exit 0,
  no errors. The genotype EM is unchanged at $\theta_T=0.459292$, log likelihood
  -23.63654; the coin EM now converges with A and B the right way round, to
  $\theta_A=0.797$, $\theta_B=0.520$. 25 cells (23 from the source plus the two
  figure cells).

- **#14 `ngs_inference_human.ipynb`** — built 2026-09-16. 99 source cells -> 123,
  with 5 quizzes and questions after 34 of 37 code cells (the 3 without are
  setup).

  **The source was two exercises in one notebook.** Cells 0-18 were a
  continuation of the previous day's *animal* mapping exercise — index the
  wildebeest BAM, call a VCF against the goat reference, filter it, view it with
  tview and mpileup, then the genotype-likelihood shiny app. That is the same
  material as the closing section of #13 `ngs_intro_animal`, so it was **dropped**
  and only the human ANGSD low-depth exercise (cells 19-98) was kept. This also
  means the exercise no longer mixes animal and human data, which is why it sits
  in the human track.

  Consequently `/course/popgen25/NGSInference/fasta` (goat reference) and
  `sams/` (the wildebeest BAM) are not needed and were not copied.

  **Data:** `data/NGSinference/` (2.1 G, 236/236 files) from
  `/course/popgen25/NGSInference/data` — 100 BAMs across 5 populations plus the
  reference and ancestral sequences. This directory was **missed by the original
  data survey**: the notebook writes `/course/popgen25` with a single path
  segment, and the survey's pattern required two.

  **Bugs fixed in the code:**
  - `grep -1 -` and `grep -v -1 -` were counting missing genotypes wrongly: `-1`
    is parsed as grep's context option, not as the pattern. Now `grep -c -- -1`.
  - a bare `print header` line in a bash cell, which is not a command and would
    error. Now a comment.
  - the association run wrote to `Results/$POP.EDAR`, where `$POP` was left over
    from an earlier loop and therefore equal to `NAM` — it silently **overwrote**
    the per-population NAM results. Now writes `Results/NAM_EAS.EDAR`.
  - the p-value cell had the test statistic `2.739244` typed in by hand from an
    old run. It now reads the statistic out of the output file.
  - a malformed markdown link to Matteo Fumagalli's original exercises.
  - doubled slashes in six figure URLs.

  **Paths** centralised per R16: `DATA`, `WORK_DIR`, `REF`, `ANC`, `CHROM`,
  `EDAR_SITE`, `EXAMPLE_SITE`, `NIND`. The work dir was `~/current_folder` and
  `~/popgen25_NGSinference`; it is now `~/ngs_inference_human`, read back from a
  dotfile by the two R cells.

  **Figures** moved into `ngs/figures/` and served from this repo instead of
  `summer2023/NGSinference/`. The two FastQC screenshots used by #12/#13 moved
  there too.

  **Typos:** yesterdays, roughty, thare, "are are called", follwoing, obversed,
  "Native amerians", accross, Fumagilli, "chin protusion", "Column knownEM if
  the estimated".

  **Attribution.** The exercise is a modified version of Matteo Fumagalli's
  low-depth NGS practical (github.com/mfumagalli/Copenhagen). The source
  notebook credited him only in a passing sentence with a broken link. It now
  carries a callout at the top and a **Credit** section at the end naming him,
  linking the original, and citing the EDAR papers (Sabeti 2007, Adhikari 2016)
  and ANGSD (Korneliussen 2014).

- **`statistics/` folder removed (2026-09-16).** Its four files moved:
  `stats_binomial.R` and `stats_normal.R` to `shiny/` because they are Shiny
  apps, and `em_algorithm.ipynb` and `haplotype_frequencies.ipynb` to a new
  `em_algorithms/` folder, since both are about building an EM algorithm from
  scratch.

- **Themes reorganised (2026-09-16).** Five exercises that build an EM algorithm
  are now grouped together, and the gene-flow material folded into admixture:

  | Moved | From | To |
  |---|---|---|
  | #19 `pca_em_human.ipynb` | `pca/` | `em_algorithms/` |
  | #24 `admixture_em_human.ipynb` | `admixture/` | `em_algorithms/` |
  | #37 `sfs_model.ipynb` | `demography/` | `em_algorithms/` |
  | #28, #29 local ancestry | `local_ancestry/` | `admixture/` |
  | #30-#34 f-stats, D-stats, admixture graphs, ChromoPainter, dating admixture | `gene_flow/` | `admixture/` |

  `local_ancestry/` and `gene_flow/` no longer exist. `em_algorithms/` now holds
  5 exercises and `admixture/` holds 11.

- **#19, #24, #37 — the three remaining EM notebooks** — built 2026-09-16.

  | | cells | quizzes | questions after code |
  |---|---|---|---|
  | `pca_em_human.ipynb` | 48 -> 68 | 0 -> 3 | 8/22 -> **22/22** |
  | `admixture_em_human.ipynb` | 36 -> 40 | 0 -> 2 | 12/14 |
  | `sfs_model.ipynb` | 32 | 5 (moved) | 9/10 |

  **All three are self-contained** — they simulate their own data and read no
  files at all, so nothing was needed from `data/` and there were no paths to
  centralise. That is unusual for advBinf exercises, which I had listed as
  `env.sh`-staged; these three are the exception.

  `pca_em_human.ipynb` needed the most work: it is the longest of the three and
  had questions after only 8 of 22 code cells, mostly missing on the function
  definitions (`emu`, `pcangsd`, `halko`, `winsvd`) where the point is to read
  the code and find the E- and M-steps. Those now all carry questions.

  `sfs_model.ipynb` already had 5 quizzes, but they loaded from
  `advBinf/exercises/quiz/`. The five JSON files were copied into
  `em_algorithms/quiz/` and the URLs repointed, so the exercise does not depend
  on the advBinf course folder. Also fixed `allleles` -> `alleles`.

  New quiz bank: `pca_em_missing`, `pca_em_lowdepth`, `pca_em_rsvd`,
  `admixture_em_model`, `admixture_em_practice` — 20 questions, all conceptual,
  covering the traps the exercises are built around (mean imputation pulling
  missing individuals to the origin, the population-frequency prior being mean
  imputation again, label switching, and why a higher likelihood at K=3 cannot
  be used to choose K).

  Still without quizzes: `em_algorithm.ipynb` and `haplotype_frequencies.ipynb`.

- **Setup questions removed (2026-09-16).** Both NGS intro notebooks asked about
  `cp -s` versus `cp` and how many files had appeared in the working folder,
  straight after the cell that links the data in. Those are questions about
  plumbing, not about the analysis, so they are gone (R17).

  Coverage after removal: `ngs_intro_human` 27/31, `ngs_intro_animal` 28/31,
  `ngs_inference_human` 34/37 — the shortfall is now entirely setup cells, which
  is what it should be.
