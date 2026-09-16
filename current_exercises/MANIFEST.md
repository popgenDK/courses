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
| 15 | `genotype_calling_imputation/genotype_calling_and_imputation_human.ipynb` | `advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb` | 2026-09-09 | yes — paths, work folder, quizzes moved in | 2 older copies (see EXERCISES.md) |
| — | `genotype_calling_imputation/quiz/*.json` (3 files) | `/course/data/current_data/popgen25_imputation/quiz/call_genotypes_quiz{1,2,3}.json` | 2025-08-05 | renamed only | none |
| 16 | `genotype_calling_imputation/imputation_human.ipynb` | `advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb` (imputation half) + `chinacourse2026/Day2_Afternoon_Genotype_Imputation.ipynb` (3 sections) | 2026-09-09 | yes — new composition, calling sections dropped | 5 older copies (see EXERCISES.md) |
| 30 | `gene_flow/f_stats_human.ipynb` | `summer2025/exercises/Day3_f_stats.ipynb` | 2025-08-05 | yes — paths, staging dirs removed, course branding | 2 older copies (see EXERCISES.md) |
| 31 | `gene_flow/gene_flow_dstat_animal.ipynb` | `kenya2026/exercises/Day3/Geneflow&Dstat.ipynb` | 2026-08-16 | yes — path cell added, paths, typos | `kenya2026/.../post_course/day3_afternoon_gene_flow_dstat.ipynb` (2026-08-25) |
| 35 | `demography/coalescence.ipynb` | `kenya2026/exercises/Day2/Coalescence_short_WoA.ipynb` | 2026-08-15 | yes — setup cell added, helper repointed | 2 older copies (see EXERCISES.md) |
| 36 | `demography/wright_fisher.ipynb` | `summer2025/exercises/Day1_morning_WrightFisherTutorial.ipynb` | 2025-08-03 | yes — helper repointed, one wording fix | none |
| 38 | `demography/sfs_animal.ipynb` | `kenya2026/exercises/Day2/SFS_WoA.ipynb` | 2026-08-15 | yes — path cell added, paths, folded-axis label | `kenya2026/.../post_course/day2_morning_sfs.ipynb` (2026-08-25) |
| 39 | `demography/psmc_demography_animal.ipynb` | `kenya2026/exercises/Day2/psmc_kenya2026.ipynb` | 2026-08-19 | yes — paths, work folder, chdir bug, quizzes moved in | `kenya2026/.../post_course/day2_afternoon_psmc.ipynb` (2026-08-25) |
| 62 | `demography/psmc_demography_human.ipynb` | `summer2025/exercises/Day5_demography.ipynb` | 2025-08-07 | yes — paths, work folder, chdir bug | 3 older copies (see EXERCISES.md) |
| 63 | `selection/sfs_fst_pbs_human.ipynb` | web page `popgen.dk/albrecht/phdcourse/html/EMBO2021sfs.html` | 2021-03-22 | yes — written as a new SoS notebook | none — new exercise |
| 43 | `relatedness_diversity/relatedness_animal.ipynb` | `kenya2026/exercises/Day5/Related.ipynb` | 2026-08-23 | yes — paths, setup cell, quizzes moved in, broken quiz JSON fixed | `kenya2026/.../post_course/day5_morning_relatedness.ipynb` (2026-08-25) |
| 64 | `relatedness_diversity/relatedness_human.ipynb` | new, written from the GWAS-intro data and analysis | 2026-09-16 | yes — new short exercise | none — new exercise |
| 46 | `gwas/gwas_intro_human.ipynb` | `novCourse2024/1GWASIntro.ipynb` | 2025-07-08 | yes — paths, staging removed, work folder | `bgi23/04.GWASintro_2023_SAIGE.ipynb` (2023-11-09) |
| — | `relatedness_diversity/quiz/*.json` (5 files) | `kenya2026/exercises/Day4/quiz_*.json` | 2026 | renamed; one fixed | none |
| — | `data/scripts/plot2dSFS.R` | `/davidData/albrecht/course/sfs/plot2dSFS.R` | 2023-07-07 | no — copied verbatim | none |
| — | `demography/quiz/psmc_*.json` (3 files) | `/davidData/users/thomas/workshop/psmc_quizzes/` | 2026 | renamed only | none |
| — | `data/scripts/simulateWF.R` | `/course/popgen25/software/simulateWF.R` | 2025 | no — copied verbatim | none |
| 19 | `em_algorithms/pca_em_human.ipynb` | `advBinf/exercises/advBinf_PCA_EM.ipynb` | 2026-09-16 | yes — 3 quizzes, 14 question blocks | none |
| 24 | `em_algorithms/admixture_em_human.ipynb` | `advBinf/exercises/advBinf_admixture_EM.ipynb` | 2026-09-14 | yes — 2 quizzes | none |
| 37 | `em_algorithms/sfs_model.ipynb` | `advBinf/exercises/advBinf_SFSmodel.ipynb` | 2026-09-16 | yes — quiz bank moved in, typo | `advBinf/exercises/SFS.md` (2024-09-20) |
| 23 | `admixture/admixture_low_depth_human.ipynb` | `advBinf/exercises/advBinf_admixture.ipynb` cells 5-70 | 2026-09-14 | yes — split, header, questions, paths | 5 older copies |
| 60 | `admixture/admixture_reference_panel_human.ipynb` | `advBinf/exercises/advBinf_admixture.ipynb` cells 71-93 | 2026-09-14 | yes — split out as its own exercise | none |
| 61 | `admixture/admixture_called_genotypes_human.ipynb` | `chinacourse2026/Day4_admix_eval_LAI.ipynb` cells 2-54 | 2026-07-30 | yes — split, header, questions, paths | none |
| 25 | `admixture/admixture_called_genotypes_animal.ipynb` | `advBinf/exercises/advBinf_admixture_bonus.ipynb` | 2026-09-14 | yes — header, questions, paths, unsafe cell removed | 6 older copies |

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
  `/course/data/current_data/popgenmsc26_exercises/linux/Exercises.zip` (verified present).
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
  -> `/course/data/current_data/kenya2026_anders`, `/course/chinacourse2026/shared` ->
  `/course/data/current_data/chinacourse2026_shared`. Both verified present. The human setup
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

- **Quizzes added to the last two EM notebooks, and a kernel bug fixed
  (2026-09-16).** `em_algorithm.ipynb` gains 3 quizzes (one per section: the two
  coins, allele frequencies from genotype likelihoods, the binomial) and
  `haplotype_frequencies.ipynb` gains 2 (the model and ambiguity; the EM itself).
  `em_algorithms/` now has 15 quiz files and every notebook in it carries
  quizzes.

  **The bug:** `jupyterquiz` is a Python package, but all four of these notebooks
  had `kernelspec: ir` — pure R. A Python quiz cell in an R notebook is executed
  as R and fails. That means the quizzes I had already added to
  `pca_em_human.ipynb` and `admixture_em_human.ipynb` would not have run.

  All four are converted to the SoS kernel that `sfs_model.ipynb` already used,
  with every cell tagged: markdown `SoS`, code `R`, quiz `Python 3 (ipykernel)`.
  Only then do mixed-language cells work.

  A reminder for the exercises still to come: **check the kernel before adding a
  quiz.** A notebook whose `kernelspec` is not `sos` cannot run a Python quiz
  cell.

- **#15 `genotype_calling_and_imputation_human.ipynb`** — copied from the 2026-09-09
  advBinf notebook. The source is a recent, well-written rewrite (5,700 words, 22
  `**Follow-up**` question blocks, a tool comparison table), so **the text, the
  structure and every command are unchanged**. 85 cells, SoS kernel kept (it mixes
  Bash, R and Python cells), committed without outputs. The work is all R4/R8/R11:

  - **R4** — `COURSE_PATH=/course/popgen25` became a single
    `DATA_PATH=/course/data/current_data/imputation`, with `SOFTWARE_PATH=$DATA_PATH/software`. The
    notebook already funnelled every path through one cell that writes an `env.sh`, which
    later cells re-`source`, so this is a two-line change and no other cell contains a
    full path. (It went via the bulk folders `popgen25_imputation` +
    `popgen25_software` first; the cleaned `data/imputation/` folder below replaced both.)
  - **R11** — the working folder `~/advBinfImputation` became
    `~/genotype_calling_imputation_human` (24 further cells updated to match).
  - **R8** — the three quiz files came out of the data folder and now sit in
    `genotype_calling_imputation/quiz/`, loaded by raw-GitHub URL like the `ngs/`
    notebooks. Renamed for meaning: `call_genotypes_quiz1/2/3.json` ->
    `imputation_calling.json`, `imputation_strategies.json`,
    `imputation_concordance.json`. The `cp -sf $DATA_PATH/quiz/*.json .` line is gone.

  **Two dangling symlinks, now fixed at the root.**
  `resources/CEU-chr20-final.b38.txt.gz` (the QUILT2 genetic map) and
  `resources/plink.chr20.GRCh38.rename.map` (the Beagle 5 map) were *symlinks* in the
  data store pointing at `../../software/...`. That resolved in the old
  `/course/popgen25/` layout but dangled in the copied one, because the sibling is
  `popgen25_software/`, not `software/`. They were the only two broken links in the
  whole 293 GB store (`find /course/data -xtype l`). Both are now **real files** in
  `data/imputation/resources/`, verified byte-identical with `cmp`, so R5 holds and the
  notebooks need no workaround — the map paths are back to their original relative
  `resources/...` form.

  **Smoke-tested** with the work folder redirected into a scratch dir: the notebook's
  own setup-check cell reports every program, jar, map, VCF and BAM present with no
  errors, and the cheap data cells run — 30 individuals in the bamlist, mean depth
  3.3-3.9x, a readable pileup, and 1,144 SNP-chip sites against 79,788 panel sites.
  The expensive cells (bcftools ~3 min, ANGSD ~3 min, QUILT2 ~2 min) were not run.

  Reads from the cleaned per-exercise folder `data/imputation/` (see below), not from
  the bulk course folders.

- **The PCA exercises were checked for duplication (2026-09-16).** #21 and #22
  turned out to be **the same exercise listed twice**:

  | Compared | Result |
  |---|---|
  | `advBinf_PCA_bonus.ipynb` vs `Day5_PCA_2.Call_genotype.ipynb` | **33 of 34 cells identical** — only the working folder differs (`~/advBinf/pca` vs `~/popgen25_pca`) |
  | `Kenya2026_PCA.ipynb` (animal) vs `advBinf_PCA_bonus.ipynb` (bonus) | 5 shared cells — genuinely different exercises |
  | `advBinf_PCA.ipynb` (human) vs `Kenya2026_PCA.ipynb` (animal) | 4 shared cells — the MDS/PCA-from-scratch intro |
  | human vs bonus | 1 shared cell |

  Tracing the lineage back, `summer2024-PCA-CalledGenotypes.ipynb` and
  `chinaCourse2025/Day4_Afternoon_PCA_bonus.ipynb` are earlier versions of the
  same wildebeest exercise. **Despite the name, none of them is about calling
  genotypes** — the content is PCAone on wildebeest data, reading in admixture
  proportions, and an IBS tree from plink distances. The "CalledGenotypes" name
  came from a course label that no longer matches the content.

  Merged into one entry under the "bonus" name.

- **Correction (2026-09-16): there are four PCA exercises, not three.** I had
  also listed `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` as superseded by
  `advBinf_PCA.ipynb`. It is not — it is a separate exercise:

  | Compared with `advBinf_PCA.ipynb` | Identical cells |
  |---|---|
  | `summer2025/Day5_PCA_1.ipynb` | 74 of 78 — genuinely a version of it |
  | `chinaCourse2025/Day4_Afternoon_PCA_main.ipynb` | 41 — an older version |
  | `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` | **23 — a different exercise** |

  #18 teaches PCA from first principles: MDS, the SVD worked by hand, then
  low-depth data. `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` instead starts
  from an **LD-pruned file of called genotypes**, runs PCAone, plots by
  population and super-population against the admixture proportions, and then
  does PC-based selection with PCAngsd. It is restored as #22
  `pca_called_genotypes_human.ipynb`.

  So the "called genotypes" name belongs to the **human** exercise. The
  wildebeest lineage carried that name in `summer2024-PCA-CalledGenotypes.ipynb`
  and `Day5_PCA_2.Call_genotype.ipynb`, but its content was never about calling
  genotypes.

  `pca/` has four exercises: human from first principles (#18), animal (#20),
  animal bonus (#21), human called genotypes (#22).

- **#16 `imputation_human.ipynb` — built as the imputation half, not as a copy of its
  recorded source (2026-09-16).** `EXERCISES.md` had it as "the separate imputation half
  of #15", but its source,
  `chinacourse2026/Day2_Afternoon_Genotype_Imputation.ipynb`, is an older and thinner
  copy of the *entire* #15 exercise: the same CEU chr20 2-5 Mb data, the same five tools
  in the same order, the same three quiz json files, the same `find_match` and
  `estimate_gt_tools` helpers, and the same NIPT bonus. 2,271 words with no follow-up
  questions, against #15's 5,704 with 22. Its headings are also duplicated ("Genotype
  calling without Imputation" appears twice) and typo'd ("pata preparation").

  A faithful copy would have put two near-identical *human* notebooks side by side,
  which R2 only allows for the human/animal split. **Asked and confirmed: build it as
  the imputation-only half.** So it is a new composition, 68 cells:

  - the imputation sections come from #15, the newest and much fuller version of the
    same material (R1)
  - the SNP and genotype **calling** sections are dropped — they are #15's job, and the
    notebook says so at the top and links to it
  - it opens with the input preparation those sections used to provide: the bamlist and
    **one ANGSD run** for the genotype likelihoods Beagle 4.1 needs. ANGSD's own
    genotype calls are kept as the **no-imputation baseline** for the comparison, which
    is the question the exercise is really asking
  - the comparison therefore has four methods, not five (`angsd`, `beagle4gl`,
    `quilt2`, `beagle5gl`); `bcftools` and the UpSet SNP-discovery plot stay in #15
  - `find_match` is defined here, since it lived in a calling cell in #15
  - quizzes 2 and 3 travel with it (`imputation_strategies`, `imputation_concordance`);
    quiz 1 is about calling and stays with #15

  **The three things only the chinacourse version had are carried over**, since that is
  the whole reason not to drop it:
  - the reference-panel and genetic-map **input formats** section (IMPUTE hap/legend,
    the map columns, bamlist, and the optional phasefile/posfile truth format), with
    "pata preparation" fixed
  - the `example.vcf.gz` walk-through, so the output fields are looked at in a real file
  - the "how about imputing the **fetus**?" question, expanded into three questions at
    the end of the NIPT bonus

  Paths, work folder (`~/imputation_human`) and quiz URLs follow #15: a single
  `DATA_PATH=/course/data/current_data/imputation` covering reads, VCFs, reference, maps and software.

- **Sample count fixed in #15 and #16.** Both said the study samples were **33** in four
  places while also saying **30** in six others. The bamlist has 30
  (`ls bams/NA*.bam | wc -l` = 30, verified). The 33 was inherited from the chinacourse
  notebook's opening line, "We will work on 33 European samples". All four now say 30.

- **PCA renamed by data type, and the theory section split out (2026-09-16).**
  Names now say which kind of data the exercise runs on:

  | # | New name | Data |
  |---|---|---|
  | 58 | `pca_mds_and_svd.ipynb` | none — a matrix typed in from the slides |
  | 18 | `pca_low_depth_human.ipynb` | genotype likelihoods (beagle) |
  | 59 | `pca_low_depth_selection_human.ipynb` | genotype likelihoods (beagle) |
  | 22 | `pca_called_genotypes_human.ipynb` | LD-pruned plink genotypes |
  | 20 | `pca_called_genotypes_animal.ipynb` | wildebeest plink genotypes |

  **`chinacourse2026/Day4_Afternoon_PCA_1.ipynb` becomes two exercises**, cut at
  its `# PC-based selection` heading: the first half is PCAone on called
  genotypes (#22), the second is `pcangsd --selection` on genotype likelihoods
  (#59). It genuinely did both, which is why it resisted classification.

  **The MDS/PCA-by-hand section becomes #58.** It was repeated in `advBinf_PCA`
  (cells 2-25), `Kenya2026_PCA` and again near the end of the China notebook.
  One copy now, in `pca/`, needing no data.

  **Found while classifying:** `Kenya2026_PCA.ipynb` carries the heading "PCA for
  low depth sequencing using PCAngsd" but contains **no PCAngsd command** — it
  runs `PCAone` on plink files throughout. The heading is left over from an
  earlier version and misdescribes the exercise. To be fixed when #20 is built.

  This takes the plan from 55 to 57 exercises.

- **Correction (2026-09-16): #21 was the same exercise as #20.** I had kept
  `advBinf_PCA_bonus.ipynb` and `Kenya2026_PCA.ipynb` apart because only 5 cells
  matched exactly. Comparing them properly rather than by exact string equality:

  | Threshold | Bonus cells matching the kenya wildebeest section |
  |---|---|
  | identical | 12 of 34 |
  | >= 0.8 similar | 19 of 34 |
  | >= 0.7 similar | 20 of 34 |

  The 5-cell figure was an artefact of comparing whole cells for exact equality;
  most of the shared cells differ only in a path or a variable name.

  The real differences: the kenya version opens with the MDS/by-hand section
  (now #58) and takes its LD-pruned file ready-made from the admixture exercise,
  while the bonus version **computes the LD pruning itself** with `PCAone --ld`,
  adjusting the LD measure for population structure, then re-runs the PCA on the
  pruned data — about 7 cells the kenya version does not have.

  Merged into #20, built from the bonus notebook as the newer and more complete
  base. `pca/` is now exactly what was asked for: two called-genotype exercises
  (#20 animal, #22 human), two low-depth ones (#18, #59), plus the theory (#58).

  **Lesson for the remaining folders:** exact-match cell comparison understates
  duplication badly. Use a similarity ratio.

- **The five PCA exercises built (2026-09-16).**

  | # | File | cells | quizzes | questions after code |
  |---|---|---|---|---|
  | 58 | `pca_mds_and_svd.ipynb` | 30 | 1 | 9/11 |
  | 18 | `pca_low_depth_human.ipynb` | 41 | 2 | 9/12 |
  | 59 | `pca_low_depth_selection_human.ipynb` | 41 | 2 | 10/11 |
  | 22 | `pca_called_genotypes_human.ipynb` | 30 | 4 | 6/8 |
  | 20 | `pca_called_genotypes_animal.ipynb` | 44 | 2 | 11/13 |

  **`advBinf_PCA.ipynb` split three ways.** Cells 2-25 became #58, 26-52 became
  #18 and 53-77 became #59. That last section is a PC-based selection scan, and
  comparing it with the China notebook's selection half showed 20 of 24 cells
  matching at 0.7 similarity or better — the same exercise again. #59 is
  therefore built from advBinf (the newer of the two) and supersedes the China
  half, which leaves the China notebook contributing only #22.

  **Data:** `data/PCA/{human_lowdepth,human_called,animal}` (461 M) from
  `/course/popgen25/pca` plus the China plink files. `/course/popgen25/pca` was
  a **third directory the original survey missed**, for the same reason as
  `NGSInference`: the notebooks write `/course/popgen25` with one path segment.

  **Fixed while building:** four quizzes inherited from the course folders had
  `indiviudals`, `indvidual` and `Who many` in their questions; `informaiton`,
  `prevous`, `keept`, `enviroment`, `Lest zoom`, "summaries the fist column" and
  "We can not perform the PCA again" (which should read "can now") in the
  notebooks. One quiz cell in #18 did `os.chdir` into the shared data folder
  before loading its quiz, which both broke R16 and would have left the notebook
  pointing at the wrong directory. The locus-zoom call in #59 hardcoded
  `~/data_folder/geneticMap/hg38`; it now reads a `GENETIC_MAP` variable set in
  the setup cell and passed to R through a dotfile, pointing at
  `data/geneticMap/hg38`.

  New quizzes: `pca_selection_qq`, `pca_selection_hit`, `pca_animal_structure`,
  `pca_ld_pruning`.

- **`data/imputation/` — cleaned per-exercise data folder (2026-09-16), 3.4 G.** Both
  #15 and #16 now read everything from this one folder, so their path cell names a
  single directory instead of two bulk course folders:

  ```
  data/imputation/
    bams/        31 low-depth CEU bams + indexes (30 study samples + the NIPT sample)   223 M
    vcfs/        reference panel, truth set, fake SNP chip, precomputed QUILT2, example  28 M
    resources/   GRCh38 reference + .fai, and the QUILT2 and Beagle 5 genetic maps      3.1 G
    software/    beagle 4.1 and 5.5 jars, and the QUILT distribution (QUILT2.R)          64 M
  ```

  Assembled with real `cp` from `popgen25_imputation/` and `popgen25_software/` (R5).
  **Zero symlinks and zero dangling links** (`find -type l` / `-xtype l`): the two
  broken map links were dereferenced into real files on the way in. Verified by `cmp`
  against the originals — the 3.1 G reference genome and both maps are byte-identical,
  and the file counts match (62 bam files, 10 vcf files, 172 QUILT files).

  Deliberately **not** copied from `popgen25_imputation/`: `quiz/` (the quiz json now
  travels with the exercise, R8), `outputs/` (7.9 M of precomputed results the
  exercises do not read) and `Day2_Imputation.ipynb` (a course notebook, not data).
  Only the software the two exercises actually invoke was taken from
  `popgen25_software/`, not the whole 1.8 G tree.

  The bulk folder `popgen25_imputation/` is untouched and still in the data-copy list;
  other exercises may need it.

  **This is why the folder is self-contained.** `popgen25_software/` was moved into
  `data/current_data/` by the reorganisation recorded below, hours after these notebooks
  were built. Because `data/imputation/software/` holds its own copy of the two Beagle
  jars and the QUILT distribution, neither notebook referenced the moved folder and both
  survived the move without an edit — verified after the fact. Pointing
  `SOFTWARE_PATH` at `popgen25_software/`, which is what the first version of #15 did,
  would have broken both.

  **Gotcha found by running it: `cp` breaks ANGSD's fasta index check.** A plain `cp` of
  the reference and its `.fai` gives both the same mtime to the second, and ANGSD wants
  the index *strictly* newer:

  ```
  -> fai index file: '...fa.fai' looks older than corresponding fastafile: '...fa'.
  -> Please reindex fasta file
  ```

  In the source folder the two were 8 hours apart, so this only appeared after the copy.
  Fixed with `touch` on the `.fai` — the file is byte-identical to the original, it just
  needed a newer timestamp. **Any future per-exercise data folder that copies a
  reference genome needs the same treatment**, or ANGSD refuses to start.

  The `.bai` and `.csi` indexes in `bams/` and `vcfs/` share their file's mtime exactly,
  because they were copied with `cp -a` from a source where they already did. samtools
  and bcftools do not enforce the comparison, and the full pipeline runs, so those were
  deliberately left alone rather than restamped.

- **Data moved into `data/current_data/` (2026-09-16).** Everything copied during
  the consolidation now sits under `current_data/`, keeping the folder names the
  exercises already use. 291 G, moved rather than copied since it is the same
  filesystem, so it was instant.

  Left in place at `/course/data/`: `1000G/` (predates this work). The imputation
  data was held back at the time and moved in once those notebooks were finished
  — see below.

  Every path in every built exercise was repointed and then checked to resolve:
  9 files, 15 paths. The two imputation notebooks were correctly left untouched,
  since their data did not move.

- **Every notebook now opens the same way (2026-09-16), per R19:** title, purpose,
  what you will do, and the data described concretely.

  **Four PCA notebooks had no introduction at all** — they opened on `## Setup`,
  because the builder put the setup cell ahead of the intro. Titles and full
  headers added.

  The data is now described from **the actual files**, not from the old prose:

  | Exercise | Data as now described |
  |---|---|
  | `pca_called_genotypes_human` | 192 individuals, 16 populations x 12, four super-populations, 317,850 LD-pruned SNPs |
  | `pca_low_depth_human` | 435 individuals, ASW/CEU/CHB/MXL/YRI with counts, beagle genotype likelihoods |
  | `pca_low_depth_selection_human` | 424 individuals, CEU/GBR/IBS/TSI with counts |
  | `pca_called_genotypes_animal` | 73 blue wildebeest, 7 localities with counts, 990,980 SNPs |
  | `ngs_intro_human` | NA19238, YRI Ibadan Nigeria, 171,880 read pairs of 100 bp, chr21 |
  | `ngs_intro_animal` | CCTauTzS_8872 blue wildebeest, 83,506 read pairs of 150 bp, mapped to a **goat** reference |
  | `ngs_inference_human` | 100 BAMs, 20 each from LWK/TSI/CHB/PEL/NAM, ~2-3x depth |
  | the five EM notebooks | simulated, with what is simulated spelled out |

  **A factual error fixed:** `ngs_inference_human` said "we will use 40 BAM files
  ... and 10 BAM files of Latinos". The folder holds **100**, 20 from each of five
  populations. The prose had drifted from the data.

  Population codes are now expanded — `CEU` is given as "Utah residents with
  northern/western European ancestry", `TSI` as "Toscani in Italy" — since a code
  alone means nothing to a student. Two things that will bite later are called out
  where they apply: the wildebeest sample sizes are very unequal (28 from one
  locality, 3 from another), and the animal NGS intro maps a wildebeest to a goat
  reference, which is explained rather than left as a surprise.

  Your two imputation notebooks already met R19 and were left untouched.

- **Imputation data moved in too (2026-09-16), now that those notebooks are
  done.** `data/current_data/imputation/` (`bams/ resources/ software/ vcfs/`,
  3.4 G) and the older bulk copy `popgen25_imputation/` (3.3 G).

  Both imputation notebooks had already been consolidated onto the single clean
  `imputation/` folder, so repointing them was one path each. Verified afterwards:
  every `/course/data/...` path in every built exercise resolves.

  `/course/data/` now holds just `current_data/` and `1000G/`.

- **#30 `f_stats_human.ipynb` and #31 `gene_flow_dstat_animal.ipynb` (2026-09-16).** The
  two halves of the f-statistics / D-statistics material, human and animal, both in
  `gene_flow/`. Names kept from the build list, which already reads as "f-stats" and
  "gene flow". Both are pure R under the SoS kernelspec; the teaching content, the
  statistics and the question blocks are unchanged.

  **#30, from the 2025-08-05 summer2025 notebook** (37 cells, `admixtools`): f2 heatmap
  and neighbour-joining tree, admixture f3, outgroup f3, f4 treeness and symmetry tests,
  and qpAdm, with four `## Task` cells left empty for the student on purpose.
  - its one bash cell existed only to `mkdir` a work folder and create the
    `~/current_folder` and `~/data_folder` symlinks that **R4 explicitly bans**. Replaced
    with an R path cell, so the notebook is now all R and reads `DATA` directly.
  - R11: dropped "CPH popgen 2025" from the practical's title.

  **#31, from the 2026-08-16 kenya2026 notebook** (13 cells): the ABBA/BABA D statistic
  computed from scratch in R on simulated CHIMP/AFR/EUR/NEA data, then f4 on real
  wildebeest data with `admixtools`.
  - it had **no path cell** — three absolute paths sat in three different cells. Added
    one at the top, matching #30 and the `ngs/` notebooks.
  - R11: "the D-statistic mentioned today" -> "above"; "what you learn today" ->
    "what you have learned"; fixed "why do use hartebeest" and the species spelling
    "Wilderbeet" -> "wildebeest".

  **Two data findings:**
  - **#30's data was missed by the original survey.** `/course/popgen25/dstats` (108 M)
    was not in the data-copy list, because the notebook reaches it through the
    `~/data_folder` symlink and so contains no literal `/course/...` path — the same
    blind spot that hid the `geneticMap` dependency. Copied to
    `data/f_stats/` (14,558/14,558 files, no symlinks).
  - **#31 depended on an instructor home directory**, `/davidData/users/thomas/workshop`
    (77 G), for a 1.5 M subfolder. `data/geneflow/` holds just the three things the
    exercise reads, 3.3 M in total, so the exercise no longer needs the 77 G tree at all.

  **Both were run end to end** (pure R, so every code cell executes), exit 0, no errors:
  - #30: 98 populations and 50,493 SNPs loaded; admixture f3 strongly negative for
    African Americans against European references (z about -31); the f4 treeness test
    non-significant in the Yoruba configuration (z = 1.64) while the other two are
    (z = 76.5 and 74.6), which is exactly the contrast its questions ask about; qpAdm
    weights and `popdrop` produced.
  - #31: 81,486 variants read, 81,259 biallelic retained, D = 0.187 — the introgression
    signal the exercise is built around — and the 42-row wildebeest f4 table with its
    plot.

- **The four admixture exercises built (2026-09-16),** mirroring the PCA folder:
  two on called genotypes, two on low depth.

  | # | File | cells | quizzes | questions after code |
  |---|---|---|---|---|
  | 23 | `admixture_low_depth_human.ipynb` | 75 | 4 | 18/21 |
  | 60 | `admixture_reference_panel_human.ipynb` | 29 | 1 | 8/10 |
  | 61 | `admixture_called_genotypes_human.ipynb` | 56 | 4 | 21/22 |
  | 25 | `admixture_called_genotypes_animal.ipynb` | 56 | 1 | 20/21 |

  **`advBinf_admixture.ipynb` was two exercises.** Cells 5-70 are NGSadmix on
  genotype likelihoods; cells 71-93 are fastNGSadmix, placing a **single**
  individual against a fixed reference panel of 7 populations. Different method,
  different question, so they are now #23 and #60.

  **Data:** `data/current_data/admixture/{human_lowdepth,human_refpanel,human_called,animal}`,
  631 M.

  **A dangerous cell removed.** The animal notebook opened with
  `mkdir -p ~/advBinf/admixture; cd ~/advBinf/admixture; rm -rf *`. Rewriting the
  hardcoded paths left `cd` with no argument immediately before `rm -rf *`, which
  would have deleted the contents of whatever directory the shell happened to be
  in. The cell is replaced — the setup cell already creates and enters the working
  folder. Worth remembering when rewriting paths in any notebook that clears its
  workspace.

  Also removed: three cells in the human notebook that read a
  `~/.sysu_admixture_paths` config file and re-did the setup, superseded by the
  standard setup cell (R16).

- **#35 `coalescence.ipynb` and #36 `wright_fisher.ipynb` (2026-09-16).** Both in
  `demography/`, both species-neutral so neither takes a `_human`/`_animal` suffix (R2).
  Pure R under the SoS kernelspec, no data files at all — everything is simulated in the
  notebook. The teaching content is unchanged.

  **They share one helper.** Both `source()`d
  `/course/popgen25/software/simulateWF.R`, which defines `WF_twoalleles`,
  `WF_manyalleles` and `track_lineages`. Per **R14** a sourced function library is not an
  exercise and lives in the data folder, so it is now `data/scripts/simulateWF.R`
  (copied verbatim, `cmp`-identical), alongside `admixFun.R`, `newPlotPlink.R` and
  `online.R`. Before this it existed only inside the 1.8 G `popgen25_software/` tree,
  which is a software directory rather than a script library.

  **#35, from the 2026-08-15 kenya2026 notebook** (24 cells): simulating a coalescence
  tree for five samples step by step — sampling which lineages coalesce, the exponential
  waiting times, and then the whole thing in a loop — followed by a backwards-in-time
  Wright-Fisher section using `track_lineages` at N = 10, 7 and 20.
  - it opened with a bare `## Environment setup` heading and **no setup cell**, while the
    helper was `source()`d from a hard-coded path much further down (cell 19). Added the
    path cell under that heading; the later cell now points at it.
  - the two figures are external images from Fernando Racimo's public `CopenhagenTutorial`
    repo. Kept, and both checked live (HTTP 200) — R11 keeps external links.

  **#36, from the 2025-08-03 summer2025 notebook** (20 cells): Wright-Fisher forwards in
  time with two alleles and with many alleles, then backwards in time with
  `track_lineages`. Four empty code cells are left for the student on purpose.
  - "Start running the R console and load the following R file" became "Run the cell
    below to load the R functions used throughout this exercise" — it is a notebook.

  **Both were run end to end** (every code cell, `set.seed(42)` for reproducibility),
  exit 0, no errors:
  - #35: the loop performs 4 coalescence events for 5 samples and ends on node 9, which
    is the correct `n-1` events and `2n-1` final node; waiting times sampled from the
    right exponential rates.
  - #36: `WF_twoalleles(5,15)` returned blue counts ending at 10, i.e. `2N` — the blue
    allele fixing, which is the outcome the first question asks students to tally.

  **Note the overlap, which is in the sources, not introduced here:** #35's
  "Wright Fisher / thinking backwards in time" section repeats #36's section 3 — the same
  three `track_lineages` calls at N = 10, 7 and 20. The build list keeps both exercises,
  so both keep the section.

- **The two local ancestry exercises built (2026-09-16).**

  | # | File | cells | questions after code |
  |---|---|---|---|
  | 28 | `local_ancestry_flare_mosaic_human.ipynb` | 108 | 28/34 |
  | 29 | `local_ancestry_hapla_human.ipynb` | 57 | 21/27 |

  **The three local-ancestry sources share nothing.** Compared cell by cell, the
  summer2025 notebook (FLARE/MOSAIC), the advBinf hapla notebook and the China
  tail have **zero** cells in common — three different tools. The plan had #28
  superseding the summer2025 notebook, which was wrong: they are different
  exercises, not versions of each other. The China tail is a 14-cell "short look"
  at hapla/fatash from precomputed files, fully covered by #29.

  **Data:** `data/current_data/local_ancestry/{flare_mosaic,hapla_cattle}`, 216 M.
  `/course/popgen25/LocalAncestry` was a **fourth directory the original survey
  missed** — same cause as the others, a single path segment after `/course`.

  The hapla exercise ends on a real **cattle** dataset (*Bos taurus*, chr25, 314
  individuals, BosTau9), which the header now says. `BosTau9.ids` is produced by
  `hapla cluster --out`, not an input, so nothing is missing.

  MOSAIC's scratch directory `~/.cache/fastfiles` is now a `FASTFILES` variable set
  in the setup cell rather than hardcoded mid-notebook.

- **#38 `sfs_animal.ipynb` (2026-09-16).** From the 2026-08-15 kenya2026 notebook, into
  `demography/`. 19 cells, SoS kernelspec — it mixes Python (the neutral 1/i expectation
  and its barplot) with R (the real spectra), which is why the kernelspec stays SoS. The
  teaching content and the six question blocks are unchanged.

  - **it had no setup cell**, and hard-coded `/course/kenya2026/harvi/sfs/inputdata/`
    *inside two different functions*, built up with `paste0`. Added a path cell after the
    title and both reads now go through `DATA`, so the directory appears once.
  - **Data: cleaned folder `data/sfs/`** (291 M, 5 files, all `cmp`-identical). The source
    `sfs/inputdata/` held exactly the 5 VCFs the exercise reads, so this is a
    straight lift out of the 293 M `kenya2026_harvi/` bulk folder.
  - **Axis label corrected.** The two wildebeest panels plot a *folded* spectrum —
    `calc_folded()` takes `pmin(alt_count, nchr - alt_count)` — but both were labelled
    "Derived allele count". A folded spectrum has no derived/ancestral distinction, which
    is the very thing Questions (5) asks about ("what information is lost when an SFS is
    folded?"). Now "Minor allele count" in both panels. The simulated panels, which are
    genuinely unfolded, keep "Derived allele count".

  **Run end to end** (all 8 R cells, including the 115 M and 176 M wildebeest VCFs),
  exit 0, no errors. The spectra come out as the exercise intends:
  - scenario a decays roughly like 1/i (0.236, 0.113, 0.077, ...) — the neutral shape
  - scenario b has a clear singleton excess (0.355) — growth
  - scenario c is extreme and ragged (0.559 singletons, sparse tail), which is what
    Questions (3) item 3 asks about ("does the curve for c look unusual?")
  - no `NA` bins anywhere, so the genotype recoding is not silently dropping
    multiallelic sites the way it would if the input had `2|0`-style codes
  - wildebeest: 159,687 segregating sites for black against 76,928 for blue, and black's
    folded spectrum *rises* after bin 1 (11384, 7604, 7454, 8320, 9257, 11056, ...) while
    blue's decays monotonically — the contrast Questions (4) is built on

  Python cells were not executed (no data, and they need the `sos` kernel to run in
  order); they are 12 lines of numpy/matplotlib with no inputs.

- **#39 `psmc_demography_animal.ipynb` and #62 `psmc_demography_human.ipynb`
  (2026-09-16).** The two PSMC exercises, animal and human, both in `demography/`. 62
  and 51 cells, SoS kernelspec (Bash plus Python for viewing the plots). The teaching
  content, the PSMC commands and the question blocks are unchanged.

  They are the **same exercise on different species**, which is why both are kept (R2):
  simulated data first, then #39 goes to three wildebeest samples and #62 to two 1000
  Genomes individuals (NA12718 CEU, NA19471 Luhya). #62 ends with a wildebeest bonus that
  #39 covers properly.

  **One shared data folder, `data/psmc/` (284 M),** rather than one per exercise, since
  they overlap almost completely — the same layout idea as `data/NGSintro/`:

  ```
  data/psmc/
    data/       simulated + 1000 Genomes psmcfa/psmc, and the wildebeest bcf for #62's bonus
    animal/     wildebeest psmcfa and precomputed psmc results (from kenya2026_rasmus)
    images/     popsize, bootstrap, 1kg_chr1, NA12718 figures
    scripts/    vcf2psmcfa.py
    software/   the psmc binary and utils/psmc_plot.pl
  ```

  Between them the two notebooks previously read from **five** places, including
  `/course/popgen25/software`, `/course/kenya2026/rasmus/psmc/data` and
  `/davidData/users/thomas/workshop/psmc_quizzes/`. Now each names one folder.

  **Three more symlinks pointing outside the data store, found and fixed.** The copied
  `popgen25_demography/data/` contained:
  - `NA12718_chr1.psmc` and `NA19471_chr1.psmc` -> `/course/popgen24/shyam/psmc_tutorial/analysis/`
    (a *different instructor's* directory), and
  - `CTauTzS_8872.chr27.filtered.bcf.gz` -> `/davidData/users/thomas/workshop/`

  All three are **load-bearing**: the two `.psmc` files are what #62 offers students who
  do not want to wait 10 minutes for PSMC, and the bcf is the input to its wildebeest
  bonus. They resolved only because those directories still exist. Dereferenced into real
  files in `data/psmc/data/`, so the folder now has **zero symlinks**.

  **A real bug in both:** the Python setup cell did `os.chdir("demography")` — a
  *relative* path, while the bash cell next to it did `cd ~; mkdir -p demography; cd
  demography`. The Python cell therefore only worked if the notebook server happened to
  start in `$HOME`, and silently put the Python cells in a different directory from the
  bash cells otherwise. Both now use
  `os.chdir(os.path.expanduser("~/<exercise name>"))`.

  **Other changes:**
  - work folder `~/demography` became `~/psmc_demography_animal` and
    `~/psmc_demography_human`. They shared one folder before, and now that both exercises
    exist side by side they would have overwritten each other's `s1_1.psmcfa`,
    `combined_coarsePattern.psmc` and plots.
  - `TOOLS_PATH` pointed at the whole 1.8 G software tree, with `PSMC=${TOOLS_PATH}/psmc/psmc`.
    It now points at `data/psmc/software`, so `PSMC=${TOOLS_PATH}/psmc` and
    `PSMC_PLOT=${TOOLS_PATH}/utils/psmc_plot.pl`.
  - the Python image cells cannot see the bash variables, so the setup cell defines
    `IMAGES_PATH` on the Python side too, and the four `mpimg.imread` calls use it.
    A comment says why the path is defined twice.
  - #39's three quizzes came out of `/davidData/users/thomas/workshop/psmc_quizzes/` into
    `demography/quiz/`, loaded by raw-GitHub URL (R8), renamed
    `psmc_quiz1_input_and_heterozygosity.json` -> `psmc_input_heterozygosity.json` and so
    on. 13 questions across the three.

  **Smoke-tested against the new folder** — the whole toolchain with real data, exit 0:
  `vcf2psmcfa.py` on `threeinds.recode.vcf.gz` produced a valid psmcfa; `psmc -p "4+5*3+4"`
  ran to completion (RS/PA records present); `psmc_plot.pl` plus `pdftoppm` produced the
  `.pdf` and `-1.png` the notebooks display; the three precomputed wildebeest `.psmc`
  files copy; the human psmcfa, the two precomputed `.psmc` and the 176 M bcf all read;
  and all four figures load in matplotlib.

  **Not run:** the full notebooks under the SoS kernel, and the real 10-minute PSMC runs
  on the 1000 Genomes and wildebeest samples. Every command in those cells is the same
  binary and options exercised above, on data verified present.

  **Trimmed to what is actually read (by request, 2026-09-16).** `data/psmc/data/` came
  across with 22 files, 10 of which neither notebook names. Those 10 were removed:
  `exercise_2.vcf.gz` (93 M, belonging to some other exercise), the chr5 and chr21
  alternatives of the 1000 Genomes VCFs and psmcfa files, and
  `1kg_2samps_chr1.recode.vcf.gz` (the with-indels version of a file the exercise only
  uses in its `_noindels` form). **211 M freed**; the folder went 410 M -> 209 M and the
  whole of `data/psmc/` 485 M -> 284 M.

  Checked before deleting, not after:
  - only the two PSMC notebooks read `data/psmc/` at all, so no other exercise could be
    affected
  - usage was matched on the **exact** filename, not a stem — a stem match wrongly counts
    `1kg_2samps_chr1.recode.vcf.gz` as used because its name is a prefix of
    `1kg_2samps_chr1_noindels.recode.vcf.gz`
  - **all 10 still exist in `popgen25_demography/data/`**, so every one is a single `cp`
    away if an exercise turns out to want it

  Re-verified after deleting: the full toolchain still runs (vcf -> psmcfa -> psmc ->
  plot -> png), all 12 files the notebooks name are present, all 7 wildebeest files, all
  4 figures, and the bcf still opens under `bcftools view`.

- **Gaps filled in the two gene-flow exercises (2026-09-16).**

  | | before | after |
  |---|---|---|
  | `f_stats_human` quizzes | 0 | 2 |
  | `f_stats_human` questions after code | 6/16 | 13/16 |
  | `gene_flow_dstat_animal` quizzes | 0 | 1 |
  | `gene_flow_dstat_animal` questions after code | 3/5 | 5/5 |

  **`f_stats_human` opened on `# setup enviroment and data`** — the setup came before
  the title, the same defect the PCA notebooks had, and with a typo. It now opens with
  the R19 header and the setup follows.

  Data descriptions taken from the files: the AADR subset is **1,646 individuals in 98
  populations, 53 modern and 45 ancient**, with the `.HO` and `_HG` suffix convention
  explained. The D-statistic exercise uses **AFR 100, EUR 100, NEA 2, CHIMP 1** for the
  simulated part and the seven blue wildebeest localities plus black wildebeest and
  hartebeest for the real part.

  The new quizzes target what these statistics actually license you to conclude: that a
  negative admixture $F_3$ is evidence while a positive one proves nothing, that qpAdm
  rejects models rather than confirming them, that a negative qpAdm weight invalidates a
  model whatever its p-value, and why the blocked jackknife standard error is what makes
  any of it interpretable.

- **#63 `sfs_fst_pbs_human.ipynb` — new exercise, built from a web page (2026-09-16).**
  Requested directly: an SoS notebook version of
  `https://www.popgen.dk/albrecht/phdcourse/html/EMBO2021sfs.html` (page updated
  2021-03-22). It was not in the repo in any form, so nothing is superseded. Put in
  `selection/`, which otherwise holds only animal exercises, so this is its human
  counterpart. 50 cells: 26 markdown, 12 Bash, 12 R.

  The content follows the page section by section — SFS in 1D and 2D from genotype
  likelihoods, pairwise Fst, a PBS sliding-window scan against a genome-wide background,
  and the chromosome 5 region. Every command, filter and option is the page's, and all of
  its questions are kept.

  **Adapted for a notebook rather than a terminal + separate R session:**
  - the page's `$ThePath=/home/albrechtsen/embo2021` layout became one
    `DATA_PATH=/course/data/current_data/sfs_fst_pbs`, written to an `env.sh` that every
    one of the 12 bash cells re-`source`s — the same robustness pattern as #15/#16, which
    matters here because the steps are deeply dependent on each other
  - `angsd` and `realSFS` are on PATH on this server, so `ANGSD`/`REAL` point at those
    rather than at a per-course `prog/angsd/` build
  - the three `x11()` calls are gone: in a notebook each plot belongs in its own cell, so
    the 2D spectra are three cells instead of one cell plus two X windows
  - "open R", "close R" and "Remember to close R" instructions dropped — the R cells are
    part of the notebook
  - `cp $ThePath/sfs/plot2dSFS.R .` dropped; the helper is `source()`d from
    `data/scripts/plot2dSFS.R` (R14), so nothing is copied into the working folder
  - the R cells `setwd()` into the working folder, since the bash `cd` does not reach them

  **Two fixes to the exercise itself:**
  - the genotype-calling bonus reused the output prefixes `yri`, `jpt` and `ceu` — the
    same as the saf runs in section 2 — so it overwrote `yri.arg`, `yri.mafs.gz` and so on
    from the earlier section. Changed to `yriGeno`/`jptGeno`/`ceuGeno`, and the R cell
    reads those. It also now passes `-anc $ANC`, which `-doMajorMinor 4` (use the
    ancestral allele) requires.
  - one of the four linked PDFs, `oulu2016/web/PBS.pdf`, is **404**. It illustrated the
    background Fst/PBS distribution, which students now generate themselves a few cells
    earlier, so the text points at their own histograms and maxima instead. The other
    three links are live (checked, HTTP 200) and kept.

  **Data: `data/sfs_fst_pbs/` (2.4 G)** copied from `/davidData/albrecht/course/sfs/` —
  the human reference and the chimp ancestral fasta with their indexes, `smallerbams/`
  (30 bams + indexes, 10 per population) and `chr5_33M_v2/` (30 bams). No symlinks.
  - the `.fai` files hit **the same ANGSD timestamp trap** recorded under R5: ANGSD
    refused to start with "Please reindex fasta file" until the copied indexes were
    `touch`ed. This is now the second time that has bitten, on unrelated data.
  - the chr5 bams have **no `.bai`**, which turned out not to matter: ANGSD reads them
    sequentially since no region is requested. Verified on a 2-bam test run, 648,443 sites
    retained.
  - `precomputed/` holds the saf files from a verification run, so the page's "if it takes
    too long, copy the results" shortcut works. The page pointed at
    `$ThePath/run/`, which does not exist here.

  **Verified by running the whole pipeline on the real data.** The answers the questions
  ask for all come out correctly:
  - 1.61 M sites retained per population; YRI has the most singletons (0.083 of the
    spectrum, against 0.054 for JPT and 0.039 for CEU) and the largest fraction of
    segregating sites (0.0050 vs 0.0035 and 0.0034), giving Ne ~23,500 against ~16,300 and
    ~16,000 — so "which population has the largest population size?" resolves to the
    African one, as it should
  - global weighted Fst: **JPT-CEU 0.085** (most closely related), YRI-CEU 0.121,
    **YRI-JPT 0.157** (most distant) — the textbook ordering the two Fst questions expect
  - the background scan produced 104 windows of 50 kb, and the chromosome 5 region 66
  - **the chr5 region is emphatically extreme**, which is the payoff of the exercise:
    `PBS_CEU` peaks at **0.490** at chr5:33,945,000, against a genome-wide background
    maximum of 0.145 — **3.4x** the most extreme background window. In the same window
    `PBS_YRI` is 0.210 and `PBS_JPT` 0.125, so "in which population is this locus under
    selection?" answers cleanly: CEU. And both CEU comparisons are elevated
    (`wFst_YRI_CEU` 0.427, `wFst_JPT_CEU` 0.433) while `wFst_YRI_JPT` is not (0.227),
    which is exactly why there are two Fst peaks and only one PBS peak.
  - `precomputed/` (106 M) holds the saf files and 1D spectra from this run, and
    `realSFS` reads them back correctly, so the shortcut cell works

  **Correction to an earlier note of mine:** I first described the chr5 peak as the
  lactase region. It is not — chr5:33.9 Mb is a **pigmentation** locus. The lactase
  example (chr2 136.5 Mb, via MCM6) belonged to the Shiny bonus, which is not included.
  The exercise text itself never named the gene, and does not now; it sends students to
  the UCSC browser to find it, which is the point of the question.

  **The page's Shiny bonus is not included (dropped by request, 2026-09-16).** Its last
  section launched a genome-wide PBS browser as a Shiny app. It is not runnable on this
  server: the app `require(rCharts)`, a GitHub-only and unmaintained package that is not
  installed, and its export path writes to a home directory that no longer exists. With
  that section gone the notebook reads **entirely** from `data/` — no path anywhere in it
  points outside the data store. The app is untouched at
  `/davidData/albrecht/course/selectionScan/` (487 M, not copied) if it is ever wanted
  back; it would need `rCharts` first.

  The notebook is therefore 50 cells (26 markdown, 12 Bash, 12 R) and ends on the
  genotype-calling bonus, which is a natural close: it is the section that answers "was it
  worth avoiding genotype calls?".

- **#43 `relatedness_animal.ipynb` (2026-09-16).** From the 2026-08-23 kenya2026 notebook,
  into `relatedness_diversity/`. 44 cells, SoS kernelspec (Bash + R + Python quiz cells).
  The pipeline on 111 Greenland reindeer, 2.8 M SNPs, 12 populations: LD pruning from
  precomputed PCAone residuals, KING kinship in plink2, PCA before and after dropping
  relatives, ADMIXTURE K=2, then admixture-aware relatedness with `relateAdmix`, and a
  KING-vs-relateAdmix comparison. Content unchanged; all five quizzes kept.

  - **R4/R11:** the notebook opened with a bare `mkdir -p ~/kenya2026/relatedness` and had
    the data path written out in full in a dozen cells. Now one path cell writing an
    `env.sh` that all 12 bash code cells re-`source`, and the work folder is
    `~/relatedness_animal`.
  - **Data: `data/relatedness/` (1.5 G),** all `cmp`-identical — the Reindeer plink
    fileset, the PCAone residuals, the precomputed ADMIXTURE K=2 output, `population.tsv`
    and the `relateAdmix` binary. The two R plotting helpers went to `data/scripts/` (R14).
  - **A pre-existing broken quiz, fixed.**
    `kenya2026/exercises/Day4/quiz_dataset_summary.json` has a trailing comma before a `}`
    and is **not valid JSON**, so the notebook's first quiz has never loaded for a student.
    Fixed in the copy. Its four answers were then checked against the data and all are
    right: 111 individuals, 2,818,432 SNPs, 12 populations, Kangerlussuaq-Sisimiut largest
    at 21.
  - removed a trailing comma in a `plot_pca_before_after()` call. R tolerates it only
    because the empty argument binds to `sample_col`, which the function never evaluates.

  **Verified by running all 11 bash cells on the real data** with an `ERR` trap: exit 0,
  trap never fired. 2,818,432 SNPs prune to 4,399; KING classifies 4,304 pairs unrelated,
  1,209 second-degree, 456 third-degree, 134 first-degree and 2 duplicate/MZ — which is
  the inflation the notebook's own question ("the values are huge, can you think why?") is
  about; 63 individuals get flagged at the 0.177 cutoff; both PCAs, ADMIXTURE and
  `relateAdmix` produce their outputs. An earlier run exited 1 from leftover state in a
  work directory I had deleted mid-run, not from the notebook.

- **#64 `relatedness_human.ipynb` — new short exercise (2026-09-16).** Requested as the
  human counterpart to #43, built on the data and the analysis from the GWAS-intro
  exercise. 25 cells: 14 markdown, 7 Bash, 4 R. Not from any existing notebook, so nothing
  is superseded.

  The GWAS-intro notebook already runs `plink --genome` as one step among many; this pulls
  that out and makes a short exercise of it: what IBD sharing is, `--genome` on the
  cohort, the Z0/Z1 plot, the PI_HAT distribution, and what to do about it.

  - **written with plain base R by request** — no `plotPlink`, no `newPlotPlink.R`
    dependency. Three plots: a Z1-against-Z0 scatter with the six textbook relationships
    marked, a `PI_HAT` histogram with the degree thresholds, and a sorted missingness plot.
  - **the "assumes no population structure" caveat is stated up front**, in a blockquote,
    with the reason (allele frequencies come from the sample, so shared ancestry looks like
    relatedness) and the note that it is fine for this single-population cohort. The
    closing questions come back to it: what the plot would look like for a mixed sample,
    and how to tell relatedness from structure.
  - **Data: `data/gwas_human/` (58 M)** — `gwa.{bed,bim,fam}` and `pheno3.txt`, extracted
    from `novo23_gwas/GWASex.tar.gz`. The exercise reads the fileset in place instead of
    copying a tarball into the home directory and untarring it, which is what the GWAS
    notebook does. The same folder will serve #46-#48 when they are built.

  **The dataset turned out to have something worth teaching, so the exercise is built
  around it.** Running the analysis first:
  - **there are no close relatives at all**: the largest `PI_HAT` is **0.148**, and only 7
    of 63,190 pairs exceed 0.08. So the honest answer to "should we remove anyone?" is no,
    and that is what makes the cohort usable for the GWAS exercises.
  - **but 54 pairs come back as `nan`**, and they are not random: every one involves the
    15 individuals who are missing **49-55%** of their genotypes. With half the genotypes
    gone in both members of a pair, too few SNPs are left where both are called for the
    method-of-moments estimate to resolve.
  - `--mind 0.2` removes exactly those 15 individuals and the NaN count goes to **zero**,
    while the largest `PI_HAT` stays 0.148 - so the conclusion does not change, which is
    itself the answer to one of the questions.

  Sections 4 and 5 of the notebook are that investigation: count the failures, find the
  individuals, measure missingness with `--missing`, plot it, and redo the estimate on the
  individuals with enough data.

- **#46 `gwas_intro_human.ipynb` (2026-09-16).** From the 2025-07-08 novCourse2024
  notebook, into `gwas/`. 32 cells, SoS kernelspec (Bash + R). The exercise is unchanged:
  an unfiltered logistic GWAS to show what goes wrong, the QQ plot, then QC — `--genome`
  for relatedness, `--cluster --mds-plot` for structure, a filtered re-run — ending on the
  chr4 signal and a LocusZoom-style plot.

  - **the `cp` + `tar -xf` staging is gone.** The original copied `GWASex.tar.gz` into
    `$HOME`, untarred it there, and wrote every plink output straight into the home
    directory. Now a path cell makes `~/gwas_intro_human` and the fileset is read in place
    from `data/gwas_human/`. Four cells dropped: the two commands and their intro
    paragraphs.
  - **R4:** `--bfile data/gwa` -> `$DATA/gwa` in five plink cells; `ls data/`,
    `head data/gwa.fam` and `tail data/gwa.fam` also repointed (they were left behind by
    the staging removal and broke the notebook at its second cell on the first run - caught
    by running it).
  - **R14:** `newPlotPlink.R` sourced from `data/scripts/` instead of
    `/course/novo23/scripts/`.
  - the R cells `setwd()` into the working folder and set `DATA`; the 7 bash code cells
    source `env.sh`.
  - the closing link pointed at `2GWASsumstats.ipynb` in the old course folder; it now
    names `gwas_sumstats_human.ipynb` (#47).
  - `plotPlink` was **kept** here. The base-R instruction was given for the new short
    relatedness exercise (#64); this notebook's plots are the source's own approach and the
    helper already lives in `data/scripts/`.

  **Verified by running it: 7 bash cells and 5 R cells, exit 0, no errors,** and the taught
  result reproduces. The unfiltered GWAS analyses 488,756 variants on 356 people (156
  cases, 200 controls) at a 0.9645 genotyping rate; after QC the top hit is
  **chr4 SNP_A-1978655 at bp 119,229,342, OR 3.91, P 5.1e-06** — the same position the
  LocusZoom cell has hard-coded, so the notebook is internally consistent.

- **What the GWAS-intro data actually contains, found while verifying #46 and #64.** The
  two exercises turn out to be looking at the same 15 individuals from opposite ends, and
  the reason is worth recording:

  - 15 of the 356 individuals are missing **49-55%** of their genotypes. They are the ones
    whose relatedness estimates come back `nan` in #64.
  - **all 15 are cases; none is a control.** Mean missingness is **6.5% in cases against
    1.2% in controls** — a five-fold difference.

  That is textbook differential missingness between case and control batches, and it is
  precisely the artifact #46 is teaching: it is why the unfiltered GWAS produces an
  inflated QQ plot, and why the QC step exists.

  **A wrinkle in #46's QC, left as taught.** Its final QC cell uses `--mind 0.55`, which
  removes exactly **1** of those 15 individuals — while the MDS cell a few cells earlier
  uses `--mind 0.2`, which removes all 15. So the filtered analysis still carries 14 cases
  with ~50% missing data. It reaches the intended answer anyway, because `--geno 0.05`
  then drops 70,801 variants including the differentially-missing ones, leaving 318,914.
  The threshold was **not changed**, since changing it changes the numbers students see;
  flagging it instead. If it was meant to be `0.055`, or to match the `0.2` used above,
  that is a one-character edit.
