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
| 3 | `statistics/stats_binomial.R` | `stat_molbio/binom.R` | 2026-01-16 | no — copied verbatim | none |
| 4 | `statistics/stats_normal.R` | `stat_molbio/normal.R` | 2026-01-16 | no — copied verbatim | none |
| 6+7 | `statistics/haplotype_frequencies.ipynb` | `advBinf/exercises/solution_haplotype_frequencies.ipynb` | 2025-09-12 | yes — rebuilt as one scaffolded notebook; print bug fixed | `advBinf/exercises/haplotype_frequencies.ipynb` (2025-09-12) |
| 12 | `ngs/ngs_intro_human.ipynb` | `chinacourse2026/Day2_Morning_NGSintro_human.ipynb` | 2026-09-14 | yes — quizzes, questions, figure, typos, paths | 5 older copies (see EXERCISES.md) |
| 13 | `ngs/ngs_intro_animal.ipynb` | `kenya2026/exercises/Day1/Kenya2026_NGSintro.ipynb` | 2026-08-17 | yes — quizzes, questions, typos, paths | 6 older copies (see EXERCISES.md) |
| — | `ngs/quiz/*.json` (10 files) | new + `kenya2024/.../quiz{1..4}.json` | 2026-09-16 | new quiz bank | kenya2024 quiz1-4 |

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

- **#13 `stats_binomial.R`, #14 `stats_normal.R`** — no quizzes, no data paths,
  no course branding, no file reads. Copied verbatim.

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
