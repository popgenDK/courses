# current_exercises — build list

Exercises are organised **by theme**, not by species. Human and animal versions
of the same analysis sit side by side in the same theme folder, distinguished by
a `_human` / `_animal` suffix on the filename.

Each entry gives the source it is taken from and **all** earlier copies it is
based on / supersedes. Rules and per-exercise recommendations live in
[`agentGuide/EXERCISE_RULES.md`](../agentGuide/EXERCISE_RULES.md) — read that first.

Dates are the last git commit touching the source file.
`[x]` = created. `[ ]` = still to do.

```
current_exercises/
  linux/  shiny/
  ngs/  genotype_calling_imputation/
  pca/  admixture/  local_ancestry/  gene_flow/
  demography/  selection/  relatedness_diversity/
  gwas/  transcriptomics/
data -> /course/data/
```

**54 exercises** in 54 entries. The numbering runs to 57 because two pairs are each
one exercise listed twice (#6+7 `haplotype_frequencies`, #21+22 `pca_bonus_animal`)
and #17 was removed from scope. Done so far: 16.

---

# linux/

### 1. [x] `intro_linux.md`
- **From:** `summer2025/BriefIntro2Linux.md` (2025-07-23)
- **Supersedes:** `summer2024/BriefIntro2Linux.md` (2024-07-04)

### 2. [x] `intro_bash_linux.md` — **converted from notebook to terminal markdown**
- **From:** `kenya2026/exercises/Day1/IntroToBash.ipynb` (2026-08-18)
- **Supersedes:** `kenya2026/exercises/post_course/day1_morning_bash_linux.ipynb` (2026-08-25)
- **Note:** the notebook form worked badly — `less`, `nano`, `top` and `man` are
  interactive and cannot run in a Jupyter cell, so the original faked them with
  canned output. As a terminal exercise they all work for real. Converted to
  `.md` with 29 fenced `bash` blocks (R15).
- **Data:** `data/popgenmsc26_exercises/linux/Exercises.zip`

---


# shiny/ — R Shiny apps

Shiny server code lives in its own folder rather than under a theme (R13).

### 9. [x] `needleman_wunsch_dna.R`
- **From:** `BSA/NW_DNA.R` (2025-09-02)
- **Supersedes:** none — only copy
- **Note:** DNA version — match/mismatch/gap scoring

### 10. [x] `needleman_wunsch_blosum50.R`
- **From:** `BSA/needleman_wunsch_shiny_app_blosum_50.r` (2025-08-30)
- **Supersedes:** none — only copy
- **Note:** protein version — BLOSUM50 substitution matrix

### 11. [x] `dotplot.R`
- **From:** `BSA/dotplotShiny.R` (2025-08-30)
- **Supersedes:** none — only copy
- **Note:** FASTA supplied through a `fileInput` upload, no fixed data path

---

### 3. [x] `stats_binomial.R`
- **From:** `stat_molbio/binom.R` (2026-01-16)
- **Supersedes:** none — only copy

### 4. [x] `stats_normal.R`
- **From:** `stat_molbio/normal.R` (2026-01-16)
- **Supersedes:** none — only copy

# em_algorithms/

Exercises that build an EM algorithm from scratch.

### 5. [x] `em_algorithm.ipynb`
- **From:** `advBinf/exercises/advBinf_EM_algorithm.ipynb` (2026-09-09)
- **Supersedes:** none — only copy

### 6+7. [x] `haplotype_frequencies.ipynb`
- **From:** `advBinf/exercises/solution_haplotype_frequencies.ipynb` (2025-09-12)
- **Supersedes:** `advBinf/exercises/haplotype_frequencies.ipynb` (2025-09-12) —
  the exercise half, which was too hard to work through
- **One notebook, not two (by request, 2026-09-16).** The old exercise/solution
  pair is replaced by a single scaffolded notebook rebuilt from the solution.
  There is no `haplotype_frequencies_solution.ipynb`.

### 19. [ ] `pca_em_human.ipynb`
- **From:** `advBinf/exercises/advBinf_PCA_EM.ipynb` (2026-09-16)
- **Supersedes:** none — new exercise (EMU / PCAngsd EM algorithms)

### 24. [ ] `admixture_em_human.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture_EM.ipynb` (2026-09-14)
- **Supersedes:** none — new exercise (EM algorithm behind ADMIXTURE)

### 37. [ ] `sfs_model.ipynb`
- **From:** `advBinf/exercises/advBinf_SFSmodel.ipynb` (2026-09-16)
- **Supersedes:** `advBinf/exercises/SFS.md` (2024-09-20)

---

# ngs/

### 12. [x] `ngs_intro_human.ipynb`
- **From:** `chinacourse2026/Day2_Morning_NGSintro_human.ipynb` (2026-09-14)
- **Supersedes:**
  - `chinaCourse2025/Day2_Morning_NGSintro_human.ipynb` (2025-08-04)
  - `summer2025/exercises/Day1_afternoon_NGSintro_human.ipynb` (2025-08-04)
  - `summer2024/exercises/NGSintro.ipynb` (2024-08-19)
  - `bgi23/NGSintro.ipynb` (2023-10-30)
  - `summer2023/IntroNGS/introNGSexercises.md` (2023-08-07)

### 13. [x] `ngs_intro_animal.ipynb`
- **From:** `kenya2026/exercises/Day1/Kenya2026_NGSintro.ipynb` (2026-08-17)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day1_afternoon_ngs_intro.ipynb` (2026-08-25)
  - `summer2025/exercises/Day1_afternoon_NGSintro_animal.ipynb` (2025-08-04)
  - `chinaCourse2025/Day2_Morning_NGSintro_animal.ipynb` (2025-08-04)
  - `kenya2024/exercises/day1_NGSintro/Day1_NGSintroV4.ipynb` (2024-08-07)
  - `kenya2024/exercises/day1_NGSintro/Day1_NGSintroV3.ipynb` (2024-07-26)
  - `kenya2024/exercises/day1_NGSintro/Day1_NGSintroV2.ipynb` (2024-07-15)

### 14. [x] `ngs_inference_human.ipynb`
- **From:** `summer2025/exercises/Day2_NGS_Inference.ipynb` (2025-08-04)
- **Supersedes:**
  - `summer2024/exercises/NGS_inference.ipynb` (2024-08-19)
  - `summer2023/NGSinference/README.md` + `solutions.md` (2023-08-05)

---

# genotype_calling_imputation/

### 15. [x] `genotype_calling_and_imputation_human.ipynb` — **merged exercise**
- **From:** `advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb` (2026-09-09)
- **Supersedes:**
  - `advBinf/exercises/genotype calling and haplotype Imputation.ipynb` (2025-09-12)
  - `advBinf/exercises/SNPandGenotypeCalling.md` (2024-09-13)
- **Paired with:** #16, the imputation half kept separately

### 16. [x] `imputation_human.ipynb` — **the imputation half of #15**
- **From:** `advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb` (2026-09-09)
  for the imputation sections, plus `chinacourse2026/Day2_Afternoon_Genotype_Imputation.ipynb`
  (2026-09-09) for the three sections only it has
- **Scope decided 2026-09-16.** The chinacourse2026 notebook is not the imputation
  half of #15 as first recorded — it is an older, thinner copy of the *whole* same
  exercise (same CEU chr20 2-5 Mb data, same five tools in the same order, same three
  quiz files, same R helpers, same NIPT bonus; 2,271 words and 0 follow-up questions
  against #15's 5,704 and 22). Building it faithfully would have produced two
  near-identical human notebooks, so #16 was built as the imputation-only half it is
  named for: the calling sections belong to #15, and #16 starts from genotype
  likelihoods.
- **Supersedes:**
  - `summer2025/exercises/Day2_Imputation.ipynb` (2025-08-05)
  - `chinaCourse2025/Day2_Afternoon_QUILT_Imputation.ipynb` (2025-07-28)
  - `bgi23/03.QUILT_Imputation_new.ipynb` (2023-11-09)
  - `bgi23/03.QUILT_Imputation.ipynb` (2023-11-09)
  - `bgi23/02.Minimac4_Imputaion.ipynb` (2023-11-09)

---

# pca/

Four exercises — two on **called genotypes** and two on **low depth
sequencing** — plus the shared theory section.

### 58. [x] `pca_mds_and_svd.ipynb` — **new: the shared theory section**
- **From:** extracted from `advBinf/exercises/advBinf_PCA.ipynb` cells 2-25 (2026-09-16)
- **Also appears in:**
  - `kenya2026/exercises/Day3/Kenya2026_PCA.ipynb` (2026-08-17) — same section
  - `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` (2026-09-14) — an MDS section near the end
- **Note:** MDS and PCA worked by hand on the small genotype matrix from the
  slides: `dist`, `cmdscale`, normalising the genotypes, the SVD by hand, the
  covariance matrix, reconstructing the data and the variance explained per PC.
  It was repeated in several notebooks, so it becomes one exercise that the
  others point back to. **Needs no data** — the matrix is typed in.

### 18. [x] `pca_low_depth_human.ipynb`
- **From:** `advBinf/exercises/advBinf_PCA.ipynb` (2026-09-16), from cell 26 on
- **Supersedes:**
  - `summer2025/exercises/Day5_PCA_1.ipynb` (2025-08-06)
  - `chinaCourse2025/Day4_Afternoon_PCA_main.ipynb` (2025-07-28)
  - `advBinf/exercises/PCA.md` (2024-09-17)
  - `summer2024/exercises/summer2024-PCA.ipynb` (2024-08-20)
- **Note:** PCAngsd on genotype likelihoods (`1000G5pops.inputgl.beagle.gz`,
  `eu1000g.small.beagle.gz`). The MDS/by-hand opening moves to #58.

### 59. [x] `pca_low_depth_selection_human.ipynb`
- **From:** `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` (2026-09-14), the
  **PC-based selection** half (from the `# PC-based selection` heading)
- **Supersedes:** none — only copy
- **Note:** `pcangsd --selection` on `eu1000g.small.beagle.gz`, then mapping the
  selection statistic back to chromosome and position.
- **Data:** `data/chinacourse2026_shared/`

### 22. [x] `pca_called_genotypes_human.ipynb`
- **From:** `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` (2026-09-14), the
  **first** half, up to the `# PC-based selection` heading
- **Supersedes:** none — only copy
- **Note:** `PCAone` on the LD-pruned called genotypes
  `human_autosomes_12pp_pcaoneLD02`, plotted by population and super-population
  against the admixture proportions.
- **Data:** `data/chinacourse2026_shared/`

### 20. [x] `pca_called_genotypes_animal.ipynb`
- **From:** `advBinf/exercises/advBinf_PCA_bonus.ipynb` (2026-09-16)
- **Supersedes:**
  - `kenya2026/exercises/Day3/Kenya2026_PCA.ipynb` (2026-08-17)
  - `kenya2026/exercises/post_course/day4_morning_pca.ipynb` (2026-08-25)
  - `summer2025/exercises/Day5_PCA_2.Call_genotype.ipynb` (2025-08-06)
  - `chinaCourse2025/Day4_Afternoon_PCA_bonus.ipynb` (2025-07-28)
  - `summer2024/exercises/summer2024-PCA-CalledGenotypes.ipynb` (2024-08-20)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_PCA-V2.ipynb` (2024-08-09)
- **Note:** `PCAone` on the wildebeest plink files, plotted against the admixture
  proportions, then LD pruning and an IBS tree.

  **This was previously listed twice, as #20 and #21.** The "bonus" notebook and
  the kenya notebook are the same exercise: 20 of the bonus notebook's 34 cells
  match the kenya one's wildebeest section, 12 of them identically. The
  differences are that the kenya version opens with the MDS/by-hand section
  (which moves to #58) and takes its LD-pruned file ready-made from the admixture
  exercise, whereas the bonus version **computes the LD pruning itself** with
  `PCAone --ld`, adjusting the LD measure for population structure, and then
  re-runs the PCA on the pruned data.

  Build from the bonus notebook, which is the newer and more complete of the two:
  it already covers the kenya content and adds the LD-pruning section. Carry over
  anything the kenya version phrases better.
- **Data:** wildebeest plink files (`blue_wildebeest_thin`, `blue_wildebeest_noLD`)
---

# admixture/

Four exercises — two on **called genotypes** and two on **low depth sequencing** —
mirroring the PCA folder, plus local ancestry and the gene-flow material.

### 23. [x] `admixture_low_depth_human.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture.ipynb` (2026-09-14), cells 5-70
- **Supersedes:**
  - `summer2025/exercises/Day3_Morning_Admixture.ipynb` (2025-08-05)
  - `advBinf/exercises/admixture.md` (2024-09-16)
  - `summer2024/exercises/admixExercise_popgen24.ipynb` (2024-08-19)
  - `bgi23/Admixture.ipynb` (2023-11-02)
  - `summer2023/InfererPopStructure/admixExercise_popgen23.ipynb` (2023-08-08)
- **Note:** NGSadmix on genotype likelihoods, then evalAdmix and the choice of K.

### 60. [x] `admixture_reference_panel_human.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture.ipynb` (2026-09-14), cells 71-93
- **Supersedes:** none — only copy
- **Note:** fastNGSadmix: the ancestry of a **single** individual against a fixed
  reference panel of 7 populations. Split out of #23, which covered two methods.

### 61. [x] `admixture_called_genotypes_human.ipynb`
- **From:** `chinacourse2026/Day4_admix_eval_LAI.ipynb` (2026-07-30), cells 2-54
- **Supersedes:** `chinaCourse2025/Day4_Morning_admixture_genotype.ipynb` — **no**, see #25
- **Note:** ADMIXTURE on LD-pruned called genotypes, convergence across seeds, and
  evalAdmix. The local-ancestry tail of the source notebook (cells 55+) is a separate
  topic and stays with #28/#29.

### 25. [x] `admixture_called_genotypes_animal.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture_bonus.ipynb` (2026-09-14)
- **Supersedes:**
  - `kenya2026/exercises/Day3/Exercises_Admixture_Kenya26_WoA.ipynb` (2026-08-21)
  - `kenya2026/exercises/post_course/day3_morning_admixture.ipynb` (2026-08-25)
  - `chinaCourse2025/Day4_Morning_admixture_genotype.ipynb` (2025-07-28)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_AdmixtureV2.ipynb` (2024-08-08)
  - `kenya2024/exercises/day3_PopulationStructure/Admixture.ipynb` (2024-08-07)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_Admixture.ipynb` (2024-07-28)
- **Note:** **#25 and #26 were the same exercise** — the kenya notebook and
  `advBinf_admixture_bonus` are 45 of 45 cells similar, 43 identical. Built from the
  advBinf one as the newer. `chinaCourse2025/Day4_Morning_admixture_genotype.ipynb` was
  wrongly listed under #23: it is **wildebeest, called genotypes**, only 2 of its 48
  cells resemble #23, and it is 72% similar to this exercise. Moved here.

### 28. [x] `local_ancestry_flare_mosaic_human.ipynb`
- **From:** `summer2025/exercises/Day4_Morning_LocalAncestry.ipynb` (2025-08-06)
- **Supersedes:** none — only copy
- **Note:** FLARE and MOSAIC on simulated admixed genomes, comparing admixture 20 vs 200
  generations ago against the known true tracts.

### 29. [x] `local_ancestry_hapla_human.ipynb`
- **From:** `advBinf/exercises/Hapla_LAI_exercise.ipynb` (2025-10-07)
- **Supersedes:** `chinacourse2026/Day4_admix_eval_LAI.ipynb` cells 55-68 (2026-07-30) —
  a 14-cell "short look" at hapla/fatash using precomputed files, which this covers fully
- **Note:** hapla cluster / admix / fatash on data simulated with `msprime`, then the same
  pipeline on a real cattle dataset (314 individuals, BosTau9 chr25).

**These three local-ancestry sources share nothing.** Compared cell by cell, the
summer2025 notebook, the advBinf hapla notebook and the China tail have **zero** cells in
common — they use different tools (FLARE/MOSAIC, hapla/fatash, and a short hapla demo).
The plan previously had #28 superseding the summer2025 notebook, which was wrong: they
are different exercises, not versions of each other.

**Retired slot: #27 `population_structure_ii_human`.** Not an admixture exercise and
not built. `bgi23/BGI2023-populationStructureII.ipynb` has three sections — *Simple
example of PCA and MDS*, *PCA for low depth sequencing using PCAngsd* and *PCAngsd and
selection* — and **53 of its 69 cells are identical to `advBinf_PCA.ipynb`**. It is a
2023 ancestor of the PCA trio and is superseded by #58, #18 and #59 in `pca/`. It was
filed under admixture because of its "population structure" title. The number is kept so
the other exercises do not shift.

### 30. [x] `f_stats_human.ipynb`
- **From:** `summer2025/exercises/Day3_f_stats.ipynb` (2025-08-05)
- **Data: MISSED BY THE ORIGINAL SURVEY.** `/course/popgen25/dstats` (108 M) was not in
  the data-copy list — the notebook reaches it through a `~/data_folder` symlink, so no
  literal `/course/...` path appears in it. Now copied to `data/f_stats/`.
- **Supersedes:**
  - `summer2024/exercises/f_stats.ipynb` (2024-08-21)
  - `summer2023/DfFstats/popgen23_f_stats.ipynb` (2023-08-09)

### 31. [x] `gene_flow_dstat_animal.ipynb`
- **From:** `kenya2026/exercises/Day3/Geneflow&Dstat.ipynb` (2026-08-16)
- **Supersedes:** `kenya2026/exercises/post_course/day3_afternoon_gene_flow_dstat.ipynb` (2026-08-25)
- **Data:** the wildebeest f2 set came from `/davidData/users/thomas/workshop`, an
  instructor home directory (77 G). Only the one 1.5 M subfolder is needed, so it and
  the two simulated-data files were copied to `data/geneflow/` (3.3 M total).

### 32. [x] `admixture_graphs_human.ipynb`
- **From:** `summer2023/DfFstats/popgen23.Admixture_Graphs_Tutorial.ipynb` (2023-08-10)
- **Supersedes:** none — only copy
- **Note:** `qpgraph` (fit a graph you specify) and `treemix` (estimate one), on
  precomputed $F_2$ statistics for 33 world populations. The exercise keeps its worked
  solutions at the end, clearly marked.
- **Data:** `data/current_data/admixture_graphs/` — 131 M, the `fdata`, `software` and
  `treemix` parts of `/course/popgen23/ben/fstats_tutorial`. The 1.9 G raw genotype
  directory is not needed: the $F_2$ statistics are precomputed.

### 33. [x] `chromopainter_finestructure_human.ipynb`
- **From:** `summer2024/exercises/ChromoPainterFineSTRUCTUREPractical.ipynb` (2024-08-21)
- **Supersedes:** none — only copy
- **Companion:** `summer2024/exercises/CopenhagenPopgenWorkshop2024_ChromoPainterFineSTRUCTUREPracticalSOLN.pdf`
- **Note:** ChromoPainter paints each genome as a mosaic of the others; fineSTRUCTURE
  clusters from the chunk counts by MCMC; then GLOBETROTTER and SOURCEFIND infer the
  ancestry of the admixed target. Compiles its software from source.
- **Data:** `data/current_data/chromopainter/` — a single 15 M tarball, unpacked by the
  setup cell.
- **Leads to** #34, which uses the ChromoPainter output from part 2.

### 34. [x] `dating_admixture_human.ipynb`
- **From:** `summer2024/exercises/DatingAdmixture.ipynb` (2024-08-21)
- **Supersedes:** none — only copy
- **Companion:** `summer2024/exercises/CopenhagenPopgenWorkshop2024_DatingAdmixturePracticalSOLN.pdf`
- **Note:** four tools on the same simulated event — ALDER, MALDER, fastGLOBETROTTER and
  MOSAIC to date it, then AdaptMix to test for selection on the admixed ancestry. The
  exercise **compiles the software from source**, so several cells take minutes.
- **Data:** `data/current_data/dating_admixture/` — a single 41 M tarball holding the data
  and the tool sources. Unpacked by the setup cell.
- **Continues from** #33 ChromoPainter/fineSTRUCTURE, which produces the input for part 2.

---

# demography/

### 35. [x] `coalescence.ipynb`
- **From:** `kenya2026/exercises/Day2/Coalescence_short_WoA.ipynb` (2026-08-15)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day2_morning_coalescence.ipynb` (2026-08-25)
  - `summer2025/exercises/Day1_morning_CoalTutorial.ipynb` (2025-08-03)

### 36. [x] `wright_fisher.ipynb`
- **From:** `summer2025/exercises/Day1_morning_WrightFisherTutorial.ipynb` (2025-08-03)
- **Supersedes:** none — only copy
- **Helper:** both #35 and #36 `source()` `simulateWF.R`, now at `data/scripts/`
  per R14. It was only inside the 1.8 G `popgen25_software/` tree before.

### 38. [x] `sfs_animal.ipynb`
- **From:** `kenya2026/exercises/Day2/SFS_WoA.ipynb` (2026-08-15)
- **Supersedes:** `kenya2026/exercises/post_course/day2_morning_sfs.ipynb` (2026-08-25)
- **Data:** the 5 VCFs it reads were copied from `kenya2026_harvi/sfs/inputdata/` to a
  cleaned `data/sfs/` folder (291 M, all `cmp`-identical).

### 39. [x] `psmc_demography_animal.ipynb`
- **From:** `kenya2026/exercises/Day2/psmc_kenya2026.ipynb` (2026-08-19)
- **Supersedes:** `kenya2026/exercises/post_course/day2_afternoon_psmc.ipynb` (2026-08-25)
- **Data:** both PSMC exercises now read one shared folder, `data/psmc/`
  (284 M: `data/` simulated + 1000 Genomes, `animal/` wildebeest, `images/`,
  `scripts/`, `software/` with the psmc binary and its utils). Trimmed to the
  files the two notebooks actually read; 10 unreferenced files (211 M) removed,
  all still available in `popgen25_demography/data/`.
- **Note:** PSMC on simulated data, then on **real wildebeest samples**. Same structure as
  #62; the two differ in which real individuals are analysed.

### 62. [x] `psmc_demography_human.ipynb`
- **From:** `summer2025/exercises/Day5_demography.ipynb` (2025-08-07)
- **Supersedes:**
  - `summer2024/exercises/summer2024-PSMC_tutorial_2024.ipynb` (2024-08-23)
  - `bgi23/PSMC_tutorial.ipynb` (2023-11-07)
  - `summer2023/DemographyInference/PSMC_tutorial.ipynb` (2023-08-10)
- **Note:** PSMC on simulated data, then on **two 1000 Genomes individuals** — NA12718, a
  CEU female of northern European ancestry, and NA19471, a Luhya female from Kenya —
  comparing how their inferred effective population sizes differ. Ends with a wildebeest
  bonus, which #39 covers properly.
- **There are two PSMC exercises, human and animal.** The plan had only the animal one,
  with `summer2025/Day5_demography.ipynb` listed as superseded by it. That was wrong: the
  two share 58% of their cells because they are the **same exercise on different species**,
  not because one is an older copy. The older notebooks really are superseded — summer2024
  is 78% similar to summer2025, summer2023 85% similar to summer2024, and bgi23 and
  summer2023 are human-only ancestors.

---

# selection/

**Five exercises**, differing in *which species*, in whether the genotypes are **called**
or the data are **low depth genotype likelihoods**, and in whether the method is
**frequency-based** or **haplotype-based**.

`summer2024/exercises/SelectionScans.ipynb` turned out to hold **three** exercises, and is
superseded piece by piece: cells 0-18 by #63, cells 19-35 by #66, cells 36-48 by #65.
Cells 49-50 are an empty "EDAR" stub. `summer2023/selectionScan/README.md` is a README,
not an exercise.

### 63. [x] `sfs_fst_pbs_human.ipynb`
- **From:** written for this set
- **Supersedes:**
  - `summer2024/exercises/SelectionScans.ipynb` cells 0-18 (2024-08-22) — the
    frequency-statistics half only. Its haplotype half is #66 and its PBS scan is #65
  - `kenya2026/exercises/Day4/SelectionScans_22nd.ipynb` cells 0-24 (2026-08-22) — the
    frequency-statistics half of that notebook, 49% similar to summer2024. Cells 25-37,
    the genome-wide PBS scan, are **not** covered by this exercise — see #65
- **Note:** SFS, Fst and PBS on human data with ANGSD. **Low depth**, genotype likelihoods.

### 66. [x] `selection_haplotype_human.ipynb`
- **From:** `summer2024/exercises/SelectionScans.ipynb` cells 19-35 (2024-08-22),
  "Exercise II: Haplotype-based methods"
- **Supersedes:** none — only copy
- **Note:** extended haplotype homozygosity — `selscan --ihs` and `--xpehh` on the lactase
  region, with normalisation within allele-frequency bins. The **only** exercise in the set
  covering haplotype-based selection; the others are all frequency-based.
- **Correction:** summer2024 was recorded as superseded by #63. #63 replaces its
  frequency-statistics half only. This haplotype half was covered by nothing.
- **Data:** `data/current_data/selection/haplotype/` (135 M) — phased CEU/YRI/CHB VCFs
  around LCT, a genetic map in cM, the selscan binary and the precomputed output.

### 65. [x] `selection_pbs_scan_human.ipynb`
- **From:** `kenya2026/exercises/Day4/SelectionScans_22nd.ipynb` cells 25-37 (2026-08-22),
  "Exercise II: Whole-genome PBS with 1000 Genomes"
- **Supersedes:** `summer2024/exercises/SelectionScans.ipynb` cells 36-48 (2024-08-22) —
  the same scan, older
- **Note:** genome-wide PBS on precomputed 50 kb windows for NAT, CHB, CEU and YRI:
  Manhattan plot, zoom into peaks, identify the gene, and check the lactase region as a
  positive control.
- **Correction:** #63 was originally recorded as superseding *all* of cells 0-37. It does
  replace the frequency-statistics half, but **not** this genome-wide scan, which nothing
  else covers. Split out as its own exercise.
- **Data:** `data/current_data/selection/pbs_scan/` (7.3 G).

### 40. [x] `selection_scans_animal.ipynb`
- **From:** `kenya2026/exercises/Day4/SelectionScans_22nd.ipynb` cells 38-50 (2026-08-22)
- **Supersedes:** `kenya2026/exercises/post_course/day4_afternoon_selection_scans.ipynb` (2026-08-25)
- **Note:** the wildebeest half of the kenya notebook — Hudson Fst and PBS on 24
  individuals in three groups. **Called genotypes** (VCF). Black wildebeest is a separate
  species and plays the outgroup role.

### 41. [x] `selection_maize.ipynb`
- **From:** `summer2025/exercises/Day4_SelectionPopGen2025.ipynb` (2025-08-06)
- **Supersedes:** none — only copy
- **Note:** **maize**, not an animal despite the old `_animal` name. Tajima's D and PBS
  from **low depth genotype likelihoods** with ANGSD. **Zero** cells in common with any
  other selection exercise.

**Also a selection scan:** [`pca/pca_low_depth_selection_human.ipynb`](#) (#59) scans along
a principal component rather than using Fst. It lives in `pca/` because it continues
directly from the PCA exercise.

**The species split, confirmed by cell comparison:** summer2024 is human; summer2025 is
**maize**; the kenya notebook is a **mix** whose human half is 49% similar to summer2024
and whose wildebeest half shares nothing with it. The old plan had #41 (`_animal`, maize)
superseding summer2024 (human), which share **0%** of their cells.

---

# relatedness_diversity/

### 42. [x] `fst_animal.ipynb`
- **From:** `kenya2026/exercises/post_course/day4_morning_fst.ipynb` (2026-08-25)
- **Supersedes:**
  - `kenya2026/exercises/Day4/Fst_Kenya2026.ipynb` (2026-08-23)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_Fst_RH.ipynb` (2024-08-09)
- **Note:** pairwise $F_{ST}$ two ways — `plink2` on **called genotypes** (95 wildebeest in
  9 groups, including black wildebeest as an outgroup species), then SAF-based estimation
  on **genotype likelihoods** (3 Greenland reindeer populations). The post-course version
  was chosen over the taught one because it reads a far smaller dataset (R7).
- **Data:** `data/current_data/fst/` — 8.0 G, of which 6.8 G is the reindeer SAF files.
  Replaced the post-course download machinery (`wget` of a zip into
  `$KENYA2026_WORK_DIR`) with a standard setup cell.

### 64. [x] `relatedness_human.ipynb` — **new, short exercise**
- **From:** written new (2026-09-16), by request, using the data and the `plink --genome`
  analysis from #46 `gwas_intro_human.ipynb`.
- **Supersedes:** none — new exercise, the human counterpart of #43.
- **Why it exists:** a GWAS treats individuals as independent. Relatives share long
  stretches of genome and so break that assumption, inflating the test statistics — which
  is why you check relatedness *before* running the association. **Leads into** `gwas/`.
- **Note:** deliberately short. IBD sharing with `plink --genome`, a Z1-vs-Z0 plot and a
  PI_HAT histogram in **base R** (no `plotPlink`), then the missingness problem the data
  turns out to have. States up front that this estimator assumes one homogeneous
  population, which holds for this cohort.
- **Data:** `data/gwas_human/` (58 M) — `gwa.{bed,bim,fam}` + `pheno3.txt`, extracted from
  `novo23_gwas/GWASex.tar.gz`. The same folder will serve #46-#48.
- **What the data shows:** no close relatives at all (largest PI_HAT 0.148, 7 of 63,190
  pairs above 0.08), but 54 pairs come back `nan` — all involving the 15 individuals
  missing 49-55% of their genotypes. `--mind 0.2` removes exactly those 15 and the NaNs
  go to zero.

### 43. [x] `relatedness_animal.ipynb`
- **From:** `kenya2026/exercises/Day5/Related.ipynb` (2026-08-23)
- **Supersedes:** `kenya2026/exercises/post_course/day5_morning_relatedness.ipynb` (2026-08-25)
- **Paired with:** #44, the merged Related&Fst form
- **Data:** `data/relatedness/` (1.5 G) — Reindeer plink fileset, the PCAone residuals,
  the precomputed ADMIXTURE K=2 output and the `relateAdmix` binary. The two R plotting
  helpers went to `data/scripts/` per R14.
- **Fixed a broken quiz:** `kenya2026/exercises/Day4/quiz_dataset_summary.json` has a
  trailing comma and is **not valid JSON**, so the notebook's first quiz has never
  loaded. Fixed in the copy.

**Retired slot: #44 `relatedness_and_fst_animal`.** Removed by request, 2026-09-16.
`kenya2026/exercises/Day4/Related&Fst.ipynb` (2026-08-22) combined the two analyses in one
notebook; they are kept separately as #42 Fst and #43 relatedness, which is how the newer
kenya2026 notebooks split them. The source stays untouched in `kenya2026/`.

### 45. [x] `heterozygosity_roh_animal.ipynb` — **post-course small-dataset version**
- **From:** `kenya2026/exercises/post_course/day5_morning_heterozygosity_roh.ipynb` (2026-08-25)
- **Supersedes:**
  - `kenya2026/exercises/Day5/Day5_GeneticDiversity.ipynb` (2026-08-23)
  - `kenya2024/exercises/day2/Day2_Inbreeding_ROH.ipynb` (2024-08-08)
  - `kenya2024/exercises/day2/Inbreeding_ROH.ipynb` (2024-08-07)
  - `kenya2024/exercises/Day1_GeneticDiversity/Day1_GeneticDiversity.ipynb` (2024-08-06)

---

# gwas/

### 46. [x] `gwas_intro_human.ipynb`
- **From:** `novCourse2024/1GWASIntro.ipynb` (2025-07-08)
- **Supersedes:** `bgi23/04.GWASintro_2023_SAIGE.ipynb` (2023-11-09)
- **Data:** `data/gwas_human/` (58 M), shared with #64 and the other GWAS exercises. The
  `cp` + `tar -xf` staging into `$HOME` was dropped; the fileset is read in place and all
  output goes to `~/gwas_intro_human`.
- **Note on its QC step:** the final QC uses `--mind 0.55`, which removes only **1** of the
  15 individuals missing ~50% of their genotypes (the MDS cell above it uses `--mind 0.2`,
  which removes all 15). Left as taught — `--geno 0.05` drops the affected SNPs anyway and
  the intended chr4 result still comes out — but see MANIFEST.

### 47. [x] `gwas_sumstats_human.ipynb`
- **From:** `novCourse2024/2GWASsumstats.ipynb` (2025-07-08)
- **Supersedes:** none — only copy
- **Data:** the 796 M Mahajan T2D summary statistics, read in place from
  `data/novo23_gwas/sumstats/`. Deliberately **not** given a cleaned per-exercise folder:
  it is a single large file already inside `data/`, and copying it again to rename it
  would waste 796 M. The original copied it into the student's home directory.

### 48. [x] `gwas_analysis_human.ipynb` — **⚠️ data not on this server**
- **From:** `chinaCourse2025/Day3_GWAS_Analysis_2025_Morning.ipynb` (2025-07-26)
- **Supersedes:** none — only copy
- **Note:** six exercises A-F: a first GWAS, the QQ plot, QC, PCA and Tracy-Widom, a
  linear mixed model with `regenie` run with and without the top 20 PCs, gene-based
  testing with SKAT/ACAT, and GWAS power. 122 cells.
- **⚠️ The data is missing.** The notebook was written for a course machine where each
  student had `European_1w` under `/home/student/<user>/GWAS/data/`. That folder does not
  exist here and the dataset is in none of the course archives — the `chinacourse2025` zip
  for this day holds only the notebook and the lecture PDFs. **Added as-is by request**, so
  the material is not lost.
- **What was still done:** all 43 data references and the conda activation were collected
  into the setup cell, so pointing `DATA` at the files is the only edit needed once they
  are found. Kernel names were normalised — the source declared `R4.4` and
  `Python 3.12 (py312)`, which do not exist on this server, so no cell would have run.
- **Overlap to check later:** Exercise E (gene-based testing) may duplicate #49
  `gene_based_testing_human`, which comes from the afternoon notebook of the same day.

### 49. [x] `gene_based_testing_human.ipynb` — **⚠️ data not on this server**
- **From:** `chinaCourse2025/Day3_Gene_Based_Testing_2025_Afternoon.ipynb` (2025-07-26)
- **Supersedes:** none — only copy
- **Note:** rare-variant gene-based testing with `regenie` — annotation, set list and mask
  files, then SKAT and ACAT — on **whole-exome data for a binary trait**
  (Charcot-Marie-Tooth disease), with a Firth-corrected null model. Ends with power
  calculations.
- **Overlap with #48, checked:** 45 of its 54 source cells repeat Exercises D, E and F of
  the morning notebook, 38 of them identically. What differs is the **data and the trait**:
  #48 runs on array genotypes with a quantitative phenotype, this one on exome data with a
  binary disease trait, which is why the null model needs Firth correction and why
  gene-based tests are used at all. Kept as a separate exercise on that basis; the header
  says so and points back to #48.
- **⚠️ The data is missing,** same as #48 — the notebook expects
  `/home/student/<user>/GWAS/data/`, which does not exist here. Added as-is by request,
  with every path collected into the setup cell.

### 50. [x] `wes_famdiab_human.ipynb` — **⚠️ data not on this server, added verbatim**
- **From:** `novCourse2024/3WESfamdiab.ipynb` (2025-07-08)
- **Supersedes:** none — only copy
- **Note:** a trio with MODY. Read the VCF by eye to work out which sample is the
  mother and which the daughter, then annotate the variants with the online
  **Ensembl VEP** tool and pick the causal one. **Run partly in a web browser.**
- **Missing:** `hnf1a.vcf` and the pedigree/VEP screenshots (`images/`), which
  were distributed with the original course and are in no archive.
- **⚠️ Data not on this server, so the notebook is a verbatim copy.** By request
  (2026-09-16), exercises whose data is missing are **not edited at all** — no
  header, no path centralisation, no added questions or quizzes. Whether they run
  does not matter; the point is that the material is not lost. Fix them when the
  data turns up.

### 51. [x] `wes_fh_human.ipynb`
- **From:** `novCourse2024/4WESfh.ipynb` (2025-07-08)
- **Supersedes:** none — only copy
- **Note:** rare recessive familial hypercholesterolemia in a family of four with
  second-cousin parents. Filter an annotated exome on gnomAD frequency and
  predicted consequence; the frequency cutoff is derived from the prevalence.
- **Data: COPIED** to `data/current_data/wes/ex02.wes.rds` (from `/course/novo23/wes/`)
- **Note:** the pedigree images (`images/1.jpg`, `images/2.jpg`) are in no archive,
  so the family is described in the header text instead.
- **Converted from `ir` to SoS** with per-cell kernels so the quiz can run (R17).

### 52. [x] `wes_diabetes_human.ipynb`
- **From:** `novCourse2024/5WESdiab.ipynb` (2025-07-08)
- **Supersedes:** none — only copy
- **Note:** four unrelated diabetes patients, each with a different damaging
  variant in a different known gene — genetic heterogeneity, and why a GWAS
  would find none of them.
- **Data: COPIED** to `data/current_data/wes/ex03.wes.rds` (from `/course/novo23/wes/`)
- **Converted from `ir` to SoS** with per-cell kernels so the quiz can run (R17).

### 53. [x] `heritability_ldscore_human.ipynb`
- **From:** `chinacourse2026/Day5_Morning_heribilty_and_ldscore.ipynb` (2026-08-07)
- **Supersedes:**
  - `chinaCourse2025/Day5_heritability_exercise.ipynb` (2025-07-26)
  - `chinaCourse2025/Day5_Afternoon_Genetic_correlation_Partitioned_Heritability.ipynb` (2025-07-26)
- **Note:** heritability estimated twice — GCTA REML on a genetic relationship
  matrix from family genotypes, and LD score regression on published summary
  statistics. The comparison is the point: the LDSC **intercept** separates
  confounding from polygenicity, which λ_GC cannot do.
- **Data: COPIED** to `data/current_data/heritability_ldscore/` — `quantfam.zip`
  and the Biobank Japan HDL/LDL summary statistics with `eas_ldscores` and
  `w_hm3.snplist` (from `/course/chinacourse2026/shared/data/`).
- **⚠️ Software gap:** `gcta64` is installed, but the **LDSC conda environment is
  not** (`/home/jonas/miniconda3` no longer exists). The two `munge_sumstats.py`
  cells will fail. Their output, the `.sumstats.gz` files, ships with the data, so
  everything downstream still runs; the notebook says so at that point.

### 54. [x] `prs_height_human.ipynb` — **⚠️ data not on this server, added verbatim**
- **From:** `chinaCourse2025/Day6_Morning1_PRS_height_pipeline.ipynb` (2025-07-26)
- **Supersedes:** none — only copy
- **Note:** three PGS Catalog scores (height, IGF-1, birth weight) harmonised
  against array genotypes and computed with `plink --score`, then compared with
  measured height and BMI.
- **Missing:** the same `European_1w` dataset as #48, plus the `PRS_data/PGS*.txt.gz`
  weights.
- **⚠️ Data not on this server, so the notebook is a verbatim copy.** By request
  (2026-09-16), exercises whose data is missing are **not edited at all** — no
  header, no path centralisation, no added questions or quizzes. Whether they run
  does not matter; the point is that the material is not lost. Fix them when the
  data turns up.

### 55. [x] `mendelian_randomization_human.ipynb`
- **From:** `chinaCourse2025/Day6_Morning2_MR.ipynb` (2025-07-26)
- **Supersedes:**
  - `bgi23/MR-exercise.ipynb` (2023-11-10)
  - `bgi23/MR.real_data.exercise.ipynb` (2023-11-10)
- **Note:** two-sample MR of BMI on coronary heart disease with `TwoSampleMR` —
  harmonisation, the MR estimators, the F-statistic, pleiotropy and heterogeneity
  tests, forest/funnel/leave-one-out plots, and the IVW estimate written out by
  hand at the end.
- **Data: COPIED** to `data/current_data/mendelian_randomization/` — the IEU
  OpenGWAS RDS files `ieu-a-2.rds` (BMI, exposure) and `ieu-a-7-out.rds` (CHD,
  outcome), **found in the bgi23 archive** at `bgi23/malthe/friday/`, not where the
  chinaCourse2025 notebook expected them. The `extract_instruments()` call that
  would fetch them over the network is kept as a comment.

### 56. [x] `proteomics_mr_human.ipynb` — **⚠️ data not on this server, added verbatim**
- **From:** `chinaCourse2025/Day6_Afternoon_Proteomics_MR.ipynb` (2025-07-26)
- **Supersedes:** none — only copy
- **Note:** three exercises — proteome-wide summary-data MR of plasma proteins on
  ischaemic stroke (cis+trans, then cis-pQTLs only), **colocalization** to tell a
  shared causal variant from two neighbouring ones, and individual-level
  proteomics MR in a UK Biobank stroke cohort.
- **Missing:** the whole `/home/student/Proteomics_MR/` tree, which is in no archive.
- **⚠️ Data not on this server, so the notebook is a verbatim copy.** By request
  (2026-09-16), exercises whose data is missing are **not edited at all** — no
  header, no path centralisation, no added questions or quizzes. Whether they run
  does not matter; the point is that the material is not lost. Fix them when the
  data turns up.

---

# Helper libraries — in the data folder, not here

R function libraries that exercises `source()` are **not exercises**. They live
in the data folder at `data/scripts/`, not under `current_exercises/` (R14).

| Done | File | Source | Date | Sourced by |
|---|---|---|---|---|
| [x] | `data/scripts/admixFun.R` | `chinaCourse2025/assets/admixFun.R` | 2025-07-26 | #25, #26 |
| [x] | `data/scripts/newPlotPlink.R` | `summer2023/InfererPopStructure/newPlotPlink.R` | 2025-07-16 | #27, #46, #47, #53 |
| [x] | `data/scripts/online.R` | `kenya2024/online.R` | 2024-08-07 | #18, #53, #54 |

All three were copied verbatim on 2026-09-16. Exercises should source them as
`data/scripts/<name>.R`.

Paths still to fix inside them:
- `newPlotPlink.R` hardcodes `/davidData/data/course/scripts/geneticMap/{hg19,hg38}`
  in four places. That tree is now copied to `data/geneticMap/`, so the four
  references need repointing.
- `online.R` writes to `~/public/albrecht/temp.pdf` and `~/public/albrecht/temp.txt`
  inside its `pdff()` helper — a personal path that should become a temp file.
- `admixFun.R` has no paths; copied clean.

# Data

Exercise data lives in **`data/current_data/`** (= `/course/data/current_data/`),
one folder per exercise or per source. `data` is a symlink in the repo root
pointing at `/course/data/`.

```
data/
  current_data/     everything the exercises read
    NGSintro/         animal/ human/ software/
    NGSinference/
    PCA/              animal/ human_called/ human_lowdepth/
    imputation/       bams/ resources/ software/ vcfs/
    geneticMap/       reference data for the locus zoom plot
    scripts/          helper R libraries that exercises source()
    BSA/
    <course folders>  raw material for the exercises not yet built
  1000G/            predates this work
```

# Data consolidation (source material)

Real `cp` copies. **Completed 2026-09-16**: 24 directories, 24 OK, 0 warnings,
0 failures, 0 skipped. Every copy verified identical file counts before and
after. Plus `BSA/` and `geneticMap/` copied separately, and `1000G/` which
predates this work.

| Done | `data/` target | Source | Size |
|---|---|---|---|
| [x] | `1000G/` | already present | 73 G |
| [x] | `BSA/` | `/davidData/data/BSA` | 1.2 G |
| [x] | `NGSintro/` | curated for exercises #12/#13 (see below) | 465 M |
| [x] | `NGSinference/` | `/course/popgen25/NGSInference/data` — for #14 | 2.1 G |
| [x] | `PCA/` | `/course/popgen25/pca` + china plink files — for #58/#18/#59/#22/#20 | 461 M |
| [x] | `geneticMap/` | `/course/scripts/geneticMap` | 480 M |
| [x] | `thomas_workshop/` | `/davidData/users/thomas/workshop` | 77 G |
| [x] | `kenyaWorkshop_anders/` | `/course/kenyaWorkshop/anders` | 76 G |
| [x] | `thomas_hmmadvbinf/` | `/davidData/users/thomas/hmmadvbinf` | 38 G |
| [x] | `bgi23_quan/` | `/course/bgi23/quan` | 27 G |
| [x] | `kenya2026/` | `/course/kenya2026/data` | 23 G |
| [x] | `popgenmsc26_exercises/` | `/course/popgenmsc26/exercises` | 17 G |
| [x] | `popgen24_anders/` | `/course/popgen24/anders` | 13 G |
| [x] | `chinacourse2026_shared/` | `/course/chinacourse2026/shared` | 6.7 G |
| [x] | `popgen25_imputation/` | `/course/popgen25/Imputation` | 3.3 G |
| [x] | `kenya2026_nuno/` | `/course/kenya2026/nuno` | 2.6 G |
| [x] | `kenya2026_anders/` | `/course/kenya2026/anders` | 2.4 G |
| [x] | `popgen23_ben/` | `/course/popgen23/ben` | 2.1 G |
| [x] | `popgen25_software/` | `/course/popgen25/software` | 1.8 G |
| [x] | `novo23_gwas/` | `/course/novo23/gwas` | 783 M |
| [x] | `advBinf_admixture_human/` | `/course/advBinf/admixtureHuman` | 746 M |
| [x] | `kenya2026_harvi/` | `/course/kenya2026/harvi` | 293 M |
| [x] | `popgen25_demography/` | `/course/popgen25/demography` | 242 M |
| [x] | `kenya2026_ida/` | `/course/kenya2026/ida` | 176 M |
| [x] | `advBinf_admixture_bonus/` | `/course/advBinf/admixtureBonus` | 173 M |
| [x] | `popgen24_cindy/` | `/course/popgen24/cindy` | 135 M |
| [x] | `kenya2026_rasmus/` | `/course/kenya2026/rasmus` | 75 M |
| [x] | `popgen24_garrett/` | `/course/popgen24/garrett` | 55 M |
| [x] | `novo23_wes/` | `/course/novo23/wes` | 46 M |
| [x] | `novo23_scripts/` | `/course/novo23/scripts` | 5.5 M |

### Cleaned, per-exercise data folders

`data/` is being reorganised so each exercise reads from a folder named after it,
rather than from a folder named after the course it happened to be taught in.
The bulk directories copied above stay for now as the source material.

Done so far:

```
data/NGSintro/
  animal/    wildebeest FASTQ pair + goat reference with bwa index   335 M
  human/     NA19238 chr21 FASTQ pair + chr21 reference with index   115 M
  software/  picard.jar                                              15 M

data/imputation/
  bams/       31 low-depth CEU bams + indexes (30 study + 1 NIPT)      223 M
  vcfs/       ref panel, truth set, fake SNP chip, QUILT2, example      28 M
  resources/  GRCh38 + .fai, QUILT2 and Beagle 5 genetic maps          3.1 G
  software/   beagle 4.1 + 5.5 jars, QUILT distribution                 64 M

data/f_stats/                                                          108 M
  ho_anc.sample_info.tsv   AADR sample metadata, 98 populations
  f2.ho_anc/               precomputed pairwise f2, 50,493 SNPs

data/geneflow/                                                         3.3 M
  hum_nea_siml.vcf.gz      simulated CHIMP/AFR/EUR/NEA data, 81,486 sites
  hum_nea_siml.tsv         sample-to-population table
  wildebeest_fstats_wildebeestref/   precomputed f2 for the wildebeest f4

data/sfs/                                                              291 M
  simld_{a,b,c}_for_sfs.vcf.gz   3 simulated scenarios, unfolded sfs
  blackwildebeest_chr1.vcf.gz             folded sfs, real data
  bluewildebeest_whitebeard_chr1.vcf.gz   folded sfs, real data

data/psmc/                                        shared by #39 and #62   284 M
  data/       simulated + 1000 Genomes psmcfa/psmc, wildebeest bcf for the bonus
  animal/     wildebeest psmcfa and precomputed psmc results
  images/     popsize, bootstrap, 1kg_chr1, NA12718 figures
  scripts/    vcf2psmcfa.py
  software/   the psmc binary and utils/psmc_plot.pl

data/sfs_fst_pbs/                                                       2.4 G
  hg19.fa.gz + hg19ancNoChr.fa.gz    human reference and chimp ancestral (+ indexes)
  smallerbams/     30 bams + indexes, 10 each CEU/JPT/YRI, reduced genome
  chr5_33M_v2/     30 bams, the 1 Mb region on chromosome 5
  precomputed/     saf files, so the angsd step can be skipped

data/relatedness/                                                      1.5 G
  Reindeer.{bed,bim,fam}   111 reindeer, 2.8 M SNPs, 12 Greenland populations
  pcaone.residuals/.mbim   precomputed PCAone residuals for the LD pruning step
  Reindeer_pruned.2.{P,Q}  precomputed ADMIXTURE K=2, in case the run is too slow
  population.tsv           population labels for the PCA plot
  relateAdmix              the relateAdmix binary

data/gwas_human/                              shared by #46-#48 and #64      58 M
  gwa.{bed,bim,fam}   356 individuals, 499,264 SNPs, 200 controls / 156 cases
  pheno3.txt          an extra phenotype used by the GWAS exercises
```

`chr21.fa.gz` was missing its `.fai`/`.gzi` index in the original course folder,
which `samtools tview` and `bcftools mpileup -f` both need. It was generated
when the clean folder was built.

`fastqc` is not shipped — it is on PATH at `/usr/bin/fastqc`.

### Data still to trace

These exercises never name their data — they source an `env.sh` or assume a
pre-staged home directory. Trace each before its data can be placed in `data/`
and its paths rewritten.

| Exercise | Indirection |
|---|---|
| #15, #18, #19, #21, #23, #24, #26 (all advBinf) | `~/advBinf/admixtureHuman/env.sh`, `~/advBinfImputation/env.sh` |
| #48, #49, #54, #55, #56 (chinaCourse2025) | `~/sysu_*`, `~/.sysu_admixture_paths` |
| #14, #22, #30, #41 (summer2025) | `~/current_folder`, `~/data_folder`, `~/popgen25_*` |
| #50 (`3WESfamdiab`) | `/home/student/USER/GWAS` |

# Build order

1. Data consolidation — **in progress**.
2. Trace the `env.sh`-staged exercises and copy their data in too.
3. Copy each remaining exercise into its theme folder under its new name.
4. Rewrite every exercise's input paths to resolve under `data/`.
5. Keep `MANIFEST.md` current.
6. Smoke-test notebooks per `agentGuide/RUN_DATA_NOTEBOOKS.md`.

# Scope

**Only exercises already in the GitHub repo are in scope** (decided 2026-09-16).
Exercises that exist on the server but are not tracked in git — for example
`/course/bsa/.../pairwise_alignment.ipynb` and `blast.ipynb` — are deliberately
**not** included. Do not add them.

**Removed from scope (2026-09-16):** #17 `phasing_shapeit_human.ipynb`, from
`bgi23/ 01.phasing.SHAPEIT_v2.ipynb` (2023-11-09). Removed by request. It was the
only SHAPEIT/phasing exercise in the list and the only copy of it, so nothing
supersedes it and nothing else in the list depends on it. Number 17 is retired
rather than reused — the numbering already skips, so later exercises keep their
numbers. `bgi23/ 01.phasing.SHAPEIT_v2.ipynb` stays untouched in `bgi23/` as
history (R6).

Its data was never copied: the notebook reads `/course/bgi23/malthe/thursday/`,
which is not in the data-copy list above, so removing it orphans nothing.
(`bgi23_quan/`, 27 G, is for the single-cell exercise, not this one.)

# Open questions

- `data/` will need cleaning so its contents line up with the exercises. Deferred
  for now; the priority is settling which exercises are in.

# Removed from scope

**#57 `scrna_seurat_human`** (from `bgi23/scRNA_Seurat_Yano.ipynb`, 2023-11-08) —
**removed by request, 2026-09-16.** Single-cell RNA-seq is not population genetics and
does not belong in this set. The `transcriptomics/` theme is gone with it; the source
notebook stays untouched in `bgi23/`.
