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
  linux/  statistics/  sequence_analysis/  shiny/
  ngs/  genotype_calling_imputation/
  pca/  admixture/  local_ancestry/  gene_flow/
  demography/  selection/  relatedness_diversity/
  gwas/  transcriptomics/
data -> /course/data/
```

**56 exercises** (numbering runs to 57; #6 and #7 are one exercise). Done so far: 12.

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

# sequence_analysis/

### 8. [ ] `motif_discovery.R` — **DEFERRED by request (2026-09-16)**
- **From:** `BSA/motif_discovery_ex.R` (2025-10-22)
- **Supersedes:** none — only copy
- **Data: COPIED** to `data/BSA/` — `motif_discovery/PUM2.top500.fa` and
  `pHMM/globins4.fasta`. The `~/work/COURSES/BIO/BSA` paths in the script are the
  workgroup share mounted into home directories.
- **Remaining blocker:** the source is an instructor working/solution draft, not
  a student-facing exercise — answers are inline as comments and lines 100-102
  (`enrichment0: # tgta / lambda0`) are not valid R. Decide whether to ship it as
  a solution script or rewrite it as an exercise.

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

### 15. [ ] `genotype_calling_and_imputation_human.ipynb` — **merged exercise**
- **From:** `advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb` (2026-09-09)
- **Supersedes:**
  - `advBinf/exercises/genotype calling and haplotype Imputation.ipynb` (2025-09-12)
  - `advBinf/exercises/SNPandGenotypeCalling.md` (2024-09-13)
- **Paired with:** #16, the imputation half kept separately

### 16. [ ] `imputation_human.ipynb` — **separate half of #15**
- **From:** `chinacourse2026/Day2_Afternoon_Genotype_Imputation.ipynb` (2026-09-09)
- **Supersedes:**
  - `summer2025/exercises/Day2_Imputation.ipynb` (2025-08-05)
  - `chinaCourse2025/Day2_Afternoon_QUILT_Imputation.ipynb` (2025-07-28)
  - `bgi23/03.QUILT_Imputation_new.ipynb` (2023-11-09)
  - `bgi23/03.QUILT_Imputation.ipynb` (2023-11-09)
  - `bgi23/02.Minimac4_Imputaion.ipynb` (2023-11-09)

### 17. [ ] `phasing_shapeit_human.ipynb`
- **From:** `bgi23/ 01.phasing.SHAPEIT_v2.ipynb` (2023-11-09) — note leading space
- **Supersedes:** none — only copy

---

# pca/

### 18. [ ] `pca_human.ipynb`
- **From:** `advBinf/exercises/advBinf_PCA.ipynb` (2026-09-16)
- **Supersedes:**
  - `chinacourse2026/Day4_Afternoon_PCA_1.ipynb` (2026-09-14)
  - `summer2025/exercises/Day5_PCA_1.ipynb` (2025-08-06)
  - `chinaCourse2025/Day4_Afternoon_PCA_main.ipynb` (2025-07-28)
  - `advBinf/exercises/PCA.md` (2024-09-17)
  - `summer2024/exercises/summer2024-PCA.ipynb` (2024-08-20)

### 19. [ ] `pca_em_human.ipynb`
- **From:** `advBinf/exercises/advBinf_PCA_EM.ipynb` (2026-09-16)
- **Supersedes:** none — new exercise (EMU / PCAngsd EM algorithms)

### 20. [ ] `pca_animal.ipynb`
- **From:** `kenya2026/exercises/Day3/Kenya2026_PCA.ipynb` (2026-08-17)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day4_morning_pca.ipynb` (2026-08-25)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_PCA-V2.ipynb` (2024-08-09)

### 21. [ ] `pca_bonus_animal.ipynb`
- **From:** `advBinf/exercises/advBinf_PCA_bonus.ipynb` (2026-09-16)
- **Supersedes:** `chinaCourse2025/Day4_Afternoon_PCA_bonus.ipynb` (2025-07-28)

### 22. [ ] `pca_called_genotypes_animal.ipynb`
- **From:** `summer2025/exercises/Day5_PCA_2.Call_genotype.ipynb` (2025-08-06)
- **Supersedes:** `summer2024/exercises/summer2024-PCA-CalledGenotypes.ipynb` (2024-08-20)

---

# admixture/

### 23. [ ] `admixture_human.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture.ipynb` (2026-09-14)
- **Supersedes:**
  - `summer2025/exercises/Day3_Morning_Admixture.ipynb` (2025-08-05)
  - `chinaCourse2025/Day4_Morning_admixture_genotype.ipynb` (2025-07-28)
  - `advBinf/exercises/admixture.md` (2024-09-16)
  - `summer2024/exercises/admixExercise_popgen24.ipynb` (2024-08-19)
  - `bgi23/Admixture.ipynb` (2023-11-02)
  - `summer2023/InfererPopStructure/admixExercise_popgen23.ipynb` (2023-08-08)

### 24. [ ] `admixture_em_human.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture_EM.ipynb` (2026-09-14)
- **Supersedes:** none — new exercise (EM algorithm behind ADMIXTURE)

### 25. [ ] `admixture_animal.ipynb`
- **From:** `kenya2026/exercises/Day3/Exercises_Admixture_Kenya26_WoA.ipynb` (2026-08-21)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day3_morning_admixture.ipynb` (2026-08-25)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_AdmixtureV2.ipynb` (2024-08-08)
  - `kenya2024/exercises/day3_PopulationStructure/Admixture.ipynb` (2024-08-07)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_Admixture.ipynb` (2024-07-28)

### 26. [ ] `admixture_bonus_animal.ipynb`
- **From:** `advBinf/exercises/advBinf_admixture_bonus.ipynb` (2026-09-14)
- **Supersedes:**
  - `summer2025/exercises/Day3_Admixture_structure_bonus.ipynb` (2025-08-05)
  - `chinaCourse2025/Day4_Morning_Admixture.bonus.ipynb` (2025-07-28)

### 27. [ ] `population_structure_ii_human.ipynb`
- **From:** `bgi23/BGI2023-populationStructureII.ipynb` (2023-11-02)
- **Supersedes:**
  - `summer2023/InfererPopStructure/popstructII2023.ipynb` (2023-08-09)
  - `summer2023/popstructureII/popstructII2023.ipynb` (2023-08-08)
  - `summer2023/popstructureII/index.md` (2023-08-08)

---

# local_ancestry/

### 28. [ ] `local_ancestry_human.ipynb`
- **From:** `chinacourse2026/Day4_admix_eval_LAI.ipynb` (2026-07-30)
- **Supersedes:** `summer2025/exercises/Day4_Morning_LocalAncestry.ipynb` (2025-08-06)

### 29. [ ] `local_ancestry_hapla_human.ipynb`
- **From:** `advBinf/exercises/Hapla_LAI_exercise.ipynb` (2025-10-07)
- **Supersedes:** none — only copy (hapla-based, distinct from #28)

---

# gene_flow/

### 30. [ ] `f_stats_human.ipynb`
- **From:** `summer2025/exercises/Day3_f_stats.ipynb` (2025-08-05)
- **Supersedes:**
  - `summer2024/exercises/f_stats.ipynb` (2024-08-21)
  - `summer2023/DfFstats/popgen23_f_stats.ipynb` (2023-08-09)

### 31. [ ] `gene_flow_dstat_animal.ipynb`
- **From:** `kenya2026/exercises/Day3/Geneflow&Dstat.ipynb` (2026-08-16)
- **Supersedes:** `kenya2026/exercises/post_course/day3_afternoon_gene_flow_dstat.ipynb` (2026-08-25)

### 32. [ ] `admixture_graphs_human.ipynb`
- **From:** `summer2023/DfFstats/popgen23.Admixture_Graphs_Tutorial.ipynb` (2023-08-10)
- **Supersedes:** none — only copy

### 33. [ ] `chromopainter_finestructure_human.ipynb`
- **From:** `summer2024/exercises/ChromoPainterFineSTRUCTUREPractical.ipynb` (2024-08-21)
- **Supersedes:** none — only copy
- **Companion:** `CopenhagenPopgenWorkshop2024_ChromoPainterFineSTRUCTUREPracticalSOLN.pdf`

### 34. [ ] `dating_admixture_human.ipynb`
- **From:** `summer2024/exercises/DatingAdmixture.ipynb` (2024-08-21)
- **Supersedes:** none — only copy
- **Companion:** `CopenhagenPopgenWorkshop2024_DatingAdmixturePracticalSOLN.pdf`

---

# demography/

### 35. [ ] `coalescence.ipynb`
- **From:** `kenya2026/exercises/Day2/Coalescence_short_WoA.ipynb` (2026-08-15)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day2_morning_coalescence.ipynb` (2026-08-25)
  - `summer2025/exercises/Day1_morning_CoalTutorial.ipynb` (2025-08-03)

### 36. [ ] `wright_fisher.ipynb`
- **From:** `summer2025/exercises/Day1_morning_WrightFisherTutorial.ipynb` (2025-08-03)
- **Supersedes:** none — only copy

### 37. [ ] `sfs_model.ipynb`
- **From:** `advBinf/exercises/advBinf_SFSmodel.ipynb` (2026-09-16)
- **Supersedes:** `advBinf/exercises/SFS.md` (2024-09-20)

### 38. [ ] `sfs_animal.ipynb`
- **From:** `kenya2026/exercises/Day2/SFS_WoA.ipynb` (2026-08-15)
- **Supersedes:** `kenya2026/exercises/post_course/day2_morning_sfs.ipynb` (2026-08-25)

### 39. [ ] `psmc_demography_animal.ipynb`
- **From:** `kenya2026/exercises/Day2/psmc_kenya2026.ipynb` (2026-08-19)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day2_afternoon_psmc.ipynb` (2026-08-25)
  - `summer2025/exercises/Day5_demography.ipynb` (2025-08-07)
  - `summer2024/exercises/summer2024-PSMC_tutorial_2024.ipynb` (2024-08-23)
  - `bgi23/PSMC_tutorial.ipynb` (2023-11-07)
  - `summer2023/DemographyInference/PSMC_tutorial.ipynb` (2023-08-10)

---

# selection/

### 40. [ ] `selection_scans_animal.ipynb`
- **From:** `kenya2026/exercises/Day4/SelectionScans_22nd.ipynb` (2026-08-22)
- **Supersedes:**
  - `kenya2026/exercises/post_course/day4_afternoon_selection_scans.ipynb` (2026-08-25)
  - `summer2023/selectionScan/README.md` (2023)

### 41. [ ] `selection_scans_popgen_animal.ipynb`
- **From:** `summer2025/exercises/Day4_SelectionPopGen2025.ipynb` (2025-08-06)
- **Supersedes:** `summer2024/exercises/SelectionScans.ipynb` (2024-08-22)
- **Note:** a different exercise from #40, not an older version of it — both kept

---

# relatedness_diversity/

### 42. [ ] `fst_animal.ipynb` — **post-course small-dataset version**
- **From:** `kenya2026/exercises/post_course/day4_morning_fst.ipynb` (2026-08-25)
- **Supersedes:**
  - `kenya2026/exercises/Day4/Fst_Kenya2026.ipynb` (2026-08-23)
  - `kenya2024/exercises/day3_PopulationStructure/Day3_Fst_RH.ipynb` (2024-08-09)
- **Paired with:** #44, the merged Related&Fst form

### 43. [ ] `relatedness_animal.ipynb`
- **From:** `kenya2026/exercises/Day5/Related.ipynb` (2026-08-23)
- **Supersedes:** `kenya2026/exercises/post_course/day5_morning_relatedness.ipynb` (2026-08-25)
- **Paired with:** #44, the merged Related&Fst form

### 44. [ ] `relatedness_and_fst_animal.ipynb` — **merged form of #42 + #43**
- **From:** `kenya2026/exercises/Day4/Related&Fst.ipynb` (2026-08-22)
- **Supersedes:** none — kept deliberately alongside the split versions

### 45. [ ] `heterozygosity_roh_animal.ipynb` — **post-course small-dataset version**
- **From:** `kenya2026/exercises/post_course/day5_morning_heterozygosity_roh.ipynb` (2026-08-25)
- **Supersedes:**
  - `kenya2026/exercises/Day5/Day5_GeneticDiversity.ipynb` (2026-08-23)
  - `kenya2024/exercises/day2/Day2_Inbreeding_ROH.ipynb` (2024-08-08)
  - `kenya2024/exercises/day2/Inbreeding_ROH.ipynb` (2024-08-07)
  - `kenya2024/exercises/Day1_GeneticDiversity/Day1_GeneticDiversity.ipynb` (2024-08-06)

---

# gwas/

### 46. [ ] `gwas_intro_human.ipynb`
- **From:** `novCourse2024/1GWASIntro.ipynb` (2025-07-08)
- **Supersedes:** `bgi23/04.GWASintro_2023_SAIGE.ipynb` (2023-11-09)

### 47. [ ] `gwas_sumstats_human.ipynb`
- **From:** `novCourse2024/2GWASsumstats.ipynb` (2025-07-08)
- **Supersedes:** none — only copy

### 48. [ ] `gwas_analysis_human.ipynb`
- **From:** `chinaCourse2025/Day3_GWAS_Analysis_2025_Morning.ipynb` (2025-07-26)
- **Supersedes:** none — only copy (distinct from #46 / #47)

### 49. [ ] `gene_based_testing_human.ipynb`
- **From:** `chinaCourse2025/Day3_Gene_Based_Testing_2025_Afternoon.ipynb` (2025-07-26)
- **Supersedes:** none — only copy

### 50. [ ] `wes_family_diabetes_human.ipynb`
- **From:** `novCourse2024/3WESfamdiab.ipynb` (2025-07-08)
- **Supersedes:** none — only copy

### 51. [ ] `wes_familial_hypercholesterolemia_human.ipynb`
- **From:** `novCourse2024/4WESfh.ipynb` (2025-07-08)
- **Supersedes:** none — only copy

### 52. [ ] `wes_diabetes_human.ipynb`
- **From:** `novCourse2024/5WESdiab.ipynb` (2025-07-08)
- **Supersedes:** none — only copy

### 53. [ ] `heritability_ldscore_human.ipynb`
- **From:** `chinacourse2026/Day5_Morning_heribilty_and_ldscore.ipynb` (2026-08-07)
- **Supersedes:**
  - `chinaCourse2025/Day5_heritability_exercise.ipynb` (2025-07-26)
  - `chinaCourse2025/Day5_Afternoon_Genetic_correlation_Partitioned_Heritability.ipynb` (2025-07-26)

### 54. [ ] `prs_height_human.ipynb`
- **From:** `chinaCourse2025/Day6_Morning1_PRS_height_pipeline.ipynb` (2025-07-26)
- **Supersedes:** none — only copy

### 55. [ ] `mendelian_randomization_human.ipynb`
- **From:** `chinaCourse2025/Day6_Morning2_MR.ipynb` (2025-07-26)
- **Supersedes:**
  - `bgi23/MR-exercise.ipynb` (2023-11-10)
  - `bgi23/MR.real_data.exercise.ipynb` (2023-11-10)

### 56. [ ] `mendelian_randomization_proteomics_human.ipynb`
- **From:** `chinaCourse2025/Day6_Afternoon_Proteomics_MR.ipynb` (2025-07-26)
- **Supersedes:** none — only copy

---

# transcriptomics/

### 57. [ ] `scrna_seurat_human.ipynb`
- **From:** `bgi23/scRNA_Seurat_Yano.ipynb` (2023-11-08)
- **Supersedes:** none — only copy
- **Companion:** `bgi23/scRNA_check_files.ipynb` (2023-11-07), `bgi23/glioblastoma.html`

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

# Data consolidation into `data/` (= `/course/data/`)

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

# Open questions

- `data/` will need cleaning so its contents line up with the exercises. Deferred
  for now; the priority is settling which exercises are in.
