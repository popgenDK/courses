# current_exercises

One current version of every exercise, organised by theme, with the human and
animal versions of an analysis side by side. Everything reads its input from the
shared data store `data/` (a symlink to `/course/data/`).

**58 of 58 exercises built.**

- **The build list, with full provenance:** [`EXERCISES.md`](EXERCISES.md)
- **What has been built and what changed in each file:** [`MANIFEST.md`](MANIFEST.md)

## How these exercises work

- **Notebooks** are [SoS](https://vatlab.github.io/sos-docs/) notebooks mixing
  bash, R and Python cells. Open them in Jupyter and pick the `SoS` kernel.
- **All paths are set in the first cell.** If the data or software moves, that
  one cell is the only thing to change — nothing further down uses a full path.
- **Quizzes and questions.** Notebooks carry `jupyterquiz` quizzes and a short
  set of questions after most code cells. Quiz files live in `<theme>/quiz/`.
- **Markdown exercises** (`.md`) are run in a terminal and carry no quizzes.
- **Shiny apps** are launched straight from R with the `source(...)` line at the
  top of each file.

## Solutions (HTML)

A solution is the notebook **executed with all its output**, exported to HTML and
saved next to the notebook as `<name>.html`.

**The links in the Solution column open the rendered page**, not its source.
GitHub serves an `.html` file as plain text rather than as a page, so every
solution link is routed through
[htmlpreview.github.io](https://htmlpreview.github.io/), which fetches the raw
file and renders it. Nothing needs to be downloaded. If a preview fails to load
— the largest render here is 4 MB and the service can time out — fall back to the
`.html` file in the repo and open it locally.

**47 of 51 are rendered**, each linked beside its exercise below. The four
without one — `gwas_analysis_human`, `gene_based_testing_human`,
`prs_height_human` and `proteomics_mr_human` — have no data on this server, so
they are expected to stay blank. `heritability_ldscore_human` renders as far as
the LD score regression, which needs a conda environment that is not installed.

The renders in the original course folders are of the *old* notebooks — before
the quizzes, questions, typo fixes and path changes — so they are deliberately
not copied here.

## Exercises

### Linux and the command line — `linux/`

**1 · [`intro_linux.md`](https://github.com/popgenDK/courses/blob/main/current_exercises/linux/intro_linux.md)**  
Linux and bash basics: navigating, editing, permissions, pipes. No data. Run in a terminal.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`summer2025/BriefIntro2Linux.md`](https://github.com/popgenDK/courses/blob/main/summer2025/BriefIntro2Linux.md) (2025-07-23)

Replaces:

- [`summer2024/BriefIntro2Linux.md`](https://github.com/popgenDK/courses/blob/main/summer2024/BriefIntro2Linux.md) (2024-07-04)

</details>

**2 · [`intro_bash_linux.md`](https://github.com/popgenDK/courses/blob/main/current_exercises/linux/intro_bash_linux.md)**  
Bash scavenger hunt: directories, text processing, archives, processes. No data. Run in a terminal.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day1/IntroToBash.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day1/IntroToBash.ipynb) (2026-08-18)

Replaces:

- [`kenya2026/exercises/post_course/day1_morning_bash_linux.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day1_morning_bash_linux.ipynb) (2026-08-25)

</details>


### Shiny apps — `shiny/`

Interactive R Shiny apps. Each is a single self-contained `.R` file, so the `source(...)` line below launches it straight from GitHub — nothing to clone, and no packages beyond the ones the app itself uses. (`shiny::runUrl()` does not apply here: it expects a zipped app directory, not a single file.)

One app also has a **run in your browser** link. That is [shinylive](https://shinylive.io), which compiles the app to WebAssembly and runs it client-side with webR — no R installation and no server, with the whole app carried in the URL. Expect 10–30 seconds on first load while the R runtime downloads. webR does not carry every CRAN package, so not every app can be published this way.

**9 · [`needleman_wunsch_dna.R`](https://github.com/popgenDK/courses/blob/main/current_exercises/shiny/needleman_wunsch_dna.R)**  
Needleman-Wunsch pairwise alignment of **DNA**, as a Shiny app. Sequences, not genotypes.  
Run in R: `source("https://raw.githubusercontent.com/popgenDK/courses/main/current_exercises/shiny/needleman_wunsch_dna.R")`  
*From [`BSA/NW_DNA.R`](https://github.com/popgenDK/courses/blob/main/BSA/NW_DNA.R) (2025-09-02) — the only copy*

**10 · [`needleman_wunsch_blosum50.R`](https://github.com/popgenDK/courses/blob/main/current_exercises/shiny/needleman_wunsch_blosum50.R)**  
Needleman-Wunsch pairwise alignment of **protein** with BLOSUM50, as a Shiny app. Sequences, not genotypes.  
Run in R: `source("https://raw.githubusercontent.com/popgenDK/courses/main/current_exercises/shiny/needleman_wunsch_blosum50.R")`  
*From [`BSA/needleman_wunsch_shiny_app_blosum_50.r`](https://github.com/popgenDK/courses/blob/main/BSA/needleman_wunsch_shiny_app_blosum_50.r) (2025-08-30) — the only copy*

**11 · [`dotplot.R`](https://github.com/popgenDK/courses/blob/main/current_exercises/shiny/dotplot.R)**  
Dot plot of two sequences from FASTA files you upload. Any species.  
Run in R: `source("https://raw.githubusercontent.com/popgenDK/courses/main/current_exercises/shiny/dotplot.R")`  
*From [`BSA/dotplotShiny.R`](https://github.com/popgenDK/courses/blob/main/BSA/dotplotShiny.R) (2025-08-30) — the only copy*

**3 · [`stats_binomial.R`](https://github.com/popgenDK/courses/blob/main/current_exercises/shiny/stats_binomial.R)**  
The binomial distribution, as an interactive Shiny app. No data.  
Run in R: `source("https://raw.githubusercontent.com/popgenDK/courses/main/current_exercises/shiny/stats_binomial.R")`  
*From [`stat_molbio/binom.R`](https://github.com/popgenDK/courses/blob/main/stat_molbio/binom.R) (2026-01-16) — the only copy*

**4 · [`stats_normal.R`](https://github.com/popgenDK/courses/blob/main/current_exercises/shiny/stats_normal.R)**  
The normal distribution, as an interactive Shiny app. No data.  
Run in R: `source("https://raw.githubusercontent.com/popgenDK/courses/main/current_exercises/shiny/stats_normal.R")` · [run in your browser][shinylive-stats_normal]  
*From [`stat_molbio/normal.R`](https://github.com/popgenDK/courses/blob/main/stat_molbio/normal.R) (2026-01-16) — the only copy*


### EM algorithms — `em_algorithms/`

Exercises that build an EM algorithm from scratch.

**5 · [`em_algorithm.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/em_algorithm.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/em_algorithm.html)  
The EM algorithm from first principles: two coins, allele frequencies from genotype likelihoods, and the binomial. Simulated data.  
*From [`advBinf/exercises/advBinf_EM_algorithm.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_EM_algorithm.ipynb) (2026-09-09) — the only copy*

**6+7 · [`haplotype_frequencies.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/haplotype_frequencies.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/haplotype_frequencies.html)  
An EM algorithm for haplotype frequencies from 2-SNP genotypes. Simulated: 1000 individuals.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`advBinf/exercises/solution_haplotype_frequencies.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/solution_haplotype_frequencies.ipynb) (2025-09-12)

Replaces:

- [`advBinf/exercises/haplotype_frequencies.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/haplotype_frequencies.ipynb) (2025-09-12) — the exercise half, which was too hard to work through

</details>

**19 · [`pca_em_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/pca_em_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/pca_em_human.html)  
The EM algorithms behind EMU and PCAngsd, and randomized SVD. Simulated: 60 individuals, 3 populations, 10,000 SNPs.  
*From [`advBinf/exercises/advBinf_PCA_EM.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_PCA_EM.ipynb) (2026-09-16) — the only copy*

**24 · [`admixture_em_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/admixture_em_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/admixture_em_human.html)  
The EM algorithm behind ADMIXTURE and NGSadmix. Simulated: 50 individuals, 2 ancestral populations.  
*From [`advBinf/exercises/advBinf_admixture_EM.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_admixture_EM.ipynb) (2026-09-14) — the only copy*

**37 · [`sfs_model.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/sfs_model.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/em_algorithms/sfs_model.html)  
Estimating the site frequency spectrum from genotype likelihoods by EM. Simulated: 10 individuals, 100,000 sites.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`advBinf/exercises/advBinf_SFSmodel.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_SFSmodel.ipynb) (2026-09-16)

Replaces:

- [`advBinf/exercises/SFS.md`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/SFS.md) (2024-09-20)

</details>


### NGS data and mapping — `ngs/`

**12 · [`ngs_intro_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/ngs/ngs_intro_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/ngs/ngs_intro_human.html)  
FASTQ &rarr; QC &rarr; bwa mapping &rarr; SAM/BAM &rarr; VCF. **Human**: NA19238 (YRI, Nigeria), chr21. Reads, not genotypes.  
<details><summary>Where it comes from — replaces 5 older copies</summary>

From [`chinacourse2026/Day2_Morning_NGSintro_human.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day2_Morning_NGSintro_human.ipynb) (2026-09-14)

Replaces:

- [`chinaCourse2025/Day2_Morning_NGSintro_human.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day2_Morning_NGSintro_human.ipynb) (2025-08-04)
- [`summer2025/exercises/Day1_afternoon_NGSintro_human.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day1_afternoon_NGSintro_human.ipynb) (2025-08-04)
- [`summer2024/exercises/NGSintro.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/NGSintro.ipynb) (2024-08-19)
- [`bgi23/NGSintro.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/NGSintro.ipynb) (2023-10-30)
- [`summer2023/IntroNGS/introNGSexercises.md`](https://github.com/popgenDK/courses/blob/main/summer2023/IntroNGS/introNGSexercises.md) (2023-08-07)

</details>

**13 · [`ngs_intro_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/ngs/ngs_intro_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/ngs/ngs_intro_animal.html)  
FASTQ &rarr; QC &rarr; bwa mapping &rarr; SAM/BAM &rarr; VCF. **Blue wildebeest**, mapped to a **goat** reference. Reads, not genotypes.  
<details><summary>Where it comes from — replaces 6 older copies</summary>

From [`kenya2026/exercises/Day1/Kenya2026_NGSintro.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day1/Kenya2026_NGSintro.ipynb) (2026-08-17)

Replaces:

- [`kenya2026/exercises/post_course/day1_afternoon_ngs_intro.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day1_afternoon_ngs_intro.ipynb) (2026-08-25)
- [`summer2025/exercises/Day1_afternoon_NGSintro_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day1_afternoon_NGSintro_animal.ipynb) (2025-08-04)
- [`chinaCourse2025/Day2_Morning_NGSintro_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day2_Morning_NGSintro_animal.ipynb) (2025-08-04)
- [`kenya2024/exercises/day1_NGSintro/Day1_NGSintroV4.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day1_NGSintro/Day1_NGSintroV4.ipynb) (2024-08-07)
- [`kenya2024/exercises/day1_NGSintro/Day1_NGSintroV3.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day1_NGSintro/Day1_NGSintroV3.ipynb) (2024-07-26)
- [`kenya2024/exercises/day1_NGSintro/Day1_NGSintroV2.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day1_NGSintro/Day1_NGSintroV2.ipynb) (2024-07-15)

</details>

**14 · [`ngs_inference_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/ngs/ngs_inference_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/ngs/ngs_inference_human.html)  
Low-depth inference with ANGSD: genotype likelihoods, calling, allele frequencies, SNP calling at EDAR. **Human**: 100 individuals from LWK, TSI, CHB, PEL, NAM. **Low depth** (~2-3x): genotype likelihoods, never called.  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`summer2025/exercises/Day2_NGS_Inference.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day2_NGS_Inference.ipynb) (2025-08-04)

Replaces:

- [`summer2024/exercises/NGS_inference.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/NGS_inference.ipynb) (2024-08-19)
- [`summer2023/NGSinference/README.md`](https://github.com/popgenDK/courses/blob/main/summer2023/NGSinference/README.md) + `solutions.md` (2023-08-05)

</details>


### Genotype calling and imputation — `genotype_calling_imputation/`

**15 · [`genotype_calling_and_imputation_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/genotype_calling_imputation/genotype_calling_and_imputation_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/genotype_calling_imputation/genotype_calling_and_imputation_human.html)  
SNP calling, genotype calling and imputation compared. **Human**: 30 CEU samples, 3 Mb of chr20. **Low depth**: the point is comparing calling with imputation.  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb) (2026-09-09)

Replaces:

- [`advBinf/exercises/genotype calling and haplotype Imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/genotype%20calling%20and%20haplotype%20Imputation.ipynb) (2025-09-12)
- [`advBinf/exercises/SNPandGenotypeCalling.md`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/SNPandGenotypeCalling.md) (2024-09-13)

</details>

**16 · [`imputation_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/genotype_calling_imputation/imputation_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/genotype_calling_imputation/imputation_human.html)  
Genotype imputation with a reference panel. **Human**: 30 CEU samples, 3 Mb of chr20. **Low depth**: imputation instead of calling.  
<details><summary>Where it comes from — replaces 5 older copies</summary>

From [`advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_genotype_calling_and_imputation.ipynb) (2026-09-09) for the imputation sections, plus [`chinacourse2026/Day2_Afternoon_Genotype_Imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day2_Afternoon_Genotype_Imputation.ipynb) (2026-09-09) for the three sections only it has

Replaces:

- [`summer2025/exercises/Day2_Imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day2_Imputation.ipynb) (2025-08-05)
- [`chinaCourse2025/Day2_Afternoon_QUILT_Imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day2_Afternoon_QUILT_Imputation.ipynb) (2025-07-28)
- [`bgi23/03.QUILT_Imputation_new.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/03.QUILT_Imputation_new.ipynb) (2023-11-09)
- [`bgi23/03.QUILT_Imputation.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/03.QUILT_Imputation.ipynb) (2023-11-09)
- [`bgi23/02.Minimac4_Imputaion.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/02.Minimac4_Imputaion.ipynb) (2023-11-09)

</details>


### PCA — `pca/`

**58 · [`pca_mds_and_svd.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_mds_and_svd.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_mds_and_svd.html)  
MDS and PCA worked by hand: distances, cmdscale, the SVD, the covariance matrix, variance explained. No data &mdash; a 5x7 matrix typed in.  
*From extracted from [`advBinf/exercises/advBinf_PCA.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_PCA.ipynb) cells 2-25 (2026-09-16) — the only copy*

**18 · [`pca_low_depth_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_low_depth_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_low_depth_human.html)  
PCAngsd on genotype likelihoods. **Human**: 435 individuals from ASW, CEU, CHB, MXL, YRI. **Low depth**: genotype likelihoods.  
<details><summary>Where it comes from — replaces 4 older copies</summary>

From [`advBinf/exercises/advBinf_PCA.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_PCA.ipynb) (2026-09-16), from cell 26 on

Replaces:

- [`summer2025/exercises/Day5_PCA_1.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day5_PCA_1.ipynb) (2025-08-06)
- [`chinaCourse2025/Day4_Afternoon_PCA_main.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day4_Afternoon_PCA_main.ipynb) (2025-07-28)
- [`advBinf/exercises/PCA.md`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/PCA.md) (2024-09-17)
- [`summer2024/exercises/summer2024-PCA.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/summer2024-PCA.ipynb) (2024-08-20)

</details>

**59 · [`pca_low_depth_selection_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_low_depth_selection_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_low_depth_selection_human.html)  
PC-based selection scan with PCAngsd. **Human**: 424 Europeans from CEU, GBR, IBS, TSI. **Low depth**: genotype likelihoods.  
*From [`chinacourse2026/Day4_Afternoon_PCA_1.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day4_Afternoon_PCA_1.ipynb) (2026-09-14), the **PC-based selection** half (from the `# PC-based selection` heading) — the only copy*

**22 · [`pca_called_genotypes_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_called_genotypes_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_called_genotypes_human.html)  
PCAone on LD-pruned called genotypes, read against admixture proportions. **Human**: 192 individuals, 16 populations. **Called genotypes** (LD-pruned plink).  
*From [`chinacourse2026/Day4_Afternoon_PCA_1.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day4_Afternoon_PCA_1.ipynb) (2026-09-14), the **first** half, up to the `# PC-based selection` heading — the only copy*

**20 · [`pca_called_genotypes_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_called_genotypes_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_called_genotypes_animal.html)  
PCAone on called genotypes, LD pruning adjusted for structure, and an IBS tree. **Blue wildebeest**: 73 individuals, 7 localities. **Called genotypes** (plink).  
<details><summary>Where it comes from — replaces 6 older copies</summary>

From [`advBinf/exercises/advBinf_PCA_bonus.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_PCA_bonus.ipynb) (2026-09-16)

Replaces:

- [`kenya2026/exercises/Day3/Kenya2026_PCA.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day3/Kenya2026_PCA.ipynb) (2026-08-17)
- [`kenya2026/exercises/post_course/day4_morning_pca.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day4_morning_pca.ipynb) (2026-08-25)
- [`summer2025/exercises/Day5_PCA_2.Call_genotype.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day5_PCA_2.Call_genotype.ipynb) (2025-08-06)
- [`chinaCourse2025/Day4_Afternoon_PCA_bonus.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day4_Afternoon_PCA_bonus.ipynb) (2025-07-28)
- [`summer2024/exercises/summer2024-PCA-CalledGenotypes.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/summer2024-PCA-CalledGenotypes.ipynb) (2024-08-20)
- [`kenya2024/exercises/day3_PopulationStructure/Day3_PCA-V2.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day3_PopulationStructure/Day3_PCA-V2.ipynb) (2024-08-09)

</details>


### Admixture, local ancestry and gene flow — `admixture/`

**23 · [`admixture_low_depth_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_low_depth_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_low_depth_human.html)  
NGSadmix on genotype likelihoods, evalAdmix and the choice of K. **Human**: 435 individuals from ASW, CEU, CHB, MXL, YRI. **Low depth**: genotype likelihoods.  
<details><summary>Where it comes from — replaces 5 older copies</summary>

From [`advBinf/exercises/advBinf_admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_admixture.ipynb) (2026-09-14), cells 5-70

Replaces:

- [`summer2025/exercises/Day3_Morning_Admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_Morning_Admixture.ipynb) (2025-08-05)
- [`advBinf/exercises/admixture.md`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/admixture.md) (2024-09-16)
- [`summer2024/exercises/admixExercise_popgen24.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/admixExercise_popgen24.ipynb) (2024-08-19)
- [`bgi23/Admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/Admixture.ipynb) (2023-11-02)
- [`summer2023/InfererPopStructure/admixExercise_popgen23.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2023/InfererPopStructure/admixExercise_popgen23.ipynb) (2023-08-08)

</details>

**60 · [`admixture_reference_panel_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_reference_panel_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_reference_panel_human.html)  
fastNGSadmix: the ancestry of a single individual against a fixed panel. **Human**: 7 reference populations, 195 individuals. **Low depth**: genotype likelihoods for one individual.  
*From [`advBinf/exercises/advBinf_admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_admixture.ipynb) (2026-09-14), cells 71-93 — the only copy*

**61 · [`admixture_called_genotypes_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_called_genotypes_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_called_genotypes_human.html)  
ADMIXTURE on called genotypes, convergence across seeds, and evalAdmix. **Human**: 192 individuals, 16 populations. **Called genotypes** (LD-pruned plink).  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`chinacourse2026/Day4_admix_eval_LAI.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day4_admix_eval_LAI.ipynb) (2026-07-30), cells 2-54

Replaces:

- [`chinaCourse2025/Day4_Morning_admixture_genotype.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day4_Morning_admixture_genotype.ipynb) — **no**, see #25

</details>

**25 · [`admixture_called_genotypes_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_called_genotypes_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_called_genotypes_animal.html)  
ADMIXTURE on called genotypes, seeds, K and model fit. **Blue wildebeest**: 73 individuals, 7 localities. **Called genotypes** (plink).  
<details><summary>Where it comes from — replaces 6 older copies</summary>

From [`advBinf/exercises/advBinf_admixture_bonus.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/advBinf_admixture_bonus.ipynb) (2026-09-14)

Replaces:

- [`kenya2026/exercises/Day3/Exercises_Admixture_Kenya26_WoA.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day3/Exercises_Admixture_Kenya26_WoA.ipynb) (2026-08-21)
- [`kenya2026/exercises/post_course/day3_morning_admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day3_morning_admixture.ipynb) (2026-08-25)
- [`chinaCourse2025/Day4_Morning_admixture_genotype.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day4_Morning_admixture_genotype.ipynb) (2025-07-28)
- [`kenya2024/exercises/day3_PopulationStructure/Day3_AdmixtureV2.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day3_PopulationStructure/Day3_AdmixtureV2.ipynb) (2024-08-08)
- [`kenya2024/exercises/day3_PopulationStructure/Admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day3_PopulationStructure/Admixture.ipynb) (2024-08-07)
- [`kenya2024/exercises/day3_PopulationStructure/Day3_Admixture.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day3_PopulationStructure/Day3_Admixture.ipynb) (2024-07-28)

</details>

**28 · [`local_ancestry_flare_mosaic_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/local_ancestry_flare_mosaic_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/local_ancestry_flare_mosaic_human.html)  
Local ancestry with FLARE and MOSAIC, 20 vs 200 generations since admixture. Simulated **human-like** admixed genomes. **Called genotypes** (phased).  
*From [`summer2025/exercises/Day4_Morning_LocalAncestry.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day4_Morning_LocalAncestry.ipynb) (2025-08-06) — the only copy*

**29 · [`local_ancestry_hapla_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/local_ancestry_hapla_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/local_ancestry_hapla_human.html)  
Local ancestry with hapla cluster/admix/fatash. Simulated **human**, then real **cattle** (chr25, 314 individuals). **Called genotypes** (phased).  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`advBinf/exercises/Hapla_LAI_exercise.ipynb`](https://github.com/popgenDK/courses/blob/main/advBinf/exercises/Hapla_LAI_exercise.ipynb) (2025-10-07)

Replaces:

- [`chinacourse2026/Day4_admix_eval_LAI.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day4_admix_eval_LAI.ipynb) cells 55-68 (2026-07-30) — a 14-cell "short look" at hapla/fatash using precomputed files, which this covers fully

</details>

**30 · [`f_stats_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/f_stats_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/f_stats_human.html)  
f2, f3 and f4 statistics and qpAdm. **Human**: 1,646 individuals in 98 populations, ancient and modern, from the AADR. **Called genotypes**, precomputed F2.  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`summer2025/exercises/Day3_f_stats.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day3_f_stats.ipynb) (2025-08-05)

Replaces:

- [`summer2024/exercises/f_stats.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/f_stats.ipynb) (2024-08-21)
- [`summer2023/DfFstats/popgen23_f_stats.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2023/DfFstats/popgen23_f_stats.ipynb) (2023-08-09)

</details>

**31 · [`gene_flow_dstat_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/gene_flow_dstat_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/gene_flow_dstat_animal.html)  
D-statistics (ABBA-BABA) built from scratch on simulated **human/Neanderthal** data, then f4 on **blue and black wildebeest**. **Called genotypes**.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day3/Geneflow&Dstat.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day3/Geneflow%26Dstat.ipynb) (2026-08-16)

Replaces:

- [`kenya2026/exercises/post_course/day3_afternoon_gene_flow_dstat.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day3_afternoon_gene_flow_dstat.ipynb) (2026-08-25)

</details>

**32 · [`admixture_graphs_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_graphs_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/admixture_graphs_human.html)  
Fitting admixture graphs with qpgraph, and estimating them with treemix. **Human**: 33 world populations. **Called genotypes**, precomputed F2.  
*From [`summer2023/DfFstats/popgen23.Admixture_Graphs_Tutorial.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2023/DfFstats/popgen23.Admixture_Graphs_Tutorial.ipynb) (2023-08-10) — the only copy*

**33 · [`chromopainter_finestructure_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/chromopainter_finestructure_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/chromopainter_finestructure_human.html)  
Chromosome painting with ChromoPainter, clustering with fineSTRUCTURE, then GLOBETROTTER and SOURCEFIND. **Human**: 16 populations, 256 individuals. **Called genotypes** (phased haplotypes + recombination map).  
*From [`summer2024/exercises/ChromoPainterFineSTRUCTUREPractical.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/ChromoPainterFineSTRUCTUREPractical.ipynb) (2024-08-21) — the only copy*

**34 · [`dating_admixture_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/dating_admixture_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/admixture/dating_admixture_human.html)  
Dating admixture with ALDER, MALDER, fastGLOBETROTTER and MOSAIC, then AdaptMix for selection. **Human**: 16 populations plus a simulated admixed group. **Called genotypes** (phased).  
*From [`summer2024/exercises/DatingAdmixture.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/DatingAdmixture.ipynb) (2024-08-21) — the only copy*


### Demography and the coalescent — `demography/`

**35 · [`coalescence.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/demography/coalescence.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/demography/coalescence.html)  
The coalescent: simulating and interpreting gene trees. Simulated data.  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`kenya2026/exercises/Day2/Coalescence_short_WoA.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day2/Coalescence_short_WoA.ipynb) (2026-08-15)

Replaces:

- [`kenya2026/exercises/post_course/day2_morning_coalescence.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day2_morning_coalescence.ipynb) (2026-08-25)
- [`summer2025/exercises/Day1_morning_CoalTutorial.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day1_morning_CoalTutorial.ipynb) (2025-08-03)

</details>

**36 · [`wright_fisher.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/demography/wright_fisher.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/demography/wright_fisher.html)  
Wright-Fisher simulations of drift. Simulated data.  
*From [`summer2025/exercises/Day1_morning_WrightFisherTutorial.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day1_morning_WrightFisherTutorial.ipynb) (2025-08-03) — the only copy*

**38 · [`sfs_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/demography/sfs_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/demography/sfs_animal.html)  
The site frequency spectrum from real sequencing data. **Blue wildebeest**. **Low depth**: SFS from genotype likelihoods.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day2/SFS_WoA.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day2/SFS_WoA.ipynb) (2026-08-15)

Replaces:

- [`kenya2026/exercises/post_course/day2_morning_sfs.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day2_morning_sfs.ipynb) (2026-08-25)

</details>

**39 · [`psmc_demography_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/demography/psmc_demography_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/demography/psmc_demography_animal.html)  
PSMC demographic history. **Blue wildebeest**. **Called genotypes**: PSMC needs high coverage.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day2/psmc_kenya2026.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day2/psmc_kenya2026.ipynb) (2026-08-19)

Replaces:

- [`kenya2026/exercises/post_course/day2_afternoon_psmc.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day2_afternoon_psmc.ipynb) (2026-08-25)

</details>

**62 · [`psmc_demography_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/demography/psmc_demography_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/demography/psmc_demography_human.html)  
PSMC demographic history. **Human**: NA12718 (CEU) and NA19471 (Luhya, Kenya). **Called genotypes**: PSMC needs high coverage.  
<details><summary>Where it comes from — replaces 3 older copies</summary>

From [`summer2025/exercises/Day5_demography.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day5_demography.ipynb) (2025-08-07)

Replaces:

- [`summer2024/exercises/summer2024-PSMC_tutorial_2024.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/summer2024-PSMC_tutorial_2024.ipynb) (2024-08-23)
- [`bgi23/PSMC_tutorial.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/PSMC_tutorial.ipynb) (2023-11-07)
- [`summer2023/DemographyInference/PSMC_tutorial.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2023/DemographyInference/PSMC_tutorial.ipynb) (2023-08-10)

</details>


### Selection — `selection/`

See also [`pca/pca_low_depth_selection_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/pca/pca_low_depth_selection_human.ipynb) — a selection scan built on the principal components rather than on Fst, kept with the PCA exercises because it continues directly from them.

**63 · [`sfs_fst_pbs_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/selection/sfs_fst_pbs_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/selection/sfs_fst_pbs_human.html)  
SFS, Fst and PBS on human data with ANGSD. **Human**: CEU and YRI. **Low depth**: genotype likelihoods (ANGSD).  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From written for this set

Replaces:

- [`summer2024/exercises/SelectionScans.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/SelectionScans.ipynb) cells 0-18 (2024-08-22) — the frequency-statistics half only. Its haplotype half is #66 and its PBS scan is #65
- [`kenya2026/exercises/Day4/SelectionScans_22nd.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day4/SelectionScans_22nd.ipynb) cells 0-24 (2026-08-22) — the frequency-statistics half of that notebook, 49% similar to summer2024. Cells 25-37, the genome-wide PBS scan, are **not** covered by this exercise — see #65

</details>

**66 · [`selection_haplotype_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_haplotype_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_haplotype_human.html)  
Extended haplotype homozygosity: selscan iHS and XP-EHH on the lactase region. **Human**: CEU 41, YRI 48, CHB 48. **Called genotypes, phased**, with a genetic map.  
*From [`summer2024/exercises/SelectionScans.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/SelectionScans.ipynb) cells 19-35 (2024-08-22), "Exercise II: Haplotype-based methods" — the only copy*

**65 · [`selection_pbs_scan_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_pbs_scan_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_pbs_scan_human.html)  
Genome-wide PBS scan: Manhattan plot, zoom into peaks, find the gene, check the lactase region. **Human**: NAT, CHB, CEU, YRI. **Called genotypes**, PBS precomputed in 50 kb windows.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day4/SelectionScans_22nd.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day4/SelectionScans_22nd.ipynb) cells 25-37 (2026-08-22), "Exercise II: Whole-genome PBS with 1000 Genomes"

Replaces:

- [`summer2024/exercises/SelectionScans.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2024/exercises/SelectionScans.ipynb) cells 36-48 (2024-08-22) — the same scan, older

</details>

**40 · [`selection_scans_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_scans_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_scans_animal.html)  
Fst and PBS selection scan on wildlife. **Wildebeest**: 24 individuals, blue and black. **Called genotypes** (VCF).  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day4/SelectionScans_22nd.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day4/SelectionScans_22nd.ipynb) cells 38-50 (2026-08-22)

Replaces:

- [`kenya2026/exercises/post_course/day4_afternoon_selection_scans.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day4_afternoon_selection_scans.ipynb) (2026-08-25)

</details>

**41 · [`selection_maize.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_maize.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/selection/selection_maize.html)  
Tajima's D and PBS for domestication and highland adaptation. **Maize** and its wild relative parviglumis, 50 samples. **Low depth**: genotype likelihoods (ANGSD/SAF).  
*From [`summer2025/exercises/Day4_SelectionPopGen2025.ipynb`](https://github.com/popgenDK/courses/blob/main/summer2025/exercises/Day4_SelectionPopGen2025.ipynb) (2025-08-06) — the only copy*


### Relatedness and genetic diversity — `relatedness_diversity/`

**42 · [`fst_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/fst_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/fst_animal.html)  
Pairwise Fst from called genotypes with plink2, then from genotype likelihoods with SAF. **Blue and black wildebeest** (95 individuals, 9 groups), then **Greenland reindeer** (3 populations). **Both**: plink on called genotypes, then SAF on genotype likelihoods.  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`kenya2026/exercises/post_course/day4_morning_fst.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day4_morning_fst.ipynb) (2026-08-25)

Replaces:

- [`kenya2026/exercises/Day4/Fst_Kenya2026.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day4/Fst_Kenya2026.ipynb) (2026-08-23)
- [`kenya2024/exercises/day3_PopulationStructure/Day3_Fst_RH.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day3_PopulationStructure/Day3_Fst_RH.ipynb) (2024-08-09)

</details>

**64 · [`relatedness_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/relatedness_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/relatedness_human.html)  
Why relatedness matters before a GWAS: IBD sharing with plink --genome, a Z1-vs-Z0 plot and a PI_HAT histogram. No close relatives here, but 54 pairs come back nan - all involving individuals missing half their genotypes. **Human**: the GWAS case/control cohort. **Called genotypes**.  
*From written new (2026-09-16), by request, using the data and the `plink --genome` analysis from #46 `gwas_intro_human.ipynb`. — the only copy*

**43 · [`relatedness_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/relatedness_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/relatedness_animal.html)  
A full relatedness pipeline: LD pruning, KING-robust kinship in plink2, a 2D SFS with ANGSD, then relateAdmix for admixture-aware relatedness. **Reindeer**. **Both**: called genotypes for KING, genotype likelihoods for the SFS.  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`kenya2026/exercises/Day5/Related.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day5/Related.ipynb) (2026-08-23)

Replaces:

- [`kenya2026/exercises/post_course/day5_morning_relatedness.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day5_morning_relatedness.ipynb) (2026-08-25)

</details>

**45 · [`heterozygosity_roh_animal.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/heterozygosity_roh_animal.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/relatedness_diversity/heterozygosity_roh_animal.html)  
Heterozygosity and runs of homozygosity. **Blue wildebeest**. **Called genotypes**.  
<details><summary>Where it comes from — replaces 4 older copies</summary>

From [`kenya2026/exercises/post_course/day5_morning_heterozygosity_roh.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/post_course/day5_morning_heterozygosity_roh.ipynb) (2026-08-25)

Replaces:

- [`kenya2026/exercises/Day5/Day5_GeneticDiversity.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2026/exercises/Day5/Day5_GeneticDiversity.ipynb) (2026-08-23)
- [`kenya2024/exercises/day2/Day2_Inbreeding_ROH.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day2/Day2_Inbreeding_ROH.ipynb) (2024-08-08)
- [`kenya2024/exercises/day2/Inbreeding_ROH.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/day2/Inbreeding_ROH.ipynb) (2024-08-07)
- [`kenya2024/exercises/Day1_GeneticDiversity/Day1_GeneticDiversity.ipynb`](https://github.com/popgenDK/courses/blob/main/kenya2024/exercises/Day1_GeneticDiversity/Day1_GeneticDiversity.ipynb) (2024-08-06)

</details>


### GWAS and human complex traits — `gwas/`

**46 · [`gwas_intro_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/gwas_intro_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/gwas_intro_human.html)  
Introduction to GWAS: association testing and QC. **Human**. **Called genotypes** (plink).  
<details><summary>Where it comes from — replaces 1 older copy</summary>

From [`novCourse2024/1GWASIntro.ipynb`](https://github.com/popgenDK/courses/blob/main/novCourse2024/1GWASIntro.ipynb) (2025-07-08)

Replaces:

- [`bgi23/04.GWASintro_2023_SAIGE.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/04.GWASintro_2023_SAIGE.ipynb) (2023-11-09)

</details>

**47 · [`gwas_sumstats_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/gwas_sumstats_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/gwas_sumstats_human.html)  
Working with GWAS summary statistics. **Human**. **Summary statistics only** — no individual genotypes.  
*From [`novCourse2024/2GWASsumstats.ipynb`](https://github.com/popgenDK/courses/blob/main/novCourse2024/2GWASsumstats.ipynb) (2025-07-08) — the only copy*

**48 · [`gwas_analysis_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/gwas_analysis_human.ipynb)** — solution _not rendered_  
A GWAS end to end: first pass, QQ plot, QC, PCA and Tracy-Widom, regenie mixed model with and without PCs, gene-based tests, power. **Human**: simulated UK Biobank, standing height. **Called genotypes**. &#9888; Data missing &mdash; do not edit further.  
*From [`chinaCourse2025/Day3_GWAS_Analysis_2025_Morning.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day3_GWAS_Analysis_2025_Morning.ipynb) (2025-07-26) — the only copy*

**49 · [`gene_based_testing_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/gene_based_testing_human.ipynb)** — solution _not rendered_  
Rare-variant gene-based testing with regenie: annotation, set list and mask files, then SKAT and ACAT. **Human**: whole-exome, Charcot-Marie-Tooth disease (binary trait). **Called genotypes**. &#9888; Data missing &mdash; do not edit further.  
*From [`chinaCourse2025/Day3_Gene_Based_Testing_2025_Afternoon.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day3_Gene_Based_Testing_2025_Afternoon.ipynb) (2025-07-26) — the only copy*

**50 · [`wes_famdiab_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/wes_famdiab_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/wes_famdiab_human.html)  
A trio with MODY: read the VCF by eye, then annotate with Ensembl VEP. **Human**: 3 individuals, *HNF1A*. **Called genotypes**. &#9888; Data missing &mdash; added verbatim, unedited.  
*From [`novCourse2024/3WESfamdiab.ipynb`](https://github.com/popgenDK/courses/blob/main/novCourse2024/3WESfamdiab.ipynb) (2025-07-08) — the only copy*

**51 · [`wes_fh_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/wes_fh_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/wes_fh_human.html)  
Rare recessive familial hypercholesterolemia: filter an exome on gnomAD frequency and predicted consequence. **Human**: a family of four with second-cousin parents. **Called genotypes** (annotated whole-exome).  
*From [`novCourse2024/4WESfh.ipynb`](https://github.com/popgenDK/courses/blob/main/novCourse2024/4WESfh.ipynb) (2025-07-08) — the only copy*

**52 · [`wes_diabetes_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/wes_diabetes_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/wes_diabetes_human.html)  
Four unrelated diabetes patients, four different genes &mdash; genetic heterogeneity, and why a GWAS finds none of them. **Human**: 4 unrelated individuals. **Called genotypes** (annotated whole-exome).  
*From [`novCourse2024/5WESdiab.ipynb`](https://github.com/popgenDK/courses/blob/main/novCourse2024/5WESdiab.ipynb) (2025-07-08) — the only copy*

**53 · [`heritability_ldscore_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/heritability_ldscore_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/heritability_ldscore_human.html)  
Heritability twice over: GCTA REML on a genetic relationship matrix, and LD score regression on summary statistics. **Human**: a simulated family cohort (**called genotypes**) and Biobank Japan HDL/LDL (**summary statistics**, East Asian).  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`chinacourse2026/Day5_Morning_heribilty_and_ldscore.ipynb`](https://github.com/popgenDK/courses/blob/main/chinacourse2026/Day5_Morning_heribilty_and_ldscore.ipynb) (2026-08-07)

Replaces:

- [`chinaCourse2025/Day5_heritability_exercise.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day5_heritability_exercise.ipynb) (2025-07-26)
- [`chinaCourse2025/Day5_Afternoon_Genetic_correlation_Partitioned_Heritability.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day5_Afternoon_Genetic_correlation_Partitioned_Heritability.ipynb) (2025-07-26)

</details>

**54 · [`prs_height_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/prs_height_human.ipynb)** — solution _not rendered_  
Three PGS Catalog scores (height, IGF-1, birth weight) computed with plink and compared with measured height and BMI. **Human**: UK Biobank, European ancestry. **Called genotypes** (array, not imputed). &#9888; Data missing &mdash; added verbatim, unedited.  
*From [`chinaCourse2025/Day6_Morning1_PRS_height_pipeline.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day6_Morning1_PRS_height_pipeline.ipynb) (2025-07-26) — the only copy*

**55 · [`mendelian_randomization_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/mendelian_randomization_human.ipynb)** — [solution](https://htmlpreview.github.io/?https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/mendelian_randomization_human.html)  
Two-sample Mendelian randomization of BMI on coronary heart disease: harmonisation, the estimators, F-statistic, pleiotropy and heterogeneity. **Human**, European: GIANT BMI and CARDIoGRAMplusC4D CHD. **Summary statistics** only.  
<details><summary>Where it comes from — replaces 2 older copies</summary>

From [`chinaCourse2025/Day6_Morning2_MR.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day6_Morning2_MR.ipynb) (2025-07-26)

Replaces:

- [`bgi23/MR-exercise.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/MR-exercise.ipynb) (2023-11-10)
- [`bgi23/MR.real_data.exercise.ipynb`](https://github.com/popgenDK/courses/blob/main/bgi23/MR.real_data.exercise.ipynb) (2023-11-10)

</details>

**56 · [`proteomics_mr_human.ipynb`](https://github.com/popgenDK/courses/blob/main/current_exercises/gwas/proteomics_mr_human.ipynb)** — solution _not rendered_  
Proteome-wide MR of plasma proteins on ischaemic stroke, cis-pQTLs, colocalization and a UK Biobank cohort. **Human**, European. **Summary statistics** plus **called genotypes**. &#9888; Data missing &mdash; added verbatim, unedited.  
*From [`chinaCourse2025/Day6_Afternoon_Proteomics_MR.ipynb`](https://github.com/popgenDK/courses/blob/main/chinaCourse2025/Day6_Afternoon_Proteomics_MR.ipynb) (2025-07-26) — the only copy*


## Status

**58 of 58 exercises built.** Numbering runs to 57 because two pairs are each one
exercise listed twice (#6+7 `haplotype_frequencies`, #21+22 `pca_bonus_animal`) and
#17 (SHAPEIT phasing) and #8 (motif discovery) were removed from scope.

| | |
|---|---|
| Data consolidation | ✅ complete — `/course/data/` holds every dataset the exercises need |
| Per-exercise data folders | most exercises have their own folder under `data/current_data/`; a few still read from the bulk course folders |
| Solution HTML | **47 of 51 rendered**; the 4 without one have no data on this server |

## Data

Exercise data lives under **`data/current_data/`** (`data` is a symlink in the
repo root pointing at `/course/data/`). Each built exercise has its own folder
there — `NGSintro/`, `NGSinference/`, `PCA/` — alongside shared support such as
`geneticMap/` and `scripts/`, the helper R libraries that exercises `source()`.
The remaining folders are raw material for exercises not yet built, and get
curated as each one is done.

Nothing outside `data/` should be read by an exercise. The symlink itself is not
in git — it resolves only on the popgen server.


[shinylive-stats_normal]: https://shinylive.io/r/app/#code=NobwRAdghgtgpmAXAAjFADugdAJTAGlQGMB7CAFzgqVQGJkBlACwEsIBPZDdZOAD3QAbEgCc2Ac2TkmcZBFEwog5ABMWAZ3JiARgFdyLMlgA6EegBkouiESbIW5exGQ5kAdwdNEp+smTqSXREiOAAKYzAmcnJ0dUQAeniRKDcscU9dPXU4YLJKCixSGHj0EnRxKgARAGl40iDs9SS4ADMmmSgVJsU2OqCRKnIAfX4cog04JvVWDnjNKHJ1IfkRRUFcCIBKU1NBFm1kkXZQ6bZ2bYgfZErWtlkAVQBJZBbRLkw9ogXDS4hdFmQAB4ALQvQT-FQABSgFXCvz89AAgh8WF8DGQpA5BHBTH4DORsdCIHBBOEwAA5BRKa4aLT7fQ-ZAAMWsRHREHUW3wO2cyHoDBYKjg2igIncnic6H06heb2kshWa1UtJ0DIx6FFsDglBEnN56kFwtFlnYgXIcL8fgNQpFIiJJItlr8f3gYiIjwgUvNEXgUEuBGQEQAsnA-chwroAAwAZm0RE23gDADclLpZABeZCRwiaOA8TORrAARk23N5TpdOVRHq9ZPUKgihAiDHIfpUopU1zgSZY3wxEZjRGjCcbyBT4IzyCLhBgbGQBeLOco+azxdLuKdcl0rurnv0ZL4o4iADVU7JQnxNrKxSpFVwIJ30IrE4Rx2n51ml3mP4WS2XN1uO7unu3pgAAjkeYAAIq6H6BjYuGYFXq8Ypgc+o5vpOhYAKwznOBYzlAfAftO-jLj+WCRiWG5+Ou8KWvyTAkG4XDIEIJCOCQLRSDIyqaKq7L3p28rIAM6i6IIizIFxPGyC0rLsnqTo9BA9qkjRlrseQADy+i1hEaiaJCwjkFyGl+EmOQigYMAACr8Dpen7gZio4JMEmmWAdEAZZIjWSwdkObpMTOWAT4KG54mSWZ5aWr5-mBXwjkhaBaERe50VeRpFy0aYFxXDcLR3P4OS+cgwjpEQpjZCIZUgi8Ck-KEbBeoQZpeleIA8giLihmyLCWbwAhiQaGLkCQyBfIIRASQssgiYq1KGXSejshui3KPVAxQP1lmhF1sW+s49UtfoAAkR0afWQKgqd5BnfWGnESdIFnYesVgTdkrnRB9F+HsmiOk6t4KB+IOrBeM6hhAOYqN5m7hasH6IzAkPIEdsPw06aVI5mOOo2BUN+pj2UbgAvvlvL0MFXqyaJGVSTJC1Uso8k2IpG7ted4MwJFHlfQMD45PZSX7RpGq5mSPNkoQd1vU2BCjnLGOBorAZy-WCtXumo4baEmxnTzOXIBT3XSU590o3zkkC1QQoiCL5oHU6EuUGSKMy99918Ar3Lq69KuNkrr2a6r2u6yz+tnSjxum-RXP3fj1uOFtdvCw5Yuxa7YQRPjnty4Tqt+7LAfQ77wdeg9Kha-OEerEoUf47HlM9cZHF0xtfErWq8cW4btJtynoKC-bg+Z8p0Nfcr0NXZ2L2V49sXPaC2RgaER3IKCAAsABUocbwA1Mgu+h9iEDiNIWBmiRkaRsbfi2CQ2QQEMy9e29Zs9QAwttlCsVp4sTJo2lj7dGZd-Bw0IOQdg6BJwREEKOUgyhMwRG0BOUcfBBBQG0B+CIh4AzsCwTglBYAbgcgcOwUcKlcEUhZjSfi9I1pZVitgvYxJQhJg-I-Z+r9CBIJoQMBsAZBBuE7JmAATPfNiJA2CLFCNwqgvDVCKnkUxHhoCVb1lLJNEgyDVaCNHOgWwJEACczdTBkx5KcDgyJ0ChH+B+f4OZSo5A-DVXyFwwBkwALpAA