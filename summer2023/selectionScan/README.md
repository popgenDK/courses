
# Exercises in selection scan

set the following

```bash
genomes=/home/rute/2023/SELECTION/REFERENCEgenome/
```

## Exercise 1. Lactase locus

Choose your favorite genome browser and find the lactase (Gene: LCT ) locus in the human genome. 
The figures in the presentation use ENSEMBL https://www.ensembl.org/index.html
You can also use UCSC genome browser https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu

### Find the lactase gene annotation in the gff file

Compare the genomic context of the lactase gene in the human with other species
Use the "Orthologs" link on the left-hand panel to find the a list of species that
you can compare against

find the lactase gene annotation by name in the gff file
```bash
zcat $genomes/Homo_sapiens.GRCh38.110.chromosome.2.gff3.gz | grep LCT
```

## Exercise 2. SNP information

Find information about the SNP rs4988235 in dbSNP
Extract the gff features that overlap with the SNP position using “bedtools intersect” 
(you will need a bed file with the SNP genomic coordinates)

I have prepared the bed file for you, check it out (remember bed files are zero-based)
```bash
more $genomes/rs4988235.bed
```

```bash
bedtools intersect  -a $genomes/Homo_sapiens.GRCh38.110.chromosome.2.gff3.gz -b $genomes/rs4988235.bed
```

## Exercise 4. Describe windows X in population B: estimate thetas and the 1dSFS

## Exercise 5. Estimate the 2dSFS for all 3 pairwise comparisons

Follow the instruction on the slides


## Exercise 6. Calculate the diversity and the PBS statistic for maize genes targeted by domestication  

you will be executing scripts but please 
take your time to read the options to better understand 
what you are doing in each step
 we will use angsd to calculate theta values, Tajima's D, Fst and PBS
more info:
http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests
http://www.popgen.dk/angsd/index.php/Fst

the paths to reference file and bamfiles are already inside the scripts
but also here in case you want to try out other options 

```bash
REF=/home/rute/2023/SELECTION/REFERENCEgenome/maizegenome-chl-mito-plasm.fasta
ANC=$REF
misc=/emily/program/albrecht/angsd/misc
```

lists of bamFiles
```bash

wild=/home/rute/2023/SELECTION/EXERCISES/wild.fileList
p750=/home/rute/2023/SELECTION/EXERCISES/750years.fileList
p2000=/home/rute/2023/SELECTION/EXERCISES/2000years.fileList
modern=/home/rute/2023/SELECTION/EXERCISES/modern.fileList
```

we will focus on individual coding regions (this is capture data)
the information on their genomic coordinates will be provide to angsd 
with region files (1-based): chr:start-end

```bash

ZAGL1=/home/rute/2023/SELECTION/COORS/ZAGL1.coors
DEHYD1A=/home/rute/2023/SELECTION/COORS/DEHYD1A.coors
SUGARY1=/home/rute/2023/SELECTION/COORS/SUGARY1.coors
AE1=/home/rute/2023/SELECTION/COORS/AE1.coors
50genes=/home/rute/2023/SELECTION/COORS/random50.coors
```

### Step 1: Finding a 'global estimate' of the SFS for each population

```bash
1-reconstruct-sfs.sh
```

### Step 2: Calculate the thetas for each site

```bash

2-calculate-thetas
```

### Step 3: calculate 2D-sfs and Fst for all population pairs

```bash

3-calculate-folded-fst
```
#### Step 4:
```bash
4-calculate-PBS
```

this will run the full analysis for gene ZAGL1
you need to provide the location of the coordinate file and a prefix for the output
```bash

./run-get-sfs-thetas-2dsfs-fst-PBS.sh $ZAGL1 ZAGL1
```
 repeat for the other genes for which we provide coordinate files

for the genes involved in the starch pathway (SUGARY1 and AE1)
plot the values of nucleotide diversity for each population per gene
check the *pestPG for the tP values
(remember to divide tP by the Nsites to get the estimate per site)

check the PBS values in the *pbs.log file
pop1=wild
pop2=2000years
pop3=750years

OPTIONAL: try to draw tree for the individual genes, negative values correspond to zero
if you use the newick format
(wild:0.284891,p2k:0,p750:0.222320);
you can use (choose radial for the depiction)
http://www.trex.uqam.ca/index.php?action=newick&project=trex

#OPTIONAL: combine the values of Tajima's D and PBS (750 and 2K) in one plot
