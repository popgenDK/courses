
# Exercises in selection scan

set the following path and make a directory to work in 

```bash
#folder with data and scripts
genomes=/course/popgen23/rute/selection

#see if scripts are there
ls $genomes/SCRIPTS/1-reconstruct-sfs.sh
ls $genomes/SCRIPTS/2-calculate-thetas
ls $genomes/SCRIPTS/3-calculate-folded-fst
ls $genomes/SCRIPTS/4-calculate-PBS
ls $genomes/SCRIPTS/run-get-sfs-thetas-2dsfs-fst-PBS.sh

## check if software is installed
which angsd
which realSFS

#make a folder and enter it. 
mkdir -p ~/selection
cd ~/selection
```

## Exercise 1. Lactase locus


Find the lactase (Gene: LCT ) locus in the human genome using the ENSEMBL browser.  https://www.ensembl.org/index.html. You might be able to use your favoriete browser e.g.  UCSC genome browser https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu but for this exercise use ENSEMBL. 


Compare the genomic context of the lactase gene in the human with other species
In the ENSEMBL browser use the "Orthologs" link on the left-hand panel to find the a list of species that you can compare against. 
 - how many primates with orthologs did you find
### Find the lactase gene annotation in the gff file

find the lactase gene annotation by name in the gff file
```bash
zcat $genomes/REFERENCEgenome/Homo_sapiens.GRCh38.110.chromosome.2.gff3.gz | grep LCT
```
 -  QUESTION
## Exercise 2. SNP information

We will try to investigate SNP rs4988235. You can find information about the SNP rs4988235 in [dbSNP](https://www.ncbi.nlm.nih.gov/snp/rs4988235). 
Extract the gff features that overlap with the SNP position using “bedtools intersect”  (you will need a bed file with the SNP genomic coordinates)
I have prepared the bed file for you, check it out (remember bed files are zero-based)
```bash
cat $genomes/REFERENCEgenome/rs4988235.bed
```
Run bedtools and extract the features in that region
```bash
bedtools intersect  -a $genomes/REFERENCEgenome/Homo_sapiens.GRCh38.110.chromosome.2.gff3.gz -b $genomes/REFERENCEgenome/rs4988235.bed
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
REF=$genome/REFERENCEgenome/maizegenome-chl-mito-plasm.fasta
ANC=$REF
misc=/emily/program/albrecht/angsd/misc
```

lists of bamFiles
```bash

wild=$genomes/EXERCISES/wild.fileList
p750=$genomes/EXERCISES/750years.fileList
p2000=$genomes/EXERCISES/2000years.fileList
modern=$genomes/EXERCISES/modern.fileList
```

we will focus on individual coding regions (this is capture data)
the information on their genomic coordinates will be provide to angsd 
with region files (1-based): chr:start-end

```bash

ZAGL1=$genomes/COORS/ZAGL1.coors
DEHYD1A=$genomes/COORS/DEHYD1A.coors
SUGARY1=$genomes/COORS/SUGARY1.coors
AE1=$genomes/COORS/AE1.coors
genes50=$genomes/COORS/random50.coors
```

### Step 1: Finding a 'global estimate' of the SFS for each population
we will use a script
<details>
<summary> 1-reconstruct-sfs.sh/summary>
REF=/home/rute/2023/SELECTION/REFERENCEgenome/maizegenome-chl-mito-plasm.fasta
ANC=$REF

id=$1 #output suffix
bamList=$2
minInd=$3
rf=$4
outSuffix=$5

out=${id}.${outSuffix} 

misc=/emily/program/albrecht/angsd/misc

angsd -bam $bamList -doSaf 1 -anc $ANC -GL 1 -P 10 -out $out -minQ 20 -minMapQ 30 -minInd $minInd -rf ${rf} 

#Obtain the maximum likelihood estimate of the SFS using the realSFS program found in the misc subfolder.
$misc/realSFS $out.saf.idx -P 10 > $out.sfs

</details>
to run in 

```bash
#set a output name
OUTNAME=WILD_ZAGL1
## run the script with the chosen file list and gene
$genomes/SCRIPTS/1-reconstruct-sfs.sh ~/selection/SAF $enomes/EXERCISES/wild.fileList 3 $ZAGL1 $OUTNAME
```



### Step 2: Calculate the thetas for each site

```bash

$genomes/2-calculate-thetas
```

### Step 3: calculate 2D-sfs and Fst for all population pairs

```bash

$genomes/3-calculate-folded-fst
```
#### Step 4:
```bash
$genomes/4-calculate-PBS
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
