# Exercises for the selection session
We will use the following programs from angsd

```bash
which angsd
which realSFS
which thetaStat

```


Set the following paths

```bash
thePath=/course/popgen23/rute/selection
REF=$thePath/REFERENCEgenome/maizegenome-chl-mito-plasm.fasta
```

BamFiles for the exercises with the maize data

```bash

wild=$thePath/EXERCISES/wild.fileList
p750=$thePath/EXERCISES/750years.fileList
p2000=$thePath/EXERCISES/2000years.fileList
modern=$thePath/EXERCISES/modern.fileList
```

We will focus on individual coding regions (this is capture data)
the information on their genomic coordinates will be provide to angsd 
with region files (1-based): chr:start-end

```bash

ZAGL1=$thePath/COORS/ZAGL1.coors
DEHYD1A=$thePath/COORS/DEHYD1A.coors
SUGARY1=$thePath/COORS/SUGARY1.coors
AE1=$thePath/COORS/AE1.coors
genes50=$thePath/COORS/random50.coors
```


 make a directory to work in and copy some scripts to it
```bash
mkdir -p ~/selection
cd ~/selection

cp $thePath/SCRIPTS/*sh .
#print the names of the scripts you copied
ls *sh
```

You will be executing scripts but please take your time to read the options to better understand what you are doing in each step.

We will use angsd to calculate theta values, Tajima's D, Fst and PBS. More information at:

http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests

http://www.popgen.dk/angsd/index.php/Fst

#### Exercises 1 and 2: 40 minutes

## Exercise 1. Describe windows X in populations A, B and C

#### Follow the instruction on the slides (starting on slide 43)

In this exercise you will estimate thetas and the 1dSFS for population B and the 2dSFS for all 3 pairwise comparisons.
#### First calculate watersons theta for population B (Slide 44)
The "number of individuals" for diploid organisms is 2N (n=6 for the example on slide 43). You can calculate 
  $$, \sum_{k=1}^{2N-1} 1/k $$  in R as
```bash
# open R for before running code
sum(1/1:5)
```
#### Calculate Tajimas theta ( $\Pi$ ) for population B (Slide 45)
For each pair of haplotype calculate the pairwise difference and fill in the table. Then calculate $\Pi$

#### calculate and plot the SFS for population B (Slide 46)
Examples for plotting the SFS for populaton A
```bash
# open R
counts <- c(6, 3, 2, 1, 0, 0, 0)
barplot(counts,names=0:6,col="hotpink",ylab="Number of sites",xlab="Number of derieved alleles")
```

####  Calculate the 2D SFS for population pairs A,C and B,C(Slide 47)
For each position in the genome find the number of derived alleles for the pair of population and increment the table by one. 




## Exercise 2.  Calculate the diversity in the maize genes from the starch pathway

Make sure you set up the paths to the bamfiles and coordinate files (see above)

The following 2 steps will be executed using the script run-get-sfs-thetas.sh. You can also execute the commands one by one if you want to follow each step independently:

- Step 1: Finding a 'global estimate' of the SFS for each population (script: 1-reconstruct-sfs.sh)

- Step 2: Calculate the thetas for each site (script: 2-calculate-thetas)

The following command lines will calculate thetas for 2 genomic regions, one including gene SUGARY1 and the other including gene AE1.

You need to provide the location of the coordinate file and a prefix for the output files as seen below if you want to try it with other genes.

```bash

./run-get-sfs-thetas.sh $SUGARY1 SUGARY1
./run-get-sfs-thetas.sh $AE1 AE1
```

For the genes involved in the starch pathway (SUGARY1 and AE1), plot the values of nucleotide diversity for each population per gene.
To do so you can pull out the nucleotide diversity (tP) values and the number of sites used to calculate it (Nsites) from the *pestPG files. Once you have that the estimate per site can be calculated by dividing tP by the Nsites. So e.g. you can extract the tP and Nsites values for the AEP1 gene for the 750 year old samples using this command:


```bash

cut -f5,14 750years.AE1.thetas.idx.pestPG
```
And then you can get the estimated nucleotide diversity by dividing the two numbers you get out (tP/Nsites). Once you have done this for the AE1 gene for the wild samples (teosinte), the 2000 year old samples, the 750 year old samples, and the modern samples (landraces) you can make the relevant plot for AE1 gene i R using the following command:
```R
vectorwithestimates = c(estWild, est2000y, est750y,estModern)
barplot(vectorwithestimates, names=c("Wild","2000y","750y","Modern"), main="Nucleotide diversity in AE1")
```
Do this for both AE1 and SUGARY1.


## Exercise 3. Calculate the diversity and the PBS statistic for maize genes targeted by domestication  

These 4 steps are run from the run-get-sfs-thetas-2dsfs-fst-PBS.sh script:

- Step 1: Finding a 'global estimate' of the SFS for each population (script: 1-reconstruct-sfs.sh)

- Step 2: Calculate the thetas for each site (script: 2-calculate-thetas)
  
- Step 3: Calculate 2D-sfs and Fst for all population pairs (script: 3-calculate-folded-fst)
  
- Step 4: Calculate PBS for genes of interest (script: 4-calculate-PBS)

The following line will run the full analysis for gene ZAGL1. You need to provide the location of the coordinate file and a prefix for the output.

```bash

./run-get-sfs-thetas-2dsfs-fst-PBS.sh $ZAGL1 ZAGL1
```

Repeat for the other genes for which we provide coordinate files (see above).

Check the PBS values in the *pbs.log file:

- pop1 = wild population 
  
- pop2 = 2000 year old population
  
- pop3 = 750 year old population

OPTIONAL: try to draw tree for the individual genes, negative values correspond to zero

If you use the newick format (the branch length is indicated after the population name)


```bash

(wild:0.284891,p2k:0,p750:0.222320);
```

you can use this tool for visualization (choose radial for the depiction):


http://www.trex.uqam.ca/index.php?action=newick&project=trex


#OPTIONAL: combine the values of Tajima's D and PBS (750 and 2K) in one plot

## Exercise 4. The lactase locus

#### Find the lactase (Gene: LCT ) locus in the human genome using the ENSEMBL browser.  

Go to the ENSEMBL browser: https://www.ensembl.org/index.html. 

(You might be able to use your favorite browser e.g. UCSC genome browser 

https://genome-euro.ucsc.edu/cgi-bin/hgGateway?redirect=manual&source=genome.ucsc.edu 

but for this exercise use ENSEMBL). 

Then search for "LCT". This will give you quite a few options. Click on the one called "LCT (Human Gene)" (the top one). Now try to explore the genomic context (e.g. which other genes are nearby).


#### Compare the genomic context of the lactase gene in the human with other species

Now based on the the above ENSEMBL page for LCT In the ENSEMBL use the "Orthologs" link on the left-hand panel to find the a list of species that you can compare against. 

How many primates with orthologs did you find?

#### Find the lactase gene annotation in the gff file

Let's also try to look up LCT in the annotation file (gff file) that comes along with the genome.

The gene id for the lactase gene is LCT.

```bash

zcat $thePath/REFERENCEgenome/Homo_sapiens.GRCh38.110.chromosome.2.gff3.gz | grep LCT
```


## Exercise 5. Obtaining information about a known SNP

We will try to investigate SNP rs4988235. You can find information about the SNP rs4988235 in

[dbSNP](https://www.ncbi.nlm.nih.gov/snp/rs4988235)


#### Extract the gff features that overlap with the SNP position using “bedtools intersect”  (you will need a bed file with the SNP genomic coordinates)

I have prepared the bed file for you, check it out (remember bed files are zero-based) by first finding the SNP position for the relevant SNP and then looking up everything that overlaps with this posotion by running these commands:

```bash

cat $thePath/REFERENCEgenome/rs4988235.bed
```

```bash

bedtools intersect  -a $thePath/REFERENCEgenome/Homo_sapiens.GRCh38.110.chromosome.2.gff3.gz -b $thePath/REFERENCEgenome/rs4988235.bed
```
- Based on the output which genes overlap with the SNP rs4988235?
