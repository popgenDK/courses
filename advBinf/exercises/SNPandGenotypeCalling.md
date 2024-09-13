# title NGS: SNP calling, genoptype calling, allele frequencies, Imputation and EM algorithm


We will work on 33 European samples where we have reduced the genome so there is very little data. 

We will be looking at the default behaviour for 3 programs

| tool | SAMtools | GATK | ANGSD |
| --- | --- | --- | --- |
| Genotype likelihood | tries to model error dependenties | simple | user specified |
| SNP caller | uses the SFS as prior | uses SFS as prior | likelihood ratio test |
| Genotype caller | MAF | ML | MAF prior |

<br>
All the programs have additional filters and additional differences. 

<contents>

### set some paths 
```
# test if software is installed
which angsd
which samtools
which bcftools
which vcftools

##must be done every time you open a new terminal
## main folder
FOLDER=/course/advBinf/calling

#beagle
BEAGLE=$FOLDER/prog/beagle.jar

#reference genome for human
REF=$FOLDER/hg19.fa.gz

# a folder with bam files
BAMFOLDER=$FOLDER/smallBam
```


### mpileup data for multiple individuals


We will require the full hg19.fa.gz reference and the bam files located in the subdir called **smallBam**/.

To allow for random access we first need to index the bam files.

```
## do not run this is done for you to same time
#for i in $BAMFOLDER/smallNA*.bam;do samtools index $i;done
```

Now you make a file list with the newly indexed bamfiles.

```
ls $BAMFOLDER/smallNA*.bam >bams.list
```

 - how many individuals do you have?

```
samtools mpileup -b bams.list  | less -S
```

 - view the data in mpileup format and identify the columns


And you should index the fasta reference
```
## do not run. This is done for you
#samtools faidx $REF
```


## VCF files
In this exercise you will try to generate the widely used **vcf** formatted files. We will use

 - angsd
 - SAMtools
 - GATK
 - vcftools

### SAMtools
   SAMtools outputs a binary version of vcf files, called bcf files. To get a the vcf equivalent you shou0ld pipe the data to bcftools which is located within the samtools/bcftools subfolder.

Samtools can do single sample genotyping/SNP calling however, it calculates the genotype likelihoods a little different than the simple GATK model

```
bcftools mpileup -f $REF -Ou -b bams.list  -r 1 | bcftools call -V indels -mv -Ov -o sam.vcf
```


The above command used Samtools to calculate genotype likelihood for chr1 and use bcftools to call variable SNP sites based on the reference genome

[[VCFfile.pdf][Overview of VCF file]]

 - Look at the sam.vcf file and identify some of the relevant columns
 -  how many variable sites does it think we have?

### GATK
  GATK was much to slow, so I have precomputed the gatk.vcf.gz file for you.
GATK uses a prior in order to call SNP sites but (not surprisingly) used the GATK genotype likelhood model

 GATK also requires some more steps. Here is the list of commands used to generate the file

```
## DO NOT RUN. SO SLOWWWWWWWWW
##first gatk doesn't support .gz compressed references
#gunzip ref/hg19.fa.gz

##GATK requirs a dictionary for the reference
## remember to change the path to the program
#java -jar CreateSequenceDictionary.jar R= ref/hg19.fa O=ref/hg19.dict
## now we can run gatk
#./gatk  -R ref/hg19.fa -T UnifiedGenotyper -I bams.list -L 1 -nt 4 -o gatk.vcf
## please zip the file again
#gzip ref/hg19.fa
```

copy the GATK output file instead
```
cp $FOLDER/gatk.vcf .
cp $FOLDER/smallChr1.geno .
```

 - Look at the gatk.vcf file, how many variable sites does it have?

### ANGSD
  
  We have another program called ANGSD, which does not output vcf files directly. But we can generate a ***mafs.gz*** (or just mafs depending the angsd version) file which contains information about estimated major/minor genomic position etc.

```
angsd -bam bams.list -r 1: -SNP_pval 0.001 -out angsd -doMaf 2 -doMajorMinor 1 -GL 1 
```

This commands can be read as: We want to run the analysis based on 'bams.list' and we limit the analysis to chromosome '1'. We are not interested in all sites, but only those sites that are variable with a likelihood ratio test with a p-value of 0.001. Our output files should be prefixed with angsd. We want to estimate the allele frequency, and that requires that we also find the major and minor allele. And we base all analysis on the Samtools model of genotype likelihoods. 

Here we use the allele frequency in order to call SNPs (as in the slides) and no reference information

The output file is here called 'angsd.mafs.gz'. 
```
chromo	position	major	minor	unknownEM	pu-EM	nInd
1	14000023	C	A	0.051607	5.401180e-04	28
1	14000202	G	A	0.270108	0.000000e+00	27
1	14000873	G	A	0.305992	0.000000e+00	33
1	14001018	T	C	0.296055	0.000000e+00	30
1	14001867	A	G	0.301093	0.000000e+00	33
1	14002342	C	T	0.042232	4.323473e-07	31
1	14002422	A	T	0.294218	0.000000e+00	29
1	14002474	T	C	0.060752	0.000000e+00	30
1	14003581	C	T	0.256778	0.000000e+00	29
```
unknownEM is the estimate allele frequency and pu-EM is the p-value, and nInd are the number of individuals with information (>0 depth)


 -  view the file, how many variable sites do we find using angsd?



### Comparing the results (SNP-discovery)
 - This requires R and the VennDiagram package.
Let us examine the difference between the three different approaches for SNP-calling. This we do by looking at the overlap. Start R.

```
##run in R
#read in data
sam <- read.table("sam.vcf")[,2]
gatk <- read.table("gatk.vcf")[,2]
angsd <- read.table("angsd.mafs.gz",he=T)[,2]
```
Our three variables now contain the positions. Now find the overlap between the three arrays.
To find the overlap between SAMtools and gatk we can do
```
##run in R
sam[sam %in% gatk]
```
This will give of the positions in sam, that exists in gatk.
We want to do this for all combinations, and visualize the different possible overlaps in call SNPs. A Venn diagram is suited for this.

```
##run in R
library(VennDiagram)
#R help function to plot vennDiagram
myVenn3 <- function(x,y,z,...){
  a1 <- length(x)
  a2 <- length(y)
  a3 <- length(z)
  n12 <- sum(x%in% y)
  n23 <- sum(y%in% z)
  n13 <- sum(x %in% z)
  n123 <- sum(x[x %in%y ] %in% z)
  cat(a1,",",a2,",",a3,"\t",n12,",",n23,",",n13,",",n123,"\n")
  draw.triple.venn(area1=a1,area2=a2,area3=a3,n12=n12,n23=n23,n13=n13,n123=n123,...)
}

myVenn3(sam,gatk,angsd,fill=1:3,category=c("SAMtools","GATK","ANGSD"))
```

keep R open


 -  How is this to be interpreted? What are the overlaps between the different approaches?

continue in R
```
##run in R
#read in known genotypes from the individuals (one is missing sorry about that)
g<-read.table("smallChr1.geno")
trueSNPs <- sapply(strsplit(names(g),"_"),function(x)x[2])

#see how many they got correct
table(sam%in%trueSNPs)
table(angsd%in%trueSNPs)
table(gatk%in%trueSNPs)
```

 - What is the strength and weaknesses of the 3 methods

## Bonus exercise
There are 5 sites that both GATK and SAMtools call as SNPs, that ANGSD doesn't find.
 -   find these 5 sites. 
 - figure out why these stand out.

### Comparing the results (Genotype-calling)
Continuing our last example. We therefore have two vcf files 'gatk.vcf' and 'sam.vcf'. Use vcftools to extract the genotypes to a more easily accessible format.

```
vcftools --vcf gatk.vcf --012 --out gatk
vcftools --vcf sam.vcf --012 --out sam
```
ignored warnings about the VCF headers.

 -  What are the output files? What do they contain? (-1 is missing)


We now call genotypes using ANGSD, and by simply calling the genotype that has the highest likelihood. We see that genotypes in the above vcf format is identified by the counts of non-reference. The default behaviour of angsd is to estimate the major/minor based on GL. But we can force it to use the reference as major by using -doMajorMinor 4 and supplying the reference (just as in done i GATK and Samtools).
```
angsd -bam bams.list -SNP_pval 0.001 -doMaf 2 -doMajorMinor 4  -r 1: \
-doGeno 2 -doPost 2  -GL 1 -out angsd2 -ref $REF
```


# Imputation

use the same individuals as the previous exercises

```
ls $BAMFOLDER/smallNA*.bam >bams.list
```

## preparing the information needed for imputation
Using angsd we will call variable sites (-SNP_pval ) based on allele frequencies (-doMaf) based on the reference allele (-ref) and the inferred alternative allele (-doMajorMinor). Based on these sites we will out put the genotype likelihoods based on the samtools model (-GL 1) in the BEAGLE software format (-doGlf 2):

```
angsd -bam bams.list  -SNP_pval 0.001 -doMaf 2 -doMajorMinor 4  -r 1: \
 -doGlf 2  -GL 1 -out angsd  -ref $REF
```

look in the BEAGLE genotype likelihood format using for example 
```
zcat  angsd.beagle.gz | less -S
```

 - why is there 3 columns for each individuals

### performing imputation on the genotype likelihoods


Run the imputation based on the genotype likelihoods 

```
java -jar $BEAGLE like=angsd.beagle.gz out=imputation
```

### comparing performance of the different methods

The 33 individuals sequenced as part of the 1000 genomes project was also part of the HapMap project. Therefore, we have accurate information about many of the genotypes from these individuals.

copy the hapmap data to your folder and uncompress it 

```
cp $FOLDER/hapmap_CEU_r23a_filteredHg19Chr1.tar.gz .
tar -xf hapmap_CEU_r23a_filteredHg19Chr1.tar.gz
```


Lets also compare with the genotypes from the previous exercise. There is a bit of data merging so don't worry to much about the R code and just paste it all

```
#open R



library(snpStats)


#read in the genotype data. Change the path to where you unpacked the genotypes
pl<-read.plink("hapmap_CEU_r23a_filteredHg19Chr1")
bim<-read.table("hapmap_CEU_r23a_filteredHg19Chr1.bim",as.is=T)


#read in the called genotypes from angsd, samtools and GATK just as before
gatk.pos <- read.table("gatk.012.pos")[,2]
sam.pos <-  read.table("sam.012.pos")[,2]
angsd <-  read.table("angsd2.geno.gz")

#keep only overlapping sites
ta<-table(c(angsd[,2],gatk.pos,sam.pos))
pos <- as.numeric(names(ta[ta==3]))
angsd<-as.data.frame(t(angsd[angsd[,2] %in% pos,][,-c(1,2)]))
gatk<-read.table("gatk.012")[,-1][gatk.pos %in% pos]
sam<-read.table("sam.012")[,-1][sam.pos %in% pos]

#find the overlab with the hapmap data
keep<-pos%in%bim[,4]
int<-which(bim[,4]%in%pos)
pl <- pl$genotypes
indNames<-rownames(pl) # individual names in HapMap
indNamesBam<-sapply(strsplit(sub("small","",basename(scan("bams.list",what="theFck"))),".m"),function(x)x[1]) #individual names in the 33 1000genome 

#convert the genotype data into integers
genoHap<-as.integer(pl[,int])
genoHap<-matrix(3-genoHap,nrow=60,ncol=length(int))
rownames(genoHap)<-indNames
genoHap[genoHap==3]<-NA

# match the strands
mafs<-read.table("angsd.mafs.gz",as.is=T,head=T) #you might have the remove the .gz
mafs<-mafs[mafs$position%in%pos,]
refStrand<-bim[int,5]==mafs$ref[keep]
refStrand<-bim[int,5]==mafs$ref[keep]
genoHap[rep(refStrand,each=60)]<-2-genoHap[rep(refStrand,each=60)]

#estimate the mean concordance rate
res <- sapply(list(sam=sam,gatk=gatk,angsd=angsd),function(x) 
mean(genoHap[indNamesBam,]==as.matrix(x[,keep]),na.rm=T))


##don't close R
```
 - the above calculations assumed that the missing genotypes are discordant.

```
# Estimate the concordance rate without missing data
angsdNA<-as.matrix(angsd)
angsdNA[angsd==-1]<-NA
samNA<-as.matrix(sam)
samNA[sam==-1]<-NA
gatkNA<-as.matrix(gatk)
gatkNA[gatkNA==-1]<-NA
resMis <- sapply(list(sam=samNA,gatk=gatkNA,angsd=angsdNA),function(x) 
  mean(is.na(x[,keep]),na.rm=T))
resNA <- 
sapply(list(sam=sam,gatk=gatk,angsd=angsd,gatkNA=gatkNA,angsdNA=angsdNA,samNA=samNA),function(x) 
  mean(genoHap[indNamesBam,]==as.matrix(x[,keep]),na.rm=T))




cat("Missing genotype ratet:\n")
print(resMis)

cat("\nDiscordance rate assuming missing is discordant:\n")
print(res)

cat("\nDiscordance rate when ignoring missing genotypes:\n")
print(resNA)
##don't close R

```


 - which of the genotyping methods performs the best in the above analysis?
 - How could the comparison of methods be improved

include the imputation results in the comparisons
```
##run in R
#read in the genotype probabilities
gprobsDat<-read.table("imputation.angsd.beagle.gz.gprobs.gz",head=T,as.is=T)
impuPos<-as.integer(sub("1_","",gprobsDat[,1])) #position

#keep only overlapping positions and convert to a matrix
gprob<-as.matrix(gprobsDat[impuPos %in% pos[keep],-c(1:3)])

#call the genotype with the highest probability
funn<-function(x)
  apply(gprob[,(x-1)*3+1:3],1,which.max)-1


imputa <- t(sapply(1:33, funn))

#get the mean concordance rate 
mean(genoHap[indNamesBam,]==imputa,na.rm=T)
```

 - why does this method perform so much better ?
 - how could you increase the accuracy even further ?
 - when will this method not increase performance?
 - How you think the performance difference correlates with sample allele frequency?



