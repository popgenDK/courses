# PCA introduction. 

Run the colab

[Colab exercises](https://colab.research.google.com/drive/1srBVxjr7D_4oVWDBO8WYdnGvYgru49gD?usp=sharing)



 PCA with admixture aware priors 

Let's try to perform PCA analysis on the same 1000 genotype genotype likelihoods that you performed admixture analysis on yesterday.
First let set the path to program and the input file
```
# Make new folder and set path to that folder
mkdir -p ~/structure2
cd ~/structure2

# NB this must be done every time you open a new terminal
ThePath=/course/popgen23/anders/popstructureII

## PCAngsd
which pcangsd

## Genotype likelihood data from previous exercise
ls $ThePath/1000G5pops.inputgl.beagle.gz 

## sample information
ls $ThePath/1000G5pops.pop.info
```


Copy the information file to your directory and view the first 10 lines. 


```
cp $ThePath/1000G5pops.pop.info .
head 1000G5pops.pop.info
```
 - See the number of individuals for each population from the sample file

```
# summaries the fist column
cut -f1 1000G5pops.pop.info | uniq -c
```

 - Count the number of lines in the genotype likelihood file

```
zcat $ThePath/1000G5pops.inputgl.beagle.gz | wc -l 
```

 - How many sites are in your genotype likelihood file?

 - What were the populations included? And how many sites?

## Run PCAngsd to perform PCA
 First let's get a list of the options in PCAngsd

```
pcangsd -h
```

Run PCAngsd on your genotype likelihood data using 5 CPU threads

```
pcangsd -b $ThePath/1000G5pops.inputgl.beagle.gz -o PCANGSD1000G -t 5
```

The program estimates the covariance matrix that can then be used for PCA. look at the output from the program

 - How many significant PCA (see MAP test in output)?

Plot the results in R

```
#open R
# Read covariance matrix estimated by PCAngsd
C <- as.matrix(read.table("~/structure2/PCANGSD1000G.cov"))

# Read population labels for each individuals
pop<-read.table("~/structure2/1000G5pops.pop.info",stringsAsFactors=T)

# Estimate the eigenvectors (principal components) from the covariance matrix
e <- eigen(C)
plot(e$vectors[,c(1,2)],col=pop[,1],xlab="PC1",ylab="PC2")
legend("left",fill=1:5,levels(pop[,1]))
## not close R
```




Compare with the estimated admixture proportions ( from previous exercise)

[plot](http://popgen.dk/albrecht/phdcourse/html/plots/BestK3.png)


```
## in R
x11()
par(mfrow=2:1)
source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")
## read and plot the output from NGSadmix from the Tuesday's exercises
pop<-read.table("/course/popgen23/ida/admixexercise/admixinput/1000G5pops.pop.info",as.is=T)
q<-read.table("/course/popgen23/ida/admixexercise/admixoutput/1000G5popsAdmixK3seed3.qopt")
ord<-orderInds(pop = pop[,1], q=q) # sort indiivduals by population and within populaoitn by admixture proportion

#plot
barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Admixture proportions",main="K=3")
text(sort(tapply(1:nrow(pop),pop[ord,1],mean)),-0.05,unique(pop[ord,1]),xpd=T) # add population labels
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)

## read for K=4
pop<-read.table("/course/popgen23/ida/admixexercise/admixinput/1000G5pops.pop.info",as.is=T)
q<-read.table("/course/popgen23/ida/admixexercise/admixoutput/1000G5popsAdmixK4seed9.qopt")
plot
barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Admixture proportions",main="K=4")
text(sort(tapply(1:nrow(pop),pop[ord,1],mean)),-0.05,unique(pop[ord,1]),xpd=T) # add population labels
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)
```

 - In the PCA plot can you identify the Mexicans with only European ancestry?
 - What about the African American with East Asian ancestry?
 - Based on the PCA would you have the same conclusion as the admixture proportions?


Try the same analysis but without estimating individual allele frequencies. This is the same as using the first iteration of the algorithm. This is similar to using mean imputation which is used by most PCA software (eigensoft, eigenstrat, smartPCA, NGStools, fastPCA, FlashPCA, plink --fast).

```
#close R if you have not done so
pcangsd -b $ThePath/1000G5pops.inputgl.beagle.gz -o PCANGSD1000G_iter0 -t 5 --iter 0
```

Plot the results in R
```
#open R
# Read covariance matrix estimated by PCAngsd
C <- as.matrix(read.table("~/structure2/PCANGSD1000G_iter0.cov"))

# Read population labels for each individuals
pop<-read.table("~/structure2/1000G5pops.pop.info",stringsAsFactors=T)

# Estimate the eigenvectors (Principal components) from the covariance matrix
e <- eigen(C)
plot(e$vectors[,1:2],col=pop[,1],xlab="PC1",ylab="PC2")
legend("top",fill=1:5,levels(pop[,1]))
#close R
```

 - Do you see any difference?
 - Would any of your conclusions change? (compared to the previous PCA plot)


## Converting a PCA into admixture proportions


Let try to use the PCA to infer admixture proportions based on the first 3 principal components. For the optimization we will use a small penalty on the admixture proportions (alpha). (PCAngsd will automatically that it needs 3 Principal components corresponding to 4 populations)
```
pcangsd -b $ThePath/1000G5pops.inputgl.beagle.gz -o PCANGSD1000G -t 5 --admix --admix_alpha 50 
```

Plot the results in R
```
#open R
# Read the admixture proportions estimated from the PCA
source("https://raw.githubusercontent.com/GenisGE/evalAdmix/master/visFuns.R")
q<-read.table("~/structure2/PCANGSD1000G.admix.4.Q")

# Read population labels for each individuals
pop<-read.table("~/structure2/1000G5pops.pop.info",stringsAsFactors=T)

## Order according to population
ord<-orderInds(pop = pop[,1], q=q) # sort indiivduals by population and within populaoitn by admixture proportion


barplot(t(q)[,ord],col=2:10,space=0,border=NA,xlab="Individuals",ylab="Admixture proportions")
text(sort(tapply(1:nrow(pop),pop[ord,1],mean)),-0.05,unique(pop[ord,1]),xpd=T)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=1.2)

## close R
```
 - how does this compare to the NGSadmix analysis?



## Bonus exercise: PCAngsd and selection

For very recent selection we can look within closely related individuals for example with in Europeans

Data:

 -   Genotype likelihoods in Beagle format
 -   ~150k random SNPs with maf > 5%
 -   Four EU populations with ~100 individuals in each
 -   whole genome sequencing
 -   depth 2-9X (1000 genome project)

| CEU | Europeans in Utah (British) |
| --- | --- |
| GBR | Great Britain |
| IBS | Iberian/Spain |
| TSI | Italien |


First let's set the paths


```
cd ~/structure2

# NB this must be done every time you open a new terminal
ThePath=/course/popgen23/anders/popstructureII


## Copy positions and sample information
cp $ThePath/eu1000g.sample.Info .

ls $ThePath/eu1000g.small.beagle.gz
```

## Explore the input data.

Take a quick look at the sample data.

First try to get an overview of the dataset by looking at the information file and making a summary using the following code:

```
# View first lines of sample info file
echo First lines in sample info file
head eu1000g.sample.Info

echo Count the number of samples from each population
cut -f 2 -d " " eu1000g.sample.Info | sed 1d| sort | uniq -c
```
 -   How many samples from each country?

Now let's have a look at the genotype likelihood (GL) file that you have created with ANGSD. It is a "beagle format" file called all.beagle.gz - and will be the input file to PCAangsd. The first line in this file is a header line and after that it contains a line for each locus with GLs. By using the unix command wc we can count the number of lines in the file:

```
gunzip -c $ThePath/eu1000g.small.beagle.gz | wc -l
```

 -   Use this to find out how many loci there are GLs for in the data set?

Next, to get an idea of what the GL file contains try (from the command line) to print the first 9 columns of the first 7 lines of the file:

```
zcat $ThePath/eu1000g.small.beagle.gz | head -n 7 | cut -f1-9 | column -t
```


## Ignore the "Broken pipe"

In general, the first three columns of a beagle file contain marker name and the two alleles, allele1 and allele2, present in the locus (in beagle A=0, C=1, G=2, T=3).

All following columns contain genotype likelihoods (three columns for each individual: first GL for homozygote for allele1, then GL for heterozygote and then GL for homozygote for allele2). Note that the GL values sum to one per site for each individuals. This is just a normalization of the genotype likelihoods in order to avoid underflow problems in the beagle software it does not mean that they are genotype probabilities.

 -   Based on this, what is the most likely genotype of Ind0 in the first locus and the locus six?


Run PCangsd with to estimate the covariance matrix while jointly estimating the individuals allele frequencies

```
pcangsd -b $ThePath/eu1000g.small.beagle.gz -o EUsmall -t 5
```

This takes around 2 min to run. The program estimates the covariance matrix that can then be used for PCA. Look at the output from the program.
    The algorithm might only need a low number of PCs to estimate the allele freuqencies. How many significant PCs (see MAP test in output)?

Now plot the results in R:

```
 ## R
 cov <- as.matrix(read.table("~/structure2/EUsmall.cov"))
 e<-eigen(cov)
 ID<-read.table("~/structure2/eu1000g.sample.Info",head=T,stringsAsFactors=T)
 plot(e$vectors[,1:2],col=ID$POP,xlab="PC1",ylab="PC2")
 legend("topleft",fill=1:4,levels(ID$POP))
 ## close R 
```

 -    Does the plot look like you expected? Which populations are close and distant to each other?

Since the European individuals in 1000G are not simple homogeneous disjoint populations it is hard to use PBS/FST or similar statistics to infer selection based on populating differences ( you will learn about these later). However, PCA offers a good description of the differences between individuals without having the define disjoint groups.

Let's try to infer selection along the genome based on the PCA
```
pcangsd -b $ThePath/eu1000g.small.beagle.gz -o EUsmall --selection --sites_save --minMaf 0 -t 5
```

The analysis takes aboud two minutes. We also need to keep track of whether a SNP is used in the analysis or not, which can be done based on the output. Create a file with the SNP location info that you will need to plot the results (the third column indicate if the site is used=1 or not=0)


```
# Create file with position and chromosome
paste <(zcat $ThePath/eu1000g.small.beagle.gz| cut -f 1 | sed 's/\_/\t/g' | sed 1d ) EUsmall.sites  > EUsmall.sites.info

head  EUsmall.sites.info 
```


Next, plot the results of the selection scan

```
# open R
library(RcppCNPy)

## function for QQplot
qqchi<-function(x,...){
lambda<-round(median(x)/qchisq(0.5,1),2)
 qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed -log10(pvalues)",xlab="Expected -log10(pvalues)",...);abline(0,1,col=2,lwd=2)
legend("topleft",paste("lambda=",lambda))
}

s <- npyLoad("EUsmall.selection.npy")
# convert test statistic to p-value
pval<-1-pchisq(s,1)
## make QQ plot to QC the test statistics
qqchi(s)
```


The above is a QQ plot of the p-values from the selection scan. If the test statistics is good them most point will follow the red line which only a few (<1%) will deviate.

 -    Did the test perform well?

Finally, let's plot the results of the scan along the genome:
```
## open R
## read positions (hg38)
p<-read.delim("~/structure2/EUsmall.sites.info",colC=c("factor","integer","integer"),head=F)
names(p)<-c("chr","pos","keep")
## make manhatten plot
plot(-log10(pval),col=p$chr[p$keep==1],xlab="Chromosomes",main="Manhattan plot")
## dont close R
```

Lest zoom in
```
# continue in R
## zoom into region
w<-range(which(pval<1e-7)) + c(-100,100)

keep<-w[1]:w[2]
plot(p$pos[keep],-log10(pval[keep]),col=p$chr[keep],xlab="HG38 Position chr2")
```

```
## see the position of the most significant SNP
cat("Base pair position of the highest peak ")
p$pos[which.max(s)]
```
See if you can make sense of the top hit:
 -   Go to the UCSC browser (https://genome.ucsc.edu/cgi-bin/hgGateway)
 -   Choose human GRCh38/hg38
 -   Search for the position of the top hit and identify the genes at that locus
