# Introduction to basic NGS data, Data processing, and formats

# Exercise 0

Log in to the server


Notes before we start:
- In coding blocks, comments are preceeded with a # and all of the code can be pasted in to the terminal
    ```bash
  # code you should run in the terminal
    ls
    ```

- **Please read the comments before running each line of code. So we know what we're doing or what we may consider.**

## Linux paths for today
Lets first set some paths. Paste the following code in the terminal

```bash
# bash command
introNGS=/course/popgen23/rute/IntroNGS
ref=$introNGS/referenceGenome/chr21.fa
realignedBam=$introNGS/2-filterBam/NA19238.q30.realign.intervals.bam
```

## Tools for today
lets see if the software is installed. Again paste the code in the terminal and look for errors. 
```bash
# bash command
which TrimmomaticPE
which fastqc
which bwa
which samtools
```

<details>
<summary>cilck here to see the output</summary>

```text
/usr/bin/TrimmomaticPE
/usr/bin/fastqc
/usr/bin/bwa
/usr/bin/samtools
```

</details>


Go to the folder called MATERIAL and view the files

```bash
cd ~/MATERIAL
ls
```



# EXERCISE 1 

 assess the quality of the fastq files with fastqc:
````bash
fastqc read1.fastq.gz
````
#Output:
#read1_fastqc.html
#read1_fastqc.zip

 - download the html to your computer and open it in a browser ( e.g. double click on it )


If you want to access the tables used for the graphics:
```bash
unzip read1_fastqc.zip
cd read1_fastqc
```
use your favorite text editor to visualize this file
```bash
gedit fastqc_data.txt
```
go back to the MATERIAL folder
```bash
cd ~/MATERIAL
```

# EXERCISE 2 
Exercise 2: Check out the Trimmomatic options:
```bash

TrimmomaticPE 
```
print the contents of the “0-rm-adapt-clean-trimmomatic” file:
```bash

cat 0-rm-adapt-clean-trimmomatic
```

note: to make a text file executable run
```bash
chmod 755 0-rm-adapt-clean-trimmomatic
```

Run the script "0-rm-adapt-clean-trimmomatic" as:
```bash
./0-rm-adapt-clean-trimmomatic read1.fastq.gz read2.fastq.gz  
```

Join the unpaired reads into a single file:

```bash
cat read1_u.fq.gz read2_u.fq.gz > read_u.fq.gz
```

And remove the redundant files:
```bash
rm read1_u.fq.gz read2_u.fq.gz 
```

# Bonus EXERCISE 3 ( optional if you have finished 1 and 2)

## fastqc stats



run fastqc on the 3 output files
```bash
fastqc read1_p.fq.gz
fastqc read2_p.fq.gz #optional 
fastqc read_u.fq.gz #optional
```
Open the html file in the browser. 

 - What is the total number of reads?
 -  - What is the lowest base quality score?
 - After which position in the read does the quality drop significantly?
 - Do the sequences contain adapter sequence?

Compare to the results from before the trimming.




# Exercise 4: mapping with bwa

#Step 1: Index the genome (DON’T do it now, it takes time).
#We will now be working with human chromosome 21 and this will be our “reference genome”.

this is the command line used for short genomes
```bash
#bwa index -a is reference.fasta
## instead of running the command we will use the already index reference
ref=$introNGS/referenceGenome/chr21.fa
```

### Step 2: map reads to chr21

first set the read group, very useful if multiple runs per sample
```bash
readGroup='@RG\tID:ID\tCN:test\tDS:description\tPL:platform\tSM:sample'
```
map the reads to the reference genome using (-t) 5 CPU threads 
```bash
bwa mem -t 5 -M -R $readGroup $ref $introNGS/NA19238.fastq.gz > NA19238.sam
```

Convert sam to bam 
```bash
samtools view -bS NA19238.sam > NA19238.bam
```
Now you should have a file named “NA19238.bam”. 

Check its size relative to the sam file.

A bam file is a binary file, we can handle it using samtools. 

Try to see the difference between the following commands:
```bash
samtools view NA19238.bam | head
samtools view -H  NA19238.bam | head
samtools view -h  NA19238.bam | head
samtools view -h  NA19238.bam > NA19238.sam
head NA19238.sam
```

 - How long are the reads? (hint: use the cigar string info or run fastqc)

```bash
head NA19238.sam
```
How many reads are mapped with quality > 30? 

```bash
samtools view -q 30 NA19238.bam | wc -l 
```

All the options with

```bash
samtools view
```

# EXERCISE 5

### Exercise 5: visualizing and filtering sam/bam files

First sort the bamfile by coordinate:
```bash

samtools sort  NA19238.bam -o NA19238.sorted.bam
```

Then index the bam file:
```bash

samtools index NA19238.sorted.bam
```

Visualize the bam file using “samtools tview”
```bash

samtools tview NA19238.sorted.bam
```

Note: chr21 does not have reference bases for the first 9719760 sites.
Type "g" to enter the position of the alignment that you want to see.
 - Type "chr21:10028350". Why do you think so many reads mapped to this position?
 - Use left/right arrows to navigate the alignment or the space key
 - Type "?" for other options
 - Browse the alignment using the arrow keys or space.
 - Color the bases by the quality scores by typing 'b' 
 -  (you can go back by typing 'm'). Do you see any patterns?
Ctrl+C to quit tview
if you want to vizualize the reference genome sequence as well:

```bash

samtools tview NA19238.sorted.bam $ref
```



save mapped reads with mapping quality >30
```bash

samtools view -b -q 30 NA19238.sorted.bam | samtools sort -o NA19238.q30.bam
```

note: paralogs are reads that map equally well to different regions of the assembly;  bwa mem assigns low quality scores to paralogs

Check the bam file again. First need to index the new bam file:
```bash
samtools index NA19238.q30.bam
```

And then open tview
```bash

samtools tview NA19238.q30.bam $ref
```

 - Type “g” and then write “chr21:10028350”. Do you still see the large stack of reads?


# EXERCISE 6 

```bash
samtools tview -p chr21:10029501 NA19238.q30.bam $ref
```

Check a region before and after realignment 
<details>
<summmary >e.g. read name SRR002987.11216777 </summary>
bam line before realignment CIGAR: 36M
#SRR002987.11216777      0       chr21   21515742        60      36M     *       0       0       TCTAGCTAGAAGACCAAACCAGTGGTAATATTCCCT      >>=>>>>>>>>?>>?F=?>A2>7?>?<<@;<<;;9<    NM:i:1  MD:Z:35C0       AS:i:35 XS:i:0  RG:Z:ID

bam line after realignment CIGAR: 32M1D4M + OC:Z:36M tag

#SRR002987.11216777      0       chr21   21515742        70      32M1D4M *       0       0       TCTAGCTAGAAGACCAAACCAGTGGTAATATTCCCT      >>=>>>>>>>>?>>?F=?>A2>7?>?<<@;<<;;9<    OC:Z:36M  
</details>
    
```bash
realignedBam=$introNGS/2-filterBam/NA19238.q30.realign.intervals.bam
samtools tview -p chr21:21515742 NA19238.q30.bam $ref
```

<details>
<summary> output </summary>
21515771  
attcccctta
......Y...
......T    
......T    
......T.    
...*....... 
...*........
</details>

```bash
samtools tview -p chr21:21515742 $realignedBam $ref
```
<details>
<summary> output </summary>

21515771  
attcccctta
...K......
...*....  
...*....  
...*..... 
...*......
...*......
</details>
M: Match (alignment column containing two letters, either different letters (mismatch) or identical.
a * denotes a deletion/missing match; more details on the pileup format at

Visualize the region containing the variant with mpileup (requires bed file)
https://en.wikipedia.org/wiki/Pileup_format

```bash
samtools mpileup -l $introNGS/variantRegion.bed --reference $ref NA19238.q30.bam
```
get depth distribution for the first 1Mb
```bash
samtools mpileup NA19238.q30.bam| head -n 1000000 |cut -f4 | sort -n | uniq -c > dep1
```

Open R  to plot depth distribution:

```R
#In R:

# Visualize the depth distribution
d<-read.table("~/MATERIAL/dep1")
barplot(d[2:20,1],names=1:19,xlab="depth") #ignores invariable sites

```

