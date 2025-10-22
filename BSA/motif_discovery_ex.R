# make exercise - 
#~/work/COURSES/BIO/BSA/pHMM/globins4.fasta
#§less globins.msa

#-------------------------------------------------------------------------------------------------

# AA slides
	#https://docs.google.com/presentation/d/1TR6eZsXPHXutIcFJQUwDjS-qjH1e_-F9_QtiIXNHzUE/edit?slide=id.g39b021e8e8b_0_758#slide=id.g39b021e8e8b_0_758

# calculate the enrichment of the sequence TGTA in ~/work/COURSES/BIO/BSA/motif_discovery/PUM2.top500.fa assuming 0,1,2 order markov models
# ignore edge effect
# assume words are independnet 

# input: ~/work/COURSES/BIO/BSA/motif_discovery/PUM2.top500.fa


#1) calculate background
#2) calculate lambda for each markov order
#3) if enrichment is significantly different from background

#-------------------------------------------------------------------------------------------------
library(stringr)
#-------------------------------------------------------------------------------------------------
input<-read.table('~/work/COURSES/BIO/BSA/motif_discovery/PUM2.top500.fa')

# head(input)
#                               V1
#1                            >1_0
#2 ATACCATTTTATTTAGAATATATTATACTTC
#3                            >1_1
#4 AAATCATTTCCAAAGCATGTGGTGTGGCCTG

torm <- seq(1, nrow(input), 2) # remove the > headers
sequences<-input[ -torm ,]

 head(sequences)
#[1] "ATACCATTTTATTTAGAATATATTATACTTC" "AAATCATTTCCAAAGCATGTGGTGTGGCCTG"
#[3] "ATTTCAAAAGCTGTGTCTCTATTAG"       "AAAAAAAACTTTTGCAAATTTTTTTATTTTT"
#[5] "ACATATGTAGATACATACCATG"          "TCACATATATACATATGTATCTTATAATCTT"

#-------------------------------------------------------------------------------------------------
all_seq <- paste(sequences, collapse = "")

total_length<-nchar(all_seq)
#[1] 13152

#sum(str_count(all_seq, "T"))
#sum(str_count(all_seq, "TG"))

#-------------------------------------------------------------------------------------------------

#0 order: P0 = p(T)p(G)p(T)p(A) 
	# frequencies of the bases
#lambda0 = p0 * length

t=sum(str_count(all_seq, "T"))/total_length  
g=sum(str_count(all_seq, "G"))/total_length  
a=sum(str_count(all_seq, "A"))/total_length  

p0 = t * g * t * a 
#lambda0 = p0 * length
lambda0 = p0 * total_length
#-------------------------------------------------------------------------------------------------

# 1 order
#1 order: p(T) * p(G|T) * p(T|G) * p(A|T)

gt = sum(str_count(all_seq, "TG"))/sum(str_count(all_seq, "T"))
tg = sum(str_count(all_seq, "GT"))/sum(str_count(all_seq, "G"))
at = sum(str_count(all_seq, "TA"))/sum(str_count(all_seq, "T"))

p1 = t * gt * tg * at

# p(G|T) = # of tg / # of t

#lambda1 = p1 * length
lambda1 = p1 * total_length


#-------------------------------------------------------------------------------------------------

# 2nd order: p(T)  * p(G|T) * p(T|TG) * p(A|TG)
#p(T| GT) = # tgt / # gt 

# first 2 nucleotides: TG
tgt = sum(str_count(all_seq, "TGT"))/sum(str_count(all_seq, "TG"))
atg = sum(str_count(all_seq, "ATG"))/sum(str_count(all_seq, "GT"))

p2 = t * gt * tgt * atg 
lambda2 = p2 * total_length

#-------------------------------------------------------------------------------------------------


# diff lambda per order

# calc enrichment

# enrichment = occurrence divided by lambda
enrichment0: # tgta / lambda0
enrichment1: # tgta / lambda1
enrichment2: # tgta / lambda2

tgta = sum(str_count(all_seq, "TGTA"))

enrichment0 = tgta / lambda0
enrichment1 = tgta / lambda1
enrichment2 = tgta / lambda2

#> enrichment0
#[1] 5.041297
#> enrichment1
#[1] 2.472306
#> enrichment2
#[1] 1.555352

# perform poisson test for stat significance 

# poison (number of occurences, lambda, lower.tail=F): right-tail probability
ppois(tgta,lambda0,lower.tail=F)
ppois(tgta,lambda1,lower.tail=F)
ppois(tgta,lambda2,lower.tail=F)

#[1] 1.283814e-114
#[1] 6.626693e-45
#[1] 1.040599e-13

# => stat sig different from the background
#cumulative distribution function (CDF) of the Poisson distribution.

# which markov model better - log likelihoods 