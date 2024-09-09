#Script that generates phenotypes and a .ped file for the Europeans from 
#1000G

#Run from /net/dumbo/home/balr/UKBB/R_package_vignette/vcfs

library("data.table")
library("stringr")
library("VariantAnnotation")
library("dplyr")
library("Matrix")

## Assume you have the following variables:
## maf_full: MAF of all v variants in gene
## annos: v*3 dataframe with 3 columns for indicators of variants belonging to the 3 masks:
##        column1 - value of 1 means variant belongs to PTV mask
##        column2 - value of 1 means variant belongs to PTV+Missense(All) mask
##        column3 - value of 1 means variant belongs to PTV+Missense(Deleterious) mask
## p_caus: vector of causal probabilities for the three annotation categories: e.g. c(0.8, 0.3, 0.01)
## n: sample size
## h: heritability of gene
## geno: n*v genotype matrix

#First read in the allele frequency file for Europeans
allele_freq_euro <- fread(cmd = str_glue("zcat NEB/chr2.1000G.NEB.european_unrel.allele.freq.gz"), 
                          header = FALSE, data.table = FALSE)

names(allele_freq_euro) <- c("CHROM", "POS", "REF", "ALT", "AN", "AC")

allele_freq_euro$AC <- as.numeric(allele_freq_euro$AC)
allele_freq_euro$AN <- as.numeric(allele_freq_euro$AN)

allele_freq_euro$AF <- allele_freq_euro$AC/allele_freq_euro$AN

#Get rid of multi-allelics and rows with NAs
allele_freq_euro <- anti_join(allele_freq_euro, allele_freq_euro[duplicated(allele_freq_euro[,2]),], by = "POS")
allele_freq_euro <- allele_freq_euro[complete.cases(allele_freq_euro),]

#Now reduce to 0 < MAF < 0.01
sum(allele_freq_euro$AF > 0.5) #151
sum(allele_freq_euro$AF > 0.99) #12
sum(allele_freq_euro$AF == 1) #12 - hence all AFs greater than 0.99 are equal to 1 so for our purposes AF = MAF

allele_freq_euro <- allele_freq_euro[allele_freq_euro$AF > 0 & allele_freq_euro$AF < 0.01,]

nrow(allele_freq_euro) #1666

#Now read in the annotation file
anno <- fread("NEB/chr2.1000G.NEB.rare_variants.annotation.txt",
              header = TRUE, data.table = FALSE)

#Match SNPs between the two data frames
anno$SNP <- paste(anno$CHROM, anno$POS, anno$REF, anno$ALT, sep = ":")
allele_freq_euro$SNP <- paste(allele_freq_euro$CHROM, allele_freq_euro$POS, allele_freq_euro$REF, allele_freq_euro$ALT, sep = ":")

snps_intersect <- intersect(anno$SNP, allele_freq_euro$SNP)
length(snps_intersect) #97 - hence 97 causal variants to work with 

anno <- anno[anno$SNP %in% snps_intersect,]
allele_freq_euro <- allele_freq_euro[allele_freq_euro$SNP %in% snps_intersect,]

#Out of curiosity check to see if any of these variants have non-zero MAF in the African VCF
allele_freq_afr <- fread(cmd = str_glue("zcat NEB/chr2.1000G.NEB.african_unrel.allele.freq.gz"), 
                          header = FALSE, data.table = FALSE)

names(allele_freq_afr) <- c("CHROM", "POS", "REF", "ALT", "AN", "AC")

allele_freq_afr$AC <- as.numeric(allele_freq_afr$AC)
allele_freq_afr$AN <- as.numeric(allele_freq_afr$AN)

allele_freq_afr$AF <- allele_freq_afr$AC/allele_freq_afr$AN

#Get rid of multi-allelics and rows with NAs
allele_freq_afr <- anti_join(allele_freq_afr, allele_freq_afr[duplicated(allele_freq_afr[,2]),], by = "POS")
allele_freq_afr <- allele_freq_afr[complete.cases(allele_freq_afr),]

allele_freq_afr$SNP <- paste(allele_freq_afr$CHROM, allele_freq_afr$POS, allele_freq_afr$REF, allele_freq_afr$ALT, sep = ":")

allele_freq_afr <- allele_freq_afr[allele_freq_afr$AF > 0 & allele_freq_afr$AF < 0.01,] #Note there is one variant with 0.99<AF<1

sum(allele_freq_afr$SNP %in% snps_intersect) #7 - hence there are some!

#Now read in the genotype matrix for the European samples
current_range <- ScanVcfParam(which = GRanges(seqnames = Rle("chr2"), ranges = IRanges(start = 151485339, end = 151734476)))
gt_current <- readGT(file = "NEB/chr2.1000G.european_unrel.NEB.vcf.gz", param = current_range)

#Note that the SNP names in gt_current are in the format chr:pos_ref/alt so change
#this in the anno and allele freq data frames in order to match 
anno$SNP_GT_FORM <- paste(paste(paste(anno$CHROM, anno$POS, sep = ":"), anno$REF, sep = "_"), anno$ALT, sep = "/")
allele_freq_euro$SNP_GT_FORM <- paste(paste(paste(allele_freq_euro$CHROM, allele_freq_euro$POS, sep = ":"), allele_freq_euro$REF, sep = "_"), allele_freq_euro$ALT, sep = "/")

gt_current <- gt_current[rownames(gt_current) %in% anno$SNP_GT_FORM,]

#Make sure they are in the same order
all.equal(rownames(gt_current), anno$SNP_GT_FORM) #TRUE
all.equal(rownames(gt_current), allele_freq_euro$SNP_GT_FORM) #TRUE

#Now change gt_current to contain only 0, 1, or 2
gt_current <- t(gt_current)
gt_current[gt_current %in% c("0/0", "0|0")] <- 0
gt_current[gt_current %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
gt_current[gt_current %in% c("1/1", "1|1")] <- 2

unique_vals <- c()
for(i in 1:ncol(gt_current)){
  unique_vals <- c(unique_vals, unique(gt_current[,i]))
}

n_row <- nrow(gt_current)
n_col <- ncol(gt_current)

geno <- suppressWarnings(Matrix::Matrix(as.numeric(gt_current), nrow = n_row, ncol = n_col, sparse = TRUE)) 
geno <- sweep(geno, 2, 2*allele_freq_euro$AF) #subtract the allele frequencies rowwise

#Now impute the missing values
is_na <- which(is.na(geno), arr.ind = TRUE)
geno[is_na] <- colMeans(geno, na.rm = TRUE)[is_na[,"col"]]

#Put in format of below code
maf_rare <- allele_freq_euro$AF
pLOF_missenseNarrow <- rep(0, nrow(anno))
pLOF_missenseNarrow[anno$LOF_num == 1] <- 1
pLOF_missenseNarrow[anno$Sum_Algs == 1] <- 1

annos <- data.frame(pLOF = anno$LOF_num,
                    pLOF_missenseBroad = rep(1, nrow(anno)),
                    pLOF_missenseNarrow = pLOF_missenseNarrow)

#Assign probability of causal to 0.8 for PTVs, 0.3 for deleterious missense, and 0.01 for non-deleterious missense
p_caus <- c(0.8, 0.3, 0.01)

#Set n and desired heritability - note that we want nh^2 approximately equal to 50 (note h = h^2 here)
n <- 503
h <- 50/503

causal1 = which(annos[,1] == 1) # find indices of PTVs
causal2 = which(annos[,3] == 1 & annos[,1] == 0) # find indices of deleterious missense variants
causal3 = which(annos[,2] == 1 & annos[,3] == 0 & annos[,1] == 0) # find indices of other missense variants

#Choose causal variants 
set.seed(1)
unif_rand = runif(length(causal1)) #Simulate random uniform numbers for PTVs
causal = causal1[which(unif_rand <= p_caus[1])] #Set PTVs as causal if corresponding uniform random number <= p_caus[1]
unif_rand = runif(length(causal2)) #Simulate random uniform numbers for deleterious missense variants
causal = c(causal, causal2[which(unif_rand <= p_caus[2])]) #Add causal del mis variants to causal set
unif_rand = runif(length(causal3)) #Simulate random uniform numbers for other missense variants
causal = c(causal, causal3[which(unif_rand <= p_caus[3])]) #Add causal other mis variants to causal set
causal = sort(causal) #Sort indices of causal variants in causal set

#Calculate effect sizes of causal variants
beta_true = rep(0, dim(annos)[1]) #Set effect sizes of all variants to 0
const = sqrt(0.5*h/sum(((log10(maf_rare[causal]))^2)*maf_rare[causal]*(1-maf_rare[causal]), na.rm = T)) #Calculate constant of proportionality
beta_true[causal] = -const * log10(maf_rare[causal]) #Set effect sizes of causal variants proportional to -log10(MAF)

#Simulate phenotype
set.seed(1)
y = geno %*% beta_true + rnorm(n, sd = sqrt(1-h))
y <- as.vector(y)

#Now create .ped files
#Read in phenotype file to get sex of individuals
df_full_ped <- fread("1kGP.3202_samples.pedigree_info.txt", 
                header = TRUE, data.table = FALSE)

df_full_ped <- df_full_ped[df_full_ped$sampleID %in% rownames(gt_current),]

#Chekc that they are in the same order
all.equal(rownames(gt_current), df_full_ped$sampleID) #Not true so rearrange
idx <- match(rownames(gt_current), df_full_ped$sampleID)
df_full_ped <- df_full_ped[idx,]
all.equal(rownames(gt_current), df_full_ped$sampleID) #TRUE now

df_ped <- data.frame(FID = rownames(gt_current),
                     IID = rownames(gt_current),
                     PATID  = rep(0, nrow(gt_current)),
                     MATID = rep(0, nrow(gt_current)),
                     Sex = df_full_ped$sex,
                     Phenotype = y)

#Now write out the .ped file 
fwrite(df_ped, 
       file = "NEB/simulated_phenotype.ped",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE, 
       col.names = TRUE)




