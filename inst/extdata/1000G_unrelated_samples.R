#Script that reads in the pedigree information from 1000G and writes out a file 
#containing only the unrelated European samples and a file containing only the 
#unrelated African samples (should have 503 and 661 samples respectively)

#Run from /net/dumbo/home/balr/UKBB/R_package_vignette/vcfs

library("data.table")
library("stringr")

df_unrel <- fread(cmd = str_glue("zgrep -v ^## 1000G_2504_high_coverage.sequence.index"),
                  header = TRUE, data.table = FALSE)

df_ped <- fread("1kGP.3202_samples.pedigree_info.txt", 
                header = TRUE, data.table = FALSE)

sum(df_unrel$SAMPLE_NAME %in% df_ped$sampleID) #2504 - hence this is the right column to use

#All European populations
pop_euro <- c("CEU", "FIN", "GBR", "IBS", "TSI")

#All African populations
pop_afr <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI")

df_unrel_euro <- df_unrel[df_unrel$POPULATION %in% pop_euro,]
df_unrel_afr <- df_unrel[df_unrel$POPULATION %in% pop_afr,]

nrow(df_unrel_euro) #503
nrow(df_unrel_afr) #661

#Hence we have the right number of samples 
euro_samps <- df_unrel_euro$SAMPLE_NAME
afr_samps <- df_unrel_afr$SAMPLE_NAME

#Now write out the two sets of samples
fileConn <- file("european_unrelateds.txt")
writeLines(euro_samps, fileConn)
close(fileConn)

fileConn <- file("african_unrelateds.txt")
writeLines(afr_samps, fileConn)
close(fileConn)











