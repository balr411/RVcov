#R script that reads in the annotation file containing the rare variants from 
#1000 Genomes chr 2 and determines which gene we will use for the vignette

library("data.table")
library("stringr")

df_rare <- fread("20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr2.recalibrated_variants.annotated.coding_rare.txt",
                 header = TRUE, data.table = FALSE)

nrow(df_rare) #68096

#Keep only SNPs with the PASS filter
df_rare <- df_rare[df_rare$FILTER == "PASS",]

nrow(df_rare) #65001

#Keep only the first occurence of each SNP
df_rare_nonDup <- df_rare[match(unique(df_rare$ID), df_rare$ID),]

nrow(df_rare_nonDup) #46439

#Find the genes that appear the most
sort(table(df_rare_nonDup$GENE), decreasing = TRUE)[1:3]

#TTN   NEB FSIP2
#Unsurprisingly TTN is the largest gene 
#Pick NEB gene instead as TTN is a very large gene (NEB has 503 variants)
#According to NCBI website (https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=4703),
#on GRCh38 the gene is between the positions 151485339 and 151734476
#See what the first and last positions are here

#First check if the data frame is sorted by position 
is.unsorted(df_rare_nonDup[df_rare_nonDup$GENE == "NEB", 2]) #FALSE

df_NEB <- df_rare_nonDup[df_rare_nonDup$GENE == "NEB",]
df_NEB$POS[1] #151485773
df_NEB$POS[nrow(df_NEB)] #151729622

#These variants both fall within the NCBI website gene boundaries so use those

#Now count the number of deleterious predictions/LOF to make sure we have enough to make masks
#For loss of function, look at the effect column and LOF column 
#For our paper we used "stop_gained", "splice_acceptor_variant", "stop_lost", "splice_donor_variant", "frameshift_variant", "start_lost" as LOF

df_NEB$LOF_num <- 0
df_NEB$LOF_num[!is.na(df_NEB$LOF)] <- 1 #Assume those with non-NA LOF entry are LOF variants (seems to be true for the one I searched - rs549794342)
sum(df_NEB$LOF_num == 1) #8

df_NEB$LOF[df_NEB$EFFECT == "frameshift_variant"] #2 variants but only the second has non-NA LOF
df_NEB$LOF[df_NEB$EFFECT == "stop_gained"] #All 3 have non-NA LOF
df_NEB$LOF[df_NEB$EFFECT == "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant"] #Only one and it has non-NA LOF
df_NEB$LOF[df_NEB$EFFECT == "splice_donor_variant&intron_variant"] #All 3 have non-NA LOF
#This gives all 8 of the LOF variants, with the only effect normally on a LOF list not having non-NA LOF (??). So I will keep this list of 8 as our pLOF variants

df_NEB$MutationTaster_pred_num = str_count(df_NEB$MutationTaster_pred, "D") + str_count(df_NEB$MutationTaster_pred, "A")
df_NEB$Polyphen2_HDIV_pred_num = str_count(df_NEB$Polyphen2_HDIV_pred, "D") + str_count(df_NEB$Polyphen2_HDIV_pred, "P")
df_NEB$Polyphen2_HVAR_pred_num = str_count(df_NEB$Polyphen2_HVAR_pred, "D") + str_count(df_NEB$Polyphen2_HVAR_pred, "P")
df_NEB$SIFT_pred_num = str_count(df_NEB$SIFT_pred, "D")

df_NEB$MutationTaster_pred_num = ifelse(df_NEB$MutationTaster_pred_num > 0, 1, 0)
df_NEB$Polyphen2_HDIV_pred_num = ifelse(df_NEB$Polyphen2_HDIV_pred_num > 0, 1, 0)
df_NEB$Polyphen2_HVAR_pred_num = ifelse(df_NEB$Polyphen2_HVAR_pred_num > 0, 1, 0)
df_NEB$SIFT_pred_num = ifelse(df_NEB$SIFT_pred_num > 0, 1, 0)
df_NEB$Sum_Algs = rowSums(df_NEB[,c("MutationTaster_pred_num", "Polyphen2_HDIV_pred_num",
                                    "Polyphen2_HVAR_pred_num", "SIFT_pred_num")])
df_NEB$Sum_Algs = ifelse(df_NEB$Sum_Algs == 4, 1, 0)

table(df_NEB$Sum_Algs) #260 0s; 184 1s

#Now write this new annotation file 
fwrite(df_NEB, "NEB/chr2.1000G.NEB.rare_variants.annotation.txt", col.names = TRUE, 
       row.names = FALSE, sep = "\t")

