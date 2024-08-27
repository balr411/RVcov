#' @title Generate a data frame with allele frequencies from VCF
#'
#' @description A function that takes a VCF file and returns a data frame containing
#' columns CHROM, POS, REF, ALT, AN, AC, AF
#'
#' @param vcf_file String containing the file path and name containing the VCF
#'  of the genotypes of the  set of individuals to be used as the reference
#'  panel.
#'
#' @importFrom stringr str_glue
#'
#' @return A data frame with columns CHROM, POS, REF, ALT, AN, AC, AF

generate_allele_frequency <- function(vcf_file){
  #Note have to add functionality for gene_start and gene_end later
  cm <- stringr::str_glue("bcftools +fill-AN-AC {vcf_file} | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AN\\t%AC\\n'") #Need to add functionality for gene_start and gene_end
  af <- system(cm, intern = TRUE)

  #Create data frame from the vector
  af_list <- strsplit(af, split = "\t")
  df <- t(as.data.frame(af_list))
  rownames(df) <- NULL
  df <- as.data.frame(df)
  colnames(df) <- c("CHROM", "POS", "REF", "ALT", "AN", "AC")
  df$POS <- as.numeric(df$POS)
  df$AN <- as.numeric(df$AN)
  df$AC <- as.numeric(df$AC)
  df$AF <- df$AC/df$AN

  #Following should work regardless if in form chr{num} or c{num}
  df$CHROM <- gsub("c", "", df$CHROM)
  df$CHROM <- as.numeric(gsub("hr", "", df$CHROM))

  return(df)
}
