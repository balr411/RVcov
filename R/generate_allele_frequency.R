#' @title Generate a data frame with allele frequencies from VCF
#'
#' @description A function that takes a VCF file and returns a data frame containing
#' columns CHROM, POS, REF, ALT, AN, AC, AF
#'
#' @param vcf_file String containing the file path and name containing the VCF
#'  of the genotypes of the  set of individuals to be used as the reference
#'  panel.
#'
#' @param chr The chromosome number on which you want to perform the analysis.
#' Must be in the same format as in the CHROM column of the VCF. This could be
#' for example c{chr_number}, chr{chr_number}, or {chr_number}. ie. c1, chr1, or 1
#' for chromosome 1.
#'
#' @param gene_start The starting base pair position of the gene on the chromosome.
#' Default = NULL. If you have large reference panel VCFs it is recommended to
#' supply gene_start and gene_end.
#'
#' @param gene_end The ending base pair position of the gene on the chromosome.
#' Default = NULL. If you have large reference panel VCFs it is recommended to
#' supply gene_start and gene_end.
#'
#' @importFrom stringr str_glue
#'
#' @return A data frame with columns CHROM, POS, REF, ALT, AN, AC, AF

generate_allele_frequency <- function(vcf_file, chr, gene_start = NULL, gene_end = NULL){

  if(!is.null(gene_start) & !is.null(gene_end)){
    range_curr <- paste0(chr, ":", as.numeric(gene_start), "-", as.numeric(gene_end))
    cm <- stringr::str_glue("tabix -h {vcf_file} {range_curr} | bcftools +fill-AN-AC | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AN\\t%AC\\n'")
  }else{
    cm <- stringr::str_glue("bcftools +fill-AN-AC {vcf_file} | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%AN\\t%AC\\n'")
  }

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
