#' @title Generate group/mask file for rare-variant aggregation testing
#'
#' @description A function that generates and writes to disk in RAREMETAL format
#'  a group/mask file for rare variant aggregation testing.
#'
#' @param anno A VEP annotation file with at minimum columns named #Uploaded_variation
#'  containing the variants in format chr:pos:ref:alt, Consequence containing the predicted
#'  consequences of the variants, SYMBOL containing the gene names that the variants
#'  belong to, and PICK which is 1 if that row is to be used and "-" otherwise
#'
#' @param allele_freq The in-sample allele frequency file with columns CHROM,
#'  POS, REF, ALT, AF, AC, and PVAL.
#'
#' @param gene The name of the gene from the annotation file to use.
#'
#' @return Nothing
#'

generate_group_file <- function(anno, allele_freq, gene){

}
