#' @title Aggregation testing using RAREMETAL and external covariance
#'
#' @description A function that perfroms various rare-variant aggregation tests
#'  given a path to a RAREMETALWORKER formatted score statistic file, a path to
#'  a VCF containing the individual-level genetic data of the individuals to be
#'  used as the external reference panel, and the path to a gzipped VEP
#'  annotation file. Option to use the two-stage approach where external
#'  reference panel is used only when the -log10(p-value) using the null
#'  covariance is above some threshold.
#'
#' @param score_stat_file String containing the file path and name containing
#'  the RAREMETALWORKER formatted score statistic file.
#'
#' @param vcf_file String containing the file path and name containing the VCF
#'  of the genotypes of the  set of individuals to be used as the reference
#'  panel.
#'
#' @param anno_file String containing the file path and name containing the
#'  VEP annotation file to create the group files. If not specified, must give
#'  argument to group_file. Default = NULL.
#'
#' @param two_stage TRUE or FALSE. Do you want to perform two-stage method
#'  where external reference panel is used only when the -log10(p-value) using
#'  the null covariance is above some threshold? Default = FALSE.
#'
#' @param two_stage_threshold What -log10(p-value) threshold would you like to
#'  use for two-stage approach? Default = 3.
#'
#' @param group_file Vector of strings containing the file paths and names of
#'  the group files to perform the rare-variant aggregation tests with. If not
#'  specified, must supply anno_file and types of group files to use. Default =
#'  NULL.
#'
#' @param pLOF TRUE or FALSE. Do you want to perform aggregation tests using
#'  pLOF variants only? Default = FALSE.
#'
#' @param pLOF_narrowMissense TRUE or FALSE. Do you want to perform aggregation
#'  tests using pLOF + missense(narrow) variants only? Default = FALSE.
#'
#' @param pLOF_broadMissense TRUE or FALSE. Do you want to perform aggregation
#'  tests using pLOF + missense(broad) variants only? Default = FALSE.
#'
#'
#' @importFrom {data.table} {fread}
#' @importFrom {stringr} {str_glue}
#'
#' @return A number.
#' @examples Use examples of data that I simulate

agg_test <- function(score_stat_file, vcf_file, anno_file = NULL,
                     two_stage = FALSE, two_stage_threshold = 3,
                     group_file = NULL, pLOF = FALSE,
                     pLOF_narrowMissense = FALSE, pLOF_broadMissense = FALSE){

  #First read in the necessary columns from the score statistic file
  if(!file.exists(score_stat_file)){
    stop("Score statistic file doesn't exist!")
  }

  allele_freq_test <- fread(cmd = str_glue("zgrep -v ^# {score_stat_file}"),
                            select = c(1, 2, 3, 4, 7, 8, 17), data.table = FALSE)

  names(allele_freq_test) <- c("CHROM", "POS", "REF", "ALT", "AF", "AC", "PVAL")

  #Get group files
  #Note that if group files are given, they do not need to be read in

  if(pLOF | pLOF_narrowMissense | pLOF_broadMissense){
    if(is.null(anno_file)){
      stop("Did not give annotation file to create the desired group files.
           If you have supplied your own group files and just want to use those,
           set pLOF, pLOF_narrowMissense, and pLOF_broadMissense to FALSE")
    }

    #Read in annotation file
    anno <- fread(cmd = paste0("zgrep -v ^## ", anno_file), header = T, data.table = F)

  }





}











