#' @title Perform two-stage approach for rare-variant aggregation testing
#'
#' @description A function that performs the two-stage approach for rare-variant
#' aggregation testing.
#'
#' @param mask_list A list where each element is a 1 column data frame of variants
#' to include in each test.
#'
#' @param allele_freq_test Data frame containing the columns "CHROM", "POS",
#' "REF", "ALT", "AF", "AC", "PVAL"
#'
#' @param two_stage_threshold -log10(p-value) threshold to use for two-stage
#' approach. Default = 3.
#'
#' @param burden Perform simple (unweighted) burden test? Must be TRUE or FALSE.
#'
#' @param wburden Perform weight burden test? Weights are from beta(MAF, 1, 25)
#' distribution. Must be TRUE or FALSE.
#'
#' @param SKAT Perform simple (unweighted) SKAT? Must be TRUE or FALSE.
#'
#' @param residual_variance The residual variance from the study
#'
#' @param sample_size The sample size from the study
#'
#' @importFrom data.table fread
#' @importFrom stringr str_glue
#' @importFrom stats pnorm
#'
#' @return TRUE/FALSE to compute covariance or not if at least one of the masks
#' passes the -log10(p-value) threshold

two_stage_test <- function(mask_list, allele_freq_test, two_stage_threshold = 3,
                           burden, wburden, SKAT, residual_variance, sample_size){

  #Create column putting the SNPs in RAREMETAL format
  allele_freq_test$SNP <- paste(allele_freq_test$CHROM, allele_freq_test$POS,
                                allele_freq_test$REF, allele_freq_test$ALT, sep = ":")

  #Loop over group files to perform two-stage approach
  to_return <- FALSE
  i <- 1
  num_masks <- length(mask_list)
  while(!to_return & i <= num_masks){
    allele_freq_test_temp <- allele_freq_test[allele_freq_test$SNP %in% mask_list[[i]][,1],]

    if(burden){
      test_stat_num <- sum(allele_freq_test_temp$U_STAT)
      test_stat_var <- sum(2*sample_size*allele_freq_test_temp$AF*(1-allele_freq_test_temp$AF))/residual_variance ## Should change this?
      test_stat <- test_stat_num/sqrt(test_stat_var)
      pval <-  2*stats::pnorm(-abs(test_stat)) #Looks like I get the right answer comparing to RAREMETAL output but do a more comprehensive look later

      if(-log10(pval) > two_stage_threshold){
        to_return <- TRUE
        break
      }
    }

    if(wburden){
      weights_curr <- dbeta(allele_freq_test_temp$AF, 1, 25)
      test_stat_num_wb <- sum(weights_curr * allele_freq_test_temp$U_STAT)
      test_stat_var_wb <- sum(weights_curr^2 * (2*sample_size*allele_freq_test_temp$AF*(1-allele_freq_test_temp$AF)))/residual_variance
      test_stat_wb <- test_stat_num_wb/sqrt(test_stat_var_wb)
      pval_wb <-  2*stats::pnorm(-abs(test_stat_wb)) #Not getting exact same but similar - need to check more comprehensively

      if(-log10(pval_wb) > two_stage_threshold){
        to_return <- TRUE
        break
      }
    }

    i <- i + 1
  }

  return(to_return)

}

