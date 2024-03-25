#' @title Covariance computation using external reference panel
#'
#' @description A function that computes covariance for rare-variant aggregation
#' testing using a VCF from an external reference panel and write it
#'
#' @param gt Sparse genotype matrix containing only the variants present in the
#' test allele frequency file (ie. bi-allelic variants w/ AF > 0).
#'
#' @param allele_freq_test The in-sample allele frequency data frame with columns
#' CHROM, POS, REF, ALT, AF, AC, U_STAT, and PVAL.
#'
#' @param n The sample size of the original study
#'
#' @param residual_variance The residual variance from the original study
#'
#' @param altCovariancePath Option to give alternative path to write covariance
#' files to. Default = NULL. If left as NULL, covariance files will be written
#' in current directory.
#'
#' @param gene The name of the gene to use in the name of the covariance file.
#'  Default = 'gene'.
#'
#' @return Nothing. Writes and bgzips covariance file in the format
#' altCovariancePath/gene.estimated.cov.gz

compute_covariance <- function(gt, allele_freq_test, n, residual_variance,
                               altCovariancePath = NULL, gene = 'gene'){

  #Check that the test allele frequency data frame and genotype matrix have the
  #same number of variants - shouldn't ever happen
  if(nrow(allele_freq_test) != ncol(gt)){
    stop("Allele frequency data frame and genotype matrix have differing number of variants.")
  }

  C_mat <- crossprod(gt)
  diag(C_mat)[diag(C_mat) == 0] <- 1
  D_g <- Diagonal(x = sqrt(diag(C_mat)))
  D_mat <- Diagonal(x = sqrt(2*n*allele_freq_test$AF*(1-allele_freq_test$AF)))
  D_g_inv <- solve(D_g)
  output_mat <- (D_mat %*% D_g_inv %*% C_mat %*% D_g_inv %*% D_mat)/(n*residual_variance)

  #Now print the matrix in RAREMETAL form
  output_mat_list <- asplit(output_mat, 2)
  N <- length(output_mat_list)
  output_mat_list_final <- lapply(1:N, function(x) paste0(c(output_mat_list[[x]][x:N], ""), collapse = ","))
  output_markers_in_window_list <- lapply(1:N, function(x) paste0(c(allele_freq_test$POS[x:N], ""), collapse = ","))

  to_write <- vector(length = 3)
  to_write[1] <- "##ProgramName=RareMetalWorker"
  to_write[2] <- "##Version=4.15.1"
  to_write[3] <- "##CHROM\tCURRENT_POS\tMARKERS_IN_WINDOW\tCOV_MATRICES"

  to_write <- c(to_write, paste(allele_freq_test$CHROM, allele_freq_test$POS,
                                output_markers_in_window_list, output_mat_list_final,
                                sep = "\t"))

  if(!is.null(altCovariancePath)){
    output_cov_files <- paste0(altCovariancePath, gene, ".estimated.cov")
  }else{
    output_cov_files <- paste0(gene, ".estimated.cov")
  }

  fileConn <- file(output_cov_files)
  writeLines(to_write, fileConn)
  close(fileConn)

  system(sprintf("bgzip -f %s", output_cov_files))
  imp0_cov_gzip <- paste0(output_cov_files, ".gz")
  system(str_glue("tabix -c \"#\" -s 1 -b 2 -e 2 {imp0_cov_gzip}"))

}

