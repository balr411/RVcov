#' @title Read in sample size and residual variance from RAREMETALWORKER score
#' stat file.
#'
#' @description Read in sample size and residual variance from RAREMETALWORKER score
#' stat file.
#'
#' @param score_stat_file String containing the file path and name containing
#'  the RAREMETALWORKER formatted score statistic file.
#'
#'  @importFrom data.table fread
#'
#' @return A list containing the residual variance and sample size of the study

ss_residual_variance <- function(score_stat_file){
  ss <- data.table::fread(cmd = paste0("zgrep ^## ", score_stat_file), sep = NULL,
              data.table = FALSE, header = FALSE)[[1]]

  sampleSize <- strsplit(ss[4], split = "=")[[1]][2]

  stop <- FALSE
  i <- 1
  while(!stop){
    if(strsplit(ss[i], split = "\t")[[1]][1] == "##Sigma_e2_Hat"){
      residual_variance <- strsplit(ss[i], split = "\t")[[1]][2]
      stop <- TRUE
    }
    i <- i + 1
  }

  to_return <- list(sampleSize = sampleSize, residualVariance = residual_variance)
  return(to_return)
}
