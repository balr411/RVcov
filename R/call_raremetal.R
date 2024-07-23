#' @title Burden testing using RAREMETAL
#'
#' @description A function that performs simple burden testing by making a call
#' to RAREMETAL
#'
#' @param mask_list Named list of vectors containing the variants to use in the
#' group files for burden testing.
#'
#' @param score_stat_file String containing the file path and name containing
#'  the RAREMETALWORKER formatted score statistic file.
#'
#' @param group_file_names Vector of strings containing the file paths and names of
#' the group files to perform the rare-variant aggregation tests with.
#'
#' @param altCovariancePath Option to give alternative path to write covariance
#' files to. Default = NULL. If left as NULL, covariance files will be written
#' in current directory.
#'
#' @param altGroupFilePath Optional alterntive path to write group files to.
#' Default = NULL.
#'
#' @param altRaremetalPath Optional alternative path to write RAREMETAL results
#' to. Default = NULL. If left as NULL, will write in current directory.
#'
#' @param gene The name of the gene to use in the name of the output file.
#'  Default = 'gene'.
#'
#' @param hwe Hardy-Weinberg equilibrium p-value cut-off. Default = 0.000001.
#'
#' @importFrom stringr str_glue
#'
#' @return Nothing. Will write the RAREMETAL results to
#' altRaremetalPath/gene.mask_name.burden.res - ??

call_raremetal <- function(mask_list, score_stat_file, group_file_names,
                           altCovariancePath = NULL, altGroupFilePath = NULL,
                           altRaremetalPath = NULL, gene = 'gene',
                           hwe = 0.000001){

  #First write temporary covFiles and scoreFiles to pass to RAREMETAL
  #scoreFiles first
  summary_files <- paste0(gene, ".RAREMETAL_scoreFiles")
  fileConn <- file(summary_files)
  writeLines(score_stat_file, fileConn)
  close(fileConn)

  #Now covariance files
  if(!is.null(altCovariancePath)){
    output_cov_files <- paste0(altCovariancePath, gene, ".estimated.cov.gz")
  }else{
    output_cov_files <- paste0(gene, ".estimated.cov.gz")
  }

  cov_files <- paste0(gene, ".RAREMETAL_covFiles")
  fileConn <- file(cov_files)
  writeLines(output_cov_files, fileConn)
  close(fileConn)

  #Now loop over mask_list calling RAREMETAL each time
  for(i in 1:length(mask_list)){
    group_file_curr_name <- names(mask_list)[i]

    #Check to see if current mask is one created by the package
    #if(group_file_curr_name %in% c("pLOF", "pLOF_narrowMissense", "pLOF_broadMissense")){
    #  if(is.null(altGroupFilePath)){
    #    group_file_curr <- paste0(gene, ".", group_file_curr_name, ".group.file")
    #  }else{
    #    group_file_curr <- paste0(altGroupFilePath, ".", gene, ".", group_file_curr_name, ".group.file")
    #  }
    #}else{
    #  group_file_curr <- paste0(gene, ".", group_file_curr_name, ".group.file")
    #}

    #Get path to group file
    group_file_curr <- group_file_names[i]

    if(is.null(altRaremetalPath)){
      full_prefix <- group_file_curr_name
    }else{
      full_prefix <- paste0(altRaremetalPath, group_file_curr_name)
    }

    cm <- str_glue("raremetal --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --burden --maf 1 --hwe {hwe} --prefix {full_prefix}")
    system(cm)
  }

  #Remove the temporary covFiles and summaryFiles
  cm <- str_glue("rm {summary_files}")
  system(cm)

  cm <- str_glue("rm {cov_files}")
  system(cm)

}


