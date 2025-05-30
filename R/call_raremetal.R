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
#' @param burden Perform simple (unweighted) burden test? Must be TRUE or FALSE.
#'
#' @param wburden Perform weight burden test? Weights are from beta(MAF, 1, 25)
#' distribution. Must be TRUE or FALSE.
#'
#' @param SKAT Perform simple (unweighted) SKAT? Must be TRUE or FALSE.
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
#' @param altRaremetalName Optional command to call RAREMETAL with. For example if
#' you have multiple RAREMETAL installations and want to use a specific one, please
#' give the path + name of the executable file for the installation you want to use.
#' Default = NULL.
#'
#' @importFrom stringr str_glue
#'
#' @return Nothing. Will write the RAREMETAL results to
#' altRaremetalPath/gene.mask_name.meta.burden.results for burden,
#' altRaremetalPath/gene.mask_name.meta.BBeta.results for weighted burden,
#' and altRaremetalPath/gene.mask_name.meta.SKAT_.results for SKAT

call_raremetal <- function(mask_list, score_stat_file, group_file_names,
                           burden, wburden, SKAT,
                           altCovariancePath = NULL, altGroupFilePath = NULL,
                           altRaremetalPath = NULL, gene = 'gene',
                           hwe = 0.000001, altRaremetalName = NULL){

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

    if(!is.null(altRaremetalName)){
      raremetal <- altRaremetalName
    }else{
      raremetal <- "raremetal"
    }

    if(burden & !wburden &!SKAT){
      cm <- str_glue("{raremetal} --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --burden --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }else if(!burden & wburden & !SKAT){
      cm <- str_glue("{raremetal} --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --BBeta --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }else if(!burden & !wburden & SKAT){
      cm <- str_glue("{raremetal} --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --SKAT --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }else if(burden & wburden & !SKAT){
      cm <- str_glue("{raremetal} --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --burden --BBeta --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }else if(burden & !wburden & SKAT){
      cm <- str_glue("{raremetal} --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --burden --SKAT --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }else if(!burden & wburden & SKAT){
      cm <- str_glue("{raremetal} --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --BBeta --SKAT --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }else{
      cm <- str_glue("raremetal --summaryFiles {summary_files} --covFiles {cov_files} --groupFile {group_file_curr} --burden --BBeta --SKAT --maf 1 --hwe {hwe} --prefix {full_prefix}")
    }

    system(cm)
  }

  #Remove the temporary covFiles and summaryFiles
  cm <- str_glue("rm {summary_files}")
  system(cm)

  cm <- str_glue("rm {cov_files}")
  system(cm)

}


