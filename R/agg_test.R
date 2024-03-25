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
#'  panel. Note that the chromosome column must be in the format c{chr_number}.
#'  ie. c1 for chromosome 1.
#'
#' @param chr The chromosome number on which you want to perform the analysis.
#' Must be in the format c{chr_number}. ie. c1 for chromosome 1.
#'
#' @param anno_file String containing the file path and name containing the
#'  VEP annotation file to create the group files. If not specified, must give
#'  argument to group_file. Default = NULL.
#'
#' @param anno Option to directly pass VEP annotation file to the function.
#' Default = NULL. Do not pass both anno_file and anno to the function as it will
#' cause an error.
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
#' @param altGroupFilePath If pLOF, pLOF_narrowMissense, or pLOF_broadMissense
#'  are TRUE, option to give alternative path to write group files to. Default =
#'  NULL. If left as NULL, group files will be written in current directory.
#'
#' @param altCovariancePath Option to give alternative path to write covariance
#' files to. Default = NULL. If left as NULL, covariance files will be written
#' in current directory.
#'
#' @param altRaremetalPath Optional alternative path to write RAREMETAL results
#' to. Default = NULL. If left as NULL, will write in current directory.
#'
#' @param mafThreshold If pLOF, pLOF_narrowMissense, or pLOF_broadMissense
#'  are TRUE, threshold for including variants in group file. Note here
#'  we assume that the alternate allele is coded as the minor allele.
#'  Default = 0.01.
#'
#' @param gene The name of the gene to use for annotation files if group files
#'  are not already supplied. Default = NULL.
#'
#' @param gene_start The starting base pair position of the gene on the chromosome.
#' Default = NULL. If you have large reference panel VCFs it is recommended to
#' supply gene_start and gene_end.
#'
#' @param gene_end The ending base pair position of the gene on the chromosome.
#' Default = NULL. If you have large reference panel VCFs it is recommended to
#' supply gene_start and gene_end.
#'
#' @param hwe Hardy-Weinberg equilibrium p-value cut-off for rare-variant
#' aggregation testing. Default = 0.000001.
#'
#' @importFrom {data.table} {fread}
#' @importFrom {stringr} {str_glue}
#'
#' @return 1 if RAREMETAL was performed using the estimated covariance successfully,
#' 0 if two-stage approach was performed successfully but covariance did not need
#' to be computed.
#'
#' @examples Use examples of data that I simulate

agg_test <- function(score_stat_file, vcf_file, chr, anno_file = NULL,
                     anno = NULL, two_stage = FALSE, two_stage_threshold = 3,
                     group_file = NULL, pLOF = FALSE,
                     pLOF_narrowMissense = FALSE, pLOF_broadMissense = FALSE,
                     altGroupFilePath = NULL, altCovariancePath = NULL,
                     altRaremetalPath = NULL, mafThreshold = 0.01, gene = NULL,
                     gene_start = NULL, gene_end = NULL, hwe = 0.000001){

  #Add checks here to make sure bcftools and RAREMETAL are installed

  #First read in the necessary columns from the score statistic file
  if(!file.exists(score_stat_file)){
    stop("Score statistic file doesn't exist!")
  }

  allele_freq_test <- fread(cmd = str_glue("zgrep -v ^# {score_stat_file}"),
                            select = c(1, 2, 3, 4, 7, 8, 14, 17), data.table = FALSE)

  names(allele_freq_test) <- c("CHROM", "POS", "REF", "ALT", "AF", "AC",
                               "U_STAT", "PVAL")

  #Now get residual variance and sample size
  ss_residual_variance_list <- ss_residual_variance(score_stat_file)
  n <- as.numeric(ss_residual_variance_list$sampleSize)
  resid_var <- as.numeric(ss_residual_variance_list$residualVariance)

  #Get group files
  #Return them as a list of data frames in addition to writing the group file to
  #make it easier to perform the two-stage approach
  if(pLOF | pLOF_narrowMissense | pLOF_broadMissense){
    if(is.null(anno_file) & is.null(anno)){
      stop("Did not give annotation file or path to annotation file needed
           to create the desired group files. If you have supplied your own
           group files and just want to use those, set pLOF, pLOF_narrowMissense,
           and pLOF_broadMissense to FALSE")
    }else if(is.null(gene)){
      stop("You must supply a gene name to create the desired group files.
           If you have supplied your own group files and just want to use those,
           please set pLOF, pLOF_narrowMissense, and pLOF_broadMissense to FALSE.")
    }else if(!is.null(anno_file) & !is.null(anno)){
      stop("You have supplied both an annotation file and a path to an annotation
           file. Please supply only one of these.")
    }else if(!is.null(anno_file) & is.null(anno)){
      #Read in annotation file if not already given
      anno <- fread(cmd = paste0("zgrep -v ^## ", anno_file), header = T, data.table = F)
    }


    #Create group files
    mask_list <- generate_group_file(anno, allele_freq_test, gene, pLOF,
                                     pLOF_narrowMissense, pLOF_broadMissense,
                                     altGroupFilePath, mafThreshold)

  }else if(is.null(group_file)){
    stop("Did not specify which group files to use and did not give any of your own")
  }else{
    #Initialize mask list to read in user-given group files
    mask_list <- list()
  }

  #If group files supplied, read those in
  if(!is.null(group_file)){
    mask_list <- c(mask_list, read_user_group_file(group_file))
  }

  #Perform two-stage approach if desired
  if(two_stage){
    calculate_covariance <- two_stage_test(mask_list, allele_freq_test,
                                           two_stage_threshold,
                                           residual_variance = resid_var,
                                           sample_size = n)
  }else{
    calculate_covariance <- TRUE
  }

  #If necessary, calculate covariance using external reference panel and perform RAREMETAL
  if(calculate_covariance){
    #Check that vcf exists
    if(!file.exists(vcf_file)){
      stop("VCF file doesn't exist!")
    }

    #First create allele frequency data frame for reference panel - note need bcftools installed here
    allele_freq_reference <- generate_allele_frequency(vcf_file, chr, gene_start, gene_end)

    #Remove multi-allelic variants in test set
    allele_freq_test <- anti_join(allele_freq_test, allele_freq_test[duplicated(allele_freq_test$POS),], by = "POS")

    #remove variants with 0 AF in test set
    allele_freq_test <- allele_freq_test[allele_freq_test$AF > 0,]

    #Keep track of the variants orignially in the VCF in case the matrix is too
    #large to read in
    #original_vars <- allele_freq_reference$POS
    #Add functionality for this later

    #Now reduce ref allele freq to only include variants found in test set
    allele_freq_reference <- allele_freq_reference[allele_freq_reference$POS %in% allele_freq_test$POS,]

    #Read in VCF file
    gt <- read_vcf(vcf_file, chr, allele_freq_test, allele_freq_reference,
                   gene_start, gene_end)

    #Compute covariance
    compute_covariance(gt, allele_freq_test, n, resid_var,
                       altCovariancePath, gene)

    #Call RAREMETAL
    call_raremetal(mask_list, score_stat_file, altCovariancePath,
                   altGroupFilePath, altRaremetalPath, gene, hwe) #Have to check what happens if RAREMETAL does not run successfully

    to_return <- 1
  }else{
    to_return <- 0
  }

  return(to_return)
}



