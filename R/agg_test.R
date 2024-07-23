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
#'  specified, must supply anno_file (or ) and types of group files to use. Default =
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
#' @importFrom data.table fread
#' @importFrom stringr str_glue
#' @importFrom dplyr anti_join
#'
#' @return 1 if RAREMETAL was performed using the estimated covariance successfully,
#' 0 if two-stage approach was performed successfully but covariance did not need
#' to be computed.
#'
#' @examples #Use examples of data that I simulate (Note this will cause devtool::check() issue)

agg_test <- function(score_stat_file, vcf_file, chr, anno_file = NULL,
                     anno = NULL, two_stage = FALSE, two_stage_threshold = 3,
                     group_file = NULL, pLOF = FALSE,
                     pLOF_narrowMissense = FALSE, pLOF_broadMissense = FALSE,
                     altGroupFilePath = NULL, altCovariancePath = NULL,
                     altRaremetalPath = NULL, mafThreshold = 0.01, gene = NULL,
                     gene_start = NULL, gene_end = NULL, hwe = 0.000001){

  #Add checks here to make sure bcftools, RAREMETAL, and tabix are installed

  #First read in the necessary columns from the score statistic file
  if(!file.exists(score_stat_file)){
    stop("Score statistic file doesn't exist!")
  }

  allele_freq_test <- fread(cmd = str_glue("zgrep -v ^# {score_stat_file}"),
                            select = c(1, 2, 3, 4, 7, 8, 14, 17), data.table = FALSE)

  names(allele_freq_test) <- c("CHROM", "POS", "REF", "ALT", "AF", "AC",
                               "U_STAT", "PVAL")

  #Delete multi-alleleic variants
  allele_freq_test <- anti_join(allele_freq_test, allele_freq_test[duplicated(allele_freq_test$POS),], by = "POS")

  #Delete variants with AF > MAF threshold and AF = 0
  allele_freq_test <- allele_freq_test[allele_freq_test$AF > 0 & allele_freq_test$AF <= mafThreshold,]

  #Now read in the reference allele frequency file and take a list of the variants
  #that are multi-allelic in the test but not the reference
  #Note need bcftools installed here
  #Check that vcf exists
  if(!file.exists(vcf_file)){
    stop("VCF file doesn't exist!")
  }

  allele_freq_reference <- generate_allele_frequency(vcf_file, chr, gene_start, gene_end)
  dup_ref <- allele_freq_reference$POS[duplicated(allele_freq_reference$POS)]

  #Keep track of the variants orignially in the VCF in case the matrix is too
  #large to read in (this is also needed when the variants in the reference panel
  #don't match the variants in the test set)
  original_snps <- paste(allele_freq_reference$CHROM,  allele_freq_reference$POS, allele_freq_reference$REF, allele_freq_reference$ALT, sep = ":")
  #Add functionality for if the matrix is too big later (functionality for when
  #the variants don't match is already there)

  if(length(dup_ref) > 0){
    allele_freq_test <- allele_freq_test[!(allele_freq_test$POS %in% dup_ref),]
    allele_freq_reference <- allele_freq_reference[!(allele_freq_reference$POS %in% dup_ref),]
  }

  if(nrow(allele_freq_test) == 0){
    stop(paste0("No variants with 0 < MAF < ", mafThreshold))
  }

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
    mask_list_obj <- generate_group_file(anno, allele_freq_test, gene, pLOF,
                                     pLOF_narrowMissense, pLOF_broadMissense,
                                     altGroupFilePath, mafThreshold)

    mask_list <- mask_list_obj[[1]]
    file_paths <- mask_list_obj[[2]]

  }else if(is.null(group_file)){
    stop("Did not specify which group files to use and did not give any of your own")
  }else{
    #Initialize mask list to read in user-given group files
    mask_list <- list()
  }

  #If group files supplied, read those in
  if(!is.null(group_file)){
    mask_list <- c(mask_list, read_user_group_file(group_file))
    if(exists("file_paths")){
      file_paths <- c(file_paths, group_file)
    }else{
      file_paths <- c(group_file)
    }

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
    #Now make sure to match based on SNP and not position (and fix allele switches)
    #First match on position
    allele_freq_reference <- allele_freq_reference[allele_freq_reference$POS %in% allele_freq_test$POS,]

    #Create copy of test allele freq that only contains the variants from the reference panel
    allele_freq_test_red_temp <- allele_freq_test[allele_freq_test$POS %in% allele_freq_reference$POS,]

    ###############################################################################
    ## The next few lines on allele switching should be put into a function of its own
    #Find allele switches
    idx_switched <- which((allele_freq_reference$REF == allele_freq_test_red_temp$ALT) & (allele_freq_reference$ALT == allele_freq_test_red_temp$REF))

    #Track the allele switches (to change them in the genotype matrix later)
    allele_freq_reference$ALLELE_SWITCH <- 0

    if(length(idx_switched) > 0){
      allele_freq_reference$ALLELE_SWITCH[idx_switched] <- 1

      #Fix the allele switches
      allele_freq_reference$REF[idx_switched] <- allele_freq_test_red_temp$REF[idx_switched]
      allele_freq_reference$ALT[idx_switched] <- allele_freq_test_red_temp$ALT[idx_switched]

      allele_freq_reference$AC[idx_switched] <- allele_freq_reference$AN[idx_switched] - allele_freq_reference$AC[idx_switched]
      allele_freq_reference$AF[idx_switched] <- 1 - allele_freq_reference$AF[idx_switched]
    }

    #Now match based on SNP
    allele_freq_test$SNP <-  paste(allele_freq_test$CHROM, allele_freq_test$POS, allele_freq_test$REF, allele_freq_test$ALT, sep = ":")
    allele_freq_reference$SNP <- paste(allele_freq_reference$CHROM,  allele_freq_reference$POS, allele_freq_reference$REF, allele_freq_reference$ALT, sep = ":")

    allele_freq_reference <- allele_freq_reference[allele_freq_reference$SNP %in% allele_freq_test$SNP,]

    #Note that here the dimension of the reference allele frequency data frame must be
    #<= dimension of the test allele frequency data
    if(nrow(allele_freq_reference) > nrow(allele_freq_test)){
      stop("Dimension of the reference allele frequency data frame is more than the
           dimension of the test allele frequency data frame after matching")
    }

    #Read in VCF file
    gt <- read_vcf(vcf_file, chr, allele_freq_test, allele_freq_reference,
                   original_snps, gene_start, gene_end)

    #Compute covariance
    compute_covariance(gt, allele_freq_test, n, resid_var,
                       altCovariancePath, gene) #currently not getting the same with InPSYght so check this

    #Call RAREMETAL
    call_raremetal(mask_list, score_stat_file, file_paths, altCovariancePath,
                   altGroupFilePath, altRaremetalPath, gene, hwe) #Have to check what happens if RAREMETAL does not run successfully

    to_return <- 1
  }else{
    to_return <- 0
  }

  return(to_return)
}



