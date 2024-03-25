#' @title Read in a VCF to a sparse genotype matrix
#'
#' @description A function that takes a VCF file and returns a genotype matrix with
#' the genotypes of the individuals, colnames corresponding to the variant positions,
#' and rownames corresponding to the ID of the individuals. Only returns genotypes
#' of variants that are bi-allelic and have allele frequency > 0 in the test set.
#' ie. the ones that are present in the allele frequency files.
#'
#' @param vcf_file String containing the file path and name containing the VCF
#'  of the genotypes of the  set of individuals to be used as the reference
#'  panel. Note that the chromosome column must be in the format c{chr_number}.
#'  ie. c1 for chromosome 1.
#'
#' @param chr The chromosome number on which you want to perform the analysis.
#' Must be in the format c{chr_number}. ie. c1 for chromosome 1.
#'
#' @param allele_freq_reference Data frame from reference panel VCF containing
#' the columns CHROM, POS, REF, ALT, AN, AC, AF
#'
#' @param gene_start The starting base pair position of the gene on the chromosome.
#' Default = NULL.
#'
#' @param gene_end The ending base pair position of the gene on the chromosome.
#' Default = NULL.
#'
#' @importFrom {VariantAnnotation} {ScanVcfParam}
#' @importFrom {VariantAnnotation} {readGT}
#' @importFrom {GenomicRanges} {GRanges}
#' @importFrom {IRanges} {IRanges}
#' @importFrom {dplyr} {anti_join}
#' @importFrom {Matrix} {Matrix}
#'
#' @return A sparse genotype matrix with the genotypes of the individuals,
#' colnames corresponding to the variant names, and rownames corresponding to
#' the IDs

read_vcf <- function(vcf_file, chr, allele_freq_test, allele_freq_reference,
                     gene_start = NULL, gene_end = NULL){

  chr_num <- as.numeric(unlist(strsplit(chr, split = "c"))[2])

  if(is.null(gene_start) & is.null(gene_end)){
    gt_current <- readGT(file = vcf_file)
  }else if((is.null(gene_start) & !is.null(gene_end)) | (!is.null(gene_start) & is.null(gene_end))){
    stop("Only one gene boundary given.")
  }else{
    current_range <- ScanVcfParam(which = GRanges(seqnames = Rle(chr_num), ranges = IRanges(start = gene_start, end = gene_end)))
    gt_current <- readGT(file = vcf_file, param = current_range)
  }

  gt_current <- t(gt_current)
  vars_gt_current <- colnames(gt_current)
  id_names <- rownames(gt_current)
  vars_gt_current_positions <- as.numeric(sapply(vars_gt_current, function(s) as.numeric(unlist(strsplit(s, split=":"))[2])))

  #Keep only single-allelic variants that are in the test set
  gt_current <- gt_current[,vars_gt_current_positions %in% allele_freq_test$POS, drop = FALSE]

  if(dim(gt_current)[2] > 0){ #Should be >= 1 variant in the gene
    allele_freq_reference_temp <- allele_freq_reference[allele_freq_reference$POS %in% vars_gt_current_positions,]
    allele_freq_test_temp <- allele_freq_test[allele_freq_test$POS %in% vars_gt_current_positions,]

    gt_current[gt_current %in% c("0/0", "0|0")] <- 0
    gt_current[gt_current %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
    gt_current[gt_current %in% c("1/1", "1|1")] <- 2

    n_row <- nrow(gt_current)
    n_col <- ncol(gt_current)

    gt_current <- suppressWarnings(Matrix(as.numeric(gt_current), nrow = n_row, ncol = n_col, sparse = TRUE)) #Note the NAs for missing values; expect warning here so suppress it

    gt_current <- sweep(gt_current, 2, 2*allele_freq_reference_temp$AF) #subtract the allele frequencies rowwise

    #Now impute the missing values
    is_na <- which(is.na(gt_current), arr.ind = TRUE)
    gt_current[is_na] <- colMeans(gt_current, na.rm = TRUE)[is_na[,"col"]]

    #If there are still NAs left, should make them 0 (these correspond to variants with no calls)
    gt_current[is_na] <- 0
  }else{
    stop("No variants in this gene!")
  }

  #Add back in the row and column names
  colnames(gt_current) <- vars_gt_current[vars_gt_current_positions %in% allele_freq_test$POS]
  rownames(gt_current) <- id_names

  return(gt_current)
}

