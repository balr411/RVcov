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
#' Must be in the format c{chr_number} or chr{chr_number}. ie. c1 or chr1 for
#' chromosome 1.
#'
#' @param allele_freq_reference Data frame from reference panel VCF containing
#' the columns CHROM, POS, REF, ALT, AN, AC, AF
#'
#' @param original_snps A vector containing the SNPs that were in the original
#' reference allele frequency file before filtering/matching to the test set
#'
#' @param gene_start The starting base pair position of the gene on the chromosome.
#' Default = NULL.
#'
#' @param gene_end The ending base pair position of the gene on the chromosome.
#' Default = NULL.
#'
#' @importFrom VariantAnnotation ScanVcfParam
#' @importFrom VariantAnnotation readGT
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors Rle
#' @importFrom dplyr anti_join
#' @importFrom Matrix Matrix
#' @importFrom Matrix as.matrix
#'
#' @return A sparse genotype matrix with the genotypes of the individuals,
#' colnames corresponding to the SNPs in format CHR:REF:ALT:POS, and rownames
#' corresponding to the participant IDs

read_vcf <- function(vcf_file, chr, allele_freq_test, allele_freq_reference,
                     original_snps, gene_start = NULL, gene_end = NULL){

  chr_num <- gsub("c", "", chr)
  chr_num <- as.numeric(gsub("hr", "", chr_num))

  if(is.null(gene_start) & is.null(gene_end)){
    gt_current <- readGT(file = vcf_file)
  }else if((is.null(gene_start) & !is.null(gene_end)) | (!is.null(gene_start) & is.null(gene_end))){
    stop("Only one gene boundary given.")
  }else{
    current_range <- ScanVcfParam(which = GRanges(seqnames = S4Vectors::Rle(chr_num), ranges = IRanges(start = gene_start, end = gene_end)))
    gt_current <- readGT(file = vcf_file, param = current_range)
  }

  gt_current <- t(gt_current)
  id_names <- rownames(gt_current)

  #Keep only single-allelic variants that are in the test set - here I think we should be doing it from the reference set since we matched based on SNP not pos
  gt_current <- gt_current[,original_snps %in% allele_freq_reference$SNP, drop = FALSE] # Note there could be no variants shared between the test and reference

  #Now switch the entries for the alleles that had switched alleles
  idx_gt_switched <- which(allele_freq_reference$ALLELE_SWITCH == 1)

  if(length(idx_gt_switched) > 0){
    gt_current[,idx_gt_switched][gt_current[,idx_gt_switched] %in% c("0/0", "0|0")] <- 2
    gt_current[,idx_gt_switched][gt_current[,idx_gt_switched] %in% c("1/1", "1|1")] <- 0
  }


  if(dim(gt_current)[2] > 0){ #Should be >= 1 variant in the gene
    allele_freq_reference_temp <- allele_freq_reference
    allele_freq_test_temp <- allele_freq_test ###

    gt_current[gt_current %in% c("0/0", "0|0")] <- 0
    gt_current[gt_current %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
    gt_current[gt_current %in% c("1/1", "1|1")] <- 2

    n_row <- nrow(gt_current)
    n_col <- ncol(gt_current)

    gt_current <- suppressWarnings(Matrix::Matrix(as.numeric(gt_current), nrow = n_row, ncol = n_col, sparse = TRUE)) #Note the NAs for missing values; expect warning here so suppress it

    gt_current <- sweep(gt_current, 2, 2*allele_freq_reference_temp$AF) #subtract the allele frequencies rowwise

    #Now impute the missing values
    is_na <- which(is.na(gt_current), arr.ind = TRUE)
    gt_current[is_na] <- colMeans(gt_current, na.rm = TRUE)[is_na[,"col"]]

    #If there are still NAs left, should make them 0 (these correspond to variants with no calls)
    gt_current[is_na] <- 0

    #Now add in the variants that weren't in the reference panel at all
    ## Need to come up with way to add the missing variants to the matrix
    gt <- matrix(0, nrow = nrow(gt_current), ncol = nrow(allele_freq_test))

    #Match columns
    idx_in_reference_panel <- match(allele_freq_reference$SNP, allele_freq_test$SNP)

    gt[,idx_in_reference_panel] <- Matrix::as.matrix(gt_current)
    gt <- Matrix::Matrix(gt, sparse = TRUE)

  }else{
    stop("No variants in this gene!")
  }

  #Check that the test allele frequency data frame and genotype matrix have the
  #same number of variants - shouldn't ever happen
  if(nrow(allele_freq_test) != ncol(gt)){
    stop("Allele frequency data frame and genotype matrix have differing number of variants.")
  }

  #Add back in the row and column names
  colnames(gt) <- allele_freq_test$SNP
  rownames(gt) <- id_names


  return(gt)
}

