#' @title Read in a VCF to a genotype matrix
#'
#' @description A function that takes a VCF file and returns a genotype matrix with
#' the genotypes of the individuals, colnames corresponding to the variant positions,
#' and rownames corresponding to the variant names.
#'
#' @param vcf_file String containing the file path and name containing the VCF
#'  of the genotypes of the  set of individuals to be used as the reference
#'  panel. Note that the chromosome column must be in the format c{chr_number}.
#'  ie. c1 for chromosome 1.
#'
#' @param chr The chromosome number on which you want to perform the analysis.
#' Must be in the format c{chr_number}. ie. c1 for chromosome 1.
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
#'
#' @return A data frame with the genotypes of the individuals and colnames
#' corresponding to the variant names.

read_vcf <- function(vcf_file, chr, gene_start = NULL, gene_end = NULL){

  chr_num <- as.numeric(unlist(strsplit(chr, split = "c"))[2])

  if(is.null(gene_start) & is.null(gene_end)){
    gt_current <- readGT(file = vcf_file)
  }else if((is.null(gene_start) & !is.null(gene_end)) | (!is.null(gene_start) & is.null(gene_end))){
    stop("Only one gene boundary given.")
  }else{
    current_range <- ScanVcfParam(which = GRanges(seqnames = Rle(chr_num), ranges = IRanges(start = gene_start, end = gene_end)))
  }




}

