#' @title Covariance computation using external reference panel
#'
#' @description A function that computes covariance for rare-variant aggregation
#' testing using a VCF from an external reference panel and write it
#'
#' @param vcf_file String containing the file path and name containing the VCF
#'  of the genotypes of the  set of individuals to be used as the reference
#'  panel.

#' @param altCovariancePath Option to give alternative path to write covariance
#' files to. Default = NULL. If left as NULL, covariance files will be written
#' in current directory.
#'
#' @param gene The name of the gene to use in the name of the covariance file.
#'  Default = 'gene'.
#'



