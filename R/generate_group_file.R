#' @title Generate group/mask file for rare-variant aggregation testing
#'
#' @description A function that generates and writes to disk in RAREMETAL format
#'  a group/mask file for rare variant aggregation testing.
#'
#' @param anno A VEP annotation file with at minimum columns named #Uploaded_variation
#'  containing the variants in format chr:pos:ref:alt, Consequence containing the predicted
#'  consequences of the variants, SYMBOL containing the gene names that the variants
#'  belong to, and PICK which is 1 if that row is to be used and "-" otherwise,
#'  a column sumAlgs which is 1 if deleterious and 0 otherwise OR LRT_pred,
#'  MutationTaster_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, and
#'  SIFT4G_pred containing LRT, MutationTaster, Polyphen2_HDIV, Polyphen2_HVAR,
#'  and SIFT4G deleteriousness predictions respectively.
#'
#' @param allele_freq The in-sample allele frequency file with columns CHROM,
#'  POS, REF, ALT, AF, AC, and PVAL.
#'
#' @param gene The name of the gene from the annotation file to use.
#'
#' @param pLOF TRUE or FALSE. Perform aggregation tests using pLOF variants only?
#'  Default = TRUE.
#'
#' @param pLOF_narrowMissense TRUE or FALSE. Perform aggregation tests using
#'  pLOF + missense(narrow) variants only? Default = TRUE.
#'
#' @param pLOF_broadMissense TRUE or FALSE. Perform aggregation tests using
#'  pLOF + missense(broad) variants only? Default = TRUE.
#'
#' @param altGroupFilePath Optional alterntive path to write group files to.
#'  Default = NULL.
#'
#' @param mafThreshold Threshold for including variants in group file. Note here
#'  we assume that the alternate allele is coded as the minor allele.
#'  Default = 0.01.
#'
#' @importFrom {dplyr} {anti_join}
#' @importFrom {stringr} {str_count}
#'
#' @return List of data frames containing the SNPs in each mask
#'

generate_group_file <- function(anno, allele_freq, gene, pLOF = TRUE,
                                pLOF_narrowMissense = TRUE,
                                pLOF_broadMissense = TRUE,
                                altGroupFilePath = NULL, mafThreshold = 0.01){

  anno <- anno[anno$PICK == 1,]

  #Delete multi-alleleic variants
  allele_freq <- anti_join(allele_freq, allele_freq[duplicated(allele_freq$POS),], by = "POS")

  #Delete variants with AF > MAF threshold and AF = 0
  allele_freq <- allele_freq[allele_freq$AF > 0 & allele_freq$AF <= mafThreshold,]

  if(nrow(allele_freq) == 0){
    stop(paste0("No variants with 0 < MAF < ", mafThreshold))
  }

  allele_freq$UploadedVariation <- paste(allele_freq$CHROM, allele_freq$POS,
                                         allele_freq$REF, allele_freq$ALT, sep = ":")


  #Match variants in annotation file to those in allele frequency file
  idx_common_anno <- which(anno$'#Uploaded_variation' %in% allele_freq$UploadedVariation)
  anno <- anno[idx_common_anno,]

  allele_freq_plof <- allele_freq_narrowMissense <- allele_freq_broadMissense <- allele_freq

  #Find pLOF variants and missense variants
  cons <- unique(anno$Consequence)
  imp_vars <- c()
  imp_vars_missense <- c()
  bad_cons <- c("stop_gained", "splice_acceptor_variant", "stop_lost", "splice_donor_variant", "frameshift_variant", "start_lost")
  if(length(cons) > 0){
    for(i in 1:length(cons)){
      temp <- unlist(strsplit(cons[i], split = ","))
      if(sum(bad_cons %in% temp) > 0){
        imp_vars <- c(imp_vars, cons[i])
      }

      if((sum(bad_cons %in% temp) == 0) & (sum(temp == "missense_variant") > 0)){
        imp_vars_missense <- c(imp_vars_missense, cons[i])
      }
    }
  }

  #Create deleteriousness prediction based on LRT, MutationTaster, Polyphen2_HDIV,
  #Polyphen2_HVAR, and SIFT4G
  if(pLOF_narrowMissense){
    if("sumAlgs" %in% names(anno)){
      warning("sumAlgs column given in VEP annotation, will use this for
              deleteriousness prediction")
    }else{
      anno$LRT_pred <- str_count(anno$LRT_pred, "D")
      anno$MutationTaster_pred <- str_count(anno$MutationTaster_pred, "D") + str_count(anno$MutationTaster_pred, "A")
      anno$Polyphen2_HDIV_pred <- str_count(anno$Polyphen2_HDIV_pred, "D") + str_count(anno$Polyphen2_HDIV_pred, "P")
      anno$Polyphen2_HVAR_pred <- str_count(anno$Polyphen2_HVAR_pred, "D") + str_count(anno$Polyphen2_HVAR_pred, "P")
      anno$SIFT4G_pred <- str_count(anno$SIFT4G_pred, "D")
      anno$LRT_pred <- ifelse(anno$LRT_pred > 0, 1, 0)
      anno$MutationTaster_pred <- ifelse(anno$MutationTaster_pred > 0, 1, 0)
      anno$Polyphen2_HDIV_pred <- ifelse(anno$Polyphen2_HDIV_pred > 0, 1, 0)
      anno$Polyphen2_HVAR_pred <- ifelse(anno$Polyphen2_HVAR_pred > 0, 1, 0)
      anno$SIFT4G_pred <- ifelse(anno$SIFT4G_pred > 0, 1, 0)
      anno$Sum_Algs <- rowSums(anno[,c("LRT_pred", "MutationTaster_pred", "Polyphen2_HDIV_pred",
                                       "Polyphen2_HVAR_pred", "SIFT4G_pred")])
      anno$Sum_Algs <- ifelse(anno$Sum_Algs == 5, 1, 0)
    }
  }

  #List to return to agg_test
  to_return <- list()

  if(pLOF){
    idx_anno_plof <- which(anno$Consequence %in% imp_vars)
    anno_plof <- anno[idx_anno_plof,]

    mask_list <- as.list(gene)

    to_add <- unique(anno_plof[anno_plof$SYMBOL == gene,]$'#Uploaded_variation')
    if(length(to_add) > 0){
      if(sum(is.na(to_add)) == 0){
        mask_list[[1]] <- c(mask_list[[1]], to_add)
      }
    }

    to_write <- vector(length = 1)
    to_write[1] <- paste(mask_list[[1]], collapse = "\t")

    if(is.null(altGroupFilePath)){
      output_plof <- paste0(gene, ".pLOF.group.file")
    }else{
      output_plof <- paste0(altGroupFilePath, ".", gene, ".pLOF.group.file")
    }

    fileConn <- file(output_plof)
    writeLines(to_write, fileConn)
    close(fileConn)

    #Create data frame to pass back to agg_test
    to_return$pLOF <- data.frame(snp = mask_list[[1]][-1])
  }

  if(pLOF_narrowMissense){
    idx_anno_narrowMissense <- c(which(anno$Consequence %in% imp_vars),
                                 which((anno$Consequence %in% imp_vars_missense) & (anno$Sum_Algs == 1)))
    idx_anno_narrowMissense <- idx_anno_narrowMissense[order(idx_anno_narrowMissense)]

    anno_narrowMissense <- anno[idx_anno_narrowMissense,]

    mask_list <- as.list(gene)

    to_add <- unique(anno_narrowMissense[anno_narrowMissense$SYMBOL == gene,]$'#Uploaded_variation')

    if(length(to_add) > 0){
      if(sum(is.na(to_add)) == 0){
        mask_list[[1]] <- c(mask_list[[1]], to_add)
      }
    }

    to_write <- vector(length = 1)
    to_write[1] <- paste(mask_list[[1]], collapse = "\t")

    if(is.null(altGroupFilePath)){
      output_narrowMissense <- paste0(gene, ".pLOF_narrowMissense.group.file")
    }else{
      output_narrowMissense <- paste0(altGroupFilePath, ".", gene, ".pLOF_narrowMissense.group.file")
    }

    fileConn <- file(output_narrowMissense)
    writeLines(to_write, fileConn)
    close(fileConn)

    #Create data frame to pass back to agg_test
    to_return$pLOF_narrowMissense <- data.frame(snp = mask_list[[1]][-1])
  }

  if(pLOF_broadMissense){
    idx_anno_broadMissense <- c(which(anno$Consequence %in% imp_vars),
                                which(anno$Consequence %in% imp_vars_missense))
    idx_anno_broadMissense <- idx_anno_broadMissense[order(idx_anno_broadMissense)]

    anno_broadMissense <- anno[idx_anno_broadMissense,]

    mask_list <- as.list(gene)

    to_add <- unique(anno_broadMissense[anno_broadMissense$SYMBOL == gene,]$'#Uploaded_variation')

    if(length(to_add) > 0){
      if(sum(is.na(to_add)) == 0){
        mask_list[[1]] <- c(mask_list[[1]], to_add)
      }
    }

    to_write <- vector(length = 1)
    to_write[1] <- paste(mask_list[[1]], collapse = "\t")

    if(is.null(altGroupFilePath)){
      output_broadMissense <- paste0(gene, ".pLOF_broadMissense.group.file")
    }else{
      output_broadMissense <- paste0(altGroupFilePath, ".", gene, ".pLOF_broadMissense.group.file")
    }

    fileConn <- file(output_broadMissense)
    writeLines(to_write, fileConn)
    close(fileConn)

    df_plof_broadMissense <- data.frame(snp = mask_list[[1]][-1])

    #Create data frame to pass back to agg_test
    to_return$pLOF_broadMissense <- data.frame(snp = mask_list[[1]][-1])
  }

  return(to_return)

}
