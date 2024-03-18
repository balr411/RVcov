#' @title Read in user-supplied group files
#'
#' @description Read in user-supplied group files
#'
#' @param group_file Vector of strings containing the file paths and names of
#'  the group files to perform the rare-variant aggregation tests with.
#'
#'  @return A list where each element is a vector of the variants in the
#'  user-supplied mask.
#'

read_user_group_file <- function(group_file){
  mask_list <- list()

  for(i in 1:length(group_file)){
    name_curr <- paste0("userMask", i)
    mask_curr <- readLines(group_file[i])
    mask_list[[i]] <- unlist(strsplit(mask_curr, split = "\t", fixed = FALSE))[-1]
    if(length(mask_list[[i]]) == 0){
      stop(paste0("Mask ", i, " has no variants. Are you sure your mask is in
                    RAREMETAL format?"))
    }
    names(mask_list)[i] <- name_curr
  }

  return(mask_list)

}
