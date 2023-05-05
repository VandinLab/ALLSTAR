LOH_load <- function() {

  LOH <- read.delim("data/BRCA_LOH.tsv")
  
  # add germline identifier
  colnames(LOH)[2:dim(LOH)[2]] <- paste0(colnames(LOH)[2:dim(LOH)[2]], "_LOH")
  
  LOH[is.na.data.frame(LOH)] <- 0
  
  return(LOH)
}