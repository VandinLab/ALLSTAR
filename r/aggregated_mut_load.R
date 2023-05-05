aggregated_mut_load <- function() {
  
  aggregated_mut <- read_delim("data/som_LOH_meth_aggregated.csv")
  
  colnames(aggregated_mut)[2:dim(aggregated_mut)[2]] <- paste(colnames(aggregated_mut)[2:dim(aggregated_mut)[2]], "_alt", sep="")
  
  return(aggregated_mut)
  
}