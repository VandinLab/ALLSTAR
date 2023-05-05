triple_negative_load <- function() {
  
  histologic_annotation <- read_excel("data/TCGA_histologic_type_annotation.xlsx", na = "NA")
  triple_negative <- histologic_annotation %>% 
    select(`CLID`, `Triple Negative Status`)
  colnames(triple_negative) <- c("TCGA_ID", "triple_negative")
  triple_negative$TCGA_ID <- substr(triple_negative$TCGA_ID, 1, nchar(triple_negative$TCGA_ID)-4)
  
  triple_negative <- triple_negative[!(triple_negative$triple_negative=="N/A" | triple_negative$triple_negative=="Not Applicable"),]
  triple_negative$triple_negative <- ifelse(triple_negative$triple_negative == "No", 1, 2) 
  
  return(triple_negative)
  
}