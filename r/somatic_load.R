somatic_load <- function() {
  
  somatic_mut <- read_table("data/BRCA_mutated_genes.maf")
  
  # clean TCGA_ID
  somatic_mut$TCGA_ID <- substr(somatic_mut$TCGA_ID, 1, 12)
  
  # remove duplicates
  somatic_mut <- somatic_mut[!duplicated(somatic_mut), ]
  
  # add somatic identifier
  somatic_mut$gene <- paste(somatic_mut$gene, "_somatic", sep="")
  # spread info on gene level
  somatic_mut$mut_presence <- 1
  somatic_mut <- somatic_mut %>% spread(gene, mut_presence, fill = 0)
  
  return(somatic_mut)
}