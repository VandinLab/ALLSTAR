clinical_load <- function() {
  
  molecular_annotation <- read_excel("data/TCGA_BRCA_histologic_features.xlsx", na = "NA")
  
  # why not taking these extra clinical features? can be useful
  clinical_extra <- molecular_annotation %>% dplyr::select(-c(`2016 Histology Annotations`, 
                                                              `Components if mixed`, 
                                                              `PAM50 and Claudin-low (CLOW) Molecular Subtype`))
  colnames(clinical_extra) <- c("TCGA_ID", "epithelial_area", "inflammation", 
                                "LCIS", "apocrine_features", "DCIS", "epithelial_tubule_formation", 
                                "LVI", "necrosis", "nuclear_pleomorphism", "fibrous_focus", "mitosis")
  
  
  # Clinical info -----------------------------------------------------------
  
  clinical <- read.delim("data/BRCA_ClinicalRawData.tsv", na.strings = "NA")
  colnames(clinical)[1] <- c("TCGA_ID")
  
  # add previous clinical info
  clinical <- left_join(clinical, clinical_extra, by = "TCGA_ID")
  
  ## group more values under the same roof:
  # race
  clinical$race[is.na(clinical$race)] <- "UNKNOWN"
  clinical$race[clinical$race == "BLACK OR AFRICAN AMERICAN"] <- "BLACK OR AFRICAN AMERICAN OR AMERICAN INDIAN"
  clinical$race[clinical$race == "AMERICAN INDIAN OR ALASKA NATIVE"] <- "BLACK OR AFRICAN AMERICAN OR AMERICAN INDIAN"
  # age at diagnosis
  clinical$age_at_diagnosis <- Hmisc::cut2(clinical$age_at_diagnosis, cuts = c(37,61))
  # Stage
  clinical$ajcc_pathologic_tumor_stage[clinical$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB", "Stage Tis")] <- "Stage I"
  clinical$ajcc_pathologic_tumor_stage[clinical$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB")] <- "Stage II"
  clinical$ajcc_pathologic_tumor_stage[clinical$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC")] <- "Stage III"
  clinical$ajcc_pathologic_tumor_stage[clinical$ajcc_pathologic_tumor_stage %in% c("Stage X")] <- "Not evaluated"
  clinical$ajcc_pathologic_tumor_stage[is.na(clinical$ajcc_pathologic_tumor_stage)] <- "Not evaluated"
  # Pathologic PT
  clinical$ajcc_tumor_pathologic_pt[clinical$ajcc_tumor_pathologic_pt %in% c("T1", "T1a", "T1b", "T1c", "TX")] <- "T1"
  clinical$ajcc_tumor_pathologic_pt[clinical$ajcc_tumor_pathologic_pt %in% c("T2", "T2a", "T2b")] <- "T2"
  clinical$ajcc_tumor_pathologic_pt[clinical$ajcc_tumor_pathologic_pt %in% c("T3", "T3a")] <- "T3"
  clinical$ajcc_tumor_pathologic_pt[clinical$ajcc_tumor_pathologic_pt %in% c("T4", "T4b", "T4d")] <- "T4"
  # Nodes PN
  clinical$ajcc_nodes_pathologic_pn[clinical$ajcc_nodes_pathologic_pn %in% c("N1", "N1a", "N1b", "N1c", "N1mi")] <- "N1"
  clinical$ajcc_nodes_pathologic_pn[clinical$ajcc_nodes_pathologic_pn %in% c("N2", "N2a")] <- "N2"
  clinical$ajcc_nodes_pathologic_pn[clinical$ajcc_nodes_pathologic_pn %in% c("N3", "N3a", "N3b", "N3c")] <- "N3"
  clinical$ajcc_nodes_pathologic_pn[clinical$ajcc_nodes_pathologic_pn %in% c("N0", "N0 (i-)", "N0 (i+)", "N0 (mol+)", "NX")] <- "N0"
  # Metastasis
  clinical$ajcc_metastasis_pathologic_pm[clinical$ajcc_metastasis_pathologic_pm %in% c("M0", "cM0 (i+)")] <- "M0"
  # Histological type
  clinical$histological_type[clinical$histological_type %in% c("Infiltrating Carcinoma NOS", 
                                                               "Medullary Carcinoma",
                                                               "Metaplastic Carcinoma",
                                                               "Mixed Histology (please specify)",
                                                               "Mucinous Carcinoma",
                                                               "Other  specify")] <- "Other (infiltrating, medullary, metaplastic, ...)"
  clinical$histological_type[is.na(clinical$histological_type)] <- "Other (infiltrating, medullary, metaplastic, ...)"
  # Menopause
  clinical$menopause_status[clinical$menopause_status %in% c("Peri (6-12 months since last menstrual period)",
                                                             "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)")] <- "Post"
  clinical$menopause_status[is.na(clinical$menopause_status)] <- "Indeterminate (neither Pre or Postmenopausal)"
  clinical$menopause_status[clinical$menopause_status == "Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement)"] <- "Pre"
  # Margin status
  clinical$margin_status[clinical$margin_status == "Close"] <- "Indeterminate"
  clinical$margin_status[is.na(clinical$margin_status)] <- "Indeterminate"
  # ER IHC
  clinical$er_status_by_ihc[is.na(clinical$er_status_by_ihc)] <- "Indeterminate"
  # PR IHC
  clinical$pr_status_by_ihc[is.na(clinical$pr_status_by_ihc)] <- "Indeterminate"
  # HER2 IHC
  clinical$her2_status_by_ihc[is.na(clinical$her2_status_by_ihc)] <- "Indeterminate"
  clinical$her2_status_by_ihc[clinical$her2_status_by_ihc == "Equivocal"] <- "Indeterminate"
  # Epithelial area
  clinical$epithelial_area[is.na(clinical$epithelial_area)] <- "<25% (Low)"
  # Inflammation
  clinical$inflammation[is.na(clinical$inflammation)] <- "Absent"
  # LCIS
  clinical$LCIS[is.na(clinical$LCIS)] <- "Absent"
  # Apocrine features
  clinical$apocrine_features[clinical$apocrine_features %in% c("1-5% (Minimum)", "6-50% (Moderate)", ">50% (Marked)")] <- "Present"
  clinical$apocrine_features[is.na(clinical$apocrine_features)] <- "Absent"
  # DCIS
  clinical$DCIS[is.na(clinical$DCIS)] <- "Absent"
  # tubule formation
  clinical$epithelial_tubule_formation[is.na(clinical$epithelial_tubule_formation)] <- "(score = 3) <10%"
  # LVI
  clinical$LVI[clinical$LVI == "Frequent"] <- "Present"
  clinical$LVI[clinical$LVI %in% c("No non-tumor tissue is present - cannot be evaluated.", 
                                   "No non-tumor tissue is present for evaluation.", 
                                   "Cannot be evaluated.")] <- "Absent"
  clinical$LVI[is.na(clinical$LVI)] <- "Absent"
  # necrosis
  clinical$necrosis[is.na(clinical$necrosis)] <- "Absent"
  # nuclear pleomorphism
  clinical$nuclear_pleomorphism[is.na(clinical$nuclear_pleomorphism)] <- "(score = 1) Small regular nuclei"
  # fibrous focus
  clinical$fibrous_focus[is.na(clinical$fibrous_focus)] <- "Absent"
  clinical$fibrous_focus[clinical$fibrous_focus == "Cannot be evaluated."] <- "Absent"
  # mitosis
  clinical$mitosis[is.na(clinical$mitosis)] <- "(score = 1) 0 to 5 per 10 HPF"
  
  # make clinical data numeric factors
  clinical_factor <- as.data.frame(lapply(clinical[,-1], factor))
  clinical_factor <- as.data.frame(lapply(clinical_factor, as.numeric))
  clinical_factor$TCGA_ID <- clinical$TCGA_ID
  
  return(clinical_factor)
  
}