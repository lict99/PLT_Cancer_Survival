# env settings ----

library(magrittr)

dir.create("00", FALSE)

# predefined objects ----

## cancer ICD codes
cancer_ICD_codes <- list(
  All_sites = "^C[0-8][0-9]|^C9[0-7]|^1[4-9][0-9]|^20[0-8]",
  E_lymphoid_haematopoietic = "^C[0-7][0-9]|^C80|^1[4-9][0-9]",
  Lymphoid_haematopoietic = "^C8[1-9]|^C9[0-6]|^20[0-8]",
  Female_breast = "^C50|^174" %>%
    set_attr("sex", "female"),
  Lung = "^C3[34]|^162",
  Colorectum = "^C1[89]|^C2[01]|^15[34]",
  Prostate = "^C61|^185" %>%
    set_attr("sex", "male"),
  Nonmelanoma_of_skin = "^C44|^173",
  Stomach = "^C16|^151",
  Head_and_neck = "^C0[0-9]|^C1[0-4]|^C3[0-2]|^14[0-9]|^16[01]",
  Liver = "^C22|^155",
  Cervix_uteri = "^C53|^180" %>%
    set_attr("sex", "female"),
  Esophagus = "^C15|^150",
  Thyroid = "^C73|^193",
  Bladder = "^C67|^188",
  Non_Hodgkin_lymphoma = "^C8[2-6]|^C96|^20[02]",
  Pancreas = "^C25|^157",
  Leukemia = "^C9[1-5]|^20[4-8]",
  Kidney = "^C6[45]|^189",
  Corpus_uteri = "^C54|^182" %>%
    set_attr("sex", "female"),
  Melanoma_of_skin = "^C43|^172",
  Ovary = "^C56|^1830" %>%
    set_attr("sex", "female"),
  Brain_nervous_system = "^C7[0-2]|^19[12]",
  Multiple_myeloma = "^C88|^C90|^203",
  Gallbladder = "^C23|^1560",
  Hodgkin_lymphoma = "^C81|^201",
  # Anus = "", # merged to colorectum
  Mesothelioma = "^C45", # no ICD-9 codes exist
  Secondary = "^C7[7-9]|^19[6-8]"
)

## cancer sites
cancer_sites <- list(
  All_sites = "All sites",
  E_lymphoid_haematopoietic =
    "All sites excluding lymphoid, haematopoietic and related tissue",
  Lymphoid_haematopoietic = "Lymphoid, haematopoietic and related tissue",
  Female_breast = "Female breast",
  Lung = "Lung",
  Colorectum = "Colorectum",
  Prostate = "Prostate",
  Nonmelanoma_of_skin = "Nonmelanoma of skin",
  Stomach = "Stomach",
  Head_and_neck = "Head and neck",
  Liver = "Liver",
  Cervix_uteri = "Cervix uteri",
  Esophagus = "Esophagus",
  Thyroid = "Thyroid",
  Bladder = "Bladder",
  Non_Hodgkin_lymphoma = "Non-Hodgkin lymphoma",
  Pancreas = "Pancreas",
  Leukemia = "Leukemia",
  Kidney = "Kidney",
  Corpus_uteri = "Corpus uteri",
  Melanoma_of_skin = "Melanoma of skin",
  Ovary = "Ovary",
  Brain_nervous_system = "Brain, nervous system",
  Multiple_myeloma = "Multiple myeloma",
  Gallbladder = "Gallbladder",
  Hodgkin_lymphoma = "Hodgkin lymphoma",
  Mesothelioma = "Mesothelioma",
  Secondary = "Secondary malignant neoplasm"
)

## cancer names (cancer types)
cancer_names <- list(
  All_sites = "Pan-cancer",
  E_lymphoid_haematopoietic = "Solid cancer",
  Lymphoid_haematopoietic = "Hematologic/lymphatic cancer",
  Female_breast = "Breast cancer",
  Lung = "Lung cancer",
  Colorectum = "Colorectal cancer",
  Prostate = "Prostate cancer",
  Nonmelanoma_of_skin = "Nonmelanoma skin cancer",
  Stomach = "Gastric cancer",
  Head_and_neck = "Head and neck cancer",
  Liver = "Liver cancer",
  Cervix_uteri = "Cervical cancer",
  Esophagus = "Esophageal cancer",
  Thyroid = "Thyroid cancer",
  Bladder = "Bladder cancer",
  Non_Hodgkin_lymphoma = "Non-Hodgkin lymphoma",
  Pancreas = "Pancreatic cancer",
  Leukemia = "Leukemia",
  Kidney = "Kidney cancer",
  Corpus_uteri = "Endometrial cancer",
  Melanoma_of_skin = "Melanoma",
  Ovary = "Ovarian cancer",
  Brain_nervous_system = "Central nervous system cancer",
  Multiple_myeloma = "Multiple myeloma",
  Gallbladder = "Gallbladder cancer",
  Hodgkin_lymphoma = "Hodgkin lymphoma",
  Mesothelioma = "Mesothelioma",
  Secondary = "Secondary cancer"
)

# data saving ----

if (
  identical(names(cancer_ICD_codes), names(cancer_names)) &&
    identical(names(cancer_ICD_codes), names(cancer_sites))
) {
  save(cancer_names, file = "00/cancer_names.RData")
  save(cancer_sites, file = "00/cancer_sites.RData")
  save(cancer_ICD_codes, file = "00/cancer_ICD_codes_with_attr.RData")
} else {
  stop("Different cancer definitions!")
}
