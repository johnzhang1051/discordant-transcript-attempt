# Installing depmap

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("depmap")

library(depmap)

expression_data <- depmap_rnai()  # Gene dependency scores

expression_data <- depmap_tpm()   # Expression data if available

# Get cell line info
cell_info <- depmap_metadata()


# Get cells with melanoma
melanoma_cells <- filter(cell_info, grepl("Melanoma", subtype_disease, ignore.case = TRUE))

#View(cell_info[cell_info$lineage == "skin",])

