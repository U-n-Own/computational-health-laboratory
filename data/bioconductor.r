# Installing packages for bioconductor, required for GEOquery package

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.18")

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))

# Install GEOquery package
BiocManager::install("GEOquery")

# Install affy package
BiocManager::install("affy", force = TRUE)