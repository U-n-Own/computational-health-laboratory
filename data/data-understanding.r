# Install and load the 'affy' package for CEL file import
# Installing packages for bioconductor, required for GEOquery package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
# Install GEOquery package
BiocManager::install("GEOquery")
# Install affy package
BiocManager::install("affy", force = TRUE)
library(affy)
# Specify the directory path containing the CEL files
cel_directory <- "./dataset/GSE235356_RAW"
# Read the CEL files and preprocess the data using the robust multi-array average (RMA) method
cel_data <- ReadAffy(celfile.path = cel_directory)
gene_expression_data <- exprs(rma(cel_data))
# Perform a log2 transformation on the gene expression data
gene_expression_data <- log2(gene_expression_data)
# Specify the Probeset names of interest
probeset_names <- c("234764_x_at", "211835_at", "1561937_x_at", "202716_at", "235305_s_at", "210538_s_at", "237461_at", "41660_at", "217892_s_at", "225822_at", "57532_at", "213489_at", "222641_s_at", "205159_at", "209012_at", "220522_at", "223709_s_at", "36129_at", "201848_s_at", "212704_at", "213622_at", "232531_at", "205666_at", "210789_x_at", "217809_at", "225291_at", "226488_at", "231131_at", "238662_at", "226098_at", "202387_at", "228217_s_at", "225553_at", "223995_at", "202613_at", "203200_s_at")
# Extract the expression data for the specified Probesets
probeset_expression_data <- gene_expression_data[probeset_names, ]
# Calculate the average expression for Down-Regulated Probesets
down_regulated <- rowMeans(probeset_expression_data[c("234764_x_at", "211835_at", "1561937_x_at", "202716_at", "235305_s_at", "210538_s_at", "237461_at", "41660_at", "217892_s_at", "225822_at", "57532_at", "213489_at", "222641_s_at", "205159_at", "209012_at", "220522_at", "223709_s_at", "36129_at", "201848_s_at", "212704_at", "213622_at", "232531_at", "205666_at", "210789_x_at"),])
# Calculate the average expression for Up-Regulated Probesets
up_regulated <- rowMeans(probeset_expression_data[c("217809_at", "225291_at", "226488_at", "231131_at", "238662_at", "226098_at", "202387_at", "228217_s_at", "225553_at", "223995_at", "202613_at", "203200_s_at"),])
# Subtract Down-Regulated from Up-Regulated
result <- up_regulated - down_regulated
# Classify samples as high-risk or low-risk MGUS based on the result
mgus_classification <- ifelse(result > 0.7, "High Risk MGUS", "Low Risk MGUS")

#Export matrix as CSV

write.csv(gene_expression_data, file = "./dataset/gene_expression_data.csv.gz", row.names = TRUE)

