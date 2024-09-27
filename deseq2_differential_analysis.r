options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("DESeq2", force = TRUE)

args <- commandArgs(trailingOnly = TRUE)

param1 <- args[1]
param2 <- args[2]
param3 <- args[3]

print(paste("Counts file:", param1)) #counts.txt
print(paste("Group info file:", param2)) #transcript_group_info.csv
print(paste("Output file:", param3))

library(DESeq2)

counts_file <- param1
group_info_file <- param2

count_data <- read.csv(param1, header = TRUE, sep = "\t", row.names = 1)
sample_info <- read.csv(param2, header = TRUE, sep = "\t", row.names = 1)

col_data <- data.frame(row.names = colnames(count_data), condition = sample_info$group)

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ condition)

dds <- DESeq(dds)

results <- results(dds)

write.csv(as.data.frame(results), param3)