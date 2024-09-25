options("repos" = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }

# BiocManager::install("limma", force = TRUE)

args <- commandArgs(trailingOnly = TRUE)

param1 <- args[1]
param2 <- args[2]
param3 <- args[3]
param4 <- args[4]

print(paste("Expression matrix file:", param1))
print(paste("Group info file:", param2))
print(paste("Output file:", param3))
print(paste("Compare group:", param4))


library(limma)

expr_matrix <- read.table(param1, header = TRUE, sep = "\t", row.names = 1)
sample_info <- read.table(param2, header = TRUE, sep = "\t", row.names = 1)

group <- factor(sample_info$group)

design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

contrast_matrix <- makeContrasts(contrasts = param4, levels = design)

fit <- lmFit(expr_matrix, design)

fit2 <- contrasts.fit(fit, contrast_matrix)

fit2 <- eBayes(fit2)

results <- topTable(fit2, adjust = "fdr", sort.by = "P", number = Inf)

write.csv(results, param3)