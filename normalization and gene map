
cat("\014")

library(affy)
library(limma)
library(simpleaffy)
library(MASS)
library(readxl)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(
  c("WGCNA", "GO.db", "impute", "AnnotationDbi"),
  ask = FALSE,
  update = FALSE
)

library(WGCNA)

library(WGCNA)


#############################
#Read Data
#############################
#BiocManager::install("hgu133plus2cdf")

setwd("C:\\Master 4\\thesis\\bioinformatics\\GEO\\GSE24982\\GSE24982_RAW")

dat<-ReadAffy()
dat


#############################
#Normalization
#############################
dat2<-rma(dat)
dat2<-justRMA()
dat2
dat.m<-exprs(dat2)

dat1<-exprs(dat)

############################
## GSE24982 - Gene Map Only
## Read CEL -> RMA -> ProbeID to GeneSymbol (GPL1355.xlsx) -> collapse -> export
############################

## 0) Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("affy", quietly = TRUE)) BiocManager::install("affy", ask = FALSE, update = FALSE)
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")

library(affy)
library(readxl)
library(WGCNA)
options(stringsAsFactors = FALSE)

############################
## 1) Paths (EDIT THESE)
############################
raw_dir  <- "C:/Master 4/thesis/bioinformatics/GEO/GSE24982/GSE24982_RAW"  # CEL files
out_dir  <- "C:/Master 4/thesis/bioinformatics/GEO/GSE24982"        # outputs
ann_file <- "C:/Master 4/thesis/bioinformatics/GEO/GSE24982/GPL1355.xlsx" # annotation: ID, Gene Symbol

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

############################
## 2) Read CEL files
############################
setwd(raw_dir)
dat <- ReadAffy()   # AffyBatch
sample_names <- sampleNames(dat)

############################
## 3) RMA normalization
############################
eset <- rma(dat)            # ExpressionSet
norm_expr <- exprs(eset)    # probes x samples
probe_ids <- rownames(norm_expr)

############################
## 4) Read & clean annotation (GPL1355.xlsx)
############################
ann <- as.data.frame(read_excel(ann_file))

# Ensure expected columns exist
stopifnot(all(c("ID", "Gene Symbol") %in% colnames(ann)))

# Clean gene symbols (remove empty / multi-annotation like "Abc1 /// Xyz2")
ann$`Gene Symbol` <- trimws(as.character(ann$`Gene Symbol`))
ann$`Gene Symbol` <- sub("\\s*///.*$", "", ann$`Gene Symbol`)
ann$`Gene Symbol` <- trimws(ann$`Gene Symbol`)

############################
## 5) Map probe -> gene symbol
############################
map <- ann[match(probe_ids, ann$ID), c("ID", "Gene Symbol"), drop = FALSE]
colnames(map) <- c("ProbeID", "GeneSymbol")

keep <- !is.na(map$GeneSymbol) & map$GeneSymbol != ""
expr_keep <- norm_expr[keep, , drop = FALSE]
genes <- map$GeneSymbol[keep]

############################
## 6) Collapse multiple probes per gene (Average)
############################
collapsed <- collapseRows(
  datET    = expr_keep,
  rowGroup = genes,
  rowID    = rownames(expr_keep),
  method   = "Average"
)

gene_expr <- collapsed$datETcollapsed   # gene x sample

############################
## 7) Build final tidy table and export
############################
final_gene_matrix <- data.frame(
  GeneSymbol = rownames(gene_expr),
  gene_expr,
  check.names = FALSE
)

# Optional: sort genes alphabetically
final_gene_matrix <- final_gene_matrix[order(final_gene_matrix$GeneSymbol), ]

write.csv(final_gene_matrix,
          file = file.path(out_dir, "GSE24982_RMA_GeneExpression.csv"),
          row.names = FALSE)

write.table(final_gene_matrix,
            file = file.path(out_dir, "GSE24982_RMA_GeneExpression.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

############################
## 8) Quick checks
############################
cat("\nSaved in:", out_dir, "\n")
cat("Samples:", ncol(norm_expr), "\n")
cat("Probes (after RMA):", nrow(norm_expr), "\n")
cat("Probes kept (mapped to genes):", nrow(expr_keep), "\n")
cat("Genes after collapse:", nrow(gene_expr), "\n")

