library(Seurat)
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(scales)

# Load the CLL object. Carefull! It contains all the cells in the PBMC
load("Data/PBMC/scRNAseq_obj_CLL.RData")
table(cll_obj$cell_types)

# Extract the CLL cells (CLL + B cells)
cll_obj <- subset(cll_obj, subset = cell_types == c("CLL", "Bcells"))

cll_obj <- NormalizeData(cll_obj)
cll_obj <- FindVariableFeatures(cll_obj)
cll_obj <- ScaleData(cll_obj)
cll_obj <- RunPCA(cll_obj)
cll_obj <- FindNeighbors(cll_obj, dims = 1:20)
cll_obj <- FindClusters(cll_obj, resolution = 0.5)
cll_obj <- RunUMAP(cll_obj, dims = 1:20)
DimPlot(cll_obj, label = TRUE)

# Extract the gene expression matrix for the CLL compartment
norm_matrix <- GetAssayData(cll_obj, layer = "data")
norm_matrix <- as.matrix(norm_matrix)

# Remove non-coding genes from the matrix
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_list <- rownames(norm_matrix)

# Query Ensembl for gene biotype
annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_list,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GENETYPE")
)

# Filter to keep only protein-coding genes
coding_genes <- annotations$SYMBOL[annotations$GENETYPE == "protein-coding"]

CLL_normalised_coding <- norm_matrix[rownames(norm_matrix) %in% coding_genes, ]

# Remove the all-zero genes
CLL_normalised_coding_0Rem <- CLL_normalised_coding[rowSums(CLL_normalised_coding) > 0, ]

write.csv(CLL_normalised_coding_0Rem, "CLL_only/CLL_normalized_matrix_coding_0Rem.csv")


##########################################################################################
##### Manipulate the seurat object #######################################################

# Temporal mapping 
DimPlot(cll_obj, reduction = "umap", group.by = "time_point") +
  ggtitle("Time point mapping") + 
  scale_color_manual(values = c("M0"  = "#66C2A5",
                                "M3"  = "#FC8D62",
                                "M27" = "#8DA0CB"))

# Mapping by the expression of particular TFs/genes
p1 <- FeaturePlot(cll_obj, "TAF7") +  ggtitle("TAF7 expression")
p2 <- FeaturePlot(cll_obj, "RBM3") +  ggtitle("RBM3 expression")

p1 | p2

# Test the significance of correlation between TF/gene expression and time
expr_df <- FetchData(cll_obj, vars = c("TAF7", "time_point")) %>%
  group_by(time_point) %>%
  summarise(mean_expr = mean(TAF7))
anova(lm(TAF7 ~ time_point, data = FetchData(cll_obj, vars = c("TAF7", "time_point"))))

###########################################################################################

# Use FindMarker function to identify the genes are more expressed in one time-point compared to the other ones
Idents(cll_obj) <- "time_point"
markers_M0  <- FindMarkers(cll_obj, ident.1 = "M0")
markers_M3  <- FindMarkers(cll_obj, ident.1 = "M3")
markers_M27 <- FindMarkers(cll_obj, ident.1 = "M27")

sig_M0  <- markers_M0  %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames()
writeLines(sig_M0, "Results/CLL_only/genes_M0_corrected.txt")
sig_M3  <- markers_M3  %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames()
writeLines(sig_M3, "Results/CLL_only/genes_M3_corrected.txt")
sig_M27 <- markers_M27 %>% filter(avg_log2FC > 0 & p_val_adj < 0.05) %>% rownames()
writeLines(sig_M27, "Results/CLL_only/genes_M27_corrected.txt")

