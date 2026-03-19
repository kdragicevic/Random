Dataset comparison top 20 marker genes example
================
Katarina Dragicevic

``` r
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# ----------------------------
# settings
# ----------------------------
set.seed(123)
n_cells_per_type <- 500

obj1 <- dataset1 #(seurat object)
obj2 <- dataset2 #(seurat object)

anno1 <- "anno1" #annotation metadata column name
anno2 <- "anno2" #annotation metadata column name

# ----------------------------
# sample up to 500 cells per cell type
# ----------------------------
sample_cells_per_type <- function(seu, anno_col, n = 500) {
  
  md <- seu@meta.data
  md$cell_id <- rownames(md)
  md$anno <- as.character(md[[anno_col]])
  
  md <- md[!is.na(md$anno), , drop = FALSE]
  
  sampled_cells <- md %>%
    dplyr::group_by(anno) %>%
    dplyr::group_modify(~{
      k <- min(nrow(.x), n)
      .x[sample(seq_len(nrow(.x)), k), ]
    }) %>%
    dplyr::ungroup() %>%
    dplyr::pull(cell_id)
  
  subset(seu, cells = sampled_cells)
}

obj1_sub <- sample_cells_per_type(obj1, anno1, n_cells_per_type)
obj2_sub <- sample_cells_per_type(obj2, anno2, n_cells_per_type)

# ----------------------------
# keep shared genes only
# ----------------------------
shared_genes <- intersect(rownames(obj1_sub), rownames(obj2_sub))
#DefaultAssay(obj1_sub) <- "RNA"
obj1_sub <- subset(obj1_sub, features = shared_genes)
obj2_sub <- subset(obj2_sub, features = shared_genes)

# ----------------------------
# marker in each dataset
# ----------------------------
DefaultAssay(obj1_sub) <- "RNA"
DefaultAssay(obj2_sub) <- "RNA"

obj1_sub <- NormalizeData(obj1_sub, normalization.method = "LogNormalize", scale.factor = 10000)
obj2_sub <- NormalizeData(obj2_sub, normalization.method = "LogNormalize", scale.factor = 10000)

Idents(obj1_sub) <- obj1_sub@meta.data[[anno1]]
Idents(obj2_sub) <- obj2_sub@meta.data[[anno2]]

markers1 <- FindAllMarkers(
  obj1_sub,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

markers2 <- FindAllMarkers(
  obj2_sub,
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)


top_n <- 20

top50_1 <- markers1 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  ungroup()

top50_2 <- markers2 %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  ungroup()

marker_genes <- union(top50_1$gene, top50_2$gene)
marker_genes <- intersect(marker_genes, intersect(rownames(obj1_sub), rownames(obj2_sub)))
```

``` r
avg3 <- AverageExpression(
  obj1_sub,
  assays = "RNA",
  features = marker_genes,
  group.by = anno1,
  slot = "data",
  verbose = FALSE
)$RNA

avg4 <- AverageExpression(
  obj2_sub,
  assays = "RNA",
  features = marker_genes,
  group.by = anno2,
  slot = "data",
  verbose = FALSE
)$RNA

avg3 <- as.matrix(avg3)
avg4 <- as.matrix(avg4)
```

``` r
common_genes <- intersect(rownames(avg3), rownames(avg4))
avg3 <- avg3[common_genes, , drop = FALSE]
avg4 <- avg4[common_genes, , drop = FALSE]

cor_mat2 <- cor(avg3, avg4, method = "pearson")

row_order <- apply(cor_mat2, 1, which.max)
cor_mat_ordered2 <- cor_mat2[order(row_order), ]
```

``` r
m <- scale(cor_mat_ordered2)
hm <- Heatmap(
  m,
  name = "Pearson_Markergene",
  cluster_rows = F,
  cluster_columns = F,
  row_title = "data1",
  column_title = "data2",
  column_names_rot = 90, col = colorRamp2(c(-2, 0, 2), c("steelblue4","white", "firebrick4"))
)
```
