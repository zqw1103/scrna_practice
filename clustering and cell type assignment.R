library(dplyr)
library(Seurat)
library(patchwork)
library(gdata)
library(ggplot2)
library(ggridges)
meta.data <- read.table("SDY997_EXP15176_celseq_meta.tsv.725704", header = TRUE)
gene_by_cell_mat<- read.table("SDY997_EXP15176_celseq_matrix_ru10_molecules.tsv.725699.gz", header= TRUE, sep="\t", row.names = 1, as.is = TRUE)
gene_by_cell_mat[is.na(gene_by_cell_mat)] <- 0
studied_cells <- CreateSeuratObject(counts = gene_by_cell_mat, project = "gene_by_cell_mat", min.cells = 3, min.features = 100)
studied_cells[["percent.mt"]] <- PercentageFeatureSet(studied_cells, pattern = "^MT-")
studied_cells <- subset(studied_cells, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 25)
studied_cells <- NormalizeData(studied_cells)

#find variable genes
studied_cells <- FindVariableFeatures(studied_cells, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(studied_cells)
all_cells_names <- colnames(studied_cells)
studied_cells <- ScaleData(studied_cells, features = all.genes)
studied_cells <- RunPCA(studied_cells, features = VariableFeatures(object = studied_cells))
ElbowPlot(studied_cells)

#clustering
studied_cells <- FindNeighbors(studied_cells, dims = 1:20)
studied_cells <- FindClusters(studied_cells, resolution = c(0.1,0.2,0.4, 0.6, 0.8, 1.0, 1.4)) # explore diff resolution
studied_cells@meta.data %>% View() #look at obj@meta.data with added column of clusters after set reslution
head(Idents(studied_cells), 5)
# Idents(studied_cells) <- "RNA_snn_res.0.8" #assign res=0.8 by assigning the identity of the clusters using the Idents() function.
#head(Idents(studied_cells), 5) 

Idents(studied_cells) <- "RNA_snn_res.0.2" #use res=0.2 for clustering and later analysis
studied_cells <- RunUMAP(studied_cells, dims = 1:20)
studied_cells <- RunTSNE(studied_cells, dims = 1:20)
DimPlot(studied_cells, reduction = "umap")
DimPlot(studied_cells, label = TRUE, reduction = "umap")
DimPlot(studied_cells, label = TRUE, reduction = "tsne")
re_0.2_cluster_cell_type <- DimPlot(studied_cells, label = TRUE,reduction = "tsne")

# print png to file
output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 're_0.2_cluster_cell_type.png') # create export of png to local
fig_file_type
png(fig_file_type)
print(re_0.2_cluster_cell_type)
dev.off()


studied_cells.markers <- FindAllMarkers(studied_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(studied_cells.markers, 8) %>% View()
studied_cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top10_tab <- studied_cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
Cells(studied_cells.markers)

#create a table combining cell name from rownames(studied_cells.c1_5@meta.data) and cluster ident from as.character(Idents(studied_cells))
df_cluster_per_cell <- data.frame(cell_name = rownames(studied_cells@meta.data), cluster = as.character(Idents(studied_cells)))
df_cluster_per_cell %>% View()
write.table(df_cluster_per_cell, file = 'clustering per cell.txt')

# access cluster per cell data after clustering and match my cluster per cell and published cluster per cell
studied_cells@active.ident %>% View()
Idents(studied_cells) %>% View()
write.table(studied_cells@active.ident, file = 'res_0.2_cell_type_per_cell.txt')
cluster_per_cell <- read.table('cluster_per_cell.736296.txt', header = TRUE)
cluster_per_cell %>% View()
write.table(cluster_per_cell, file = 'cell_type_per_cell_published.txt')
i <- match(rownames(studied_cells@meta.data), cluster_per_cell$cell_name)

df_cluster_per_cell$cell_name %>% View()
cluster_per_cell$cell_name %>% View()
i %>% View()
write.table(i, file = 'matched_cell.txt')
matched_cluster_removed_NA = cluster_per_cell$cluster[na.omit(i)]    # why i have cells that were not assigned clusters in your file
matched_cluster_removed_NA %>% View()
matched_df_cluster_per_cell <- data.frame(my_cell_name = rownames(studied_cells@meta.data)[na.omit(i)], pub_cluster = matched_cluster_removed_NA )
matched_df_cluster_per_cell %>% View()
write.table(matched_df_cluster_per_cell, file = 'matched_my_cell_vs_published_cluster.txt')
# included NA in matched cell name to create df with my cluster
matched_cluster_all = cluster_per_cell$cluster[i]
matched_df_cluster_per_cell_all <- data.frame(my_cell_name = rownames(studied_cells@meta.data)[i], pub_cluster = matched_cluster_all, my_cluster = as.character(Idents(studied_cells)) )
matched_df_cluster_per_cell_all %>% View()
write.table(matched_df_cluster_per_cell_all, file = 'matched_my_cell_vs_published_cluster_vs_my_cluster.txt')


top10_tab %>% View()
studied_cells@meta.data %>% View()
write.table(top10_tab, file = 'res_0.2_cell_type_top_10.csv')
DoHeatmap(studied_cells, features = top10_tab$gene) 

top20_tab <- studied_cells.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC) 
top20_tab %>% View()
write.table(top20_tab, file = 'res_0.2_cell_type_top_20.csv')

# install package for cell type assignment
install.packages("clustermole")
library(clustermole)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")
install.packages("GSVA")
install.packages("htmltools")

cluster_cell_type.func <- function(i) {
  cluster_top_10_gene_tab <- studied_cells.markers %>% group_by(cluster) %>% filter(cluster == i) %>% top_n(n = 10, wt = avg_logFC) %>% select(gene)
  my_overlaps = clustermole_overlaps(genes = pull(cluster_top_10_gene_tab[, 2 ], gene), species = "hs")
  write.table(as.data.frame(my_overlaps), file = glue("res_0.2_cluster_cell_type {i}.csv"))
  
}


cluster_cell_type.func(10)

# scCATCH cluster cell type assignment
install.packages(pkgs = 'devtools')
devtools::install_github('ZJUFanLab/scCATCH')
library(scCATCH)
clu_markers <- findmarkergenes(object =  studied_cells, species = 'Human', cluster = 'All', match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
clu_markers$clu_markers %>% View()
clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Human',
                   cancer = NULL,
                   tissue = 'Kidney')

library(glue)
clu_ann[ , ] %>% View()
write.table(clu_ann[ , ], file = "re_0.2_cluster_cell_type.csv")















