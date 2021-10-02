library(dplyr)
library(Seurat)
library(patchwork)
library(gdata)
library(ggplot2)
library(ggridges)
library(tidyr)
meta.data <- read.table("SDY997_EXP15176_celseq_meta.tsv.725704", header = TRUE)
meta.data %>% View()
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

# Subset cluster1,5 and sub-cluster
studied_cells.c1_5 <- SubsetData(studied_cells, subset.name = "RNA_snn_res.0.2", accept.value = c(1,5))
studied_cells.c1_5$RNA_snn_res.1.4 <- NULL
studied_cells.c1_5$RNA_snn_res.1 <- NULL
studied_cells.c1_5$RNA_snn_res.0.8 <- NULL
studied_cells.c1_5$RNA_snn_res.0.6 <- NULL
studied_cells.c1_5$RNA_snn_res.0.4 <- NULL
studied_cells.c1_5$RNA_snn_res.0.1 <- NULL

studied_cells.c1_5@meta.data %>% View() #only cluster 1 & 5

# Find variable feature
studied_cells.c1_5 <- FindVariableFeatures(studied_cells.c1_5, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(studied_cells.c1_5), 10)
plot1 <- VariableFeaturePlot(studied_cells.c1_5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale data
all.genes <- rownames(studied_cells.c1_5)
studied_cells.c1_5 <- ScaleData(studied_cells.c1_5, features = all.genes)

# linear Dimreduc and PCA plot
studied_cells.c1_5 <- RunPCA(studied_cells.c1_5, features = VariableFeatures(object = studied_cells.c1_5))
DimPlot(studied_cells.c1_5, reduction = "pca")

#using meta.data as ref for PC and apply matched pc1 and pc2 and type
meta.data <- read.table("SDY997_EXP15176_celseq_meta.tsv.725704", header = TRUE)
cell_name_meta <- meta.data[["cell_name"]]
cell_name_get <- colnames(studied_cells.c1_5)
j <- match(cell_name_get, cell_name_meta)
matched_type_cluster_1_5 <- meta.data$type[j]

pca_embeddings_mat <-  studied_cells[["pca"]]@cell.embeddings
pc1_vec <- pca_embeddings_mat[,1]
pc2_vec <- pca_embeddings_mat[,2]
df_with_matched_type_cluster_1_5 <- data.frame(pc1 = pc1_vec[j], pc2 = pc2_vec[j], type = as.factor(matched_type_cluster_1_5))

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 'pc1_vs_pc2_cluster_1_5_type_color.png')
fig_file_type
png(fig_file_type)

plot_match_cell_type <- ggplot(df_with_matched_type_cluster_1_5, aes(x = pc1, y = pc2, color = type)) + 
  geom_point(size = 1) + 
  labs(x = 'pc1', y = 'pc2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(plot_match_cell_type)
dev.off()

# determine Dim
ElbowPlot(studied_cells.c1_5, ndims=40)

# Cluster the cells
studied_cells.c1_5 <- FindNeighbors(studied_cells.c1_5, dims = 1:30, k.param = 5)
studied_cells.c1_5 <- FindClusters(studied_cells.c1_5, resolution = 0.3)
head(Idents(studied_cells.c1_5), 5)

# only umap on cluster 1,5
studied_cells.c1_5 <- RunUMAP(studied_cells.c1_5, dim = 1:30)
studied_cells.c1_5 <- RunTSNE(studied_cells.c1_5, dims = 1:30)
DimPlot(studied_cells.c1_5, label = TRUE, reduction = "umap")
DimPlot(studied_cells.c1_5, label = TRUE, reduction = "tsne")
re_0.3_cluster_cell_type <- DimPlot(studied_cells, label = TRUE,reduction = "tsne")
DoHeatmap(studied_cells.c1_5,features = VariableFeatures(object = studied_cells.c1_5))

#save Seurat obj after tSNE
saveRDS(studied_cells.c1_5, file = "subcluster_of_myeloid_cells.rds")

#read Seurat RDS
studied_cells.c1_5 <- readRDS(file = "subcluster_of_myeloid_cells.rds")

#Label only cluster 1,5
studied_cells.c1_5$sub_cluster1_5 <- as.character(Idents(studied_cells.c1_5))
studied_cells.c1_5@meta.data %>% View()
sort(studied_cells.c1_5@active.ident) %>% View()
Idents(studied_cells.c1_5) %>% View()
write.table(sort(studied_cells.c1_5@active.ident), file = 'subclustering of 1&5.txt')

rownames(studied_cells.c1_5@meta.data) %>% View()
as.character(studied_cells.c1_5@active.ident) == Idents(studied_cells.c1_5)

write.table(studied_cells.c1_5@active.ident, file = 'res_0.3_cluster_1_5_cell_type_per_cell.txt')
cluster_per_cell <- read.table('cluster_per_cell.736296.txt', header = TRUE)
cluster_per_cell %>% View()
write.table(cluster_per_cell, file = 'cell_type_per_cell_published.txt')

#remove urine cells using regex pattern match 
#regmatches(colnames(studied_cells.c1_5),gregexpr("U[0-9]{3}.{9}", colnames(studied_cells.c1_5)))
#cell_name_no_urine_cell <- gsub("U[0-9]{3}.{9}", '', colnames(studied_cells.c1_5))
cell_name_no_urine_cell %>% View()
cluster_per_cell$cell_name %>% View()

# matched cell name to create df with my cluster
i <- match(colnames(studied_cells.c1_5), cluster_per_cell$cell_name)
df <- data.frame(my = colnames(studied_cells.c1_5), pub = cluster_per_cell$cell_name[i], 
                 pub_cluster = cluster_per_cell$cluster[i], my_cluster =studied_cells.c1_5$sub_cluster1_5)
df %>% View()

df_no_urine_cells <- df[-grep("U[0-9]{3}.{9}", df$my), ] 
df_no_urine_cells %>% View() # 485 rows
write.table(df_no_urine_cells, file = 'cell_name_vs_cluster_no_urine_cells_df.txt')

# remane Idents
studied_cells.c1_5 <- RenameIdents(object = studied_cells.c1_5,
                                   '0' = 'CM0',
                                   '1' = 'CM3',
                                   '2' = 'CM1',
                                   '3' = 'CM4 + CM1',
                                   '4' = 'CM0',
                                   '5' = 'CM2',
                                   '6' = 'CM4',
                                   '7' = 'CB2b',
                                   '8' = 'CM1 + CM2',
                                   '9' = 'CM3',
                                   '10' = 'CM1')
cluster_1_5_labeled <- DimPlot(studied_cells.c1_5, label = TRUE, reduction = "tsne")
DimPlot(studied_cells.c1_5, label = TRUE, reduction = "tsne")

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 'cluster_1_5_labeled.png')
fig_file_type
png(fig_file_type)

print(cluster_1_5_labeled)
dev.off()


#create a reference tsne plot with published cell type vs tsne coord #485 rows
tsne_embeddings_mat <- studied_cells.c1_5[['tsne']]@cell.embeddings
studied_cells.c1_5[['tsne']]@cell.embeddings %>% View()
x <- match(rownames(studied_cells.c1_5[['tsne']]@cell.embeddings), cluster_per_cell$cell_name)
tsne_1_vec <- tsne_embeddings_mat[,1]
tsne_1_vec_cell_name <- data.frame(tsne_1 = tsne_1_vec, my =rownames(studied_cells.c1_5[['tsne']]@cell.embeddings), pub = cluster_per_cell$cell_name[x], cell_type = as.factor(cluster_per_cell$cluster[x]))
tsne_1_vec_cell_name %>% View()
tsne_2_vec <- tsne_embeddings_mat[,2]
cell_type = as.factor(cluster_per_cell$cluster[x])
rownames(studied_cells.c1_5[['tsne']]@cell.embeddings) == colnames(studied_cells.c1_5)
df_tsne_subcluster_myeloid <- data.frame(tsne_1 = tsne_1_vec, tsne_2 = tsne_2_vec, my =rownames(studied_cells.c1_5[['tsne']]@cell.embeddings), pub = cluster_per_cell$cell_name[x], cell_type = as.factor(cluster_per_cell$cluster[x]))
df_tsne_subcluster_myeloid <- drop_na(df_tsne_subcluster_myeloid)
df_tsne_subcluster_myeloid %>% View()
write.table(df_tsne_subcluster_myeloid, file = 'df_tsne_subcluster_myeloid.txt')

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 'tsne_subcluster_myeloid.png')
fig_file_type
png(fig_file_type)

plot_tsne_subcluster_myeloid <- ggplot(df_tsne_subcluster_myeloid, aes(x = tsne_1, y = tsne_2, color = cell_type)) + 
  geom_point(size = 1) + 
  labs(x = 'tsne1', y = 'tsne2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(plot_tsne_subcluster_myeloid)
dev.off()

##Valani gene_vs_cluster_mat.var_genes
library('sos')

valani_mat <- read.table(file.path("valani_paper_gene_vs_cluster_mat.var_genes.txt"), header = TRUE)
valani_mat %>% View()
head(valani_mat, 5)

my_gene = rownames(studied_cells.c1_5$RNA@data)
my_gene
villani_gene = rownames(valani_mat)
villani_gene

common_genes <- intersect(my_gene, villani_gene)
common_genes
glimpse(common_genes)

studied_cells.c1_5$RNA@data[common_genes, ] %>% View()

length(colnames(studied_cells.c1_5$RNA@data))

# calculate cor btwn vallani cluster and my cluster and assign max cor of vallani cluster to my cluster
cluster_with_max_corr <- c()
for(cell_idx in 1:length(colnames(studied_cells.c1_5$RNA@data))){
  my_gene_vec <- studied_cells.c1_5$RNA@data[common_genes, cell_idx]
  corr_per_cluster <- c()
for(cluster_idx in 1:length(colnames(valani_mat))){
  valani_vec <- valani_mat[common_genes, cluster_idx]
  corr_per_cluster[cluster_idx] <- cor(my_gene_vec, valani_vec)
  
} 
  sort_results <- sort.int(corr_per_cluster, decreasing = TRUE, index.return = TRUE)
cluster_with_max_corr[cell_idx] <- sort_results$ix[1]
}
cluster_with_max_corr 

# add new assigned vallani cluster to metadata
studied_cells.c1_5[['vallani_cluster']] <-  cluster_with_max_corr 

#plot tSNE
png("tsne_subcluster_colored_by_vallani_cluster_myeloid.png")
TSNEPlot(studied_cells.c1_5, group.by='vallani_cluster')
dev.off()

studied_cells.c1_5@meta.data


?sort.int

plot_tsne_subcluster_myeloid <- ggplot(df_tsne_subcluster_myeloid, aes(x = tsne_1, y = tsne_2, color = cluster_with_max_corr)) + 
  geom_point(size = 1) + 
  labs(x = 'tsne1', y = 'tsne2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(plot_tsne_subcluster_myeloid)
dev.off()

rownames(studied_cells.c1_5[['RNA']]@scale.data)
  
  
  
  
  
  