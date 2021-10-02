meta.data <- read.table("SDY997_EXP15176_celseq_meta.tsv.725704")
head(meta.data, 10)
tail(meta.data, 10)
colnames(meta.data)
?Read10X
 
library(dplyr)
library(Seurat)
library(patchwork)
library(gdata)
library(ggplot2)
library(ggridges)
library(tidyr)
library(tibble)
gene_by_cell_mat<- read.table("SDY997_EXP15176_celseq_matrix_ru10_molecules.tsv.725699.gz", header= TRUE, sep="\t", row.names = 1, as.is = TRUE)
gene_by_cell_mat[is.na(gene_by_cell_mat)] <- 0
studied_cells <- CreateSeuratObject(counts = gene_by_cell_mat, project = "gene_by_cell_mat", min.cells = 3, min.features = 100)
gene_by_cell_mat
studied_cells
gene_by_cell_mat.data[c("UBE2J2", "VWA1", "CDK11A"), 1:20]
studied_cells[["percent.mt"]] <- PercentageFeatureSet(studied_cells, pattern = "^MT-")
head(studied_cells@meta.data,5)
?VlnPlot
VlnPlot(studied_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.25)
?FeatureScatter
plot1 <- FeatureScatter(studied_cells, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = rep("red", 24), pt.size = 0.1)
plot2 <- FeatureScatter(studied_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = rep("red", 24), pt.size = 0.1)
plot1
plot2

#add a column in Seurat obect metadata(object@meta.data)   #following are same
studied_cells$log10GenesPerUMI <- log10(studied_cells$nFeature_RNA) / log10(studied_cells$nCount_RNA)  
studied_cells[["log10GenesPerUMI"]] <- log10(studied_cells$nFeature_RNA) / log10(studied_cells$nCount_RNA)
head(studied_cells@meta.data,5)

#studied_cells@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

#studied_cells@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, fill = "red")) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#Assays(studied_cells)
#GetAssay(object = studied_cells, assay = "RNA")
#head(x = rownames(x = studied_cells[["RNA"]]))
#head(x = colnames(x = studied_cells[["RNA"]]))
#GetAssayData(object = studied_cells, assay = "RNA", slot = "data")

#counts_per_cell <- Matrix::colSums(studied_cells)
#hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
#counts_per_gene <- Matrix::rowSums(studied_cells)
#hist(log10(counts_per_gene+1), main='counts per gene', col='wheat')
#genes_per_cell <- Matrix::colSums(gene_by_cell_mat>0)
#hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')

# filter high QC
studied_cells <- subset(studied_cells, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 25)

#normalize data
studied_cells <- NormalizeData(studied_cells)


#find variable genes
studied_cells <- FindVariableFeatures(studied_cells, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(studied_cells), 10)
plot1 <- VariableFeaturePlot(studied_cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

all.genes <- rownames(studied_cells)
all_cells_names <- colnames(studied_cells)
studied_cells <- ScaleData(studied_cells, features = all.genes)


studied_cells <- RunPCA(studied_cells, features = VariableFeatures(object = studied_cells))
studied_cells[["pca"]]
print(studied_cells[["pca"]], dims = 1:5, nfeatures = 15)
VizDimLoadings(studied_cells, dims = 1:5, reduction = "pca", nfeatures = 15)
DimPlot(studied_cells, reduction = "pca")

#studied_cells <- SCTransform(studied_cells, vars.to.regress = "percent.mt", verbose = FALSE)
#studied_cells <- RunPCA(studied_cells, verbose = FALSE)
#DimPlot(studied_cells, label = TRUE)

? as.data.frame
studied_cells@assays
DimHeatmap(studied_cells, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(studied_cells)

#Elbow plot: quantitative approach

# Determine percent of variation associated with each PC
pct <- studied_cells[["pca"]]@stdev / sum(studied_cells[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

#clustering steps
studied_cells <- FindNeighbors(studied_cells, dims = 1:20)
studied_cells <- FindClusters(studied_cells, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4)) # explore diff resolution
studied_cells@meta.data %>% View() #look at obj@meta.data with added column of clusters after set reslution
head(Idents(studied_cells), 5)
Idents(studied_cells) <- "RNA_snn_res.0.6" #assign res=0.6 by assigning the identity of the clusters using the Idents() function.
head(Idents(studied_cells), 5) 

studied_cells <- RunUMAP(studied_cells, dims = 1:20)
DimPlot(studied_cells, reduction = "umap") + label()
Idents(studied_cells) <- "RNA_snn_res.0.8"
DimPlot(studied_cells, label = TRUE, reduction = "umap")

Idents(studied_cells) <- "RNA_snn_res.0.6"
studied_cells <- RunTSNE(studied_cells, dims = 1:20)
DimPlot(studied_cells, reduction = "tsne")
Idents(studied_cells) <- "RNA_snn_res.0.8"
DimPlot(studied_cells, reduction = "tsne")

saveRDS(studied_cells, file = "AMP_tutorial.rds")


#DE find cluster biomarkers
cluster1.markers <- FindMarkers(studied_cells, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster2.markers <- FindMarkers(studied_cells, ident.1 = 2, ident.2 = c(3, 4), min.pct = 0.25)
rownames(head(cluster2.markers, n = 10))
rownames(cluster2.markers)

studied_cells.markers <- FindAllMarkers(studied_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(studied_cells.markers, 8) %>% View()
studied_cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


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


# add gene annotation file. if 'gene' is rownames, using rownames_to_column(var='gene)
annotations <- read.csv("annotation.csv")
cluster_ann_markers <- studied_cells.markers %>% left_join(y = unique(annotations[, c("gene_name", "description")]), by = c("gene" = "gene_name"))
cluster_ann_markers %>% View()

#assigning cell types to clusters using "clustermole"
cluster0_top_10_gene_tab <- studied_cells.markers %>% group_by(cluster) %>% filter(cluster == "0") %>% top_n(n = 10, wt = avg_logFC) %>% select(gene) # with conditions to create a table of 2 colums(cluster and gene)
cluster0_top_10_gene <- cluster0_top_10_gene_tab[, 2 ]
pull(cluster0_top_10_gene, gene) # extract values of gene column
my_overlaps = clustermole_overlaps(genes = pull(cluster0_top_10_gene, gene), species = "hs") %>% View()
my_overlaps

#visualize cluster 0
VlnPlot(studied_cells, features = pull(cluster0_top_10_gene, gene))
FeaturePlot(studied_cells, features = pull(cluster0_top_10_gene, gene))

top10_tab <- studied_cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) 
DoHeatmap(studied_cells, features = top10_tab$gene) 

# assign cluster_cell_type function by input of cluster number and write to file
cluster_cell_type.func <- function(i) {
  library(glue)
  cluster_top_10_gene_tab <- studied_cells.markers %>% group_by(cluster) %>% filter(cluster == i) %>% top_n(n = 10, wt = avg_logFC) %>% select(gene)
    my_overlaps = clustermole_overlaps(genes = pull(cluster_top_10_gene_tab[, 2 ], gene), species = "hs")
    write.table(as.data.frame(my_overlaps), file = glue("cluster_cell_type {i}.csv"))

}


cluster_cell_type.func(18)


cluster_top_10_gene_tab <- studied_cells.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) %>% select(gene)
cluster_top_10_gene_tab %>% View()
cluster_top_10_gene <- cluster_top_10_gene_tab[, 2 ]
cluster_top_10_gene %>% View()
pull(cluster_top_10_gene, gene)
my_overlaps = clustermole_overlaps(genes = pull(cluster_top_10_gene, gene), species = "hs") %>% View()


#pct <- studied_cells[["pca"]]@stdev / sum(studied_cells[["pca"]]@stdev) * 100
#pct

#colnames(studied_cells[["pca"]]@cell.embeddings)

#as_tibble(
  #studied_cells[[c("nCount_RNA","nFeature_RNA","percent.mt")]],
  #rownames="Cell.Barcode"
#) -> qc.metrics

#qc.metrics

#studied_cells[["pca"]]@cell.embeddings[,] %>%
  #as_tibble(rownames = "Cell.Barcode") -> barcode
#barcode

#studied_cells[["pca"]]@feature.loadings[,] %>%
  #as_tibble(rownames = "Feature.name") -> feature
#feature

#studied_cells@reductions[['pca']]
#studied_cells[["pca"]]

#generate a figure saved to local after ggplot of pca colored by batch
pca_embeddings_mat <-  studied_cells[["pca"]]@cell.embeddings
pca_embeddings_mat
pc1_vec <- pca_embeddings_mat[,1]
pc1_vec
pc2_vec <- pca_embeddings_mat[,2]
colnames(studied_cells@meta.data)
batch_per_cell <- studied_cells@meta.data$orig.ident
batch_per_cell
df_to_show <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, batch = as.factor(batch_per_cell))
df_to_show %>% View()

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file <- paste0(output_dir, 'pc1_vs_pc2.png')
fig_file
png(fig_file)
p <- ggplot(df_to_show, aes(x = pc1, y = pc2, color = batch)) + 
  geom_point(size = 1) + 
  labs(x = 'pc1', y = 'pc2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(p)
dev.off()

?SubsetData

studied_cells[['RNA']]@data

features <- rownames(studied_cells[["pca"]]@feature.loadings)
features




# match cell name btw metadata and Seurat pca data and plot with colored by condintion 
cell_name_metadata <- meta.data[["cell_name"]]
cell_name_to_get <- colnames(studied_cells)

i <- match(cell_name_to_get, cell_name_metadata)
matched_condition <- meta.data$disease[i]

pca_embeddings_mat <-  studied_cells[["pca"]]@cell.embeddings
pca_embeddings_mat
pc1_vec <- pca_embeddings_mat[,1]
length(pc1_vec)
pc2_vec <- pca_embeddings_mat[,2]
df_with_matched_cell_name <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, condition = as.factor(matched_condition))

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_condition <- paste0(output_dir, 'pc1_vs_pc2_condition_color.png')
fig_file_condition
png(fig_file_condition)

plot_match_cell_name <- ggplot(df_with_matched_cell_name, aes(x = pc1, y = pc2, color = condition)) + 
  geom_point(size = 1) + 
  labs(x = 'pc1', y = 'pc2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(plot_match_cell_name)
dev.off()

#match cell name btw metadata and Seurat pca data and plot with colored by type
cell_name_metadata <- meta.data[["cell_name"]] #access cell name from metadata
cell_name_to_get <- colnames(studied_cells)

# cell_name_to_get_bool <- cell_name_metadata %in% cell_name_to_get # no need in correct code. use match()
# match_cell_name <- cell_name_metadata[cell_name_to_get_bool] # no need in correct code. use match()

#correct here as instructed
i <- match(cell_name_to_get, cell_name_metadata) #use match() to get exact order of cell type as cell name by position
matched_type <- meta.data$type[i]
matched_type

pca_embeddings_mat <-  studied_cells[["pca"]]@cell.embeddings # access pcs cell name
pc1_vec <- pca_embeddings_mat[,1] # access pc1 coord based on cell name
pc2_vec <- pca_embeddings_mat[,2] # access pc2 coord based on cell name
df_with_matched_type <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, type = as.factor(matched_type)) # df with factored type

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 'pc1_vs_pc2_type_color.png') # create export of png to local
fig_file_type
png(fig_file_type)

plot_match_cell_type <- ggplot(df_with_matched_type, aes(x = pc1, y = pc2, color = type)) + 
  geom_point(size = 1) + 
  labs(x = 'pc1', y = 'pc2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(plot_match_cell_type)
dev.off()

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
studied_cells.c1_5 <- FindClusters(studied_cells.c1_5, resolution = 0.1)
head(Idents(studied_cells.c1_5), 5)

# only umap on cluster 1,5
studied_cells.c1_5 <- RunUMAP(studied_cells.c1_5, dim = 1:30)
studied_cells.c1_5 <- RunTSNE(studied_cells.c1_5, dims = 1:30)
DimPlot(studied_cells.c1_5, label = TRUE, reduction = "umap")
DimPlot(studied_cells.c1_5, label = TRUE, reduction = "tsne")
re_0.1_cluster_cell_type <- DimPlot(studied_cells, label = TRUE,reduction = "tsne")
DoHeatmap(studied_cells.c1_5,features = VariableFeatures(object = studied_cells.c1_5))

#Label only cluster 1,5
studied_cells.c1_5$sub_cluster1_5 <- as.character(Idents(studied_cells.c1_5))
studied_cells.c1_5@meta.data %>% View()
sort(studied_cells.c1_5@active.ident) %>% View()
Idents(studied_cells.c1_5) %>% View()
write.table(sort(studied_cells.c1_5@active.ident), file = 'subclustering of 1&5.txt')

rownames(studied_cells.c1_5@meta.data) %>% View()
as.character(studied_cells.c1_5@active.ident) == Idents(studied_cells.c1_5)

write.table(studied_cells.c1_5@active.ident, file = 'res_0.1_cluster_1_5_cell_type_per_cell.txt')
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
                                   '0' = 'most CM0 + most CM1 + CM4 + CT4',
                                   '1' = 'CM1',
                                   '2' = 'most CM2 + most CM3 + CM4',
                                   '3' = 'CB2b',
                                   '4' = 'CM3')
cluster_1_5_labeled <- DimPlot(studied_cells.c1_5, label = TRUE, reduction = "tsne")
DimPlot(studied_cells.c1_5, label = TRUE, reduction = "tsne")

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 'cluster_1_5_labeled.png')
fig_file_type
png(fig_file_type)

print(cluster_1_5_labeled)
dev.off()


# headder = True to read.table()
meta.data <- read.table("SDY997_EXP15176_celseq_meta.tsv.725704", header = TRUE)

# all cell name/bar code
rownames(studied_cells[["pca"]]@cell.embeddings)
colnames(studied_cells)

#all gene name
rownames(studied_cells) %>% View()
rownames(studied_cells[["pca"]]@feature.loadings)

#access Seurat metadata(nFeature,nCount,percentage.mito))
studied_cells@meta.data
#studied_cells$variablename could be added to metadata column

#access various meta data
studied_cells@meta.data$
  
studied_cells@reductions

studied_cells@assays$RNA

studied_cells$RNA@scale.data

studied_cells$RNA@counts

#acess DimReduc meta data table
studied_cells$pca@feature.loadings
studied_cells$tsne@cell.embeddings
studied_cells$umap@

  #access DimReduc
studied_cells[["pca"]]@cell.embeddings
studied_cells[["pca"]]@feature.loadings
studied_cells[["pca"]]@stdev

# create a df after pca, with pc coordinates for custom plot, set parameter as as.fator for color in ggplot
pca_embeddings_mat <-  studied_cells[["pca"]]@cell.embeddings
pca_embeddings_mat
pc1_vec <- pca_embeddings_mat[,1]
length(pc1_vec)
pc2_vec <- pca_embeddings_mat[,2]
df_with_matched_cell_name <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, condition = as.factor(matched_condition))

#access batch name/batch name + cell name  
studied_cells@meta.data$orig.ident
Idents(studied_cells)
rownames(studied_cells@meta.data)

colnames(studied_cells) = rownames(studied_cells@meta.data)

#assay slot:sacle.data object$RNA@scale.data shows normalized counts mtx with column as cell names, rownames as gene names
studied_cells$RNA@scale.data
str(studied_cells[["RNA"]]@scale.data)
rownames(studied_cells$RNA@scale.data)
cellnames(studied_cells$RNA@scale.data)

#extract specific columns from a table
meta.data %>% select(1, 4) %>% View()

#use match() to get exact order of cell type as cell name by position
i <- match(cell_name_to_get, cell_name_metadata)
matched_type_correct <- meta.data$type[i]

# clustering steps
studied_cells <- FindNeighbors(studied_cells, dims = 1:20)
studied_cells <- FindClusters(studied_cells, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4)) # explore diff resolution
studied_cells@meta.data %>% View() #look at obj@meta.data with added column of clusters after set reslution
head(Idents(studied_cells), 5)
Idents(studied_cells) <- "RNA_snn_res.0.6" #assign res=0.6 by assigning the identity of the clusters using the Idents() function.
head(Idents(studied_cells), 5) 

studied_cells <- RunUMAP(studied_cells, dims = 1:20)
DimPlot(studied_cells, reduction = "umap")
Idents(studied_cells) <- "RNA_snn_res.0.8"
DimPlot(studied_cells, reduction = "umap")

Idents(studied_cells) <- "RNA_snn_res.0.6"
studied_cells <- RunTSNE(studied_cells, dims = 1:20)
DimPlot(studied_cells, reduction = "tsne")
Idents(studied_cells) <- "RNA_snn_res.0.8"
DimPlot(studied_cells, reduction = "tsne")

#install "clustermole" for cell type id
install.packages("clustermole")
library(clustermole)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSVA")
install.packages("GSVA")
install.packages("htmltools")
my_genes = rownames(head(cluster2.markers, n = 10))
my_overlaps = clustermole_overlaps(genes = my_genes, species = "hs") %>% View()
my_overlaps

# assign cluster_cell_type function by iput of cluster number
cluster_cell_type.func <- function(i) {
  cluster_top_10_gene_tab <- studied_cells.markers %>% group_by(cluster) %>% filter(cluster == i) %>% top_n(n = 10, wt = avg_logFC) %>% select(gene)
  my_overlaps = clustermole_overlaps(genes = pull(cluster_top_10_gene_tab[, 2 ], gene), species = "hs")
  print(my_overlaps)
  
}


cluster_cell_type.func(10)

# if remove whole column, paste(), strsplit(). sapply
studied_cells.c0$RNA_snn_res.1.4 <- NULL
Idents(studied_cells.c0)
studied_cells.c0$sub_cluster0
studied_cells.c0$sub_cluster0[Cells(studied_cells.c0)] <- paste("c0",Idents(studied_cells.c0), sep='_')
strsplit(studied_cells.c0$sub_cluster0[Cells(studied_cells.c0)], '_')
sapply(strsplit(studied_cells.c0$sub_cluster0[Cells(studied_cells.c0)], '_'), '[[', 1)
DimPlot(studied_cells.c0, group.by = "sub_cluster0")

# Generate a new column called sub_cluster in the metadata
studied_cells$sub_cluster <- as.character(Idents(studied_cells))
studied_cells@meta.data %>% View()

# Change the information of cells containing sub-cluster information (incoporated to overall umap)
studied_cells$sub_cluster[Cells(studied_cells.c0)] <- paste("c0",Idents(studied_cells.c0))
DimPlot(studied_cells, group.by = "sub_cluster")

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


clu_ann[ , ] %>% View()
write.table(clu_ann[ , ], file = "re_0.2_cluster_cell_type.csv")

# string split
strsplit()

#regex
grep()
df_no_urine_cells <- df[-grep("U[0-9]{3}.{9}", df$my), ] #remove rows with specified pattern

#drop NA row
drop_na()

# subset specified rows and columns
small_counts[1:3, c("Sample_1", "Sample_3")]

# Sum the counts for each sample
sample_sums = apply(small_counts, MARGIN = 2, sum)
print(sample_sums)
#Sample_1 Sample_2 Sample_3 Sample_4 

#For each gene in ResultsTable_small$SYMBOL, does it appear in my_genes?
my_genes <- c("Smad7", "Wif1", "Fam102b", "Tppp3")
#You can use match to get the same subset, but this time make sure they are in the same order as in my_genes. This can be useful if you are trying to merge together data from two sources.

match(my_genes, ResultsTable_small$SYMBOL)
#[1] 20  1 33 37
#This returns the index for where each of the genes in my_list appears in ResultsTable_small$SYMBOL. We can then use this to subset the columns from ResultsTable_small.

#ResultsTable_small[match(my_genes, ResultsTable_small$SYMBOL), ]
#ENTREZID  SYMBOL     logFC  AveExpr         t      P.Value    adj.P.Val
#20    17131   Smad7  1.972771 6.717519  14.14348 6.493642e-09 5.131276e-06
#1     24117    Wif1  1.819943 2.975545  20.10780 1.063770e-10 1.016240e-06
#33   329739 Fam102b -1.520291 4.188130 -12.75357 2.120968e-08 1.015751e-05
#37    67971   Tppp3 -2.653398 4.908163 -12.22845 3.416616e-08 1.419445e-05
#8411     7942     8887     8098 

#find a which package a function is from
install.packages("sos")
library(sos)
???plot_grid()

# stringr package: cover base packages, str_split() could return df w/o using apply func(following are same output)
gene_id <- c('ENG000.2', 'ENG000.3', 'ENG000.4')
str_split(gene_id, '[.]', simplify = T)[,1] #df and extrat column 1
unlist(lapply(str_split(gene_id, '[.]'), function(x){x[1]}))
sapply(str_split(gene_id, '[.]'), '[[', 1)
#[1] "ENG000" "ENG000" "ENG000"

#str_replace_all()
library(stringr)
x <- "what kind of cheese do you like? [wonder] Swiss [groan] [laugh]"
str_replace_all(x, "\\[.*?\\]", "")
# "what kind of cheese do you like?  Swiss  "

#aggregate(X, by = list(g=), FUN)
data <- data.frame(
  s1=c(1, 2, 3),
  s2=c(4, 5, 6),
  s3=c(7, 8, 9),
  row.names = c("p1", "p2", "p3")
)

anno <- data.frame(
  g=c("g1", "g2", "g1"),
  row.names = c("p1", "p2", "p3")
)

new_data <- aggregate(
  data, by = list(g=anno$g), FUN = max
)#new col'g'(g1, g2) after aggregate

# 修改行名
rownames(new_data) <- new_data$g # rownames = g1, g2
new_data$g <- NULL #delte col 'g'

