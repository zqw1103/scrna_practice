meta.data <- read.table("SDY997_EXP15176_celseq_meta.tsv.725704", header = TRUE)
gene_by_cell_mat<- read.table("SDY997_EXP15176_celseq_matrix_ru10_molecules.tsv.725699.gz", header= TRUE, sep="\t", row.names = 1, as.is = TRUE)
gene_by_cell_mat[is.na(gene_by_cell_mat)] <- 0
studied_cells <- CreateSeuratObject(counts = gene_by_cell_mat, project = "gene_by_cell_mat", min.cells = 3, min.features = 100)
studied_cells[["percent.mt"]] <- PercentageFeatureSet(studied_cells, pattern = "^MT-")
studied_cells <- subset(studied_cells, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 25)
studied_cells <- NormalizeData(studied_cells)
studied_cells <- FindVariableFeatures(studied_cells, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(studied_cells)
all_cells_names <- colnames(studied_cells)
studied_cells <- ScaleData(studied_cells, features = all.genes)
studied_cells <- RunPCA(studied_cells, features = VariableFeatures(object = studied_cells))

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
