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

cell_name_metadata <- meta.data[["cell_name"]]
cell_name_to_get <- colnames(studied_cells)
cell_name_to_get_bool <- cell_name_metadata %in% cell_name_to_get # no need in correct code. use match()
match_cell_name <- cell_name_metadata[cell_name_to_get_bool] # no need in correct code. use match()

#correct here as instructed
i <- match(cell_name_to_get, cell_name_metadata)
matched_type_correct <- meta.data$type[i]
matched_type_correct

#original code
matched_metadata <- subset(meta.data, cell_name[1:length(meta.data$cell_name)] %in% match_cell_name)
matched_type_incorrect <- matched_metadata$type[1:length(matched_metadata$type)]
matched_type_incorrect

#logic test but seems same
matched_type_correct == matched_type_incorrect


pca_embeddings_mat <-  studied_cells[["pca"]]@cell.embeddings
pc1_vec <- pca_embeddings_mat[,1]
pc2_vec <- pca_embeddings_mat[,2]
df_with_matched_type <- data.frame(pc1 = pc1_vec, pc2 = pc2_vec, type = as.factor(matched_type_correct))

output_dir = '/Users/zunqiuwang/Desktop/'
fig_file_type <- paste0(output_dir, 'pc1_vs_pc2_coorect_type_color.png')
fig_file_type
png(fig_file_type)

plot_match_cell_type <- ggplot(df_with_matched_type, aes(x = pc1, y = pc2, color = type)) + 
  geom_point(size = 1) + 
  labs(x = 'pc1', y = 'pc2') +
  theme(axis.title.x = element_text(face="bold", size=20), 
        axis.title.y = element_text(face="bold", size=20))
print(plot_match_cell_type)
dev.off()
