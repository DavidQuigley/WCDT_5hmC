###########################################################################################################
###########################################################################################################
# Figure 2A,E clustering
# 2022-05-02
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1/"

fn_load_5hmc_environment <- createFn("/data/hmc_220428.RData")
###########################################################################################################
# Load environment
load(fn_load_5hmc_environment)
lapply(libraries_to_load, library, character.only = T)

###########################################################################################################
###########################################################################################################
# Figure 2A, dimension reduction
###########################################################################################################
###########################################################################################################
# Prep data
hmc_wcdt_tissue <- counts_wcdt_tissue_tpm[,c(samples_mcrpc_uniq_all_data, samples_wcdt_tissue_nt)]
hmc_benign <- counts_benign_prostate_tpm
colnames(hmc_benign) <- gsub("_5hmc","",colnames(hmc_benign))
hmc_localized <- counts_localized_tpm
full_mat <- as.matrix(cbind(hmc_wcdt_tissue,hmc_localized,hmc_benign))
full_mat <- log2(full_mat + 1)

all(rownames(full_mat)==rownames(ensembl2sym))
pc_mat <- full_mat[ensembl2sym$type=="protein_coding",]

q_split <- 0.9
vars <- apply(pc_mat,1,FUN = function(x) { var(x, na.rm = T)})
top_vars <- vars > quantile(vars, q_split, na.rm = T)
all(!is.na(top_vars))
mat <- pc_mat[top_vars,]

mat <- t(scale(t(mat)))

# Annotation
anno_df <- data.frame(sample = colnames(mat))
rownames(anno_df) <- anno_df$sample
anno_df$stage <- NA
anno_df$stage[grepl("DTB",anno_df$sample)] <- "mCRPC"
anno_df$stage[grepl("NT",anno_df$sample)] <- "NT"
anno_df$stage[grepl("CPCG",anno_df$sample)] <- "Localized"
anno_df$stage[grepl("BSG",anno_df$sample)] <- "Benign"

all(rownames(anno_df)==colnames(mat))

# PCA
res.pca <- prcomp(t(mat), scale = T)

theDFpc <- as.data.frame(res.pca$x)
all(rownames(theDFpc)==colnames(mat))
all(rownames(theDFpc)==rownames(anno_df))
theDFpc <- cbind(theDFpc,anno_df)

p <- ggplot(theDFpc, aes(x = PC1, y = PC2, shape = stage, col = stage)) +
  geom_point(size =3) + 
  scale_color_manual(values = as.vector(cols_disease_stage[c("Benign Prostate","Localized PCa","mCRPC","Adjacent tissue")])) + 
  scale_shape_manual(values = c(19,17,18,15)) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.text = element_text(size = 14),
        legend.title = element_blank())
p

pdf(file = paste0(dir_out_figures, "/Figure_2A.pdf"), width = 5, height = 5)
print(p)
dev.off()


###########################################################################################################
###########################################################################################################
# Figure 2E, hierarchical clustering
###########################################################################################################
###########################################################################################################
## Add more annotations
samples_plot <- samples_mcrpc_uniq_all_data

mat <- as.matrix(log2(counts_wcdt_tissue_tpm[,samples_plot]+1))
all(rownames(mat)==gencode.v28$gene_ensembl)
mat <- mat[gencode.v28$gene_type == "protein_coding",]

vars <- apply(mat,1,FUN = function(x) { var(x, na.rm = T)})
top_10 <- vars > quantile(vars,0.9,na.rm = T)
all(!is.na(top_10))
mat <- mat[top_10,]
mat <- t(scale(t(mat)))

d <- dist(t(mat), method = "euclidean") 
fit <- hclust(d, method="ward.D2")
plot(fit) # display dendogram
hmc_tissue_clusters <- cutree(fit, k=3) 
rect.hclust(fit, k=3, border="red")

annoDF <- metadata_wcdt_genomics
annoDF$select_muts <- FALSE
annoDF$select_muts[annoDF$CHD1] <- "CDH1"
annoDF$select_muts[annoDF$SPOP] <- "SPOP"
annoDF$select_muts[annoDF$IDH1] <- "IDH1"
annoDF$select_muts[annoDF$BRAF] <- "BRAF"
annoDF$select_muts[annoDF$TET2] <- "TET2"
annoDF$select_muts[annoDF$DNMT3B] <- "DNMT3B"

columns_keep <- c("site","AR_amp","TMPRSS2_ERG","RB1_2loss","PTEN_2loss","MYC_amp","TP53_2loss","CDK12_2loss","select_muts","CMP")
annoDF <- annoDF[colnames(mat),columns_keep]
all(rownames(annoDF)==colnames(mat))
all(names(hmc_tissue_clusters)==rownames(annoDF))
annoDF$cluster_5hmc <- hmc_tissue_clusters

columns_order <- c("cluster_5hmc","site","AR_amp","TMPRSS2_ERG","RB1_2loss","PTEN_2loss","MYC_amp","TP53_2loss","CDK12_2loss","select_muts","CMP")
annoDF <- annoDF[,columns_order]

table(annoDF$CMP,annoDF$cluster_5hmc)
fisher.test(table(annoDF$CMP,annoDF$cluster_5hmc)) #6.9x10^-12

ha <- HeatmapAnnotation(df=annoDF,
                        gap = unit(0.05, "in"),
                        border  = TRUE,
                        gp = gpar(col = "black",lwd = 0.5),
                        col = list(cluster_5hmc = cols_5hmc_clusters,
                                   select_muts = cols_select_mutations,
                                   PTEN_2loss = cols_true_false,
                                   CMP = cols_true_false,
                                   RB1_2loss = cols_true_false,
                                   CDK12_2loss = cols_true_false,
                                   MYC_amp = cols_true_false,
                                   AR_amp = cols_true_false,
                                   TP53_2loss = cols_true_false,
                                   TMPRSS2_ERG = cols_true_false,
                                   site = cols_met_site))

hm <- Heatmap(mat, 
              name = "5hmC levels",
              cluster_columns = function(x) hclust(dist(x),method = "ward.D2"),
              column_dend_height = unit(1, "in"),
              column_split = annoDF$cluster_5hmc,
              column_title = NULL,
              show_column_names = F,
              cluster_rows = function(x) hclust(dist(x),method = "ward.D2"),
              show_row_names = F,
              show_row_dend = F,
              top_annotation = ha,
              heatmap_height = unit(8, units = "in"))


pdf(paste0(dir_out_figures,"/Figure_2E.pdf"), height = 8.5, width = 12)
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left")
dev.off()