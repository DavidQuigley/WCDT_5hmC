###########################################################################################################
# Figure 1C,D,E, gene expression model
# 2022-04-28
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1/"

fn_load_5hmc_environment <- createFn("/data/hmc_220428.RData")

fn_correlation_table_out <- createFn("results/correlation_table.txt")
###########################################################################################################
# Load environment
load(fn_load_5hmc_environment)
lapply(libraries_to_load,require,character.only = T)

###########################################################################################################
###########################################################################################################
# Analysis
###########################################################################################################
###########################################################################################################

###########################################################################################################
# Model gene expresssion
###########################################################################################################
# Model GEX based on 5hmC
# 93 mCRPC tissue samples with all data

matrix_hmc <- as.matrix(counts_wcdt_tissue_tpm[,samples_mcrpc_uniq_all_data])
matrix_gex <- as.matrix(tpm[,samples_mcrpc_uniq_all_data])
matrix_cna <- as.matrix(matrix_CNA_gencode28[,samples_mcrpc_uniq_all_data])
matrix_pm <- as.matrix(promoter_methylation[,samples_mcrpc_uniq_all_data])

all(rownames(matrix_hmc)==rownames(matrix_gex))
all(rownames(matrix_hmc)==rownames(matrix_cna))
all(rownames(matrix_hmc) == gencode.v28$gene_ensembl)

all(rownames(matrix_hmc)==rownames(matrix_pm))
gencode.v28[gencode.v28$gene_ensembl %in% rownames(matrix_hmc)[! rownames(matrix_hmc) %in% rownames(matrix_pm)],]
# Methylation data is missing 37 genes that are on the mitochondrial chromosome, keep as NA. 

all(colnames(matrix_hmc)==colnames(matrix_gex))
all(colnames(matrix_hmc)==colnames(matrix_cna))
all(colnames(matrix_hmc)==colnames(matrix_pm))

samples_all <- colnames(matrix_hmc)

###########################################################################################################
# This is based on the data tables used in the Cell 2018 paper.
gof_names = c(
  'activating_missense',
  'activating_sv')
lof_names = c(
  'inactivating_missense',
  'nonsense',
  'inactivating_germline',
  'inactivating_sv')
sv_activating <- c("fusion","tandem")
sv_inactivating <- c("breakend","deletion","inversion")

###########################################################################################################
# gex_correlation <- as.data.frame(matrix(NA,nrow=nrow(matrix_hmc),ncol = 13))
# 
# colnames(gex_correlation) <- c("gene","ensembl_id",
#                                "cor_5hmc","p_5hmc",
#                                "cor_cna", "p_cna",
#                                "cor_met", "p_met",
#                                "lm_5hmc_adj","p_5hmc_adj",
#                                "adj_r_sq_no_5hmc","adj_r_sq_5hmc",
#                                "p_improvement")
# 
# for(i in 1:nrow(gex_correlation)){
# 
#   gene <- rownames(matrix_hmc)[i]
#   gene_symbol <- gencode.v28[gencode.v28$gene_ensembl==gene,"gene_name"]
# 
#   gex_correlation[i,"ensembl_id"] <- gene
#   gex_correlation[i,"gene"] <- gene_symbol
# 
#   if(i %% 1000 == 0){
#     print(paste("Done with",i,"genes!"))
#   }
# 
#   gex <- log2(as.numeric(matrix_gex[gene,samples_all,drop=T])+1)
#   cna <- as.numeric(matrix_cna[gene,samples_all])
#   hmc <- log2(as.numeric(matrix_hmc[gene,samples_all,drop=T])+1)
# 
#   # Force missing genes to be NA
#   if(!gene %in% rownames(matrix_pm)){
#     met <- rep(NA,length(samples_all))
#   }else{
#     met <- as.numeric(promoter_methylation[gene,samples_all])
#   }
# 
#   # Create z-scores to have comparable coefficients. Does not affect rank-based correlations.
#   hmc_scale <- scale(hmc)
#   gex_scale <- scale(gex)
# 
#   # Prepare SNV and SV data
#   
#   if(gene_symbol %in% setdiff(unique(c(curated_fs$symbol, curated_missense$gene, curated_sv$fiveprime,curated_sv$threeprime)),"RP11-356O9.1")){
# 
#     all_eff = allele_effect(gene_symbol)$alleles
#     all_eff = all_eff[samples_all, ]
# 
#     lof = apply(all_eff[,lof_names],1,sum)
#     gof = apply(all_eff[,gof_names],1,sum)
# 
# 
#   }else if(gene_symbol %in% rownames(matrix_missense)){
# 
#    has_missense <- matrix_missense[gene_symbol,samples_all]!="." # Using Strelka calls
#    has_activating_sv <- matrix_SV[gene_symbol,samples_all] %in% sv_activating # Using lumpy calls
# 
#    has_inactivating_somatic_snv <-  matrix_inactive[gene_symbol,samples_all]!="."
#    has_inactivating_germline <- matrix_germline[gene_symbol,samples_all] != "."
#    has_inactivating_sv <- matrix_SV[gene_symbol,samples_all] %in% sv_inactivating
# 
#    lof <- has_inactivating_germline + has_inactivating_somatic_snv + has_activating_sv
#    gof <- has_missense + has_activating_sv
# 
#   # Force missing to be 0.
#   }else{
#     lof <- gof <- rep(0,length(samples_all))
#     names(gof) <- names(lof) <- samples_all
#   }
# 
#   # Set genes with missing data or invariable 0 to NA
#   if( (sum(gex!=0)) < 4 | (sum(cna!=0)) < 4 | (sum(hmc!=0)) < 4 | (sum(met!=0, na.rm=T)) < 4 | (sum(!is.na(met)) < 10)){
#     gex_correlation[i,3:ncol(gex_correlation)] <- NA
#     next()
#   }
# 
#   cor_5hmc <- cor(hmc_scale,gex_scale,method = "spearman")
#   p_5hmc <- cor.test(hmc_scale,gex_scale,method = "spearman")$p.value
#   cor_cna <- cor(cna,gex_scale,method="spearman")
#   p_cna <- cor.test(cna,gex_scale,method = "spearman")$p.value
#   cor_met <- cor(met,gex_scale,method="spearman", use = "na.or.complete")
#   p_met <- cor.test(met,gex_scale,method="spearman")$p.value
# 
#   lm_fit <- summary(lm(gex_scale~hmc_scale+cna+met+gof+lof))
#   lm_5hmc_adj <- lm_fit$coefficients["hmc_scale","Estimate"]
#   p_5hmc_adj <- lm_fit$coefficients["hmc_scale","Pr(>|t|)"]
#  
#   model1 <- lm(gex_scale~met+cna+lof+gof)
#   model2 <- lm(gex_scale~met+cna+lof+gof+hmc_scale)
# 
#   anova_fit <- anova(model1,model2)
# 
#   p_improvement <- anova_fit$`Pr(>F)`[2]
#   
#   adj_r_sq_no_5hmc <- summary(model1)$adj.r.squared
#   adj_r_sq_5hmc <- summary(model2)$adj.r.squared
# 
#   gex_correlation[i,"cor_5hmc"] <- cor_5hmc
#   gex_correlation[i,"p_5hmc"] <- p_5hmc
#   gex_correlation[i,"cor_cna"] <- cor_cna
#   gex_correlation[i,"p_cna"] <- p_cna
#   gex_correlation[i,"cor_met"] <- cor_met
#   gex_correlation[i,"p_met"] <- p_met
#   gex_correlation[i,"lm_5hmc_adj"] <- lm_5hmc_adj
#   gex_correlation[i,"p_5hmc_adj"] <- p_5hmc_adj
#   gex_correlation[i,"p_improvement"] <- p_improvement
#   gex_correlation[i,"adj_r_sq_no_5hmc"] <- adj_r_sq_no_5hmc
#   gex_correlation[i,"adj_r_sq_5hmc"] <- adj_r_sq_5hmc
# }
# 
# 
# gex_correlation$adj_fit_improvement <- gex_correlation$adj_r_sq_5hmc-gex_correlation$adj_r_sq_no_5hmc
# gex_correlation$adj_fit_improvement_perc <- (gex_correlation$adj_r_sq_5hmc-gex_correlation$adj_r_sq_no_5hmc)/gex_correlation$adj_r_sq_no_5hmc

# Output correlation tablee
# write.table(gex_correlation, fn_correlation_table_out)
gex_correlation <- read.table(fn_correlation_table_out, stringsAsFactors = F)
############################################################################################################
############################################################################################################
# Plotting
############################################################################################################
############################################################################################################

############################################################################################################
# Boxplot for model fit
############################################################################################################
# Gene subsets
all(gex_correlation$ensembl_id == gencode.v28$gene_ensembl)
gex_correlation$ensembl_base <- gsub("\\..+$","",gex_correlation$ensembl_id)

gex_correlation$is.protein.coding <- gex_correlation$ensembl_id %in% gencode.v28$gene_ensembl[gencode.v28$gene_type=="protein_coding"]
gex_correlation$is.prostate.specific <- gex_correlation$ensembl_base %in% genelist_pca_specific_overexpr$ensembl_base 
gex_correlation$is.ar.response <- gex_correlation$gene %in% hallmark_pathways[["HALLMARK_ANDROGEN_RESPONSE"]] 

df_boxplot <- rbind(data.frame(adj_rsq = gex_correlation$adj_r_sq_no_5hmc, genes =  "All", model = "no_hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_5hmc, genes =  "All", model = "hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_no_5hmc[gex_correlation$is.protein.coding], genes =  "Protein coding", model = "no_hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_5hmc[gex_correlation$is.protein.coding], genes =  "Protein coding", model = "hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_no_5hmc[gex_correlation$is.prostate.specific], genes =  "PCa specific", model = "no_hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_5hmc[gex_correlation$is.prostate.specific], genes =  "PCa specific", model = "hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_no_5hmc[gex_correlation$is.ar.response], genes =  "AR response", model = "no_hmc"),
                    data.frame(adj_rsq = gex_correlation$adj_r_sq_5hmc[gex_correlation$is.ar.response], genes =  "AR response", model = "hmc"))

p <- ggplot(df_boxplot, aes(x = genes, y = adj_rsq, fill = model)) +
  geom_boxplot() +
  theme_minimal() + 
  scale_fill_manual(values = c("gray","dodgerblue"), name = "Model", labels = c("CN + PM + SNV + SV","CN + PM + SNV + SV + 5hmC")) +
  theme(axis.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size  = 16),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = c(0.65,0.08),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(title = "" , x = "", y = expression("Adjusted R"^"2"))
p

#pdf(paste0(dir_out_figures,"/Figure_1D.pdf"), width = 5, height = 7)
print(p)
#dev.off()

# Extract stats
median(gex_correlation[gex_correlation$is.protein.coding,"adj_r_sq_no_5hmc"], na.rm = T) # 9.72%
median(gex_correlation[gex_correlation$is.protein.coding,"adj_r_sq_5hmc"], na.rm = T) # 21.79%
# 12.07 improvement in median between groups
wilcox.test(gex_correlation[gex_correlation$is.protein.coding,"adj_r_sq_no_5hmc"],gex_correlation[gex_correlation$is.protein.coding,"adj_r_sq_5hmc"])$p.value # 0

median(gex_correlation[gex_correlation$is.ar.response,"adj_r_sq_no_5hmc"], na.rm = T) #14.8
median(gex_correlation[gex_correlation$is.ar.response,"adj_r_sq_5hmc"], na.rm = T) # 42.1
wilcox.test(gex_correlation[gex_correlation$is.ar.response,"adj_r_sq_no_5hmc"],gex_correlation[gex_correlation$is.ar.response,"adj_r_sq_5hmc"])$p.value #1.53217e-18

median(gex_correlation[gex_correlation$is.prostate.specific,"adj_r_sq_no_5hmc"], na.rm = T) # 0.118
median(gex_correlation[gex_correlation$is.prostate.specific,"adj_r_sq_5hmc"], na.rm = T) # 0.3519
wilcox.test(gex_correlation[gex_correlation$is.prostate.specific,"adj_r_sq_no_5hmc"],gex_correlation[gex_correlation$is.prostate.specific,"adj_r_sq_5hmc"])$p.value #1.063743e-17

### Plot for all Hallmark pathways

for(i in 1:length(hallmark_pathways)){
  thePathway <- names(hallmark_pathways)[i]
  theGenes <- hallmark_pathways[[i]]
  HMC <-  gex_correlation[gex_correlation$gene %in% theGenes,"adj_r_sq_5hmc"]
  noHMC <- gex_correlation[gex_correlation$gene %in% theGenes,"adj_r_sq_no_5hmc"]
  theDF <- rbind(data.frame(pathway = thePathway,
                            value = HMC,
                            type = "HMC"),
                 data.frame(pathway = thePathway,
                            value = noHMC,
                            type = "noHMC")
  )
  if(i == 1){
    improvement_all <-  theDF
  }else{
    improvement_all <- rbind(improvement_all,
                             theDF)
  }
}

improvement_all$type <- factor(improvement_all$type,levels = c("noHMC","HMC"))

theOrder <- improvement_all %>%
  dplyr::filter(type == "HMC") %>%
  dplyr::group_by(pathway) %>%
  dplyr::summarise(median(value, na.rm=T))

improvement_all$pathway <- factor(improvement_all$pathway, levels = levels(reorder(theOrder$pathway,theOrder$`median(value, na.rm = T)`)))

p <- ggplot(improvement_all, aes(x = pathway, y = value, fill = type)) +
  geom_boxplot() +
  theme_minimal() + 
  scale_fill_manual(values = c("gray","dodgerblue"), name = "Model", labels = c("CN + PM + SNV + SV","CN + PM + SNV + SV + 5hmC")) +
  theme(axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size  = 16),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = c(0.25,0.08),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  labs(title = "" , x = "", y = expression("Adjusted R"^"2"))
p

#pdf(paste0(dir_out_figures,"/Figure_S3B.pdf"), width = 15, height = 15)
print(p)
#dev.off()


############################################################################################################
# Volcano plot 
############################################################################################################
# 
select_genes <- c("AR")

gex_correlation_pc <- gex_correlation[gex_correlation$is.protein.coding,]
gex_correlation_pc <- gex_correlation_pc[!(is.na(gex_correlation_pc$lm_5hmc_adj) | is.na(gex_correlation_pc$p_5hmc_adj)),] #1369 with missing
gex_correlation_pc$select_genes <- gex_correlation_pc$gene %in% select_genes 
gex_correlation_pc$is.ar.response <- gex_correlation_pc$gene %in% c(hallmark_pathways[["HALLMARK_ANDROGEN_RESPONSE"]],"AR")

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

gex_correlation_pc$density <- get_density(gex_correlation_pc$lm_5hmc_adj, -log10(gex_correlation_pc$p_5hmc_adj), n = 1000)
myPalette <- colorRampPalette(rev(brewer.pal(11, "PuOr")))

pval <- signif(wilcox.test(gex_correlation_pc$lm_5hmc_adj[gex_correlation_pc$is.ar.response],gex_correlation_pc$lm_5hmc_adj[!gex_correlation_pc$is.ar.response])$p.value,2)

p1 <- ggplot(gex_correlation_pc, aes(x=lm_5hmc_adj, y=-log10(p_5hmc_adj), label = gene)) +
  geom_point(aes(color = density), alpha = 0.1, data = gex_correlation_pc[!gex_correlation_pc$is.ar.response,]) +
  geom_point(color = "black", alpha = 1, data = gex_correlation_pc[gex_correlation_pc$is.ar.response,], aes(fill = "Androgen Response gene")) +
  scale_fill_discrete(name = "") + 
  scale_color_viridis(option = "plasma", name = "Density") +
  geom_text_repel(data = subset(gex_correlation_pc, select_genes), max.overlaps = 100, size = 10) +
#  xlim(-0.6,1.3) +
#  ylim(0,22) +
  labs(title = "", 
       x = "Scaled 5hmC coeffecient", y = "-log10(p improvement of model fit)") +
  annotate("text", x = -0.3, y = 20, label = paste0("AR response vs.\n other protein coding\np = ", pval), size = 6) +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = c(0.3,0.6),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16))

plist = list()
plist[[1]] <- p1

gex_correlation_pc$is.ar.response <- factor(gex_correlation_pc$is.ar.response, levels = c(FALSE,TRUE),labels = c("Other","AR response"))

test = wilcox.test(gex_correlation_pc$lm_5hmc_adj[gex_correlation_pc$is.ar.response == "AR response"],
                   gex_correlation_pc$lm_5hmc_adj[gex_correlation_pc$is.ar.response == "Other"])

p2 <- ggplot(gex_correlation_pc,aes(x=is.ar.response,y=lm_5hmc_adj)) +
  xlab('') + 
  ylab('') +
  geom_boxplot() +
  theme_classic() + 
  theme(axis.text.y=element_text(size = 12),
        axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.line=element_blank()) + 
  coord_flip()

plist[[2]] <- p2

p <- ggarrange(plotlist=plist,heights=c(4,1), ncol = 1, nrow = 2, align = "v")

#pdf(paste0(dir_out_figures,"/Figure 1E.pdf"), width = 7, height = 9)
print(p)
#dev.off()


############################################################################################################
# Calculate how many genes have a significant improvment
############################################################################################################
thePs <- gex_correlation_pc$p_5hmc_adj
table(thePs < 0.05, useNA = "ifany")
theQs <- p.adjust(thePs, method = "fdr")
table(theQs < 0.05, useNA = "ifany") # 10637/18532 = 57%

############################################################################################################
# Density plot of gene-wise correlations 
############################################################################################################
# Plot comparing 5hmC and promoter methylation
df_density <- data.frame(corr = c(gex_correlation_pc$cor_5hmc,
                                  gex_correlation_pc$cor_met,
                                  gex_correlation_pc$cor_cna),
                         type = c(rep("hmc",nrow(gex_correlation_pc)),
                                  rep("met",nrow(gex_correlation_pc)),
                                  rep("cna",nrow(gex_correlation_pc)))
                         )

df_density$type <- factor(df_density$type, levels = c("met","hmc","cna"))

peak_met <- median(df_density$corr[df_density$type=="met"],na.rm=T)
peak_hmc <- median(df_density$corr[df_density$type=="hmc"],na.rm=T)
peak_cna <- median(df_density$corr[df_density$type=="cna"],na.rm=T)

# Test for  difference, reverse methylation
wilcox.test(df_density$corr[df_density$type == "hmc"],df_density$corr[df_density$type == "cna"]) # < 2.2e-16
wilcox.test(df_density$corr[df_density$type == "hmc"],-df_density$corr[df_density$type == "met"]) # < 2.2e-16
wilcox.test(df_density$corr[df_density$type == "cna"],-df_density$corr[df_density$type == "met"]) # < 2.2e-16


p <- ggplot(df_density, aes(x = corr, fill = type))+
  geom_density(alpha = 0.4)+
  geom_vline(xintercept = peak_cna, color = cols_data_types["CNA"], alpha = 0.5, linetype = "dashed", size = 2) +
  geom_vline(xintercept = peak_met, color = cols_data_types["PM"], alpha = 0.5, linetype = "dashed", size = 2) +
  geom_vline(xintercept = peak_hmc, color = cols_data_types["HMC"], alpha = 0.5, linetype = "dashed", size = 2) +
  labs(title = "",
       y = "Density", x = "Spearman's Rho")+ 
  scale_fill_manual(values  = c(cols_data_types[["PM"]], cols_data_types[["HMC"]],cols_data_types[["CNA"]]), 
                    name = "", labels = c("Promoter methylation","5hmC gene body counts","Copy number") ) +
  theme_minimal() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        legend.position = c(0.2,0.8),
        legend.text = element_text(size = 14)) +
  annotate("text", x = -0.63, y = 1, label = paste0("P (all comparisons) < 2.2e-16"), size = 4) 
p

pdf(paste0(dir_out_figures,"/Figure_1C.pdf"), width = 8, height = 5)
print(p)
dev.off()



###########################################################################################################
# GSEA for correlation strength
###########################################################################################################
genes_rank <- gex_correlation_pc[,"cor_5hmc"]
names(genes_rank) <- gex_correlation_pc$gene
genes_rank <- genes_rank[!is.na(genes_rank)]

set.seed(12)
fgseaRes <- fgseaMultilevel(pathways = hallmark_pathways, stats = genes_rank, maxSize = 500, nproc = 6)

plotDF <- as.data.frame(fgseaRes)
rownames(plotDF) <- plotDF$pathway

thePlot <- ggplot(plotDF, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill = padj < 0.05)) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "",fill= paste0("FDR < ",0.05)) +
  scale_fill_manual(values = c("slategray3","forestgreen"),breaks = c("FALSE","TRUE")) +
  theme_few() +
  theme(legend.position = c(0.8,0.3))

print(thePlot)

pdf(paste0(dir_out_figures,"/Figure_S3.pdf"),width = 10, height = 10)
print(thePlot)
dev.off()
