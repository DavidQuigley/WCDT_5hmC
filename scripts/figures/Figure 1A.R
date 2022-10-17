###########################################################################################################
# Figure 1A, locations for HMRs and 5hmC peaks.
# 2025-05-05
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce/cancer_research_revision_1/"

fn_load_5hmc_environment <- createFn("/data/hmc_220428.RData")
###########################################################################################################
###########################################################################################################
# Load environment
load(fn_load_5hmc_environment)
lapply(libraries_to_load, library, character.only = T)
###########################################################################################################
###########################################################################################################
# 5hmC peaks 
###########################################################################################################
###########################################################################################################
# Annotate location
peakAnnoList_wcdt_tissue <- lapply(peaks_wcdt_tissue[samples_mcrpc_uniq_all_data], annotatePeak, tssRegion = c(-2000,500), TxDb=gencode.v28.txdb) 
###########################################################################################################
# Plot peak location
locationData_peaks_wcdt_tissue <- data.frame(matrix(NA,nrow=10, ncol = 0))
rownames(locationData_peaks_wcdt_tissue) <- peakAnnoList_wcdt_tissue[[1]]@annoStat[,"Feature"]

for(i in 1:length(peakAnnoList_wcdt_tissue)){
  sampleid <- names(peakAnnoList_wcdt_tissue)[i]
  locationData_peaks_wcdt_tissue[,sampleid] <- peakAnnoList_wcdt_tissue[[i]]@annoStat[,2]
}

locationData <- data.frame(matrix(NA,nrow=10, ncol = 0))
rownames(locationData) <- peakAnnoList_wcdt_tissue[[1]]@annoStat[,"Feature"]

locationData[,"mcrpc"] <- rowMeans(locationData_peaks_wcdt_tissue)

df <- locationData
df$feature <- rownames(df)

df$feature_plot <- c("Promoter","Promoter","5'UTR", "3'UTR","Exon","Exon", "Intron", "Intron", "Downstream","Distal Intergenic")

df <- aggregate(df$mcrpc, by = list(feature = df$feature_plot), FUN = sum)
colnames(df)[2] <- "value"

df <- df %>%
  arrange(desc(feature)) %>%
  mutate(lab.ypos = cumsum(value) - 0.5*value)
df

df[df$feature=="5'UTR","lab.ypos"] <- 95
df[df$feature=="3'UTR","lab.ypos"] <- 100
df[df$feature=="Downstream (<=300)","lab.ypos"] <- 83

mycols <- rev(brewer.pal(7, name = "Set2"))

p <- ggplot(df, aes(x = 2, y = value, fill = feature)) + 
  geom_bar(stat = "identity", color = "black") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab.ypos, label = feature), color = "black", size = 5, fontface = "bold") +
  scale_fill_manual(values = mycols) +
  theme_void() +
  xlim(0.5, 2.5) +
  labs(fill = "Genomic Feature")
p

#png(paste0(dir_out_figures,"/Figure_1A_hmc_peaks.png"),height = 7, width = 7, units = "in", res = 300)
print(p)
#dev.off()

# Print percentages
df # gene body (excluding promoter overlap) 59%, promoter 18%

###########################################################################################################
###########################################################################################################
# HMR location
###########################################################################################################
###########################################################################################################
hmrAnnoList_wcdt_tissue <- lapply(hmr[samples_mcrpc_uniq_all_data], annotatePeak, tssRegion = c(-2000,500), TxDb=gencode.v28.txdb) 
###########################################################################################################
# Plot
locationData_hmr_wcdt_tissue <- data.frame(matrix(NA,nrow=10, ncol = 0))
rownames(locationData_hmr_wcdt_tissue) <- hmrAnnoList_wcdt_tissue[[1]]@annoStat[,"Feature"]

for(i in 1:length(hmrAnnoList_wcdt_tissue)){
  sampleid <- names(hmrAnnoList_wcdt_tissue)[i]
  locationData_hmr_wcdt_tissue[,sampleid] <-hmrAnnoList_wcdt_tissue[[i]]@annoStat[,2]
}

locationData_hmr <- data.frame(matrix(NA,nrow=10, ncol = 0))
rownames(locationData_hmr) <- hmrAnnoList_wcdt_tissue[[1]]@annoStat[,"Feature"]

locationData_hmr[,"mcrpc"] <- rowMeans(locationData_hmr_wcdt_tissue)

df_hmr <- locationData_hmr
df_hmr$feature <- rownames(df_hmr)

df_hmr$feature_plot <- c("Promoter","Promoter","5'UTR", "3'UTR","Exon","Exon", "Intron", "Intron", "Downstream","Distal Intergenic")

df_hmr <- aggregate(df_hmr$mcrpc, by = list(feature = df_hmr$feature_plot), FUN = sum)
colnames(df_hmr)[2] <- "value"

df_hmr <- df_hmr %>%
  arrange(desc(feature)) %>%
  mutate(lab.ypos = cumsum(value) - 0.5*value)
df_hmr

df_hmr[df_hmr$feature=="5'UTR","lab.ypos"] <- 95
df_hmr[df_hmr$feature=="3'UTR","lab.ypos"] <- 100
df_hmr[df_hmr$feature=="Downstream (<=300)","lab.ypos"] <- 83

mycols <- rev(brewer.pal(7, name = "Set2"))

p <- ggplot(df_hmr, aes(x = 2, y = value, fill = feature)) + 
  geom_bar(stat = "identity", color = "black") +
  coord_polar(theta = "y", start = 0) +
  geom_text(aes(y = lab.ypos, label = feature), color = "black", size = 5, fontface = "bold") +
  scale_fill_manual(values = mycols) +
  theme_void() +
  xlim(0.5, 2.5) +
  labs(fill = "Genomic Feature")
p

png(paste0(dir_out_figures,"/Figure_1A_hmr_peaks.png"),height = 7, width = 7, units = "in", res = 300)
print(p)
dev.off()

# Print percentages 
df_hmr # gene body (excluding promoter overlap) 27%  promoter 47% 

###########################################################################################################
###########################################################################################################
# Test differences in locations between 5hmC peaks and HMRs
###########################################################################################################
###########################################################################################################
# Genebody
genebody_hmr <- colSums(locationData_hmr_wcdt_tissue[c("1st Exon","Other Exon","1st Intron","Other Intron"),samples_mcrpc_uniq_all_data])
genebody_hmc <- colSums(locationData_peaks_wcdt_tissue[c("1st Exon","Other Exon","1st Intron","Other Intron"),samples_mcrpc_uniq_all_data])

wilcox.test(genebody_hmr,genebody_hmc)$p.value # p-value = 5.066495e-32

# Promtoer
promoter_hmr <- colSums(locationData_hmr_wcdt_tissue[c("Promoter (<=1kb)","Promoter (1-2kb)"),samples_mcrpc_uniq_all_data])
promoter_hmc <- colSums(locationData_peaks_wcdt_tissue[c("Promoter (<=1kb)","Promoter (1-2kb)"),samples_mcrpc_uniq_all_data])

wilcox.test(promoter_hmr,promoter_hmc)$p.value # p-value = 5.066495e-32


