###########################################################################################################
# Figure 5E
# 2022-01-06
# Martin Sjoestroem
###########################################################################################################
###########################################################################################################
# Directories and filenames
createFn <- function(fn, dirBase = dir_base) {paste0(dirBase, fn)}

dir_base <- "/home/martin/mnt/data1/projects/WCDT_5hmc_2019/publication_reproduce"

fn_load_5hmc_environment <- createFn("/data/hmc_211221.RData")
###########################################################################################################
# Load environment
###########################################################################################################
load(fn_load_5hmc_environment)
lapply(libraries_to_load,require,character.only = T)
###########################################################################################################
###########################################################################################################
# Plot first 10, sorted by p-val
nTop <- 10

theDF <- homer_res_erg[order(homer_res_erg$P.value)[1:nTop],]
theDF$Motif <- gsub("/.+$","",theDF$Motif.Name)
table(duplicated(theDF$Motif))

p <- ggplot(theDF, aes(x = reorder(Motif, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = "identity", width = 0.5) +
  coord_flip() + 
  theme_classic() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
p

pdf(file = paste0(dir_out_figures,"/Figure_5E.pdf"), width = 7, height = nTop)
print(p)
dev.off()
    
  

