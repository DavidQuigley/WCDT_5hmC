regionType = c("Hypomethylated regions","H3K27me3","H3K27ac, in promoter","H3K27ac, outside promoter","H3K4me3, in promoter","CpG islands","Genebody")
# View(cbind(names(data_ngs_plot_enrich),regionType)) # Names match

### Functions
#######

getMultiFeatureData <- function(dataList, featureName, samplesToUse, nBins = 101){
  theList = dataList[[featureName]][samplesToUse]
  theArray = abind(theList,along=3)
  theMeans = apply(theArray, c(1,2), FUN = function(x) mean(x,na.rm = T))
  dataToPlot = as.data.frame(theMeans)
  dataToPlot$x = c(1:nBins)
  return(dataToPlot)
}

plotMultiFeature <- function(dataToPlot, nLevels = 3, levelNames = c("Low","Medium","High"), 
                             regionType = "Genebody", startPos = 21, endPos = 81,
                             xLab = "Genomic position", xBreaks = c(0,21,41,61,81,101), xLabels = c("-5kb", "TSS", "33%", "67%","TES","+5kb"), 
                             legendXpos = 0.1, legendYpos = 0.8, legendName = "Expression level",
                             yLab = "Log2(Fold change over control)", nColors = 3, plotHline = T,
                             savePlot = F, outDir = NULL, fileName = NULL, fileW = 7, fileH = 5)
{
  theTitle = regionType
  theColors = c(viridis_pal(option = "viridis")(nColors)[-nColors],'gold4')
  p <- ggplot(dataToPlot, aes(x = x)) + 
    geom_vline(xintercept = startPos, linetype = "dotted", color = "red", size = 1, alpha = 0.5) +
    geom_vline(xintercept = endPos, linetype = "dotted", color = "red", size = 1, alpha = 0.5) +
    labs(title = theTitle, y = yLab, x = xLab,size=0.5) +
    scale_x_continuous(breaks = xBreaks, labels = xLabels) +
    theme_few()
  if(plotHline){
    p <- p + geom_hline(yintercept = 0, linetype = "dotted", size = 1, alpha = 0.5, color = "gray")
  }
  if(nLevels == 3){
    p <- p + geom_spline(aes(y = dataToPlot[,1], color = theColors[3]), size = 1) +
      geom_spline(aes(y = dataToPlot[,2], color = theColors[2]), size = 1) +
      geom_spline(aes(y = dataToPlot[,3], color = theColors[1]), size = 1) +
      scale_color_identity(name = legendName, breaks = theColors, labels = rev(levelNames), guide = "legend") +
      theme(legend.position = c(legendXpos,legendYpos), legend.background = element_blank(),legend.text = element_text(size=6),plot.tag.position=c(1.02, 0.5),plot.margin = unit(c(1, 1, 1, 1), "cm"))
  }
  
  if(nLevels == 5){
    p <- p + geom_spline(aes(y = dataToPlot[,1], color = theColors[5]), size = 2) +
      geom_spline(aes(y = dataToPlot[,2], color = theColors[4]), size = 2) +
      geom_spline(aes(y = dataToPlot[,3], color = theColors[3]), size = 2) +
      geom_spline(aes(y = dataToPlot[,4], color = theColors[2]), size = 2) +
      geom_spline(aes(y = dataToPlot[,5], color = theColors[1]), size = 2) +
      scale_color_identity(name = legendName, breaks = theColors, labels = rev(levelNames), guide = "legend") +
      theme(legend.position = c(legendXpos,legendYpos), legend.background = element_blank(),legend.text = element_text(size=6),plot.tag.position=c(1.02, 0.5),plot.margin = unit(c(1, 1, 1, 1), "cm"))
  }
  print(p)}

#######
Genebody_plot <- plotMultiFeature(dataToPlot = getMultiFeatureData(dataList = data_ngs_plot_enrich,
                                                                   samplesToUse = samples_mcrpc_uniq_all_data, 
                                                                   featureName = "GENEBODY"),
                                  nLevels = 5, nColors = 5,
                                  levelNames = c("0-20%","21-40%","41-60%","61-80%","81-100%"),
                                  yLab = "5hmC enrichment",
                                  xLab = "Gene body",
                                  legendName = "Expression Rank",
                                  regionType = "") +
  theme(plot.title = element_text(size=30,face='bold'),
        axis.text = element_text(size=20),
        axis.title = element_text(size=30),
        legend.title =element_text(size=20),
        legend.position = 'right',
        legend.text = element_text(size=18),
        plot.tag.position=c(1.02, 0.5),
        plot.margin = unit(c(1, 1, 1, 1), "cm"))

Genebody_plot

plot_list = list()
for (i in 1:6) {
  p = plotMultiFeature(dataToPlot = getMultiFeatureData(dataList = data_ngs_plot_enrich, featureName = names(data_ngs_plot_enrich)[i]),
                       nLevels = 3, nColors = 3, yLab = "5hmC enrichment",
                       startPos = 41, endPos = 61,regionType=regionType[i],
                       xLabels = c("-5kb","-2.5kb","Start","End","+2.5kb","+5kb")) +
    theme(legend.position='bottom',plot.title = element_text(size=18),axis.text = element_text(size=15),axis.title = element_text(size=12),legend.title =element_text(size=17),
           legend.text = element_text(size=15),plot.tag.position=c(1.02, 0.5),plot.margin = unit(c(1, 1, 1, 1), "cm"))
  plot_list[[i]] = p
}




tmp <- do.call('ggarrange',c(plot_list[6],plot_list[5],plot_list[3],plot_list[4],plot_list[1],plot_list[2], ncol = 2,nrow=3,common.legend = T,labels='auto'))

pdf(file=paste0(dir_out_figures,"/Figure_S2.pdf"),paper = 'special',onefile = T,width = 11, height = 14,pointsize = 10,colormodel = 'cmyk')
tmp
dev.off()

pdf(file=paste0(dir_out_figures,"/Figure_1B.pdf"),paper = 'special',onefile = T,width = 12, height = 7,pointsize = 10,colormodel = 'cmyk')
Genebody_plot
dev.off()
        

