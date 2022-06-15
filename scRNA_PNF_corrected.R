###########################################################################
###########################################################################
### SINGLE-CELL RNA-SEQ ANALYSIS OF 3 PNF SAMPLES
###########################################################################
###########################################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## SOME OF THE ANALYSES SHOWN MAY HAVE BEEN DISCONTINUED, 
## AND ARE THEREFORE NOT INCLUDED IN THE FDP MANUSCRIPT
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(DropletUtils)
library(scater)
library(scuttle)
library(gridExtra)
library(ggplot2)
library(scran)
library(cowplot)
library(bluster)
library(RColorBrewer)
library(EnsDb.Hsapiens.v75)
library(batchelor)
library(limma)
library(sva)
library(Seurat)
library(ggsignif)
library(scDblFinder)
library(BiocSingular)
library(edgeR)
library(celldex)
library(SingleR)
library(pheatmap)
library(stringr)
library(TSCAN)
library(zellkonverter)
library(velociraptor)

###########################################################################
## LOADING
###########################################################################
sessionInfo()
setwd("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF")
sce_PNF19 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/run_count_scRNA_PNF19/outs/filtered_feature_bc_matrix")
sce_PNF20 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/run_count_scRNA_PNF20/outs/filtered_feature_bc_matrix")
sce_PNF23 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/run_count_scRNA_PNF23/outs/filtered_feature_bc_matrix")

sce_PNF19$batch <- "PNF19"
sce_PNF20$batch <- "PNF20"
sce_PNF23$batch <- "PNF23"

rowData(sce_PNF19) <- rowData(sce_PNF20)
rowData(sce_PNF19) <- rowData(sce_PNF23)

combined <- cbind(sce_PNF19, sce_PNF20, sce_PNF23)

combined$batch <- factor(combined$batch)

###########################################################################
## QUALITY CONTROL
###########################################################################

gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")  #Identifying mitochondrial transcripts

is.mito <- rownames(combined) %in% gns$gene_id
sum(is.mito)

df <- perCellQCMetrics(combined, subsets=list(Mito=is.mito)) #Generate the dataframe


reasons <- perCellQCFilters(df, sub.fields=c("subsets_Mito_percent")) #We make use of this function to identify outliers for our QC metrics
colSums(as.matrix(reasons)) #We obtain the number of cells to discard


combined <- combined[,!reasons$discard]  #We discard low-quality cells from further analysis

###########################################################################
## NORMALIZATION
###########################################################################

set.seed(1101)
lib.sf <- librarySizeFactors(combined)  #compute the library size factors
summary(lib.sf)
hist(log10(lib.sf), xlab="Log10[Size factor]", col='grey80')

combined <- logNormCounts(combined, size.factors = lib.sf)  #normalize the counts
assayNames(combined)

rownames(combined) <- elementMetadata(combined)$Symbol  #rename row names to gene symbols

###########################################################################
## FEATURE SELECTION
###########################################################################

set.seed(11101)
dec.pois.filtered.norm <- modelGeneVarByPoisson(combined)
dec.pois.filtered.norm <- dec.pois.filtered.norm[order(dec.pois.filtered.norm$bio, decreasing=TRUE),] #order
dec.pois.filtered.norm #We will use this for downstream analysis
citation("scran")

chosen <- getTopHVGs(dec.pois.filtered.norm, prop=0.2)  #get the 20% most important genes
str(chosen)

combined.hvg <- combined[chosen,]  #new sce object containing only HVGs
dim(combined.hvg)

###########################################################################
## BATCH CORRECTION
###########################################################################

# assay(combined.hvg, "corrected") <- removeBatchEffect(logcounts(combined.hvg), batch=combined.hvg$batch)
# assay(combined.hvg, "corrected") <- ComBat(logcounts(combined.hvg), batch = combined.hvg$batch)
# combined.hvg <- runPCA(combined.hvg, ncomponents=30, exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
# set.seed(1010)
# combined.hvg <- runTSNE(combined.hvg, dimred="PCA", perplexity=80)
# plotReducedDim(combined.hvg, "TSNE", colour_by = "batch")   # need to correct


set.seed(1001100)
combined.mnn <- fastMNN(logcounts(combined.hvg), batch = combined.hvg$batch)
reducedDim(combined.hvg, "corrected") <- reducedDim(combined.mnn, "corrected")

###########################################################################
## DIMENSIONALITY REDUCTION
###########################################################################

# combined.hvg <- runPCA(combined.hvg, ncomponents=30, exprs_values="corrected", BSPARAM=BiocSingular::ExactParam())

set.seed(1010)
combined.hvg <- runTSNE(combined.hvg, dimred="corrected", perplexity=80)
plotReducedDim(combined.hvg, "TSNE", colour_by = "batch") + theme(axis.title = element_text(size = 30),legend.text = element_text(size = 30), legend.title =  element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))


set.seed(1234)
combined.hvg <- runUMAP(combined.hvg, dimred="corrected")
plotReducedDim(combined.hvg, "UMAP", colour_by = "batch")

############################################################################################
############################################################################################
## ELIMINATING STRESSED CELLS (>0.35) - QUALITY CONTROL - 2nd STEP
############################################################################################
############################################################################################

sce.seurat <- as.Seurat(combined.hvg, counts = "counts", data = "logcounts")

genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/ER_genes.txt", sep = "\t", header = FALSE)[,1]

genes <- genes[genes %in% rownames(combined.hvg)]
genes <- unique(genes)
set.seed(12345)
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")

combined.hvg$signature <- sce.seurat$signature1


combined.hvg$stressed <- combined.hvg$signature > 0.35
sce_unstressed <- combined.hvg[,!combined.hvg$stressed]
combined_unstressed <- combined[,!combined.hvg$stressed]

violin <- data.frame(combined.hvg$batch, combined.hvg$signature)
ggplot(data = violin, mapping = aes(x = combined.hvg.batch, y = combined.hvg.signature)) + geom_violin(aes(fill = combined.hvg.batch)) +
  scale_fill_manual(values = c("#719fcf", "#ffa04d", "#69c15e"))+ geom_signif(comparisons = list(c("PNF19", "PNF20")),  map_signif_level=FALSE, test = "t.test", size = 2, textsize = 10) +
  geom_signif(comparisons = list(c("PNF19", "PNF23")) ,y_position = 1.2, map_signif_level=FALSE, test = "t.test", size = 2, textsize = 10) + labs(x="Sample", y="Score", fill = "Sample") +
  theme_linedraw() + geom_hline(aes(yintercept = 0.35), color = "red", size = 2) + annotate(geom = "text", x=1.5, y=0.4, label ="0.35", color = "red", size = 10) + theme(legend.position = "blank",text = element_text(size = 40), axis.title.x = element_blank()) 
  

plotReducedDim(combined.hvg, "TSNE", colour_by = "stressed") + scale_colour_manual(values = c("grey70", "red")) +  theme_classic() + theme(legend.position = "bottom", legend.justification = "center",axis.title  = element_text(size=30), legend.text = element_text(size=30), legend.title = element_blank(), title = element_text(size = 30))  + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))


cut_df <- data.frame(combined.hvg$batch, combined.hvg$stressed)
ggplot(cut_df, mapping = aes(combined.hvg.batch, fill = combined.hvg.stressed)) + geom_bar(position = position_stack(reverse = TRUE)) + 
  geom_text(aes(label = ..count..), stat = "count", position = position_stack(reverse=TRUE), size = 15) + labs(x = "Sample", y= "Number of cells", fill = "Stressed") + scale_fill_manual(values = c("chartreuse3", "red3"))+ theme_classic() + theme(text = element_text(size = 40), axis.title.x = element_blank(), legend.title = element_blank(), legend.position = "bottom", legend.justification = "center")


set.seed(100)
dbl.dens <- computeDoubletDensity(sce_unstressed, d=ncol(reducedDim(sce_unstressed, "corrected")))
summary(dbl.dens)

sce_unstressed$DoubletScore <- dbl.dens
set.seed(1010)
plotTSNE(sce_unstressed, colour_by="DoubletScore") + theme(text = element_text(size = 30), legend.key.size = unit(1, "cm"))

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call", p=0.01)

summary(dbl.calls)
plotColData(sce_unstressed, x="batch", y="DoubletScore", colour_by=I(dbl.calls)) + theme(axis.text  = element_text(size = 20), text = element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))

sce_unstressed <- sce_unstressed[,dbl.calls=="singlet"]
combined_unstressed <- combined_unstressed[,dbl.calls=="singlet"]
###########################################################################
## BATCH CORRECTION
###########################################################################

# assay(sce_unstressed.hvg, "corrected") <- removeBatchEffect(logcounts(sce_unstressed.hvg), batch=sce_unstressed.hvg$batch)
# assay(sce_unstressed.hvg, "corrected") <- ComBat(logcounts(sce_unstressed.hvg), batch = sce_unstressed$batch)
# sce_unstressed <- runPCA(sce_unstressed, ncomponents=30, exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
# set.seed(1010)
# sce_unstressed <- runTSNE(sce_unstressed, dimred="PCA", perplexity=80)
# plotReducedDim(sce_unstressed, "TSNE", colour_by = "batch")  # need to correct

?fastMNN
set.seed(1001100)
sce_unstressed.mnn <- fastMNN(logcounts(sce_unstressed), batch = sce_unstressed$batch)
reducedDim(sce_unstressed, "corrected") <- reducedDim(sce_unstressed.mnn, "corrected")

###########################################################################
## DIMENSIONALITY REDUCTION
###########################################################################

# sce_unstressed <- runPCA(sce_unstressed, ncomponents=30, exprs_values="corrected", BSPARAM=BiocSingular::ExactParam())

set.seed(1010)  #001001101101
sce_unstressed <- runTSNE(sce_unstressed, dimred="corrected", perplexity=120)
plotReducedDim(sce_unstressed, "TSNE", colour_by = "batch") + theme(axis.title = element_text(size = 30),legend.text = element_text(size = 30), legend.title =  element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))


set.seed(1234)
sce_unstressed <- runUMAP(sce_unstressed, dimred="corrected")
plotReducedDim(sce_unstressed, "UMAP", colour_by = "batch")

###########################################################################
## CLUSTERING
###########################################################################

set.seed(11111)
clust2s <- clusterCells(sce_unstressed, use.dimred="corrected", BLUSPARAM=TwoStepParam(first=KmeansParam(centers=2500),
                                                                                     second=NNGraphParam(k=35, type = "jaccard", cluster.fun = "infomap")))

table(clust2s)


colLabels(sce_unstressed) <- clust2s

p1 <- plotReducedDim(sce_unstressed, "TSNE", colour_by="label", text_by = "label", text_size = 10) + theme(axis.title = element_text(size = 30),legend.text = element_text(size = 30), legend.title =  element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))
p1

plotReducedDim(sce_unstressed, "UMAP", colour_by="label", text_by = "label")

###########################################################################
## EXPRESSION PLOTS
###########################################################################
reducedDim(combined_unstressed, "TSNE") <- reducedDim(sce_unstressed, "TSNE")
combined_unstressed$label <- sce_unstressed$label

gene <- "S100B"     #plug in any gene
p2 <- plotReducedDim(combined_unstressed, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =6) +
  labs(title = paste0(gene, " expression per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()

p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)

genes <- c("STMN1", "EHBP1", "GAS7")

plots <- list()
i <- 1
for (gene in genes){
  p2 <- plotReducedDim(sce_unstressed, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =10) +
     theme_classic() + theme(axis.title  = element_text(size=20), legend.text = element_text(size=30), legend.key.size = unit(1.5, "cm"), legend.title = element_text(size=30))
  
  p2 <- p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)
  plots[[i]] <- p2
  i <- i+1
}
grid.arrange(grobs = plots, ncol = 3)

grid.arrange(p1,p2, ncol = 3)
###########################################################################
## SIGNATURE SCORING
###########################################################################


sce.seurat <- as.Seurat(sce_unstressed, counts = "counts", data = "logcounts")

genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/up_day30.txt", sep = "\t", header = FALSE)[,1]

genes <- genes[genes %in% rownames(sce_unstressed)]
genes <- unique(genes)
genes <- genes[1:length(genes)-1]
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")


sce_unstressed$signature <- sce.seurat$signature1
signature_plot <- plotReducedDim(sce_unstressed, "TSNE", colour_by = "signature", text_by = "label", text_size = 10)
p4 <- signature_plot + scale_colour_gradient2(low = "blue", mid = "#F7F7F7", high = "red",  midpoint = 0) 
p3 <- p4 + labs(color = "Scores") + theme(axis.title  = element_text(size=20), legend.text = element_text(size=30), legend.key.size = unit(1.5, "cm"), legend.title = element_text(size=30))
p1+p2+p3+p4
grid.arrange(p1,p2,p3, nrow=1)
# GO:0001837 and GO:0007409 are interesting GO:0048762  GO:0048699  --> cluster 9 goes both ways
###########################################################################
## DE ANALYSIS
###########################################################################

summed <- aggregateAcrossCells(sce_unstressed, 
                               id=colData(sce_unstressed)[,c("celltype.main", "label")])
summed





y <- DGEList(counts(summed), samples=colData(summed))
y

discarded <- summed$ncells < 10
y <- y[,!discarded]
summary(discarded)


y <- calcNormFactors(y)

y


y$samples$group <- "others"
# y$samples$celltype[y$samples$label %in% c(13,11,10,9,7)] <- "SchwannControl"
# y$samples$celltype[y$samples$label %in% c(2) & y$samples$batch == "PNF19"] <- "SchwannCase"
# y$samples$celltype[y$samples$label %in% c(16) & ] <- "double1"
# y$samples$celltype[y$samples$label %in% c(16)] <- "double1"
y$samples$group[y$samples$label %in% c(1) ] <- "case"
y$samples$group[y$samples$label %in% c(2) ] <- "control"

y$samples$group


label <- y$samples$label
batch <- y$samples$batch
group <- y$samples$group


plotMDS(y, col = as.numeric(batch))

mm <- model.matrix(~ 0 + group)

mm
voom <- voom(y, mm, plot = T)


fit <- lmFit(voom, mm)
head(coef(fit))


contr.matrix <- makeContrasts(
  SCBvsSCA = groupcase-groupcontrol,
  levels = colnames(coef(fit)))


contr.matrix

tmp <- contrasts.fit(fit, contr.matrix)

tmp <- eBayes(tmp)

top.table <- topTable(tmp, sort.by = "logFC", n = Inf)
top.table <- top.table[top.table$logFC>0,]
head(top.table, 50)

sink("test_genes.txt")
for (i in 1:50){
  cat(rownames(head(top.table, 50))[i])
  cat("\n")
}
sink()
rownames(head(top.table, 20))


###########################################################################
## MARKER DETECTION
###########################################################################
sce_unstressed$cdh <- 2
sce_unstressed[,assay(sce_unstressed, "counts")["CDH19",]>0]$cdh <- 1

assay(sce_unstressed, "counts")["CDH19",]
################################################################
m.out <- scoreMarkers(sce_unstressed, sce_unstressed$cdh, block=sce_unstressed$batch)

chosen <- m.out[["1"]]
detect.only <- chosen[,grepl("logFC.detected", colnames(chosen))]
ordered <- detect.only[order(detect.only$mean.logFC.detected,decreasing=TRUE),]

ordered[1:10,]
rownames(ordered[1:100,])
plotGroupedHeatmap(sce_unstressed, features=rownames(ordered[1:20,]), group="label", 
                   center=TRUE, zlim=c(-3, 3))


plotExpression(combined.hvg, features=head(rownames(ordered)), 
               x="label", colour_by="batch")


  sink("test_genes.txt")
for (i in 1:30){
  cat(rownames(ordered[1:30,])[i])
  cat("\n")
}
sink()


###########################################################################
## CELL TYPE ASSIGNMENT
###########################################################################


ref <- HumanPrimaryCellAtlasData()

pred <- SingleR(test=sce_unstressed, ref=ref, labels=ref$label.fine)
table(pred$labels)

tab <- table(Assigned=pred$pruned.labels, Cluster=colLabels(sce_unstressed))
# 
# c <- 12
# sort((tab[,c] / sum(tab[,c]))*100)  # obtain the percentages for each type
# # sum(tab[,c])  # number of cells in the cluster

pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101)) # Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.

sce_unstressed$celltype <- pred$first.labels

table(sce_unstressed$celltype[sce_unstressed$label ==9])



sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Fibroblasts.*", "Fibroblasts")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Endothelial_cells.*", "Endothelial cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "B_cell.*", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "T_cell.*", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "NK.*", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Tissue_stem_cells.*", "Tissue stem cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "iPS_cells.*", "iPS cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Macrophage.*", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Monocyte.*", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "DC.*", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Smooth_muscle.*", "Smooth muscle cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "MEP", "Immune cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Chondrocytes.*", "Chondrocytes")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Neuron.*", "Schwann cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Neuro.*", "Neuroepithelial cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Epi.*", "Epithelial cells")
sce_unstressed$celltype <- str_replace_all(sce_unstressed$celltype, "Keratinocytes.*", "Epithelial cells")
types <- rownames(table(sce_unstressed$celltype)[table(sce_unstressed$celltype)>50])

sce_unstressed$celltype[!sce_unstressed$celltype %in% types] <- "others"

sce_unstressed$celltype.main <- sce_unstressed$celltype
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Schwann cells", "Schwann cell")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Neuroepithelial cells", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Immune cells", "Immune cells")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Fibroblasts", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Tissue stem cells", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Smooth muscle cells", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Chondrocytes", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Epithelial cells", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Astrocyte.*", "Schwann cell")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Fibroblasts", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "Osteoblasts", "Mesenchymal")
sce_unstressed$celltype.main <- str_replace_all(sce_unstressed$celltype.main, "MSC", "Mesenchymal")


plotReducedDim(sce_unstressed, "TSNE", colour_by="celltype", text_by = "label", text_size = 10) + scale_color_manual(name = "Cell Type", values = c("Immune cells" = "red", "Schwann cells" = "forestgreen", "Endothelial cells" = "orange",
                                                                                                                             "Fibroblasts" = "firebrick4", "Smooth muscle cells" = "peru", "Tissue stem cells" = "lightpink", 
                                                                                                                             "MSC" = "purple", "iPS cells" = "turquoise", "Astrocyte:Embryonic_stem_cell-derived" = "blue4", 
                                                                                                                             "others" = "azure4", "Epithelial cells" = "yellow", "Neuroepithelial cells" = "green", "Chondrocytes" = "magenta")) + theme(legend.position = "bottom", legend.justification = "center", legend.text = element_text(size = 12),axis.title  = element_text(size = 20), legend.title = element_text(size = 15)) +
   guides(colour = guide_legend(override.aes = list(shape = 15, size = 6))) 



plotReducedDim(sce_unstressed, "TSNE", colour_by="celltype.main", text_by = "label", text_size = 10) + scale_color_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey")) +
  theme(text = element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10))) + theme(legend.position = "blank")


celltype_prop <- prop.table(as.array(table(sce_unstressed$celltype, sce_unstressed$label)), margin = 2)*100

df <- data.frame(sce_unstressed$batch, sce_unstressed$celltype.main)

ggplot(df, mapping = aes(x=sce_unstressed.batch, fill=sce_unstressed.celltype.main)) + geom_bar(position = "fill") + labs(x = "Sample of origin", y= "Proportion") +
  theme_minimal() + theme(legend.position = "bottom", text = element_text(size = 30))  + scale_fill_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey")) + theme(legend.position = "blank", axis.title.x = element_blank())
                                                                        
###########################################################################
## DOUBLET DETECTION
###########################################################################

set.seed(100)
dbl.dens <- computeDoubletDensity(sce_unstressed, d=ncol(reducedDim(sce_unstressed, "corrected")))
summary(dbl.dens)

sce_unstressed$DoubletScore <- dbl.dens
set.seed(1010)
plotTSNE(sce_unstressed, colour_by="DoubletScore", text_by = "label", text_colour = "red")

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call", p=0.01)

summary(dbl.calls)
plotColData(sce_unstressed, x="label", y="DoubletScore", colour_by=I(dbl.calls))



###########################################################################
## TRAJECTORY ANALYSIS
###########################################################################

sce.schwann <- sce_unstressed[,sce_unstressed$label %in% c(9, 11, 3, 5, 1, 4)]
plotReducedDim(sce.schwann, "TSNE", colour_by = "label", text_by = "label")


by.cluster <- aggregateAcrossCells(sce.schwann, ids=colLabels(sce.schwann))

centroids <- reducedDim(by.cluster, "corrected")


mst <- createClusterMST(centroids, clusters=NULL)  #generate the mst
mst

line.data <- reportEdges(by.cluster, mst=mst, clusters=NULL, use.dimred="TSNE")

plotTSNE(sce.schwann, colour_by="label") + 
  geom_line(data=line.data, mapping=aes(x=dim1, y=dim2, group=edge))

map.tscan <- mapCellsToEdges(sce.schwann, mst=mst, use.dimred="corrected")


tscan.pseudo <- orderCells(map.tscan, mst, start = 9) #forcing 2 to be the start gives the same result
head(tscan.pseudo)

common.pseudo <- averagePseudotime(tscan.pseudo) 
plotTSNE(sce.schwann, colour_by=I(common.pseudo), 
         text_by="label", text_colour="red") +
  geom_line(data=line.data, mapping=aes(x=dim1, y=dim2, group=edge)) +
  labs(title = "Trajectory plot for Schwann component")

# CHARACTERIZING TRAJECTORIES

pseudo <- testPseudotime(sce.schwann, pseudotime=tscan.pseudo[,2])[[1]]

sorted <- pseudo[order(pseudo$logFC, decreasing =TRUE),] ## CDK2AP2, S100A4, GADD45B 

sorted[1:10,]
gene <- "ID2"
p2 <- plotReducedDim(sce.schwann, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =6) +
  labs(title = paste0(gene, " expression per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()

p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)

sce.schwann$TSCAN.first <- pathStat(tscan.pseudo)[,1]
sce.schwann$TSCAN.second <- pathStat(tscan.pseudo)[,2]

best <- head(rownames(sorted), 4)
plotExpression(sce.schwann, features=best,
               x="TSCAN.second", colour_by="label")

sorted2 <- pseudo[order(pseudo$logFC, decreasing =FALSE),]
sorted2[1:10,]
best2 <- head(rownames(sorted2), 4)
plotExpression(sce.schwann, features=best2,
               x="TSCAN.second", colour_by="label")


on.first.path <- !is.na(sce.schwann$TSCAN.first)
plotHeatmap(sce.schwann[,on.first.path], order_columns_by="TSCAN.first", 
            colour_columns_by="label", features=head(rownames(sorted2), 50),
            center=TRUE)

head(rownames(sorted2), 50)

sink("test_genes.txt")
for (i in 1:50){
  cat(head(rownames(sorted2), 50)[i])
  cat("\n")
}
sink()


###########################################################################
## RNA VELOCITY
###########################################################################


spliced_19 <- readH5AD("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/data/PNF19/counts_filtered/PNF19.h5ad", reader = "R")
spliced_20 <- readH5AD("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/data/PNF20/counts_filtered/PNF20.h5ad", reader = "R")
spliced_23 <- readH5AD("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/data/PNF23/counts_filtered/PNF23.h5ad", reader = "R")

spliced_19$barcode <- paste0(spliced_19$bcs, "-1")
spliced_20$barcode <- paste0(spliced_20$bcs, "-1")
spliced_23$barcode <- paste0(spliced_23$bcs, "-1")

rownames(spliced_19) <- rowData(spliced_19)[,1]
rownames(spliced_20) <- rowData(spliced_20)[,1]
rownames(spliced_23) <- rowData(spliced_23)[,1]

# spliced_19 <- spliced_19[,spliced_19$barcode %in% sce_PNF19$Barcode]
# sce_PNF19 <- sce_PNF19[,sce_PNF19$Barcode %in% spliced_19$barcode]

# assay(sce_PNF19, "spliced") <- assay(spliced_19, "spliced")


spliced_19$batch <- "PNF19"
spliced_20$batch <- "PNF20"
spliced_23$batch <- "PNF23"

rowData(spliced_19) <- rowData(spliced_20)
rowData(spliced_19) <- rowData(spliced_23)

spliced <- cbind(spliced_19, spliced_20, spliced_23)

spliced$batch <- factor(spliced$batch)

colnames(sce_unstressed) <- paste0(sce_unstressed$Barcode, ".", sce_unstressed$batch)
colnames(spliced) <- paste0(spliced$barcode, ".", spliced$batch)

spliced <- spliced[,colnames(spliced) %in% colnames(sce_unstressed)]
spliced <- spliced[rownames(spliced) %in% rowData(sce_unstressed)[,1],]


reducedDim(spliced, "corrected") <- reducedDim(sce_unstressed[,colnames(sce_unstressed) %in% colnames(spliced)], "corrected")

reducedDim(spliced, "TSNE") <- reducedDim(sce_unstressed[,colnames(sce_unstressed) %in% colnames(spliced)], "TSNE")
plotReducedDim(spliced, "TSNE", colour_by = "batch")

spliced$label <- sce_unstressed[,colnames(sce_unstressed) %in% colnames(spliced)]$label
plotReducedDim(spliced, "TSNE", colour_by = "label", text_by = "label")


set.seed(1101)
velo.out <- scvelo(spliced, assay.X="spliced", use.dimred="corrected")
velo.out

spliced$pseudotime <- velo.out$velocity_pseudotime

# Also embedding the velocity vectors, for some verisimilitude.
embedded <- embedVelocity(reducedDim(spliced, "TSNE"), velo.out)
grid.df <- gridVectors(reducedDim(spliced, "TSNE"), embedded, resolution=50)

plotTSNE(spliced, colour_by="pseudotime", text_by = "label", text_colour = "red", point_alpha=0.15) +
  geom_segment(data=grid.df, 
               mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), 
                 arrow=arrow(length=unit(0.05, "inches"), type="closed")) + scale_color_viridis_c(option = "plasma") + labs(colour = "Pseudotime")


spliced$celltype.main <- sce_unstressed[,colnames(sce_unstressed) %in% colnames(spliced)]$celltype.main
plotTSNE(spliced, colour_by="celltype.main", text_by = "label", text_colour = "red", point_alpha=0.5, text_size = 10) +
  geom_segment(data=grid.df, 
               mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), 
               arrow=arrow(length=unit(0.05, "inches"), type="closed")) + scale_color_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey")) + theme(legend.position = "right", text = element_text(size = 30))

#########################################################################################
#########################################################################################

reducedDim(combined_unstressed, "TSNE") <- reducedDim(sce_unstressed, "TSNE")
combined_unstressed$label <- sce_unstressed$label


df <- data.frame(counts = assay(combined_unstressed, "counts")["NF1",], cluster = sce_unstressed$celltype.main)
ggplot(df, mapping = aes(x = cluster, y = counts)) + geom_boxplot()
means <- aggregate(counts ~  cluster, df, mean)
means$counts <- round(means$counts, digits = 2)
ggplot(df, mapping = aes(x = cluster, y = counts, fill = cluster)) + geom_boxplot() + coord_cartesian(ylim = c(0,1)) + labs(x = "CDH19 and PRRX1", y = "Number of reads", fill = "CDH19 and PRRX1") + stat_summary(fun = "mean")+
  geom_text(data = means, aes(label = paste0("Mean = ", counts), y = counts - 1000), size = 7) 
