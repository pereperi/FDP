###########################################################################
###########################################################################
### SINGLE-CELL RNA-SEQ ANALYSIS OF 3D NEUROFIBROMA MODELS + PRIMARY TUMORS
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
library(scDblFinder)
library(BiocSingular)
library(Seurat)
library(ggsignif)
library(celldex)
library(SingleR)
library(pheatmap)
library(stringr)
library(edgeR)
###########################################################################
## LOADING
###########################################################################
sessionInfo()
setwd("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF")
sce_PNF19 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/run_count_scRNA_PNF19/outs/filtered_feature_bc_matrix")
sce_PNF20 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/run_count_scRNA_PNF20/outs/filtered_feature_bc_matrix")
sce_PNF23 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/preprocessed_data/run_count_scRNA_PNF23/outs/filtered_feature_bc_matrix")
sce_3MM <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_spheroids/preprocessed_data/run_count_scRNA_3MM/outs/filtered_feature_bc_matrix")
sce_FIPS <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_spheroids/preprocessed_data/run_count_scRNA_FIPS/outs/filtered_feature_bc_matrix")
sce_D12 <- read10xCounts("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_spheroids/preprocessed_data/run_count_scRNA_D12/outs/filtered_feature_bc_matrix")



sce_PNF19$batch <- "PNF19"
sce_PNF20$batch <- "PNF20"
sce_PNF23$batch <- "PNF23"
sce_3MM$batch <- "3MM"
sce_D12$batch <- "D12"
sce_FIPS$batch <- "FIPS"

rowData(sce_PNF19) <- rowData(sce_PNF20)
rowData(sce_PNF19) <- rowData(sce_PNF23)
rowData(sce_PNF19) <- rowData(sce_3MM)
rowData(sce_PNF19) <- rowData(sce_D12)
rowData(sce_PNF19) <- rowData(sce_FIPS)



###########################################################################
## QUALITY CONTROL
###########################################################################

gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")  #Identifying mitochondrial transcripts

all.sce <- list(
  sce_D12 = sce_D12,
  sce_FIPS = sce_FIPS,
  sce_PNF19 = sce_PNF19,
  sce_PNF20 = sce_PNF20,
  sce_PNF23 = sce_PNF23
)
stats <- reasons <- list()
for (n in names(all.sce)) {
  current <- all.sce[[n]]
  is.mito <- rownames(current) %in% gns$gene_id
  stats[[n]] <- perCellQCMetrics(current, subsets=list(Mito=is.mito))
  reasons[[n]] <- perCellQCFilters(stats[[n]], sum.field = "sum", sub.fields=c("subsets_Mito_percent"))
  all.sce[[n]] <- current[,!reasons[[n]]$discard]
}

is.mito <- rownames(sce_3MM) %in% gns$gene_id
sum(is.mito)

df <- perCellQCMetrics(sce_3MM, subsets=list(Mito=is.mito)) #Generate the dataframe

reasons <- perCellQCFilters(df, sub.fields=c("subsets_Mito_percent")) #We make use of this function to identify outliers for our QC metrics
reasons$low_lib_size <- df$sum < 2000 #we manually set the threshold given it presents a bimodal distribution
reasons$discard <- reasons$low_lib_size | reasons$low_n_features | reasons$high_subsets_Mito_percent
colSums(as.matrix(reasons)) #We obtain the number of cells to discard

sce_3MM <- sce_3MM[,!reasons$discard]  #We discard low-quality cells from further analysis

merge <- cbind(all.sce$sce_PNF19, all.sce$sce_PNF20, all.sce$sce_PNF23, all.sce$sce_D12, all.sce$sce_FIPS, sce_3MM)

merge$batch <- factor(merge$batch)


###########################################################################
## NORMALIZATION
###########################################################################

set.seed(1101)
lib.sf <- librarySizeFactors(merge)  #compute the library size factors
summary(lib.sf)
hist(log10(lib.sf), xlab="Log10[Size factor]", col='grey80')

merge <- logNormCounts(merge, size.factors = lib.sf)  #normalize the counts
assayNames(merge)

rownames(merge) <- elementMetadata(merge)$Symbol  #rename row names to gene symbols

###########################################################################
## FEATURE SELECTION
###########################################################################

set.seed(11101)
dec.pois.filtered.norm <- modelGeneVarByPoisson(merge)
dec.pois.filtered.norm <- dec.pois.filtered.norm[order(dec.pois.filtered.norm$bio, decreasing=TRUE),] #order
dec.pois.filtered.norm #We will use this for downstream analysis

chosen <- getTopHVGs(dec.pois.filtered.norm, prop=0.2)  #get the 20% most important genes
str(chosen)

merge.hvg <- merge[chosen,]  #new sce object containing only HVGs
dim(merge.hvg)

###########################################################################
## BATCH CORRECTION
###########################################################################

# assay(merge.hvg, "corrected") <- removeBatchEffect(logcounts(merge.hvg), batch=merge.hvg$batch)
# assay(merge.hvg, "corrected") <- ComBat(logcounts(merge.hvg), batch = merge.hvg$batch)
# merge.hvg <- runPCA(merge.hvg, ncomponents=30, exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
# set.seed(1010)
# merge.hvg <- runTSNE(merge.hvg, dimred="PCA", perplexity=80)
# plotReducedDim(merge.hvg, "TSNE", colour_by = "batch")  # need to correct


set.seed(1001100)
merge.mnn <- fastMNN(logcounts(merge.hvg), batch = merge.hvg$batch)
reducedDim(merge.hvg, "corrected") <- reducedDim(merge.mnn, "corrected")

###########################################################################
## DIMENSIONALITY REDUCTION
###########################################################################

set.seed(1010)
merge.hvg <- runTSNE(merge.hvg, dimred="corrected", perplexity=80)
plotReducedDim(merge.hvg, "TSNE", colour_by = "batch")


set.seed(1234)
merge.hvg <- runUMAP(merge.hvg, dimred="corrected")
plotReducedDim(merge.hvg, "UMAP", colour_by = "batch")

###########################################################################
## QUALITY CONTROL - 2nd STEP
###########################################################################


set.seed(100)
dbl.dens <- computeDoubletDensity(merge.hvg, d=ncol(reducedDim(merge.hvg, "corrected")))
summary(dbl.dens)

merge.hvg$DoubletScore <- dbl.dens
set.seed(1010)
plotTSNE(merge.hvg, colour_by="DoubletScore", text_colour = "red")

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call", p=0.01)

summary(dbl.calls)
plotColData(merge.hvg, x="batch", y="DoubletScore", colour_by=I(dbl.calls))

merge.hvg <- merge.hvg[,dbl.calls=="singlet"]



sce.seurat <- as.Seurat(merge.hvg, counts = "counts", data = "logcounts")

genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/ER_genes.txt", sep = "\t", header = FALSE)[,1]

genes <- genes[genes %in% rownames(merge.hvg)]
genes <- unique(genes)
set.seed(12345)
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")

merge.hvg$signature <- sce.seurat$signature1
signature_plot <- plotReducedDim(merge.hvg, "TSNE", colour_by = "signature")
p4 <- signature_plot + scale_colour_gradient2(low = "blue", mid = "#F7F7F7", high = "red",  midpoint = 0) 
p4 <- p4 + labs(color = "Scores")
p4

merge.hvg$stressed <- merge.hvg$signature > 0.35
merge.hvg_stressed <- merge.hvg
merge.hvg <- merge.hvg[,!merge.hvg$stressed]


violin <- data.frame(merge.hvg_stressed$batch, merge.hvg_stressed$signature)
ggplot(data = violin, mapping = aes(x = merge.hvg_stressed.batch, y = merge.hvg_stressed.signature)) + geom_violin(aes(fill = merge.hvg_stressed.batch)) + labs(x="Sample", y="Score of response to ER stress", fill = "Sample") +
  theme_linedraw() + geom_hline(aes(yintercept = 0.35), color = "red") + annotate(geom = "text", x=1.5, y=0.4, label ="Threshold = 0.35", color = "red")

plotReducedDim(merge.hvg_stressed, "TSNE", colour_by = "stressed") + scale_colour_manual(values = c("grey70", "red")) + labs(color = "Stressed (Score>0.35)")


cut_df <- data.frame(merge.hvg_stressed$batch, merge.hvg_stressed$stressed)
ggplot(cut_df, mapping = aes(merge.hvg_stressed.batch, fill = merge.hvg_stressed.stressed)) + geom_bar(position = position_stack(reverse = TRUE)) + 
  geom_text(aes(label = ..count..), stat = "count", position = position_stack(reverse=TRUE), ) + labs(x = "Sample", y= "Number of cells", fill = "Stressed (Score>0.35)") + scale_fill_manual(values = c("chartreuse3", "red3"))


###########################################################################
## BATCH CORRECTION
###########################################################################

# assay(merge.hvg, "corrected") <- removeBatchEffect(logcounts(merge.hvg), batch=merge.hvg$batch)
# assay(merge.hvg, "corrected") <- ComBat(logcounts(merge.hvg), batch = merge.hvg$batch)
# merge.hvg <- runPCA(merge.hvg, ncomponents=30, exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
# set.seed(1010)
# merge.hvg <- runTSNE(merge.hvg, dimred="PCA", perplexity=80)
# plotReducedDim(merge.hvg, "TSNE", colour_by = "batch")  # need to correct


set.seed(1001100)
merge.mnn <- fastMNN(logcounts(merge.hvg), batch = merge.hvg$batch)
reducedDim(merge.hvg, "corrected") <- reducedDim(merge.mnn, "corrected")

###########################################################################
## DIMENSIONALITY REDUCTION
###########################################################################

set.seed(101010)
merge.hvg <- runTSNE(merge.hvg, dimred="corrected", perplexity=120)
plotReducedDim(merge.hvg, "TSNE", colour_by = "batch")

merge.hvg$type[merge.hvg$batch %in% c("3MM", "D12", "FIPS")] <- "Spheres" 
merge.hvg$type[merge.hvg$batch %in% c("PNF19", "PNF20", "PNF23")] <- "Primary Tumor"
plotReducedDim(merge.hvg, "TSNE", colour_by = "type") + theme(axis.title = element_text(size = 30),legend.text = element_text(size = 30), legend.title =  element_blank(), legend.position = "bottom", legend.justification = "center") + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))

set.seed(1234)
merge.hvg <- runUMAP(merge.hvg, dimred="corrected")
plotReducedDim(merge.hvg, "TSNE", colour_by = "batch") + theme(axis.title = element_text(size = 30),legend.text = element_text(size = 30), legend.title =  element_blank(), legend.position = "bottom") + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10))) 


###########################################################################
## CLUSTERING
###########################################################################

set.seed(11111)
clust2s <- clusterCells(merge.hvg, use.dimred="corrected", BLUSPARAM=TwoStepParam(first=KmeansParam(centers=2500),
                                                                                       second=NNGraphParam(k=35, type = "jaccard", cluster.fun = "infomap")))

table(clust2s)


colLabels(merge.hvg) <- clust2s

p1 <- plotReducedDim(merge.hvg, "TSNE", colour_by="label", text_by = "label", text_size = 10) + theme(axis.title = element_text(size = 30),legend.text = element_text(size = 30), legend.title =  element_blank(), legend.position = "blank") + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10)))
p1

plotReducedDim(merge.hvg, "UMAP", colour_by="type", text_by = "label")

###########################################################################
## EXPRESSION PLOTS
###########################################################################

  gene <- "PTPRN"     #plug in any gene
p2 <- plotReducedDim(merge.hvg, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =6) +
  labs(title = paste0(gene, " expression per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()

p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)


genes <- c("SOX10","CDH19","S100B", "TWIST1", "PRRX1")

plots <- list()
i <- 1
for (gene in genes){
  p2 <- plotReducedDim(merge.hvg, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts") +
    theme_classic() + theme(axis.title  = element_text(size=20), legend.text = element_text(size=30), legend.key.size = unit(1.5, "cm"), legend.title = element_text(size=30))
  
  p2 <- p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)
  plots[[i]] <- p2
  i <- i+1
}
grid.arrange(grobs = plots, ncol = 3)

###########################################################################
## CELL TYPE ASSIGNMENT
###########################################################################


ref <- HumanPrimaryCellAtlasData()


pred <- SingleR(test=merge.hvg, ref=ref, labels=ref$label.fine)
table(pred$labels)

tab <- table(Assigned=pred$pruned.labels, Cluster=colLabels(merge.hvg))
# 
# c <- 12
# sort((tab[,c] / sum(tab[,c]))*100)  # obtain the percentages for each type
# # sum(tab[,c])  # number of cells in the cluster


pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101)) # Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.

merge.hvg$celltype <- pred$labels




merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Fibroblasts.*", "Fibroblasts")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Endothelial_cells.*", "Endothelial cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "B_cell.*", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "T_cell.*", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "NK.*", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Tissue_stem_cells.*", "Tissue stem cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "iPS_cells.*", "iPS cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Macrophage.*", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Monocyte.*", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "DC.*", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Smooth_muscle.*", "Smooth muscle cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "MEP", "Immune cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Chondrocytes.*", "Chondrocytes")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Neuron.*", "Schwann cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Neuro.*", "Neuroepithelial cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Epi.*", "Epithelial cells")
merge.hvg$celltype <- str_replace_all(merge.hvg$celltype, "Keratinocytes.*", "Epithelial cells")
types <- rownames(table(merge.hvg$celltype)[table(merge.hvg$celltype)>50])

merge.hvg$celltype[!merge.hvg$celltype %in% types] <- "others"

merge.hvg$celltype.main <- merge.hvg$celltype
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Schwann cells", "Schwann cell")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Neuroepithelial cells", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Immune cells", "Immune cells")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Fibroblasts", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Endothelial cells", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Tissue stem cells", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Smooth muscle cells", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Chondrocytes", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Epithelial cells", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Astrocyte.*", "Schwann cell")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Fibroblasts", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "Osteoblasts", "Mesenchymal")
merge.hvg$celltype.main <- str_replace_all(merge.hvg$celltype.main, "MSC", "Mesenchymal")


plotReducedDim(merge.hvg, "TSNE", colour_by="celltype", text_by = "label") + scale_color_manual(name = "Cell Type", values = c("Immune cells" = "red", "Schwann cells" = "forestgreen", "Endothelial cells" = "orange",
                                                                                                                                    "Fibroblasts" = "firebrick4", "Smooth muscle cells" = "peru", "Tissue stem cells" = "lightpink", 
                                                                                                                                    "MSC" = "purple", "iPS cells" = "turquoise", "Astrocyte:Embryonic_stem_cell-derived" = "blue4", 
                                                                                                                                    "others" = "azure4", "Epithelial cells" = "yellow", "Neuroepithelial cells" = "green", "Chondrocytes" = "magenta"))


plotReducedDim(merge.hvg, "TSNE", colour_by="celltype.main") + scale_color_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey")) +
  theme(text = element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10))) + theme(legend.position = "blank")                                               

celltype_prop <- prop.table(as.array(table(merge.hvg$celltype, merge.hvg$label)), margin = 2)*100

df <- data.frame(merge.hvg$label, merge.hvg$celltype)

ggplot(df, mapping = aes(x=merge.hvg.label, fill=merge.hvg.celltype)) + geom_bar(position = "fill") + labs(x = "Cluster", y= "Proportion") +
  theme_minimal() + theme(legend.position = "bottom")  + scale_fill_manual(name = "Cell Type", values = c("Immune cells" = "red", "Schwann cells" = "forestgreen", "Endothelial cells" = "orange", "Fibroblasts" = "firebrick4", "Smooth muscle cells" = "peru", "Tissue stem cells" = "lightpink", 
                                                                                                          "MSC" = "purple", "iPS cells" = "turquoise", "Astrocyte:Embryonic_stem_cell-derived" = "blue4", "others" = "azure4", "Epithelial cells" = "yellow", "Neuroepithelial cells" = "green", "Chondrocytes" = "magenta"))

###########################################################################
## MARKER DETECTION
###########################################################################

m.out <- scoreMarkers(merge.hvg, colLabels(merge.hvg), block=merge.hvg$batch)

chosen <- m.out[["11"]]
detect.only <- chosen[,grepl("logFC.detected", colnames(chosen))]
ordered <- detect.only[order(detect.only$mean.logFC.detected,decreasing=TRUE),]

ordered[1:10,]
rownames(ordered[1:30,])
plotGroupedHeatmap(merge.hvg, features=rownames(ordered[1:20,]), group="label", 
                   center=TRUE, zlim=c(-3, 3))


plotExpression(merge.hvg, features=head(rownames(ordered)), 
               x="label", colour_by="batch")


sink("test_genes.txt")
for (i in 1:30){
  cat(rownames(ordered[1:30,])[i])
  cat("\n")
}
sink()                                                                                                         


###########################################################################
## SIGNATURE SCORING
###########################################################################


sce.seurat <- as.Seurat(merge.hvg, counts = "counts", data = "logcounts")

genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/up_fb.txt", sep = "\t", header = FALSE)[,1]

genes <- genes[genes %in% rownames(merge.hvg)]
genes <- unique(genes)
genes <- genes[1:length(genes)-1]
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")


merge.hvg$signature <- sce.seurat$signature1
signature_plot <- plotReducedDim(merge.hvg[,merge.hvg$signature < 1.5], "TSNE", colour_by = "signature")
p4 <- signature_plot + scale_colour_gradient2(low = "blue", mid = "#F7F7F7", high = "red",  midpoint = 0) 
p4 <- p4 + labs(color = "Scores") + theme(axis.title  = element_text(size=20), legend.text = element_text(size=30), legend.key.size = unit(1.5, "cm"), legend.title = element_text(size=30))
p4                                                                                               

grid.arrange(p1,p2,p3, nrow=1)
###########################################################################
## DE ANALYSIS
###########################################################################

summed <- aggregateAcrossCells(merge.hvg, 
                               id=colData(merge.hvg)[,c("batch", "label")])
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
y$samples$group[y$samples$label %in% c(11)] <- "case"
y$samples$group[y$samples$label %in% c(9)] <- "control"

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

top.table <- topTable(tmp, sort.by = "P", n = Inf)
top.table <- top.table[top.table$logFC>0,]
head(top.table, 50)

sink("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/test_genes.txt")
for (i in 1:50){
  cat(rownames(head(top.table, 50))[i])
  cat("\n")
}
sink()
rownames(head(top.table, 20))


###########################################################################
## CELL CYCLE (TAKES LONG)
###########################################################################

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
rownames(merge) <- elementMetadata(merge)$ID
assignments <- cyclone(merge, hs.pairs)  # we want to use the whole set of genes for prediction (takes a while)
table(assignments$phases)
merge$cycle <- assignments$phases

merge <- merge[, paste0(merge$Barcode, merge$batch) %in% paste0(merge.hvg$Barcode, merge.hvg$batch)]
merge.hvg$cycle <- merge$cycle

cycle <- plotReducedDim(merge.hvg, "TSNE", colour_by = "cycle", text_by = "label")


cycle + scale_color_manual(values=c(G2M="red", S="green", G1 = "grey70")) +  
  labs(title = "Predicted cell cycle phase per cell", color = "Phase")

cycle_table <- table(merge.hvg$batch, merge.hvg$cycle)
for (x in 1:3){
  cycle_table[x,] <- cycle_table[x,] / sum(cycle_table[x,]) * 100
}
cycle_table


