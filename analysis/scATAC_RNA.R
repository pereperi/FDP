###########################################################################
###########################################################################
### MULTIOME SINGLE-CELL ATAC + GENE EXPRESSION ANALYSIS OF 3 PNF SAMPLES
###########################################################################
###########################################################################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
## SOME OF THE ANALYSES SHOWN MAY HAVE BEEN DISCONTINUED, 
## AND ARE THEREFORE NOT INCLUDED IN THE FDP MANUSCRIPT
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
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
library(stringr)
library(scDblFinder)
library(BiocSingular)
library(Seurat)
library(ggsignif)
library(celldex)
library(SingleR)
library(pheatmap)
library(stringr)
library(dplyr)
library(regioneR)

############################################################################################
############################################################################################
## LOADING THE DATA
############################################################################################
############################################################################################

set.seed(1234)

counts <- Read10X_h5("/imppc/labs/eslab/share/PerePericot/Projects/scRNAseq_scATAC_PNF/preprocessed_data/PNF_aggr/outs/filtered_feature_bc_matrix.h5")

fragpath <- "/imppc/labs/eslab/share/PerePericot/Projects/scRNAseq_scATAC_PNF/preprocessed_data/PNF_aggr/outs/atac_fragments.tsv.gz"

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

multi <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

rna <- as.SingleCellExperiment(multi) #we are going to process the rna assay separately with bioconductor, to be consistent with our previous analysis

############################################################################################
############################################################################################
############################################################################################
## GENE EXPRESSION DATA PROCESSING (BIOC)
############################################################################################
############################################################################################
############################################################################################



rna$batch <- colnames(rna)
rna$batch <- str_replace_all(rna$batch, ".*1", "PNF19")
rna$batch <- str_replace_all(rna$batch, ".*2", "PNF20")
rna$batch <- str_replace_all(rna$batch, ".*3", "PNF23")

rna$batch <- factor(rna$batch)
###########################################################################
## QUALITY CONTROL
###########################################################################

gns <- genes(EnsDb.Hsapiens.v75, filter = ~ seq_name == "MT")  #Identifying mitochondrial transcripts

is.mito <- rownames(rna) %in% gns$symbol
sum(is.mito)

df <- perCellQCMetrics(rna, subsets=list(Mito=is.mito)) #Generate the dataframe


reasons <- perCellQCFilters(df, sub.fields=c("subsets_Mito_percent")) #We make use of this function to identify outliers for our QC metrics

reasons$low_lib_size <- df$sum < 1000 
reasons$discard <- reasons$low_lib_size | reasons$low_n_features | reasons$high_subsets_Mito_percent
colSums(as.matrix(reasons)) #We obtain the number of cells to discard


rna <- rna[,!reasons$discard]  #We discard low-quality cells from further analysis


###########################################################################
## NORMALIZATION
###########################################################################

set.seed(1101)
lib.sf <- librarySizeFactors(rna)  #compute the library size factors
summary(lib.sf)
hist(log10(lib.sf), xlab="Log10[Size factor]", col='grey80')

rna <- logNormCounts(rna, size.factors = lib.sf)  #normalize the counts
assayNames(rna)


###########################################################################
## FEATURE SELECTION
###########################################################################

set.seed(11101)
dec.pois.filtered.norm <- modelGeneVarByPoisson(rna)
dec.pois.filtered.norm <- dec.pois.filtered.norm[order(dec.pois.filtered.norm$bio, decreasing=TRUE),] #order
dec.pois.filtered.norm #We will use this for downstream analysis

chosen <- getTopHVGs(dec.pois.filtered.norm, prop=0.2)  #get the 20% most important genes
str(chosen)

rna.hvg <- rna[chosen,]  #new sce object containing only HVGs
dim(rna.hvg)

###########################################################################
## BATCH CORRECTION
###########################################################################

# assay(rna.hvg, "corrected") <- removeBatchEffect(logcounts(rna.hvg), batch=rna.hvg$batch)
# assay(rna.hvg, "corrected") <- ComBat(logcounts(rna.hvg), batch = rna.hvg$batch)
# rna.hvg <- runPCA(rna.hvg, ncomponents=30, exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
# ?runPCA
# rna.hvg <- runPCA(rna.hvg, ncomponents = 30, exprs_values = "logcounts")
# set.seed(1010)
# rna.hvg <- runTSNE(rna.hvg, dimred="PCA", perplexity=80)
# plotReducedDim(rna.hvg, "TSNE", colour_by = "batch")  # need to correct


set.seed(1001100)
rna.mnn <- fastMNN(logcounts(rna.hvg), batch = rna.hvg$batch)
reducedDim(rna.hvg, "corrected") <- reducedDim(rna.mnn, "corrected")

###########################################################################
## DIMENSIONALITY REDUCTION
###########################################################################

set.seed(1010)
rna.hvg <- runTSNE(rna.hvg, dimred="corrected", perplexity=80)
plotReducedDim(rna.hvg, "TSNE", colour_by = "batch")


set.seed(1234)
rna.hvg <- runUMAP(rna.hvg, dimred="corrected")
plotReducedDim(rna.hvg, "UMAP", colour_by = "batch")

###########################################################################
## QUALITY CONTROL - 2nd STEP
###########################################################################


set.seed(100)
dbl.dens <- computeDoubletDensity(rna.hvg, d=ncol(reducedDim(rna.hvg, "corrected")))
summary(dbl.dens)

rna.hvg$DoubletScore <- dbl.dens
set.seed(1010)
plotTSNE(rna.hvg, colour_by="DoubletScore", text_colour = "red")

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call", p=0.01)

summary(dbl.calls)
plotColData(rna.hvg, x="batch", y="DoubletScore", colour_by=I(dbl.calls))

rna.hvg <- rna.hvg[,dbl.calls=="singlet"]



sce.seurat <- as.Seurat(rna.hvg, counts = "counts", data = "logcounts")

genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/ER_genes.txt", sep = "\t", header = FALSE)[,1]

genes <- genes[genes %in% rownames(rna.hvg)]
genes <- unique(genes)
set.seed(12345)
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")

rna.hvg$signature <- sce.seurat$signature1
signature_plot <- plotReducedDim(rna.hvg, "TSNE", colour_by = "signature")
p4 <- signature_plot + scale_colour_gradient2(low = "blue", mid = "#F7F7F7", high = "red",  midpoint = 0) 
p4 <- p4 + labs(color = "Scores")
p4

rna.hvg$stressed <- rna.hvg$signature > 0.2
rna.hvg_stressed <- rna.hvg
rna.hvg <- rna.hvg[,!rna.hvg$stressed]


violin <- data.frame(rna.hvg_stressed$batch, rna.hvg_stressed$signature)
ggplot(data = violin, mapping = aes(x = rna.hvg_stressed.batch, y = rna.hvg_stressed.signature)) + geom_violin(aes(fill = rna.hvg_stressed.batch)) + labs(x="Sample", y="Score of response to ER stress", fill = "Sample") +
  theme_linedraw() + geom_hline(aes(yintercept = 0.2), color = "red") + annotate(geom = "text", x=1.5, y=0.4, label ="Threshold = 0.2", color = "red")

plotReducedDim(rna.hvg_stressed, "TSNE", colour_by = "stressed") + scale_colour_manual(values = c("grey70", "red")) + labs(color = "Stressed (Score>0.35)")


cut_df <- data.frame(rna.hvg_stressed$batch, rna.hvg_stressed$stressed)
ggplot(cut_df, mapping = aes(rna.hvg_stressed.batch, fill = rna.hvg_stressed.stressed)) + geom_bar(position = position_stack(reverse = TRUE)) + 
  geom_text(aes(label = ..count..), stat = "count", position = position_stack(reverse=TRUE), ) + labs(x = "Sample", y= "Number of cells", fill = "Stressed (Score>0.2)") + scale_fill_manual(values = c("chartreuse3", "red3"))


###########################################################################
## BATCH CORRECTION
###########################################################################

# assay(rna.hvg, "corrected") <- removeBatchEffect(logcounts(rna.hvg), batch=rna.hvg$batch)
# assay(rna.hvg, "corrected") <- ComBat(logcounts(rna.hvg), batch = rna.hvg$batch)
# rna.hvg <- runPCA(rna.hvg, ncomponents=30, exprs_values="logcounts", BSPARAM=BiocSingular::ExactParam())
# set.seed(1010)
# rna.hvg <- runTSNE(rna.hvg, dimred="PCA", perplexity=80)
# plotReducedDim(rna.hvg, "TSNE", colour_by = "batch")  # need to correct


set.seed(1001100)
rna.mnn <- fastMNN(logcounts(rna.hvg), batch = rna.hvg$batch)
reducedDim(rna.hvg, "corrected") <- reducedDim(rna.mnn, "corrected")

###########################################################################
## DIMENSIONALITY REDUCTION
###########################################################################

set.seed(101010)
rna.hvg <- runTSNE(rna.hvg, dimred="corrected", perplexity=120)
plotReducedDim(rna.hvg, "TSNE", colour_by = "batch")


set.seed(1234)
rna.hvg <- runUMAP(rna.hvg, dimred="corrected")
plotReducedDim(rna.hvg, "UMAP", colour_by = "batch")


###########################################################################
## CLUSTERING
###########################################################################

set.seed(11111)
clust2s <- clusterCells(rna.hvg, use.dimred="corrected", BLUSPARAM=TwoStepParam(first=KmeansParam(centers=2500),
                                                                                second=NNGraphParam(k=35, type = "jaccard", cluster.fun = "infomap")))

table(clust2s)


colLabels(rna.hvg) <- clust2s

p1 <- plotReducedDim(rna.hvg, "TSNE", colour_by="label", text_by = "label")
p1

plotReducedDim(rna.hvg, "UMAP", colour_by="label", text_by = "label")


###########################################################################
## EXPRESSION PLOTS
###########################################################################
reducedDim(rna, "TSNE") <- reducedDim(rna.hvg, "TSNE")
rna$label <- rna.hvg$label


gene <- "CDH19"     #plug in any gene
p2 <- plotReducedDim(rna.hvg, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =6) +
  labs(title = paste0(gene, " expression per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()

p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)


genes <- c("S100B", "CDH19", "SOX10", "TWIST1", "PRRX1")

plots <- list()
i <- 1
for (gene in genes){
  p2 <- plotReducedDim(rna.hvg, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =6) +
    labs(title = paste0(gene, " expression per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()
  
  p2 <- p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)
  plots[[i]] <- p2
  i <- i+1
}
grid.arrange(grobs = plots, ncol = 3)


###########################################################################
## CELL TYPE ASSIGNMENT
###########################################################################


ref <- HumanPrimaryCellAtlasData()


pred <- SingleR(test=rna.hvg, ref=ref, labels=ref$label.fine)
table(pred$labels)

tab <- table(Assigned=pred$pruned.labels, Cluster=colLabels(rna.hvg))
# 
# c <- 12
# sort((tab[,c] / sum(tab[,c]))*100)  # obtain the percentages for each type
# # sum(tab[,c])  # number of cells in the cluster


pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101)) # Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.

rna.hvg$celltype <- pred$labels



rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Fibroblasts.*", "Fibroblasts")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Endothelial_cells.*", "Endothelial cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "B_cell.*", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "T_cell.*", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "NK.*", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Tissue_stem_cells.*", "Tissue stem cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "iPS_cells.*", "iPS cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Macrophage.*", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Monocyte.*", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "DC.*", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Smooth_muscle.*", "Smooth muscle cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "MEP", "Immune cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Chondrocytes.*", "Chondrocytes")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Neuron.*", "Schwann cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Neuro.*", "Neuroepithelial cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Epi.*", "Epithelial cells")
rna.hvg$celltype <- str_replace_all(rna.hvg$celltype, "Keratinocytes.*", "Epithelial cells")
types <- rownames(table(rna.hvg$celltype)[table(rna.hvg$celltype)>50])

rna.hvg$celltype[!rna.hvg$celltype %in% types] <- "others"

rna.hvg$celltype.main <- rna.hvg$celltype
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Schwann cells", "Schwann cell")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Neuroepithelial cells", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Immune cells", "Immune cells")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Fibroblasts", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Tissue stem cells", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Smooth muscle cells", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Chondrocytes", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Epithelial cells", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Astrocyte.*", "Schwann cell")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Fibroblasts", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "Osteoblasts", "Mesenchymal")
rna.hvg$celltype.main <- str_replace_all(rna.hvg$celltype.main, "MSC", "Mesenchymal")


plotReducedDim(rna.hvg, "TSNE", colour_by="celltype", text_by = "label") + scale_color_manual(name = "Cell Type", values = c("Immune cells" = "red", "Schwann cells" = "forestgreen", "Endothelial cells" = "orange",
                                                                                                                             "Fibroblasts" = "firebrick4", "Smooth muscle cells" = "peru", "Tissue stem cells" = "lightpink", 
                                                                                                                             "MSC" = "purple", "iPS cells" = "turquoise", "Astrocyte:Embryonic_stem_cell-derived" = "blue4", 
                                                                                                                             "others" = "azure4", "Epithelial cells" = "yellow", "Neuroepithelial cells" = "green", "Chondrocytes" = "magenta"))


plotReducedDim(rna.hvg, "TSNE", colour_by="celltype.main", text_by = "label") + scale_color_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey"))


celltype_prop <- prop.table(as.array(table(rna.hvg$celltype, rna.hvg$label)), margin = 2)*100

df <- data.frame(rna.hvg$label, rna.hvg$celltype)

ggplot(df, mapping = aes(x=rna.hvg.label, fill=rna.hvg.celltype)) + geom_bar(position = "fill") + labs(x = "Cluster", y= "Proportion") +
  theme_minimal() + theme(legend.position = "bottom")  + scale_fill_manual(name = "Cell Type", values = c("Immune cells" = "red", "Schwann cells" = "forestgreen", "Endothelial cells" = "orange", "Fibroblasts" = "firebrick4", "Smooth muscle cells" = "peru", "Tissue stem cells" = "lightpink", "MSC" = "purple", "iPS cells" = "turquoise", "Astrocyte:Embryonic_stem_cell-derived" = "blue4","others" = "azure4", "Epithelial cells" = "yellow", "Neuroepithelial cells" = "green", "Chondrocytes" = "magenta")) 
                                                                                                          

############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

# GENERATION OF SOME PLOTS

genes <- c("SOX10", "CDH19", "S100B", "TWIST1","PRRX1")

plots <- list()
i <- 1
for (gene in genes){
  p2 <- plotReducedDim(rna.hvg, dimred = "TSNE", colour_by =gene, by_exprs_values = "logcounts", text_by = "label", text_colour = "black", text_size =10) +
    theme_classic() + theme(axis.title  = element_text(size=20), legend.text = element_text(size=30), legend.key.size = unit(1.5, "cm"), legend.title = element_text(size=30))
  
  p2 <- p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)
  plots[[i]] <- p2
  i <- i+1
}
p2 <- grid.arrange(grobs = plots, ncol = 2)

grid.arrange(p1,p2, nrow = 2)




plotReducedDim(rna.hvg, "TSNE", colour_by="celltype.main", text_by = "label", text_size = 10) + scale_color_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey")) +
  theme(text = element_text(size = 30)) + guides(colour = guide_legend(override.aes = list(shape = 16, size = 10))) + theme(legend.position = "blank")

df <- data.frame(batch = sce_unstressed$batch, component = sce_unstressed$celltype.main, exp = "Conventional")
df2 <- data.frame(batch = rna.hvg$batch,component = rna.hvg$celltype.main, exp = "Multiome")
df <- rbind(df, df2)
ggplot(df, mapping = aes(x=exp, fill=component)) + geom_bar(position = "fill") + labs(x = "Sample of origin", y= "Proportion") + facet_grid(cols = vars(batch)) +
  theme_minimal() + theme(legend.position = "bottom", text = element_text(size = 30))  + scale_fill_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal" = "orange", "Schwann cell" = "#6fc364", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey")) + theme(legend.position = "blank", axis.title.x = element_blank())

############################################################################################
############################################################################################
## ATAC DATA PROCESSING AND TRANSFER OF METADATA
############################################################################################
############################################################################################

# create ATAC assay and add it to the object
multi[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

multi <- multi[,colnames(multi) %in% colnames(rna.hvg)]

multi@meta.data$batch <- rna.hvg$batch
multi@meta.data$celltype <- rna.hvg$celltype
multi@meta.data$celltype.main <- rna.hvg$celltype.main
multi@meta.data$label_rna <- rna.hvg$label


multi[["corrected"]] <- CreateDimReducObject(embeddings = reducedDim(rna.hvg, "corrected"), key = "corrected_", assay = DefaultAssay(multi))

DefaultAssay(multi) <- "ATAC"

multi <- NucleosomeSignal(multi)
multi <- TSSEnrichment(multi)

VlnPlot(
  object = multi,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

# filter out low quality cells (ATAC QC)
multi <- subset(
  x = multi,
  subset = nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
multi

############################################################################################
############################################################################################
## PEAK CALLING WITH MACS2
############################################################################################
############################################################################################


# call peaks using MACS2
set.seed(1010)
peaks <- CallPeaks(multi, macs2.path = "/software/debian-10/bin/macs2")

  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
set.seed(111)
macs2_counts <- FeatureMatrix(
  fragments = Fragments(multi),
  features = peaks,
  cells = colnames(multi)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
multi[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

############################################################################################
############################################################################################
## DIMENSIONALITY REDUCTION (NOT USED FOR ATAC)
############################################################################################
############################################################################################


DefaultAssay(multi) <- "peaks"
multi <- FindTopFeatures(multi, min.cutoff = 5)
set.seed(1101)
multi <- RunTFIDF(multi)
set.seed(11020)
multi <- RunSVD(multi)

DepthCor(multi) # first component captures sequencing depth rather than biological variation, so we remove it from downstream analysis

# build a joint neighbor graph using both assays

set.seed(1101)
multi <- FindMultiModalNeighbors(
  object = multi,
  reduction.list = list("corrected", "lsi"), 
  dims.list = list(1:50, 2:50),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
set.seed(1010)
multi <- RunUMAP(
  object = multi,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(multi, label = TRUE, repel = TRUE, reduction = "umap") # joint dimensionality recudtion visualisation

set.seed(1234)
multi <- RunUMAP(object = multi, reduction = 'lsi', dims = 2:50, reduction.name = "umap_lsi") #using ATAC


############################################################################################
############################################################################################
## CONVERT TO SingleCellExperiment FOR BETTER MANIPULATION
############################################################################################
############################################################################################

multi_sce <- as.SingleCellExperiment(multi)

multi_sce <- swapAltExp(multi_sce, "RNA")

set.seed(1234)
multi_sce <- runUMAP(multi_sce, dimred="CORRECTED", name = "umap_rna")
plotReducedDim(multi_sce, "umap_rna", colour_by = "celltype.main")

set.seed(11111)
clust2s <- clusterCells(multi_sce, use.dimred="CORRECTED", BLUSPARAM=TwoStepParam(first=KmeansParam(centers=1000),
                                                                                second=NNGraphParam(k=35, type = "jaccard", cluster.fun = "infomap")))

table(clust2s)


multi_sce$label_rna <- clust2s
plotReducedDim(multi_sce, "umap_rna", colour_by = "label_rna")

multi_sce <- swapAltExp(multi_sce, "RNA")
plotReducedDim(multi_sce, "UMAP", colour_by="batch") + labs(x = "UMAP 1", y = "UMAP 2") 
multi_sce <- swapAltExp(multi_sce, "peaks")
plotReducedDim(multi_sce, "UMAP_LSI", colour_by="celltype.main") + labs(x = "UMAP 1", y = "UMAP 2") + scale_color_manual(name = "Cellular components", values = c("Immune cells" = "red", "Mesenchymal component" = "orange", "Schwann cell component" = "green", "Endothelial cells" = "pink", "iPS cells" = "turquoise", "others" = "grey"))

############################################################################################
############################################################################################
## EXPRESSION AND ACCESSIBILITY PLOTS
############################################################################################
############################################################################################

multi_sce <- swapAltExp(multi_sce, "RNA")
assay(multi_sce, "logcounts") <- assay(rna, "logcounts")

gene <- "SOX10"     #plug in any gene
p2 <- plotReducedDim(multi_sce, dimred = "umap_rna", colour_by =gene, by_exprs_values = "logcounts", text_by = "label_rna", text_colour = "black", text_size =6) +
  labs(title = paste0(gene, " expression per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()

p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)

multi_sce <- swapAltExp(multi_sce, "peaks")
reducedDim(multi_sce, "umap_rna") <- reducedDim(altExp(multi_sce,"RNA"), "umap_rna")
gene <- "chr9-129760957-129763749"     #plug in any region
p2 <- plotReducedDim(multi_sce, dimred = "umap_rna", colour_by =gene, by_exprs_values = "logcounts", text_by = "label_rna", text_colour = "black", text_size =6) +
  labs(title = paste0(gene, " accessibility per cell")) + theme(plot.title = element_text(size=15)) + theme_classic()

p2 + scale_color_gradientn(colours = colorRampPalette(c("grey70", "orange3", "firebrick", "firebrick", "red", "red"))(10)) +  labs(color = gene)

# SOME RELEVANT PEAKS:
# S100B: chr21-46604748-46605651
# CDH19: chr18-66603601-66604573
# SOX10. chr22-37985406-37986360
# chr13-113410049-113411935
# PRRX1: chr1-170729009-170730027 
# EGFR: chr7-55132597-55134204
############################################################################################
############################################################################################
## PREDICTING LINKS  (TAKES ~12h FOR WHOLE GENOME)
############################################################################################
############################################################################################

DefaultAssay(multi) <- "peaks"

# first compute the GC content for each peak
multi <- RegionStats(multi, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
multi <- LinkPeaks(
  object = multi,
  peak.assay = "peaks",
  expression.assay = "RNA"
)


############################################################################################
############################################################################################
## ASSIGNING CELLS TO ROADMAP STAGE (NC, Early, Mature, or Fb)
############################################################################################
############################################################################################

sce.seurat <- as.Seurat(rna.hvg, counts = "counts", data = "logcounts")
sce.seurat <- sce.seurat[, colnames(sce.seurat) %in% colnames(multi)]
genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/up_fb.txt", sep = "\t", header = FALSE)[,1]


genes <- genes[genes %in% rownames(rna.hvg)]
genes <- unique(genes)
set.seed(12345)
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")

multi_sce$signature <- sce.seurat$signature1
signature_plot <- plotReducedDim(multi_sce, "umap_rna", colour_by = "signature")
p4 <- signature_plot + scale_colour_gradient2(low = "blue", mid = "#F7F7F7", high = "red",  midpoint = 0) 
p4 <- p4 + labs(color = "Scores")
p4


multi_sce$signature_fb <- multi_sce$signature / norm(as.matrix(multi_sce$signature), type = "2") #do the same for the other signatures

df <- data.frame("NC" = multi_sce$signature_nc, "Early" = multi_sce$signature_7, "Mature" = multi_sce$signature_30, "Fb" = multi_sce$signature_fb,"Unassigned" = 0.015)

boxplot(df)

multi@meta.data$stage <- colnames(df)[max.col(df,ties.method="first")]
multi_sce$stage <- colnames(df)[max.col(df,ties.method="first")]
summary(as.factor(multi_sce$stage))

multi$stage[multi_sce$label_rna != 3 & multi_sce$stage == "NC"] <- "Unassigned"
multi_sce$stage[multi_sce$label_rna != 3 & multi_sce$stage == "NC"] <- "Unassigned"
summary(as.factor(multi_sce$stage))
plotReducedDim(multi_sce, "umap_rna", colour_by = "stage")

############################################################################################
############################################################################################
## COVERAGE PLOTS (INTERACTIVE)
############################################################################################
############################################################################################


gene <- "SOX10" #plug in any gene

multi_sce <- swapAltExp(multi_sce, "RNA")
multi_sce$expression <- assay(multi_sce, "logcounts")[gene,]

quartile <- ntile(multi_sce$expression[multi_sce$expression != 0], 5)
multi_sce$expression[multi_sce$expression != 0] <- quartile
multi_sce$expression[multi_sce$expression == 0] <- "No expression"
multi_sce$expression[multi_sce$expression == 1] <- "Low expression"
multi_sce$expression[multi_sce$expression == 2] <- "Low expression"
multi_sce$expression[multi_sce$expression == 3] <- "Low expression"
multi_sce$expression[multi_sce$expression == 4] <- "Low expression"
multi_sce$expression[multi_sce$expression == 5] <- "High expression"

multi@meta.data$expression <- multi_sce$expression
multi@meta.data$label_rna <- multi_sce$label_rna



Idents(multi) <- multi$stage
multi@assays$RNA@counts <- assay(rna, "logcounts")
CoverageBrowser(multi, region = "STMN1") #interactive plots 


############################################################################################
############################################################################################
## DIFFERENTIAL ACCESSIBILITY ANALYSIS
############################################################################################
############################################################################################

DefaultAssay(multi) <- 'peaks'

da_peaks <- FindMarkers(
  object = multi,
  ident.1 = "Schwann cell component",
  ident.2 = "Mesenchymal component",
  group.by = "celltype.main",
  min.pct = 0.05,
  test.use = 'roc'
)
# chr6-101391291-101394251  for SCs (cluster 1)
# chr13-113410049-113411935, chr9-129760957-129763749 for MC (cluster 7)

ClosestFeature(multi, regions = rownames(head(da_peaks, 200)))
da_peaks[da_peaks$avg_diff<0,]
da_peaks[da_peaks$myAUC>0.4,]


CoverageBrowser(multi, region = "chr20-58016858-58018164") #interactive plots

multi_sce$double <- multi@assays$peaks@counts["chr13-113410049-113411935",] > 0 & multi@assays$peaks@counts["chr6-101391291-101394251",] > 0
plotReducedDim(multi_sce, "umap_rna", colour_by = "double") + scale_color_manual(values = c("grey", "red"))

pp <- toGRanges(rownames(multi@assays$peaks@counts))
subsetByOverlaps(pp, toGRanges("chr1-170660000-170662000"))



##########################################################################################
##########################################################################################

sce.seurat <- as.Seurat(rna.hvg, counts = "counts", data = "logcounts")
sce.seurat <- sce.seurat[, colnames(sce.seurat) %in% colnames(multi)]
genes <- read.delim("/imppc/labs/eslab/share/PerePericot/Projects/scRNA_PNF/results/up_fb.txt", sep = "\t", header = FALSE)[,1]

genes <- toupper(x = genes)
genes <- genes[genes %in% rownames(rna.hvg)]
genes <- unique(genes)
set.seed(12345)
sce.seurat <- AddModuleScore(sce.seurat,
                             features = list(genes),
                             name="signature")

multi_sce$signature <- sce.seurat$signature1
signature_plot <- plotReducedDim(multi_sce, "umap_rna", colour_by = "signature")
p4 <- signature_plot + scale_colour_gradient2(low = "blue", mid = "#F7F7F7", high = "red",  midpoint = 0) 
p4 <- p4 + labs(color = "Scores")
p4

sum(multi$stage == "NC" & multi$celltype.main == "Mesenchymal component")
