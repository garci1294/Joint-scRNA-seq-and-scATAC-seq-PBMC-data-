library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(ggplot2)
library(SeuratDisk)
library(patchwork)

set.seed(1234)


#==============================================================================================================
#reference

@UNPUBLISHED{signac,
  title    = "Multimodal single-cell chromatin analysis with Signac",
  author   = "Stuart, Tim and Srivastava, Avi and Lareau, Caleb and Satija,
              Rahul",
  journal  = "bioRxiv",
  pages    = "2020.11.09.373613",
  month    =  nov,
  year     =  2020,
  url      = "https://www.biorxiv.org/content/10.1101/2020.11.09.373613v1",
  language = "en"
}

##   Cole Trapnell and Davide Cacchiarelli et al (2014): The dynamics
##   and regulators of cell fate decisions are revealed by
##   pseudo-temporal ordering of single cells. Nature Biotechnology
##
##   Xiaojie Qiu, Andrew Hill, Cole Trapnell et al (2017):
##   Single-cell mRNA quantification and differential analysis with
##   Census. Nature Methods
##
##   Xiaojie Qiu, Cole Trapnell et al (2017): Reverse graph embedding
##   resolves complex single-cell developmental trajectories. BioRxiv

#==============================================================================================================
# load data and generate objects

# load the RNA and ATAC data
counts <- Read10X_h5("/Users/garci624/Desktop/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
fragpath <- "/Users/garci624/Desktop/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"

# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
)

Annotation(pbmc[["ATAC"]]) <- annotation

# create fragment object
fragments <- CreateFragmentObject(
  path = fragpath,
  cells = colnames(pbmc),
  validate.fragments = FALSE
)

# fraction counts in regions from ATAC data
pbmc$blacklist_fraction <- FractionCountsInRegion(
  object = pbmc,
  assay = 'ATAC',
  regions = blacklist_hg38
)

Fragments(pbmc[["ATAC"]]) <- fragments

DefaultAssay(pbmc) <- "ATAC"

# obtain nucleosome signal and TSS enrichements
pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc)

#==============================================================================================================
# generate violin plots

VlnPlot(
  object = pbmc,
  features = "nCount_RNA",
  cols = "turquoise",
  ncol = 1,
  pt.size = 0.1
)

VlnPlot(
  object = pbmc,
  features = "nCount_ATAC",
  cols = "turquoise",
  ncol = 1,
  pt.size = 0.1
)

VlnPlot(
  object = pbmc,
  features = "TSS.enrichment",
  cols = "turquoise",
  ncol = 1,
  pt.size = 0.1
)

VlnPlot(
  object = pbmc,
  features = "nucleosome_signal",
  cols = "turquoise",
  ncol = 1,
  pt.size = 0.1
)

VlnPlot(
  object = pbmc,
  features = "blacklist_fraction",
  cols = "turquoise",
  ncol = 1,
  pt.size = 0.1
)

# filter out low quality cells
pbmc <- subset(
  x = pbmc,
  subset = nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    nCount_ATAC < 100000 &
    nCount_ATAC > 1000 &
    TSS.enrichment > 1 &
    nucleosome_signal < 2 &
    blacklist_fraction < 0.015
)

#==============================================================================================================
# calling peaks

# call peaks using MACS2
peaks <- CallPeaks(pbmc, macs2.path = "/Users/garci624/opt/miniconda3/envs/py2/bin/macs2")

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pbmc),
  features = peaks,
  cells = colnames(pbmc)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pbmc[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
)

Annotation(pbmc[["ATAC"]]) <- annotation
#==============================================================================================================
# RNA and peaks assay analysis

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc)
pbmc <- RunPCA(pbmc)

DefaultAssay(pbmc) <- "peaks"
pbmc <- FindTopFeatures(pbmc, min.cutoff = 5)
pbmc <- RunTFIDF(pbmc)
pbmc <- RunSVD(pbmc)
#==============================================================================================================
# cell annotation

# load PBMC reference
reference <- LoadH5Seurat("/Users/garci624/Desktop/pbmc_multimodal.h5seurat")

DefaultAssay(pbmc) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype.l2,
  weight.reduction = pbmc[['pca']],
  dims = 1:50
)

pbmc <- AddMetaData(
  object = pbmc,
  metadata = predictions
)

# set the cell identities to the cell type predictions
Idents(pbmc) <- "predicted.id"

# set a reasonable order for cell types to be displayed when plotting
levels(pbmc) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT","CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", 
                  "NK", "NK_CD56bright","NK Proliferating", "gdT","Treg", "B naive", 
                  "B intermediate", "B memory", "Plasmablast","CD14 Mono", "CD16 Mono", 
                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")
#==============================================================================================================
# single cell RNA-ATAC multiomics analysis

# build a joint neighbor graph using both assays
pbmc <- FindMultiModalNeighbors(
  object = pbmc,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
pbmc <- RunUMAP(
  object = pbmc,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(pbmc, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

#==============================================================================================================
# Signal analysis

DefaultAssay(pbmc) <- "peaks"

cov_plot <- CoveragePlot(
  object = pbmc,
  region = "chr1-40189344-40252549",
  annotation = FALSE,
  peaks = FALSE
)
cov_plot

#==============================================================================================================
# MOTIF analysis

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
pbmc <- AddMotifs(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "CD4 Naive",
  ident.2 = "B memory",
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

# test enrichment
enriched.motifs <- FindMotifs(
  object = pbmc,
  features = top.da.peak
)

MotifPlot(
  object = pbmc,
  motifs = head(rownames(enriched.motifs))
)

#==============================================================================================================
# FOOTPRINT analysis

# gather the footprinting information for sets of motifs
pbmc <- Footprint(
  object = pbmc,
  motif.name = c("NRF1", "ZBTB14", "EGR3", "KLF15", "EGR1", "KLF14"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(pbmc, features = "NRF1")
p2 + patchwork::plot_layout(ncol = 1)

p3 <- PlotFootprint(pbmc, features = "ZBTB14")
p3 + patchwork::plot_layout(ncol = 1)

p4 <- PlotFootprint(pbmc, features = "EGR3")
p4 + patchwork::plot_layout(ncol =1)

p5 <- PlotFootprint(pbmc, features = "KLF15")
p5 + patchwork::plot_layout(ncol = 1)

p6 <- PlotFootprint(pbmc, features = "EGR1")
p6 + patchwork::plot_layout(ncol = 1)

p7 <- PlotFootprint(pbmc, features = "KLF14")
p7 + patchwork::plot_layout(ncol =1)
#==============================================================================================================
# PSEUDOTIME analysis

DefaultAssay(pbmc) <- "RNA"

pbmc.cds <- as.cell_data_set(pbmc)
pbmc.cds <- cluster_cells(cds = pbmc.cds, reduction_method = "UMAP")
pbmc.cds <- learn_graph(pbmc.cds, use_partition = TRUE)

plot_cells(pbmc.cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)

pbmc.cds <- order_cells(pbmc.cds)
plot_cells(pbmc.cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)

#==============================================================================================================


