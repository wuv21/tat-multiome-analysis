library(ArchR)
suppressMessages(library(Seurat))
suppressMessages(library(glue))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ggrastr))

dyn.load("/usr/local/hdf5/lib/libhdf5_hl.so.310")
library(hdf5r)

set.seed(42)

addArchRThreads(threads = 12)

source("scripts/helperUtils.R")
source("scripts/graphs.R")

source("scripts/plotBrowserTrackCustom.R")
tmpFun <- get("plotBrowserTrack", envir = asNamespace("ArchR"))
environment(plotBrowserTrack2) <- environment(tmpFun)
attributes(plotBrowserTrack2) <- attributes(tmpFun)


COLORS <- list(
  treatment = c(
    "18hr Tat" = "#1b9e77",
    "18hr Luc" = "#d95f02",
    "24hr Tat" = "#7570b3",
    "24hr Luc" = "#e7298a",
    "72hr Tat" = "#66a61e",
    "72hr Luc" = "#666666",
    "16hr P/I" = "#e6ab02",
    "24hr P/I" = "#a6761d"),
  deg = c("Up in Tat tx" = "#d95f02", "Up in Luc tx" = "#1b9e77", "n/a" = "#cccccc")
)

################################################################################
# Downstream preprocessing
################################################################################

samples_late <- list.files("data/cellranger_late_out/")
samples_early <- list.files("data/cellranger_early_out/")
samples <- c(samples_late, samples_early)

proj <- loadArchRProject(
  path = "dogma_qcFilt_all_v1",
  showLogo = FALSE)


spMatRna <- lapply(samples, function(x) {
  message(x)
  
  if (x %in% c("D2_10-6_Tat_2", "D3_10-6_Luc_2", "D3_10-6_Tat_2", "NA_10-6_PmaIono")) {
    crOutDir <- "cellranger_late_out"
  } else {
    crOutDir <- "cellranger_early_out"
  }
  
  mat <- Read10X_h5(glue("data/{crOutDir}/{x}/outs/filtered_feature_bc_matrix.h5"))
  mat2 <- mat$`Gene Expression`
  
  colnames(mat2) <- paste0(x, "#", colnames(mat2))
  return(mat2)
})

spMatRnaFinal <- BiocGenerics::Reduce(cbind, spMatRna)
spMatRnaFinal <- spMatRnaFinal[, proj$cellNames]

rnaAssay <- CreateAssayObject(counts = spMatRnaFinal)
seu <- CreateSeuratObject(rnaAssay)
detach("package:hdf5r")

seu <- NormalizeData(seu, assay = "RNA")

transferToSeu <- data.frame(
  sample = proj$Sample,
  gex_UMI = proj$Gex_nUMI
)
rownames(transferToSeu) <- proj$cellNames
seu <- AddMetaData(seu, transferToSeu)

proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  saveIterations = FALSE,
  name = "LSI_ATAC",
  iterations = 3,
  firstSelection = "Top",
  varFeatures = 20000,
  dimsToUse = 1:30,
  sampleCellsPre = 15000,
  threads = 12,
  clusterParams = list(
    resolution = 2, 
    sampleCells = 15000,
    maxClusters = 25
  )
)

proj <- addIterativeLSI(
  ArchRProj = proj, 
  clusterParams = list(
    resolution = 2, 
    sampleCells = 15000,
    maxClusters = 25
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  sampleCellsPre = 15000,
  firstSelection = "Var",
  binarize = FALSE,
  iterations = 3,
  name = "LSI_RNA"
)

# add harmony based on experiment date
proj$experimentDate <- ifelse(proj$Sample %in% c("D2_10-6_Tat_2", "D3_10-6_Luc_2", "D3_10-6_Tat_2", "NA_10-6_PmaIono"), "late", "early")


proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  name = "LSI_ATAC_harmonized",
  groupBy = "experimentDate"
)

proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_RNA",
  name = "LSI_RNA_harmonized",
  groupBy = "experimentDate"
)

proj <- addCombinedDims(
  proj,
  reducedDims = c("LSI_ATAC_harmonized", "LSI_RNA_harmonized"),
  name = "LSI_Combined_harmonized")

proj@reducedDims$LsiCleanForSeurat <- proj@reducedDims$LSI_Combined_harmonized

# taken from ArchR MultiModal.R to calculate variance by column
# will sort these by variance explained in case FindNeighbors or FindClusters 
# are influenced by order (since FindClusters will throw an error if 
# the colnames are not labeled in sequential order)
cV <- apply(proj@reducedDims$LsiCleanForSeurat$matRD, 2, var) 

proj@reducedDims$LsiCleanForSeurat$matRD <- 
  proj@reducedDims$LsiCleanForSeurat$matRD[, order(cV, decreasing=TRUE)]

colnames(proj@reducedDims$LsiCleanForSeurat$matRD) <-
  paste0("tmp", c(1: ncol(proj@reducedDims$LsiCleanForSeurat$matRD)))

proj <- addUMAP(
  proj,
  reducedDims = "LsiCleanForSeurat",
  name = "UMAP_Combined",
  minDist = 0.1,
  nNeighbors = 20,
  force = TRUE)

projFinal <- saveArchRProject(proj, outputDirectory = "dogma_umap_all_v1", load = TRUE)

save.image("rdata/all_postUmap.rdata")


################################################################################
# Start here 
################################################################################
load("rdata/all_postUmap.rdata")

projFinal$experimentDate <- ifelse(projFinal$Sample %in% c("D2_10-6_Tat_2", "D3_10-6_Luc_2", "D3_10-6_Tat_2", "NA_10-6_PmaIono"), "late", "early")

treatmentData <- data.frame(sample = projFinal$Sample, expDate = projFinal$experimentDate) %>%
  dplyr::mutate(sample = case_when(
    sample == "18_luc" ~ "18_10-6_Luc",
    sample == "18_tat" ~ "18_10-6_Tat",
    sample == "24_luc" ~ "24_10-6_Luc",
    sample == "24_tat" ~ "24_10-6_Tat",
    sample == "24_pmaI" ~ "24_10-6_PmaIono",
    TRUE ~ sample)) %>%
  tidyr::separate(sample, sep = "_", into = c("time", "jlat", "treatment", "extra"), remove = FALSE) %>%
  dplyr::mutate(treatment = gsub("PmaIono", "P/I", treatment)) %>%
  dplyr::mutate(newTime = case_when(
    expDate == "late" & treatment == "P/I" ~ 16,
    time == "D2" ~ 48,
    time == "D3" ~ 72,
    TRUE ~ as.numeric(time)))

projFinal$treatment <- treatmentData$treatment
projFinal$time <- treatmentData$newTime
projFinal$combinedTreatment <- paste0(treatmentData$newTime, "hr ", treatmentData$treatment)

################################################################################
# Fig2B - umap
################################################################################
umapDf <- getEmbedding(proj, embedding = "UMAP_Combined")
stopifnot(sum(rownames(umapDf) == proj$cellNames) == length(proj$cellNames))

umapDf$time <- proj$time
umapDf$treatment <- proj$treatmnet 
umapDf$combinedTreatment <- proj$combinedTreatment

umapDf <- umapDf %>%
  dplyr::filter(time != 48) %>%
  dplyr::mutate(combinedTreatment = factor(combinedTreatment, levels = names(COLORS$treatment)))

# randomize
umapDf <- dplyr::slice(umapDf, sample(1:n()))
umapTheme <- theme(
  legend.position = "right",
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.title = element_blank(),
  axis.title = element_text(size = BASEPTFONTSIZE),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line = element_blank(),
  legend.box.margin = margin(t = 0, b = -10),
  plot.margin = margin(l = -50),
  plot.background = element_blank(),
  legend.key.size = unit(BASEPTFONTSIZE, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA))

figA <- umapDf %>% 
  ggplot(aes(x = `LsiCleanForSeurat#UMAP_Dimension_1`, y = `LsiCleanForSeurat#UMAP_Dimension_2`, color = combinedTreatment)) +
  geom_hline(yintercept = -Inf, linewidth = 12/14) +
  geom_vline(xintercept = -Inf, linewidth = 12/14) +
  rasterize(geom_point(alpha = 0.3, size = 0.05), dpi = 300) +
  guides(colour = guide_legend(override.aes = list(size = 1), ncol = 1)) +
  scale_y_continuous(limits = c(-9, 15)) +
  scale_x_continuous(limits = c(-9, 15)) +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = COLORS$treatment) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic() +
  umapTheme


################################################################################
# Fig2b - atac track
################################################################################
fig2Track <- plotBrowserTrack2(
  proj,
  geneSymbol = "HIV",
  groupBy = "combinedTreatment",
  useGroups = names(COLORS$treatment),
  plotSummary = c("bulkTrack", "geneTrack"),
  highlight = makeGRangesFromDataFrame(data.frame(seqname = "chr9", start = 136468584, end = 136478783)),
  upstream = 2500,
  downstream = 12500,
  minCells = 1000,
  sizes = c(3, 0))

fig2TrackFinal <- fig2Track$HIV$bulktrack +
  labs(y = "Bulk ATAC Signal\n(Normalized)") +
  scale_fill_manual(values = COLORS$treatment) +
  scale_color_manual(values = COLORS$treatment) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = BASEPTFONTSIZE),
    axis.text = element_text(size = BASEPTAXISFONTSIZE),
    strip.text = element_text(size = BASEPTFONTSIZE, hjust = 0),
    plot.margin = margin(b = 0, l = 10),
    strip.background = element_rect(color = "#00000000", fill = "#00000000"),
  )


################################################################################
# Fig2C - vln plot of rna
################################################################################
sum(proj$cellNames == Cells(seu))
proj

seu$combinedTreatment <- proj$combinedTreatment

fig2CDf <- data.frame(HIV = seu[["RNA"]]@data["HIV",], combinedTreatment = seu$combinedTreatment)

fig2C <- fig2CDf %>%
  dplyr::filter(combinedTreatment != "48hr Tat") %>%
  dplyr::mutate(combinedTreatment = factor(combinedTreatment, levels = names(COLORS$treatment))) %>%
  ggplot(aes(x = combinedTreatment, y = HIV)) + 
  geom_hline(yintercept = 0, linewidth = 12/14) +
  geom_vline(xintercept = -Inf, linewidth = 12/14) +
  geom_jitter(height = 0, size = 0.05, alpha = 0.1) +
  geom_violin(aes(fill = combinedTreatment), trim = TRUE, scale = "width", alpha = 0.8) + 
  scale_fill_manual(values = COLORS$treatment) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(y = "HIV RNA expression") +
  theme_classic() +
  theme(legend.position = "blank",
    axis.text = element_text(size = BASEPTAXISFONTSIZE),
    axis.text.x = element_text(size = BASEPTAXISFONTSIZE, angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_text(size = BASEPTFONTSIZE),
    plot.title = element_blank())

################################################################################
# Fig2d - deg analysis/volcano
################################################################################
seu$time <- proj$time
seu$treatment <- proj$treatment 
seu$combinedTreatment <- proj$combinedTreatment

seu <- ScaleData(seu, assay = "RNA")
seu <- FindVariableFeatures(seu, assay = "RNA", nfeatures = 2000)
seu <- RunPCA(seu, assay = "RNA")

Idents(seu) <- "combinedTreatment"
deg_lnps <- FindMarkers(
  seu,
  ident.1 = "72hr Luc",
  ident.2 = "72hr Tat",
  logfc.threshold = 0,
  min.pct = 0.2,
  slot = "data",
  test.use = "wilcox")

fig2d_df <- deg_lnps %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::mutate(pValAdj2 = p.adjust(p_val, method = "holm")) %>%
  dplyr::mutate(pValAdj2 = ifelse(pValAdj2 == 0, .Machine$double.xmin, pValAdj2)) %>%
  dplyr::mutate(importance = case_when(
    pValAdj2 < 0.01 & avg_log2FC >= 0.5 ~ "Up in Luc tx",
    pValAdj2 < 0.01 & avg_log2FC <= -0.5 ~ "Up in Tat tx",
    TRUE ~ "n/a"
  )) %>%
  dplyr::group_by(importance) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::mutate(group_rank = row_number()) %>%
  dplyr::mutate(label = ifelse(importance != "n/a" & group_rank < 10, gene, "")) 

fig2d_label <- fig2d_df %>%
  dplyr::group_by(importance) %>%
  dplyr::slice_sample(n = 150)


fig2d <- ggplot(fig2d_df, aes(x = avg_log2FC, y = -log10(p_val_adj), fill = importance)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#00000030") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#00000030") +
  geom_hline(yintercept = 0, linewidth = 12/14) +
  geom_vline(xintercept = -Inf, linewidth = 12/14) +
  geom_point(pch = 21, color = "#000000", alpha = 0.3, size = 0.5) +
  guides(fill = guide_legend(override.aes = list(alpha = 1, size = 1))) +
  ggrepel::geom_text_repel(data = fig2d_label,
    aes(label = label),
    max.overlaps = 300,
    segment.alpha = 0.3,
    force = 5,
    min.segment.length = 0,
    segment.size = 0.5,
    size = 1.5) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(-4, 2), expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = COLORS$deg) +
  theme(
    legend.text = element_text(size = BASEPTFONTSIZE),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(b = 0, t = 5, r = 0, l = 0),
    panel.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
    legend.box.margin = margin(t = -3, b = 0, r = 0, l = 0),
    legend.margin = margin(t = -3, b = 0, r = 0, l = 0),
    legend.background = element_rect(color = "#00000000", fill = "#00000000"),
    legend.box.background = element_rect(color = "#00000000", fill = "#00000000"),
    legend.key.size = unit(BASEPTFONTSIZE, 'points'),
    axis.title = element_text(size = BASEPTFONTSIZE),
    axis.text = element_text(color = "#000000", size = BASEPTAXISFONTSIZE - 2),
    axis.line = element_blank()) +
  labs(x = "Average log2 Fold Change",
    y = "-log10(adjusted p)")


################################################################################
# fig2e 
################################################################################
proApotosisGenes <- c("DIABLO", "ENDOG", "AIFM1", "BAK1", "BAX", "PMAIP1", "BID", "BAD",
  "BIK", "BMF", "BCL2L11", "BBC3", "FADD", "TRADD", "DISC1", "FASLG", "TNFSF10")

figDotPlotApop <- DotPlot(seu, features = proApotosisGenes, group.by = "combinedTreatment") +
  coord_flip() +
  scale_y_discrete(limits = names(COLORS$treatment)) +
  labs(x = "Pro-apoptosis genes") +
  theme(
    legend.key.size = unit(0.5, "lines"),
    text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.line = element_line(color = "#000000", linewidth = 12/14),
    panel.grid.major.y = element_line(color = "#cfcfcf50"),
    legend.box.margin = margin(t = 10, b = 10),
    axis.title.x = element_blank())

FindMarkers(
  seu,
  ident.1 = "18hr Luc",
  ident.2 = "18hr Tat",
  logfc.threshold = 0.2,
  features = proApotosisGenes,
  min.pct = 0.2,
  slot = "data",
  test.use = "wilcox")


################################################################################
# fig2...
################################################################################
layout <- c(
  area(1, 1, 2, 3), # umap
  area(1, 4, 3, 9), # tracks
  area(3, 1, 4, 3), # vln plot
  area(5, 1, 6, 3), # deg volcano 72 hrs
  area(4, 4, 6, 9) # dot plot
)


final_p <- patchwork::free(figA, side = "trb") + 
  patchwork::free(fig2TrackFinal, side = "lr") +
  patchwork::free(fig2C, side = "trb") +
  patchwork::free(fig2d, side = "trb") +
  figDotPlotApop +
  plot_layout(design = layout) + 
  plot_annotation(tag_levels = list(c("A", "B", "", "C", "D"))) &
  theme(plot.tag = element_text(size = 10, color = "#000000"))

savePlot(plot = final_p, fn = "fig2", devices = c("png", "pdf"), gwidth = 8, gheight = 5)

