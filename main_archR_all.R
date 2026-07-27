library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(ArchR)
library(patchwork)

set.seed(21)
addArchRThreads(8)

samples_late <- list.files("data/cellranger_late_out/")
fragmentFiles_late <- glue("data/cellranger_late_out/{samples_late}/outs/atac_fragments.tsv.gz")

samples_early <- list.files("data/cellranger_early_out/")
fragmentFiles_early <- glue("data/cellranger_early_out/{samples_early}/outs/atac_fragments.tsv.gz")

geneAnnot_10.6 <- readRDS("make_bsgenome/rds/10-6_geneAnnot.rds")
genomeAnnot_10.6 <- readRDS("make_bsgenome/rds/10-6_genomeAnnot.rds")
chrNamesDiscard_10.6 <- readRDS("make_bsgenome/discard_10.6.rds")

arrowfiles <- createArrowFiles(
  inputFiles = c(fragmentFiles_late, fragmentFiles_early),
  sampleNames = c(samples_late, samples_early),
  outputNames = c(paste0("origArrow/", samples_late), paste0("origArrow/", samples_early)),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard_10.6,
  geneAnnotation = geneAnnot_10.6,
  genomeAnnotation = genomeAnnot_10.6,
  subThreading = FALSE,
  threads = 8,
  force = FALSE,
  addGeneScoreMat = TRUE)

# setting it lower to reduce RAM usage
addArchRThreads(threads = 12)

initProjDir <- "dogma_init_all_v1"
qcFiltProjDir <- "dogma_qcFilt_all_v1"
if (!dir.exists(qcFiltProjDir)) {
  proj <- ArchRProject(
    ArrowFiles = arrowfiles,
    outputDirectory = initProjDir,
    geneAnnotation = geneAnnot_10.6,
    genomeAnnotation = genomeAnnot_10.6,
    copyArrows = TRUE,
    showLogo = FALSE)
  
  # filter using AMULET
  multCells <- lapply(c(samples_late, samples_early), function(x) {
    if (x %in% c("D2_10-6_Tat_2", "D3_10-6_Luc_2", "D3_10-6_Tat_2", "NA_10-6_PmaIono")) {
      amuletOutDir <- "amulet_late_out"
    } else {
      amuletOutDir <- "amulet_early_out"
    }
    
    
    amuletOut <- unlist(read.table(paste0("data/", amuletOutDir, "/", x, "/MultipletBarcodes_01.txt"), header = FALSE))
    amuletOut <- paste0(x, "#", amuletOut)
    
    return(amuletOut)
  })
  multCells <- unlist(multCells)
  
  singlets <- setdiff(proj$cellNames, multCells)
  projQcFilter <- proj[singlets, ]
  
  ##############################################################################
  # add RNA into ArchRProj
  ##############################################################################
  seRNA <- import10xFeatureMatrix(
    input = c(
      paste0("data/cellranger_late_out/", samples_late, "/outs/filtered_feature_bc_matrix.h5"),
      paste0("data/cellranger_early_out/", samples_early, "/outs/filtered_feature_bc_matrix.h5")),
    names = c(samples_late, samples_early)
  )
  
  # intersection of barcodes in atac and rna
  multiBcs <- intersect(projQcFilter$cellNames, colnames(seRNA))
  projQcFilter <- subsetArchRProject(
    cells = multiBcs,
    ArchRProj = projQcFilter,
    outputDirectory = "dogma_subsetMatched_all_v1",
  )
  
  # this step already adds in Mito and Ribo automatically :)
  projQcFilter <- addGeneExpressionMatrix(
    input = projQcFilter,
    seRNA = seRNA,
    force = TRUE,
    strictMatch = TRUE,
    excludeChr = c(chrNamesDiscard_10.6, "chrNA"))
  
  
  tmp <- getCellColData(ArchRProj = projQcFilter, select = c("TSSEnrichment", "Sample", "Gex_nGenes", "Gex_nUMI", "Gex_MitoRatio"))
  
  ggplot(as.data.frame(tmp), aes(x = Sample, y = TSSEnrichment)) +
    geom_violin() +
    # geom_hline(yintercept = 500) +
    # geom_hline(yintercept = 6000) +
    theme_classic()
  
  
  # qc filter
  idxPass <- BiocGenerics::which(
    projQcFilter$TSSEnrichment >= 10 &
      log10(projQcFilter$nFrags) >= 3.5 &
      projQcFilter$Gex_nGenes > 500 &
      projQcFilter$Gex_nGenes < 6000 &
      projQcFilter$Gex_MitoRatio < 0.1)
  
  cellsPass <- projQcFilter$cellNames[idxPass]
  projQcFilter <- projQcFilter[cellsPass, ]
  
  rm(tmp)
  gc()
  
  # 
  # ##############################################################################
  # # add information to ArchRProj and metadata
  # ##############################################################################

  
  # save archr project
  projUnhashed <- saveArchRProject(
    ArchRProj = projQcFilter,
    outputDirectory = qcFiltProjDir,
    load = TRUE
  )
  
} else {
  projUnhashed <- loadArchRProject(
    path = qcFiltProjDir,
    showLogo = FALSE
  )
}
