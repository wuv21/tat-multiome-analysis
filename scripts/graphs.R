################################################################################
# COLOR SCHEMES
################################################################################
COLORS <- list(
)

BASEPTAXISFONTSIZE <- 7
BASEPTFONTSIZE <- 7
BASEFONTSIZE <- BASEPTFONTSIZE / ggplot2:::.pt

pValSymnum <- function(x, showNs = TRUE) {
  tmp <- sapply(x, function(y) {
    if (is.na(y)) {
      return(NA)
    }
    
    if (y < 0.001) {
      return("***")
    } else if (y >= 0.001 & y < 0.01) {
      return("**")
    } else if (y >= 0.01 & y < 0.05) {
      return("*")
    } else {
      return(ifelse(showNs, "ns", ""))
    }
  })
  
  return(tmp)
}

################################################################################
# Generate UMAP df from ArchRProject
################################################################################
generateUmapDfFromArchR <- function(
  proj,
  cluster = "Clusters",
  customSort = FALSE,
  donorColumn = NULL,
  colorLabelCluster = NULL,
  embedding = "Harmony_UMAP_Combined_v2") {
  
  umapFromArchr <- getEmbedding(proj, embedding = embedding)
  
  # check to make sure embedding order is same as project order
  if (sum(rownames(umapFromArchr) == proj$cellNames) != length(proj$cellNames)) {
    message("Embedding rownames from getEmbedding() not in the same order as ArchR project's cell names. Will fix now.")
    
    # sort embedding df by proj$cellNames
    umapFromArchr <- umapFromArchr[proj$cellNames, ]
  }
  
  df <- data.frame(
    cbc = rownames(umapFromArchr),
    x = umapFromArchr[, 1],
    y = umapFromArchr[, 2],
    sample = proj$Sample)
  

  if (is.vector(cluster) && length(cluster) > 1) {
    tmp <- getCellColData(proj, select = cluster)
    tmp <- tidyr::unite(as.data.frame(tmp), col = "tmp", sep = ": ")
    df$cluster <- tmp$tmp
    
  } else {
    df$cluster <- getCellColData(proj, select = cluster)[, 1]
    
    if (customSort) {
      df$cluster <- factor(df$cluster, levels = customSortAnnotation(df$cluster))
    }
  }
  
  if (!is.null(colorLabelCluster)) {
    df$colorLabelCluster <- getCellColData(proj, select = colorLabelCluster)[, 1]
  }
  
  if (!is.null(donorColumn)) {
    df$donor <- getCellColData(proj, select = donorColumn)[, 1]
  }
  
  return(df)
}


generateNonUmapDfFromArchR <- function(
  proj,
  groupingVar = NULL,
  groupingSort = NULL,
  splittingVar = NULL,
  splittingSort = NULL,
  fillVar = NULL,
  colorVar = NULL
  ) {
  
  df <- data.frame(cbc = proj$cellNames)
  
  if (is.vector(groupingVar) && length(groupingVar) > 1) {
    tmp <- getCellColData(proj, select = groupingVar)
    tmp <- tidyr::unite(as.data.frame(tmp), col = "tmp", sep = ": ")
    
    stopifnot(sum(rownames(tmp) == df$cbc) == nrow(df))
    
    df$groupingVar <- tmp$tmp
    
  } else {
    tmp <- getCellColData(proj, select = groupingVar)
    stopifnot(sum(rownames(tmp) == df$cbc) == nrow(df))
    
    df$groupingVar <- tmp[, 1]
    
    if (!is.null(groupingSort)) {
      df$groupingVar <- factor(df$groupingVar, levels = groupingSort)
    }
  }
  
  if (!is.null(splittingVar)) {
    df$splittingVar <- getCellColData(proj, select = splittingVar)[, 1]
  }
  
  if (!is.null(fillVar)) {
    df$fillVar <- getCellColData(proj, select = fillVar)[, 1]
  }
  
  if (!is.null(colorVar)) {
    df$colorVar <- getCellColData(proj, select = colorVar)[, 1]
  }
  
  return(df)
}

################################################################################
# UMAP plot themes
################################################################################
umapTheme <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.title = element_blank(),
  axis.title = element_text(size = BASEPTFONTSIZE),
  axis.text = element_text(size = BASEPTAXISFONTSIZE),
  legend.spacing.x = unit(BASEPTFONTSIZE / 2, 'points'),
  # legend.spacing.y = unit(BASEPTFONTSIZE / 4, 'points'),
  plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5),
  plot.title.position = "panel",
  legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA),
  panel.background = element_rect(fill = "transparent", colour = NA),
  plot.background = element_rect(fill = "transparent", colour = NA))


################################################################################
# Plot umap
################################################################################
plotUmap <- function(
  proj,
  fn,
  width = 4,
  height = 3.5,
  legendNRow = 5,
  returnPlot = FALSE,
  cellFilter = NULL,
  title = NULL,
  devices = c("png"),
  ggtheme = umapTheme,
  colorBy = "Clusters",
  colorLabelBy = colorBy,
  colorLabel = colorBy,
  colorScheme = NULL,
  rasterize = TRUE,
  bringToTop = FALSE,
  customClusterSort = FALSE,
  propInLegend = FALSE,
  propDigits = 1,
  embedding = "UMAP") {
  
  df <- generateUmapDfFromArchR(
    proj,
    cluster = colorBy,
    embedding = embedding,
    colorLabelCluster = colorLabelBy,
    customSort = customClusterSort)
  
  if (!is.null(cellFilter)) {
    df <- df[df$cbc %in% cellFilter, ]
  }
  
  
  if (bringToTop == -1) {
    df <- df %>%
      arrange(desc(cluster))
  } else if (bringToTop) {
    df <- df %>%
      arrange(cluster)
  }
  
  if (length(colorBy) == 1 && colorBy == "haystackOut") {
    df <- df %>%
      dplyr::mutate(cluster = ifelse(cluster, "HIV+", "HIV-"))
  }
  
  if (propInLegend) {
    # need to carry over the mapped color scheme
    
    if (!is.null(colorScheme) & sum(class(colorScheme) == "Scale") == 0) {
      colorSchemeTmp <- df %>%
        group_by(cluster) %>%
        dplyr::mutate(cluster2 = glue("{cluster} ({round(n() / nrow(.) * 100, digits = propDigits)}%)")) %>%
        select(cluster, cluster2) %>%
        distinct()
      
      
      colorSchemeTmp <- as.data.frame(colorSchemeTmp)
      
      rownames(colorSchemeTmp) <- colorSchemeTmp$cluster
      colorSchemeTmp$color <- colorScheme[colorSchemeTmp$cluster]
      
      colorScheme <- colorSchemeTmp$color
      names(colorScheme) <- colorSchemeTmp$cluster2
      
      colorScheme <- scale_color_manual(values = colorScheme)
    }
    
    df <- df %>%
      group_by(cluster) %>%
      dplyr::mutate(cluster = glue("{cluster} ({round(n() / nrow(.) * 100, digits = propDigits)}%)"))
  }
  
  if (customClusterSort) {
    # TODO fix because cluster has more info...
    
    df$cluster <- factor(df$cluster, levels = customSortAnnotation(df$cluster))
  }
  
  p1 <- ggplot(df, aes(x = x, y = y))
  
  if (rasterize) {
    p1 <- p1 + rasterize(geom_point(alpha = 0.6, aes(color = cluster), size = 0.25), dpi = 300)
  } else {
    p1 <- p1 + geom_point(alpha = 0.8, aes(color = cluster), size = 0.25)
  }
  
  p1 <- p1 +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      color = colorLabel) +
    theme_classic() +
    ggtheme
  
  if (!is.null(title)) {
    p1 <- p1 + labs(title = title)
  }
  
  if (!is.numeric(df$cluster)) {
    p1 <- p1 +
      guides(colour = guide_legend(override.aes = list(size = 4), nrow = legendNRow))
  }
  
  if (is.null(colorScheme) && (colorBy == "haystackOut" | colorBy == "hivRNA")) {
    p1 <- p1 + scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR))
  } else if (!is.null(colorScheme) && sum(class(colorScheme) == "Scale") > 0) {
    p1 <- p1 + colorScheme
  } else if(!is.null(colorScheme) && sum(class(colorScheme) == "Scale") == 0) {
    p1 <- p1 + scale_color_manual(values = colorScheme)
  }
  
  if (!is.null(colorBy) & !is.null(colorLabelBy)) {
    clusterLabelUmapPos <- df %>% 
      group_by(colorLabelCluster) %>% 
      dplyr::summarize(
        topX = quantile(x, c(.6)),
        topY = quantile(y, c(.6)),
        bottomX = quantile(x, c(.4)),
        bottomY = quantile(y, c(.4)),
        x = (topX + bottomX) / 2,
        y = (topY + bottomY) / 2)
    
    p1 <- p1 + 
      ggrepel::geom_label_repel(data = clusterLabelUmapPos,
        aes(x = x, y = y, label = colorLabelCluster),
        label.size = 0.05,
        force = 300,
        max.time = 5,
        max.iter = 25000,
        size = BASEFONTSIZE,
        seed = 21,
        min.segment.length = 0.25,
        segment.color = "#BBBBBB")
  }
  
  if (!is.null(fn)) {
    savePlot(plot = p1, fn = fn, devices = devices, gheight = height, gwidth = width)
  }
  
  if (returnPlot) {
    return(p1)
  }
}


cleanUpTrackAndSave <- function(
  archrTrack,
  fn,
  gheight = 2.25,
  gwidth = 4,
  devices = c("png", "rds")
) {
  # plot the raw figure
  p <- wrap_plots(archrTrack, ncol = 1, heights = c(3,2,1))
  savePlot(p, fn = paste0(fn, "_orig"), devices = "png", gheight = gheight, gwidth = gwidth)
  
  # clean up theme
  cleanUpTheme <- theme(
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    text = element_text(family = "Arial", size = 6),
    axis.text = element_text(family = "Arial", size = 6),
    panel.border = element_blank())
  
  # clean up
  archrTrack$bulktrack <- archrTrack$bulktrack +
    labs(y = "Normalized signal") +
    cleanUpTheme
  
  archrTrack$sctrack <- archrTrack$sctrack +
    labs(y = "Binarized signal") +
    cleanUpTheme +
    scale_color_manual(values = c(HIVNEGCOLOR, HIVPOSCOLOR))
  
  archrTrack$genetrack <- archrTrack$genetrack +
    cleanUpTheme +
    labs(y = "Gene") +
    scale_color_manual(values = c("blue", "orange"))
  
  savePlot(archrTrack, customSavePlot = wrap_plots(archrTrack, ncol = 1, heights = c(2,1.5,1.5)),
    fn = fn, devices = devices, gheight = gheight, gwidth = gwidth)
}
