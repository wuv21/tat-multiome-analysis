################################################################################
# generic functions
################################################################################
savePlot <- function(plot, fn, devices, gheight, gwidth, rdsPlot = NULL, scale = 1, customSavePlot = NULL) {
  if (!is.vector(devices)) {
    devices <- c(devices)
  }
  
  for (d in devices) {
    gfn <- glue("outs/{d}/{fn}.{d}")
    
    if (d == "rds" & !is.null(rdsPlot)) {
      saveRDS(rdsPlot, gfn)
    } else if (d == "rds") {
      saveRDS(plot, gfn)
    } else if (!is.null(customSavePlot)) {
      ggsave(gfn, plot = customSavePlot, dpi = "retina", device = d, width = gwidth, height = gheight, scale = scale)
    } else {
      ggsave(gfn, plot = plot, dpi = "retina", device = d, width = gwidth, height = gheight, scale = scale)  
    }
  }
}


customSortAnnotation <- function(x, reverse = FALSE, ignoreClusterPrefix = TRUE) {
  priority <- c("MAIT", "Tfh", "Naive", "Tcm", "Tem", "Activated", "Treg")
  
  origNames <- x
  if (ignoreClusterPrefix) {
    x <- gsub("C\\d+: ", "", x)
    
    mapNames <- origNames
    names(mapNames) <- x
  }
  
  x <- sort(unique(x))
  
  cd4Indices <- grepl("^CD4", x)
  
  oldCd4 <- x[cd4Indices]
  newCd4 <- c()
  
  for (p in priority) {
    pIndices <- grepl(p, oldCd4)
    newCd4 <- append(newCd4, oldCd4[pIndices])
    oldCd4 <- oldCd4[!pIndices]
  }
  
  if (length(oldCd4) > 0) {
    newCd4 <- append(newCd4, oldCd4)
  }
  
  finalSort <- c(newCd4, x[!cd4Indices])
  
  if (ignoreClusterPrefix) {
    finalSort <- mapNames[finalSort]
  }
  
  return(finalSort)
}