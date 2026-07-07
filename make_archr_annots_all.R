library(GenomicFeatures)
library(ArchR)
library(dplyr)
library(BSgenome.Hg38.Jlat106Integrated)

txdb_new_10.6 <- makeTxDbFromGFF(
  file = "make_bsgenome/gtf/hg38_jlat10-6_integrated.renamed.gtf",
  format = "auto",
  organism = "Homo sapiens",
  circ_seqs = "chrM"
)


prepareGeneGRanges <- function(gtf_fn, txdb) {
  geneSymbols <- read.table(gtf_fn, sep = "\t")
  geneSymbolsClean <- geneSymbols %>%
    filter(V3 == "gene") %>%
    mutate(
      gene_id = stringr::str_match(V9, "(gene_id )(ENSG\\d+)")[, 3],
      symbol = stringr::str_match(V9, "(gene_name )([^\\s;]+)")[, 3]
    ) %>%
    distinct()
  
  # note that hiv gene is given an arbitrary ensembl ID of ENSG99999999999
  geneSymbolDict <- geneSymbolsClean$symbol
  names(geneSymbolDict) <- geneSymbolsClean$gene_id
  
  actualGenes <- genes(txdb)
  finalSymbols <- geneSymbolDict[actualGenes$gene_id]
  names(finalSymbols) <- NULL
  
  # https://github.com/GreenleafLab/ArchR/issues/422
  actualGenes <- GRanges(symbol = finalSymbols, actualGenes)
  
  return(actualGenes)
}

keepBSgenomeSequences <- function(genome, seqnames) {
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  
  return(genome)
}


cleanGenome <- function(bsgenome, discardfn) {
  chrNames <- seqnames(bsgenome)
  chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM)", chrNames)]
  chrNamesKeep <- paste0("chr", c(1:22, "X", "Y"))
  
  genome <- bsgenome
  genome <- keepBSgenomeSequences(genome, chrNamesKeep)
  
  saveRDS(chrNamesDiscard, discardfn)
  return(list("genome" = genome, "discard" = chrNamesDiscard))
}

prepareDenylist <- function(bed) {
  denylist <- read.table(bed, header = FALSE, sep = "\t")
  denylistGr <- GRanges(
    seqnames = denylist$V1,
    ranges = IRanges(
      start = denylist$V2,
      end = denylist$V3))
  
  return(denylistGr)
}

prepareChromSizes <- function(genome) {
  new_seqinfo <- seqinfo(genome)
  chromSizes <- GRanges(
    seqnames = seqnames(new_seqinfo),
    ranges = IRanges(
      start = 1,
      end = seqlengths(new_seqinfo)
    )
  )
  
  return(chromSizes)
}

saveAllRds <- function(txdb, genes, genome, chromsizes, denylist, discard, name) {
  geneAnnot <- createGeneAnnotation(
    genes = genes,
    exons = exons(txdb),
    TSS = promoters(txdb, upstream = 0, downstream = 1),
    annoStyle = "ENSEMBL")
  
  saveRDS(geneAnnot, paste0("make_bsgenome/rds/", name, "_geneAnnot.rds"))
  
  genomeAnnot <- createGenomeAnnotation(
    genome = genome,
    chromSizes = chromsizes,
    blacklist = denylist,
    filter = TRUE,
    filterChr = discard
  )
  
  saveRDS(genomeAnnot, paste0("make_bsgenome/rds/", name, "_genomeAnnot.rds"))
}


genes_10.6 <- prepareGeneGRanges("make_bsgenome/gtf/hg38_jlat10-6_integrated.renamed.gtf", txdb_new_10.6)
genome_10.6 <- cleanGenome(BSgenome.Hg38.Jlat106Integrated, "make_bsgenome/discard_10.6.rds")
denylist_10.6 <- prepareDenylist("make_bsgenome/denylist/10.6.bed")
chromsizes_10.6 <- prepareChromSizes(genome_10.6$genome)
saveAllRds(txdb_new_10.6, genes_10.6, genome_10.6$genome, chromsizes_10.6, denylist_10.6, genome_10.6$discard, name = "10-6")


