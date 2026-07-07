suppressPackageStartupMessages(library(BSgenome))

forgeBSgenomeDataPkg(
  x = "make_bsgenome/BSgenome.Hg38.Jlat106Integrated.seed",
  seqs_srcdir = "make_bsgenome/seqs",
  destdir = "make_bsgenome",
  verbose = TRUE)








# /opt/R/4.1.1/bin/R CMD build BSgenome.Hg38.Jlat106Integrated
# /opt/R/4.1.1/bin/R CMD check BSgenome.Hg38.Jlat106Integrated_1.0.0.tar.gz
# install.packages("make_bsgenome/BSgenome.Hg38.Jlat106Integrated_1.0.0.tar.gz", repos = NULL, type = "source")