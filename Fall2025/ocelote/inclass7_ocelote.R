BiocManager::install(c("EnsDb.Hsapiens.v86","chromVAR",
                  "JASPAR2020","motifmatchr",
                  "BSgenome.Hsapiens.UCSC.hg38",
                       "glmGamPoi"), update=F)
remotes::install_github("immunogenomics/presto",upgrade=F)
