remotes::install_version("systemfonts",version="1.2.3",
    repos = "http://cran.us.r-project.org",upgrade=F)
remotes::install_version("textshaping",version="1.0.3",
    repos = "http://cran.us.r-project.org",upgrade=F)
remotes::install_version("ragg",version="1.4.0",
    repos = "http://cran.us.r-project.org",upgrade=F)
install.packages(c("harmony","tidyverse","viridis","pheatmap"),
    repos = "http://cran.us.r-project.org")
remotes::install_github("SingleR-inc/SingleR@RELEASE_3_13",upgrade=F)
remotes::install_version("matrixStats", version="1.0.0",
    repos = "http://cran.us.r-project.org",upgrade=F)
remotes::install_github("satijalab/seurat-data",upgrade=F)

install.packages("BiocManager", repos = "https://cloud.r-project.org")
BiocManager::install(c("scRNAseq", "celldex", "scran"), update=F)
#And one last package for fun
remotes::install_version("dbplyr", version="2.3.4",
    repos = "https://cloud.r-project.org",upgrade=F)
