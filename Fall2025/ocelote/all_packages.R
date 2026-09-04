#IC2
install.packages("remotes", repos="http://cran.us.r-project.org")
#These ones use remotes to install an older version of a package
remotes::install_version("Matrix", version = "1.6-5",
    repos = "http://cran.us.r-project.org")
remotes::install_version("curl", version = "5.2.3",
repos = "http://cran.us.r-project.org")
remotes::install_version("openssl", version = "2.3.2",
    repos = "http://cran.us.r-project.org")
remotes::install_version("ggplot2", version = "3.5.2",
    repos = "http://cran.us.r-project.org")
install.packages(c("png","Seurat"), repos="http://cran.us.r-project.org",
	upgrade=F)
#This last one uses remotes to install a package from GitHub
remotes::install_github("mojaveazure/seurat-disk",upgrade=F)
remotes::install_github("immunogenomics/presto",upgrade=F)

#IC3
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
    
#IC4
remotes::install_version("Cairo", version="1.6-2",
    repos = "http://cran.us.r-project.org")
remotes::install_version("terra", version="1.8-70",
    repos = "http://cran.us.r-project.org")
BiocManager::install(c('Rhdf5lib', 'rhdf5filters', 'lme4', 'batchelor', 'HDF5Array',
                       'ggrastr'), update=F)
install.packages(c("statmod","R.utils"), repos = "http://cran.us.r-project.org")
remotes::install_github(c("r-spatial/sf", "cole-trapnell-lab/monocle3@v1.3.1",
                           "satijalab/seurat-wrappers"),upgrade=F)

#IC6
BiocManager::install("biovizBase", update=F)
install.packages(c("patchwork","Signac"),repos="http://cran.us.r-project.org")

#IC7
BiocManager::install(c("EnsDb.Hsapiens.v86","chromVAR",
                  "JASPAR2020","motifmatchr",
                  "BSgenome.Hsapiens.UCSC.hg38",
                       "glmGamPoi"), update=F)
remotes::install_github("immunogenomics/presto",upgrade=F)
