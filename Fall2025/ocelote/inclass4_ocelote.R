remotes::install_version("Cairo", version="1.6-2",
    repos = "http://cran.us.r-project.org")
remotes::install_version("terra", version="1.8-70",
    repos = "http://cran.us.r-project.org")
BiocManager::install(c('Rhdf5lib', 'rhdf5filters', 'lme4', 'batchelor', 'HDF5Array',
                       'ggrastr'), update=F)
install.packages(c("statmod","R.utils"), repos = "http://cran.us.r-project.org")
remotes::install_github(c("r-spatial/sf", "cole-trapnell-lab/monocle3@v1.3.1",
                           "satijalab/seurat-wrappers"),upgrade=F)
