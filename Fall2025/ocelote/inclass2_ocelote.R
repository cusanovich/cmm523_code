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