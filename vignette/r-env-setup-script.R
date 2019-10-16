#' 
#' This script sets up the R environment in the 
#' Google Colab virtual maschine. Mainly it 
#' downloads the precompiled R packages and 
#' installs needed Ubuntu packages.
#'

#' R package cache (precompiled packages based on colab env)
message("Download R package cache")
system("wget --continue https://drive.google.com/open?id=11FiOzQf7YFESIhYG4tFaL6odbnXDAHop")

#' Unpack cache locally
message("Unzipping R package cache")
system("tar -xzf CTD_0.0.0.9000.tar.gz -C /", intern=TRUE)

#' Set correct library path
.libPaths(Sys.getenv("R_LIBS_USER"))


# Options to make plots smaller
options(repr.plot.width=4, repr.plot.height=4)
