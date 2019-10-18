#' 
#' This script sets up the R environment in the 
#' Google Colab virtual maschine. Mainly it 
#' downloads the precompiled R packages and 
#' installs needed Ubuntu packages.
#'

#' R package cache (precompiled packages based on colab env)
message("Download R package caches")
system("wget --continue https://drive.google.com/uc?export=download&id=1HjI0ee2ImcQvL9OiJyKIwkARqLTkDHe6")

#' Unpack cache locally
message("Unzipping R package cache")
system("tar -xzf r_binaries.tar.gz -C /", intern=TRUE)

#' Set correct library path
.libPaths(Sys.getenv("R_LIBS_USER"))


# Options to make plots smaller
options(repr.plot.width=4, repr.plot.height=4)
