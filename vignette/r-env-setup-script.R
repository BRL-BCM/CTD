#' 
#' This script sets up the R environment in the 
#' Google Colab virtual maschine. Mainly it 
#' downloads the precompiled R packages and 
#' installs needed Ubuntu packages.
#'

#' Function to download bigger files from google drive automatically.
downloadGDriveFile <- function(id, out){
    system(paste0(
        'wget --continue --load-cookies /tmp/cookies.txt ',
            '"https://docs.google.com/uc?export=download&confirm=$(wget ',
            '--quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate ',
            "'https://docs.google.com/uc?export=download&id=", id, "' -O- | ",
            "sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=", id, '" ',
            "-O ", out, " && rm -rf /tmp/cookies.txt"), intern=TRUE)
}

#' R package cache (precompiled packages based on colab env)
message("Download R package cache")
downloadGDriveFile(out="CTD_0.0.0.9000.tar.gz", id="11FiOzQf7YFESIhYG4tFaL6odbnXDAHop")

#' Unpack cache locally
message("Unzipping R package cache")
system("tar -xzf CTD_0.0.0.9000.tar.gz -C /", intern=TRUE)

#' Set correct library path
.libPaths(Sys.getenv("R_LIBS_USER"))


# Options to make plots smaller
options(repr.plot.width=4, repr.plot.height=4)
