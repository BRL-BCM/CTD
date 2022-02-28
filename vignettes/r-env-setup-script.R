#' 
#' This script sets up the R environment in the 
#' Google Colab virtual maschine. Mainly it 
#' downloads the precompiled R packages and 
#' installs needed Ubuntu packages.
#'

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
downloadGDriveFile(out="r_binaries.tar.gz", id="166AfpzFwfWU8gtu6su120RYhmETgSZrn")

system("ls")

#' Unpack cache locally
message("Unzipping R package cache")
system("tar -xzf r_binaries.tar.gz -C / ", intern=TRUE)

#' Set correct library path
.libPaths(Sys.getenv("R_LIBS_USER"))


# Options to make plots smaller
options(repr.plot.width=10, repr.plot.height=10)
