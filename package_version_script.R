# Package Versions Loaded/Installed
###################################

cat(
  "Bioconductor:", as.character(packageVersion("BiocManager")),
  "\ndplyr:", as.character(packageVersion("dplyr")),
  "\nggplot2:", as.character(packageVersion("ggplot2")),
  "\npatchwork:", as.character(packageVersion("patchwork")),
  "\nknitr:", as.character(packageVersion("knitr")),
  "\nstringr:", as.character(packageVersion("stringr")),
  "\nkableExtra:", as.character(packageVersion("kableExtra")),
  "\nmosaic:", as.character(packageVersion("mosaic"))
  )

##Quick environment sanity check
cat("\n=== Session Info ===\n")
cat("R:", R.version.string, "\n")
cat("Bioc:", as.character(BiocManager::version()), "\n")
cat("Library paths:\n")
print(.libPaths())
cat("\ncuratedMetagenomicData:", as.character(packageVersion("curatedMetagenomicData")), "\n\n")