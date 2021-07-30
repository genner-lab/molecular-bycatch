#!/usr/bin/env Rscript

# libs
suppressMessages({
    library("here")
    library("rdryad")
})

# get dryad data
files <- dryad_download(dois="10.5061/dryad.b8f6s44")

# copy
file.copy(c(files$"10.5061/dryad.b8f6s44"[1],files$"10.5061/dryad.b8f6s44"[2]), here("temp/data"))

# report
writeLines("\nData obtained from Dryad\n")
