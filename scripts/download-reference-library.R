#!/usr/bin/env Rscript

# libs
suppressMessages({
    library("here")
    library("tidyverse")
    library("ape")
})

# load remote references and scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")

# get version
gb.version <- reflib.orig %>% filter(source=="GENBANK") %>% distinct(genbankVersion) %>% pull(genbankVersion)

# write out combined
reflib.orig %>% 
    write_csv(file=here(paste0("meta-fish-pipe/assets/meta-fish-lib-v",gb.version,".csv")))

# report
writeLines(paste0("\nReference library saved to 'meta-fish-pipe/assets/meta-fish-lib-v",gb.version,".csv'\n"))
