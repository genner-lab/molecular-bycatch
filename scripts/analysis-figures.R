#!/usr/bin/env Rscript

#' ---
#' Title: "Estuarine molecular bycatch as a landscape-wide biomonitoring tool"
#' Author: "Stefano Mariani, Lynsey R. Harper, Rupert A. Collins, Charles Baillie, Owen S. Wangensteen, Allan D. McDevitt, Morton Heddell-Cowie and Martin J. Genner"
#' Date: "9th August 2021"
#' ---
#' 
#' Three estuaries across the UK were sampled in 2016 and 2017:
#' - Esk in Autumn 2016 and Spring 2017
#' - Tees in Autumn 2016, Spring 2017, and Autumn 2017
#' - Tweed in Spring 2017 and Autumn 2017
#' 
#' eDNA samples were screened with the MiFish (Miya et al. 2015) or 
#' Tele02 (Taberlet et al. 2018) primers to target the 12S region in
#' fish and examine seasonal variation in fish communities for each 
#' estuary. However, 12S primers generally detect other vertebrates
#' as well as fish as this region is highly conserved across taxa.
#' 
#' A number of bird and mammal taxa were detected alongside fish, which
#' emphasises the potential of estuaries as collectors for catchment-scale
#' biodiversity, similar to ponds (Harper et al. 2019) and rivers (Sales
#' et al. 2020).
#' 
#' We will compare alpha and beta diversity between estuaries, and 
#' analyse seasonal variation in bird and mammal assemblages.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. 
#' 

## load
source(here::here("scripts/check-resids-function.R"))
library("here")
library("reshape2")
library("ggpubr")
library("scales")
library("car")
library("MASS")
library("FSA")
library("lawstat")
library("vegan")
library("betapart")
library("iNEXT")
library("tidyverse")

#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()

#' create output dir
dir.create(here("results/figures"),recursive=TRUE)

#' ---
#' 
#' ## 1) Metabarcoding dataset
#' 
#' The output from OBItools is not particularly suited to applying
#' statistical analysis in R. Files must be manipulated for downstream 
#' analyses.
#' 

## Import metabarcoding data
EA.raw <- read.csv(here("assets/ea-estuary-raw-data.csv"), header = TRUE)

## Remove redundant columns
EA.reduced <- EA.raw %>% 
        dplyr::select(-c(phylum, order, family, genus, nASVs, 
                         sintaxBSabund, sintaxBSmean, nucsAbund, 
                         totalReads))

## Store taxonomy information in separate dataframe
taxonomy <- EA.reduced %>% dplyr::select(c(class, species))

## Remove class column
EA.reduced <- EA.reduced %>% dplyr::select(-class)

## Transpose dataframe while maintaining wide format
EA.wide <- setNames(data.frame(t(EA.reduced[,-1])), EA.reduced[,1])
EA.wide <- rownames_to_column(EA.wide, "sample_ID")

## Split sample_ID column into new columns
EA.temp <- EA.wide %>% 
        separate(sample_ID, c("primer_set","library","site","bio_rep",
                              "pcr_rep","serial"), "_")

## Split site column into new columns then split newly created location
## column into a further two columns named site and number. Add leading
## 0s to numbers less than 10 in the number column. Recombine the site
## and number columns into a single column called location. Finally, 
## remove uninformative columns.
EA.master <- EA.temp %>% 
        separate(site, c("sample_type","year","month","day","estuary","location"), "[.]") %>%
        separate(location, c("site", "number"), "(?<=[A-Z])(?=[0-9])") %>%
        mutate(number = str_pad(number, 2, side = "left", pad = 0)) %>%
        mutate(number = replace_na(number, "")) %>%
        unite("location", c(site, number), sep = "", remove = FALSE) %>%
        dplyr::select(-c(day, serial, site, number)) %>%
        mutate(sample_type = replace(sample_type, bio_rep == "FieldBlank", "FieldBlank"))

## Add column containing unique identifier for each sample to EA.master
EA.master <- EA.master %>%
        mutate(identifier = 1:n()) %>%
        relocate(identifier, .after = pcr_rep)

## Tidy the master data:
## Rename primers as "Teleo02" and "MiFish-U"
EA.master$primer_set <- gsub("tele02", "Tele02", EA.master$primer_set)
EA.master$primer_set <- gsub("mifish.u", "MiFish-U", EA.master$primer_set)

## Rename libraries as "Library 1" and "Library 2"
EA.master$library <- gsub("lib1", "Library 1", EA.master$library)
EA.master$library <- gsub("lib2", "Library 2", EA.master$library)

## Replace "EA" with "eDNA" in sample_type column
EA.master$sample_type <- gsub("EA", "eDNA", EA.master$sample_type)

## Create new column with season that eDNA samples were collected. Move 
## new column to come after year column in master dataframe. Create new 
## column containing season and year, then remove year, season and month
## columns.
EA.master <- EA.master %>%
        mutate(season = if_else(month %in% c("05","06"), "Spring", 
                                if_else(month %in% c("09","10"), "Autumn",
                                        ""))) %>% 
        relocate(season, .after = year) %>%
        unite("season_year", c(season, year), sep = " ", remove = FALSE) %>%
        dplyr::select(-c(year, season, month))

## Combine newly created sample metadata with original sample IDs for
## downsteam diversity analyses using vegan
metadata <- data.frame(EA.wide$sample_ID, EA.master[,1:8])
metadata <- rename(metadata, sample_ID = EA.wide.sample_ID)

## Convert master dataframe to long format by melting
## Species detected and read counts will be condensed into two columns
EA.dat <- melt(EA.master, id = c("primer_set","library","sample_type", 
                                 "season_year","estuary","location", 
                                 "bio_rep","pcr_rep","identifier"))

## Rename new columns and only keep taxa with read counts for each 
## sample
EA.dat <- EA.dat %>% 
        rename(species = variable, reads = value) %>%
        filter(reads > 0) %>%
        droplevels()

## Merge taxonomy with EA.dat
EA.dat <- merge(EA.dat, taxonomy, by = "species", all.x = TRUE)

## Reorder and rename columns
EA.dat <- EA.dat %>% 
        relocate(c(class, species), .before = reads) %>%
        rename(group = class)

## Replace class names with common group names
EA.dat$group <- gsub("Actinopterygii", "Ray-finned fishes", EA.dat$group)
EA.dat$group <- gsub("Amphibia", "Amphibians", EA.dat$group)
EA.dat$group <- gsub("Aves", "Birds", EA.dat$group)
EA.dat$group <- gsub("Cephalaspidomorphi", "Jawless fishes", EA.dat$group)
EA.dat$group <- gsub("Elasmobranchii", "Elasmobranchs", EA.dat$group)
EA.dat$group <- gsub("Mammalia", "Mammals", EA.dat$group)
EA.dat$group <- gsub("Anthozoa|Asteroidea|Bivalvia|Gastropoda|Insecta|Malacostraca|Maxillopoda", "Invertebrates", EA.dat$group)

## Remove red porgy (Pagrus pagrus) as this is a marine fish that was 
## detected in seafood mislabelling projects, and witch (Glyptocephalus 
## cynoglossus) as this species was only detected in negative process 
## controls
EA.dat <- EA.dat %>%
        filter(!grepl("Glyptocephalus cynoglossus|Pagrus pagrus", species)) %>%
        droplevels()

## Calculate the total number of species and number of reads belonging
## to each group
EA.group.richness <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        group_by(group, species) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(group) %>%
        rename(richness = n) 

EA.group.reads <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        group_by(group) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() 

## Calculate the total number of species and number of reads belonging
## to each group for each estuary
EA.est.richness <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        group_by(estuary, group, species) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(estuary, group) %>%
        rename(richness = n) 

EA.est.reads <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        group_by(estuary, group) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup()

## Calculate the total number of species and number of reads belonging
## to each group for each season within each estuary
EA.season.richness <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        group_by(season_year, estuary, group, species) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(season_year, estuary, group) %>%
        rename(richness = n) 

EA.season.reads <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        group_by(season_year, estuary, group) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() 

## Calculate the total number of samples for each season in each estuary
EA.season.size <- EA.dat %>%
        filter(grepl("eDNA", sample_type)) %>%
        dplyr::select(season_year, estuary, location, bio_rep) %>%
        distinct() %>%
        group_by(season_year, estuary) %>%
        count(bio_rep) %>%
        ungroup()

## Merge EA.season.richness and EA.season.reads
EA.groups <- merge(EA.season.richness, EA.season.reads, 
                   by = c("season_year", "estuary", "group"), 
                   all.x = TRUE)

## Reorder factor levels of season_year and group
EA.groups <-  EA.groups %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf)) %>%
        mutate(group = fct_relevel(group,
                                   "Invertebrates","Amphibians","Birds",
                                   "Mammals","Elasmobranchs","Jawless fishes",
                                   "Ray-finned fishes"))

## Plot number of species detected for each group in each season by 
## estuary
p1a <- ggplot(EA.groups,
             aes(x = season_year, y = richness, 
                 group = group, fill = group)) + 
        geom_bar(stat = "identity", 
                 width = 0.5, colour = "black") +
        scale_fill_manual(name = "Group",
                          values = c("grey90","purple","gold",
                                     "violetred","darkcyan",
                                     "azure","navy")) +
        scale_y_continuous(limits = c(0, 120), breaks = seq(0,120,20)) +
        labs(title = "(a)", x = "Season", y = "Taxon richness") +
        theme_bw() +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(1.5, "lines")) +
        facet_grid(. ~ estuary, scales = "free", space = "free")
p1a

## Plot number of reads for each group in each season by estuary
p1b <- ggplot(EA.groups,
              aes(x = season_year, y = reads, 
                  group = group, fill = group)) + 
        geom_bar(stat = "identity", 
                 width = 0.5, colour = "black") +
        scale_fill_manual(name = "Group",
                          values = c("grey90","purple","gold",
                                     "violetred","darkcyan",
                                     "azure","navy")) +
        scale_y_continuous(limits = c(0, 700000), 
                           breaks = seq(0,700000,100000),
                           labels = comma) +
        labs(title = "(b)", x = "Season", y = "Read counts") +
        theme_bw() +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(1.5, "lines")) +
        facet_grid(. ~ estuary, scales = "free", space = "free")
p1b


## Summarise plots
g1 <- ggarrange(p1a, p1b, ncol = 1, nrow = 2,
                common.legend = TRUE, legend = "right",
                align = "hv")

ggsave(filename=here("results/figures/FigS1_raw_data_summaries.png"), 
       plot = g1, width = 13, height = 13, dpi = 300, units = "in")


## Now we will examine the molecular by-catch of using fish metabarcoding 
## primers on water from estuaries. Subset EA.dat for non-fish detections.
EA.nonfish <- EA.dat %>%
        filter(!grepl("Elasmobranchs|Invertebrates|Jawless fishes|Ray-finned fishes", 
                      group)) %>%
        droplevels()

## Rename vertebrate groups
EA.nonfish$group <- gsub("Amphibians", "Amphibian", EA.nonfish$group)
EA.nonfish$group <- gsub("Birds", "Bird", EA.nonfish$group)
EA.nonfish$group <- gsub("Mammals", "Mammal", EA.nonfish$group)



#' --- 
#' 
#' ## 2) Refine dataset
#' 
#' Now, we need to further refine the metabarcoding dataset by removing
#' spurious taxa.
#' 
#' 1. NBN Atlas was used to check species occurrence records and ensure 
#'    they matched sampling locations. This is also a good source for 
#'    checking current taxonomy.    
#' 2. Where species distribution did not overlap with sampling locations,
#'    the most common sequence for that taxa was used as input for an 
#'    NCBI BLAST search to confirm species identity.
#' 3. If the NCBI BLAST search returned multiple hits for different taxa
#'    and no single taxon possessed >98% identity to the query sequence,
#'    then we reassigned sequences as the species most likely to occur at
#'    the sampling locations.
#' 4. Where two or more species returned by NCBI BLAST were likely to 
#'    occur at the sampling locations, we reassigned sequences to genus
#'    or family level.
#'    
#' The original assignments and the proposed reassignments were stored
#' in a separate file that will be imported and matched against taxa
#' in the EA.dat dataframe.
#' 

## Import file containing original and new assignments
reassignments <- read.csv(here("assets/reassignments.csv"), header = TRUE)

## Correct misassignments by merging EA.dat and reasssignments dataframes
EA.correct <- merge(EA.nonfish, reassignments, 
                    by = "species", all.x = TRUE)

## Remove original assignments, rename column of corrected assignments,
## and move new columns to come before read counts. Convert all columns 
## except reads from character to factor.
EA.correct <- EA.correct %>% 
        dplyr::select(-species) %>%
        rename(sci_name = reassignment) %>%
        relocate(c(family, sci_name, common_name), .before = reads) %>% 
        mutate_if(is.character, as.factor)

## Now aggregate read counts for each species as some species detected 
## in each sample may have been duplicated during reassignment
estuary.dat <- EA.correct %>%
        group_by(primer_set, library, sample_type, season_year, 
                 estuary, location, bio_rep, pcr_rep, identifier,
                 group, family, sci_name, common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup()

## Calculate total reads for each PCR replicate, each biological replicate,
## and each sampling location. The latter two will be needed to calculate 
## the proportional read counts for taxa remaining in each biological 
## replicate and sampling location after removal of false positives and
## contaminants. Using the total reads for each PCR replicate, calculate
## the proportional read counts for taxa detected in each PCR replicate.
estuary.dat <- estuary.dat %>%
        group_by(primer_set, library, estuary, location, bio_rep, 
                 pcr_rep, identifier) %>%
        mutate(total_reads_by_pcr_rep = sum(reads)) %>%
        ungroup() %>%
        group_by(primer_set, library, estuary, location, bio_rep) %>%
        mutate(total_reads_by_bio_rep = sum(reads)) %>%
        ungroup() %>%
        group_by(estuary, location) %>%
        mutate(total_reads_by_location = sum(reads)) %>%
        ungroup() %>%
        mutate(PRC = reads/total_reads_by_pcr_rep)

        

#' --- 
#' 
#' ## 3) Contamination
#' 
#' Two species that occur in the UK are clear false positives as they
#' have highly restricted distributions and do not occur at the sampling 
#' locations.
#' 
#' - Edible dormouse (Glis glis) -> max. read frequency = 0.018292683
#' - Greater white-toothed shrew (Crocidura russula) -> max. read 
#'   frequency = 0.510040161
#' 
#' We will remove these species from the metabarcoding dataset, then 
#' examine how much contamination remains in the process controls.
#' 

## Remove Glis glis and Crocidura russula from dataset
spp.removed <- estuary.dat %>% 
        filter(!grepl("Glis glis|Crocidura russula", sci_name)) %>%
        dplyr::select(-c(total_reads_by_bio_rep, total_reads_by_location)) %>%
        droplevels()

## Create separate dataframes for eDNA samples and controls
samples <- spp.removed %>%
        filter(sample_type == "eDNA") %>%
        droplevels()

controls <- spp.removed %>% 
        filter(sample_type != "eDNA") %>%
        droplevels()

## Convert all columns except reads from character to factor, replace NA
## values with empty cells, then remove redundant columns from controls 
## dataframe.
controls <- controls %>% 
        mutate_if(is.factor, as.character) %>%
        replace(is.na(.), "") %>%
        mutate(season_year = gsub("NA", "", season_year)) %>%
        mutate(location = gsub("NA", "", location)) %>%
        unite("control_ID", c(sample_type, season_year, estuary, location), sep = " ", remove = FALSE) %>%
        mutate(control_ID = gsub(" ", "-", control_ID)) %>%
        mutate(control_ID = gsub("----", "", control_ID)) %>%
        dplyr::select(-c(season_year, estuary, location, bio_rep, 
                         pcr_rep, identifier)) %>%
        mutate(common_name = fct_relevel(common_name, "Human","Cow","Pig",
                                         "Sheep", "Gulls", "Field vole",
                                         "Water vole","European mole"))

## Add unique identifier to controls with the same ID so that each
## control will be plotted
controls$control_ID <- ave(controls$control_ID, controls$control_ID, 
                           FUN = function(i) paste0(i, '-', seq_along(i)))

## Add space in between blank types and reorder factor levels
controls <- controls %>%
        mutate(sample_type = gsub("Blank", " Blank", 
                                  controls$sample_type)) %>%
        mutate(sample_type = fct_relevel(sample_type,
                                         "Field Blank", 
                                         "Extraction Blank",
                                         "PCR Blank",
                                         "Well Blank"))

## Plot contamination found in process controls
hm1 <- ggplot(controls, aes(x = control_ID, 
                            y = fct_rev(as_factor(common_name)),
                            fill = PRC)) +
        geom_tile(colour = "black") +
        scale_fill_gradient(name = "Proportional\nread counts", 
                            limits = c(0, 1),
                            breaks = c(0, 0.25, 0.50, 0.75, 1),
                            low = "white", high = "red",
                            guide = guide_colourbar(frame.colour = "black",
                                                    ticks.colour = "black")) +
        labs(x = "Process controls", y = "Taxon") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              strip.text.y = element_text(angle = 360),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_blank(),
              axis.text.y = element_text(colour = "black"),
              axis.ticks.x = element_blank(),
              text = element_text(size = 20),
              legend.key.size = unit(1, "lines")) + 
        facet_grid(sample_type ~ library, scales = "free", space = "free")
hm1

ggsave(filename=here("results/figures/FigS2_Tele02_contamination.png"), 
       plot = hm1, width = 10, height = 8, dpi = 300, units = "in")



#' ---
#' 
#' ## 4) Remove potential false positives
#' 
#' Most contamination in the process controls is from cow (Bos taurus), 
#' sheep (Ovis aries), pig (Sus scrofa domesticus), and human (Homo 
#' sapiens). These species will be removed alongside other domestic 
#' species, including dog (Canis lupus familiaris), chicken (Gallus gallus
#' domesticus) and horse (Equus caballus ferus).
#' 
#' However, four wild species that could occur at the sampling locations, 
#' water vole (Arvicola amphibius), field vole (Microtus agrestis), mole
#' (Talpa europaea) and gulls (Larus spp.), were also detected in field
#' blanks, extraction blanks and PCR negative controls as contaminants.
#' 
#' Tissue samples from water vole, field vole and mole, along with several 
#' other mammal species, was sequenced for another study (Sales et al. 2020) 
#' before more stringent workspace separation and decontamination procedures 
#' were put in place. Therefore, these are all potential false positives
#' via contamination.
#' 
#' The proportional read counts for these contaminants are too high to use 
#' as sequence thresholds. Therefore, we will have to remove these species 
#' from eDNA samples instead. This will potentially result in false negatives 
#' for these species in the study estuaries. We will produce summary plots
#' of eDNA detections for each estuary with and without potential 
#' contaminants (i.e. species not detected in negative process controls 
#' but whose tissue was sequenced) for comparison.
#' 
#' The list of species constituting potential contaminants were stored 
#' in a separate file that will  be imported and matched against taxa in 
#' the samples dataframe.
#' 

#-----------------------------#
# INC. POTENTIAL CONTAMINANTS #
#-----------------------------#

## Remove human, cow, pig, sheep, horse, dog, chicken, water vole, field
## vole, European mole, and gulls from eDNA samples
all.rep <- samples %>% 
        filter(!grepl("Human|Cow|Sheep|Pig|Horse|Dog|Chicken|Water vole|Field vole|European mole|Gulls", 
                      common_name)) %>%
        dplyr::select(-identifier) %>%
        droplevels()

## Read counts for some species in some samples are exceptionally low 
## (<10 reads) and could be the result of sequencing error, tag jumping,
## or cross-contamination. We cannot apply a threshold based on actual
## read counts to remove these detections due to variation in total read
## counts across samples, i.e. some samples have <100 reads so 10 reads
## would constitute an abundant taxon in these samples. Instead, we will
## apply a threshold based on read count frequency, where all detections
## less than 1% of the total reads for each sample are removed.

## Apply false positive threshold
all.rep <- all.rep %>%
        group_by(primer_set, library, estuary, location, bio_rep, pcr_rep) %>%
        filter(PRC >= 0.01) %>%
        ungroup() %>%
        droplevels()

## Export as .csv file
write.csv(all.rep, 
          here("results/figures/EA_estuary_bio_pcr_rep_data_filtered_Salford-spp-inc_20210809.csv"), 
          row.names = FALSE)

## Calculate the total number of species belonging to each group in the
## refined dataset
RD.group.richness <- all.rep %>%
        group_by(group, sci_name, common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(group) %>%
        rename(richness = n)


## Now pool data for PCR replicates belonging to the same biological 
## replicate for each sample
pooled.pcr <- all.rep %>%
        group_by(primer_set, library, season_year, estuary, 
                 location, bio_rep, group, family, sci_name,
                 common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup()

## Add total and proportional read count columns for each biological
## replicate
pooled.pcr <- pooled.pcr %>%
        left_join(dplyr::select(estuary.dat, -c(sample_type, pcr_rep, 
                                                identifier, reads, 
                                                total_reads_by_pcr_rep, 
                                                total_reads_by_location, 
                                                PRC)), 
                  by = c("primer_set", "library", 
                         "season_year", "estuary", 
                         "location", "bio_rep", "group", 
                         "family", "sci_name", "common_name")) %>%
        distinct() %>%
        mutate(PRC = reads/total_reads_by_bio_rep) %>%
        droplevels()

## Export as .csv file
write.csv(pooled.pcr, 
          here("results/figures/EA_estuary_bio_rep_data_filtered_Salford-spp-inc_20210809.csv"), 
          row.names = FALSE)

## Calculate the number of samples from each season in each estuary 
## remaining in the refined dataset
RD.sample.size <- pooled.pcr %>%
        dplyr::select(season_year, estuary, location, bio_rep) %>%
        distinct() %>%
        group_by(season_year, estuary) %>%
        count(bio_rep) %>%
        ungroup()


## Now pool data for biological replicates belonging to the same sampling
## location
pooled.bio <- pooled.pcr %>%
        group_by(season_year, estuary, location, group, family, sci_name,
                 common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup()

## Remove total and proportional read count columns for each biological 
## replicate and add total and proportional read count columns for each 
## sampling location
pooled.bio <- pooled.bio %>%
        left_join(dplyr::select(estuary.dat, -c(primer_set, library, 
                                                sample_type, bio_rep, 
                                                pcr_rep, identifier, reads, 
                                                total_reads_by_pcr_rep, 
                                                total_reads_by_bio_rep, 
                                                PRC)), 
                  by = c("season_year", "estuary", "location", "group", 
                         "family", "sci_name", "common_name")) %>%
        distinct() %>%
        mutate(PRC = reads/total_reads_by_location) %>%
        droplevels()

## Export as .csv file
write.csv(pooled.bio, 
          here("results/figures/EA_estuary_sampling_location_data_filtered_Salford-spp-inc_20210809.csv"), 
          row.names = FALSE)


#-----------------------------#
# EXC. POTENTIAL CONTAMINANTS #
#-----------------------------#

## Import list of species whose tissue was sequenced at the University of
## Salford
Salford.spp <- read.csv(here("assets/salford-species-to-remove.csv"), header = TRUE)

## Remove human, cow, pig, sheep, horse, dog, chicken, water vole, field
## vole, European mole, and gulls from eDNA samples
all.rep.nc <- samples %>% 
        filter(!grepl("Human|Cow|Sheep|Pig|Horse|Dog|Chicken|Water vole|Field vole|European mole|Gulls", 
                      common_name)) %>%
        dplyr::select(-identifier) %>%
        droplevels()

## Remove taxa found in Salford.spp from all.rep
all.rep.nc <- all.rep.nc %>%
        anti_join(Salford.spp) %>%
        droplevels()

## Read counts for some species in some samples are exceptionally low 
## (<10 reads) and could be the result of sequencing error, tag jumping,
## or cross-contamination. We cannot apply a threshold based on actual
## read counts to remove these detections due to variation in total read
## counts across samples, i.e. some samples have <100 reads so 10 reads
## would constitute an abundant taxon in these samples. Instead, we will
## apply a threshold based on read count frequency, where all detections
## less than 1% of the total reads for each sample are removed.

## Apply false positive threshold
all.rep.nc <- all.rep.nc %>%
        group_by(primer_set, library, estuary, location, bio_rep, pcr_rep) %>%
        filter(PRC >= 0.01) %>%
        ungroup() %>%
        droplevels()

## Export as .csv file
write.csv(all.rep.nc, 
          here("results/figures/EA_estuary_bio_pcr_rep_data_filtered_Salford-spp-exc_20210809.csv"), 
          row.names = FALSE)

## Calculate the total number of species belonging to each group in the
## refined dataset
RD.nc.group.richness <- all.rep.nc %>%
        group_by(group, sci_name, common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(group) %>%
        rename(richness = n)


## Now pool data for PCR replicates belonging to the same biological 
## replicate for each sample
pooled.pcr.nc <- all.rep.nc %>%
        group_by(primer_set, library, season_year, estuary, 
                 location, bio_rep, group, family, sci_name,
                 common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup()

## Add total and proportional read count columns for each biological
## replicate
pooled.pcr.nc <- pooled.pcr.nc %>%
        left_join(dplyr::select(estuary.dat, -c(sample_type, pcr_rep, 
                                                identifier, reads, 
                                                total_reads_by_pcr_rep, 
                                                total_reads_by_location, 
                                                PRC)), 
                  by = c("primer_set", "library", 
                         "season_year", "estuary", 
                         "location", "bio_rep", "group", 
                         "family", "sci_name", "common_name")) %>%
        distinct() %>%
        mutate(PRC = reads/total_reads_by_bio_rep) %>%
        droplevels()

## Export as .csv file
write.csv(pooled.pcr.nc, 
          here("results/figures/EA_estuary_bio_rep_data_filtered_Salford-spp-exc_20210809.csv"), 
          row.names = FALSE)

## Calculate the number of samples from each season in each estuary 
## remaining in the refined dataset
RD.nc.sample.size <- pooled.pcr.nc %>%
        dplyr::select(season_year, estuary, location, bio_rep) %>%
        distinct() %>%
        group_by(season_year, estuary) %>%
        count(bio_rep) %>%
        ungroup()


## Now pool data for biological replicates belonging to the same sampling
## location
pooled.bio.nc <- pooled.pcr.nc %>%
        group_by(season_year, estuary, location, group, family, sci_name,
                 common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup()

## Remove total and proportional read count columns for each biological 
## replicate and add total and proportional read count columns for each 
## sampling location
pooled.bio.nc <- pooled.bio.nc %>%
        left_join(dplyr::select(estuary.dat, -c(primer_set, library, 
                                                sample_type, bio_rep, 
                                                pcr_rep, identifier, reads, 
                                                total_reads_by_pcr_rep, 
                                                total_reads_by_bio_rep, 
                                                PRC)), 
                  by = c("season_year", "estuary", "location", "group", 
                         "family", "sci_name", "common_name")) %>%
        distinct() %>%
        mutate(PRC = reads/total_reads_by_location) %>%
        droplevels()

## Export as .csv file
write.csv(pooled.bio.nc, 
          here("results/figures/EA_estuary_sampling_location_data_filtered_Salford-spp-exc_20210809.csv"), 
          row.names = FALSE)



#' ---
#'
#' ## 6) Basic summaries
#'
#' Summarise data in terms of proportional read counts and detection 
#' rate across samples with and without potential contaminants removed.
#' 

#-----------------------------#
# INC. POTENTIAL CONTAMINANTS #
#-----------------------------#

## Create dataframe containing richness of individual biological 
## replicates
bio.rep.richness <- pooled.pcr %>%
        count(season_year, estuary, location, bio_rep) %>%
        rename(richness = n) %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf))
        
## Plot richness of individual biological replicates
p2 <- ggplot(bio.rep.richness,
             aes(x = location, y = richness, 
                 group = bio_rep, fill = bio_rep)) + 
        geom_bar(stat = "identity", 
                 position = position_dodge(preserve = "single"), 
                 width = 0.5, colour = "black") +
        scale_fill_manual(name = "Biological \nreplicate",
                          values = c("black","grey60","grey90")) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(x = "Sampling location", y = "Taxon richness") +
        theme_bw() +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(1, "lines")) +
        facet_grid(estuary ~ season_year, scales = "free", space = "free")
p2


## Create dataframe containing detection rate across biological replicates 
## for each taxa
detection.rate <- pooled.pcr %>%
        group_by(season_year, estuary, location, sci_name, common_name) %>%
        count(bio_rep) %>%
        ungroup() %>%
        rename(no_bio_rep = n) %>%
        group_by(season_year, estuary, location) %>%
        mutate(total_bio_rep = n_distinct(bio_rep)) %>%
        ungroup() %>%
        mutate(SO = no_bio_rep/total_bio_rep) %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf))

## Plot number of biological replicates each taxa was detected in, i.e.
## sample occupancy
b1 <- ggplot(detection.rate, 
             aes(x = location, y = fct_rev(common_name), 
                 size = no_bio_rep)) + 
        geom_point() + 
        scale_size_continuous(name = "Number of \nbiological \nreplicates",
                              range = c(1, 5),
                              breaks = c(1, 2)) + 
        labs(x = "Sampling location", y = "Taxon") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(1, "lines")) + 
        facet_grid(estuary ~ season_year, scales = "free", space = "free")
b1


## Create vectors containing terrestrial bird and aquatic mammal species
terrestrial.bird <- c("Collared dove","Common blackbird",
                      "Common pheasant","Common starling",
                      "Corvids","Grey partridge",
                      "House sparrow","Red-legged partridge",
                      "Rock dove","Willow ptarmigan",
                      "Wood pigeon")
aquatic.mammal <- c("Eurasian beaver","Grey seal",
                    "Harbour porpoise","Harbour seal")

## Create vectors containing IUCN category for each species
rspb.red <- c("Common starling","Grey partridge","House sparrow",
              "Lapwing","Whimbrel","White-fronted goose")
rspb.amber <- c("Black-headed gull","Dunlin","Eurasian spoonbill",
                "Greylag goose","Guillemot","Mute swan",
                "Oystercatcher","Pintail","Redshank",
                "Shelduck","Whooper swan","Willow ptarmigan")
rspb.green <- c("Collared dove","Common blackbird","Common moorhen",
                "Corvids","Golden plover","Goosander",
                "Grey heron","Rock dove","Wood pigeon")
rspb.intro <- c("Canada goose","Common pheasant","Mandarin duck",
                "Muscovy duck","Red-legged partridge")
not.assessed <- c("Brown rat","European rabbit")
LC <- c("Bank vole","Daubenton's bat","European badger","Grey seal",
        "Harbour seal","Mustelids","Red deer","Smooth newt")
VU <- c("Harbour porpoise")
EN <- c("Eurasian beaver")

## Sort pooled.bio dataframe by group then species, then reset factor 
## levels of species according to new order of dataframe. Create new 
## column specifying whether bird and mammal species are aquatic or 
## terrestrial, and create new column containing red list category for
## each species according to RSPB (birds) and the Mammal Society 
## (mammals).
spp.summary <- pooled.bio %>%
        arrange(group, common_name) %>%   
        mutate(common_name = factor(common_name, unique(common_name))) %>%
        mutate(common_name = fct_rev(common_name)) %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf)) %>%
        mutate(subgroup = if_else(group == "Bird" & common_name %in% terrestrial.bird, "Terrestrial bird", 
                                  if_else(group == "Bird" & !common_name %in% terrestrial.bird, "Aquatic bird",
                                          if_else(group == "Mammal" & common_name %in% aquatic.mammal, "Aquatic mammal",
                                                  if_else(group == "Mammal" & !common_name %in% aquatic.mammal, "Terrestrial mammal",
                                                          "Amphibian"))))) %>%
        relocate(subgroup, .after = group) %>%
        mutate(subgroup = fct_relevel(subgroup, 
                                      c("Amphibian", 
                                        "Aquatic bird","Terrestrial bird",
                                        "Aquatic mammal","Terrestrial mammal"))) %>%
        mutate(category = if_else(common_name %in% rspb.red, "RSPB Red List", 
                                  if_else(common_name %in% rspb.amber, "RSPB Amber List",
                                          if_else(common_name %in% rspb.green, "RSPB Green List",
                                                  if_else(common_name %in% rspb.intro, "RSPB Introduced List",
                                                          if_else(common_name %in% not.assessed, "Not assessed",
                                                                  if_else(common_name %in% LC, "Least Concern",
                                                                          if_else(common_name %in% VU, "Vulnerable",
                                                                                  if_else(common_name %in% EN, "Endangered",
                                                                                          "Captive"))))))))) %>%
        relocate(category, .after = common_name)

## Calculate the total number of aquatic and terrestrial bird and mammal 
## species
spp.richness <- spp.summary %>%
        group_by(subgroup, sci_name, common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(subgroup) %>%
        rename(richness = n) 

## Plot species detected in different estuaries across different seasons
b2 <- ggplot(spp.summary, 
             aes(x = location, y = common_name, 
                 fill = subgroup, size = PRC)) + 
        geom_point(pch = 21) +
        scale_fill_manual(name = "Group",
                          values = c("#fff44f",
                                     "deepskyblue","lightskyblue1",
                                     "hotpink", "lightpink")) +
        scale_size_continuous(name = "Proportional \nread counts",
                              range = c(3, 10),
                              limits = c(0, 1),
                              breaks = c(0, 0.25, 0.5, 0.75, 1),
                              labels = c("> 0","0.25","0.50","0.75","1.00")) + 
        guides(fill = guide_legend(override.aes = list(size=10))) +
        labs(x = "Sampling location", y = "Taxon") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(2, "lines")) +
        facet_grid(estuary ~ season_year, scales = "free", space = "free")
b2

ggsave(filename=here("results/figures/FigS4_detection_site_season_Salford-spp-inc.svg"), 
       plot = b2, width = 20, height = 25, dpi = 300, units = "in")



#-----------------------------#
# EXC. POTENTIAL CONTAMINANTS #
#-----------------------------#

## Create dataframe containing richness of individual biological 
## replicates
bio.rep.richness.nc <- pooled.pcr.nc %>%
        count(season_year, estuary, location, bio_rep) %>%
        rename(richness = n) %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf))

## Plot richness of individual biological replicates
p3 <- ggplot(bio.rep.richness.nc,
             aes(x = location, y = richness, 
                 group = bio_rep, fill = bio_rep)) + 
        geom_bar(stat = "identity", 
                 position = position_dodge(preserve = "single"), 
                 width = 0.5, colour = "black") +
        scale_fill_manual(name = "Biological \nreplicate",
                          values = c("black","grey60","grey90")) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(x = "Sampling location", y = "Taxon richness") +
        theme_bw() +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(1, "lines")) +
        facet_grid(estuary ~ season_year, scales = "free", space = "free")
p3

ggsave(filename=here("results/figures/FigS3_bio_rep_richness.png"), 
       plot = p3, width = 15, height = 10, dpi = 300, units = "in")


## Create dataframe containing detection rate across biological replicates 
## for each taxa
detection.rate.nc <- pooled.pcr.nc %>%
        group_by(season_year, estuary, location, sci_name, common_name) %>%
        count(bio_rep) %>%
        ungroup() %>%
        rename(no_bio_rep = n) %>%
        group_by(season_year, estuary, location) %>%
        mutate(total_bio_rep = n_distinct(bio_rep)) %>%
        ungroup() %>%
        mutate(SO = no_bio_rep/total_bio_rep) %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf))

## Plot number of biological replicates each taxa was detected in, i.e.
## sample occupancy
b3 <- ggplot(detection.rate.nc, 
             aes(x = location, y = fct_rev(common_name), 
                 size = no_bio_rep)) + 
        geom_point() + 
        scale_size_continuous(name = "Number of \nbiological \nreplicates",
                              range = c(1, 5),
                              breaks = c(1, 2)) + 
        labs(x = "Sampling location", y = "Taxon") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(1, "lines")) + 
        facet_grid(estuary ~ season_year, scales = "free", space = "free")
b3


## Create vectors containing terrestrial bird and aquatic mammal species
terrestrial.bird <- c("Collared dove","Common blackbird",
                      "Common pheasant","Common starling",
                      "Corvids","Grey partridge",
                      "House sparrow","Red-legged partridge",
                      "Rock dove","Willow ptarmigan",
                      "Wood pigeon")
aquatic.mammal <- c("Grey seal","Harbour porpoise","Harbour seal")

## Create vectors containing IUCN category for each species
rspb.red <- c("Common starling","Grey partridge","House sparrow",
              "Lapwing","Whimbrel","White-fronted goose")
rspb.amber <- c("Black-headed gull","Dunlin","Eurasian spoonbill",
                "Greylag goose","Guillemot","Mute swan",
                "Oystercatcher","Pintail","Redshank",
                "Shelduck","Whooper swan","Willow ptarmigan")
rspb.green <- c("Collared dove","Common blackbird","Common moorhen",
                "Corvids","Golden plover","Goosander",
                "Grey heron","Rock dove","Wood pigeon")
rspb.intro <- c("Canada goose","Common pheasant","Mandarin duck",
                "Muscovy duck","Red-legged partridge")
not.assessed <- c("European rabbit")
LC <- c("Daubenton's bat","Grey seal","Harbour seal")
VU <- c("Harbour porpoise")

## Sort pooled.bio dataframe by group then species, then reset factor 
## levels of species according to new order of dataframe. Create new 
## column specifying whether bird and mammal species are aquatic or 
## terrestrial, and create new column containing red list category for
## each species according to RSPB (birds) and the Mammal Society 
## (mammals).
spp.summary.nc <- pooled.bio.nc %>%
        arrange(group, common_name) %>%   
        mutate(common_name = factor(common_name, unique(common_name))) %>%
        mutate(common_name = fct_rev(common_name)) %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf)) %>%
        mutate(subgroup = if_else(group == "Bird" & common_name %in% terrestrial.bird, "Terrestrial bird", 
                                  if_else(group == "Bird" & !common_name %in% terrestrial.bird, "Aquatic bird",
                                          if_else(group == "Mammal" & common_name %in% aquatic.mammal, "Aquatic mammal",
                                                  if_else(group == "Mammal" & !common_name %in% aquatic.mammal, "Terrestrial mammal",
                                                          "Other"))))) %>%
        relocate(subgroup, .after = group) %>%
        mutate(subgroup = fct_relevel(subgroup, 
                                      c("Aquatic bird","Terrestrial bird",
                                        "Aquatic mammal","Terrestrial mammal"))) %>%
        mutate(category = if_else(common_name %in% rspb.red, "RSPB Red List", 
                                  if_else(common_name %in% rspb.amber, "RSPB Amber List",
                                          if_else(common_name %in% rspb.green, "RSPB Green List",
                                                  if_else(common_name %in% rspb.intro, "RSPB Introduced List",
                                                          if_else(common_name %in% not.assessed, "Not assessed",
                                                                  if_else(common_name %in% LC, "Least Concern",
                                                                          if_else(common_name %in% VU, "Vulnerable",
                                                                                  "Captive")))))))) %>%
        relocate(category, .after = common_name)

## Calculate the total number of aquatic and terrestrial bird and mammal 
## species
spp.richness.nc <- spp.summary.nc %>%
        group_by(subgroup, sci_name, common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        count(subgroup) %>%
        rename(richness = n) 

## Plot species detected in different estuaries across different seasons
b4 <- ggplot(spp.summary.nc, 
             aes(x = location, y = common_name, 
                 fill = subgroup, size = PRC)) + 
        geom_point(pch = 21) +
        scale_fill_manual(name = "Group",
                          values = c("deepskyblue","lightskyblue1",
                                     "hotpink", "lightpink")) +
        scale_size_continuous(name = "Proportional \nread counts",
                              range = c(3, 10),
                              limits = c(0, 1),
                              breaks = c(0, 0.25, 0.5, 0.75, 1),
                              labels = c("> 0","0.25","0.50","0.75","1.00")) + 
        guides(fill = guide_legend(override.aes = list(size=10))) +
        labs(x = "Sampling location", y = "Taxon") + 
        theme_bw() + 
        theme(panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black", angle = 60, hjust = 1),
              axis.text.y = element_text(colour = "black"),
              strip.text.y = element_text(angle = 360),
              text = element_text(size = 20),
              legend.key.size = unit(2, "lines")) +
        facet_grid(estuary ~ season_year, scales = "free", space = "free")
b4

ggsave(filename=here("results/figures/Fig1_detection_site_season.svg"), 
       plot = b4, width = 20, height = 22, dpi = 300, units = "in")


## Now that false positives and contaminants have been removed, and the
## data have been pooled by sampling location, we can perform diversity
## analyses.



#' ---
#'
#' ## 7) Alpha and beta diversity 
#'
#' Compare species richness and community dissimilarity between estuaries
#' across seasons using 'vegan' and 'betapart'. First, we will examine
#' diversity in each estuary, then we will create datasets for each 
#' estuary to examine seasonal changes in diversity.
#' 

######################
# ESTUARY COMPARISON #
######################

#-----------------#
# ALPHA DIVERSITY #
#-----------------#

## Check if data are normally distributed
## Histogram:
hist(bio.rep.richness.nc$richness)  # asymmetrical

## Quantile-quantile plot (Q-Q plot):
qqnorm(bio.rep.richness.nc$richness)
qqline(bio.rep.richness.nc$richness, lty=2)  # observations don't tail away much

## Quantitative tests for normal distribution
## Shapiro-Wilk test for normality:
shapiro.test(bio.rep.richness.nc$richness)  # P = 2.432e-05

## Kolmogorov-Smirnov test:
ks.test(bio.rep.richness.nc$richness, pnorm)  # P < 2.2e-16

## Check whether any common data transformations improve data distribution
## Natural log:
hist(log(bio.rep.richness.nc$richness))
shapiro.test(log(bio.rep.richness.nc$richness))  # P = 0.000277

## log10:
hist(log10(bio.rep.richness.nc$richness))
shapiro.test(log10(bio.rep.richness.nc$richness))  # P = 0.000277

## Exponential:
hist(exp(bio.rep.richness.nc$richness))
shapiro.test(exp(bio.rep.richness.nc$richness))  # P < 2.2e-16

## logit:
hist(log(bio.rep.richness.nc$richness/(bio.rep.richness.nc$richness)))
#shapiro.test(log(bio.rep.richness.nc$richness/(bio.rep.richness.nc$richness)))  # NA

## Square root:
hist(sqrt(bio.rep.richness.nc$richness))
shapiro.test(sqrt(bio.rep.richness.nc$richness))  # P = 0.0009333

## Reciprocal:
hist(1/bio.rep.richness.nc$richness)
shapiro.test(1/(bio.rep.richness.nc$richness))  # P = 2.295e-08

## Taxon richness data is not normally distributed and transformations 
## do not resolve this issue. Compare variance in taxon richness between 
## different estuaries and seasons.
levene.test(bio.rep.richness.nc$richness,
            bio.rep.richness.nc$estuary,
            location = "mean")  # P = 0.03576

## It is unlikely that the data will conform to all the assumptions of a
## one-way ANOVA, but we will run the model and assess the residuals
anova <- lm(richness ~ estuary, data = bio.rep.richness.nc)
summary(anova)  # summaries for each factor level
summary.aov(anova)  # summary for factor as a whole
model.tables(aov(anova), "means", se = TRUE)  # mean values for each factor level
TukeyHSD(aov(anova))  # Tukey post-hoc test
plot(TukeyHSD(aov(anova)))  # visualise pairwise comparisons

# Get standardised residuals for model validation
sresid <- resid(anova, type="pearson")

## Check assumption of normal distribution
## Most models are robust to slight deviations from normality in the 
## residuals
hist(sresid)
lines(density(sresid,adjust=1))
qqnorm(sresid, cex=1.8, pch=20)
qqline(sresid, lty=2, lwd=2)
shapiro.test(sresid) # P = 0.0202
## Not normally distributed

## Check assumption of homogeneity of variances
plot(sresid ~ anova$fitted.values)
plot(sresid ~ bio.rep.richness.nc$estuary)
## No heteroscedasticity present

## Check assumption of no collinearity
## Only one variable being modelled so no collinearity present.

## Check dataset does not contain serial auto-correlation
## This can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(anova)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation

## Check model not biased by unduly influential observations
influence <- influence.measures(anova)
summary(influence)
CD <- cooks.distance(anova)
plot(CD ~ sresid)
## Cook's distance <1 for so observations not exerting strong influence 
## on model parameters

## Several model assumptions are violated and data cannot be transformed to
## conform to a normal distribution, thus a one-way ANOVA is not applicable. 
## In these scenarios, switching to a Generalised Linear Model (GLM) and 
## changing the combination of error family and link-function terms to
## achieve normally distributed residuals is recommended. A Poisson 
regression <- glm(richness ~ estuary, 
                  family = poisson(link="log"), 
                  data = bio.rep.richness.nc)
summary(regression)
anova(regression, test = "Chi")
drop1(regression, test = "Chi")
TukeyHSD(aov(regression))

## Check model meets GLM assumptions
## Test for overdispersion
83.365/66
1-pchisq(83.365, df = 66)  # overdispersed

## Try negative binomial distribution to counter overdispersion
regression <- glm.nb(richness ~ estuary, data = bio.rep.richness.nc)
summary(regression)
anova(regression, test = "Chi")
drop1(regression, test = "Chi")
TukeyHSD(aov(regression))

## Check model meets GLM assumptions
## Test for overdispersion
66.276/66
1-pchisq(66.276, df = 66)  # not overdispersed

## Plot the fitted data against the observed data
plot(bio.rep.richness.nc$richness ~ fitted(regression))

## Perform model validation checks to ensure model is good fit to data 
## and making reliable predictions
## Assumption 1: residuals are normally distributed
sresid <- resid(regression, type = "pearson")
hist(sresid)
lines(density(sresid,adjust = 1))
qqnorm(sresid, cex = 1.8, pch = 20)
qqline(sresid, lty = 2, lwd = 2)
chkres(regression)
shapiro.test(sresid) # P = 0.005519
## Some deviation from normality as residuals are not normally distributed
## therefore model is unreliable.

## Assumption 2: no heteroscedascity
## Plot standardised residuals against each independent variable to 
## identify any source of heterogeneity i.e. independent variable that is 
## non-linearly associated with y
plot(sresid ~ regression$fitted.values, pch = 20, cex = 2, cex.lab = 1.5)
plot(sresid ~ bio.rep.richness.nc$estuary, pch = 20, cex = 2, cex.lab = 1.5) 
## No heteroscedascity

## Assumption 3: no collinearity
## Only one variable being modelled so no collinearity present.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
durbinWatsonTest(regression)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation

## Assumption 5: model not biased by unduly influential observations
influence <- influence.measures(regression)
summary(influence)
CD <- cooks.distance(regression)
plot(CD ~ sresid)
## Cook's distance < 1 so observations not exerting strong influence on 
## model parameters

## GLM did not improve distribution of the residuals. We will use the 
## non-parametric Kruskal-Wallis test instead (NB: this tests for 
## differences between the median values of different groups, not the 
## mean values).
kruskal.test(richness ~ estuary, data = bio.rep.richness.nc)

## Dunn's test of multiple comparisons
dunnTest(richness ~ estuary, 
         data = bio.rep.richness.nc,
         method = "bh")

## KW test indicates significant difference in species richness between
## estuaries. Plot species richness in each estuary.
p4 <- ggplot(bio.rep.richness.nc,
             aes(x = estuary, y = richness)) + 
        geom_jitter(aes(colour = estuary), cex = 5, width = 0.3, alpha = 0.7) + 
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        annotate("text", 
                 x = c("ESK", "TEES", "TWEED"), y = 15, 
                 label = c("a","b","a"), cex=10) + 
        scale_colour_manual(name = "Estuary",
                            values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(title = expression(bold("(a)"~alpha~-"diversity")),
             x = "Estuary", y = "Taxon richness") +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_line(colour = "white"),
              panel.grid.minor = element_line(colour = "white"), 
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              text = element_text(size = 20),
              legend.key = element_blank(),
              legend.position = "none")
p4


#----------------------------------#
# RAREFACTION/EXTRAPOLATION CURVES #
#----------------------------------#

## Look at whether differences in richness are influenced by sampling
## effort by performing rarefaction and extrapolation curves using the
## iNEXT package.

## Make a copy of the pooled.pcr dataframe where reads are summed by
## season_year, estuary, location, biological replicate, and species
estuary.beta <- pooled.pcr.nc %>%
        group_by(season_year, estuary, location, bio_rep, sci_name,
                 common_name) %>%
        summarise(reads = sum(reads), .groups = "keep") %>%
        ungroup() %>%
        unite("sample_ID", c(season_year, estuary, location, bio_rep), 
              sep = "-", remove = FALSE) %>%
        mutate(sample_ID = gsub(" ", "-", sample_ID))

## Create metadata for beta diversity analyses
estuary.metadata <- estuary.beta %>% 
        dplyr::select(-c(sci_name, common_name, reads)) %>%
        distinct() %>%
        mutate(season_year = fct_relevel(season_year, 
                                         "Autumn 2017", 
                                         after = Inf))

## Now convert estuary.beta into wide format for beta diversity analyses
estuary.beta <- estuary.beta %>%
        dplyr::select(-c(season_year, estuary, location, bio_rep, 
                         common_name)) %>%
        pivot_wider(names_from = sci_name, 
                    values_from = reads,
                    values_fill = 0) %>%
        column_to_rownames("sample_ID")

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
estuary.beta <- estuary.beta[!sapply(estuary.beta, function(x) all(x == 0))]
estuary.beta <- estuary.beta[!apply(estuary.beta == 0, 1, all),]

## Convert read count data to presence-absence
estuary.beta[estuary.beta > 0] <- 1

## Create vectors of sample IDs associated with each estuary.
ESK <- droplevels(subset(estuary.metadata, select = "sample_ID", estuary == "ESK"))
TEES <- droplevels(subset(estuary.metadata, select = "sample_ID", estuary == "TEES"))
TWEED <- droplevels(subset(estuary.metadata, select = "sample_ID", estuary == "TWEED"))

## Create copies of estuary.beta for each estuary
ESK.beta <- estuary.beta[rownames(estuary.beta) %in% ESK$sample_ID,]
TEES.beta <- estuary.beta[rownames(estuary.beta) %in% TEES$sample_ID,]
TWEED.beta <- estuary.beta[rownames(estuary.beta) %in% TWEED$sample_ID,]

## Calculate incidence frequency for each species detected in samples 
## from each estuary
ESK.freq <- colSums(ESK.beta)
TEES.freq <- colSums(TEES.beta)
TWEED.freq <- colSums(TWEED.beta)

## Add sample size to beginning of each richness vector. This is because 
## the first entry of each list for iNEXT must be the total number of 
## sampling units, followed by the species incidence frequencies.
ESK.freq <- append(ESK.freq, 33, after = 0)
TEES.freq <- append(TEES.freq, 16, after = 0)
TWEED.freq <- append(TWEED.freq, 20, after = 0)

## Make list of incidence frequency data for each estuary
estuary.richness <- list(ESK.freq, TEES.freq, TWEED.freq)
names(estuary.richness) <- c("ESK", "TEES", "TWEED")

## Run iNEXT function with 50 samples for each estuary:
est.re.50 <- iNEXT(estuary.richness, q = 0, datatype = "incidence_freq", 
                    endpoint = 50, knots = 10, se = TRUE, conf = 0.95,
                    nboot = 1000)
est.re.50

## Sample-size-based R/E curves
p5a <- ggiNEXT(est.re.50, type = 1) + 
        scale_colour_manual(values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_fill_manual(values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_shape_manual(values = c(19,19,19)) +
        scale_x_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
        scale_y_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
        labs(title = "\n(b) Sample size-based R/E curve",
             y = "Taxon diversity") +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, colour = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              legend.position = "bottom",
              legend.box = "vertical",
              legend.key = element_blank(),
              text = element_text(size = 20))
p5a

## Sample completeness curves
p5b <- ggiNEXT(est.re.50, type=2) +
        scale_colour_manual(values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_fill_manual(values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_shape_manual(values = c(19,19,19)) +
        scale_x_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
        scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
        labs(title = "Sample completeness curve") +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, colour = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.key = element_blank(),
              text = element_text(size = 20))
p5b

## Coverage-based R/E curves
p5c <- ggiNEXT(est.re.50, type=3) + 
        scale_colour_manual(values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_fill_manual(values = c("grey30","goldenrod2","dodgerblue3")) +
        scale_shape_manual(values = c(19,19,19)) +
        scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
        scale_y_continuous(limits = c(0,50), breaks = seq(0,50,10)) +
        labs(title = "Coverage-based R/E curve",
             y = "Taxon diversity") +
        theme(panel.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, colour = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              legend.position = "bottom",
              legend.box = "horizontal",
              legend.key = element_blank(),
              text = element_text(size = 20))
p5c

## Apply the estimateD() function to obtain diversity estimates of order 
## q = 0, 1, 2 for any particular level of sample size (base="size") or
## any specified level of sample coverage (base="coverage")
estimateD(estuary.richness, datatype="incidence_freq", base="size",
          level=NULL, conf=0.95)
estimateD(estuary.richness, datatype="incidence_freq", base="coverage", 
          level=0.95, conf=0.95)


#----------------#
# BETA DIVERSITY #
#----------------#

## Remove samples that are extreme outliers: Autumn-2017-TEES-SEI06-B
estuary.beta <- estuary.beta[which(!grepl("Autumn-2017-TEES-SEI06-B", 
                                          rownames(estuary.beta))),]

## Also remove these outlier samples from metadata
estuary.metadata <- estuary.metadata[which(!grepl("Autumn-2017-TEES-SEI06-B",
                                                  estuary.metadata$sample_ID)),]
rownames(estuary.metadata) <- NULL

## Beta diversity across ESK samples
ESK.multi <- beta.multi(ESK.beta, index.family="jaccard")
print(ESK.multi)

## Beta diversity across TEES samples
TEES.multi <- beta.multi(TEES.beta, index.family="jaccard")
print(TEES.multi)

## Beta diversity across TWEED samples
TWEED.multi <- beta.multi(TWEED.beta, index.family="jaccard")
print(TWEED.multi)

## The majority of total beta diversity arises from taxon turnover
## rather than nestedness for samples from each estuaries.

## Pairwise between-site values of each component of beta diversity
estuary.dist <- beta.pair(estuary.beta, index.family="jaccard")


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
estuary.bd.turn <- betadisper(estuary.dist$beta.jtu, estuary.metadata$estuary)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
mod.turn <- with(estuary.metadata, estuary.bd.turn)
mod.turn

## Compute mean distance to centroid per group
tapply(estuary.bd.turn$distances, estuary.metadata$estuary, mean)

## Compute variance per group
tapply(estuary.bd.turn$distances, estuary.metadata$estuary, var)

## Ordination plot of distances to centroid
plot(estuary.bd.turn)

## Boxplot of distances to centroid
boxplot(estuary.bd.turn, xlab="estuary", xaxt="n", bty="n")
axis(side=1, at=c(1:3), labels=c("ESK","TEES","TWEED"))

## Plots indicate that there is some difference in multivariate 
## dispersion of turnover partition between estuaries. Statistically 
## check whether variance is different between estuaries using standard 
## parametric anova or permutation tests.
anova(estuary.bd.turn)     # Significant difference between estuaries
permutest(estuary.bd.turn) # Significant difference between estuaries

## Analyse pairwise differences between groups (estuaryators) using 
## parametric Tukey's HSD test.
TukeyHSD(estuary.bd.turn)  # Significant difference between groups

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
estuary.comm.turn <- metaMDS(estuary.dist$beta.jtu, 
                             dist = "jaccard", 
                             k = 2,
                             maxit = 999,
                             trymax = 1000,
                             noshare = TRUE,
                             wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
estuary.comm.turn$stress
stressplot(estuary.comm.turn)

## Plot site scores as text
ordiplot(estuary.comm.turn, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
estuary.NMDS1 <- estuary.comm.turn$points[,1]
estuary.NMDS2 <- estuary.comm.turn$points[,2]
estuary.turn.NMDS <- data.frame(NMDS1=estuary.NMDS1, 
                             NMDS2=estuary.NMDS2,
                             estuary = estuary.metadata$estuary)

## Check data
head(estuary.turn.NMDS)

## Plot data frame
p6a <- ggplot(estuary.turn.NMDS, 
              aes(x = NMDS1, y = NMDS2, colour = estuary)) + 
        geom_point(cex = 5, alpha = 0.3) + 
        stat_ellipse() + 
        scale_colour_manual(name="Estuary",
                            values=c("grey30","goldenrod2","dodgerblue3")) +
        scale_x_continuous(limits=c(-0.6,0.8), breaks=seq(-0.6,0.8,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.6,0.6), breaks=seq(-0.6,0.6,0.3),
                                labels=scales::number_format(accuracy = 0.1)) + 
        labs(title = expression(bold("(c)"~beta~-"diversity")), 
             subtitle = "(i) Turnover", 
             x = "NMDS1", y = "NMDS2") + 
        annotate("text", x = 0.7, y = 0.6, 
                 label = "stress == 0.131", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.key = element_blank(),
              legend.position = "none")
p6a


## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
estuary.turn.adonis <- adonis(estuary.dist$beta.jtu ~ estuary, 
                              estuary.metadata)

## Inspect results:
estuary.turn.adonis

## Result is significant. There is a substantial difference in taxon 
## replacement (i.e. turnover) between estuaries. Therefore, taxa
## in one estuary are substituted by another species in a different
## estuary.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
estuary.bd.nest <- betadisper(estuary.dist$beta.jne, estuary.metadata$estuary)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
estuary.nest <- with(estuary.metadata, estuary.bd.nest)
estuary.nest

## Compute mean distance to centroid per group
tapply(estuary.bd.nest$distances, estuary.metadata$estuary, mean)

## Compute variance per group
tapply(estuary.bd.nest$distances, estuary.metadata$estuary, var)

## Ordination plot of distances to centroid
plot(estuary.bd.nest)

## Boxplot of distances to centroid
boxplot(estuary.bd.nest, xlab="estuary", xaxt="n", bty="n")
axis(side=1, at=c(1:3), labels=c("ESK","TEES","TWEED"))

## Plots indicate that there is no difference in multivariate dispersion
## between estuaries. Statistically check whether variance is different 
## between estuaries using standard parametric anova or permutation tests.
anova(estuary.bd.nest)     # No significant difference between estuaries
permutest(estuary.bd.nest) # No significant difference between estuaries

## Analyse pairwise differences between groups using parametric Tukey's HSD 
## test.
TukeyHSD(estuary.bd.nest)  # No significant difference between estuaries

## Ordination of beta diversity partitioned by nestedness:
estuary.comm.nest <- metaMDS(estuary.dist$beta.jne, 
                             dist = "jaccard",
                             k = 2,
                             maxit = 999,
                             trymax = 1000,
                             noshare = TRUE,
                             wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
estuary.comm.nest$stress
stressplot(estuary.comm.nest)

## plot site scores as text
ordiplot(estuary.comm.nest, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
estuary.NMDS1 <- estuary.comm.nest$points[,1]
estuary.NMDS2 <- estuary.comm.nest$points[,2]
estuary.nest.NMDS <- data.frame(NMDS1=estuary.NMDS1, 
                                NMDS2=estuary.NMDS2,
                                estuary = estuary.metadata$estuary)

## Check data
head(estuary.nest.NMDS)

## Plot data frame
p6b <- ggplot(estuary.nest.NMDS, 
              aes(x = NMDS1, y = NMDS2, colour = estuary)) + 
        geom_point(cex = 5, alpha = 0.3) + 
        stat_ellipse() + 
        scale_colour_manual(name="Estuary",
                            values=c("grey30","goldenrod2","dodgerblue3")) + 
        scale_x_continuous(limits=c(-0.8,0.6), breaks=seq(-0.8,0.6,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.6,0.6), breaks=seq(-0.6,0.6,0.3),
                           labels=scales::number_format(accuracy = 0.1)) + 
        labs(title = "", 
             subtitle = "(ii) Nestedness-resultant", 
             x = "NMDS1", y = "NMDS2") + 
        annotate("text", x = 0.5, y = 0.5, 
                 label="stress == 0.203", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.key = element_blank(),
              legend.position = "none")
p6b

## Statistically check difference in nestedness of communities
## Look at PermANOVA using adonis(), considered more robust than anosim()
estuary.nest.adonis <- adonis(estuary.dist$beta.jne ~ estuary, 
                              estuary.metadata)

## Inspect results
## no summary() or plot() diets included
estuary.nest.adonis

## Result is not significant. There is no substantial difference in taxon
## loss or gain (i.e. nestedness) between estuaries.


## 3. TOTAL BETA DIVERSITY
estuary.bd.total <- betadisper(estuary.dist$beta.jac, estuary.metadata$estuary)

## Check homogeneity of multivariate dispersions. Groups being tested should 
## have the sample multivariate spread to conform to assumptions of PERMANOVA
mod.total <- with(estuary.metadata, estuary.bd.total)
mod.total

## Compute mean distance to centroid per group
tapply(estuary.bd.total$distances, estuary.metadata$estuary, mean)

## Compute variance per group
tapply(estuary.bd.total$distances, estuary.metadata$estuary, var)

## Ordination plot of distance to centroids
plot(estuary.bd.total)

## Boxplot of distance to centroids
boxplot(estuary.bd.total, xlab="estuary", xaxt="n", bty="n")
axis(side=1, at=c(1:3), labels=c("ESK","TEES","TWEED"))

## Plots indicates that there is some difference in multivariate dispersion 
## between estuaries. Statistically check whether variance is different 
## between estuaries using standard parametric anova or permutation tests.
anova(estuary.bd.total)     # Significant difference between estuaries
permutest(estuary.bd.total) # Significant difference between estuaries

## Analyse pairwise differences between groups using parametric Tukey's 
## HSD test.
TukeyHSD(estuary.bd.total)  # Significant difference between estuaries

## Ordination of total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
estuary.comm.total <- metaMDS(estuary.dist$beta.jac, 
                              dist = "jaccard",
                              k = 2,
                              maxit = 999,
                              trymax = 1000,
                              noshare = TRUE,
                              wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
estuary.comm.total$stress
stressplot(estuary.comm.total)

## plot site scores as text
ordiplot(estuary.comm.total, display = "sites", type = "text")

## Build data frame with NMDS coordinates and metadata
estuary.NMDS1 <- estuary.comm.total$points[,1]
estuary.NMDS2 <- estuary.comm.total$points[,2]
estuary.total.NMDS <- data.frame(NMDS1 = estuary.NMDS1,
                                 NMDS2 = estuary.NMDS2,
                                 estuary = estuary.metadata$estuary)

## Check data
head(estuary.total.NMDS)

## Plot data frame
p6c <- ggplot(estuary.total.NMDS, 
              aes(x = NMDS1, y = NMDS2, colour = estuary)) + 
        geom_point(cex = 5, alpha = 0.3) + 
        stat_ellipse() + 
        scale_x_continuous(limits=c(-0.4,0.8), breaks=seq(-0.4,0.8,0.2),
                           labels=scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits=c(-0.5,0.5), breaks=seq(-0.5,0.5,0.25),
                           labels=scales::number_format(accuracy = 0.01)) + 
        labs(title = "", 
             subtitle = expression(bold("(iii) Total"~beta~-"diversity")), 
             x = "NMDS1", y = "NMDS2") + 
        annotate("text", x = 0.7, y = 0.5, 
                 label = "stress == 0.149", 
                 cex = 5, parse = TRUE) + 
        scale_colour_manual(name="Estuary",
                            values=c("grey30","goldenrod2","dodgerblue3")) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.key = element_blank(),
              legend.position = "none")
p6c

## Statistically check difference in total beta diversity
## Look at PermANOVA using adonis(), considered more robust than anosim()
estuary.total.adonis <- adonis(estuary.dist$beta.jac ~ estuary,
                               estuary.metadata)

## Inspect results
## no summary() or plot() diets included
estuary.total.adonis

## Again result is significant. There is substantial variation in overall 
## community composition between estuaries.

## Plot all diversity results
g2 <- ggarrange(ggarrange(p4, p5a, nrow = 2, ncol = 1),
                ggarrange(p6a, p6b, p6c, nrow = 3, ncol = 1, 
                          align = "hv", common.legend = TRUE,
                          legend = "bottom"), 
                nrow = 1, ncol = 2)

ggsave(filename=here("results/figures/Fig2_estuary_diversity.png"), 
       plot = g2, width = 15, height = 15, dpi = 300, units = "in")



#######################
# SEASONAL COMPARISON #
#######################

#-----------------#
# ALPHA DIVERSITY #
#-----------------#

## Create dataframe of taxon richness for each estuary
ESK.richness <- filter(bio.rep.richness.nc, estuary == "ESK")
TEES.richness <- filter(bio.rep.richness.nc, estuary == "TEES")
TWEED.richness <- filter(bio.rep.richness.nc, estuary == "TWEED")

## Use the non-parametric Kruskal-Wallis test and Dunn's test of multiple 
## comparisons to examine seasonal variation in species richness for each 
## estuary
kruskal.test(richness ~ season_year, data = ESK.richness)
dunnTest(richness ~ season_year,
         data = ESK.richness,
         method = "bh")

kruskal.test(richness ~ season_year, data = TEES.richness)
dunnTest(richness ~ season_year,
         data = TEES.richness,
         method = "bh")

kruskal.test(richness ~ season_year, data = TWEED.richness)
dunnTest(richness ~ season_year,
         data = TWEED.richness,
         method = "none")

## KW test indicates that seasonal variation only occurs in the ESK
## estuary.

## Create text labels denoting significance
p7.text <- tribble(~season_year, ~estuary, ~richness, ~label,
                   "Autumn 2016", "ESK", 15, "a",
                   "Spring 2017", "ESK", 15, "b",
                   "Autumn 2017", "ESK", 15, "a",
                   "Autumn 2016", "TEES", 15, "a",
                   "Spring 2017", "TEES", 15, "a",
                   "Autumn 2017", "TEES", 15, "a",
                   "Spring 2017", "TWEED", 15, "a",
                   "Autumn 2017", "TWEED", 15, "a")

## Plot seasonal variation in species richness in each estuary
p7 <- ggplot(bio.rep.richness.nc,
             aes(x = season_year, y = richness)) + 
        geom_jitter(aes(colour = season_year, fill = season_year, 
                        shape = season_year), 
                    cex = 5, width = 0.2, alpha = 0.7) + 
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(x = "Sampling event", y = "Taxon richness") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black", angle=60, hjust=1),
              axis.text.y = element_text(colour="black"),
              strip.text.y = element_text(angle=360),
              text = element_text(size=20),
              legend.position = "none") +
        facet_grid(. ~ estuary, scales="free") +
        geom_text(data = p7.text, 
                  aes(x = season_year, y = richness,  label = label),
                  cex = 10)
p7


#-----------------------------#
# BETA DIVERSITY: ESK ESTUARY #
#-----------------------------#

## Subset metadata for ESK estuary
ESK.metadata <- filter(estuary.metadata, estuary == "ESK")

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
ESK.beta <- ESK.beta[!sapply(ESK.beta, function(x) all(x == 0))]
ESK.beta <- ESK.beta[!apply(ESK.beta == 0, 1, all),]

## Create vectors of sample IDs associated with each season
Aut16 <- droplevels(subset(ESK.metadata, select = "sample_ID", season_year == "Autumn 2016"))
Spr17 <- droplevels(subset(ESK.metadata, select = "sample_ID", season_year == "Spring 2017"))
Aut17 <- droplevels(subset(ESK.metadata, select = "sample_ID", season_year == "Autumn 2017"))

## Subset ESK.beta for samples from each season
Aut16.beta <- ESK.beta[rownames(ESK.beta) %in% Aut16$sample_ID,]
Spr17.beta <- ESK.beta[rownames(ESK.beta) %in% Spr17$sample_ID,]
Aut17.beta <- ESK.beta[rownames(ESK.beta) %in% Aut17$sample_ID,]

## Beta diversity across Autumn 2016 samples
Aut16.multi <- beta.multi(Aut16.beta, index.family="jaccard")
print(Aut16.multi)

## Beta diversity across Spring 2017 samples
Spr17.multi <- beta.multi(Spr17.beta, index.family="jaccard")
print(Spr17.multi)

## Beta diversity across Autumn 2017 samples
Aut17.multi <- beta.multi(Aut17.beta, index.family="jaccard")
print(Aut17.multi)

## The majority of total beta diversity arises from taxon turnover
## rather than nestedness for ESK samples from all seasons.

## Pairwise between-site values of each component of beta diversity
ESK.dist <- beta.pair(ESK.beta, index.family="jaccard")


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
ESK.bd.turn <- betadisper(ESK.dist$beta.jtu, ESK.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
ESK.turn <- with(ESK.metadata, ESK.bd.turn)
ESK.turn

## Compute mean distance to centroid per group
tapply(ESK.bd.turn$distances, ESK.metadata$season_year, mean)

## Compute variance per group
tapply(ESK.bd.turn$distances, ESK.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(ESK.bd.turn)

## Boxplot of distances to centroid
boxplot(ESK.bd.turn, xlab = "Sampling event", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:3), labels = c("Autumn 2016","Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions
## between seasons. Statistically check whether turnover is different 
## between seasons using standard parametric anova or permutation tests.
anova(ESK.bd.turn)     # No significant difference between seasons
permutest(ESK.bd.turn) # No significant difference between seasons

## Analyse pairwise differences between groups (seasons) using 
## parametric Tukey's HSD test.
TukeyHSD(ESK.bd.turn)  # No significant difference between seasons

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
ESK.comm.turn <- metaMDS(ESK.dist$beta.jtu, 
                         distance = "jaccard", 
                         k = 2,
                         maxit = 999,
                         trymax = 1000,
                         noshare = TRUE,
                         wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
ESK.comm.turn$stress
stressplot(ESK.comm.turn)

## Plot site scores as text
ordiplot(ESK.comm.turn, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
ESK.NMDS1 <- ESK.comm.turn$points[,1]
ESK.NMDS2 <- ESK.comm.turn$points[,2]
ESK.turn.NMDS <- data.frame(NMDS1 = ESK.NMDS1, 
                            NMDS2 = ESK.NMDS2,
                            Season = ESK.metadata$season_year)

## Check data
head(ESK.turn.NMDS)

## Plot data frame
p8a <- ggplot(ESK.turn.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.65, 1), breaks = seq(-0.6, 1, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "(a) ESK",
             subtitle = "(i) Turnover", 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.85, y = 0.6, 
                 label = "stress == 0.093", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "bottom",
              legend.key = element_blank())
p8a

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
ESK.turn.adonis <- adonis(ESK.dist$beta.jtu ~ season_year, ESK.metadata)

## Inspect results:
ESK.turn.adonis

## Result is significant. There is a substantial difference in species 
## replacement (i.e. turnover) between seasons within the ESK estuary. 
## Therefore, species in one season are substituted by species in a 
## different season.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
ESK.bd.nest <- betadisper(ESK.dist$beta.jne, ESK.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
ESK.nest <- with(ESK.metadata, ESK.bd.nest)
ESK.nest

## Compute mean distance to centroid per group
tapply(ESK.bd.nest$distances, ESK.metadata$season_year, mean)

## Compute variance per group
tapply(ESK.bd.nest$distances, ESK.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(ESK.bd.nest)

## Boxplot of distances to centroid
boxplot(ESK.bd.nest, xlab = "Season", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:3), labels = c("Autumn 2016","Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions
## between seasons. Statistically check whether nestedness is different 
## between seasons using standard parametric anova or permutation tests.
anova(ESK.bd.nest)     # No significant difference between seasons
permutest(ESK.bd.nest) # No significant difference between seasons

## Analyse pairwise differences between groups (sites) using 
## parametric Tukey's HSD test.
TukeyHSD(ESK.bd.nest)  # No significant difference between seasons

## Ordination of beta diversity partitioned by nestedness:
## The metaMDS function automatically transforms data and checks solution
## robustness
ESK.comm.nest <- metaMDS(ESK.dist$beta.jne, 
                         distance = "jaccard", 
                         k = 2,
                         maxit = 999,
                         trymax = 1000,
                         noshare = TRUE,
                         wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
ESK.comm.nest$stress
stressplot(ESK.comm.nest)

## Plot site scores as text
ordiplot(ESK.comm.nest, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
ESK.NMDS1 <- ESK.comm.nest$points[,1]
ESK.NMDS2 <- ESK.comm.nest$points[,2]
ESK.nest.NMDS <- data.frame(NMDS1 = ESK.NMDS1, 
                            NMDS2 = ESK.NMDS2,
                            Season = ESK.metadata$season_year)

## Check data
head(ESK.nest.NMDS)

## Plot data frame
p8b <- ggplot(ESK.nest.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = "(ii) Nestedness-resultant", 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.5, y = 0.6, 
                 label = "stress == 0.187", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p8b

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
ESK.nest.adonis <- adonis(ESK.dist$beta.jne ~ season_year, ESK.metadata)

## Inspect results:
ESK.nest.adonis

## Result is significant. There is a difference in species loss or 
## gain (i.e. nestedness) between sites.


## 3. TOTAL BETA DIVERSITY
ESK.bd.total <- betadisper(ESK.dist$beta.jac, ESK.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
ESK.total <- with(ESK.metadata, ESK.bd.total)
ESK.total

## Compute mean distance to centroid per group
tapply(ESK.bd.total$distances, ESK.metadata$season_year, mean)

## Compute variance per group
tapply(ESK.bd.total$distances, ESK.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(ESK.bd.total)

## Boxplot of distances to centroid
boxplot(ESK.bd.total, xlab = "Season", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:3), labels = c("Autumn 2016","Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions 
## between seasons. Statistically check whether variance is different 
## between seasons using standard parametric anova or permutation tests.
anova(ESK.bd.total)     # No significant difference between seasons
permutest(ESK.bd.total) # No significant difference between seasons

## Analyse pairwise differences between groups (seasons) using 
## parametric Tukey's HSD test.
TukeyHSD(ESK.bd.total)  # No significant difference between seasons

## Ordination of beta diversity partitioned by total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
ESK.comm.total <- metaMDS(ESK.dist$beta.jac, 
                          distance = "jaccard", 
                          k = 2,
                          maxit = 999,
                          trymax = 1000,
                          noshare = TRUE,
                          wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
ESK.comm.total$stress
stressplot(ESK.comm.total)

## Plot site scores as text
ordiplot(ESK.comm.total, display = "sites", type = "text", cex=0.5)

## Build data frame with NMDS coordinates and metadata
ESK.NMDS1 <- ESK.comm.total$points[,1]
ESK.NMDS2 <- ESK.comm.total$points[,2]
ESK.total.NMDS <- data.frame(NMDS1 = ESK.NMDS1, 
                             NMDS2 = ESK.NMDS2,
                             Season = ESK.metadata$season_year)

## Check data
head(ESK.total.NMDS)

## Plot data frame
p8c <- ggplot(ESK.total.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = expression(bold("(iii) Total"~beta~-"diversity")), 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.7, y = 0.6, 
                 label = "stress == 0.119", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p8c

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
ESK.total.adonis <- adonis(ESK.dist$beta.jac ~ season_year, ESK.metadata)

## Inspect results:
ESK.total.adonis

## Again result is significant. There is substantial variation in overall 
## community composition of samples from different seasons within estuaries.


#------------------------------#
# BETA DIVERSITY: TEES ESTUARY #
#------------------------------#

## Subset metadata for ESK estuary
TEES.metadata <- filter(estuary.metadata, estuary == "TEES")

## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
TEES.beta <- TEES.beta[!sapply(TEES.beta, function(x) all(x == 0))]
TEES.beta <- TEES.beta[!apply(TEES.beta == 0, 1, all),]

## Remove samples that are extreme outliers: Autumn-2017-TEES-SEI06-B
TEES.beta <- TEES.beta[which(!grepl("Autumn-2017-TEES-SEI06-B", 
                                    rownames(TEES.beta))),]

## Also remove these outlier samples from metadata
TEES.metadata <- TEES.metadata[which(!grepl("Autumn-2017-TEES-SEI06-B",
                                            TEES.metadata$sample_ID)),]
rownames(TEES.metadata) <- NULL

## Create vectors of sample IDs associated with each season
Aut16 <- droplevels(subset(TEES.metadata, select = "sample_ID", season_year == "Autumn 2016"))
Spr17 <- droplevels(subset(TEES.metadata, select = "sample_ID", season_year == "Spring 2017"))
Aut17 <- droplevels(subset(TEES.metadata, select = "sample_ID", season_year == "Autumn 2017"))

## Subset TEES.beta for samples from each season
Aut16.beta <- TEES.beta[rownames(TEES.beta) %in% Aut16$sample_ID,]
Spr17.beta <- TEES.beta[rownames(TEES.beta) %in% Spr17$sample_ID,]
Aut17.beta <- TEES.beta[rownames(TEES.beta) %in% Aut17$sample_ID,]

## Beta diversity across August 2016 samples
Aut16.multi <- beta.multi(Aut16.beta, index.family="jaccard")
print(Aut16.multi)

## Beta diversity across Spring 2017 samples
Spr17.multi <- beta.multi(Spr17.beta, index.family="jaccard")
print(Spr17.multi)

## Beta diversity across August 2017 samples
Aut17.multi <- beta.multi(Aut17.beta, index.family="jaccard")
print(Aut17.multi)

## The majority of total beta diversity arises from taxon turnover
## rather than nestedness for TEES samples from all seasons.

## Pairwise between-site values of each component of beta diversity
TEES.dist <- beta.pair(TEES.beta, index.family="jaccard")


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
TEES.bd.turn <- betadisper(TEES.dist$beta.jtu, TEES.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
TEES.turn <- with(TEES.metadata, TEES.bd.turn)
TEES.turn

## Compute mean distance to centroid per group
tapply(TEES.bd.turn$distances, TEES.metadata$season_year, mean)

## Compute variance per group
tapply(TEES.bd.turn$distances, TEES.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(TEES.bd.turn)

## Boxplot of distances to centroid
boxplot(TEES.bd.turn, xlab = "Sampling event", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:3), labels = c("Autumn 2016","Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions
## between seasons. Statistically check whether turnover is different 
## between seasons using standard parametric anova or permutation tests.
anova(TEES.bd.turn)     # Significant difference between seasons
permutest(TEES.bd.turn) # Significant difference between seasons

## Analyse pairwise differences between groups (seasons) using 
## parametric Tukey's HSD test.
TukeyHSD(TEES.bd.turn)  # Significant difference between seasons

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
TEES.comm.turn <- metaMDS(TEES.dist$beta.jtu, 
                          distance = "jaccard", 
                          k = 2,
                          maxit = 999,
                          trymax = 1000,
                          noshare = TRUE,
                          wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
TEES.comm.turn$stress
stressplot(TEES.comm.turn)

## Plot site scores as text
ordiplot(TEES.comm.turn, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
TEES.NMDS1 <- TEES.comm.turn$points[,1]
TEES.NMDS2 <- TEES.comm.turn$points[,2]
TEES.turn.NMDS <- data.frame(NMDS1 = TEES.NMDS1, 
                             NMDS2 = TEES.NMDS2,
                             Season = TEES.metadata$season_year)

## Check data
head(TEES.turn.NMDS)

## Plot data frame
p9a <- ggplot(TEES.turn.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.6, 1), breaks = seq(-0.6, 1, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.8, 0.6), breaks = seq(-0.8, 0.6, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "(b) TEES",
             subtitle = "(i) Turnover", 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.9, y = 0.6, 
                 label = "stress < 0.001", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p9a

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
TEES.turn.adonis <- adonis(TEES.dist$beta.jtu ~ season_year, TEES.metadata)

## Inspect results:
TEES.turn.adonis

## Result is significant. There is a substantial difference in species 
## replacement (i.e. turnover) between seasons within the TEES estuary. 
## Therefore, species in one season are substituted by species in a 
## different season.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
TEES.bd.nest <- betadisper(TEES.dist$beta.jne, TEES.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
TEES.nest <- with(TEES.metadata, TEES.bd.nest)
TEES.nest

## Compute mean distance to centroid per group
tapply(TEES.bd.nest$distances, TEES.metadata$season_year, mean)

## Compute variance per group
tapply(TEES.bd.nest$distances, TEES.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(TEES.bd.nest)

## Boxplot of distances to centroid
boxplot(TEES.bd.nest, xlab = "Season", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:3), labels = c("Autumn 2016","Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions
## between seasons. Statistically check whether nestedness is different 
## between seasons using standard parametric anova or permutation tests.
anova(TEES.bd.nest)     # No significant difference between seasons
permutest(TEES.bd.nest) # No significant difference between seasons

## Analyse pairwise differences between groups (sites) using 
## parametric Tukey's HSD test.
TukeyHSD(TEES.bd.nest)  # No significant difference between seasons

## Ordination of beta diversity partitioned by nestedness:
## The metaMDS function automatically transforms data and checks solution
## robustness
TEES.comm.nest <- metaMDS(TEES.dist$beta.jne, 
                          distance = "jaccard", 
                          k = 2,
                          maxit = 999,
                          trymax = 1000,
                          noshare = TRUE,
                          wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
TEES.comm.nest$stress
stressplot(TEES.comm.nest)

## Plot site scores as text
ordiplot(TEES.comm.nest, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
TEES.NMDS1 <- TEES.comm.nest$points[,1]
TEES.NMDS2 <- TEES.comm.nest$points[,2]
TEES.nest.NMDS <- data.frame(NMDS1 = TEES.NMDS1, 
                             NMDS2 = TEES.NMDS2,
                             Season = TEES.metadata$season_year)

## Check data
head(TEES.nest.NMDS)

## Plot data frame
p9b <- ggplot(TEES.nest.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25),
                           labels = scales::number_format(accuracy = 0.01)) + 
        scale_y_continuous(limits = c(-1.05, 1), breaks = seq(-1, 1, 0.5),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = "(ii) Nestedness-resultant", 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.85, y = 1, 
                 label = "stress == 0.063", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p9b

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
TEES.nest.adonis <- adonis(TEES.dist$beta.jne ~ season_year, TEES.metadata)

## Inspect results:
TEES.nest.adonis

## Result is not significant. There is no difference in species loss or 
## gain (i.e. nestedness) between seasons.


## 3. TOTAL BETA DIVERSITY
TEES.bd.total <- betadisper(TEES.dist$beta.jac, TEES.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
TEES.total <- with(TEES.metadata, TEES.bd.total)
TEES.total

## Compute mean distance to centroid per group
tapply(TEES.bd.total$distances, TEES.metadata$season_year, mean)

## Compute variance per group
tapply(TEES.bd.total$distances, TEES.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(TEES.bd.total)

## Boxplot of distances to centroid
boxplot(TEES.bd.total, xlab = "Season", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:3), labels = c("Autumn 2016","Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions 
## between seasons. Statistically check whether variance is different 
## between seasons using standard parametric anova or permutation tests.
anova(TEES.bd.total)     # No significant difference between seasons
permutest(TEES.bd.total) # No significant difference between seasons

## Analyse pairwise differences between groups (seasons) using 
## parametric Tukey's HSD test.
TukeyHSD(TEES.bd.total)  # No significant difference between seasons

## Ordination of beta diversity partitioned by total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
TEES.comm.total <- metaMDS(TEES.dist$beta.jac, 
                           distance = "jaccard", 
                           k = 2,
                           maxit = 999,
                           trymax = 1000,
                           noshare = TRUE,
                           wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
TEES.comm.total$stress
stressplot(TEES.comm.total)

## Plot site scores as text
ordiplot(TEES.comm.total, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
TEES.NMDS1 <- TEES.comm.total$points[,1]
TEES.NMDS2 <- TEES.comm.total$points[,2]
TEES.total.NMDS <- data.frame(NMDS1 = TEES.NMDS1, 
                              NMDS2 = TEES.NMDS2,
                              Season = TEES.metadata$season_year)

## Check data
head(TEES.total.NMDS)

## Plot data frame
p9c <- ggplot(TEES.total.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.8, 1.2), breaks = seq(-0.8, 1.2, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = expression(bold("(iii) Total"~beta~-"diversity")), 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 1, y = 1, 
                 label = "stress < 0.001", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p9c

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
TEES.total.adonis <- adonis(TEES.dist$beta.jac ~ season_year, TEES.metadata)

## Inspect results:
TEES.total.adonis

## Again result is significant. There is substantial variation in overall 
## community composition of samples from different seasons within estuaries.


#-------------------------------#
# BETA DIVERSITY: TWEED ESTUARY #
#-------------------------------#

## Subset metadata for ESK estuary
TWEED.metadata <- estuary.metadata %>%
        filter(estuary == "TWEED") %>%
        droplevels()
        
## Remove empty samples and taxonomic assignments as these are problematic 
## for vegan
TWEED.beta <- TWEED.beta[!sapply(TWEED.beta, function(x) all(x == 0))]
TWEED.beta <- TWEED.beta[!apply(TWEED.beta == 0, 1, all),]

## Create vectors of sample IDs associated with each season
Spr17 <- droplevels(subset(TWEED.metadata, select = "sample_ID", season_year == "Spring 2017"))
Aut17 <- droplevels(subset(TWEED.metadata, select = "sample_ID", season_year == "Autumn 2017"))

## Subset TWEED.beta for samples from each season
Spr17.beta <- TWEED.beta[rownames(TWEED.beta) %in% Spr17$sample_ID,]
Aut17.beta <- TWEED.beta[rownames(TWEED.beta) %in% Aut17$sample_ID,]

## Beta diversity across TWEED Spring 2017 samples
Spr17.multi <- beta.multi(Spr17.beta, index.family="jaccard")
print(Spr17.multi)

## Beta diversity across TWEED Autumn 2017 samples
Aut17.multi <- beta.multi(Aut17.beta, index.family="jaccard")
print(Aut17.multi)

## The majority of total beta diversity arises from taxon turnover
## rather than nestedness for TWEED samples from all seasons.

## Pairwise between-site values of each component of beta diversity
TWEED.dist <- beta.pair(TWEED.beta, index.family="jaccard")


## 1. TURNOVER PARTITION (COMPLETE CHANGE IN COMMUNITIES)
TWEED.bd.turn <- betadisper(TWEED.dist$beta.jtu, TWEED.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
TWEED.turn <- with(TWEED.metadata, TWEED.bd.turn)
TWEED.turn

## Compute mean distance to centroid per group
tapply(TWEED.bd.turn$distances, TWEED.metadata$season_year, mean)

## Compute variance per group
tapply(TWEED.bd.turn$distances, TWEED.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(TWEED.bd.turn)

## Boxplot of distances to centroid
boxplot(TWEED.bd.turn, xlab = "Sampling event", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:2), labels = c("Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions
## between seasons. Statistically check whether turnover is different 
## between seasons using standard parametric anova or permutation tests.
anova(TWEED.bd.turn)     # No significant difference between seasons
permutest(TWEED.bd.turn) # No significant difference between seasons

## Analyse pairwise differences between groups (seasons) using 
## parametric Tukey's HSD test.
TukeyHSD(TWEED.bd.turn)  # No significant difference between seasons

## Ordination of beta diversity partitioned by turnover:
## The metaMDS function automatically transforms data and checks solution
## robustness
TWEED.comm.turn <- metaMDS(TWEED.dist$beta.jtu, 
                           distance = "jaccard", 
                           k = 2,
                           maxit = 999,
                           trymax = 1000,
                           noshare = TRUE,
                           wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
TWEED.comm.turn$stress
stressplot(TWEED.comm.turn)

## Plot site scores as text
ordiplot(TWEED.comm.turn, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
TWEED.NMDS1 <- TWEED.comm.turn$points[,1]
TWEED.NMDS2 <- TWEED.comm.turn$points[,2]
TWEED.turn.NMDS <- data.frame(NMDS1 = TWEED.NMDS1, 
                              NMDS2 = TWEED.NMDS2,
                              Season = TWEED.metadata$season_year)

## Check data
head(TWEED.turn.NMDS)

## Plot data frame
p10a <- ggplot(TWEED.turn.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(24,22)) +
        scale_colour_manual(values = c("purple","cyan4")) +
        scale_fill_manual(values = c("purple","cyan4")) +
        scale_x_continuous(limits = c(-0.6, 0.8), breaks = seq(-0.6, 0.8, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "(c) TWEED",
             subtitle = "(i) Turnover", 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.65, y = 0.6, 
                 label = "stress == 0.149", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p10a

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
TWEED.turn.adonis <- adonis(TWEED.dist$beta.jtu ~ season_year, TWEED.metadata)

## Inspect results:
TWEED.turn.adonis

## Result is significant. There is a substantial difference in species 
## replacement (i.e. turnover) between seasons within the TWEED estuary. 
## Therefore, species in one season are substituted by species in a 
## different season.


## 2. NESTEDNESS PARTITION (SUBSET OF CHANGE IN WIDER COMMUNITY)
TWEED.bd.nest <- betadisper(TWEED.dist$beta.jne, TWEED.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
TWEED.nest <- with(TWEED.metadata, TWEED.bd.nest)
TWEED.nest

## Compute mean distance to centroid per group
tapply(TWEED.bd.nest$distances, TWEED.metadata$season_year, mean)

## Compute variance per group
tapply(TWEED.bd.nest$distances, TWEED.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(TWEED.bd.nest)

## Boxplot of distances to centroid
boxplot(TWEED.bd.nest, xlab = "Season", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:2), labels = c("Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions
## between seasons. Statistically check whether nestedness is different 
## between seasons using standard parametric anova or permutation tests.
anova(TWEED.bd.nest)     # Significant difference between seasons
permutest(TWEED.bd.nest) # Significant difference between seasons

## Analyse pairwise differences between groups (sites) using 
## parametric Tukey's HSD test.
TukeyHSD(TWEED.bd.nest)  # Significant difference between seasons

## Ordination of beta diversity partitioned by nestedness:
## The metaMDS function automatically transforms data and checks solution
## robustness
TWEED.comm.nest <- metaMDS(TWEED.dist$beta.jne, 
                           distance = "jaccard", 
                           k = 2,
                           maxit = 999,
                           trymax = 1000,
                           noshare = TRUE,
                           wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
TWEED.comm.nest$stress
stressplot(TWEED.comm.nest)

## Plot site scores as text
ordiplot(TWEED.comm.nest, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
TWEED.NMDS1 <- TWEED.comm.nest$points[,1]
TWEED.NMDS2 <- TWEED.comm.nest$points[,2]
TWEED.nest.NMDS <- data.frame(NMDS1 = TWEED.NMDS1, 
                              NMDS2 = TWEED.NMDS2,
                              Season = TWEED.metadata$season_year)

## Check data
head(TWEED.nest.NMDS)

## Plot data frame
p10b <- ggplot(TWEED.nest.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(24,22)) +
        scale_colour_manual(values = c("purple","cyan4")) +
        scale_fill_manual(values = c("purple","cyan4")) +
        scale_x_continuous(limits = c(-0.8, 1), breaks = seq(-0.8, 1, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.7, 0.5), breaks = seq(-0.7, 0.5, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = "(ii) Nestedness-resultant", 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.85, y = 0.5, 
                 label = "stress = 0.190", 
                 cex = 5) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p10b

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
TWEED.nest.adonis <- adonis(TWEED.dist$beta.jne ~ season_year, TWEED.metadata)

## Inspect results:
TWEED.nest.adonis

## Result is not significant. There is no difference in species loss or 
## gain (i.e. nestedness) between seasons.


## 3. TOTAL BETA DIVERSITY
TWEED.bd.total <- betadisper(TWEED.dist$beta.jac, TWEED.metadata$season_year)

## Check homogeneity of multivariate dispersions. Groups being tested 
## should have the sample multivariate spread to conform to assumptions 
## of PERMANOVA.
TWEED.total <- with(TWEED.metadata, TWEED.bd.total)
TWEED.total

## Compute mean distance to centroid per group
tapply(TWEED.bd.total$distances, TWEED.metadata$season_year, mean)

## Compute variance per group
tapply(TWEED.bd.total$distances, TWEED.metadata$season_year, var)

## Ordination plot of distances to centroid
plot(TWEED.bd.total)

## Boxplot of distances to centroid
boxplot(TWEED.bd.total, xlab = "Season", xaxt = "n", bty = "n")
axis(side = 1, at = c(1:2), labels = c("Spring 2017","Autumn 2017"))

## Plots indicate that there is some difference in multivariate dispersions 
## between seasons. Statistically check whether variance is different 
## between seasons using standard parametric anova or permutation tests.
anova(TWEED.bd.total)     # No significant difference between seasons
permutest(TWEED.bd.total) # No significant difference between seasons

## Analyse pairwise differences between groups (seasons) using 
## parametric Tukey's HSD test.
TukeyHSD(TWEED.bd.total)  # No significant difference between seasons

## Ordination of beta diversity partitioned by total beta diversity:
## The metaMDS function automatically transforms data and checks solution
## robustness
TWEED.comm.total <- metaMDS(TWEED.dist$beta.jac, 
                            distance = "jaccard", 
                            k = 2,
                            maxit = 999,
                            trymax = 1000,
                            noshare = TRUE,
                            wascores = TRUE)

## Assess goodness of ordination fit (stress plot)
TWEED.comm.total$stress
stressplot(TWEED.comm.total)

## Plot site scores as text
ordiplot(TWEED.comm.total, display = "sites", type = "text", cex = 0.5)

## Build data frame with NMDS coordinates and metadata
TWEED.NMDS1 <- TWEED.comm.total$points[,1]
TWEED.NMDS2 <- TWEED.comm.total$points[,2]
TWEED.total.NMDS <- data.frame(NMDS1 = TWEED.NMDS1, 
                               NMDS2 = TWEED.NMDS2,
                               Season = TWEED.metadata$season_year)

## Check data
head(TWEED.total.NMDS)

## Plot data frame
p10c <- ggplot(TWEED.total.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(24,22)) +
        scale_colour_manual(values = c("purple","cyan4")) +
        scale_fill_manual(values = c("purple","cyan4")) +
        scale_x_continuous(limits = c(-0.6, 0.8), breaks = seq(-0.6, 0.8, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = expression(bold("(iii) Total"~beta~-"diversity")), 
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.7, y = 0.6, 
                 label = "stress == 0.153", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p10c

## Statistically check difference in spatial turnover of communities
## Look at PERMANOVA using adonis(), considered more robust than anosim()
TWEED.total.adonis <- adonis(TWEED.dist$beta.jac ~ season_year, TWEED.metadata)

## Inspect results:
TWEED.total.adonis

## Again result is significant. There is substantial variation in overall 
## community composition of samples from different seasons within estuaries.


#---------#
# SUMMARY #
#---------#

## Plot all beta diversity components for seasonal community composition
g4 <- ggarrange(p8a, p8b, p8c,
                p9a, p9b, p9c,
                p10a, p10b, p10c,
                nrow = 3, ncol = 3,
                align = "hv",
                common.legend = TRUE, 
                legend = "bottom")

ggsave(filename=here("results/figures/FigS5_seasonal_beta_diversity_all_components.png"), 
       plot = g4, width = 20, height = 15, dpi = 300, units = "in")


## Re-plot seasonal richness and community composition in each estuary 
## and produce summary plot showing alpha and beta diversity
p11a <- ggplot(bio.rep.richness.nc,
              aes(x = season_year, y = richness)) + 
        geom_jitter(aes(colour = season_year, fill = season_year, 
                        shape = season_year), 
                    cex = 5, width = 0.2, alpha = 0.7) + 
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        scale_shape_manual(name = "Season",
                           values = c(22,24,22)) +
        scale_colour_manual(name = "Season",
                            values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(name = "Season",
                          values = c("darkorange","purple","cyan4")) +
        scale_y_continuous(limits = c(0, 15)) +
        labs(title = expression(bold("(a)"~alpha~-"diversity")),
             x = "Sampling event", y = "Taxon richness") +
        theme_bw() +
        theme(panel.grid.major = element_line(colour="white"),
              panel.grid.minor = element_line(colour="white"), 
              axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour="black", angle=60, hjust=1),
              axis.text.y = element_text(colour="black"),
              strip.text.y = element_text(angle=360),
              text = element_text(size=20),
              legend.position = "none",
              legend.key = element_blank()) +
        facet_grid(. ~ estuary, scales="free_x") +
        geom_text(data = p7.text, 
                  aes(x = season_year, y = richness,  label = label),
                  cex = 10)
p11a

p11b <- ggplot(ESK.total.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.4, 0.6), breaks = seq(-0.4, 0.6, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = expression(bold("(b) Total"~beta~-"diversity")),
             subtitle = "(i) ESK",
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.7, y = 0.6, 
                 label = "stress == 0.119", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p11b

p11c <- ggplot(TEES.total.NMDS, 
              aes(x = NMDS1, y = NMDS2, 
                  colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(22,24,22)) +
        scale_colour_manual(values = c("darkorange","purple","cyan4")) +
        scale_fill_manual(values = c("darkorange","purple","cyan4")) +
        scale_x_continuous(limits = c(-0.8, 1.2), breaks = seq(-0.8, 1.2, 0.4),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-1.1, 1), breaks = seq(-1, 1, 0.5),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = "(ii) TEES",
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 1, y = 1, 
                 label = "stress < 0.001", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p11c

p11d <- ggplot(TWEED.total.NMDS, 
               aes(x = NMDS1, y = NMDS2, 
                   colour = Season, fill = Season, shape = Season)) + 
        geom_point(cex = 5, alpha = 0.5) + 
        stat_ellipse() +
        scale_shape_manual(values = c(24,22)) +
        scale_colour_manual(values = c("purple","cyan4")) +
        scale_fill_manual(values = c("purple","cyan4")) +
        scale_x_continuous(limits = c(-0.6, 0.8), breaks = seq(-0.6, 0.8, 0.2),
                           labels = scales::number_format(accuracy = 0.1)) + 
        scale_y_continuous(limits = c(-0.6, 0.6), breaks = seq(-0.6, 0.6, 0.3),
                           labels = scales::number_format(accuracy = 0.1)) + 
        labs(title = "",
             subtitle = "(iii) TWEED",
             x = "NMDS1",y = "NMDS2") + 
        annotate("text", x = 0.65, y = 0.6, 
                 label = "stress == 0.153", 
                 cex = 5, parse = TRUE) + 
        theme(panel.background = element_rect(fill = "white"),
              axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid"),
              axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
              axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
              axis.text.x = element_text(colour = "black"),
              axis.text.y = element_text(colour = "black"),
              plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
              plot.subtitle = element_text(face = "bold", hjust = 0, color = "black", margin = unit(c(2, 0, 0, 0), "mm")),
              text = element_text(size = 20),
              legend.position = "none")
p11d


## Summarise plots
g5 <- ggarrange(p11a,
                ggarrange(p11b, p11c, p11d, 
                          nrow = 1, ncol = 3, 
                          align = "hv"), 
                nrow = 2, ncol = 1,
                common.legend = TRUE, 
                legend = "bottom")

ggsave(filename=here("results/figures/Fig3_seasonal_alpha_beta_diversity.png"), 
       plot = g5, width = 18, height = 15, dpi = 300, units = "in")

writeLines("\n...Analysis completed!\n...")
#################
# END OF SCRIPT #
#################
