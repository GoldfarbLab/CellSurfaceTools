library(tidyverse)
library(here)
library(homologene)
library(data.table)
library(devtools)

#############################################
## Create Gene-level table
#############################################

# read human surface/plasma
load(here('data/human_surface_and_plasma_membrane_geneLevel.rda'))
human_surface_plasma_geneLevel <- as_tibble(human_surface_and_plasma_membrane_geneLevel)
human_surface_plasma_geneLevel <- mutate(human_surface_plasma_geneLevel, GeneID = as.integer(GeneID))

load(here('data/human_surface_and_plasma_membrane_protLevel.rda'))
human_surface_plasma_protLevel <- as_tibble(human_surface_and_plasma_membrane_protLevel)
human_surface_plasma_protLevel <- mutate(human_surface_plasma_protLevel, GeneID = as.integer(GeneID))

# read hamster surface
uniprot_cell_surface <- read_tsv(here('data-raw/hamster_cell_surface.txt'))
uniprot_cell_surface <- select(uniprot_cell_surface, "Entry", "Status")
uniprot_cell_surface <- distinct(uniprot_cell_surface)

# read hamster plasma
uniprot_plasma_membrane <- read_tsv(here('data-raw/hamster_plasma_membrane.txt'))
uniprot_plasma_membrane <- select(uniprot_plasma_membrane, "Entry", "Status")
uniprot_plasma_membrane <- distinct(uniprot_plasma_membrane)

# read hamster id mapping
uniprot_id_mapping <- read_tsv(here('data-raw/UP_hamster.idmapping'), col_names = c("UniProt", "ID Type", "ID"))

# create table of hamster genes
hamster_gene_IDS <- filter(uniprot_id_mapping, `ID Type` == "GeneID")
hamster_gene_names <- filter(uniprot_id_mapping, `ID Type` == "Gene_Name")
hamster_genes <- left_join(hamster_gene_names, hamster_gene_IDS, by=c("UniProt"))
hamster_genes <- select(hamster_genes, UniProt, ID.y, ID.x)
names(hamster_genes) <- c("UniProt", "GeneID", "Gene Name")
hamster_genes <- mutate(hamster_genes, GeneID = as.integer(GeneID))
hamster_genes <- distinct(hamster_genes)

# align GO
hamster_genes <- left_join(hamster_genes, uniprot_cell_surface, by=c("UniProt" = "Entry"))
hamster_genes$`GO cell surface` <- !is.na(hamster_genes$Status)
hamster_genes <- select(hamster_genes, -Status)

hamster_genes <- left_join(hamster_genes, uniprot_plasma_membrane, by=c("UniProt" = "Entry"))
hamster_genes$`GO plasma membrane` <- !is.na(hamster_genes$Status)
hamster_genes <- select(hamster_genes, -Status)

# use gene name if no homolog that mapped to a surface/plasma protein
hamster_genes$`Gene Name - human` = toupper(hamster_genes$`Gene Name`)
gene_name_homologs <- left_join(hamster_genes, human_surface_plasma_geneLevel, by = c("Gene Name - human" = "Gene Name"), suffix = c("", " - human"))
# collapse duplicate hamster genes
gene_name_homologs_geneLevel <- as_tibble(setDT(gene_name_homologs)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `Gene Name`])
# collapse duplicate hamster proteins
gene_name_homologs_protLevel <- as_tibble(setDT(gene_name_homologs)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `UniProt`])

convertMulticolumnToSingle <- function(data) {
  data$`GO plasma membrane - human` <- ifelse(is.na(data$`GO plasma membrane - human`) | !str_detect(data$`GO plasma membrane - human`, "TRUE"), F, T)
  data$`GO cell surface - human` <- ifelse(is.na(data$`GO cell surface - human`) | !str_detect(data$`GO cell surface - human`, "TRUE"), F, T)
  data$`GO plasma membrane` <- ifelse(is.na(data$`GO plasma membrane`) | !str_detect(data$`GO plasma membrane`, "TRUE"), F, T)
  data$`GO cell surface` <- ifelse(is.na(data$`GO cell surface`) | !str_detect(data$`GO cell surface`, "TRUE"), F, T)
  data$`Surfaceome` <- ifelse(is.na(data$`Surfaceome`) | !str_detect(data$`Surfaceome`, "TRUE"), F, T)
  data$`CRISPR` <- ifelse(is.na(data$`CRISPR`)  | !str_detect(data$`CRISPR`, "TRUE"), F, T)
  # Compute number of cell surface evidence
  data$`Num cell surface evidence` <- as.integer(data$`CSPA category` != "") +
    as.integer(data$Surfaceome) + 
    as.integer(data$CRISPR) + 
    as.integer(data$`GO cell surface`) + 
    as.integer(data$`GO cell surface - human`)
  # Compute number of plasma membrane evidence
  data$`Num plasma membrane evidence` <- as.integer(data$`GO plasma membrane - human`) + 
    as.integer(data$`GO plasma membrane`)
  return(data)
}

# merge results
hamster_geneLevel <- convertMulticolumnToSingle(gene_name_homologs_geneLevel)
# merge results
hamster_protLevel <- convertMulticolumnToSingle(gene_name_homologs_protLevel)


# Cell Surface and Plasma membrane
hamster_surface_and_plasma_membrane_geneLevel <- filter(hamster_geneLevel, hamster_geneLevel$`Num cell surface evidence` > 0 | hamster_geneLevel$`GO plasma membrane`==T)
hamster_surface_and_plasma_membrane_protLevel <- filter(hamster_protLevel, hamster_protLevel$`Num cell surface evidence` > 0 | hamster_protLevel$`GO plasma membrane`==T)

#write_tsv(hamster_surface_and_plasma_membrane_geneLevel, here("data/hamster_surface_and_plasma_membrane_genes.txt"))
#write_tsv(hamster_surface_and_plasma_membrane_protLevel, here("data/hamster_surface_and_plasma_membrane_proteins.txt"))
usethis::use_data(hamster_surface_and_plasma_membrane_geneLevel, overwrite = TRUE)
usethis::use_data(hamster_surface_and_plasma_membrane_protLevel, overwrite = TRUE)

# Cell Surface
hamster_surface_geneLevel <- filter(hamster_geneLevel, hamster_geneLevel$`Num cell surface evidence` > 0)
hamster_surface_protLevel <- filter(hamster_protLevel, hamster_protLevel$`Num cell surface evidence` > 0)

#write_tsv(hamster_surface_geneLevel, here("data/hamster_surface_genes.txt"))
#write_tsv(hamster_surface_protLevel, here("data/hamster_surface_proteins.txt"))
usethis::use_data(hamster_surface_geneLevel, overwrite = TRUE)
usethis::use_data(hamster_surface_protLevel, overwrite = TRUE)

# Plasma membrane
hamster_plasma_membrane_geneLevel <- filter(hamster_geneLevel, hamster_geneLevel$`GO plasma membrane`==T)
hamster_plasma_membrane_protLevel <- filter(hamster_protLevel, hamster_protLevel$`GO plasma membrane`==T)

#write_tsv(hamster_plasma_membrane_geneLevel, here("data/hamster_plasma_membrane_genes.txt"))
#write_tsv(hamster_plasma_membrane_protLevel, here("data/hamster_plasma_membrane_proteins.txt"))
usethis::use_data(hamster_plasma_membrane_geneLevel, overwrite = TRUE)
usethis::use_data(hamster_plasma_membrane_protLevel, overwrite = TRUE)


