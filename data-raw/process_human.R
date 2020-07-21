library(tidyverse)
library(here)
library(devtools)

#############################################
## Create Gene-level table
#############################################

# Read in data
uniprot_id_mapping <- read_tsv(here('data-raw/UP_human_reviewed.idmapping'))

aebersold <- read_tsv(here('data-raw/human_cellsurface_aebersold.txt'))
aebersold <- select(aebersold, "CSPA category", "ID_link")
aebersold <- distinct(aebersold)

PNAS <- read_tsv(here('data-raw/human_cellsurface_PNAS.txt'))
PNAS <- select(PNAS, "Surfaceome Label", "UniProt accession")
PNAS <- distinct(PNAS)

uniprot_cell_surface <- read_tsv(here('data-raw/human_cellsurface_uniprot.tsv'))
uniprot_cell_surface <- select(uniprot_cell_surface, "Entry", "Status")
uniprot_cell_surface <- distinct(uniprot_cell_surface)

uniprot_plasma_membrane <- read_tsv(here('data-raw/human_plasmamembrane_uniprot.tsv'))
uniprot_plasma_membrane <- select(uniprot_plasma_membrane, "Entry", "Status")
uniprot_plasma_membrane <- distinct(uniprot_plasma_membrane)

CRISPR_library <- read_tsv(here('data-raw/ThermoFisher_LentiArray CRISPR Library_Human Cell Surface_Gene List.txt'))
CRISPR_library <- select(CRISPR_library, "gene_symbol", "gene_id")
CRISPR_library <- distinct(CRISPR_library)

# Create table with gene names and gene_ids
human_gene_IDS <- filter(uniprot_id_mapping, `ID Type` == "GeneID")
human_gene_names <- filter(uniprot_id_mapping, `ID Type` == "Gene_Name")
human_genes <- left_join(human_gene_names, human_gene_IDS, by=c("UniProt"))
human_genes <- select(human_genes, UniProt, ID.y, ID.x)
names(human_genes) <- c("UniProt", "GeneID", "Gene Name")
human_genes <- mutate(human_genes, GeneID = as.integer(GeneID))
human_genes <- distinct(human_genes)

# Align aebersold
human_genes <- left_join(human_genes, aebersold, by=c("UniProt" = "ID_link"))

# Align PNAS Surfaceome
human_genes <- left_join(human_genes, PNAS, by=c("UniProt" = "UniProt accession"))
human_genes$Surfaceome <- !is.na(human_genes$`Surfaceome Label`)
human_genes <- select(human_genes, -`Surfaceome Label`)

# Align CRISPR
human_genes <- left_join(human_genes, CRISPR_library, by=c("GeneID" = "gene_id"))
human_genes$CRISPR <- !is.na(human_genes$gene_symbol)
human_genes <- select(human_genes, -gene_symbol)

# Align GO cell surface
human_genes <- left_join(human_genes, uniprot_cell_surface, by=c("UniProt" = "Entry"))
human_genes$`GO cell surface` <- !is.na(human_genes$Status)
human_genes <- select(human_genes, -Status)

# Align GO plasma membrane
human_genes <- left_join(human_genes, uniprot_plasma_membrane, by=c("UniProt" = "Entry"))
human_genes$`GO plasma membrane` <- !is.na(human_genes$Status)
human_genes <- select(human_genes, -Status)


convertMulticolumnToSingle <- function(data) {
  data$`GO plasma membrane` <- ifelse(is.na(data$`GO plasma membrane`) | !str_detect(data$`GO plasma membrane`, "TRUE"), F, T)
  data$`GO cell surface` <- ifelse(is.na(data$`GO cell surface`) | !str_detect(data$`GO cell surface`, "TRUE"), F, T)
  data$`Surfaceome` <- ifelse(is.na(data$`Surfaceome`) | !str_detect(data$`Surfaceome`, "TRUE"), F, T)
  data$`CRISPR` <- ifelse(is.na(data$`CRISPR`)  | !str_detect(data$`CRISPR`, "TRUE"), F, T)
  # Compute number of cell surface evidence
  data$`Num cell surface evidence` <- as.integer(data$`CSPA category` != "") +
    as.integer(data$Surfaceome) + 
    as.integer(data$CRISPR) + 
    as.integer(data$`GO cell surface`)
  return(data)
}

# collapse to unique gene names
human_geneLevel <- as_tibble(setDT(human_genes)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `Gene Name`])
human_geneLevel <- convertMulticolumnToSingle(human_geneLevel)

# collapse to unique proteins
human_protLevel <- as_tibble(setDT(human_genes)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `UniProt`])
human_protLevel <- convertMulticolumnToSingle(human_protLevel)


# Cell Surface and Plasma membrane
human_surface_and_plasma_membrane_geneLevel <- filter(human_geneLevel, human_geneLevel$`Num cell surface evidence` > 0 | human_geneLevel$`GO plasma membrane`==T)
human_surface_and_plasma_membrane_protLevel <- filter(human_protLevel, human_protLevel$`Num cell surface evidence` > 0 | human_protLevel$`GO plasma membrane`==T)

write_tsv(human_surface_and_plasma_membrane_geneLevel, here("data/human_surface_and_plasma_membrane_genes.txt"))
write_tsv(human_surface_and_plasma_membrane_protLevel, here("data/human_surface_and_plasma_membrane_proteins.txt"))

# Cell Surface
human_surface_geneLevel <- filter(human_geneLevel, human_geneLevel$`Num cell surface evidence` > 0)
human_surface_protLevel <- filter(human_protLevel, human_protLevel$`Num cell surface evidence` > 0)

write_tsv(human_surface_geneLevel, here("data/human_surface_genes.txt"))
write_tsv(human_surface_protLevel, here("data/human_surface_proteins.txt"))
usethis::use_data(human_surface_geneLevel, overwrite = TRUE)
usethis::use_data(human_surface_protLevel, overwrite = TRUE)

# Plasma membrane
human_plasma_membrane_geneLevel <- filter(human_geneLevel, human_geneLevel$`GO plasma membrane`==T)
human_plasma_membrane_protLevel <- filter(human_protLevel, human_protLevel$`GO plasma membrane`==T)

write_tsv(human_plasma_membrane_geneLevel, here("data/human_plasma_membrane_genes.txt"))
write_tsv(human_plasma_membrane_protLevel, here("data/human_plasma_membrane_proteins.txt"))
usethis::use_data(human_plasma_membrane_geneLevel, overwrite = TRUE)
usethis::use_data(human_plasma_membrane_protLevel, overwrite = TRUE)

# CRISPR
human_CRISPR_geneLevel <- filter(human_geneLevel, human_geneLevel$`CRISPR` == T)
human_CRISPR_protLevel <- filter(human_protLevel, human_protLevel$`CRISPR` == T)

write_tsv(human_CRISPR_geneLevel, here("data/human_CRISPR_genes.txt"))
write_tsv(human_CRISPR_protLevel, here("data/human_CRISPR_proteins.txt"))
usethis::use_data(human_CRISPR_geneLevel, overwrite = TRUE)
usethis::use_data(human_CRISPR_protLevel, overwrite = TRUE)



