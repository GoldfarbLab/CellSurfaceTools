library(tidyverse)
library(here)
library(homologene)
library(data.table)
library(devtools)

#############################################
## Create Gene-level table
#############################################

# read human surface/plasma
human_surface_plasma_geneLevel <- read_tsv(here('data/human_surface_and_plasma_membrane_genes.txt'))
human_surface_plasma_geneLevel <- mutate(human_surface_plasma_geneLevel, GeneID = as.integer(GeneID))

human_surface_plasma_protLevel <- read_tsv(here('data/human_surface_and_plasma_membrane_proteins.txt'))
human_surface_plasma_protLevel <- mutate(human_surface_plasma_protLevel, GeneID = as.integer(GeneID))

# read monkey surface
uniprot_cell_surface <- read_tsv(here('data-raw/monkey_cell_surface.txt'))
uniprot_cell_surface <- select(uniprot_cell_surface, "Entry", "Status")
uniprot_cell_surface <- distinct(uniprot_cell_surface)

# read monkey plasma
uniprot_plasma_membrane <- read_tsv(here('data-raw/monkey_plasma_membrane.txt'))
uniprot_plasma_membrane <- select(uniprot_plasma_membrane, "Entry", "Status")
uniprot_plasma_membrane <- distinct(uniprot_plasma_membrane)

# read monkey id mapping
uniprot_id_mapping <- read_tsv(here('data-raw/UP_monkey.idmapping'), col_names = c("UniProt", "ID Type", "ID"))

# create table of monkey genes
monkey_gene_IDS <- filter(uniprot_id_mapping, `ID Type` == "GeneID")
monkey_gene_names <- filter(uniprot_id_mapping, `ID Type` == "Gene_Name")
monkey_genes <- left_join(monkey_gene_names, monkey_gene_IDS, by=c("UniProt"))
monkey_genes <- select(monkey_genes, UniProt, ID.y, ID.x)
names(monkey_genes) <- c("UniProt", "GeneID", "Gene Name")
monkey_genes <- mutate(monkey_genes, GeneID = as.integer(GeneID))
monkey_genes <- distinct(monkey_genes)

# align GO
monkey_genes <- left_join(monkey_genes, uniprot_cell_surface, by=c("UniProt" = "Entry"))
monkey_genes$`GO cell surface` <- !is.na(monkey_genes$Status)
monkey_genes <- select(monkey_genes, -Status)

monkey_genes <- left_join(monkey_genes, uniprot_plasma_membrane, by=c("UniProt" = "Entry"))
monkey_genes$`GO plasma membrane` <- !is.na(monkey_genes$Status)
monkey_genes <- select(monkey_genes, -Status)

# map monkey genes to human homologs
homologs <- autoTranslate(monkey_genes$GeneID, human_surface_plasma_geneLevel$`Gene Name`, possibleTargets=c(9606))
# align to human annotations
homologs <- left_join(homologs, human_surface_plasma_geneLevel, by = c("9606_ID" = "GeneID"))


# join with monkey genes
monkey_genes_homologene <- left_join(monkey_genes, homologs, by=c("GeneID" = "9544_ID"), suffix = c("", " - human"))
# keep only the ones that actually mapped and are surface/plasma
monkey_genes_homologene <- filter(monkey_genes_homologene, !is.na(`9606_ID`))
monkey_genes_homologene <- rename(monkey_genes_homologene, "GeneID - human" = "9606_ID")
# filter for surface/plasma
monkey_genes_homologene <- select(filter(monkey_genes_homologene, `UniProt` != ""), -`9544`, -`9606` )
# collapse duplicate monkey genes
monkey_genes_homologene_geneLevel <- as_tibble(setDT(monkey_genes_homologene)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `Gene Name`])
# collapse duplicate monkey proteins
monkey_genes_homologene_protLevel <- as_tibble(setDT(monkey_genes_homologene)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `UniProt`])

# use gene name if no homolog that mapped to a surface/plasma protein
gene_name_homologs <- left_join(monkey_genes, homologs, by=c("GeneID" = "9544_ID"), suffix = c("", " - human"))
gene_name_homologs <- select(filter(gene_name_homologs, is.na(`9606_ID`)), UniProt, GeneID, "Gene Name", "GO cell surface", "GO plasma membrane")
gene_name_homologs$`Gene Name - human` = toupper(gene_name_homologs$`Gene Name`)
gene_name_homologs <- left_join(gene_name_homologs, human_surface_plasma_geneLevel, by = c("Gene Name - human" = "Gene Name"), suffix = c("", " - human"))
# collapse duplicate monkey genes
gene_name_homologs_geneLevel <- as_tibble(setDT(gene_name_homologs)[, lapply(.SD, function(x) paste(unique(na.omit(x)), collapse = ";")), by = `Gene Name`])
# collapse duplicate monkey proteins
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
monkey_geneLevel <- rbind(monkey_genes_homologene_geneLevel, gene_name_homologs_geneLevel)
monkey_geneLevel <- convertMulticolumnToSingle(monkey_geneLevel)
# merge results
monkey_protLevel <- rbind(monkey_genes_homologene_protLevel, gene_name_homologs_protLevel)
monkey_protLevel <- convertMulticolumnToSingle(monkey_protLevel)


# Cell Surface and Plasma membrane
monkey_surface_and_plasma_membrane_geneLevel <- filter(monkey_geneLevel, monkey_geneLevel$`Num cell surface evidence` > 0 | monkey_geneLevel$`GO plasma membrane`==T)
monkey_surface_and_plasma_membrane_protLevel <- filter(monkey_protLevel, monkey_protLevel$`Num cell surface evidence` > 0 | monkey_protLevel$`GO plasma membrane`==T)

write_tsv(monkey_surface_and_plasma_membrane_geneLevel, here("data/monkey_surface_and_plasma_membrane_genes.txt"))
write_tsv(monkey_surface_and_plasma_membrane_protLevel, here("data/monkey_surface_and_plasma_membrane_proteins.txt"))
usethis::use_data(monkey_surface_and_plasma_membrane_geneLevel, overwrite = TRUE)
usethis::use_data(monkey_surface_and_plasma_membrane_protLevel, overwrite = TRUE)

# Cell Surface
monkey_surface_geneLevel <- filter(monkey_geneLevel, monkey_geneLevel$`Num cell surface evidence` > 0)
monkey_surface_protLevel <- filter(monkey_protLevel, monkey_protLevel$`Num cell surface evidence` > 0)

write_tsv(monkey_surface_geneLevel, here("data/monkey_surface_genes.txt"))
write_tsv(monkey_surface_protLevel, here("data/monkey_surface_proteins.txt"))
usethis::use_data(monkey_surface_geneLevel, overwrite = TRUE)
usethis::use_data(monkey_surface_protLevel, overwrite = TRUE)

# Plasma membrane
monkey_plasma_membrane_geneLevel <- filter(monkey_geneLevel, monkey_geneLevel$`GO plasma membrane`==T)
monkey_plasma_membrane_protLevel <- filter(monkey_protLevel, monkey_protLevel$`GO plasma membrane`==T)

write_tsv(monkey_plasma_membrane_geneLevel, here("data/monkey_plasma_membrane_genes.txt"))
write_tsv(monkey_plasma_membrane_protLevel, here("data/monkey_plasma_membrane_proteins.txt"))
usethis::use_data(monkey_plasma_membrane_geneLevel, overwrite = TRUE)
usethis::use_data(monkey_plasma_membrane_protLevel, overwrite = TRUE)


