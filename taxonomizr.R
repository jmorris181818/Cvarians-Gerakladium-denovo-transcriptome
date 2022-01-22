#### packages ####

# install.packages("taxonomizr")
library(taxonomizr)
library(dplyr)
library(tidyr)
library(tibble)


#### database download ####

# NOTE: this will require a lot of hard drive space, bandwidth, and time to process all the data from NCBI
prepareDatabase('accessionTaxa.sql')


#### data import and taxonomy generation ####

nomatch <- read.csv(file = "../allblast.br",
                         header = F,
                         sep = "\t",
                         stringsAsFactors = F)

# creating a lookup table with contig and taxa IDs generated from blast
nomatch %>%
  dplyr::select(V1, V6) -> taxaID

# generating full taxonomy from taxa IDs
taxonomy <- getTaxonomy(taxaID$V6,'accessionTaxa.sql')

# generating a dataframe of contig ID to full taxonomy
condenseTaxa(taxonomy, groupings = nomatch$V1) -> contig_taxa

# exporting full table for reference
contig_taxa <- as.data.frame(contig_taxa)
write.table(contig_taxa, file="nomatch_taxa.txt",sep="\t", quote = FALSE)


#### subsetting host/symbiont matches ####

symbiont <- subset(contig_taxa, class == "Dinophyceae")
# 75 contigs
symbiont %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") -> symbiont_id

host <- subset(contig_taxa, phylum == "Annelida" | phylum == "Arthopoda" | phylum == "Chordata" | phylum == "Cnidaria" | phylum == "Mollusca" | phylum == "Porifera")
# 507 contigs
host %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") -> host_id

# exporting host/symbiont lookup tables
write.table(symbiont_id$id, file = "nomatch_symbiont.txt",sep="\t", quote = FALSE, row.names = F, col.names = F)
write.table(host_id$id, file = "nomatch_host.txt",sep="\t", quote = FALSE, row.names = F, col.names = F)


# Now return to Cvarians_transcriptome_hostsym_separation_README.txt to finish transcriptome assemblies