library(dplyr)
library(tidyr)

# Load post-MED dataset
postmed <- read.csv("~/Documents/Zoox PCR/20250603T195144_cabeler1/post_med_seqs/596_20250604T195602_DBV_20250605T195712.seqs.absolute.abund_and_meta.csv")
postmed <- postmed[-34,]
postmed <- postmed %>% subset(sample_name != 'HBL4')

# Annotate Morph and Location
postmed <- postmed %>%
  mutate(
    Location = case_when(
      grepl("^K", sample_name) ~ "Kaua`i",
      grepl("^O", sample_name) ~ "O`ahu",
      grepl("^H", sample_name) ~ "Big Island",
      TRUE ~ "Unknown"
    ),
    Morph = case_when(
      grepl("BR", sample_name) ~ "Brown",
      grepl("BL", sample_name) ~ "Blue",
      TRUE ~ "Unknown"
    )
  )

# After adding Location and Morph
metadata_cols <- c("sample_uid", "sample_name", "fastq_fwd_file_name", "fastq_fwd_sha256_file_hash",
                   "fastq_rev_file_name", "fastq_rev_sha256_file_hash", "data_set_uid", "data_set_name",
                   "raw_contigs", "post_qc_absolute_seqs", "post_qc_unique_seqs",
                   "post_taxa_id_absolute_symbiodiniaceae_seqs", "post_taxa_id_unique_symbiodiniaceae_seqs",
                   "size_screening_violation_absolute", "size_screening_violation_unique",
                   "post_taxa_id_absolute_non_symbiodiniaceae_seqs", "post_taxa_id_unique_non_symbiodiniaceae_seqs",
                   "post_med_absolute", "post_med_unique", "sample_type", "host_phylum", "host_class",
                   "host_order", "host_family", "host_genus", "host_species", "collection_latitude",
                   "collection_longitude", "collection_date", "collection_depth",
                   "Location", "Morph")  # <- Add these here

type_cols <- setdiff(names(postmed), metadata_cols) #returns difference, ie cols not in metadata_cols

# Convert to long format for type abundance
postmed_long <- postmed %>%
  pivot_longer(cols = all_of(type_cols), names_to = "type", values_to = "abundance") %>%
  filter(abundance > 0)

# Count unique post-MED types
n_distinct(postmed_long$type)
length(unique(postmed_long$type))

# By Morph
types_by_morph <- postmed_long %>%
  group_by(Morph) %>%
  summarise(unique = n_distinct(type))

# By Location
types_by_location <- postmed_long %>%
  group_by(Location) %>%
  summarise(unique = n_distinct(type))

postmed_long %>%
  group_by(Morph) %>%
  summarise(unique_types = list(unique(type)))

# Find which types are present in each island
kauai_types <- type_cols[colSums(postmed[postmed$Location == "Kaua`i", type_cols] > 0) > 0]
oahu_types <- type_cols[colSums(postmed[postmed$Location == "O`ahu", type_cols] > 0) > 0]
hawaii_types <- type_cols[colSums(postmed[postmed$Location == 'Big Island', type_cols] > 0) > 0]

shared_all <- Reduce(intersect, list(kauai_types, oahu_types, hawaii_types))
shared_b_o <- intersect(oahu_types, hawaii_types)
shared_k_b <- intersect(kauai_types, hawaii_types)
shared_k_o <- intersect(kauai_types, oahu_types)

length(shared_all) #4
length(shared_b_o) #7
length(shared_k_b) #7
length(shared_k_o) #10

# Compare
only_hawaii <- setdiff(hawaii_types, union(kauai_types, oahu_types))
length(only_hawaii) #37
only_kauai <- setdiff(kauai_types, union(hawaii_types, oahu_types))
length(only_kauai) #46
only_oahu <- setdiff(oahu_types, union(hawaii_types, kauai_types))
length(only_oahu) #46

# Define function to extract types for a given location and morph
get_types_by_group <- function(df, location, morph, type_cols) {
  subset_df <- df %>% filter(Location == location, Morph == morph)
  present <- colSums(subset_df[, type_cols] > 0) > 0
  type_cols[present]
}

# Example usage:
kauai_blue_types <- get_types_by_group(postmed, "Kaua`i", "Blue", type_cols)
kauai_brown_types <- get_types_by_group(postmed, "Kaua`i", "Brown", type_cols)
oahu_blue_types <- get_types_by_group(postmed, "O`ahu", "Blue", type_cols)
oahu_brown_types <- get_types_by_group(postmed, "O`ahu", "Brown", type_cols)

length(kauai_blue_types) #45
length(kauai_brown_types) #18
length(oahu_blue_types) #51
length(oahu_brown_types) #15
length(hawaii_types) #47 #hawaii blue
n_distinct(hawaii_types)

# Compare sets
intersect(kauai_blue_types, oahu_blue_types)
setdiff(kauai_blue_types, oahu_blue_types)
setdiff(kauai_brown_types, oahu_brown_types)

total <- postmed %>% summarise(sum(post_taxa_id_absolute_symbiodiniaceae_seqs))

total_by_morph <- postmed %>% group_by(Morph) %>% 
  summarise(sum(post_taxa_id_absolute_symbiodiniaceae_seqs))

avg_unique_morph <- postmed %>% group_by(Morph) %>%
  summarise(mean(post_taxa_id_unique_symbiodiniaceae_seqs))


n_distinct(only_kauai)
n_distinct(only_oahu)
n_distinct(oahu_types)
n_distinct(kauai_types)


# Split by sample
types_by_sample <- split(postmed_long$type, postmed_long$sample_name)

# Find profiles present in all samples
shared_types_all <- Reduce(intersect, types_by_sample)

# View result
print(shared_types_all)
length(shared_types_all)

# ITS2 types only in Kaua`i Blue
only_kbl <- setdiff(kauai_blue_types, union(union(union(oahu_blue_types, oahu_brown_types), hawaii_types), kauai_brown_types))
length(only_kbl)

# ITS2 types only in Kaua`i Brown
only_kbr <- setdiff(kauai_brown_types, union(union(union(oahu_blue_types, oahu_brown_types), hawaii_types), kauai_blue_types))
length(only_kbr)

only_obl <- setdiff(oahu_blue_types, union(union(union(kauai_brown_types, oahu_brown_types), hawaii_types), kauai_blue_types))
length(only_obl)

only_obr <- setdiff(oahu_brown_types, union(union(union(kauai_brown_types, oahu_blue_types), hawaii_types), kauai_blue_types))
length(only_obr)

blue <- postmed_long %>% filter(Morph == 'Blue') %>% pull(type) %>%
  unique()
length(blue) #129

brown <- postmed_long %>% filter(Morph == 'Brown') %>% pull(type) %>%
  unique()
length(brown) #30

brown_only <- setdiff(brown, blue)
length(brown_only)
blue_only <- setdiff(blue, brown)
length(blue_only)
shared_morphs <- intersect(blue, brown)
