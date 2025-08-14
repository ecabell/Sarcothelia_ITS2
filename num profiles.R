library(dplyr)
library(tidyr)

# Load your dataset
# Read the full file
raw <- read.table("~/Documents/Zoox PCR/20250603T195144_cabeler1/its2_type_profiles/596_20250604T195602_DBV_20250605T195712.profiles.absolute.abund_and_meta.txt",
                  header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

# Subset the portion with the ITS2 data (starting from row 7 onwards)
its2_data <- raw[7:nrow(raw), ]

# Assign first row as column names for ITS2 types
colnames(its2_data) <- c("RowID", "Sample", as.character(unlist(its2_data[1, 3:ncol(its2_data)])))

# Remove the first row (used as column names)
its2_data <- its2_data[-1, ]

# Remove RowID column and convert remaining numeric columns
its2_data <- its2_data %>%
  select(-RowID) %>%
  mutate(across(-Sample, as.numeric))

# Extract morph and location from sample names
its2_data <- its2_data %>%
  mutate(
    Location = case_when(
      grepl("^K", Sample) ~ "Kaua`i",
      grepl("^O", Sample) ~ "O`ahu",
      grepl("^H", Sample) ~ "Big Island",
      TRUE ~ "Unknown"
    ),
    Morph = case_when(
      grepl("BR", Sample) ~ "Brown",
      grepl("BL", Sample) ~ "Blue",
      TRUE ~ "Unknown"
    )
  )

# Filter out samples with unknown morph or location (i.e., likely controls)
its2_data <- its2_data %>%
  filter(Morph != "Unknown", Location != "Unknown")

# Pivot to long format for ggplot
data_long <- its2_data %>%
  pivot_longer(
    cols = -c(Sample, Location, Morph),
    names_to = "ITS2_type",
    values_to = "Abundance"
  )

# Calculate relative abundance within each sample
data_long <- data_long %>%
  group_by(Sample) %>%
  mutate(Total = sum(Abundance, na.rm = TRUE),
         RelAbundance = Abundance / Total) %>%
  ungroup()

# Summarize by Morph and Location
summary_data <- data_long %>%
  group_by(Location, Morph, ITS2_type) %>%
  summarize(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

# Keep only profiles with non-zero relative abundance
df_present <- summary_data %>% filter(MeanRelAbundance > 0)

all_profiles <- unique(df_present$ITS2_type)
length(all_profiles) #14

# Extract ITS2 profiles per island
kauai_profiles <- unique(df_present$ITS2_type[df_present$Location == "Kaua`i"])
oahu_profiles <- unique(df_present$ITS2_type[df_present$Location == "O`ahu"])
big_island_profiles <- unique(df_present$ITS2_type[df_present$Location == "Big Island"])

length(kauai_profiles) #5
length(oahu_profiles) #7
length(big_island_profiles) #4

# Shared profiles across all three islands
shared_all <- Reduce(intersect, list(kauai_profiles, oahu_profiles, big_island_profiles))
shared_k_o <- intersect(oahu_profiles, kauai_profiles)
shared_b_o <- intersect(oahu_profiles, big_island_profiles)
shared_k_b <- intersect(kauai_profiles, big_island_profiles)

# Output
length(shared_all)  # Number of shared profiles
shared_all           # Profile names

only_oahu <- setdiff(oahu_profiles, union(big_island_profiles, kauai_profiles))
only_kauai <- setdiff(kauai_profiles, union(big_island_profiles, oahu_profiles))
only_bi <- setdiff(big_island_profiles, union(kauai_profiles, oahu_profiles))

length(only_oahu)
length(only_kauai)
length(only_bi)

# Create vectors for each group
kauai_blue  <- df_present %>% filter(Location == "Kaua`i", Morph == "Blue") %>% pull(ITS2_type) %>% unique()
kauai_brown <- df_present %>% filter(Location == "Kaua`i", Morph == "Brown") %>% pull(ITS2_type) %>% unique()
oahu_blue   <- df_present %>% filter(Location == "O`ahu", Morph == "Blue") %>% pull(ITS2_type) %>% unique()
oahu_brown  <- df_present %>% filter(Location == "O`ahu", Morph == "Brown") %>% pull(ITS2_type) %>% unique()
bi_blue     <- df_present %>% filter(Location == "Big Island", Morph == "Blue") %>% pull(ITS2_type) %>% unique()

length(kauai_blue) #4
length(kauai_brown) #1
length(oahu_blue) #5
length(oahu_brown) #2
length(bi_blue) #4


only_kbl <- setdiff(kauai_blue, union(union(union(oahu_blue, oahu_brown), bi_blue), kauai_brown))
length(only_kbl)

only_kbr <- setdiff(kauai_brown, union(union(union(oahu_blue, oahu_brown), bi_blue), kauai_blue))
length(only_kbr)

only_obl <- setdiff(oahu_blue, union(union(union(kauai_brown, oahu_brown), bi_blue), kauai_blue))
length(only_obl)

only_obr <- setdiff(oahu_brown, union(union(union(kauai_brown, oahu_blue), bi_blue), kauai_blue))
length(only_obr)

blue <- df_present %>% filter(Morph == 'Blue') %>% pull(ITS2_type) %>%
  unique()
length(blue) #13

brown <- df_present %>% filter(Morph == 'Brown') %>% pull(ITS2_type) %>%
  unique()
length(brown) #3

brown_only <- setdiff(brown, blue)
length(brown_only)
blue_only <- setdiff(blue, brown)
length(blue_only)
shared_morphs <- intersect(blue, brown)

data_long <- data_long %>% filter(RelAbundance > 0)

# Split by sample
profiles_by_sample <- split(data_long$ITS2_type, data_long$Sample)

# Find profiles present in all samples
shared_profiles_all <- Reduce(intersect, profiles_by_sample)

# View result
print(shared_profiles_all)
length(shared_profiles_all)
