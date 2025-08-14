library(dplyr)
library(tidyr)
library(ggplot2)

# Read the file
raw <- read.table("~/Documents/Zoox PCR/20250603T195144_cabeler1/post_med_seqs/596_20250604T195602_DBV_20250605T195712.seqs.relative.abund_and_meta.txt", header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

# Identify the start of ITS2 data columns
# Assuming ITS2 types start from the column after "collection_depth"
its2_start_col <- which(colnames(raw) == "collection_depth") + 1

# Select relevant columns
its2_data <- raw %>%
  select(sample_name, all_of(colnames(raw)[its2_start_col:ncol(raw)]))

# Convert ITS2 columns to numeric
its2_data <- its2_data %>%
  mutate(across(-sample_name, as.numeric))

# Extract Location and Morph
its2_data<- its2_data %>%
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

# Filter out unknowns
its2_data <- its2_data %>%
  filter(Morph != "Unknown", Location != "Unknown")

# Pivot longer
data_long <- its2_data %>%
  pivot_longer(
    cols = -c(sample_name, Location, Morph),
    names_to = "ITS2_type",
    values_to = "Abundance"
  )

# Calculate relative abundance
data_long <- data_long %>%
  group_by(sample_name) %>%
  mutate(Total = sum(Abundance, na.rm = TRUE),
         RelAbundance = Abundance / Total) %>%
  ungroup()

big_island <- data_long %>% subset(Location == 'Big Island')

# Summarize
summary_data <- data_long %>%
  group_by(Location, Morph, ITS2_type) %>%
  summarize(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

# Plot
ggplot(summary_data, aes(x = "", y = MeanRelAbundance, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_grid(Morph ~ Location) +
  theme_void() +
  labs(title = "Mean Relative Abundance of ITS2 Types by Morph and Location", fill = "ITS2 Type")

