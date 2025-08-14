library(tidyverse)
library(viridis)
library(dplyr)

# Read the full file
raw <- read.table("~/Downloads/20250603T195144_cabeler1/its2_type_profiles/596_20250604T195602_DBV_20250605T195712.profiles.relative.abund_and_meta.txt",
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

# install if needed
# install.packages("ggstatsplot")

library(ggstatsplot)

data_long$Abundance <- as.numeric(data_long$Abundance)  # Ensure numeric
# Create a new integer count column (e.g., by scaling and rounding)
data_long <- data_long %>%
  mutate(IntAbundance = round(Abundance * 1000)) %>%
  filter(!is.na(IntAbundance), IntAbundance > 0)

grouped_ggpiestats(
  data = data_long,
  x = ITS2_type,
  counts = IntAbundance,
  grouping.var = Location,
  label = "percentage",
  label.repel = TRUE,
  ggplot.component = ggplot2::scale_fill_viridis_d(option = "cividis"),
  ggtheme = theme_minimal(),
  results.subtitle = TRUE,
  title.text = "ITS2 Type Distributions by Location",
  plotgrid.args = list(nrow = 1)
)

grouped_ggpiestats(
  data = data_long,
  x = ITS2_type,
  counts = IntAbundance,
  grouping.var = Morph,
  label = "percentage",
  label.repel = TRUE,
  ggplot.component = ggplot2::scale_fill_viridis_d(option = "cividis"),
  ggtheme = theme_minimal(),
  results.subtitle = TRUE,
  title.text = "ITS2 Type Distributions by Morph",
  plotgrid.args = list(nrow = 1)
)

# Prepare data: convert proportions to pseudo-counts
data_counts <- data_long %>%
  mutate(IntCount = round(RelAbundance * 1000)) %>%
  filter(IntCount > 0)  # remove zero counts

# Use ggbarstats with simulated count data
ggbarstats(
  data = data_counts,
  x = ITS2_type,
  y = Morph,              # y must be categorical, so switch with Location
  weight = IntCount,         # this simulates the count data ggbarstats needs
  label = "both",
  facet.vars = Location,   # Facet by Location
  results.subtitle = TRUE,
  ggtheme = theme_minimal(),
  ggplot.component = ggplot2::scale_fill_viridis_d(option = "cividis")
)

