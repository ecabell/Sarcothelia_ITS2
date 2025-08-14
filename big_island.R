
library(dplyr)
library(tidyr)
library(ggplot2)

# Read the file
raw <- read.table("~/Downloads/20250603T195144_cabeler1/post_med_seqs/596_20250604T195602_DBV_20250605T195712.seqs.relative.abund_and_meta.txt", header = TRUE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE)

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

filtered_data <- big_island%>%
  filter(RelAbundance > 0) %>%
  group_by(sample_name) %>%
  arrange(desc(RelAbundance), .by_group = TRUE) %>%
  mutate(ITS2_type = factor(ITS2_type, levels = unique(ITS2_type))) %>%
  ungroup()


label <- filtered_data %>% group_by(sample_name) %>%
  mutate(perc = RelAbundance/ sum(RelAbundance)) %>% 
  mutate(labels = scales::percent(perc)) %>% 
  arrange(desc(ITS2_type)) %>% ## arrange in the order of the legend
  mutate(text_y = cumsum(RelAbundance) - RelAbundance/2) ### ca

label %>%
  ggplot(aes(x = "", y = RelAbundance, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~ sample_name) +
  geom_label_repel(aes(label = paste0(ITS2_type, ": ", round(RelAbundance * 100, 1), "%")),
             nudge_x = 1,size = 3) +
  labs(title = "ITS2 Type Distribution per Sample") +
  scale_fill_viridis_d() +
  theme_void()


label %>%
  ggplot(aes(x = sample_name, y = RelAbundance, fill = ITS2_type)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  labs(title = "ITS2 Type Distribution per Sample", x = "Sample", y = "Relative Abundance", fill = "ITS2 Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(colorspace)

n_types <- length(unique(filtered_data$ITS2_type))
palette_colors <- qualitative_hcl(n_types, palette = "Dark 3")  # or "Set 2", "Harmonic"

# Set ITS2_type as factor with consistent ordering
filtered_data <- filtered_data %>%
  mutate(ITS2_type = factor(ITS2_type, levels = unique(ITS2_type)))

# Pie chart
ggplot(filtered_data, aes(x = "", y = RelAbundance, fill = ITS2_type)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  facet_wrap(~ sample_name) +
  theme_void() +
  geom_text(aes(y = ypos, label = label), size = 3, check_overlap = TRUE) +
  scale_fill_manual(values = palette_colors) +
  labs(title = "ITS2 Type Distribution per Sample")
