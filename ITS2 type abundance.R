# Load libraries
library(tidyverse)
library(viridis)
library(dplyr)
library(ggrepel)

# Read the full file
raw <- read.table("~/Documents/Zoox PCR/20250603T195144_cabeler1/its2_type_profiles/596_20250604T195602_DBV_20250605T195712.profiles.relative.abund_and_meta.txt",
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

its2_data<- its2_data %>% slice(1:(n() - 2))# Extract morph and location from sample names
its2_data <- its2_data %>%
  mutate(
    Location = case_when(
      grepl("^K", Sample) ~ "Kaua`i",
      grepl("^O", Sample) ~ "O`ahu",
      grepl("^H", Sample) ~ "Hawai`i Island",
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
    names_to = "ITS2_Profile",
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
  group_by(Location, Morph, ITS2_Profile) %>%
  summarize(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE),
            Abundance = mean(Abundance), .groups = "drop")

# Plot stacked bar chart
ggplot(summary_data, aes(x = interaction(Location, Morph), y = Abundance, fill = ITS2_Profile)) +
  geom_bar(stat = "identity") +
  labs(x = "Location-Morph", y = "Relative Abundance", fill = "ITS2 Profile") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d() +
  geom_text(aes(label = round(MeanRelAbundance, 2)), 
            position = position_stack(vjust = 0.5), 
            color = "white")

# Plot stacked bar chart
ggplot(summary_data, aes(x = Morph, y = Abundance, fill = ITS2_Profile)) +
  geom_bar(stat = "identity") +
  labs(x = "Morph", y = "Relative Abundance", fill = "ITS2 Profile") +
  facet_wrap(~Location) +
  scale_fill_viridis_d() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot pie chart
ggplot(summary_data, aes(x = "", y = MeanRelAbundance, fill = ITS2_Profile)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_wrap(~Location + Morph) +
  labs(title = "ITS2 Profile Relative Abundance by Morph and Location",
       fill = "ITS2 Profile") +
  scale_fill_viridis_d() +
  theme_void()   # Removes axes and background for pie charts

# Plot stacked bar chart
ggplot(summary_data, aes(x = Location, y = MeanRelAbundance, fill = ITS2_Profile)) +
  geom_bar(stat = "identity") +
  labs(x = "Location", y = "Relative Abundance", fill = "ITS2 Profile") +
  facet_wrap(~Morph) +
  scale_fill_viridis_d(option = "cividis") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

summary_data$Morph <- factor(summary_data$Morph, levels = c('Blue', 'Brown'))
summary_data$Location <- factor(summary_data$Location, levels = c("Hawai`i Island", 'Kaua`i', 'O`ahu'))

# Plot pie chart
ggplot(summary_data, aes(x = "", y = MeanRelAbundance, fill = ITS2_Profile)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(x = "Location", y = "Relative Abundance", fill = "ITS2 Profile") +
  facet_grid(rows = vars(Morph), cols = vars(Location), switch = "y") +
  scale_fill_viridis_d(option = "mako") +
  labs(title = "Symbiodiniaceae ITS2 Profile Relative Abundance by Morph and Location",
       fill = "ITS2 Profile") +
  theme_void() + 
  theme(strip.text = element_text(size = 12), # ← Change size here
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)), 
        legend.position = "right")  #center title and add space below it 

# Plot pie chart
ggplot(summary_data, aes(x = "", y = MeanRelAbundance, fill = ITS2_Profile)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  labs(x = "Location", y = "Relative Abundance", fill = "ITS2 Profile") +
  facet_grid(rows = vars(Morph), cols = vars(Location), switch = "y") +
  scale_fill_viridis_d() +
  labs(title = "Symbiodiniaceae ITS2 Profile Relative Abundance by Morph and Location",
       fill = "ITS2 Profile") +
  theme_void() + 
  theme(strip.text = element_text(size = 12), # ← Change size here
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)), 
        legend.position = "right")  #center title and add space below it 

# Add group column based on ITS2_type
data_long <- data_long %>%
  mutate(Group = case_when(
    grepl("^A", ITS2_Profile) ~ "Symbiodinium",
    grepl("^B", ITS2_Profile) ~ "Breviolum",
    grepl("^C", ITS2_Profile) ~ "Cladocopium",
    grepl("^F", ITS2_Profile) ~ "Fugacium",
    TRUE ~ "Other"
  ))

# Summarize by group
summary_grouped <- data_long %>%
  group_by(Location, Morph, Group) %>%
  summarize(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

ggplot(summary_grouped, aes(x = Morph, y = MeanRelAbundance, fill = Group)) +
  geom_bar(stat = "identity") +
  labs(x = "Morph", y = "Relative Abundance", fill = "ITS2 Type") +
  facet_wrap(~Location, drop = TRUE) +
  scale_fill_viridis_d()+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Sum abundances by group instead of averaging
summary_grouped <- data_long %>%
  group_by(Location, Morph, Group) %>%
  summarise(SumRelAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  group_by(Location, Morph) %>%
  mutate(NormalizedAbundance = SumRelAbundance / sum(SumRelAbundance)) %>%
  ungroup()

summary_grouped$Group <- factor(summary_grouped$Group, levels = c('Symbiodinium', 'Breviolum', 'Cladocopium', 'Fugacium'))

# Plot corrected pie chart
ggplot(summary_grouped, aes(x = "", y = NormalizedAbundance, fill = Group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_grid(Morph ~ Location) +
  theme_void() +
  scale_fill_viridis_d()+
  labs(title = "Symbiodiniaceae Community Composition by Morph and Location ",
       fill = "Genus")

summary_grouped <- summary_grouped %>%
  group_by(Location, Morph) %>%
  mutate(
    Fraction = NormalizedAbundance / sum(NormalizedAbundance),
    ymax = cumsum(Fraction),
    ymin = lag(ymax, default = 0),
    labelPosition = (ymax + ymin) / 2,
    labelText = paste0("italic('", Group, "')*'  ", round(Fraction * 100, 1), "%'")
  ) %>%
  ungroup()

ggplot(summary_grouped, aes(x = "", y = NormalizedAbundance, fill = Group)) +
  geom_col(width = 1, colour = "white", size = 0.1, alpha = 0.85) +
  coord_polar("y") +
  geom_text_repel(data = summary_grouped %>% filter(Fraction > 0),
                  aes(y = labelPosition, label = labelText, x = 1.1), size = 4, parse = TRUE,
                  nudge_x = 1,
                  direction = "y",
                  segment.size = 0.05,
                  show.legend = FALSE)+
  facet_grid(rows = vars(Morph), cols = vars(Location), switch = "y") +
  theme_void() +
  scale_fill_viridis_d(
   # option = "cividis",
    labels = c(
      "Cladocopium" = expression(italic("Cladocopium")),
      "Breviolum" = expression(italic("Breviolum")),
      "Symbiodinium" = expression(italic("Symbiodinium")),
      "Fugacium" = expression(italic("Fugacium"))
      
    ) ) +
  labs(title = "Symbiodiniaceae Community Composition (Genus)\nby Morph and Location",
       fill  = "Genus") +
  theme(strip.text = element_text(size = 13), # ← Change size here
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5, size = 14, margin = margin(b = 10)), 
        legend.position = "right")  #center title and add space below it 

# Plot bar chart
ggplot(summary_grouped, aes(x = "", y = NormalizedAbundance, fill = Group)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(Morph ~ Location) +
  scale_fill_viridis_d()+
  theme_minimal()+
  labs(title = "Symbiodiniaceae Community Composition by Morph and Location ",
       fill = "Genus")

summary_data$ITS2_Profile <- factor(summary_data$ITS2_Profile, levels = unique(summary_data$ITS2_Profile))
# ─────────────────────────────────────────────────────────────
# Build a facet-by-facet pie-chart data set
# ─────────────────────────────────────────────────────────────
summary_pie <- summary_data %>%                       # your original table
  group_by(Location, Morph) %>%                       # one pie per facet
  arrange(desc(ITS2_Profile), .by_group = TRUE) %>%   # stable slice order
  mutate(
    Fraction = MeanRelAbundance / sum(MeanRelAbundance),   # slice height
    ymax      = cumsum(Fraction),                          # top of slice
    ymin      = lag(ymax, default = 0),                    # bottom of slice
    mid       = (ymax + ymin) / 2,                         # label position
    label     = paste0("'", ITS2_Profile, "'  ",
                       round(Fraction * 100, 1), "%")      # pretty label
  ) %>%
  ungroup()

# ─────────────────────────────────────────────────────────────
# Draw the pie charts with correctly centred labels
# ─────────────────────────────────────────────────────────────
ggplot(summary_pie,
       aes(x = 1, y = Fraction, fill = ITS2_Profile)) +
  geom_col(width = 1, colour = "white", size = 0.1, alpha = 0.85) +
  geom_text_repel(data =summary_pie %>% filter(Fraction > 0.007), 
                  aes(y = mid, label = label, x = 1.1),
                  parse = FALSE,      # allows the italicised/quoted bits
                  size  = 4,
                  segment.size = 0.05,
                  nudge_x = 1,
                  colour = "black") +
  coord_polar(theta = "y") +
  facet_grid(rows = vars(Morph),
             cols  = vars(Location),
             switch = "y") +
  scale_fill_viridis_d(option = "mako") +
  labs(title = "Symbiodiniaceae Community Composition (ITS2 Profiles)\nby Morph and Location",
       fill  = "ITS2 Profile") +
  theme_void() +
  theme(strip.text      = element_text(size = 13),
        strip.placement = "outside",
        plot.title      = element_text(hjust = .5, size = 14,
                                       margin = margin(b = 10)),
        legend.position = "right")

