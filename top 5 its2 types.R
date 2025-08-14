# Load libraries
library(tidyverse)
library(viridis)
library(dplyr)

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

top5_all <- data_long %>%
  group_by(ITS2_type) %>%
  summarize(TotalAbundance = sum(RelAbundance, na.rm = TRUE)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 5)
top5_all$percent <- round((top5_all$TotalAbundance / 27.65297) * 100,)
print("Top 5 ITS2 types overall:")
print(top5_all)
print(sum(top5_all$TotalAbundance))

top5_by_morph <- data_long %>%
  group_by(Morph, ITS2_type) %>%
  summarize(TotalAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Morph, desc(TotalAbundance)) %>%
  group_by(Morph) %>%
  slice_head(n = 5)

top5_by_morph <- top5_by_morph %>% group_by(Morph) %>% mutate(percent = round(TotalAbundance/sum(TotalAbundance) * 100), )

print("Top 5 ITS2 types by Morph:")
print(top5_by_morph)

top5_by_location <- data_long %>%
  group_by(Location, ITS2_type) %>%
  summarize(TotalAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Location, desc(TotalAbundance)) %>%
  group_by(Location) %>%
  slice_head(n = 5)

top5_by_location <- top5_by_location %>% group_by(Location) %>% mutate(percent = round(TotalAbundance/sum(TotalAbundance) * 100), )

print("Top 5 ITS2 types by Location:")
print(top5_by_location)

top5_by_location_morph <- data_long %>%
  group_by(Location, Morph, ITS2_type) %>%
  summarize(TotalAbundance = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Location, desc(TotalAbundance)) %>%
  group_by(Morph, Location) %>%
  slice_head(n = 5)
top5_by_location_morph <- top5_by_location_morph %>% group_by(Morph, Location) %>% mutate(percent = round(TotalAbundance/sum(TotalAbundance) * 100), )


library(kableExtra) #a table lib, for example
kable(top5_all, format = 'html') %>% # or 'simple'
  kable_classic_2(full_width = F)

kable(list(top5_all, top5_by_morph, top5_by_location), caption = c('Total', 'By Morph', 'By Island'), format = 'html') %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

top5_all %>%
  kable("html", align = 'clc', caption = "Top 5 ITS2 Profiles Overall") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), position = 'left')

top5_by_morph %>%
  kable("html", align = 'clc', caption = "Top 5 ITS2 Profiles by Morph") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), position = 'right')

top5_by_location %>%
  kable("html", align = 'clc', caption = 'Top 5 ITS2 Profiles by Location') %>%
  kable_styling(full_width = F, position = "right")


t1 <- top5_all %>%
  kable("html", align = 'clc', caption = "Top 5 ITS2 Profiles Overall") %>%
  kable_styling(full_width = F, position = "float_left")

t2 <- top5_by_morph %>%
  kable("html", align = 'clc', caption = "Top 5 ITS2 Profiles by Morph") %>%
  kable_styling(full_width = F, position = "float_right")

t3 <- top5_by_location %>%
  kable("html", align = 'clc', caption = 'Top 5 ITS2 Profiles by Location') %>%
  kable_styling(full_width = F, position = "right")

t1 <- reactable(top5_all, columns = list(
 # Morph = colDef(name = 'Morphotype'), 
  ITS2_type = colDef(name = 'ITS Profile'), 
  TotalAbundance = colDef(name = 'Total Abundance'), 
  percent = colDef(name = 'Percent of Total')),
  bordered = TRUE, striped = TRUE, highlight = TRUE)


t2 <- reactable(top5_by_location, columns = list(
  # Morph = colDef(name = 'Morphotype'), 
  ITS2_type = colDef(name = 'ITS Profile'), 
  TotalAbundance = colDef(name = 'Total Abundance'), 
  percent = colDef(name = 'Percent of Total')
  ),
  bordered = TRUE, striped = TRUE, highlight = TRUE)

t3 <- reactable(top5_by_morph, columns = list(
  Morph = colDef(name = 'Morphotype'), 
  ITS2_type = colDef(name = 'ITS Profile'), 
  TotalAbundance = colDef(name = 'Total Abundance'), 
  percent = colDef(name = 'Percent of Total')
  ),
  bordered = TRUE, striped = TRUE, highlight = TRUE,
  defaultPageSize = 10)

table <- reactable(top5_by_location_morph, columns = list(
  Morph = colDef(name = 'Morphotype'), 
  ITS2_type = colDef(name = 'ITS Profile'), 
  TotalAbundance = colDef(name = 'Total Abundance'), 
  percent = colDef(name = 'Percent of Total')
),
bordered = TRUE, striped = TRUE, highlight = TRUE,
defaultPageSize = 10)


## top 5 for each morph at each location
html <- "table.html"
library(htmlwidgets)
saveWidget(table, html)
library(webshot)
webshot(html, "rtableSnapshot.png")

# top 5 overall
html <- "t1.html"
library(htmlwidgets)
saveWidget(t1, html)
library(webshot)
webshot(html, "t1.png")

# top 5 by location
html <- "t2.html"
library(htmlwidgets)
saveWidget(t2, html)
library(webshot)
webshot(html, "t2.png")

# top 5 by morph 
html <- "t3.html"
library(htmlwidgets)
saveWidget(t3, html)
library(webshot)
webshot(html, "t3.png")
