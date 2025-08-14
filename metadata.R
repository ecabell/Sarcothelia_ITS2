library(readr)
library(dplyr)
library(tidyr)
library(stringr)

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
metadata <- its2_data %>%
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
metadata <- metadata %>%
  filter(Morph != "Unknown", Location != "Unknown")

metadata <- metadata %>% subset(Sample != "HBL4")

colnames(metadata)[1] <- "Sample_Name"
