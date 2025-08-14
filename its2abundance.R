library(dplyr)

# Load your data
file_path <- "~/Documents/Zoox PCR/20250603T195144_cabeler1/post_med_seqs/596_20250604T195602_DBV_20250605T195712.seqs.absolute.abund_and_meta.csv"
data <- read.csv(file_path, check.names = FALSE)

# Subset only the sequence abundance columns (assuming they start from column 34 onward)
abund_data <- data[, 34:ncol(data)]

# Ensure numeric
abund_data <- abund_data %>% mutate(across(everything(), as.numeric))

# Total abundance
total_abundance <- sum(abund_data, na.rm = TRUE)

# A3 columns and abundance
a3_cols <- grep("A3", colnames(abund_data), value = TRUE)
a3_abundance <- sum(abund_data[, a3_cols], na.rm = TRUE)
a3_percent <- (a3_abundance / total_abundance) * 100

# A-clade columns and abundance
a_cols <- grep("\\bA\\b|_A|^A", colnames(abund_data), value = TRUE)
a_abundance <- sum(abund_data[, a_cols], na.rm = TRUE)
a_percent <- (a_abundance / total_abundance) * 100

# Non-A3 A-clade columns and abundance
non_a3_a_cols <- setdiff(a_cols, a3_cols)
non_a3_a_abundance <- sum(abund_data[, non_a3_a_cols], na.rm = TRUE)
non_a3_a_percent <- (non_a3_a_abundance / total_abundance) * 100

# Top 5 non-A3 A-clade types
non_a3_a_sums <- colSums(abund_data[, non_a3_a_cols], na.rm = TRUE)
top5_non_a3_a <- sort(non_a3_a_sums, decreasing = TRUE)[1:5]
top5_non_a3_a_abundance <- sum(top5_non_a3_a)
top5_non_a3_a_percent <- (top5_non_a3_a_abundance / total_abundance) * 100

# Remaining non-A3 A-clade abundance
remaining_non_a3_a_abundance <- non_a3_a_abundance - top5_non_a3_a_abundance
remaining_non_a3_a_percent <- (remaining_non_a3_a_abundance / total_abundance) * 100

# Non-clade A abundance
non_a_abundance <- total_abundance - a_abundance
non_a_percent <- (non_a_abundance / total_abundance) * 100

# Output
cat("Total abundance:", total_abundance, "\n")
cat("A3 abundance (%):", round(a3_percent, 2), "\n")
cat("A clade abundance (%):", round(a_percent, 2), "\n")
cat("Non-A3 A clade abundance (%):", round(non_a3_a_percent, 2), "\n")
cat("Top 5 Non-A3 A clade types and percent:\n")
print(round((top5_non_a3_a / total_abundance) * 100, 2))
cat("Remaining Non-A3 A clade abundance (%):", round(remaining_non_a3_a_percent, 2), "\n")
cat("Non-clade A abundance (%):", round(non_a_percent, 2), "\n")

