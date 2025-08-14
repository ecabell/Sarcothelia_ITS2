# Install ggmsa if you haven't already
# install.packages("ggmsa")

library(ggmsa)

# Assuming your aligned FASTA file is named "aligned_sequences.fasta"
# You might need to adjust the file path depending on your setup
aligned_fasta <- "~/Documents/Zoox PCR/20250603T195144_cabeler1/between_profile_distances/A/clade_A_seqs.aligned.fasta"

# Create the heatmap
ggmsa(aligned_fasta, color = "Chemistry_AA")  #  Customize color as needed
#A visual exploration tool for multiple sequence alignment and associated data