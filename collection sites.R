# Load libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(ggrepel)

# Prepare your data
coords <- data.frame(
  Location = c("North Shore (Blue)", "Kailua (Brown)",
               "Koloa (Blue)", "Anini (Brown)",
               "South Kona (Blue)"),
  Lat = c(21.62982, 21.396089, 21.88431, 22.22545904, 19.26685),
  Lon = c(-158.07432, -157.721375, -159.48763, -159.4387606, -155.89490),
  Type = c("Blue","Brown","Blue","Brown","Blue")
)
# Add sample size column
coords$n <- c("n = 7", "n = 7", "n = 7", "n = 7", "n = 5")

# Convert to spatial
pts <- st_as_sf(coords, coords = c("Lon","Lat"), crs = 4326)

# Get base map
hi <- ne_states(country = "United States of America", returnclass = "sf") %>%
  filter(name == "Hawaii")

# Plot
p_no_n <- ggplot() +
  geom_sf(data = hi, fill = "#f4e7d3", color = "black", alpha = 0.7) +
  geom_sf(data = pts, aes(color = Type, shape = Type), size = 3) +
  geom_text_repel(
    data = coords,
    aes(x = Lon, y = Lat, label = Location),
    size = 4,
    nudge_x = c(0.1, 0.1, 0.1, 0.1, -0.05),  # right for 4, left for 1
    nudge_y = c(0.2, 0, -0.05, 0.1, 0)     # slight up/down nudges
  ) +
  scale_color_manual(values = c("Blue" = "blue", "Brown" = "brown")) +
  scale_shape_manual(values = c("Blue" = 16, "Brown" = 17)) +
  coord_sf(xlim = c(-160.5, -154.5), ylim = c(18.5, 23), expand = FALSE) +
  labs(title = expression("Collection Sites of " * italic("S. edmondsoni") * " Samples"),
       color = "Morphotype",
       shape = "Morphotype", 
       x = 'Longitude', 
       y = 'Latitude') +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#e0f7fa", color = NA),
    plot.background = element_rect(fill = "#e0f7fa", color = NA),  # <- this makes outer background white
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10), 
    legend.position = "right",
    legend.key = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.6)
  )

# Plot
p <- ggplot() +
  geom_sf(data = hi, fill = "#f4e7d3", color = "black", alpha = 0.7) +
  geom_sf(data = pts, aes(color = Type, shape = Type), size = 3) +
  
  # Main label
  geom_text_repel(
    data = coords,
    aes(x = Lon, y = Lat, label = Location),
    size = 4,
    nudge_x = c(0.1, 0.1, 0.1, 0.1, -0.05),
    nudge_y = c(0.2, 0, -0.05, 0.1, 0)
  ) +
  
  # Sample size label (smaller, just beneath)
  geom_text_repel(
    data = coords,
    aes(x = Lon, y = Lat, label = n),
    size = 3,
    nudge_x = c(0.15, 0.35, 0.15, 0.14, -.6),
    nudge_y = c(0.3, 0.13, -0.25, 0.24, -0.15), 
    segment.color = NA
  ) +
  
  scale_color_manual(values = c("Blue" = "blue", "Brown" = "brown")) +
  scale_shape_manual(values = c("Blue" = 16, "Brown" = 17)) +
  coord_sf(xlim = c(-160.5, -154.5), ylim = c(18.5, 23), expand = FALSE) +
  labs(
    title = expression("Collection Sites of " * italic("S. edmondsoni") * " Samples"),
    color = "Morphotype",
    shape = "Morphotype", 
    x = 'Longitude', 
    y = 'Latitude'
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "#e0f7fa", color = NA),
    plot.background = element_rect(fill = "#e0f7fa", color = NA),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.margin = margin(t = 10, r = 25, b = 10, l = 10), 
    legend.position = "right",
    legend.key = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.6)
  )

# Save figure
ggsave("hawaii_sample_sites.png", p, width = 8, height = 6, dpi = 500)
ggsave("hawaii_sample_sites_no_n.png", p_no_n, width = 8, height = 6, dpi = 500)
