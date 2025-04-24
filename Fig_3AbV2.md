library(tidyverse)
library(readxl)
library(ggplot2)
library(sf)
library(ggrepel)

# Fichier contenant ZONE et latitude
umap_data <- read_tsv("UMAP_zones_latitude.tsv")

# Fichier Excel contenant la longitude
merged_data <- read_excel("Merged_Lobster_Data.xlsx") %>%
  rename(FID = `Sample ID 2`)

# Jointure par FID
full_data <- left_join(umap_data, merged_data, by = "FID")

# Regrouper par latitude, longitude et zone, puis compter les individus
localisation_summary <- full_data %>%
  group_by(LATITUDE, Long.x, ZONE) %>%
  summarise(n = n(), .groups = "drop")

# Thème personnalisé avec tailles agrandies
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# Créer la carte
ggplot(localisation_summary, aes(x = Long.x, y = LATITUDE)) +
  borders("world", fill = "gray90", colour = "gray50") +
  geom_point(aes(color = ZONE, size = n), alpha = 0.7) +
  scale_size_continuous(name = "Nombre d'individus") +
  scale_color_manual(
    name = "Zone",
    values = c("Nord" = "blue", "Sud" = "red")
  ) +
  coord_cartesian(xlim = c(-66, -51), ylim = c(41, 52)) +
  labs(
    title = "Carte des échantillons de homards - Atlantique canadien",
    x = "Longitude", y = "Latitude"
  ) +
  custom_theme
