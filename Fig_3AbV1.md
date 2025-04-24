library(tidyverse)
library(readxl)
library(ggplot2)
library(sf)
library(ggrepel)
library(patchwork)

# Charger les données
umap_data <- read_tsv("UMAP_zones_latitude.tsv")
merged_data <- read_excel("Merged_Lobster_Data.xlsx") %>%
  rename(FID = `Sample ID 2`)

# Fusionner les données
full_data <- left_join(umap_data, merged_data, by = "FID")

# Regrouper et compter par site
localisation_summary <- full_data %>%
  group_by(LATITUDE, Long.x, ZONE) %>%
  summarise(n = n(), .groups = "drop")

# Définir un thème personnalisé avec des tailles de texte plus grandes
custom_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 21, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 17),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  )

# Fonction pour créer une carte stylisée
base_map <- function(data, zone_label) {
  ggplot(data, aes(x = Long.x, y = LATITUDE)) +
    borders("world", fill = "gray90", colour = "gray50") +
    geom_point(aes(size = n), color = ifelse(zone_label == "Nord", "blue", "red"), alpha = 0.7) +
    scale_size_continuous(name = "Nombre d'individus") +
    coord_cartesian(xlim = c(-66, -51), ylim = c(41, 52)) +
    labs(
      title = paste("Carte des échantillons - Zone", zone_label),
      x = "Longitude", y = "Latitude"
    ) +
    custom_theme
}

# Créer les cartes Nord et Sud
map_nord <- base_map(filter(localisation_summary, ZONE == "Nord"), "Nord")
map_sud  <- base_map(filter(localisation_summary, ZONE == "Sud"), "Sud")

# Afficher côte à côte
map_nord + map_sud
