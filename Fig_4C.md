# Charger les bibliothèques nécessaires
library(dplyr)
library(readxl)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)

# 1. Lire les scores UMAP + Latitude + Zone (généré précédemment)
umap_data <- read.table("UMAP_zones_latitude.tsv", header = TRUE, sep = "\t")

# 2. Lire le fichier Excel avec les noms des populations
# Attention : adapte le nom de l’onglet si ce n’est pas le premier
metadata <- read_excel("Merged_Lobster_Data.xlsx", sheet = 1)

# Vérifie les noms de colonnes
colnames(metadata)

# 3. Renommer la colonne contenant les FID si nécessaire
metadata <- metadata %>%
  rename(FID = `Sample ID 2`)  # adapte ce nom si différent

# Supposons que la colonne avec le nom de la population s'appelle "Population"
# Et que "Longitude" et "Latitude" sont aussi bien nommées dans le fichier Excel

# 1. Fusionner les données
merged_data <- left_join(umap_data, metadata, by = "FID")

# 2. Filtrer les données avec coordonnées valides et nom de population non manquant
merged_data <- merged_data %>%
  filter(!is.na(Long.x), !is.na(Lat.x), !is.na(`Population ID`))

# 3. Créer les datasets Nord/Sud
nord_data <- merged_data %>% filter(ZONE == "Nord")
sud_data  <- merged_data %>% filter(ZONE == "Sud")

# 4. Pour chaque groupe, créer un seul point par population avec moyenne des coordonnées
nord_pop <- nord_data %>%
  group_by(`Population ID`) %>%
  summarise(
    Latitude = mean(Lat.x),
    Longitude = mean(Long.x)
  )

sud_pop <- sud_data %>%
  group_by(`Population ID`) %>%
  summarise(
    Latitude = mean(Lat.x),
    Longitude = mean(Long.x)
  )

# 5. Charger carte
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

# 6. Fonction pour faire une carte avec une seule fois chaque population
plot_group_map <- function(data, title) {
  ggplot() +
    geom_sf(data = world, fill = "gray90", color = "gray40") +
    geom_point(data = data, aes(x = Longitude, y = Latitude), color = "blue", size = 2.5) +
    geom_text_repel(data = data, aes(x = Longitude, y = Latitude, label = `Population ID`), size = 3.2) +
    coord_sf(xlim = c(-70, -50), ylim = c(40, 52), expand = FALSE) +
    theme_minimal() +
    labs(title = title, x = "Longitude", y = "Latitude")
}

# 7. Générer les cartes
map_nord <- plot_group_map(nord_pop, "North Group")
map_sud  <- plot_group_map(sud_pop, "South Group")

# 8. Afficher côte à côte
library(gridExtra)
grid.arrange(map_nord, map_sud, ncol = 2)
