# Installer les packages nécessaires si ce n’est pas encore fait
```
install.packages(c("tidyverse", "readxl", "sf", "rnaturalearth", "rnaturalearthdata", "ggspatial", "ggrepel", "viridis"))
```
# Charger les bibliothèques
```
library(tidyverse)
library(readxl)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(ggrepel)
library(viridis)
```
# Lire les données
```
data <- read_excel("Data_lobster.xlsx")
```
# Convertir en objet spatial
```
data_sf <- st_as_sf(data, coords = c("Long", "Lat"), crs = 4326)
```
# Charger les fonds de carte naturels
```
land <- ne_countries(scale = "medium", returnclass = "sf")
ocean <- ne_download(scale = 50, type = "ocean", category = "physical", returnclass = "sf")
```
# Création de la carte
# Créer la carte
```
ggplot() +
  # Fond océan
  geom_sf(data = ocean, fill = "lightblue", color = NA) +
  # Terre
  geom_sf(data = land, fill = "cornsilk", color = "gray60", size = 0.3) +
  # Cercle vide (juste le contour noir)
  geom_sf(data = data_sf, aes(size = `Nb`), 
          shape = 21, fill = NA, color = "black", stroke = 0.8) +
  # Noms des populations
  geom_text_repel(data = data_sf,
                  aes(geometry = geometry, label = Population_ID),
                  stat = "sf_coordinates",
                  size = 4,
                  fontface = "bold",
                  box.padding = 0.3,
                  max.overlaps = Inf,
                  color = "black") +
  # Ajustement des limites de carte
  coord_sf(xlim = range(st_coordinates(data_sf)[,1]) + c(-1, 1),
           ylim = range(st_coordinates(data_sf)[,2]) + c(-1, 1),
           expand = FALSE) +
  # Légende
  scale_size_continuous(name = "Number of samples", range = c(3, 10)) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "lightblue"),
        legend.position = "right") +
  labs(title = "Lobster population map",
       x = "Longitude", y = "Latitude")
```
