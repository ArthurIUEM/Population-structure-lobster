library(tidyverse)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(here)

# Chemins d'accès
env_file <- "Filtered_ACP_Lobster_with_Lat_Env3.txt"
zone_file <- "UMAP_zones_latitude.tsv"
pop_file <- "Merged_Lobster_Data.xlsx"

# Charger les données
df <- read_delim(env_file, delim = "\t")
zones <- read_tsv(zone_file) %>% select(FID, ZONE)
pop_data <- read_excel(pop_file) %>% select(FID, `Population ID`)

# Fusionner
df <- df %>%
  left_join(zones, by = "FID") %>%
  left_join(pop_data, by = "FID") %>%
  filter(!is.na(ZONE), !is.na(`Population ID`)) %>%
  mutate(ZONE = case_when(
    ZONE == "Nord" ~ "North",
    ZONE == "Sud"  ~ "South",
    TRUE ~ ZONE
  ))

# Ajouter les régions
df <- df %>%
  mutate(Region = case_when(
    `Population ID` %in% c("NUSPS", "NUSBB", "NUSWR") ~ "Prince Edward Island",
    `Population ID` %in% c("NSBC", "NSBDL", "NSF", "NSH", "NSLP") ~ "South-East Nova Scotia",
    `Population ID` %in% c("HB", "PS", "CC", "CBSPB", "CBTB") ~ "Newfoundland",
    `Population ID` %in% c("BOFSMC", "NSC") ~ "Bay of Fundy",
    `Population ID` %in% c("GB", "NC", "NSNWP") ~ "Gulf of Maine",
    TRUE ~ NA_character_
  ))

# Identifier les colonnes environnementales
env_cols <- df %>%
  select(-FID, -ZONE, -Region, -`Population ID`, -starts_with("PC"), -N_SITES, -Longitude) %>%
  colnames()

# Moyennes par population
df_env_grouped <- df %>%
  group_by(`Population ID`, Region) %>%
  summarise(across(all_of(env_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

# Ordre personnalisé des populations par région
df_env_grouped <- df_env_grouped %>%
  arrange(factor(Region, levels = c("Prince Edward Island", "South-East Nova Scotia", 
                                    "Newfoundland", "Bay of Fundy", "Gulf of Maine"))) %>%
  column_to_rownames("Population ID")

# Centrage-réduction
df_scaled <- scale(df_env_grouped[, env_cols])

# Annotation Nord/Sud
zone_info <- df %>%
  select(`Population ID`, ZONE) %>%
  distinct() %>%
  group_by(`Population ID`) %>%
  slice(1) %>%  # Prend la première occurrence
  ungroup() %>%
  column_to_rownames("Population ID")

# Gaps entre groupes pour créer effet de regroupement visuel
region_vector <- df_env_grouped$Region
region_labels <- unique(region_vector)
gaps_row <- cumsum(rle(region_vector)$lengths)

# Couleurs Nord/Sud
zone_colors <- list(ZONE = c("North" = "#1f78b4", "South" = "#e31a1c"))

# Palette de couleurs pour la heatmap
palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100)

# Heatmap
pheatmap(df_scaled,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         annotation_row = zone_info,
         annotation_colors = zone_colors,
         gaps_row = gaps_row,
         labels_row = rownames(df_scaled),
         color = palette,
         fontsize = 10,
         angle_col = 45,
         border_color = NA,
         main = "Influence of environmental variables by population and region")
