# Charger les bibliothèques
```
library(dplyr)
```
# Lire le fichier UMAP
```
df <- read.table("UMAP_with_clusters_assigned.txt", header = TRUE, sep = "\t")
```
# Colonnes à trimmer
```
cols_to_trim <- c("PC1", "PC2", "PC3", "PC4", "UMAP1", "UMAP2")
```
# Appliquer le trimming 5% supérieurs pour chaque colonne
```
for (col in cols_to_trim) {
  upper <- quantile(df[[col]], 0.95, na.rm = TRUE)
  lower <- quantile(df[[col]], 0.05, na.rm = TRUE)
  df <- df %>% filter(.data[[col]] >= lower & .data[[col]] <= upper)
}
```
# Sauvegarder le fichier filtré
```
write.table(df, file = "UMAP_with_clusters_trimmed.txt", quote = FALSE, row.names = FALSE, sep = "\t")
```
# Charger les packages
```
library(ggplot2)
library(patchwork)
```
# Charger le fichier filtré après trimming
```
df <- read.table("UMAP_with_clusters_trimmed.txt", header = TRUE, sep = "\t")
```
# S'assurer que LATITUDE est bien numérique
```
df$LATITUDE <- as.numeric(as.character(df$LATITUDE))
```
# Définir la palette de couleurs : blue → green → yellow → orange → red
```
color_palette <- c("blue", "green", "yellow", "orange", "red")
```
# Créer les graphiques
```
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "PC1 vs PC2", color = "Latitude") +
  theme_minimal()
```
```
p2 <- ggplot(df, aes(x = PC2, y = PC3, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "PC2 vs PC3", color = "Latitude") +
  theme_minimal()
```
```
p3 <- ggplot(df, aes(x = PC3, y = PC4, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "PC3 vs PC4", color = "Latitude") +
  theme_minimal()
```
```
p4 <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "UMAP1 vs UMAP2", color = "Latitude") +
  theme_minimal()
```
# Affichage en panel 2x2
```
panel <- (p1 | p2) / (p3 | p4)
print(panel)
```
