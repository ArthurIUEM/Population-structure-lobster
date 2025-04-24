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
# Définir la palette de couleurs
```
color_palette <- c("red", "orange", "yellow", "green", "blue")
```
# Faire le graphique
```
ggplot(df, aes(x = UMAP1, y = UMAP2, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(
    title = "UMAP1 vs UMAP2",
    color = "Latitude",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 21, face = "bold"),
    axis.title = element_text(size = 17),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  )
```
