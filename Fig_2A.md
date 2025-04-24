# Charger les packages
```
library(ggplot2)
library(patchwork)
```
# Charger les donnees
```
df <- read.table("UMAP_with_clusters_trimmed.txt", header = TRUE, sep = "\t")
df$LATITUDE <- as.numeric(as.character(df$LATITUDE))
```
# Faire la palette de coueleurs
```
color_palette <- c("red", "orange", "yellow", "green", "blue")
```
# Définir le thème commun
```
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
```
# Fixer les limites pour harmoniser les axes (facultatif mais utile)
```
xlims <- range(c(df$PC1, df$PC2, df$PC3, df$PC4), na.rm = TRUE)
ylims <- xlims  # pour un rapport 1:1
```
# Faire les graphiques
```
p1 <- ggplot(df, aes(x = PC1, y = PC2, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "PC1 vs PC2", color = "Latitude", x = "PC1", y = "PC2") +
  custom_theme +
  coord_fixed() +
  xlim(xlims) + ylim(ylims)
```
```
p2 <- ggplot(df, aes(x = PC2, y = PC3, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "PC2 vs PC3", color = "Latitude", x = "PC2", y = "PC3") +
  custom_theme +
  coord_fixed() +
  xlim(xlims) + ylim(ylims)
```
```
p3 <- ggplot(df, aes(x = PC3, y = PC4, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "PC3 vs PC4", color = "Latitude", x = "PC3", y = "PC4") +
  custom_theme +
  coord_fixed() +
  xlim(xlims) + ylim(ylims)
```
# Affichage en panel
```
panel <- (p1 | p2 | p3)
print(panel)
```
