# Charger les bibliothèques nécessaires
library(ggplot2)
library(dplyr)
library(gridExtra)
library(umap)

# Charger les données ACP et profondeur de séquençage
acp <- read.table("Lobster_PCA_chr3.eigenvec", header=FALSE)
depth <- read.table("Lobster1MB.idepth", header=TRUE)

colnames(acp) <- c("FID", "IID", paste0("PC", 1:(ncol(acp)-2)))
if (!"FID" %in% colnames(depth)) colnames(depth)[1] <- "FID"

acp <- acp %>% select(-IID) %>% merge(depth, by="FID")

for (i in 1:9) {
  upper_pc <- quantile(acp[[paste0("PC", i)]], 0.95)
  acp <- acp %>% filter(acp[[paste0("PC", i)]] <= upper_pc)
}

fam <- read.table("Lobster1MB.fam", header=FALSE)
if (ncol(fam) < 6) stop("Erreur : Le fichier .fam ne contient pas suffisamment de colonnes.")
colnames(fam)[1:7] <- c("FID", "IID", "PID", "MID", "SEX", "KO", "OK")
fam$LATITUDE <- fam$OK

acp <- merge(acp, fam[, c("FID", "LATITUDE")], by="FID", all.x=TRUE)

eigenval <- scan("Lobster_PCA.eigenval")
variance <- eigenval / sum(eigenval) * 100
variance_labels <- sprintf("%.1f%%", variance)

color_palette <- c("red", "orange", "yellow", "green", "blue")

# Définir un thème commun
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 21, face = "bold"),
    axis.title = element_text(size = 17),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13)
  )

# Créer la liste de graphes PCA
plot_list <- list()
for (i in 1:3) {
  x_col <- paste0("PC", i)
  y_col <- paste0("PC", i+1)
  
  min_x <- min(acp[[x_col]], na.rm = TRUE)
  max_x <- max(acp[[x_col]], na.rm = TRUE)
  min_y <- min(acp[[y_col]], na.rm = TRUE)
  max_y <- max(acp[[y_col]], na.rm = TRUE)
  min_axis <- min(min_x, min_y)
  max_axis <- max(max_x, max_y)
  
  p <- ggplot(acp, aes_string(x = x_col, y = y_col, color = "LATITUDE")) +
    geom_point(size = 2) +
    scale_color_gradientn(colors = color_palette) +
    coord_fixed() +
    xlim(min_axis, max_axis) + ylim(min_axis, max_axis) +
    custom_theme +
    labs(
      title = paste0(x_col, " vs ", y_col),
      x = paste0(x_col, " (", variance_labels[i], ")"),
      y = paste0(y_col, " (", variance_labels[i+1], ")"),
      color = "Latitude"
    )
  
  plot_list[[i]] <- p
}

# Calculer UMAP
umap_input <- acp %>% select(starts_with("PC"))
set.seed(123)
umap_result <- umap(umap_input)
acp$UMAP1 <- umap_result$layout[,1]
acp$UMAP2 <- umap_result$layout[,2]

# Graphe UMAP
p_umap <- ggplot(acp, aes(x = UMAP1, y = UMAP2, color = LATITUDE)) +
  geom_point(size = 2) +
  scale_color_gradientn(colors = color_palette) +
  custom_theme +
  labs(
    title = "UMAP1 vs UMAP2",
    x = "UMAP1",
    y = "UMAP2",
    color = "Latitude"
  )

plot_list[[4]] <- p_umap

# Affichage 2x2
do.call(grid.arrange, c(plot_list, ncol = 2))

# Définir les zones Nord/Sud selon UMAP2
acp <- acp %>%
  mutate(ZONE = ifelse(UMAP2 < 0, "Nord", "Sud"))

# Sélection et export
acp_export <- acp %>%
  select(FID, UMAP1, UMAP2, LATITUDE, ZONE)

write.table(acp_export, file = "UMAP_zones_latitude.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
