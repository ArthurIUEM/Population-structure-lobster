# Chargement des librairies
library(ggplot2)
library(RColorBrewer)
library(dplyr)

# 1. Lecture des fichiers
geno_pca <- read.table("Lobster_top50K_PCA.eigenvec", header = FALSE)
colnames(geno_pca)[1:2] <- c("FID", "IID")

env <- read.delim("Filtered_ACP_Lobster_with_Lat_Env3.txt", header = TRUE)
zones <- read.delim("UMAP_zones_latitude.tsv", header = TRUE)

# 2. Fusion des fichiers par FID
df <- merge(geno_pca, env, by = "FID")
df <- merge(df, zones, by = "FID")

# 3. Scaling des PC environnementales (colonnes V3 à V6 = PC1 à PC4 environnementales)
for (i in 3:6) {
  df[[paste0("V", i)]] <- scale(df[[paste0("V", i)]])
}

# 4. Fonction pour tracer et sauvegarder les figures
plot_env_vs_geno_PC <- function(df, pc_index) {
  p <- ggplot(df, aes_string(x = paste0("PC", pc_index), 
                             y = paste0("V", pc_index + 2), 
                             color = "ZONE")) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      x = paste0("Environmental PC", pc_index, " (scaled)"),
      y = paste0("Genetic PC", pc_index),
      title = paste0("Environmental PC", pc_index, " vs Genetic PC", pc_index, " (top 50K)")
    ) +
    theme_minimal(base_size = 14) +
    scale_color_brewer(palette = "Dark2")
  
  # Sauvegarde
  ggsave(filename = paste0("Env_vs_Genetic_PC", pc_index, ".png"), plot = p, width = 8, height = 6)
  
  return(p)
}

# 5. Génération et sauvegarde des graphiques pour PC1 à PC4
plot_env_vs_geno_PC(df, 1)
plot_env_vs_geno_PC(df, 2)
plot_env_vs_geno_PC(df, 3)
plot_env_vs_geno_PC(df, 4)

