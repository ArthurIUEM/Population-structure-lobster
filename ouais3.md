# Lire le fichier avec le bon séparateur (essaie "\t" ou " " selon ton fichier)
env_all <- read.table("Filtered_ACP_Lobster_with_Lat_Env3.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# Affiche les noms de colonnes pour identifier celles utiles
colnames(env_all)

env_vars <- c("Bathymetry", "Max.temperature", "Mean.temperature", "Min.temperature",
              "Max.dissolution", "Mean.dissolution", "Min.dissolution",
              "Max.salinity", "Mean.salinity", "Min.salinity")

env_data <- env_all[, env_vars]

env_pca <- prcomp(env_data, scale. = TRUE)

# Afficher les résultats
summary(env_pca)

# Graphe des individus
biplot(env_pca)

library(factoextra)

# Plot des individus dans l’espace environnemental
fviz_pca_ind(env_pca, 
             geom.ind = "point", 
             col.ind = "black", 
             addEllipses = TRUE, 
             repel = TRUE)

library(ggplot2)
library(factoextra)

# Lire le fichier complet
env_all <- read.table("Filtered_ACP_Lobster_with_Lat_Env3.txt", 
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE)

# Vérifier les noms de colonnes
colnames(env_all)

# Extraire uniquement les colonnes environnementales (adapter les noms si nécessaire)
env_vars <- c("Bathymetry", "Max.temperature", "Mean.temperature", "Min.temperature",
              "Max.dissolution", "Mean.dissolution", "Min.dissolution",
              "Max.salinity", "Mean.salinity", "Min.salinity")
env_data <- env_all[, env_vars]

# ACP sur les variables environnementales
env_pca <- prcomp(env_data, scale. = TRUE)

# Extraire les coordonnées individuelles sur les deux premiers axes
pca_ind <- as.data.frame(env_pca$x[, 1:4])

# Ajouter la latitude depuis les données initiales
pca_ind$Latitude <- env_all$Latitude  # adapte le nom si nécessaire

ggplot(pca_ind, aes(x = PC1, y = PC2, color = Latitude)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +  # palette continue jolie
  labs(title = "ACP environnementale des individus",
       x = "PC1", y = "PC2", color = "Latitude") +
  theme_minimal()

ggplot(pca_ind, aes(x = PC2, y = PC3, color = Latitude)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(title = "ACP environnementale : PC2 vs PC3",
       x = "PC2", y = "PC3", color = "Latitude") +
  theme_minimal()

ggplot(pca_ind, aes(x = PC3, y = PC4, color = Latitude)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(title = "ACP environnementale : PC3 vs PC4",
       x = "PC3", y = "PC4", color = "Latitude") +
  theme_minimal()


df_plot <- data.frame(
  PC1_env = env_pca$x[, 1],
  PC3_gen = env_all$PC3,      # Remplace par le vrai nom si nécessaire
  Latitude = env_all$Latitude # Pour la couleur éventuellement
)

ggplot(df_plot, aes(x = PC1_env, y = PC3_gen, color = Latitude)) +
  geom_point(size = 2) +
  scale_color_viridis_c() +
  labs(title = "PC1 environnementale vs PC3 génétique",
       x = "PC1 environnementale",
       y = "PC3 génétique",
       color = "Latitude") +
  theme_minimal()
