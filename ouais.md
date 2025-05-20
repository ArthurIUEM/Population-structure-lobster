### read in your gds file
gdsfmt::showfile.gds(closeall=TRUE)
filename.gds <- "Lobster_top10K_fst.gds"
gdsin = SeqArray::seqOpen(filename.gds, readonly = TRUE)
### get it ready for rda
### (only run once)
dos.all = SeqVarTools::altDosage(gdsin)
SeqArray::seqResetFilter(gdsin)
### set up different genomic input datasets, with mean impute of missing since rda cannot have missing data
impute_mean <- function(x) {
  mean_value <- mean(x, na.rm = TRUE)
  ifelse(is.na(x), mean_value, x)
}
library(dplyr)
dos.all.input = SeqVarTools::altDosage(gdsin) %>%
  as.data.frame() %>%
  mutate(across(everything(), impute_mean))
rownames(dos.all.input) <- SeqArray::seqGetData(gdsin, "sample.id")

## Visualisation 

### Packages nécessaires
library(tidyverse)
library(vegan)

### Lire les variables environnementales
cols_env <- c(
  "FID", "Bathymetry", 
  "Max.temperature", "Mean.temperature", "Min.temperature",
  "Max.dissolution", "Mean.dissolution", "Min.dissolution",
  "Max.salinity", "Mean.salinity", "Min.salinity"
)

env <- read.delim("Filtered_ACP_Lobster_with_Lat_Env3.txt", header = TRUE)

### Garder seulement les colonnes d'intérêt
env_sub <- env[, cols_env]

### Vérifier l'ordre des individus (les FID doivent être les rownames de dos.all.input)
env_ordered <- env_sub[match(rownames(dos.all.input), env_sub$FID), ]
### RDA : on enlève la colonne "FID" pour ne garder que les variables explicatives
env_vars <- env_ordered %>% select(-FID)
### Lire les FID
env_sub <- env[, cols_env]   # contient FID
rownames(env_sub) <- as.character(env_sub$FID)

### Reordonner pour matcher dos.all.input
env_ordered <- env_sub[match(rownames(dos.all.input), rownames(env_sub)), ]


### Retirer la colonne FID car maintenant elle est dans les rownames
env_vars <- env_ordered %>% select(-FID)
### Ajouter la latitude séparément
env_lat <- env$Latitude[match(env_ordered$FID, env$FID)]

### Imputation par la moyenne
impute_mean <- function(x) {
  mean(x, na.rm = TRUE)
}

# Imputation
env_vars_imputed <- env_vars %>%
  mutate(across(everything(), ~ ifelse(is.na(.), impute_mean(.), .))) %>%
  mutate(Latitude = env_lat)

# Scaling (centrage-réduction)
env_vars_scaled <- as.data.frame(scale(env_vars_imputed))

### Synchroniser la matrice génotypique avec les mêmes lignes
geno_clean <- dos.all.input[match(env_ordered$FID, rownames(dos.all.input)), ]

### Trouver les lignes sans NA dans les deux jeux de données
complete_rows <- complete.cases(env_vars_scaled) & complete.cases(geno_clean)

### Appliquer ce filtre aux deux matrices
env_vars_scaled <- env_vars_scaled[complete_rows, ]
geno_clean <- geno_clean[complete_rows, ]

### Charger dplyr si ce n'est pas déjà fait
library(dplyr)

### Transformer en matrice si nécessaire
geno_clean_matrix <- as.matrix(geno_clean)

### Imputer la moyenne colonne par colonne
for (i in 1:ncol(geno_clean_matrix)) {
  missing_idx <- is.na(geno_clean_matrix[, i])
  if (any(missing_idx)) {
    col_mean <- mean(geno_clean_matrix[, i], na.rm = TRUE)
    geno_clean_matrix[missing_idx, i] <- col_mean
  }
}

### Vérification
sum(is.na(geno_clean_matrix))  # doit maintenant renvoyer 0

### Optionnel : renommer proprement
geno_clean_imputed <- geno_clean_matrix

rda_result <- rda(geno_clean_imputed ~ ., data = env_vars_scaled)


### Résumé et tests
summary(rda_result)
plot(rda_result, scaling = 2)

### Préparer les scores des individus
sc_si <- scores(rda_result, display = "sites", scaling = 2)

### Définir les axes à tracer
rdax <- 1
rday <- 2
perc <- round(100 * summary(rda_result)$cont$importance[2, c(rdax, rday)], 1)

### Synchroniser les lignes de données génotypiques et environnementales
complete_rows <- complete.cases(env_vars) & complete.cases(dos.all.input)  # Retirer les lignes avec des NA dans les deux
env_vars_clean <- env_vars[complete_rows, ]
geno_clean <- dos.all.input[complete_rows, ]

### Récupérer et synchroniser correctement la latitude
latitude_vector <- env$Latitude
names(latitude_vector) <- rownames(env_vars_clean)

### Appliquer la palette de couleurs pour la latitude
colors <- colorRampPalette(c("red", "green", "blue"))(100)
norm_vals <- as.numeric(cut(latitude_vector[rownames(sc_si)], breaks = 100))
point_cols <- colors[norm_vals]

### RDA et graphique
rda_result <- rda(geno_clean ~ ., data = env_vars_clean)

### Préparer les scores des individus
sc_si <- scores(rda_result, display = "sites", scaling = 2)
sc_arrows <- scores(rda_result, display = "bp", scaling = 2)
### Définir les axes à tracer
rdax <- 1
rday <- 2
eig_vals <- rda_result$CCA$eig
perc <- round(100 * eig_vals[c(rdax, rday)] / sum(eig_vals), 1)

### Extraire les scores des individus et des flèches
sc_si <- scores(rda_result, display = "sites", scaling = 2)
sc_bp <- scores(rda_result, display = "bp", scaling = 2)

### Choisir les deux axes
xaxis <- rdax
yaxis <- rday

### Définir les limites manuellement
x_range <- range(sc_si[, xaxis], sc_bp[, xaxis]) * 1.2
y_range <- range(sc_si[, yaxis], sc_bp[, yaxis]) * 1.2

### Triplot manuel
plot(rda_result,
     scaling = 2,
     type = "none",
     frame = FALSE,
     xlim = c(-8, 8),
     ylim = c(-10, 10),
     main = "Triplot RDA - scaling 2 (Top 10K)",
     xlab = paste0("RDA", rdax, " (", perc[1], "%)"),
     ylab = paste0("RDA", rday, " (", perc[2], "%)"),
     cex.axis = 1.5,
     cex.lab = 1.5)

### Points pour les individus colorés selon latitude
points(sc_si,
       cex = 2,
       col = point_cols,
       pch = 20)

### Flèches des variables explicatives
text(rda_result, 
     scaling = 2, 
     display = "bp", 
     col = "black", 
     cex = 1.2,          # taille du texte (ex: 1.2 plus petit que 2)
     arrow.mul = 5)    # longueur des flèches (0.5 = moitié de la taille par défaut)


### Scores des SNPs (variables dépendantes)
sc_species <- scores(rda_result, display = "species", scaling = 2)

### Exemple : loadings pour RDA1
loadings_rda1 <- sc_species[, 1]  # ou [, 2] pour RDA2

loadings_df <- data.frame(
  SNP = rownames(sc_species),
  Loading = loadings_rda1
)

threshold <- 2 * sd(loadings_df$Loading)
loadings_df$Significant <- abs(loadings_df$Loading) > threshold

library(ggplot2)

ggplot(loadings_df, aes(x = 1:nrow(loadings_df), y = Loading)) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = c(-threshold, threshold), linetype = "dashed", color = "blue") +
  labs(title = "RDA1 Loadings - Manhattan Plot (Top 10K)",
       x = "SNP Index",
       y = "Loading on RDA1") +
  theme_minimal()

### Exemple : loadings pour RDA1
loadings_rda2 <- sc_species[, 2]  # ou [, 2] pour RDA2

loadings_df <- data.frame(
  SNP = rownames(sc_species),
  Loading = loadings_rda2
)

threshold <- 2 * sd(loadings_df$Loading)
loadings_df$Significant <- abs(loadings_df$Loading) > threshold

library(ggplot2)

ggplot(loadings_df, aes(x = 1:nrow(loadings_df), y = Loading)) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = c(-threshold, threshold), linetype = "dashed", color = "blue") +
  labs(title = "RDA2 Loadings - Manhattan Plot (Top 10K)",
       x = "SNP Index",
       y = "Loading on RDA2") +
  theme_minimal()

library(qqman)

### Charger le fichier BIM
bim <- read.table("Lobster_top10K_fst.bim", header = FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")

### Vérifier que le nombre de SNPs est identique à celui des loadings
stopifnot(nrow(bim) == nrow(sc_species))

### RDA1
loadings_rda1 <- abs(sc_species[, 1])  # Valeurs absolues des loadings

manhattan_df_rda1 <- data.frame(
  CHR = bim$CHR,
  BP = bim$BP,
  SNP = paste0("SNP", 1:nrow(bim)),  # Noms génériques si besoin
  P = loadings_rda1                 # PAS de log10 ici
)

manhattan(
  manhattan_df_rda1,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  logp = FALSE,                     # IMPORTANT : ne pas appliquer de log10
  genomewideline = 2 * sd(loadings_rda1),  # seuil de sur/sous-expression
  suggestiveline = FALSE,
  ylim = c(0, 0.20),
  main = "Manhattan Plot - RDA1 (abs loadings) (Top 10K)",
  ylab = "abs(loading)"
)

### RDA2
loadings_rda2 <- abs(sc_species[, 2])

manhattan_df_rda2 <- data.frame(
  CHR = bim$CHR,
  BP = bim$BP,
  SNP = paste0("SNP", 1:nrow(bim)),
  P = loadings_rda2
)

manhattan(
  manhattan_df_rda2,
  chr = "CHR",
  bp = "BP",
  snp = "SNP",
  p = "P",
  logp = FALSE,
  genomewideline = 2 * sd(loadings_rda2),
  suggestiveline = FALSE,
  ylim = c(0, 0.1),
  main = "Manhattan Plot - RDA2 (abs loadings) (Top 10K)",
  ylab = "abs(loading)"
)


# Charger les packages nécessaires
library(ggplot2)
library(patchwork)

# 1. ACP environnementale
env_pca <- prcomp(env_vars_scaled, center = FALSE, scale. = FALSE)
env_scores <- as.data.frame(env_pca$x[, 1:4])
colnames(env_scores) <- paste0("PC", 1:4, "_env")
env_scores$FID <- rownames(env_vars_scaled)

# 2. Charger les scores PCA génétiques
geno_pca <- read.table("Lobster_top10K_PCA.eigenvec", header = FALSE)
colnames(geno_pca)[1:6] <- c("FID", "IID", paste0("PC", 1:4, "_geno"))

# 3. Charger les zones
zones <- read.table("UMAP_zones_latitude.tsv", header = TRUE, sep = "\t")
# Vérifier que la colonne ZONE et FID existent
stopifnot("FID" %in% colnames(zones), "ZONE" %in% colnames(zones))

# 4. Fusionner les données
merged_df <- merge(env_scores, geno_pca[, c("FID", paste0("PC", 1:4, "_geno"))], by = "FID")
merged_df <- merge(merged_df, zones[, c("FID", "ZONE")], by = "FID")

# 5. Graphiques comparant PC1-PC4 environnement vs génétique
plot_list <- lapply(1:4, function(i) {
  ggplot(merged_df, aes_string(x = paste0("PC", i, "_env"), y = paste0("PC", i, "_geno"), color = "ZONE")) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    scale_color_manual(values = c("Nord" = "#1b9e77", "Sud" = "#d95f02")) +
    labs(
      title = paste("PC", i, ": Environnement vs Génétique"),
      x = paste0("PC", i, " environnement"),
      y = paste0("PC", i, " génétique"),
      color = "Zone"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
})

# 6. Affichage en grille 2x2
(plot_list[[1]] | plot_list[[2]]) / (plot_list[[3]] | plot_list[[4]])


# Charger les packages nécessaires
library(tidyverse)

# 1. Lire les variables environnementales
cols_env <- c(
  "FID", "Bathymetry", 
  "Max.temperature", "Mean.temperature", "Min.temperature",
  "Max.dissolution", "Mean.dissolution", "Min.dissolution",
  "Max.salinity", "Mean.salinity", "Min.salinity"
)
env <- read.delim("Filtered_ACP_Lobster_with_Lat_Env3.txt", header = TRUE)

# 2. Sélectionner les variables d'intérêt + imputation des NA
env_vars <- env[, cols_env] %>%
  drop_na(FID) %>%
  column_to_rownames("FID") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))

# 3. Centrage-réduction
env_scaled <- scale(env_vars)

# 4. ACP
env_pca <- prcomp(env_scaled, center = FALSE, scale. = FALSE)

# 5. Résumé de l'ACP
summary(env_pca)

# 6. Screeplot (variance expliquée)
plot(env_pca, type = "l", main = "Screeplot - ACP Environnementale")

# 7. Visualisation des individus
library(ggplot2)

scores <- as.data.frame(env_pca$x)
scores$FID <- rownames(scores)

ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(
    title = "ACP des variables environnementales",
    x = paste0("PC1 (", round(100 * summary(env_pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(env_pca)$importance[2, 2], 1), "%)")
  ) +
  theme_minimal()







