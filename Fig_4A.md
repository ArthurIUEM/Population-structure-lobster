# Sans coloration 
# 🔹 Chargement des packages
library(data.table)
library(vegan)

# 🔹 Paramètres
threshold <- 1e-5
pcs <- 1:4

# 🔹 Extraction des SNPs significatifs après filtrage sur CHROM
for (pc in pcs) {
  file <- paste0("PC", pc, "_assoc.PC", pc, ".glm.linear")
  
  data <- read.table(file, header = TRUE)
  colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTEDA1_FREQ",
                      "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
  
  # S'assurer que CHROM est bien un caractère
  data$CHROM <- as.character(data$CHROM)
  
  # 🔥 Filtrer : Retirer les lignes où CHROM == "24712526"
  data <- data[data$CHROM != "24712526", ]
  
  # 🔹 Extraction des SNPs sous le seuil
  outliers <- data[data$P < threshold, "ID"]
  
  # 🔹 Sauvegarde
  write.table(outliers, paste0("snps_PC", pc, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# 🔹 Chargement du fichier fam pour récupérer les latitudes
fam <- fread("Lobster1MB.fam", header = FALSE)

# Supposons que la latitude soit dans la 6e colonne (à adapter si nécessaire)
# Vérifie si besoin avec head(fam)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX","Ok", "Latitude", "Longitude", "NAFO Zone")

lat_data <- fam[, .(ID = IID, Latitude)]


# 🔹 Fonction pour exécuter la RDA
run_rda <- function(raw_prefix, env_data, pc_label) {
  message(paste0("\n🔎 Analyse RDA pour ", pc_label))
  
  geno <- fread(paste0(raw_prefix, ".raw"))
  
  # Nettoyage des colonnes inutiles
  geno <- geno[, !c("FID", "MAT", "PAT", "SEX", "PHENOTYPE"), with = FALSE]
  setnames(geno, "IID", "ID")
  setnames(env_data, "FID", "ID")
  
  # Identifier les individus communs
  common_ids <- intersect(geno$ID, env_data$ID)
  geno <- geno[ID %in% common_ids]
  env_data <- env_data[ID %in% common_ids]
  setkey(geno, ID)
  setkey(env_data, ID)
  
  # Supprimer les lignes avec NA dans env_data
  env_data <- na.omit(env_data)
  geno <- geno[ID %in% env_data$ID]
  
  geno_final <- geno[, !"ID"]
  env_final <- env_data[, !"ID"]
  
  # Conversion en numérique
  geno_final <- geno_final[, lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x))))]
  
  # 🔥 Supprimer les SNPs avec >5% de NA
  geno_final <- geno_final[, which(colMeans(is.na(geno_final)) < 0.05), with = FALSE]
  
  # 🔥 Imputation des valeurs manquantes restantes
  geno_final <- geno_final[, lapply(.SD, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))]
  
  # 🔹 Nettoyage de l'environnement
  env_final <- env_final[, which(sapply(env_final, is.numeric)), with = FALSE]
  env_final <- env_final[, which(colSums(is.na(env_final)) < nrow(env_final)), with = FALSE]
  env_final <- env_final[, which(sapply(env_final, sd, na.rm = TRUE) > 0), with = FALSE]
  
  # 🔹 Standardisation
  env_final <- as.data.table(scale(env_final))
  
  # Vérification finale
  if (any(is.na(geno_final)) || any(is.na(env_final))) stop("❌ Données incomplètes après nettoyage.")
  
  # 🔥 Lancement RDA
  message("🚀 Lancement de la RDA...")
  rda_model <- rda(geno_final ~ ., data = env_final)
  print(summary(rda_model))
  
  return(rda_model)
}

# 🔹 Chargement des données environnementales
env_data_full <- fread("Filtered_ACP_Lobster_with_Lat_Env2.txt")
cols_env <- c("FID", "Bathymetrie_moy", "Temperature_max", "Temperature_moy", "Temperature_min",
              "Dissolution_max", "Dissolution_moy", "Dissolution_min",
              "Salinity_max", "Salinity_moy", "Salinity_min")
env_data <- env_data_full[, ..cols_env]

# 🔹 Lancement des modèles RDA pour PC2, PC3, PC4
rda_PC2 <- run_rda("PC2_outliers", copy(env_data), "PC2")
rda_PC3 <- run_rda("PC3_outliers", copy(env_data), "PC3")
rda_PC4 <- run_rda("PC4_outliers", copy(env_data), "PC4")

# 🔹 Fonction pour afficher joliment les RDA
plot_rda_with_var <- function(rda_model, rda_label) {
  # Extraire la variance expliquée
  var_exp <- summary(rda_model)$constrained$importance[2, ]
  x_lab <- paste0("RDA1 (", round(var_exp[1] * 100, 1), "%)")
  y_lab <- paste0("RDA2 (", round(var_exp[2] * 100, 1), "%)")
  
  # 🔹 Triplot général
  plot(rda_model, type = "n", main = paste("RDA -", rda_label, "without sex chr"), xlab = x_lab, ylab = y_lab)
  points(rda_model, display = "sites", col = "grey30", pch = 16, cex = 0.8)
  text(rda_model, display = "bp", col = "blue", cex = 1)
  text(rda_model, display = "species", col = "red", cex = 0.6)
  
  # 🔹 Plot sites uniquement
  plot(rda_model, type = "n", main = paste("RDA - Indiv -", rda_label), xlab = x_lab, ylab = y_lab)
  points(rda_model, display = "sites", col = "black", pch = 16, cex = 0.8)
  
  # 🔹 Plot environnement uniquement
  plot(rda_model, type = "n", main = paste("RDA - Environnement -", rda_label), xlab = x_lab, ylab = y_lab)
  text(rda_model, display = "bp", col = "blue", cex = 1)
  
  # 🔹 Plot SNPs uniquement
  plot(rda_model, type = "n", main = paste("RDA - SNPs -", rda_label), xlab = x_lab, ylab = y_lab)
  text(rda_model, display = "species", col = "red", cex = 0.5)
}

# 🔹 Affichage des résultats
plot_rda_with_var(rda_PC2, "PC2")
plot_rda_with_var(rda_PC3, "PC3")
plot_rda_with_var(rda_PC4, "PC4")


# Colore par latitude
# ------------------- Chargement des packages -------------------
library(data.table)
library(vegan)

# ------------------- Paramètres -------------------
threshold <- 1e-5
pcs <- 1:4

# ------------------- Extraction des SNPs significatifs par PC -------------------
for (pc in pcs) {
  file <- paste0("PC", pc, "_assoc.PC", pc, ".glm.linear")
  data <- fread(file)
  
  # Renommer les colonnes
  setnames(data, c("CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTEDA1_FREQ", "A1_FREQ",
                   "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE"))
  
  # Retirer lignes où CHROM == 24712526
  data <- data[CHROM != 24712526]
  
  # Garder seulement SNPs significatifs
  outliers <- data[P < threshold, ID]
  
  # Sauvegarder
  fwrite(data.table(outliers), paste0("snps_PC", pc, ".txt"), col.names = FALSE)
}

# ------------------- Chargement des données environnementales -------------------
env_data_full <- fread("Filtered_ACP_Lobster_with_Lat_Env2.txt")

# Colonnes environnementales sélectionnées
cols_env <- c("FID", "Bathymetrie_moy", "Temperature_max", "Temperature_moy", "Temperature_min",
              "Dissolution_max", "Dissolution_moy", "Dissolution_min",
              "Salinity_max", "Salinity_moy", "Salinity_min")

env_data <- env_data_full[, ..cols_env]

# Récupération de la latitude
lat_data <- env_data_full[, .(ID = FID, LATITUDE)]

# ------------------- Fonction pour lancer la RDA -------------------
run_rda <- function(raw_prefix, env_data, pc_label) {
  message(paste0("\n🔎 Analyse RDA pour ", pc_label))
  
  # Charger génotypes
  geno <- fread(paste0(raw_prefix, ".raw"))
  geno <- geno[, !c("FID", "MAT", "PAT", "SEX", "PHENOTYPE"), with = FALSE]
  setnames(geno, "IID", "ID")
  setnames(env_data, "FID", "ID")
  
  # Gérer les individus communs
  common_ids <- intersect(geno$ID, env_data$ID)
  geno <- geno[ID %in% common_ids]
  env_data <- env_data[ID %in% common_ids]
  setkey(geno, ID)
  setkey(env_data, ID)
  
  # Nettoyage des NAs
  env_data <- na.omit(env_data)
  geno <- geno[ID %in% env_data$ID]
  
  geno_final <- geno[, !"ID"]
  env_final <- env_data[, !"ID"]
  
  # Convertir génotypes en numériques
  geno_final <- geno_final[, lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x))))]
  
  # Enlever SNPs avec >5% de NA
  geno_final <- geno_final[, which(colMeans(is.na(geno_final)) < 0.05), with = FALSE]
  
  # Imputation des NA
  geno_final <- geno_final[, lapply(.SD, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))]
  
  # Sélectionner seulement variables environnementales numériques valides
  env_final <- env_final[, which(sapply(env_final, is.numeric)), with = FALSE]
  env_final <- env_final[, which(colSums(is.na(env_final)) < nrow(env_final)), with = FALSE]
  env_final <- env_final[, which(sapply(env_final, sd, na.rm = TRUE) > 0), with = FALSE]
  
  # Standardiser environnement
  env_final <- as.data.table(scale(env_final))
  
  # ➡️ Corriger pour rownames
  geno_final_df <- as.data.frame(geno_final)
  rownames(geno_final_df) <- geno$ID
  
  env_final_df <- as.data.frame(env_final)
  rownames(env_final_df) <- env_data$ID
  
  # Vérifications
  if (any(is.na(geno_final_df)) || any(is.na(env_final_df))) stop("❌ Données incomplètes.")
  
  # Lancer la RDA
  message("🚀 Lancement de la RDA...")
  rda_model <- rda(geno_final_df ~ ., data = env_final_df)
  print(summary(rda_model))
  
  return(rda_model)
}

# ------------------- Fonction pour afficher la RDA -------------------
plot_rda_with_var <- function(rda_model, rda_label, lat_data) {
  message("🎨 Plot RDA - ", rda_label)
  
  # Variance expliquée
  var_exp <- summary(rda_model)$constrained$importance[2, ]
  x_lab <- paste0("RDA1 (", round(var_exp[1] * 100, 1), "%)")
  y_lab <- paste0("RDA2 (", round(var_exp[2] * 100, 1), "%)")
  
  # Scores des sites (individus)
  site_scores <- scores(rda_model, display = "sites")
  site_ids <- rownames(site_scores)
  
  # Garder latitudes correspondantes
  lat_data_use <- lat_data[ID %in% site_ids]
  lat_data_use <- lat_data_use[match(site_ids, lat_data_use$ID)]
  
  # Récupérer latitudes
  latitudes <- as.numeric(lat_data_use$LATITUDE)
  
  # Nettoyer
  valid_latitudes <- !is.na(latitudes) & is.finite(latitudes)
  latitudes <- latitudes[valid_latitudes]
  site_scores <- site_scores[valid_latitudes, , drop = FALSE]
  
  if (length(latitudes) == 0) stop("❌ Problème : aucune latitude valide après nettoyage.")
  
  # ➡️ Créer un gradient continu
  color_palette <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))(100)
  
  # ➡️ Normaliser latitudes entre 1 et 100
  lat_min <- min(latitudes)
  lat_max <- max(latitudes)
  color_indices <- round((latitudes - lat_min) / (lat_max - lat_min) * 99) + 1
  
  # ➡️ Plot Sites avec gradient
  plot(rda_model, type = "n", main = paste("RDA - Sites par Latitude -", rda_label), xlab = x_lab, ylab = y_lab)
  points(site_scores[,1], site_scores[,2], pch = 16, cex = 0.8, col = color_palette[color_indices])
  
  # ➡️ Ajouter une légende en dégradé
  legend_gradient <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"))(5)
  legend("topright", legend = round(seq(lat_min, lat_max, length.out = 5), 2), fill = legend_gradient, title = "Latitude", cex = 0.8)
  
  # ➡️ Plot Variables
  plot(rda_model, type = "n", main = paste("RDA - Variables Environnementales -", rda_label), xlab = x_lab, ylab = y_lab)
  text(rda_model, display = "bp", col = "blue", cex = 1)
  
  # ➡️ Plot SNPs
  plot(rda_model, type = "n", main = paste("RDA - SNPs -", rda_label), xlab = x_lab, ylab = y_lab)
  text(rda_model, display = "species", col = "red", cex = 0.5)
}


# ------------------- Lancer les RDA pour les PCs -------------------
rda_PC2 <- run_rda("PC2_outliers", copy(env_data), "PC2")
rda_PC3 <- run_rda("PC3_outliers", copy(env_data), "PC3")
rda_PC4 <- run_rda("PC4_outliers", copy(env_data), "PC4")

# ------------------- Plots -------------------
plot_rda_with_var(rda_PC2, "PC2", lat_data)
plot_rda_with_var(rda_PC3, "PC3", lat_data)
plot_rda_with_var(rda_PC4, "PC4", lat_data)

# Pour le groupe nord et sud
# 🔹 Chargement des packages
library(data.table)
library(vegan)

# 🔹 Paramètres
threshold <- 1e-5
pcs <- 1:4

# 🔹 Extraction des SNPs significatifs après filtrage sur CHROM
for (pc in pcs) {
  file <- paste0("PC", pc, "_assoc.PC", pc, ".glm.linear")
  
  data <- read.table(file, header = TRUE)
  colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTEDA1_FREQ",
                      "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
  
  data$CHROM <- as.character(data$CHROM)
  data <- data[data$CHROM != "24712526", ]
  
  outliers <- data[data$P < threshold, "ID"]
  
  write.table(outliers, paste0("snps_PC", pc, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# 🔹 Chargement du fichier fam pour les latitudes
fam <- fread("Lobster1MB.fam", header = FALSE)
colnames(fam) <- c("FID", "IID", "PAT", "MAT", "SEX","Ok", "Latitude", "Longitude", "NAFO_Zone")
lat_data <- fam[, .(ID = IID, Latitude)]

# 🔹 Chargement des données de zone (Nord/Sud)
zones <- fread("UMAP_zones_latitude.tsv")
colnames(zones)[1] <- "ID"

# 🔹 Chargement des données environnementales
env_data_full <- fread("Filtered_ACP_Lobster_with_Lat_Env2.txt")
cols_env <- c("FID", "Bathymetrie_moy", "Temperature_max", "Temperature_moy", "Temperature_min",
              "Dissolution_max", "Dissolution_moy", "Dissolution_min",
              "Salinity_max", "Salinity_moy", "Salinity_min")
env_data <- env_data_full[, ..cols_env]

# 🔹 Fusion des données environnementales et des zones
env_data <- merge(env_data, zones[, .(ID, ZONE)], by.x = "FID", by.y = "ID")

# 🔹 Fonction pour exécuter la RDA
run_rda <- function(raw_prefix, env_data, pc_label) {
  message(paste0("\n🔎 Analyse RDA pour ", pc_label))
  
  geno <- fread(paste0(raw_prefix, ".raw"))
  geno <- geno[, !c("FID", "MAT", "PAT", "SEX", "PHENOTYPE"), with = FALSE]
  setnames(geno, "IID", "ID")
  setnames(env_data, "FID", "ID")
  
  common_ids <- intersect(geno$ID, env_data$ID)
  geno <- geno[ID %in% common_ids]
  env_data <- env_data[ID %in% common_ids]
  setkey(geno, ID)
  setkey(env_data, ID)
  
  env_data <- na.omit(env_data)
  geno <- geno[ID %in% env_data$ID]
  
  geno_final <- geno[, !"ID"]
  env_final <- env_data[, !"ID"]
  
  geno_final <- geno_final[, lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x))))]
  geno_final <- geno_final[, which(colMeans(is.na(geno_final)) < 0.05), with = FALSE]
  geno_final <- geno_final[, lapply(.SD, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))]
  
  env_final <- env_final[, which(sapply(env_final, is.numeric)), with = FALSE]
  env_final <- env_final[, which(colSums(is.na(env_final)) < nrow(env_final)), with = FALSE]
  env_final <- env_final[, which(sapply(env_final, sd, na.rm = TRUE) > 0), with = FALSE]
  env_final <- as.data.table(as.data.frame(scale(env_final)))
  
  if (any(is.na(geno_final)) || any(is.na(env_final))) stop("❌ Données incomplètes après nettoyage.")
  
  message("🚀 Lancement de la RDA...")
  rda_model <- rda(geno_final ~ ., data = env_final)
  print(summary(rda_model))
  
  return(rda_model)
}

# 🔹 Fonction d'affichage avec choix d'axes
plot_rda_axes <- function(rda_model, rda_label, axes = c(1, 2)) {
  var_exp <- summary(rda_model)$constrained$importance[2, ]
  x_lab <- paste0("RDA", axes[1], " (", round(var_exp[axes[1]] * 100, 1), "%)")
  y_lab <- paste0("RDA", axes[2], " (", round(var_exp[axes[2]] * 100, 1), "%)")
  
  plot(rda_model, type = "n", choices = axes, main = paste("RDA -", rda_label, "axes", axes[1], "vs", axes[2]),
       xlab = x_lab, ylab = y_lab)
  points(rda_model, display = "sites", choices = axes, col = "grey30", pch = 16, cex = 0.8)
  text(rda_model, display = "bp", choices = axes, col = "blue", cex = 1)
  text(rda_model, display = "species", choices = axes, col = "red", cex = 0.6)
}

# 🔹 Lancement pour chaque PC et chaque zone
zones_list <- c("Nord", "Sud")
pcs_to_run <- 2:4

for (zone in zones_list) {
  for (pc in pcs_to_run) {
    env_sub <- env_data[ZONE == zone]
    rda_result <- run_rda(paste0("PC", pc, "_outliers"), copy(env_sub), paste0("PC", pc, " - ", zone))
    
    # Affichages multiples
    plot_rda_axes(rda_result, paste0("PC", pc, " - ", zone), axes = c(1, 2))
    plot_rda_axes(rda_result, paste0("PC", pc, " - ", zone), axes = c(2, 3))
    plot_rda_axes(rda_result, paste0("PC", pc, " - ", zone), axes = c(3, 4))
  }
}

