# Code de base pour la base de donnee complete

# Chargement des fichiers d'association et extraction des SNPs significatifs
threshold <- 1e-5
pcs <- 1:4

for (pc in pcs) {
  file <- paste0("PC", pc, "_assoc.PC", pc, ".glm.linear")
  data <- read.table(file, header = TRUE)
  colnames(data) <- c("CHROM", "POS", "ID", "REF", "ALT", "PROVISIONAL_REF", "A1", "OMITTEDA1_FREQ", 
                      "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
  outliers <- data[data$P < threshold, "ID"]
  write.table(outliers, paste0("snps_PC", pc, ".txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

library(data.table)
library(vegan)

run_rda <- function(raw_prefix, env_data, pc_label) {
  message(paste0("\n🔎 Analyse RDA pour ", pc_label))
  
  geno <- fread(paste0(raw_prefix, ".raw"))
  geno <- geno[, !c("FID", "MAT", "PAT", "SEX", "PHENOTYPE"), with = FALSE]
  setnames(geno, "IID", "ID")
  setnames(env_data, "FID", "ID")
  
  # Identifier et réaligner les individus communs
  common_ids <- intersect(geno$ID, env_data$ID)
  geno <- geno[ID %in% common_ids]
  env_data <- env_data[ID %in% common_ids]
  setkey(geno, ID)
  setkey(env_data, ID)
  
  # Nettoyage des données
  env_data <- na.omit(env_data)
  geno <- geno[ID %in% env_data$ID]
  
  geno_final <- geno[, !"ID"]
  env_final <- env_data[, !"ID"]
  
  # Conversion des génotypes en numériques si besoin
  geno_final <- geno_final[, lapply(.SD, function(x) suppressWarnings(as.numeric(as.character(x))))]
  
  # Suppression des SNPs avec >5% de NA
  geno_final <- geno_final[, which(colMeans(is.na(geno_final)) < 0.05), with = FALSE]
  
  # Imputation des valeurs manquantes restantes
  geno_final <- geno_final[, lapply(.SD, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))]
  
  # Garder uniquement les colonnes numériques et non constantes dans l'environnement
  env_final <- env_final[, which(sapply(env_final, is.numeric)), with = FALSE]
  env_final <- env_final[, which(colSums(is.na(env_final)) < nrow(env_final)), with = FALSE]
  env_final <- env_final[, which(sapply(env_final, sd, na.rm = TRUE) > 0), with = FALSE]
  
  env_final <- as.data.table(scale(env_final))
  
  # Vérifications
  if (any(is.na(geno_final)) || any(is.na(env_final))) stop("❌ Données incomplètes.")
  
  # RDA
  message("🚀 Lancement de la RDA...")
  rda_model <- rda(geno_final ~ ., data = env_final)
  print(summary(rda_model))
  return(rda_model)
}

# Chargement des données environnementales
env_data_full <- fread("Filtered_ACP_Lobster_with_Lat_Env.txt")
cols_env <- c("FID", "Bathymetrie_moy", "Temperature_max", "Temperature_moy", "Temperature_min", 
              "Dissolution_max", "Dissolution_moy", "Dissolution_min", 
              "Salinity_max", "Salinity_moy", "Salinity_min")
env_data <- env_data_full[, ..cols_env]

# Lancement des modèles RDA pour les PCs d'intérêt
rda_PC2 <- run_rda("PC2_outliers", copy(env_data), "PC2")
rda_PC3 <- run_rda("PC3_outliers", copy(env_data), "PC3")
rda_PC4 <- run_rda("PC4_outliers", copy(env_data), "PC4")

var_exp <- summary(run_rda)$constrained$importance["Proportion Explained", ]

plot_rda_with_var <- function(rda_model, rda_label) {
  # Extraire la variance expliquée
  var_exp <- summary(rda_model)$constrained$importance[2, ]
  x_lab <- paste0("RDA1 (", round(var_exp[1] * 100, 1), "%)")
  y_lab <- paste0("RDA2 (", round(var_exp[2] * 100, 1), "%)")
  
  # Triplot global
  plot(rda_model, type = "n", main = paste("RDA -", rda_label), xlab = x_lab, ylab = y_lab)
  points(rda_model, display = "sites", col = "grey30", pch = 16, cex = 0.8)
  text(rda_model, display = "bp", col = "blue", cex = 1)
  text(rda_model, display = "species", col = "red", cex = 0.6)
  
  # Sites
  plot(rda_model, type = "n", main = paste("RDA - Indiv -", rda_label), xlab = x_lab, ylab = y_lab)
  points(rda_model, display = "sites", col = "black", pch = 16, cex = 0.8)
  
  # Environnement
  plot(rda_model, type = "n", main = paste("RDA - Environnement -", rda_label), xlab = x_lab, ylab = y_lab)
  text(rda_model, display = "bp", col = "blue", cex = 1)
  
  # SNPs
  plot(rda_model, type = "n", main = paste("RDA - SNPs -", rda_label), xlab = x_lab, ylab = y_lab)
  text(rda_model, display = "species", col = "red", cex = 0.5)
}

# Afficher les résultats pour PC2, PC3 et PC4
plot_rda_with_var(rda_PC2, "PC2")
plot_rda_with_var(rda_PC3, "PC3")
plot_rda_with_var(rda_PC4, "PC4")
