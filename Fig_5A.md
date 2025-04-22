library(readr)
library(dplyr)

# Charger les FST
fst_data <- read_tsv("lobster_fst.fst")

# Garder les 150000 SNPs avec FST les plus élevés
top_fst <- fst_data %>%
  arrange(desc(FST)) %>%
  slice(1:150000)

# Sauvegarder la liste des SNPs à garder
write_tsv(top_fst %>% select(SNP), "top150k_snps.txt", col_names = FALSE)

# Filtrer les SNPs
./plink --bfile Lobster_no_024712526 \
      --extract top10k_snps.txt \
      --make-bed \
      --out Lobster_top10k –allow-extra-chr

# Créer le fichier .raw pour dosage
./plink --bfile Lobster_top10k \
      --recode A \
      --out Lobster_top10k –allow-extra-chr

# 📦 Chargement des packages
library(data.table)
library(dplyr)
library(vegan)

# 🔧 Fonction pour retirer les colonnes constantes
remove_constant_columns <- function(df) {
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}

# 🔧 Fonction principale pour faire une RDA sur un sous-ensemble
do_rda <- function(data, titre = "RDA") {
  # 🧬 Sélection génotypique
  geno <- data %>% select(where(is.numeric), -FID)
  
  # 🌿 Variables environnementales autorisées
  env_var_names <- c(
    "Bathymetrie_moy",
    "Temperature_max", "Temperature_moy", "Temperature_min",
    "Dissolution_max", "Dissolution_moy", "Dissolution_min",
    "Salinity_max", "Salinity_moy", "Salinity_min"
  )
  
  env_vars <- data %>%
    select(all_of(env_var_names)) %>%
    remove_constant_columns() %>%
    mutate(across(everything(), as.numeric))
  
  # 🧭 RDA
  rda_res <- rda(geno ~ ., data = env_vars)
  
  # 🎯 Variance expliquée
  eig_vals <- summary(rda_res)$cont$importance[2, 1:2] * 100
  xlab_txt <- paste0("RDA1 (", round(eig_vals[1], 1), "%)")
  ylab_txt <- paste0("RDA2 (", round(eig_vals[2], 1), "%)")
  
  # 📊 Plot
  plot(rda_res, scaling = 2,
       main = titre,
       xlim = if (titre == "RDA - Nord") c(-40, 40) else c(-40, 40),
       ylim = if (titre == "RDA - Nord") c(-20, 20) else c(-20, 20),
       xlab = xlab_txt,
       ylab = ylab_txt)
  
  return(rda_res)
}

# ================================
# 1. Chargement des données
# ================================

# 📂 Fichier .raw
genos <- fread("Lobster_top10k.raw", sep = " ")

# 🧼 Nettoyage des colonnes
genos.dose <- genos %>% select(-PAT, -MAT, -SEX, -PHENOTYPE)

# 🧩 Imputation (valeur allèle la plus fréquente par SNP)
genos.dose.imp <- genos.dose %>%
  select(-FID) %>%
  as.data.frame() %>%
  apply(2, function(x) {
    x[is.na(x)] <- as.numeric(names(which.max(table(x))))
    return(x)
  }) %>%
  as.data.frame()

# Réattacher les FID
genos.dose.imp$FID <- genos$FID

# 🌿 Données environnementales
env <- fread("Filtered_ACP_Lobster_with_Lat_Env_clean.tsv")

# 🌍 Zones Nord/Sud
zones <- fread("UMAP_zones_latitude.tsv")  # Doit contenir FID et Zone

# 🧷 Fusion
data_all <- genos.dose.imp %>%
  inner_join(env, by = "FID") %>%
  inner_join(zones, by = "FID")

# ================================
# 2. Séparation par Zone
# ================================

data_nord <- data_all %>% filter(ZONE.x == "Nord")
data_sud  <- data_all %>% filter(ZONE.x == "Sud")

# ================================
# 3. RDA
# ================================

rda_nord <- do_rda(data_nord, titre = "RDA - Nord")
rda_sud  <- do_rda(data_sud, titre = "RDA - Sud")


# 📦 Chargement des packages
library(data.table)
library(dplyr)
library(vegan)

# 🔧 Fonction pour retirer les colonnes constantes
remove_constant_columns <- function(df) {
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}

# 🔧 Fonction principale pour faire une RDA sur un sous-ensemble
do_rda <- function(data, titre = "RDA") {
  # 🧬 Sélection génotypique
  geno <- data %>% select(where(is.numeric), -FID)
  
  # 🌿 Variables environnementales autorisées
  env_var_names <- c(
    "Bathymetrie_moy",
    "Temperature_max", "Temperature_moy", "Temperature_min",
    "Dissolution_max", "Dissolution_moy", "Dissolution_min",
    "Salinity_max", "Salinity_moy", "Salinity_min"
  )
  
  env_vars <- data %>%
    select(all_of(env_var_names)) %>%
    remove_constant_columns() %>%
    mutate(across(everything(), as.numeric))
  
  # 🧭 RDA
  rda_res <- rda(geno ~ ., data = env_vars)
  
  # 📍 Scores des individus
  sites_scores <- scores(rda_res, display = "sites", scaling = 2)
  plot_data <- as.data.frame(sites_scores)
  plot_data$LATITUDE <- data$LATITUDE.x

  # 🎨 Couleur personnalisée
  color_palette <- c("blue", "green", "yellow", "orange", "red")
  p <- ggplot() +
    geom_point(data = plot_data, aes(x = RDA1, y = RDA2, color = LATITUDE), size = 2, alpha = 0.8) +
    scale_color_gradientn(colors = color_palette) +
    labs(title = titre, x = "RDA1", y = "RDA2", color = "Latitude") +
    theme_minimal(base_size = 14)

  # 🔍 Limites personnalisées
  if (titre == "RDA - Nord") {
    p <- p + xlim(-40, 20) + ylim(-20, 25)
  } else {
    p <- p + xlim(-40, 35) + ylim(-20, 30)
  }

  print(p)
  return(rda_res)
}
# ================================
# 1. Chargement des données
# ================================

# 📂 Fichier .raw
genos <- fread("Lobster_top10k.raw", sep = " ")

# 🧼 Nettoyage des colonnes
genos.dose <- genos %>% select(-PAT, -MAT, -SEX, -PHENOTYPE)

# 🧩 Imputation (valeur allèle la plus fréquente par SNP)
genos.dose.imp <- genos.dose %>%
  select(-FID) %>%
  as.data.frame() %>%
  apply(2, function(x) {
    x[is.na(x)] <- as.numeric(names(which.max(table(x))))
    return(x)
  }) %>%
  as.data.frame()

# Réattacher les FID
genos.dose.imp$FID <- genos$FID

# 🌿 Données environnementales
env <- fread("Filtered_ACP_Lobster_with_Lat_Env_clean.tsv")

# 🌍 Zones Nord/Sud
zones <- fread("UMAP_zones_latitude.tsv")  # Doit contenir FID et Zone

# 🧷 Fusion
data_all <- genos.dose.imp %>%
  inner_join(env, by = "FID") %>%
  inner_join(zones, by = "FID")

# ================================
# 2. Séparation par Zone
# ================================

data_nord <- data_all %>% filter(ZONE.x == "Nord")
data_sud  <- data_all %>% filter(ZONE.x == "Sud")

# ================================
# 3. RDA
# ================================

rda_nord <- do_rda(data_nord, titre = "RDA - Nord")
rda_sud  <- do_rda(data_sud, titre = "RDA - Sud")
