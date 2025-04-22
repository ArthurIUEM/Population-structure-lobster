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
  
  # 📊 Plot sans variance expliquée
  plot(rda_res, scaling = 2,
       main = titre,
       xlim = c(-40, 40),
       ylim = c(-20, 20))
  
  return(rda_res)
}

# ================================
# 1. Chargement des données
# ================================

# 📂 Fichier .raw
genos <- fread("Lobster_top10k.raw", sep = " ")

# 🧼 Nettoyage des colonnes inutiles
genos.dose <- genos %>% select(-PAT, -MAT, -SEX, -PHENOTYPE)

# 🧩 Imputation (remplacement des NA par l’allèle le plus fréquent par SNP)
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

# 🧷 Fusion de toutes les données
data_all <- genos.dose.imp %>%
  inner_join(env, by = "FID") %>%
  inner_join(zones, by = "FID")

# ================================
# 2. Séparation par Zone
# ================================

data_nord <- data_all %>% filter(ZONE.x == "Nord")
data_sud  <- data_all %>% filter(ZONE.x == "Sud")

# ================================
# 3. RDA par groupe
# ================================

rda_nord <- do_rda(data_nord, titre = "RDA - Nord")
rda_sud  <- do_rda(data_sud,  titre = "RDA - Sud")

# ================================
# 4. RDA sur tout le jeu de données
# ================================

rda_total <- do_rda(data_all, titre = "RDA - Total")



# 📦 Chargement des packages
library(data.table)
library(dplyr)
library(vegan)
library(ggplot2)

# 🔧 Fonction pour retirer les colonnes constantes
remove_constant_columns <- function(df) {
  df[, sapply(df, function(x) length(unique(x)) > 1)]
}

# 🎨 Fonction pour plot RDA avec ggplot2 et gradient latitude
plot_rda_by_lat <- function(rda_res, data, titre = "RDA") {
  scores_ind <- scores(rda_res, display = "sites", scaling = 2)
  df_plot <- as.data.frame(scores_ind)
  df_plot$LATITUDE.x <- data$LATITUDE.x
  
  ggplot(df_plot, aes(x = RDA1, y = RDA2, color = LATITUDE.x)) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_gradientn(colors = c("blue", "green", "yellow", "orange", "red")) +
    theme_minimal() +
    labs(title = titre, color = "Latitude")
}

# 🔧 Fonction principale pour faire une RDA et la retourner
do_rda <- function(data) {
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
  rda(geno ~ ., data = env_vars)
}

# ================================
# 1. Chargement des données
# ================================

# 📂 Fichier .raw
genos <- fread("Lobster_top10k.raw", sep = " ")
genos.dose <- genos %>% select(-PAT, -MAT, -SEX, -PHENOTYPE)

# 🧩 Imputation
genos.dose.imp <- genos.dose %>%
  select(-FID) %>%
  as.data.frame() %>%
  apply(2, function(x) {
    x[is.na(x)] <- as.numeric(names(which.max(table(x))))
    return(x)
  }) %>%
  as.data.frame()
genos.dose.imp$FID <- genos$FID

# 🌿 Environnement et zones
env <- fread("Filtered_ACP_Lobster_with_Lat_Env_clean.tsv")
zones <- fread("UMAP_zones_latitude.tsv")  # contient FID et Zone

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
# 3. RDA + Plot
# ================================

# Nord
rda_nord <- do_rda(data_nord)
plot_rda_by_lat(rda_nord, data_nord, "RDA - Nord")

# Sud
rda_sud <- do_rda(data_sud)
plot_rda_by_lat(rda_sud, data_sud, "RDA - Sud")

# Total
rda_total <- do_rda(data_all)
plot_rda_by_lat(rda_total, data_all, "RDA - Total")
