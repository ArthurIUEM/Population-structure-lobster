# Extraction des données environnementales

## Charger les bibliothèques
library(tidyverse)
library(data.table)
library(readxl)
library(sdmpredictors)
library(raster)

## Lire le fichier Excel avec FID, latitude et longitude
env_points <- read_excel("Data for environmental data.xlsx")

## Harmoniser les noms de colonnes
colnames(env_points) <- tolower(colnames(env_points))
names(env_points)[names(env_points) %in% c("lat", "latitude")] <- "lat"
names(env_points)[names(env_points) %in% c("lon", "longitude")] <- "lon"

## Éviter les timeouts réseau
options(timeout = max(1000, getOption("timeout")))

## Télécharger les couches environnementales du fond marin
environment.bottom <- load_layers(c("BO2_tempmax_bdmean",
                                    "BO2_tempmean_bdmean",
                                    "BO2_tempmin_bdmean",
                                    "BO2_dissoxmax_bdmean",
                                    "BO2_dissoxmean_bdmean",
                                    "BO2_dissoxmin_bdmean",
                                    "BO2_salinitymax_bdmean",
                                    "BO2_salinitymean_bdmean",
                                    "BO2_salinitymin_bdmean"))

## Télécharger la couche de bathymétrie
bathymetry <- load_layers("BO_bathymean")

## Extraire les coordonnées
coordinates <- env_points[, c("lon", "lat")]

## Extraire les valeurs environnementales
env_extracted <- extract(environment.bottom, coordinates)
depths <- extract(bathymetry, coordinates)

## Créer le tableau final
Env <- data.frame(FID = env_points$fid,
                  lat = env_points$lat,
                  lon = env_points$lon,
                  depth = depths,
                  env_extracted)

## Sauvegarder en TSV
write.table(Env,
            "KessCod2020_Env.tsv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

# Top 10K

## Filtre  

### Lire le fichier FST
fst <- read.table("lobster_fst.fst", header = TRUE)

### Vérifier les premières lignes
head(fst)

### Trier par FST décroissant
fst_top <- fst[order(-fst$FST), ]

### Garder les 10 000 SNPs les plus différenciés
top10k_snps <- fst_top[1:10000, "SNP"]

### Sauvegarder dans un fichier texte pour PLINK
write.table(top10k_snps, "top10k_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Creation des fichiers .bim, .bed et .fam 

 ./plink --bfile Lobster_no_024712526 \
      --extract top10k_snps.txt \
      --make-bed \
      --out Lobster_top10k_fst --allow-extra-chr

## Preparation des fichiers pour la RDA

### big SNPR way 
library(bigsnpr)
plink_file <- "Lobster_top10k_fst"  # no extension
bed <- snp_readBed(paste0(plink_file, ".bed"))
obj.bigSNP <- snp_attach(paste0(plink_file, ".rds"))
G   <- obj.bigSNP$genotypes  # genotype matrix
fam <- obj.bigSNP$fam        # FAM file data
map <- obj.bigSNP$map        # BIM file data
#### OR TRY SEQARRAY
library(SeqArray)
library(SeqVarTools)
### make a gds obejct from you rbed file
seqBED2GDS(
  bed.fn = "Lobster_top10k_fst.bed",   # Path to .bed file
  bim.fn = "Lobster_top10k_fst.bim",   # Path to .bim file
  fam.fn = "Lobster_top10k_fst.fam",   # Path to .fam file
  out.gdsfn = "Lobster_top10k_fst.gds"    # Output GDS file, call it whatever you wwant
)
### read in your gds file
gdsfmt::showfile.gds(closeall=TRUE)
filename.gds <- "Lobster_top10k_fst.gds"
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

env_vars_clean <- env_vars %>%
  mutate(across(everything(), ~ ifelse(is.na(.), impute_mean(.), .))) %>%
  mutate(Latitude = env_lat)

### Synchroniser la matrice génotypique avec les mêmes lignes
geno_clean <- dos.all.input[match(env_ordered$FID, rownames(dos.all.input)), ]

### Trouver les lignes sans NA dans les deux jeux de données
complete_rows <- complete.cases(env_vars_clean) & complete.cases(geno_clean)

### Appliquer ce filtre aux deux matrices
env_vars_clean <- env_vars_clean[complete_rows, ]
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

rda_result <- rda(geno_clean_imputed ~ ., data = env_vars_clean)

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
colors <- colorRampPalette(c("blue", "green", "red"))(100)
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
     ylim = c(-9, 9),
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
bim <- read.table("Lobster_top50k_fst.bim", header = FALSE)
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
  ylim = c(0, 0.2),
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
  ylim = c(0, 0.2),
  main = "Manhattan Plot - RDA2 (abs loadings) (Top 10K)",
  ylab = "abs(loading)"
)

      
# Top 25K

## Filtre 

### Lire le fichier FST
fst <- read.table("lobster_fst.fst", header = TRUE)

### Vérifier les premières lignes
head(fst)

### Trier par FST décroissant
fst_top <- fst[order(-fst$FST), ]

### Garder les 25 000 SNPs les plus différenciés
top10k_snps <- fst_top[1:25000, "SNP"]

### Sauvegarder dans un fichier texte pour PLINK
write.table(top10k_snps, "top25k_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Creation des fichiers .bim, .bed et .fam 

 ./plink --bfile Lobster_no_024712526 \
      --extract top25k_snps.txt \
      --make-bed \
      --out Lobster_top25k_fst --allow-extra-chr

## Preparation des fichiers pour la RDA

### big SNPR way 

library(bigsnpr)
plink_file <- "Lobster_top25k_fst"  # no extension
bed <- snp_readBed(paste0(plink_file, ".bed"))
obj.bigSNP <- snp_attach(paste0(plink_file, ".rds"))
G   <- obj.bigSNP$genotypes  # genotype matrix
fam <- obj.bigSNP$fam        # FAM file data
map <- obj.bigSNP$map        # BIM file data

#### OR TRY SEQARRAY

library(SeqArray)
library(SeqVarTools)

### make a gds obejct from you rbed file

seqBED2GDS(
  bed.fn = "Lobster_top25k_fst.bed",   # Path to .bed file
  bim.fn = "Lobster_top25k_fst.bim",   # Path to .bim file
  fam.fn = "Lobster_top25k_fst.fam",   # Path to .fam file
  out.gdsfn = "Lobster_top25k_fst.gds"    # Output GDS file, call it whatever you wwant
)

### read in your gds file

gdsfmt::showfile.gds(closeall=TRUE)
filename.gds <- "Lobster_top25k_fst.gds"
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

env_vars_clean <- env_vars %>%
  mutate(across(everything(), ~ ifelse(is.na(.), impute_mean(.), .))) %>%
  mutate(Latitude = env_lat)

### Synchroniser la matrice génotypique avec les mêmes lignes
geno_clean <- dos.all.input[match(env_ordered$FID, rownames(dos.all.input)), ]

### Trouver les lignes sans NA dans les deux jeux de données
complete_rows <- complete.cases(env_vars_clean) & complete.cases(geno_clean)

### Appliquer ce filtre aux deux matrices
env_vars_clean <- env_vars_clean[complete_rows, ]
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

rda_result <- rda(geno_clean_imputed ~ ., data = env_vars_clean)

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
colors <- colorRampPalette(c("blue", "green", "red"))(100)
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
     ylim = c(-9, 9),
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
bim <- read.table("Lobster_top50k_fst.bim", header = FALSE)
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
  ylim = c(0, 0.2),
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
  ylim = c(0, 0.2),
  main = "Manhattan Plot - RDA2 (abs loadings) (Top 10K)",
  ylab = "abs(loading)"
)

# Top 50K

## Filtre

### Lire le fichier FST
fst <- read.table("lobster_fst.fst", header = TRUE)

### Vérifier les premières lignes
head(fst)

### Trier par FST décroissant
fst_top <- fst[order(-fst$FST), ]

### Garder les 50 000 SNPs les plus différenciés
top10k_snps <- fst_top[1:50000, "SNP"]

### Sauvegarder dans un fichier texte pour PLINK
write.table(top10k_snps, "top50k_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Creation des fichiers .bim, .bed et .fam 

 ./plink --bfile Lobster_no_024712526 \
      --extract top50k_snps.txt \
      --make-bed \
      --out Lobster_top50k_fst --allow-extra-chr

## Preparation des fichiers pour la RDA

### big SNPR way 
library(bigsnpr)
plink_file <- "Lobster_top50k_fst"  # no extension
bed <- snp_readBed(paste0(plink_file, ".bed"))
obj.bigSNP <- snp_attach(paste0(plink_file, ".rds"))
G   <- obj.bigSNP$genotypes  # genotype matrix
fam <- obj.bigSNP$fam        # FAM file data
map <- obj.bigSNP$map        # BIM file data
#### OR TRY SEQARRAY
library(SeqArray)
library(SeqVarTools)
### make a gds obejct from you rbed file
seqBED2GDS(
  bed.fn = "Lobster_top50k_fst.bed",   # Path to .bed file
  bim.fn = "Lobster_top50k_fst.bim",   # Path to .bim file
  fam.fn = "Lobster_top50k_fst.fam",   # Path to .fam file
  out.gdsfn = "Lobster_top50k_fst.gds"    # Output GDS file, call it whatever you wwant
)
### read in your gds file
gdsfmt::showfile.gds(closeall=TRUE)
filename.gds <- "Lobster_top50k_fst.gds"
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

env_vars_clean <- env_vars %>%
  mutate(across(everything(), ~ ifelse(is.na(.), impute_mean(.), .))) %>%
  mutate(Latitude = env_lat)

### Synchroniser la matrice génotypique avec les mêmes lignes
geno_clean <- dos.all.input[match(env_ordered$FID, rownames(dos.all.input)), ]

### Trouver les lignes sans NA dans les deux jeux de données
complete_rows <- complete.cases(env_vars_clean) & complete.cases(geno_clean)

### Appliquer ce filtre aux deux matrices
env_vars_clean <- env_vars_clean[complete_rows, ]
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

rda_result <- rda(geno_clean_imputed ~ ., data = env_vars_clean)

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
colors <- colorRampPalette(c("blue", "green", "red"))(100)
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
     ylim = c(-9, 9),
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
bim <- read.table("Lobster_top50k_fst.bim", header = FALSE)
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
  ylim = c(0, 0.2),
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
  ylim = c(0, 0.2),
  main = "Manhattan Plot - RDA2 (abs loadings) (Top 10K)",
  ylab = "abs(loading)"
)
