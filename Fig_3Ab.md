# R
## Lire le fichier .bim
```
bim <- read.table("Lobster_no_024712526.bim", header = FALSE)
```
## Vérification du nombre total de SNPs disponibles
```
cat("Nombre total de SNPs dans le fichier :", nrow(bim), "\n")
```
## Tirer 50 000 SNPs au hasard (sans remplacement)
```
set.seed(42)  # pour rendre les tirages reproductibles
subset_snps <- sample(bim$V2, 50000)
```
## Sauvegarder la liste dans un fichier texte
```
write.table(subset_snps, file = "snps_subset.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
```
# Bash
## Faire les fichiers .bed, .bim et .fam pour le sous-ensemble
```
./plink --bfile Lobster_no_024712526 --extract snps_subset.txt --make-bed --out Lobster_50000SNPs --allow-extra-chr
```
## Lancer ADMiXTURE
```
admixture Lobster_50000SNPs.bed K
```
# R
## Charger les librairies
```
library(ggplot2)
library(dplyr)
library(tidyr)
library(readxl)
```
## PARAMÈTRES
```
basename <- "Lobster_5000SNPs"        # préfixe des fichiers PLINK/ADMIXTURE
K_values <- 2:10                      # valeurs de K à tester
fam_file <- paste0(basename, ".fam") # fichier .fam
meta_file <- "Merged_Lobster_Data.xlsx"  # fichier avec Sample_ID et Population_ID
```
## IMPORT DES DONNÉES
### Lire le .fam
```
fam <- read.table(fam_file)
colnames(fam)[2] <- "Sample_ID"
```
### Lire les métadonnées (attention à l'encodage des colonnes)
```
meta <- read_excel(meta_file)
meta <- meta %>% rename(Sample_ID = `Sample ID 2`)
```
### Harmoniser les IDs (si besoin)
```
fam$Sample_ID <- as.character(fam$Sample_ID)
meta$Sample_ID <- as.character(meta$Sample_ID)
```
### Fusion
```
fam_meta <- left_join(fam, meta, by = "Sample_ID")
```
### Créer un dossier de sortie
```
dir.create("ADMIXTURE_plots", showWarnings = FALSE)
```
## BOUCLE SUR LES K
```
for (K in K_values) {
  
  Qfile <- paste0(basename, ".", K, ".Q")
  if (!file.exists(Qfile)) {
    cat("Fichier manquant pour K =", K, "\n")
    next
  }
  
  Q <- read.table(Qfile)
  colnames(Q)[1:K] <- paste0("Cluster", 1:K)
  
  Q$Sample_ID <- fam_meta$Sample_ID
  Q$Population <- fam_meta$`Population ID`
  
  # Format long pour ggplot
  Q_long <- Q %>%
    pivot_longer(cols = starts_with("Cluster"), 
                 names_to = "Cluster", values_to = "Ancestry")
  
  Q_long$Sample_ID <- factor(Q_long$Sample_ID, levels = fam_meta$Sample_ID)
```  
  ## Créer le graphique
  ```
  p <- ggplot(Q_long, aes(x = Sample_ID, y = Ancestry, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1) +
    facet_grid(~ Population, scales = "free_x", space = "free_x") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      strip.text = element_text(size = 6)
    ) +
    labs(title = paste("Structure génétique – K =", K),
         x = NULL, y = "Proportion d'ancestralité")
  ```
  ## Sauvegarder
  ```
  ggsave(filename = paste0("ADMIXTURE_plots/admixture_K", K, ".png"),
         plot = p, width = 10, height = 4 + K * 0.2, dpi = 300)
  
  cat("Graphique pour K =", K, "enregistré.\n")
}
```
# Graphique par population
## Charger les librairies
```
library(tidyverse)
library(readxl)
```
## Paramètres 
```
K <- 9
qfile <- paste0("lobster_top5000.", K, ".Q")
famfile <- "lobster_top5000.fam"
meta_xlsx <- "Merged_Lobster_Data.xlsx"
```
## Charger les fichiers ADMIXTURE
```
Q <- read.table(qfile)
colnames(Q) <- paste0("Cluster", 1:K)
fam <- read.table(famfile)
colnames(fam)[1] <- "FID"
```
## Charger les métadonnées
```
meta <- read_excel(meta_xlsx)
```
### S'assurer que les ID sont bien comparables
```
meta <- meta %>%
  mutate(FID = as.character(`Sample ID 2`))
```
### Ajouter les ID au Q matrix
```
Q <- Q %>%
  mutate(FID = fam$FID) %>%
  left_join(meta[, c("FID", "Population ID")], by = "FID") %>%
  rename(Population = `Population ID`) %>%
  arrange(Population) %>%
  mutate(Indiv = factor(FID, levels = FID))
```
## Long format pour ggplot
```
Q_long <- Q %>%
  pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Ancestry")
```
## ADMIXTURE plot par Population ID 
```
ggplot(Q_long, aes(x = Indiv, y = Ancestry, fill = Cluster)) +
  geom_bar(stat = "identity", width = 1) +
  facet_grid(~ Population, scales = "free_x", space = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = paste("ADMIXTURE plot (K =", K, ") par population"),
       x = "Individus",
       y = "Proportion d’ancestralité") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(hjust = 0.5))
```
