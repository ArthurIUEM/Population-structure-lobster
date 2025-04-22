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

# Générer un fichier avec uniquement les 150k meilleurs SNPs
./plink --bfile Lobster1MB --extract top150k_snps.txt --make-bed --out Lobster1MB_topFST --allow-extra-chr

# Sous-échantillonnage Nord/Sud
awk 'NR>1 && $5 == "Nord" {print $1, $1}' UMAP_zones_latitude.tsv > nord_ids.txt
awk 'NR>1 && $5 == "Sud"  {print $1, $1}' UMAP_zones_latitude.tsv > sud_ids.txt

plink2 --bfile Lobster1MB_topFST --keep nord_ids.txt --make-bed --out Lobster_Nord_topFST
plink2 --bfile Lobster1MB_topFST --keep sud_ids.txt --make-bed --out Lobster_Sud_topFST

# Recode en .raw pour l'analyse RDA
plink2 --bfile Lobster_Nord_topFST --recode A --out Lobster_Nord_topFST
plink2 --bfile Lobster_Sud_topFST --recode A --out Lobster_Sud_topFST

