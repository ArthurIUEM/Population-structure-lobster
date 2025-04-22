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


