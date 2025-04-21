library(tidyverse)

# Charger les FST
fst_data <- read.table("lobster_fst.fst", header = TRUE)

# Trier les SNPs par FST décroissant
fst_sorted <- fst_data %>% arrange(desc(FST))

# Garder le top 1%
n_top <- ceiling(0.01 * nrow(fst_sorted))
top_snps <- fst_sorted[1:n_top, "SNP"]  # Remplace "SNP" si le nom de colonne est différent

# Sauvegarder dans un fichier
write.table(top_snps, "top1pct_snps.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

./plink --bfile Lobster1MB --extract top1pct_snps.txt --recodeA --out top1pct_genotypes --allow-extra-chr
