# Extraire la zone et créer un fichier FID ZONE
cut -f1,5 UMAP_zones_latitude.tsv > lobster_zones.txt

# Créer les fichiers north_keep.txt et south_keep.txt au format FID IID
awk '$2 == "Nord" {print $1, $1}' lobster_zones.txt > north_keep.txt
awk '$2 == "Sud"  {print $1, $1}' lobster_zones.txt > south_keep.txt

awk '{print $1, $2, 1}' north_keep.txt > clusters.txt
awk '{print $1, $2, 2}' south_keep.txt >> clusters.txt

./plink --bfile Lobster_no_024712526 \
        --within clusters.txt \
        --fst \
        --allow-extra-chr \
        --out lobster_fst

library(ggplot2)
library(dplyr)

# Charger les données FST
fst <- read.table("lobster_fst.fst", header = TRUE)
fst_clean <- fst %>%
  na.omit() %>%
  mutate(CHR = as.factor(CHR))

# Générer la position cumulée
chr_info <- fst_clean %>%
  group_by(CHR) %>%
  summarize(chr_len = max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(CHR, tot)

fst_clean <- left_join(fst_clean, chr_info, by = "CHR") %>%
  mutate(BPcum = POS + tot)

# Obtenir le top 5000 SNPs
top_fst <- fst_clean %>%
  arrange(desc(FST)) %>%
  slice_head(n = 5000)

# Identifier les loci à visualiser
loci_with_top_snps <- unique(top_fst$CHR)

# Créer un dossier pour sauvegarder les plots (optionnel)
dir.create("plots_by_locus", showWarnings = FALSE)

# Générer un plot pour chaque locus contenant un top SNP
for (chr in loci_with_top_snps) {
  
  locus_data <- fst_clean %>% filter(CHR == chr)
  top_locus <- top_fst %>% filter(CHR == chr)
  
  p <- ggplot(locus_data, aes(x = POS, y = FST)) +
    geom_point(color = "grey70", size = 0.5) +
    geom_point(data = top_locus, aes(x = POS, y = FST), color = "red", size = 1) +
    theme_minimal() +
    labs(title = paste("Manhattan plot for locus", chr),
         x = "Position",
         y = expression(F[ST])) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0("plots_by_locus/Manhattan_", chr, ".png"),
         plot = p, width = 8, height = 4, dpi = 300)
}

library(ggplot2)
library(dplyr)

# Charger les données FST
fst <- read.table("lobster_fst.fst", header = TRUE)

# Nettoyer
fst_clean <- fst %>%
  na.omit() %>%
  mutate(CHR = as.factor(CHR))

# Sélection du top 5000 SNPs
top_fst <- fst_clean %>%
  arrange(desc(FST)) %>%
  slice_head(n = 5000)

# Panel avec un subplot par locus
ggplot(top_fst, aes(x = POS, y = FST)) +
  geom_point(color = "red", size = 0.6) +
  facet_wrap(~ CHR, scales = "free_x") +
  theme_minimal() +
  labs(title = "Top 5000 FST SNPs (par locus)",
       x = "Position (bp)",
       y = expression(F[ST])) +
  theme(strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(hjust = 0.5))

# Dans R
write.table(top_fst$SNP, "top5000_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
