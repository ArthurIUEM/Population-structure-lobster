# Charger les bibliothèques
library(ggplot2)
library(dplyr)
library(readr)
library(gridExtra)

# Lire le fichier brut (sans entête)
raw_data <- read_table2("PCA_results.eigenvec.var", col_names = FALSE)

# Ajouter les noms de colonnes
colnames(raw_data) <- c("CHROM", "ID", "MAJ", "NONMAJ", paste0("PC", 1:10))

# Convertir les colonnes PC1 à PC10 en numérique (au cas où il y a du texte dedans)
for (i in 1:10) {
  raw_data[[paste0("PC", i)]] <- as.numeric(raw_data[[paste0("PC", i)]])
}

# Supprimer les lignes avec des NA dans PC2
pc_loading <- raw_data %>% filter(!is.na(PC2))

# Ajouter une position locale par chromosome
pc_loading <- pc_loading %>%
  group_by(CHROM) %>%
  mutate(Pos_cum = row_number()) %>%
  ungroup()

# Calculer le décalage cumulatif des chromosomes pour positionner sur l’axe x
chr_offset <- pc_loading %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(Pos_cum)) %>%
  mutate(offset = cumsum(chr_len) - chr_len) %>%
  select(CHROM, offset)

# Ajouter la position cumulée
pc_loading <- pc_loading %>%
  inner_join(chr_offset, by = "CHROM") %>%
  mutate(pos_cum = Pos_cum + offset)

# Points de repère pour l'axe x (centres des chromosomes)
axis_df <- pc_loading %>%
  group_by(CHROM) %>%
  summarise(center = mean(pos_cum))

# Manhattan plot global (PC2 brut)
p1 <- ggplot(pc_loading, aes(x = pos_cum, y = abs(PC2), color = as.factor(CHROM))) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = rep(c("gray20", "skyblue"), length.out = length(unique(pc_loading$CHROM)))) +
  scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$center) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 18, face = "bold")
  ) +
  labs(
    x = "Contigs",
    y = "|PC2 loading|",
    title = "(A) Manhattan plot of PC2 loadings"
  )

# Liste des loci d'intérêt
highlight_positions <- c(5435, 8796, 9929, 19067, 21423, 21972, 24398, 25413, 29272, 35994, 44179)

# Filtrer le chromosome sexuel
chr_sex <- "NW_024712526.1"
sex_chr_data <- pc_loading %>%
  filter(CHROM == chr_sex) %>%
  mutate(highlight = ifelse(Pos_cum %in% highlight_positions, "highlight", "normal"))

# Zoom sur le chromosome sexuel avec loci en rouge
# Zoom sur le chromosome sexuel avec loci d'intérêt bien visibles
p2 <- ggplot(sex_chr_data, aes(x = Pos_cum, y = PC2)) +
  # D'abord les points "normaux" en gris
  geom_point(data = subset(sex_chr_data, highlight == "normal"),
             color = "gray40", alpha = 0.6, size = 1.5) +
  # Ensuite les points d'intérêt en rouge (passent au premier plan)
  geom_point(data = subset(sex_chr_data, highlight == "highlight"),
             color = "red", alpha = 0.9, size = 2.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "(B) Zoom in on interest contig (NW_024712526.1) with sex loci (in red) identify by Benestan et al. (2017)",
    x = "Position on the contig",
    y = "PC2 loading"
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "none"
  )


# Afficher les deux graphiques en panel
grid.arrange(p1, p2, ncol = 1)
