# CHARGER LES PACKAGES

library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)

# === 1. CHARGER LES DONNÉES POUR LES GRAPHIQUES ACP + UMAP ===

df <- read.table("UMAP_with_clusters_trimmed.txt", header = TRUE, sep = "\t")
df$LATITUDE <- as.numeric(as.character(df$LATITUDE))

# Palette de couleurs pour la latitude
color_palette <- c("red", "orange", "yellow", "green", "blue")

# Thème personnalisé
custom_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

# Limites communes pour les axes
xlims <- range(c(df$PC1, df$PC2, df$PC3, df$PC4), na.rm = TRUE)
ylims <- xlims  # pour coord_fixed

# === 2. GRAPHES ACP ===

p1 <- ggplot(df, aes(x = PC1, y = PC2, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "(A) PC1 vs PC2", x = "PC1", y = "PC2", color = "Latitude") +
  custom_theme +
  coord_fixed() +
  xlim(xlims) + ylim(ylims)

p2 <- ggplot(df, aes(x = PC2, y = PC3, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "(B) PC2 vs PC3", x = "PC2", y = "PC3", color = "Latitude") +
  custom_theme +
  coord_fixed() +
  xlim(xlims) + ylim(ylims)

p3 <- ggplot(df, aes(x = PC3, y = PC4, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "(C) PC3 vs PC4", x = "PC3", y = "PC4", color = "Latitude") +
  custom_theme +
  coord_fixed() +
  xlim(xlims) + ylim(ylims)

# === 3. UMAP ===

p4 <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = LATITUDE)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_gradientn(colors = color_palette) +
  labs(title = "(D) UMAP1 vs UMAP2", x = "UMAP1", y = "UMAP2", color = "Latitude") +
  custom_theme +
  coord_fixed()

# === 4. MANHATTAN PLOT PC2 ===

# Lire les données PC2 loadings
raw_data <- read_table2("PCA_results.eigenvec.var", col_names = FALSE)
colnames(raw_data) <- c("CHROM", "ID", "MAJ", "NONMAJ", paste0("PC", 1:10))
for (i in 1:10) raw_data[[paste0("PC", i)]] <- as.numeric(raw_data[[paste0("PC", i)]])

pc_loading <- raw_data %>%
  filter(!is.na(PC2)) %>%
  group_by(CHROM) %>%
  mutate(Pos_cum = row_number()) %>%
  ungroup()

# Décalage cumulé pour position sur X
chr_offset <- pc_loading %>%
  group_by(CHROM) %>%
  summarise(chr_len = max(Pos_cum)) %>%
  mutate(offset = cumsum(chr_len) - chr_len)

pc_loading <- pc_loading %>%
  inner_join(chr_offset, by = "CHROM") %>%
  mutate(pos_cum = Pos_cum + offset)

axis_df <- pc_loading %>%
  group_by(CHROM) %>%
  summarise(center = mean(pos_cum))

# Manhattan global
p5 <- ggplot(pc_loading, aes(x = pos_cum, y = abs(PC2), color = as.factor(CHROM))) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = rep(c("gray20", "skyblue"), length.out = length(unique(pc_loading$CHROM)))) +
  scale_x_continuous(label = axis_df$CHROM, breaks = axis_df$center) +
  labs(title = "(E) Manhattan plot of PC2 loadings", x = "Contigs", y = "|PC2 loading|") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold")
  )

# === 5. ZOOM SEX CONTIG ===

highlight_positions <- c(5435, 8796, 9929, 19067, 21423, 21972, 24398, 25413, 29272, 35994, 44179)
chr_sex <- "NW_024712526.1"

sex_chr_data <- pc_loading %>%
  filter(CHROM == chr_sex) %>%
  mutate(highlight = ifelse(Pos_cum %in% highlight_positions, "highlight", "normal"))

p6 <- ggplot(sex_chr_data, aes(x = Pos_cum, y = PC2)) +
  geom_point(data = subset(sex_chr_data, highlight == "normal"), color = "gray40", alpha = 0.6, size = 1.5) +
  geom_point(data = subset(sex_chr_data, highlight == "highlight"), color = "red", alpha = 0.9, size = 2.5) +
  labs(
    title = "(F) Zoom on sex-linked contig (NW_024712526.1)",
    x = "Position on contig",
    y = "PC2 loading"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

# === 6. ASSEMBLAGE FINAL AVEC PATCHWORK ===

final_panel <- (
  (p1 | p2 | p3 | p4) /     # Ligne 1 : les 3 ACP
      p5   /          # Ligne 2 : UMAP + Manhattan PC2
     p6                   # Ligne 3 : Zoom contig sexuel
) + plot_layout(heights = c(1, 1, 1))  # Ajuste la hauteur relative

# AFFICHAGE
print(final_panel)
  
