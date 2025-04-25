# Créer le fichier de groupes
```
awk 'BEGIN{OFS="\t"} NR>1 {print $1, $1, $5}' UMAP_zones_latitude.tsv > lobster_clusters.txt
```
# Calculer le FST avec PLINK
```
./plink --bfile Lobster_no_024712526 \
      --fst \
      --within lobster_clusters.txt \
      --out lobster_fst --allow-extra-chr
```
# Script R pour le Manhattan plot

## Charger les librairies
```
library(data.table)
library(ggplot2)
```
## Charger les fichiers
```
fst <- fread("lobster_fst.fst", header=TRUE)
bim <- fread("Lobster_no_024712526.bim", header=FALSE)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
```
## Merge pour récupérer les positions
```
fst_data <- merge(fst, bim[, .(SNP, BP)], by = "SNP")
```
## Nettoyage
```
fst_data[, CHR := as.factor(CHR)]
fst_data[, CHR := as.numeric(as.character(CHR))]
fst_data <- fst_data[!is.na(CHR)]
```
## Ordre et position cumulée
```
fst_data <- fst_data[order(CHR, BP)]
chrom_cumlen <- fst_data[, .(maxBP = max(BP)), by = CHR]
chrom_cumlen[, offset := cumsum(shift(maxBP, fill=0))]
fst_data <- merge(fst_data, chrom_cumlen[, .(CHR, offset)], by = "CHR")
fst_data[, pos_cum := BP + offset]
```
## Position des étiquettes
```
axis_df <- fst_data[, .(center = mean(pos_cum)), by = CHR]
```
## Plot
```
ggplot(fst_data, aes(x = pos_cum, y = FST, color = as.factor(CHR))) +
  geom_point(size = 0.6) +
  scale_color_manual(values = rep(c("black", "grey60"), length.out=length(unique(fst_data$CHR.x)))) +
  scale_x_continuous(label = axis_df$CHR.x, breaks = axis_df$center) +
  labs(x = "Chromosome", y = "FST", title = "Manhattan plot of FST (Nord vs Sud)") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)  # <--- texte en diagonale
  )

```
# Sortir le tableau
unique_chr_table <- axis_df[order(center), .(CHR = CHR)]
fwrite(unique_chr_table, "liste_chromosomes_ordonnees.tsv", sep = "\t")
