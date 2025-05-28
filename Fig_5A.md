# R
## Charger les données
```
df <- read.table("UMAP_zones_latitude.tsv", header = TRUE)
```
## Assumer que FID = IID (comme souvent dans PLINK)
```
df$IID <- df$FID
```
## Séparer les groupes
```
nord <- df[df$ZONE == "Nord", c("FID", "IID")]
sud <- df[df$ZONE == "Sud", c("FID", "IID")]
```
## Sauvegarder les fichiers
```
write.table(nord, "nord.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(sud, "sud.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
```
# Bash
## Pour le groupe Nord
```
./plink --bfile Lobster_no_024712526 --keep nord.txt --make-bed --out GONE_NORD --allow-extra-chr
```
## Pour le groupe Sud
```
./plink --bfile Lobster_no_024712526 --keep sud.txt --make-bed --out GONE_SUD --allow-extra-chr
```
## Nord
```
./plink --bfile GONE_NORD --recode --out GONE_NORD --allow-extra-chr
```
## Sud
```
./plink --bfile GONE_SUD --recode --out GONE_SUD --allow-extra-chr
```
## Lancer GONE2 pour le nord
```
./gone2 -r 1.1 Lobster_Nord_selected.tped
```
## Lancer GONE2 pour le sud
```
./gone2 -r 1.1 Lobster_Sud_selected.tped
```
# Lire les fichiers de sortie GONE2
nord <- read.table("Lobster_Nord_selected.tped_GONE2_Ne", header = TRUE)
sud  <- read.table("Lobster_Sud_selected.tped_GONE2_Ne", header = TRUE)

# Conversion Générations → Années
nord$Year <- 972 + 7 * nord$Generation
sud$Year  <- 972 + 7 * sud$Generation

# Déterminer la plage des années (commune pour cohérence visuelle)
year_range <- range(c(nord$Year, sud$Year))

# Afficher deux graphiques côte à côte
par(mfrow = c(1, 2), mar = c(7, 4, 4, 2))  # marge basse augmentée pour texte incliné

# Graphique Nord
plot(nord$Year, nord$Ne_diploids, type = "l", col = "steelblue",
     xlab = "Calendar year", ylab = "Effective population size (Ne)",
     main = "North zone", ylim = c(0, 52000), lwd = 2, xlim = year_range, xaxt = "n")

years_to_show_nord <- nord$Year[nord$Year %% 3 == 0]
axis(1, at = years_to_show_nord, labels = FALSE)
text(x = years_to_show_nord, y = par("usr")[3] - 1000,
     labels = years_to_show_nord, srt = 45, adj = 1, xpd = TRUE, cex = 0.8)

# Lignes verticales pour le nord
abline(v = 1455, col = "steelblue", lty = 2, lwd = 1.5)
abline(v = 1868, col = "steelblue", lty = 2, lwd = 1.5)

# Graphique Sud
plot(sud$Year, sud$Ne_diploids, type = "l", col = "darkorange",
     xlab = "Calendar year", ylab = "Effective population size (Ne)",
     main = "South zone", ylim = c(0, 400000), lwd = 2, xlim = year_range, xaxt = "n")

years_to_show_sud <- sud$Year[sud$Year %% 3 == 0]
axis(1, at = years_to_show_sud, labels = FALSE)
text(x = years_to_show_sud, y = par("usr")[3] - 10000,
     labels = years_to_show_sud, srt = 45, adj = 1, xpd = TRUE, cex = 0.8)

# Lignes verticales pour le sud
abline(v = 1413, col = "darkorange", lty = 2, lwd = 1.5)
abline(v = 1945, col = "darkorange", lty = 2, lwd = 1.5)

# Remettre le layout par défaut
par(mfrow = c(1, 1))
