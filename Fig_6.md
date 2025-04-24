# Charger les données
df <- read.table("UMAP_zones_latitude.tsv", header = TRUE)

# Assumer que FID = IID (comme souvent dans PLINK)
df$IID <- df$FID

# Séparer les groupes
nord <- df[df$ZONE == "Nord", c("FID", "IID")]
sud <- df[df$ZONE == "Sud", c("FID", "IID")]

# Sauvegarder les fichiers
write.table(nord, "nord.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(sud, "sud.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Pour le groupe Nord
./plink --bfile Lobster_no_024712526 --keep nord.txt --make-bed --out GONE_NORD --allow-extra-chr

# Pour le groupe Sud
./plink --bfile Lobster_no_024712526 --keep sud.txt --make-bed --out GONE_SUD --allow-extra-chr

# Nord
./plink --bfile GONE_NORD --recode --out GONE_NORD --allow-extra-chr

# Sud
./plink --bfile GONE_SUD --recode --out GONE_SUD --allow-extra-chr

# Lire les fichiers de sortie GONE2
nord <- read.table("Lobster_Nord_selected.tped_GONE2_Ne", header = TRUE)
sud  <- read.table("Lobster_Sud_selected.tped_GONE2_Ne", header = TRUE)

# Conversion de génération → année
nord$Year <- 2029 - 7 * nord$Generation
sud$Year  <- 2029 - 7 * sud$Generation

# Tracer la courbe pour la zone Nord
plot(nord$Year, nord$Ne_diploids, type = "l", col = "darksteelblueorange",
     xlab = "Calendar year", ylab = "Effective population size (Ne)",
     main = "Historical Ne estimates South zone",
     ylim = c(0, 400000), lwd = 2, xlim = c(min(nord$Year, sud$Year), 2022))

# Ajouter la courbe Sud
lines(sud$Year, sud$Ne_diploids, col = "darkorange", lwd = 2)

# Ajouter une légende
legend("topright", legend = c("North Zone", "South Zone"),
       col = c("steelblue", "darkorange"), lty = 1, lwd = 2)

plot(sud$Year, sud$Ne_diploids, type = "l", col = "darkorange",
     xlab = "Calendar year", ylab = "Effective population size (Ne)",
     main = "Historical Ne estimates South zone",
     ylim = c(0, 400000), lwd = 2, xlim = c(min(nord$Year, sud$Year), 2022))
