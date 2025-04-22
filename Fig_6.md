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
