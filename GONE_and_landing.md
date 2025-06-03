# === 1. Chargement des librairies ===
library(tidyverse)
library(readxl)
library(gridExtra)

# === 2. Chargement et préparation des landings ===

# Sud et Nord combinés à partir des fichiers d'origine
sGSL <- read.csv("LandingsByLFA_1892_2023.csv") %>%
  select(Year, Land_mt) %>%
  mutate(Population = "North")

mar_his <- read.csv("MaritimesHistoricLandings.csv") %>%
  mutate(Population = ifelse(LFA == 27, "North", "South")) %>%
  rename(Year = YEAR) %>%
  group_by(Year, Population) %>%
  summarize(Land_mt = sum(land_mt, na.rm = TRUE), .groups = "drop")

mar_cont <- read.csv("Maritimes_Annual_Landings.csv") %>%
  pivot_longer(LFA27:LFA38, names_to = "LFA", values_to = "land_mt") %>%
  mutate(Population = ifelse(LFA == "LFA27", "North", "South")) %>%
  group_by(Year, Population) %>%
  summarize(Land_mt = sum(land_mt, na.rm = TRUE), .groups = "drop")

nl <- read_xls("Landings_NL.xls", sheet = "Sheet1") %>%
  mutate(Population = "North")

# Création séparée pour le Nord et le Sud
land_north <- bind_rows(nl, mar_his, mar_cont, sGSL) %>%
  filter(Population == "North") %>%
  group_by(Year) %>%
  summarize(Landings = sum(Land_mt, na.rm = TRUE), .groups = "drop") %>%
  filter(Year >= 1892)

land_south <- bind_rows(mar_his, mar_cont) %>%
  filter(Population == "South") %>%
  group_by(Year) %>%
  summarize(Landings = sum(Land_mt, na.rm = TRUE), .groups = "drop")

# === 3. Chargement des données Ne (GONE) ===

ne_north <- read.table("Lobster_Nord_selected.tped_GONE2_Ne", header = TRUE) %>%
  mutate(Year = 972 + 7 * Generation)

ne_south <- read.table("Lobster_Sud_selected.tped_GONE2_Ne", header = TRUE) %>%
  mutate(Year = 972 + 7 * Generation)

# === 4. Création des graphiques ===

g1 <- ggplot(land_north, aes(x = Year, y = Landings)) +
  geom_line(color = "steelblue", size = 1) +
  labs(title = "(A) Estimated Landings - North", y = "Estimated Landings (t)", x = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

g2 <- ggplot(land_south, aes(x = Year, y = Landings)) +
  geom_line(color = "darkorange", size = 1) +
  labs(title = "(B) Estimated Landings - South", y = "Estimated Landings (t)", x = "") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

g3 <- ggplot(ne_north, aes(x = Year, y = Ne_diploids)) +
  geom_line(color = "steelblue", size = 1) +
  labs(title = "(C) Effective population size (Ne) - North", y = "Ne", x = "Year") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

g4 <- ggplot(ne_south, aes(x = Year, y = Ne_diploids)) +
  geom_line(color = "darkorange", size = 1) +
  labs(title = "(D) Effective population size (Ne) - South", y = "Ne", x = "Year") +
  theme_bw() +
  theme(plot.title = element_text(face = "bold"))

# === 5. Affichage en panel 2x2 ===
grid.arrange(g1, g2, g3, g4, ncol = 2, nrow = 2)
