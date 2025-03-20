#Install packages 
install.packages(c("arsenal","tidyverse"))
install.packages("data.table")
install.packages("stats")
install.packages("ggplot2")
install.packages("plotly")
install.packages("readxl")
install.packages("openxlsx")
install.packages("dplyr")
#Jedes mal, wenn ich Paket nützen möchte 
library("arsenal")
library("tidyverse")
library("data.table")
library("stats")
library("ggplot2")
library("plotly")
library('readxl')
library("openxlsx")
library("dpylr")
#Datensatz laden 
df <- read_excel(file.choose())
#Strukturelle Übersicht anzeigen lassen 
str(df)
summary(df)
view(df)
#Logarithmus rückgängig machen mit 2^Zellinhalt
df <- df %>%
  mutate(across(6:ncol(.), ~ 2^.x))
view(df)
#Neue Excel Tabelle mit Rohdaten exportieren: 
write.xlsx(df, "V3_undo log2.xlsx", overwrite = TRUE)
#Mittelwerte aus zwei technischen replikaten
# Sicherstellen, dass ab Spalte 6 alles numerisch ist
df[6:ncol(df)] <- lapply(df[6:ncol(df)], function(x) as.numeric(as.character(x)))

# Mittelwerte für jede 2er-Gruppe berechnen
for (i in 6:ncol(df)) {  
  n_rows <- nrow(df)
  n_rows_even <- n_rows - (n_rows %% 2)  # Falls ungerade, kürzen
  
  # Sicherstellen, dass nur numerische Werte für die Berechnung verwendet werden
  column_data <- as.numeric(df[[i]][1:n_rows_even])  # Verwende df[[i]] für den Vektor
  
  # Berechnung der Mittelwerte aus den 2er-Gruppen
  mittelwerte <- rowMeans(matrix(column_data, ncol = 2, byrow = TRUE), na.rm = TRUE)
  
  # Neuer Spaltenname für den Mittelwert
  new_col_name <- paste0(colnames(df)[i], "_mean")
  
  # Neue Spalte für die Mittelwerte einfügen
  df[[new_col_name]] <- NA
  df[1:n_rows_even, new_col_name] <- rep(mittelwerte, each = 2)
}

# Ergebnis anzeigen
head(df)
#mittelwerte anzeigen lassen 
colnames(df)
View(df[, grepl("_mean$", colnames(df))])  # Zeigt nur die Mittelwert-Spalten
#speicher als neue Tabelle
library("openxlsx")
write.xlsx(df, "V3 means double", overwrite = TRUE)
#Nur noch mit der Mittelwert-Tabelle arbeiten
df_means <- df[, c(1:5, grep("_mean$", colnames(df)))]
View(df_means)
#Speichern
write.xlsx(df_means, "V2_mittelwerte_tabelle.xlsx", overwrite = TRUE)
#Lösche alle Zeilen, die bei Technical Replicate 2 stehen haben 
df_means <- df_means[df_means[[4]] != 2, ]
View(df_means)

#Speichere neue Tabelle
write.xlsx(df_means, "V3_mittelwerte_bereinigt.xlsx", overwrite = TRUE)
#Berechne den Mittelwerte der Expressionswerte für alle ALS pro Protein

#Normalisierte ALS werte: ALS1/Mean alle ALS für protein 
#Nur noch Mittelwerte bereinigt ALS reinladen
df_ALS <- read_excel(file.choose())
View(df_ALS)
# Berechnung des Mittelwerts für jede Spalte ab Spalte 6
column_means <- colMeans(df_ALS[, 6:ncol(df_ALS)], na.rm = TRUE)
# Normalisierung: Jede Zelle durch den Spaltenmittelwert teilen
df_ALS_normalized <- df_ALS
df_ALS_normalized[, 6:ncol(df_ALS)] <- sweep(df_ALS[, 6:ncol(df_ALS)], 2, column_means, "/")
View(df_ALS_normalized)
# Speichern der normalisierten Daten als Excel-Datei
write.xlsx(df_ALS_normalized, "V3_Normalisierte_Werte_ALS.xlsx", overwrite = TRUE)
# Anzeige der ersten Zeilen zur Überprüfung
head(df_ALS_normalized)


#Normalisiere CTRL: Mean CTRL 1/ Mean all CTRL 
#Lade andere Tabelle 
library('readxl')
df_CTRL <- read_excel(file.choose())
View(df_means)
#Normalisierung 
# Laden des Pakets für den Export
library(openxlsx)

# Berechnung des Mittelwerts für jede Spalte ab Spalte 6
column_means_ctrl <- colMeans(df_CTRL[, 6:ncol(df_CTRL)], na.rm = TRUE)

# Normalisierung: Jede Zelle durch den Spaltenmittelwert teilen
df_CTRL_normalized <- df_CTRL
df_CTRL_normalized[, 6:ncol(df_CTRL)] <- sweep(df_CTRL[, 6:ncol(df_CTRL)], 2, column_means_ctrl, "/")

# Speichern der normalisierten Daten als Excel-Datei
write.xlsx(df_CTRL_normalized, "Normalisierte_Werte_CTRL.xlsx", overwrite = TRUE)

# Anzeige der ersten Zeilen zur Überprüfung
head(df_CTRL_normalized)

#Scatterplot for ranked protein abundance in health 
#1.Mittelwert pro Protein berechnen 
#Lade Tabelle nicht normalisierte Einzelmittelwerte pro CTRL sample
df_rpa_CTRL <- read_excel(file.choose())
View(df_rpa_CTRL)
# Berechnung des Mittelwerts pro Protein für die Kontrollen ab Spalte 6
mean_protein_ctrl <- colMeans(df_CTRL[, 6:ncol(df_CTRL)], na.rm = TRUE)
View(mean_protein_ctrl)

# Erstelle ein DataFrame mit Spaltennamen und Mittelwerten
df_ranked_ctrl <- data.frame(Protein = names(mean_protein_ctrl), Abundance = mean_protein_ctrl)
# Sortieren der Proteine nach Abundance (absteigend)
df_ranked_ctrl <- df_ranked_ctrl[order(-df_ranked_ctrl$Abundance), ]
# Ergebnis anzeigen
View(df_ranked_ctrl)
# Speichern als Excel
library(openxlsx)
write.xlsx(df_ranked_ctrl, "Ranked_Protein_Abundance_CTRL.xlsx", overwrite = TRUE)
# Log10-Transformation der Abundance-Werte
df_ranked_ctrl$Log10_Abundance <- log10(df_ranked_ctrl$Abundance)
# Ergebnis anzeigen
View(df_ranked_ctrl)
# Optional: Speichern als Excel
library(openxlsx)
write.xlsx(df_ranked_ctrl, "Ranked_Protein_Abundance_CTRL_log10.xlsx", overwrite = TRUE)

#Proteine meiner Signatur, die hervorgehoben werden sollen 
# Liste der Proteine, die hervorgehoben werden sollen 
highlighted_proteins_group1 <- c("CRYM_mean", "CAPZA2_mean", "ALDH16A1_mean","PFKL_mean","SERPINC1_mean","HP_mean") #Gruppe 1 (rot)
highlighted_proteins_group2 <- c("LTF_mean","LCN1_mean","ALB_mean","LYZ_mean","IGHA1_mean") #Gruppe 2 (neon)
highlighted_proteins_group3 <- c("CHI3L2_mean") #Gruppe 3 (türkis)
# Sortieren der Proteine nach Abundance (aufsteigend) 
df_ranked_ctrl <- df_ranked_ctrl[order(df_ranked_ctrl$Abundance), ]
# Neue Spalte für die Farbmarkierung erstellen
# Neue Spalte für die Farbmarkierung und Gruppierung erstellen
df_ranked_ctrl$Highlight <- ifelse(df_ranked_ctrl$Protein %in% highlighted_proteins_group1, "group1", 
                                   ifelse(df_ranked_ctrl$Protein %in% highlighted_proteins_group2, "group2", 
                                          ifelse(df_ranked_ctrl$Protein %in% highlighted_proteins_group3, "group3", "default")))
#Graph ploten
library(ggplot2)
library(openxlsx)
library(ggrepel)
# Rang für jedes Protein erstellen
df_ranked_ctrl$Rank <- seq_len(nrow(df_ranked_ctrl))
# Ranked Protein Abundance Plot
ggplot(df_ranked_ctrl, aes(x = Rank, y = Log10_Abundance, color = Highlight, size = Highlight)) +
  geom_point(alpha = 0.6) +  # Punktgröße nach Highlight
  
  # Beschriftung für Gruppe 1 (rot)
  geom_text_repel(
    data = subset(df_ranked_ctrl, Highlight == "group1"),
    aes(label = Protein), 
    color = "red",  
    size = 4,  
    box.padding = 0.5,  
    point.padding = 0.2  
  ) +
  
  # Beschriftung für Gruppe 2 (neon)
  geom_text_repel(
    data = subset(df_ranked_ctrl, Highlight == "group2"),
    aes(label = Protein), 
    color = "limegreen",  
    size = 4,  
    box.padding = 0.5,  
    point.padding = 0.2  
  ) +
  
  # Beschriftung für Gruppe 3 (türkis)
  geom_text_repel(
    data = subset(df_ranked_ctrl, Highlight == "group3"),
    aes(label = Protein), 
    color = "turquoise",  
    size = 4,  
    box.padding = 0.5,  
    point.padding = 0.2  
  ) +
  
  theme_minimal() +
  labs(title = "Ranked Protein Abundance (CTRL)", 
       x = "Protein Rank", 
       y = "Log10 Abundance") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black"),   
    axis.ticks = element_line(color = "black"),  
    axis.ticks.length = unit(0.2, "cm"),          
    panel.grid = element_blank() 
  ) +
  
  scale_color_manual(values = c("group1" = "#FF0000",  
                                "group2" = "#32CD32",  
                                "group3" = "#40E0D0",  
                                "default" = "#555555"))  +
  
  scale_size_manual(values = c(
    "group1" = 4,
    "group2" = 4,
    "group3" = 4,
    "default" = 3
  ))  # Punktgrößen je nach Highlight


#Jetzt das ganze für ALS und Kontrollen
#Erstmal dataframe für ALS machen
#Lade Tabelle nicht normalisierte Einzelmittelwerte pro ALS sample
df_rpa_ALS <- read_excel(file.choose())
View(df_rpa_ALS)
# Berechnung des Mittelwerts pro Protein für die Kontrollen ab Spalte 6
mean_protein_ALS <- colMeans(df_ALS[, 6:ncol(df_ALS)], na.rm = TRUE)
View(mean_protein_ALS)
#Dataframe für ALS erstellen
# Erstelle ein DataFrame mit Spaltennamen und Mittelwerten
df_ranked_ALS <- data.frame(Protein = names(mean_protein_ALS), Abundance = mean_protein_ALS)
# Sortieren der Proteine nach Abundance (absteigend)
df_ranked_ALS <- df_ranked_ALS[order(-df_ranked_ALS$Abundance), ]
# Ergebnis anzeigen
View(df_ranked_ALS)
# Speichern als Excel
library(openxlsx)
write.xlsx(df_ranked_ALS, "Ranked_Protein_Abundance_ALS.xlsx", overwrite = TRUE)
# Log10-Transformation der Abundance-Werte
df_ranked_ALS$Log10_Abundance <- log10(df_ranked_ALS$Abundance)
# Ergebnis anzeigen
View(df_ranked_ALS)
# Optional: Speichern als Excel
library(openxlsx)
write.xlsx(df_ranked_ALS, "Ranked_Protein_Abundance_ALS_log10.xlsx", overwrite = TRUE)
#Proteine meiner Signatur, die hervorgehoben werden sollen 
# Liste der Proteine, die hervorgehoben werden sollen 
#Proteins, die gehighlighted werden sollen
highlighted_proteins_ctrl_group1 <- c("CRYM_mean", "CAPZA2_mean", "ALDH16A1_mean","PFKL_mean","SERPINC1_mean","HP_mean") #Gruppe 1 (rot)
highlighted_proteins_ctrl_group2 <- c("LTF_mean","LCN1_mean","ALB_mean","LYZ_mean","IGHA1_mean") #Gruppe 2 (neon)
highlighted_proteins_ctrl_group3 <- c("CHI3L2_mean") #Gruppe 3 (türkis)

highlighted_proteins_ALS_group1 <- c("CRYM_mean", "CAPZA2_mean", "ALDH16A1_mean","PFKL_mean","SERPINC1_mean","HP_mean") #Gruppe 1 (rot)
highlighted_proteins_ALS_group2 <- c("LTF_mean","LCN1_mean","ALB_mean","LYZ_mean","IGHA1_mean") #Gruppe 2 (neon)
highlighted_proteins_ALS_group3 <- c("CHI3L2_mean") #Gruppe 3 (türkis)
df_combined$HighlightGroup <- ifelse(df_combined$Protein %in% highlighted_proteins_ctrl_group1, "group1_ctrl", 
                                     ifelse(df_combined$Protein %in% highlighted_proteins_ctrl_group2, "group2_ctrl", 
                                            ifelse(df_combined$Protein %in% highlighted_proteins_ctrl_group3, "group3_ctrl", 
                                                   ifelse(df_combined$Protein %in% highlighted_proteins_ALS_group1, "group1_als", 
                                                          ifelse(df_combined$Protein %in% highlighted_proteins_ALS_group2, "group2_als", 
                                                                 ifelse(df_combined$Protein %in% highlighted_proteins_ALS_group3, "group3_als", "default"))))))
# Sortieren der Proteine nach Abundance (aufsteigend) 
df_ranked_ALS <- df_ranked_ALS[order(df_ranked_ALS$Abundance), ]
df_ranked_ALS$Rank <- seq_len(nrow(df_ranked_ctrl))
View(df_ranked_ALS)
# Neue Spalte für die Farbmarkierung und Gruppierung erstellen
df_ranked_ALS$Highlight <- ifelse(df_ranked_ALS$Protein %in% highlighted_proteins_group1, "group1", 
                                   ifelse(df_ranked_ALS$Protein %in% highlighted_proteins_group2, "group2", 
                                          ifelse(df_ranked_ALS$Protein %in% highlighted_proteins_group3, "group3", "default")))
# Verknüpfen der beiden DataFrames (CTRL und ALS)
df_ranked_ctrl$Group <- "CTRL"  # Gruppe für die CTRL-Daten
df_ranked_ALS$Group <- "ALS"    # Gruppe für die ALS-Daten
# Kombiniere beide DataFrames (d.h. df_ranked_ctrl und df_ranked_ALS)
df_combined <- rbind(df_ranked_ctrl, df_ranked_ALS)
#Wegen Fehler: # Überprüfen der Spaltennamen von df_ranked_ctrl und df_ranked_ALS
colnames(df_ranked_ctrl)
colnames(df_ranked_ALS)
# Kombiniere beide DataFrames (d.h. df_ranked_ctrl und df_ranked_ALS)
df_combined <- rbind(df_ranked_ctrl, df_ranked_ALS)
View(df_combined)
#Speichern als Excel
write.xlsx(df_combined, "Ranked_Protein_Abundance_ALS+CTRL.xlsx", overwrite = TRUE)
# Rang für jedes Protein in beiden Gruppen erstellen
df_combined$Rank <- rep(seq_len(nrow(df_combined[df_combined$Group == "CTRL", ])), 2)


#Plotte 
# Ranked Protein Abundance Plot für ALS vs. CTRL
# Neuen Plot für df_combined erstellen
# Plot erstellen
ggplot(df_combined, aes(x = Rank, y = Log10_Abundance, color = Highlight)) +
  geom_point(alpha = 0.5) +  # Punkte für jedes Protein
  
  
  # Beschriftung für Gruppe 1 CTRL (grau)
  geom_text_repel(
    data = subset(df_combined, Highlight == "group1_ctrl"),  # Nur Gruppe 1 in CTRL beschriften
    aes(label = Protein), 
    color = "#7F7F7F",  # Grau für Gruppe 1
    size = 4,  # Schriftgröße der Labels
    box.padding = 0.5,  # Abstand der Labels von den Punkten
    point.padding = 0.2  # Abstand zwischen Punkt und Label
  ) +
  
  # Beschriftung für Gruppe 2 CTRL (blau)
  geom_text_repel(
    data = subset(df_combined, Highlight == "group2_ctrl"),  # Nur Gruppe 2 in CTRL beschriften
    aes(label = Protein), 
    color = "blue",  # Blau für Gruppe 2
    size = 4,  # Schriftgröße der Labels
    box.padding = 0.5,  # Abstand der Labels von den Punkten
    point.padding = 0.2  # Abstand zwischen Punkt und Label
  ) +
  
  # Beschriftung für Gruppe 3 CTRL (türkis)
  geom_text_repel(
    data = subset(df_combined, Highlight == "group3_ctrl"),  # Nur Gruppe 3 in CTRL beschriften
    aes(label = Protein), 
    color = "#40E0D0",  # Türkis für Gruppe 3
    size = 4,  # Schriftgröße der Labels
    box.padding = 0.5,  # Abstand der Labels von den Punkten
    point.padding = 0.2  # Abstand zwischen Punkt und Label
  ) +
  
  # Beschriftung für Gruppe 1 ALS (rot)
  geom_text_repel(
    data = subset(df_combined, Highlight == "group1_ALS"),  # Nur Gruppe 1 in ALS beschriften
    aes(label = Protein), 
    color = "#FF4500",  # Orangerot für Gruppe 1 in ALS
    size = 4,  # Schriftgröße der Labels
    box.padding = 0.5,  # Abstand der Labels von den Punkten
    point.padding = 0.2  # Abstand zwischen Punkt und Label
  ) +
  
  # Beschriftung für Gruppe 2 ALS (neon-grün)
  geom_text_repel(
    data = subset(df_combined, Highlight == "group2_ALS"),  # Nur Gruppe 2 in ALS beschriften
    aes(label = Protein), 
    color = "#32CD32",  # Neon-Grün für Gruppe 2 in ALS
    size = 4,  # Schriftgröße der Labels
    box.padding = 0.5,  # Abstand der Labels von den Punkten
    point.padding = 0.2  # Abstand zwischen Punkt und Label
  ) +
  
  # Beschriftung für Gruppe 3 ALS (türkis)
  geom_text_repel(
    data = subset(df_combined, Highlight == "group3_ALS"),  # Nur Gruppe 3 in ALS beschriften
    aes(label = Protein), 
    color = "#40E0D0",  # Türkis für Gruppe 3 in ALS
    size = 4,  # Schriftgröße der Labels
    box.padding = 0.5,  # Abstand der Labels von den Punkten
    point.padding = 0.2  # Abstand zwischen Punkt und Label
  ) +
  
  theme_minimal() +
  labs(title = "Ranked Protein Abundance (ALS vs CTRL)", 
       x = "Protein Rank", 
       y = "Log10 Abundance") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black"),   # X- und Y-Achse schwarz
    axis.ticks = element_line(color = "black"),  # Achsen-Ticks schwarz
    axis.ticks.length = unit(0.2, "cm"),         # Länge der Ticks anpassen
    panel.grid = element_blank()                 # Entfernt das Gitter im Hintergrund
  ) +
  scale_color_manual(values = c("group1_ctrl" = "#FF0000",  # Rot für Gruppe 1 in CTRL
                                "group2_ctrl" = "#32CD32",  # Neon-Grün für Gruppe 2 in CTRL
                                "group3_ctrl" = "#40E0D0",  # Türkis für Gruppe 3 in CTRL
                                "group1_ALS" = "#FF4500",   # Orangerot für Gruppe 1 in ALS
                                "group2_ALS" = "#32CD32",   # Neon-Grün für Gruppe 2 in ALS
                                "group3_ALS" = "#40E0D0",   # Türkis für Gruppe 3 in ALS
                                "default" = "#555555"))     # Grau für alle anderen






#Nur für ALS 
#Scatterplot for ranked protein abundance in ALS
#Erstmal dataframe für ALS machen
#Lade Tabelle nicht normalisierte Einzelmittelwerte pro ALS sample
library(openxlsx)
library(readxl)
df_rpa_ALS <- read_excel(file.choose())
View(df_rpa_ALS)
# Berechnung des Mittelwerts pro Protein für die Kontrollen ab Spalte 6
mean_protein_ALS <- colMeans(df_rpa_ALS[, 6:ncol(df_rpa_ALS)], na.rm = TRUE)
View(mean_protein_ALS)
#Dataframe für ALS erstellen
# Erstelle ein DataFrame mit Spaltennamen und Mittelwerten
df_ranked_ALS <- data.frame(Protein = names(mean_protein_ALS), Abundance = mean_protein_ALS)
# Sortieren der Proteine nach Abundance (absteigend)
df_ranked_ALS <- df_ranked_ALS[order(-df_ranked_ALS$Abundance), ]
# Ergebnis anzeigen
View(df_ranked_ALS)
# Speichern als Excel
library(openxlsx)
write.xlsx(df_ranked_ALS, "Ranked_Protein_Abundance_ALS.xlsx", overwrite = TRUE)
# Log10-Transformation der Abundance-Werte
df_ranked_ALS$Log10_Abundance <- log10(df_ranked_ALS$Abundance)
# Ergebnis anzeigen
View(df_ranked_ALS)
# Optional: Speichern als Excel
write.xlsx(df_ranked_ALS, "Ranked_Protein_Abundance_ALS_log10.xlsx", overwrite = TRUE)
#Proteine meiner Signatur, die hervorgehoben werden sollen 
# Liste der Proteine, die hervorgehoben werden sollen 
highlighted_proteins_group1 <- c("CRYM_mean", "CAPZA2_mean", "ALDH16A1_mean","PFKL_mean","SERPINC1_mean","HP_mean") #Gruppe 1 (rot)
highlighted_proteins_group2 <- c("LTF_mean","LCN1_mean","ALB_mean","LYZ_mean","IGHA1_mean") #Gruppe 2 (neon)
highlighted_proteins_group3 <- c("CHI3L2_mean") #Gruppe 3 (türkis)
# Sortieren der Proteine nach Abundance (aufsteigend) 
df_ranked_ALS <- df_ranked_ALS[order(df_ranked_ALS$Abundance), ]
# Neue Spalte für die Farbmarkierung erstellen
# Neue Spalte für die Farbmarkierung und Gruppierung erstellen
df_ranked_ALS$Highlight <- ifelse(df_ranked_ALS$Protein %in% highlighted_proteins_group1, "group1", 
                                   ifelse(df_ranked_ALS$Protein %in% highlighted_proteins_group2, "group2", 
                                          ifelse(df_ranked_ALS$Protein %in% highlighted_proteins_group3, "group3", "default")))
#Graph ploten
library(ggplot2)
library(openxlsx)
library(ggrepel)
# Rang für jedes Protein erstellen
df_ranked_ALS$Rank <- seq_len(nrow(df_ranked_ALS))
#Plot des Graphs
ggplot(df_ranked_ALS, aes(x = Rank, y = Log10_Abundance, color = Highlight, size = Highlight)) +
  geom_point(alpha = 0.6) +  # Point size is mapped to Highlight
  # Replace `size` with `linewidth` for lines
  geom_text_repel(
    data = subset(df_ranked_ALS, Highlight == "group1"),
    aes(label = Protein),
    color = "red",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2
  ) +
  geom_text_repel(
    data = subset(df_ranked_ALS, Highlight == "group2"),
    aes(label = Protein),
    color = "limegreen",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2
  ) +
  geom_text_repel(
    data = subset(df_ranked_ALS, Highlight == "group3"),
    aes(label = Protein),
    color = "turquoise",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2
  ) +
  theme_minimal() +
  labs(title = "Ranked Protein Abundance (ALS)", 
       x = "Protein Rank", 
       y = "Log10 Abundance") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = c("group1" = "#FF0000",  
                                "group2" = "#32CD32", 
                                "group3" = "#40E0D0",  
                                "default" = "#555555")) +
  scale_size_manual(values = c(
    "group1" = 4,
    "group2" = 4,
    "group3" = 4,
    "default" = 3
  ))


#Differentielle Expression Mean ALS/Mean CTRL 
#Lade Tabelle mit Mittelwerte
library(openxlsx)
df_rpa_ALS_CTRL <- read_excel(file.choose())
View(df_rpa_ALS_CTRL)
#Benennen der relevanten Zeilen 
colnames(df_rpa_ALS_CTRL)[c(1,3)] <- c("Protein", "Log10_Abundance")
#Nach Log10-Abundanz sortieren (absteigend)
df_ranked_ALS_CTRL <- df_rpa_ALS_CTRL[order(-df_rpa_ALS_CTRL$Log10_Abundance), ]
#Rank vergeben
df_ranked_ALS_CTRL$Rank <- seq_len(nrow(df_ranked_ALS_CTRL))
df_ranked_ALS_CTRL$Rank <- rank(df_ranked_ALS_CTRL$Log10_Abundance, ties.method = "first")

View(df_ranked_ALS_CTRL)
#Als Excel Speichern
write.xlsx(df_ranked_ALS_CTRL, "Ranked_Protein_Abundance_ALS_CTRL_log10.xlsx", overwrite = TRUE)
# Bibliotheken laden
library(ggplot2)
library(ggrepel)

# Liste der spezifisch hervorzuhebenden Proteine
highlighted_proteins <- c("CRYM", "CAPZA2", "ALDH16A1", "PFKL", "SERPINC1", "HP")

# Neue Spalte für die Farbmarkierung erstellen
df_ranked_ALS_CTRL$Highlight <- ifelse(df_ranked_ALS_CTRL$Protein %in% highlighted_proteins, "highlighted",
                                       ifelse(df_ranked_ALS_CTRL$Log10_Abundance > 0, "positive",
                                              ifelse(df_ranked_ALS_CTRL$Log10_Abundance < 0, "negative", "neutral")))

# Top 10 Proteine mit den höchsten positiven Log10_Abundance-Werten
top10_positive <- df_ranked_ALS_CTRL %>%
  filter(Log10_Abundance > 0) %>%
  slice_max(Log10_Abundance, n = 10)


# Top 10 Proteine mit den niedrigsten negativen Log10_Abundance-Werten
top10_negative <- df_ranked_ALS_CTRL %>%
  filter(Log10_Abundance < 0) %>%
  slice_min(Log10_Abundance, n = 10)

# Neue Spalte für die Farbmarkierung: Standardmäßig grau
df_ranked_ALS_CTRL$Highlight <- ifelse(df_ranked_ALS_CTRL$Protein %in% highlighted_proteins, "highlighted", "neutral")
# Top 10 positive (rot) und Top 10 negative (blau) einfärben
df_ranked_ALS_CTRL$Highlight[df_ranked_ALS_CTRL$Protein %in% top10_positive$Protein] <- "positive"
df_ranked_ALS_CTRL$Highlight[df_ranked_ALS_CTRL$Protein %in% top10_negative$Protein] <- "negative"
# Plot erstellen
library(ggplot2)
library(ggrepel)
library(dplyr)

ggplot(df_ranked_ALS_CTRL, aes(x = Rank, y = Log10_Abundance, color = Highlight)) +
  geom_point(aes(size = Highlight), alpha = 0.6) +  # Punktgröße nach Highlight anpassen
 
  
  # Beschriftung für spezifisch hervorgehobene Proteine
  geom_text_repel(
    data = subset(df_ranked_ALS_CTRL, Highlight == "highlighted"),
    aes(label = Protein),
    color = "green",
    size = 4,  # Fehlende Größe ergänzt
    box.padding = 0.5,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  # Beschriftung der Top 10 positiven Proteine
  geom_text_repel(
    data = top10_positive,
    aes(label = Protein),
    color = "red",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  # Beschriftung der Top 10 negativen Proteine
  geom_text_repel(
    data = top10_negative,
    aes(label = Protein),
    color = "blue",
    size = 4,
    box.padding = 0.5,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  
  theme_minimal() +
  labs(title = "Ranked Protein Abundance (CTRL + ALS)",
       x = "Overall Abundance",
       y = "Log10 Abundance") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid = element_blank()
  ) +
  scale_color_manual(values = c(
    "highlighted" = "green",
    "positive" = "red",
    "negative" = "blue",
    "neutral" = "#555555"
  )) +
  
  scale_size_manual(values = c(
    "highlighted" = 4,
    "positive" = 4,
    "negative" = 4,
    "neutral" = 3
  )) 
