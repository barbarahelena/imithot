## Data cleaning IMITHOT
## b.j.verhaar@amsterdamumc.nl

### Libraries
library(tidyverse)

### Script
#### Opening clinical data
df <- readRDS('data/IMITHOT_clinical_long.RDS')
dim(df)
head(df)

df <- df |> mutate(
  Group = as.factor(Group),
  Timepoint = case_when(Timepoint == "2" ~ 1,
                    Timepoint == "4" ~ 2,
                    Timepoint == "5" ~ 3,
                    Timepoint == "6" ~ 4,
                    Timepoint == "7" ~ 5
                  ),
  Timepoint = as.factor(Timepoint),
  FUtime = case_when(Timepoint == "1" ~ 0,
                    Timepoint == "2" ~ 3,
                    Timepoint == "3" ~ 6,
                    Timepoint == "4" ~ 12,
                    Timepoint == "5" ~ 24
                  ),
  ID = str_c(Study, "_", ID),
  Hashimoto = case_when(
    Group %in% c("Allogenic", "Autologous") ~ "Hashimoto",
    Group == "Healthy" ~ "Control"
  ),
  Hashimoto = as.factor(Hashimoto),,
  Donor = case_when(str_detect(ID, "1100") ~ TRUE, .default = FALSE),
) |>
  mutate(FUtime = as.numeric(FUtime)) |> 
  select(-Study)
head(df)

# IMITHOT_ALT023 should have been IMITHOT_023
df$ID[which(df$ID == "IMITHOT_ALT023")] <- "IMITHOT_023"
df$SampleID[which(df$SampleID == "IMITHOT_ALT023_2")] <- "IMITHOT_023_2"
df$Group[which(df$SampleID == "IMITHOT_023_2")] <- "Allogenic"
df$Hashimoto[which(df$SampleID == "IMITHOT_023_2")] <- "Control"
df$Donor[which(df$SampleID == "IMITHOT_023_2")] <- FALSE

# Make HC001_4 a baseline sample
df$SampleID[which(df$SampleID == "IMITHOT_HC001_4")] <- "IMITHOT_HC001_2"
df$Timepoint[which(df$SampleID == "IMITHOT_HC001_2")] <- "2"
df$FUtime[which(df$SampleID == "IMITHOT_HC001_2")] <- 0

# Add IMITHOT_008_6 and IMITHOT_014_6
df <- df %>% add_row(ID = "IMITHOT_008", SampleID = "IMITHOT_008_6", Hashimoto = "Hashimoto",
                          Group = "Allogenic", FUtime = 12, Timepoint = "4", Donor = FALSE) |>
              add_row(ID = "IMITHOT_014", SampleID = "IMITHOT_014_6", Hashimoto = "Hashimoto",
                          Group = "Autologous", FUtime = 12, Timepoint = "4", Donor = FALSE)

# Remove follow-up healthy samples, these samples were collected under AB Tx
hc_antibiotics <- df |> filter(Group == "Healthy") |> filter(FUtime == 12) |> pull(SampleID)
df <- df |> filter(!SampleID %in% hc_antibiotics)

head(df)
summary(df$Group)
summary(df$Timepoint)
saveRDS(df, "data/imithot_meta.RDS")
