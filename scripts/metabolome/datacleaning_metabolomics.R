## Data cleaning metabolomics IMITHOT
## b.j.verhaar@amsterdamumc.nl

### Libraries
library(tidyverse)

### Script
#### Opening metabolomics data, batch-corrected
df <- readxl::read_xlsx("data/metabolomics_imithot_table.xlsx")
head(df)
dfmeta <- readxl::read_xlsx("data/metabolomics_imithot_sampledata.xlsx")
head(dfmeta)
info <- readxl::read_xlsx("data/metabolomics_imithot_info.xlsx")
head(info)
meta <- readRDS("data/imithot_meta.RDS")

# Fix IDs
df <- df |> 
  mutate(
    SampleID2 = dfmeta$CLIENT_IDENTIFIER[match(dfmeta$SampleID, SampleID)],
    SampleID = str_replace_all(str_remove(SampleID2, "-\\d{2}/\\d{2}/\\d{4}$"), "-", "_"),
    SampleID = str_remove(str_replace(SampleID, "IMI", "IMITHOT"), "V"),
    SampleID = str_replace(SampleID, "HC", "IMITHOT_HC"),
    SampleID = str_replace(SampleID, "IMITHOT016", "IMITHOT_016"),
    visitdate = str_extract(SampleID2, "\\d{2}/\\d{2}/\\d{4}$"),
    visitdate = dmy(visitdate),
  ) |> 
  select(-SampleID2) |> 
  select(SampleID, visitdate, everything()) |> 
  filter(!is.na(visitdate)) # all non-imithot samples have no date
head(df)
names(df)
dim(df)

df$SampleID

# Fix metabolite names
colnames(df)[which(str_detect(colnames(df), "X"))] <- info$PLOT_NAME[match(str_c("X", info$MetID), colnames(df)[which(str_detect(colnames(df), "X"))])]
names(df)

# Remove all metabolites with more than 10% zero
countna <- function(x) sum(x == 0)
imp_list <- sapply(df, countna)
sum(imp_list > 10) # zero for more than 10 samples of 103 total (>10%) -> 470 metabolites
imp <- names(imp_list[imp_list > 10])
df <- df |> select(-any_of(imp))
df <- as.data.frame(df)
dim(df) # 1088 metabolites left

# Imputation with half the lowest value
df_imputed <- df %>% mutate(across(where(is.numeric), ~ ifelse(is.na(.x) | .x == 0,
                         min(.x[.x > 0], na.rm = TRUE) / 2, .x)))

# Normalize and transform
df_norm <- df_imputed %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x, na.rm = TRUE))) |> 
  mutate(across(where(is.numeric), ~ .x * median(rowSums(df_imputed |> select(where(is.numeric)), na.rm = TRUE))))
df_log <- df_norm %>% mutate(across(where(is.numeric), ~ log10(.x)))
df_scaled <- df_log %>% mutate(across(where(is.numeric), scale))

# Remove all unnamed and xenobiotic metabolites
mets_X <- df |> select(starts_with("X")) |> names()
mets_xeno <- info |> filter(SUPER_PATHWAY == "Xenobiotics") |> pull(PLOT_NAME)
removemet <- c(mets_X, mets_xeno)
df_scaled_filtered <- df_scaled |> select(-any_of(removemet))
dim(df_scaled_filtered)
head(df_scaled_filtered)[1:5,1:5]
df_scaled_filtered <- df_scaled_filtered |> inner_join(meta)
df_scaled_filtered$FUtime

df_scaled_filtered$SampleID
tail(df_scaled_filtered)[,1:7]
df_scaled_filtered[,1:7]

# Checks
hist(df_imputed$`S-1-pyrroline-5-carboxylate`)
hist(df_log$`S-1-pyrroline-5-carboxylate`)
hist(df_scaled_filtered$`S-1-pyrroline-5-carboxylate`)

# Save file
saveRDS(df_scaled_filtered, "data/metabolomics_clean.RDS")
