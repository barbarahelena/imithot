# ML model looking at differences Hashimoto-healthy: XGBoost input

library(tidyverse)
library(ggpubr)
rm(list=ls())

# make data for machine learning XGB classification models

# writes input data files for XGB models as tab-delimited 
# subject ids and feature ids are written as separate tab-delimited files
# write X data / predictors
write_data <- function(x, data_path){
    x <- as.matrix(x)
    if(any(is.na(x))){
        cat('There are missing values in the input data!\n')
    }
    write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write y / predicted outcome
write_y <- function(x, name_y, data_path){
    if(missing(name_y)){
        cat('\n\nYou need to provide a name for the y data file!\n')
    }
    if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
        cat('\nThe file name is not compatible with XGBeast!\n' )
    }
    if(any(is.na(x))){
        cat('\nThere are missing values in the outcome data!\n')
    }
    write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}

## Open dataframe
df <- readRDS('data/imithot_meta.RDS')
df <- df |> filter(FUtime == 0) |> 
        mutate(Hashimoto = case_when(Hashimoto == "Hashimoto" ~ 1, Hashimoto == "Control" ~ 0))
mb <- readRDS('data/microbiome_filtered_pt_hc.RDS')
mb <- mb[rownames(mb) %in% df$SampleID, ]
mb <- as.data.frame(mb)
dim(mb)
dim(df)

clindf <- df[match(rownames(mb), df$SampleID), ] # put IDs in order of mb
all(clindf$SampleID == rownames(mb)) # TRUE
clindf$SampleID; rownames(mb) # check if these are the same

path <- 'hashimoto'
dir.create(path, showWarnings = FALSE); dir.create("hashimoto/input_data", showWarnings = FALSE)
head(mb) # check if IDs are rownames, data present and numerics, colnames are species
write_data(mb, file.path(path, 'input_data'))
y <- as.data.frame(clindf$Hashimoto)
y # classification model - check if these are 0 and 1
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
