## Data cleaning microbiome IMITHOT
## b.j.verhaar@amsterdamumc.nl

## Load libraries
library(tidyverse)
library(rio)
library(forcats)

# Function for plots
theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.5, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

# Function to add Group to the missing samples
add_metadata <- function(mbdata, metadata, metadata_cols) {
  # Filter the original data to get the metadata columns
  metadata <- meta %>%
    filter(Timepoint == "1") %>%
    select(ID, all_of(metadata_cols)) %>%
    filter(!duplicated(ID))
  # Merge the metadata with the microbiome data, redefine metadata columns
  mbdata_with_metadat <- mbdata %>% select(-Group) %>%
    left_join(., metadata, by = "ID")
  return(mbdata_with_metadat)
}

impute_missing <- function(x) {
  n <- length(x)
  if (n < 3) return(x)  # If there are less than 3 elements, return as is
  
  for (i in 2:(n - 1)) {
    if (is.na(x[i]) & !is.na(x[i - 1]) & !is.na(x[i + 1])) {
      x[i] <- mean(c(x[i - 1], x[i + 1]), na.rm = TRUE)
    } else if (is.na(x[i]) & !is.na(x[i - 1])){
      x[i] <- x[i - 1] + (x[i - 1] - x[i - 2])
      if(x[i] < 0) x[i] <- x[i - 1]
    }
  }
  return(x)
}

# Metadata and shotgun
meta <- readRDS("data/imithot_meta.RDS")
sg <- readRDS("data/shotgun_abundance.RDS")
rownames(sg)[which(rownames(sg) == "IMITHOT_ALT023_2")] <- "IMITHOT_023_2"
meta2 <- meta |> filter(SampleID %in% rownames(sg))
head(sg)[1:5,1:5]
dim(sg)
rownames(sg)

# Pruning of species - IMITHOT patients
imi <- meta2 |> filter(Group != "Healthy") |> pull(SampleID)
sg_imithot <- sg[which(rownames(sg) %in% imi),]
threshold <- 0.001 # 0.1% threshold
min_samples <- 0.3 * nrow(sg_imithot) # in 10% of samples
species_filter <- apply(sg_imithot, 2, function(x) sum(x >= threshold) >= min_samples)
sg_filt1 <- sg_imithot[, species_filter]
dim(sg_filt1) # 115 species and 81 samples

saveRDS(sg_filt1, "data/microbiome_filtered_pt.RDS")
write.csv(sg_filt1, "data/microbiome_filtered_pt.csv")

# Pruning of species - IMITHOT patients and healthy controls
threshold <- 0.001 # 0.1% threshold
min_samples <- 0.3 * nrow(sg_imithot) # in 10% of samples
species_filter <- apply(sg, 2, function(x) sum(x >= threshold) >= min_samples)
sg_filt2 <- sg[, species_filter]
dim(sg_filt2) # 135 species and 111 samples

saveRDS(sg_filt2, "data/microbiome_filtered_pt_hc.RDS")
write.csv(sg_filt2, "data/microbiome_filtered_pt_hc.csv")

# Plot a bar plot of available samples
## In total per timepoint per group
sample_counts <- meta2 %>%
  group_by(Group, FUtime) %>%
  summarise(count = n()) %>%
  ungroup()
sample_counts
ggplot(sample_counts, aes(x = FUtime, y = count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of samples",
       x = "Follow-up time (months)",
       y = "Number of samples") +
  ggsci::scale_fill_nejm() +
  scale_x_continuous(breaks = c(0, 3, 6, 12, 24)) +
  theme_Publication() +
  theme(legend.title = element_blank())
ggsave("results/microbiome/sample_counts_barplot.pdf", width = 10, height = 6)

## Per patient
sample_counts <- meta %>%
  group_by(ID, FUtime) %>%
  summarise(count = n(), Group = Group) %>%
  ungroup()
ggplot(sample_counts, aes(x = as.factor(FUtime), y = count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of samples",
       x = "Follow-up time (months)",
       y = "Number of samples") +
  ggsci::scale_fill_nejm() +
  # scale_x_continuous(breaks = c(0, 3, 6, 12, 24)) +
  facet_wrap(~ID) +
  theme_Publication() +
  theme(legend.title = element_blank())
ggsave("results/microbiome/sample_counts_barplot_id.pdf", width = 15, height = 15)

## Imputing data for TCAM - only IMITHOT patients
## First make a complete dataset: microbiome and metadata
sg_filt1 <- as.data.frame(sg_filt1)
sg_filt1$SampleID <- rownames(sg_filt1)
meta3 <- meta2 |> filter(Group != "Healthy") |> droplevels()
all_levels <- levels(meta3$Timepoint)
all_levels
table(meta3$FUtime)
dftot <- left_join(meta3, sg_filt1, by = "SampleID") # merge metadata with abundance table

## Add missing timepoints to table and fill in metadata for these samples
metadata_cols <- c("Group") # to add to missing samples below
dftot_compl <- dftot %>% 
  group_by(ID) %>% 
  complete(Timepoint = all_levels) %>% # complete levels to make missing timepoints explicit
  add_metadata(dftot_compl2, metadata_cols) %>% # see function at the top
  mutate(FUtime = case_when(Timepoint == "1" ~ 0,
                    Timepoint == "2" ~ 3,
                    Timepoint == "3" ~ 6,
                    Timepoint == "4" ~ 12,
                    Timepoint == "5" ~ 24
                  )) %>% # redefine for missing samples
  arrange(FUtime) %>%
  mutate(Timepoint = fct_inorder(as.factor(Timepoint)), # so that it sorts nicely later
         ID = fct_inorder(as.factor(ID))) |> 
  droplevels()

microbiome_cols <- names(dftot_compl)[c(5:(ncol(dftot_compl)-1))] # last three cols are metadata
head(microbiome_cols); tail(microbiome_cols)

# Pivot to long format and imputation
dftot_long <- dftot_compl |> 
  pivot_longer(cols = all_of(microbiome_cols), names_to = "microbe", values_to = "abundance") |> 
  arrange(ID, FUtime, microbe) |> 
  mutate(Timepoint = fct_inorder(Timepoint), 
         ID = fct_inorder(as.factor(ID)),ID) |> 
  group_by(ID, microbe) |> 
  arrange(Timepoint) %>%
  mutate(abundance = impute_missing(abundance)) %>%
  ungroup()

# And pivot back to wide format
df_imputed_wide <- dftot_long %>%
  pivot_wider(names_from = microbe, values_from = abundance) %>%
  select(SampleID, ID, Timepoint, FUtime, Group, all_of(microbiome_cols))

dim(df_imputed_wide) # 84 samples (2 samples added, by filling in the missing timepoint at 12 months)
table(df_imputed_wide$Timepoint) # now each timepoint has 21 samples
table(df_imputed_wide$Group, df_imputed_wide$FUtime) # 9 autologous and 12 allogenic

## Per patient
sample_counts <- df_imputed_wide %>%
  filter(Timepoint != "2") |> 
  group_by(ID, FUtime) %>%
  summarise(count = n(), Group = Group) %>%
  ungroup()
ggplot(sample_counts, aes(x = as.factor(FUtime), y = count, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of samples",
       x = "Follow-up time (months)",
       y = "Number of samples") +
  ggsci::scale_fill_nejm() +
  # scale_x_continuous(breaks = c(0, 3, 6, 12, 24)) +
  facet_wrap(~ID) +
  theme_Publication() +
  theme(legend.title = element_blank())
ggsave("results/microbiome/imp_sample_counts_barplot_id.pdf", width = 15, height = 15)

## Save the filtered data (no mice with > 2 missing timepoints, ignoring wk 18)
saveRDS(df_imputed_wide, "data/filtered_imputed_microbiome_intervention.RDS")
write.csv(df_imputed_wide, "data/filtered_imputed_microbiome_intervention.csv", row.names = FALSE)
