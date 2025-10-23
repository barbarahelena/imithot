# Beta diversity plots IMITHOT
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Librarieslibrary(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ape)
library(vegan)

theme_Publication <- function(base_size=14, base_family="sans") {
        library(grid)
        library(ggthemes)
        library(stringr)
        (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
        size = rel(1.0), hjust = 0.5),
        text = element_text(),
        panel.background = element_rect(colour = NA, fill = NA),
        plot.background = element_rect(colour = NA, fill = NA),
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
        legend.spacing = unit(0, "cm"),
        # legend.title = element_text(face="italic"),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        strip.text = element_text(face="bold"),
        plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
}

## Load data
mb <- readRDS("data/shotgun_cleaned.RDS")
df <- readRDS("data/imithot_meta.RDS")

#### Bray-Curtis distance ####
df_intervention <- df |> filter(Group %in% c("Autologous", "Allogenic"))
mb2 <- mb[rownames(mb) %in% df_intervention$SampleID,]
head(mb2)[1:5,1:5]
rowSums(mb2)

bray <- vegan::vegdist(mb2, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$SampleID <- rownames(dbray)
dbray <- left_join(dbray, df_intervention, by = 'SampleID') # add metadata / covariates

#### Bray-Curtis per time point ####
braypertimepoint <- function(timepoint, df = dbray, tab = mb2) {
    set.seed(14)
    dfsel <- df %>% filter(FUtime == timepoint)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$SampleID,], method = 'bray')
    return(adonis2(bray ~ Group, data = dfsel))
}

braypertimepoint(timepoint = 0, tab = mb2)
baseline <- c(braypertimepoint(timepoint = 0)['Pr(>F)'][[1]][1], "0")
halfyear <- c(braypertimepoint(timepoint = 6)['Pr(>F)'][[1]][1], "6")
fullyear <- c(braypertimepoint(timepoint = 12)['Pr(>F)'][[1]][1], "12")
twoyear <- c(braypertimepoint(timepoint = 24)['Pr(>F)'][[1]][1], "24")

res <- rbind(baseline, halfyear, fullyear, twoyear)
colnames(res) <- c("pvalue", "FUtime")
res <- as.data.frame(res)
res$FUtime_fct <- factor(res$FUtime, levels = c("0", "6", "12", "24"))
res$FUtime_fct <- fct_recode(res$FUtime, "0 months" = "0", "6 months" = "6", 
        "12 months" = "12", "24 months" = "24")
res$pvalue <- as.numeric(res$pvalue)
dbray$FUtime_fct <- fct_recode(as.character(dbray$FUtime), 
    "0 months" = "0", "6 months" = "6", "12 months" = "12", "24 months" = "24")
dbray$FUtime_fct <- fct_relevel(dbray$FUtime_fct, "6 months", after = 1L)

## Plots
(braycurt <- dbray %>% 
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Group), size = 3, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(8)[c(8,6)]) +
                    scale_fill_manual(values = pal_nejm()(8)[c(8,6)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = Group, fill = Group), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("PERMANOVA p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
                    facet_wrap(~ FUtime_fct))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_group.pdf", width = 8, height = 8)


(braycurt <- dbray %>% 
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = FUtime_fct), size = 3, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(8)[c(3:8)]) +
                    scale_fill_manual(values = pal_nejm()(8)[c(3:8)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = FUtime_fct, 
                        fill = FUtime_fct), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    # geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("PERMANOVA p = ", round(pvalue, 3))),
                    #           hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
                    facet_wrap(~ Group))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_time.pdf", width = 8, height = 8)


# Healthy vs Hashimoto
df_baseline <- df |> filter(FUtime == 0)
mb3 <- mb[rownames(mb) %in% df_baseline$SampleID,]
head(mb3)[1:5,1:5]

bray <- vegan::vegdist(mb3, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Relative_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$SampleID <- rownames(dbray)
dbray <- left_join(dbray, df_baseline, by = 'SampleID') # add metadata / covariates

set.seed(14)
dfsel <- dbray %>% filter(FUtime == 0)
bray <- vegan::vegdist(mb3[rownames(mb3) %in% dfsel$SampleID,], method = 'bray')
hashimoto <- adonis2(bray ~ Hashimoto, data = dfsel)
hashimotores <- hashimoto['Pr(>F)'][[1]][1]
hashimotores <- as.data.frame(hashimotores)
colnames(hashimotores) <- c("pvalue")

(braycurt <- dbray %>%
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Hashimoto), size = 3, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(8)[c(7,6)]) +
                    scale_fill_manual(values = pal_nejm()(8)[c(7,6)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = Hashimoto, 
                        fill = Hashimoto), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    geom_text(data = hashimotores, aes(x = Inf, y = Inf, label = paste0("PERMANOVA p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_hashimoto.pdf", width = 5, height = 5)

# Donors vs Hashimoto
df_donors <- df |> filter((!(Group == "Healthy" & Donor == FALSE)) & Group != "Autologous")
mb4 <- mb[rownames(mb) %in% df_donors$SampleID,]

bray <- vegan::vegdist(mb4, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$SampleID <- rownames(dbray)
dbray <- left_join(dbray, df_donors, by = 'SampleID') # add metadata / covariates

set.seed(14)
braypertimepoint <- function(timepoint, df = dbray, tab = mb4) {
    set.seed(14)
    dfsel <- df |> filter(Donor == TRUE | (FUtime == timepoint))
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$SampleID,], method = 'bray')
    return(adonis2(bray ~ Donor, data = dfsel))
}

braypertimepoint(timepoint = 0)
baseline <- c(braypertimepoint(timepoint = 0)['Pr(>F)'][[1]][1], "0")
halfyear <- c(braypertimepoint(timepoint = 6)['Pr(>F)'][[1]][1], "6")
fullyear <- c(braypertimepoint(timepoint = 12)['Pr(>F)'][[1]][1], "12")
twoyear <- c(braypertimepoint(timepoint = 24)['Pr(>F)'][[1]][1], "24")

res <- rbind(baseline, halfyear, fullyear, twoyear)
colnames(res) <- c("pvalue", "FUtime")
res <- as.data.frame(res)
res$FUtime_fct <- factor(res$FUtime, levels = c("0", "6", "12", "24"))
res$FUtime_fct <- fct_recode(res$FUtime, "0 months" = "0", "6 months" = "6", 
        "12 months" = "12", "24 months" = "24")
res$pvalue <- as.numeric(res$pvalue)
dbray$FUtime_fct <- fct_recode(as.character(dbray$FUtime), 
    "0 months" = "0", "6 months" = "6", "12 months" = "12", "24 months" = "24")
dbray$FUtime_fct <- fct_relevel(dbray$FUtime_fct, "6 months", after = 1L)

# Duplicate records of donors to have data at the sep time points
ids_to_repeat <- dbray |> filter(Donor == TRUE) |>  pull(SampleID)
df_repeat6 <- dbray |> filter(SampleID %in% ids_to_repeat & FUtime == 0) |> 
  mutate(Timepoint = "3", FUtime = 6, FUtime_fct = "6 months")
df_repeat12 <- dbray |> filter(SampleID %in% ids_to_repeat & FUtime == 0) |> 
  mutate(Timepoint = "4", FUtime = 12, FUtime_fct = "12 months")
df_repeat24 <- dbray |> filter(SampleID %in% ids_to_repeat & FUtime == 0) |> 
  mutate(Timepoint = "5", FUtime = 24, FUtime_fct = "24 months")
df_new <- bind_rows(dbray, df_repeat6, df_repeat12, df_repeat24) |> 
    mutate(Donor = case_when(Donor == TRUE ~ "Donor", .default = "Recipient"),
            FUtime_fct = as.factor(FUtime_fct),
            FUtime_fct = fct_relevel(FUtime_fct, "6 months", after = 1L)) 

(braycurt <- df_new %>%
                ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = Donor), size = 3, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(8)[c(7,6)]) +
                    scale_fill_manual(values = pal_nejm()(8)[c(7,6)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = Donor, 
                        fill = Donor), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    facet_wrap(~ FUtime_fct) +
                    geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("PERMANOVA p = ", round(pvalue, 3))),
                              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_donorrecipient.pdf", width = 8, height = 8)
