# Beta diversity plots IMITHOT
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(doParallel)
library(ape)
registerDoParallel(8)

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
mb <- readRDS("data/microbiome_filtered_pt.RDS")
df <- readRDS("data/imithot_meta.RDS")
head(mb)[1:5,1:5]
rowSums(mb)

#### Bray-Curtis distance ####
bray <- vegan::vegdist(mb, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance_bray <- pcoord$values$Rel_corr_eig * 100
head(expl_variance_bray)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$SampleID <- rownames(dbray)
dbray <- left_join(dbray, df, by = 'SampleID') # add metadata / covariates

#### Bray-Curtis per time point ####
braypertimepoint <- function(timepoint, df = dbray, tab = mb) {
    set.seed(14)
    dfsel <- dbray %>% filter(FUtime == timepoint)
    bray <- vegan::vegdist(tab[rownames(tab) %in% dfsel$SampleID,], method = 'bray')
    return(adonis2(bray ~ Group, data = dfsel))
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
res$FUtime_fct <- fct_recode(res$FUtime, "0 months" = "0", "6 monhts" = "6", "12 months" = "12", "24 months" = "24")
res$pvalue <- as.numeric(res$pvalue)
dbray$FUtime_fct <- fct_recode(as.character(dbray$FUtime), 
    "0 months" = "0", "6 monhts" = "6", "12 months" = "12", "24 monhts" = "24")

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
                    geom_point(aes(color = as.factor(FUtime)), size = 3, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance_bray[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance_bray[2], digits = 1),'%)')) +
                    scale_color_manual(values = pal_nejm()(8)[c(3:8)]) +
                    scale_fill_manual(values = pal_nejm()(8)[c(3:8)], guide = "none") +
                    theme_Publication() +
                    labs(color = "", fill = "", title = "PCoA Bray-Curtis distance") +
                    stat_ellipse(geom = "polygon", aes(color = as.factor(FUtime), fill = as.factor(FUtime)), type = "norm",
                        alpha = 0.1, linewidth = 1.0) +
                    theme(legend.position = "top") +
                    # geom_text(data = res, aes(x = Inf, y = Inf, label = paste0("PERMANOVA p = ", round(pvalue, 3))),
                    #           hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
                    facet_wrap(~ Group))
ggsave(braycurt, filename = "results/microbiome/betadiversity/PCoA_BrayCurtis_time.pdf", width = 8, height = 8)
