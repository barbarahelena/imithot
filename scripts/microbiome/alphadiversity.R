# Alpha diversity plots IMITHOT
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(ggsci)

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
                axis.text = element_text(size = rel(0.9)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(1.0, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

## Load data
mb <- readRDS("data/shotgun_cleaned.RDS")
mb <- as.data.frame(mb)
mb <- mb |> mutate(across(everything(), function (x) x*100))
df <- readRDS("data/imithot_meta.RDS")
head(mb)[1:5,1:5]
rowSums(mb)

# Diversity metrics between 4 groups
## Shannon plot
shannon <- vegan::diversity(mb, index = 'shannon')
df_shan <- data.frame(SampleID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df, by = "SampleID")
(plshan <- ggplot(data = df_shan |> filter(Group != "Healthy"), aes(x = as.factor(FUtime), 
                        y = shannon, fill = Group,
                        group = interaction(as.factor(FUtime), Group))) +
    scale_fill_manual(values = pal_nejm()(8)[c(8,6)]) +
    stat_compare_means(aes(group = interaction(as.factor(FUtime), Group)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Shannon index", y = "Shannon index", x="FU time (months)") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon_time.pdf", width = 7, height = 5)

(plshan <- ggplot(data = df_shan |> filter(FUtime == 0), 
                    aes(x = Hashimoto, y = shannon, fill = Hashimoto)) +
    scale_fill_manual(values = pal_nejm()(8)[c(7,6)], guide = "none") +
    stat_compare_means(method = "wilcox.test", label = "p.format") + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Shannon index", y = "Shannon index", x="") +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon_baseline.pdf", width = 5, height = 5)

## Simpsons
simpson <- vegan::diversity(mb, index = 'simpson')
df_simp <- data.frame(SampleID = names(simpson), simpson = simpson)
df_simp <- left_join(df_simp, df, by = "SampleID")
comp <- list(c("0", "6"), c("0", "12"), c("0", "24"))
(plsimp <- ggplot(data = df_simp |> filter(Group != "Healthy"), aes(x = as.factor(FUtime), 
                        y = simpson, fill = Group,
                        group = interaction(as.factor(FUtime), Group))) +
    scale_fill_manual(values = pal_nejm()(8)[c(8,6)]) +
    stat_compare_means(comparisons = comp, paired = TRUE,
        label = "p.signif", hide.ns = TRUE) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Simpsons index", y = "Simpsons index", x="FU time (months)") +
    facet_wrap(~ Group) +
    theme_Publication())
#ggsave(plshan, filename = "results/microbiome/alphadiversity/shannon.svg", width = 12, height = 5)
ggsave(plsimp, filename = "results/microbiome/alphadiversity/simpson.pdf", width = 12, height = 5)

## Species richness
specrich <- specnumber(mb)
dfspec <- data.frame(SampleID = names(specrich), richness = specrich)
dfspec <- left_join(dfspec, df, by = "SampleID") |> filter(Group != "Healthy")
(plrich <- ggplot(data = dfspec, 
                 aes(x = as.factor(FUtime), 
                        y = richness, fill = Group,
                        group = interaction(as.factor(FUtime), Group))) +
    scale_fill_manual(values = pal_nejm()(8)[c(8,6)]) +
    geom_boxplot(outlier.shape = NA) +
    # stat_compare_means(aes(group = interaction(as.factor(FUtime), Group)), 
    #         method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Species richness", y = "Number of species", x="FU time (months)") +
    facet_wrap(~ Group) +
    theme_Publication())
(plrich <- ggplot(data = dfspec, 
                 aes(x = as.factor(FUtime), 
                        y = richness, fill = Group,
                        group = interaction(as.factor(FUtime), Group))) +
    scale_fill_manual(values = pal_nejm()(8)[c(8,6)]) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(aes(group = interaction(as.factor(FUtime), Group)), 
            method = "wilcox.test", label = "p.signif", hide.ns = TRUE) + 
    geom_jitter(position = position_dodge(0.75)) +
    labs(title = "Species richness", y = "Number of species", x="FU time (months)") +
    theme_Publication())
ggsave(plrich, filename = "results/microbiome/alphadiversity/richness.pdf", width = 4, height = 5)
#ggsave(plrich, filename = "results/microbiome/alphadiversity/richness.svg", width = 4, height = 5)