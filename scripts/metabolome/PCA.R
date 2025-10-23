# PCA metabolomics
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

#### Libraries ####
library(tidyverse)
library(ggsci)
library(ggpubr)

##### Functions ####
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.line.y = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 

#### Opening data files ####
meta <- readRDS("data/imithot_meta.RDS")
dim(meta)
met <- readRDS("data/metabolomics_clean.RDS")
head(met)
dim(met)
met$SampleID
met[,1:5]
meta <- meta |> filter(SampleID %in% met$SampleID)
dim(meta)

#### PCA plot metabolites ####
head(met)[1:5,1:5]
dim(met)
metmat <- met |> select(3:(ncol(met)-4)) |> as.matrix()
rownames(metmat) <- met$SampleID
set.seed(112)
pca <- prcomp(metmat, center = FALSE, scale = TRUE)
df <- as.data.frame(pca$x[, 1:5])
summary(pca)
pc1_ev <- round(summary(pca)$importance[2,1] * 100, 2)
pc2_ev <- round(summary(pca)$importance[2,2] * 100, 2)
head(df)
df$SampleID <- rownames(df)
df$Group <- meta$Group[match(df$SampleID, meta$SampleID)]
df$FUtime <- meta$FUtime[match(df$SampleID, meta$SampleID)]
df$Hashimoto <- meta$Hashimoto[match(df$SampleID, meta$SampleID)]
df$Donor <- meta$Donor[match(df$SampleID, meta$SampleID)]
df$ID <-meta$ID[match(df$SampleID, meta$SampleID)]
df <- df %>% arrange(Group)
df <- df |> filter(!is.na(Group)) |> 
  mutate(
    FUtime = str_c(FUtime, " months"),
    FUtime = as.factor(FUtime),
    FUtime = fct_relevel(FUtime, "6 months", after = 1L),
    Donor = case_when(Donor == TRUE ~ "Donor", .default = "Recipient")
  )
head(df); tail(df)

(pca1 <- ggplot(df |> filter(Group != "Healthy"), 
          aes(x=PC1, y=PC2, color=Group)) +
        stat_ellipse(geom = "polygon", aes(fill = Group, color = Group),
                     alpha = 0.1) +
        geom_point(size = 3) +
        ggtitle('PCA plasma metabolites') +
        scale_color_manual(values = pal_nejm()(8)[c(8,6)]) +
        scale_fill_manual(values = pal_nejm()(8)[c(8,6)]) +
        labs(x=str_c('PC1 ', pc1_ev, '%'), y=str_c('PC2 ', pc2_ev, '%')) +
        theme_Publication() +
        facet_wrap(~ as.factor(FUtime)) +
        theme(legend.title = element_blank()))
ggsave("results/metabolomics/pca_time.pdf", width = 10, height = 8)

(pca2 <- ggplot(df |> filter(FUtime == "0 months"), 
          aes(x=PC1, y=PC2, color=Hashimoto)) +
        stat_ellipse(geom = "polygon", aes(fill = Hashimoto, color = Hashimoto),
                     alpha = 0.1) +
        geom_point(size = 3) +
        ggtitle('Baseline differences') +
        scale_color_manual(values = pal_nejm()(8)[c(7,6)]) +
        scale_fill_manual(values = pal_nejm()(8)[c(7,6)]) +
        labs(x=str_c('PC1 ', pc1_ev, '%'), y=str_c('PC2 ', pc2_ev, '%')) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("results/metabolomics/pca_baseline.pdf", width = 7, height = 6)
