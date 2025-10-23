# Linear mixed for shotgun IMITHOT
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(afex)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold", size = rel(0.7))
        ))
    
} 

# Data
## Open phyloseq object
sg <- readRDS("data/microbiome_filtered_pt.RDS")
head(sg)
sg <- as.data.frame(sg)
taxno <- ncol(sg)
sgrel_nolog <- sg |> mutate(across(everything(), function (x) x * 100))
sgrel_nolog$SampleID <- rownames(sgrel_nolog)
sgrel <- sg |> mutate(across(everything(), function (x) log10((x * 100) + 0.01)))
head(sgrel)
sgrel$SampleID <- rownames(sgrel)

# Metadata
clindata <- readRDS("data/imithot_meta.RDS")
clindf <- clindata |> filter(SampleID %in% rownames(sgrel))
dim(clindf)

df_tot <- left_join(sgrel, clindf, by = "SampleID")

statres <- c()
for(i in c(1:taxno)) {
    df_tot1 <- df_tot %>% filter(FUtime %in% c(0, 6))
    df_tot1$tax <- df_tot1[,i]
    taxname <- colnames(df_tot1)[i]
    model1 <- lmer(tax ~ Group*FUtime + (1|ID), 
                      data = df_tot1)
    res <- summary(model1)
    confint_model1 <- confint(model1)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model1[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model1[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(taxname, group1 = 0, group2 = 6, pval, 
                          sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
    
    df_tot2 <- df_tot %>% filter(FUtime %in% c(0, 12))
    df_tot2$tax <- df_tot2[,i]
    taxname <- colnames(df_tot2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = df_tot2)
    res <- summary(model2)
    confint_model2 <- confint(model2)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model2[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model2[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(taxname, group1 = 0, group2 = 12, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
  
  df_tot2 <- df_tot %>% filter(FUtime %in% c(0, 24))
    df_tot2$tax <- df_tot2[,i]
    taxname <- colnames(df_tot2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = df_tot2)
    res <- summary(model2)
    confint_model2 <- confint(model2)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model2[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model2[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(taxname, group1 = 0, group2 = 24, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
}

statres <- as.data.frame(statres)
statres <- statres |> 
    arrange(pval, group1) |> 
    mutate(qval = p.adjust(pval, method = "fdr"),
           group1 = as.numeric(group1),
           group2 = as.numeric(group2))
head(statres, n = 20)
maxsig <- statres %>% filter(pval <= 0.05) %>% filter(!duplicated(taxname))
maxsig24 <- statres |> filter(pval <=0.05 & group2 == 24)

df_tot <- left_join(sgrel_nolog, clindf, by = "SampleID")
plist <- list()
for(i in 1:nrow(maxsig)){
    taxname2 <- maxsig$tax[i]
    pval <- maxsig$pval[i]
    taxnumber <- maxsig$taxname[i]
    df_tot$tax <- df_tot[,maxsig$taxname[i]]
    df_means <- df_tot %>% group_by(Group, FUtime) %>% 
        summarise(mean = mean(tax), sd = sd(tax), n = length(tax))
    res_lmm <- statres %>% filter(taxname == taxname2) %>% dplyr::select(-taxname) %>% filter(sig != "")
    asvmax <- log10(max(df_tot$tax)*1.1)
    asvmax <- ifelse(asvmax < 0, asvmax*-1, asvmax)
    posy <- if (nrow(res_lmm) == 1) {
                0.65 * asvmax
              } else if (nrow(res_lmm) == 2) {
                c(0.70, 0.65) * asvmax
              } else if (nrow(res_lmm) == 3) {
                c(0.75, 0.65, 0.55) * asvmax
              }
    print(posy)
    pl2 <- ggplot() +
        geom_line(data = df_tot, aes(x = FUtime, y = tax,
                                    color = Group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_tot, aes(x = FUtime, y = tax,
                                     color = Group, group = Group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = FUtime, y = mean, 
                                          color = Group, group = Group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = FUtime, y = mean, 
                                           color = Group, group = Group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = FUtime,
                          color = Group), width=0.1) +
        stat_pvalue_manual(res_lmm, y.position = posy, label = "{sig}", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_manual(values = c(pal_nejm()(8)[c(8,6)])) + 
        scale_y_log10() +
        scale_x_continuous(breaks = c(0,6,12,24)) +
        theme_Publication() +
        labs(x = "Months", y = "log10(relative abundance)", title = taxname2,
             color = "") 
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:11],
          nrow = 3, ncol = 4))
ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.pdf", width = 14, height = 12)
ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.png", width = 14, height = 12)
ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.svg", width = 14, height = 12)
write.csv2(statres, "results/microbiome/lmer/lmm_results.csv")

## Subset of highest ranked predictors ML
topranked <- rio::import("hashimoto/output_XGB_class_hashimoto_2025_10_23__16-46-23/feature_importance.txt")
topranked <- topranked |> arrange(-RelFeatImp) |> slice(1:20)
dftop <- df_tot1 |> select((ncol(df_tot1)-7):ncol(df_tot1), any_of(topranked$FeatName))

(featimppl <- ggplot(topranked |> mutate(FeatName = fct_rev(fct_inorder(FeatName)))) +
    geom_bar(aes(x = FeatName, y = RelFeatImp), 
        stat = "identity", fill = pal_nejm()(7)[6]) +
    labs(x = "", y = "Relative importance (%)") +
    coord_flip() +
    theme_Publication())
ggsave("results/microbiome/featimp_hashimoto.pdf", width = 8, height = 7)

statres <- c()
for(i in c(8:38)) {
    df_tot1 <- df_tot1 %>% filter(FUtime %in% c(0, 6))
    df_tot1$tax <- df_tot1[,i]
    taxname <- colnames(df_tot1)[i]
    model1 <- lmer(tax ~ Group*FUtime + (1|ID), 
                      data = df_tot1)
    res <- summary(model1)
    confint_model1 <- confint(model1)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model1[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model1[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(taxname, group1 = 0, group2 = 6, pval, 
                          sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
    
    df_tot2 <- df_tot %>% filter(FUtime %in% c(0, 12))
    df_tot2$tax <- df_tot2[,i]
    taxname <- colnames(df_tot2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = df_tot2)
    res <- summary(model2)
    confint_model2 <- confint(model2)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model2[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model2[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(taxname, group1 = 0, group2 = 12, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
  
  df_tot2 <- df_tot %>% filter(FUtime %in% c(0, 24))
    df_tot2$tax <- df_tot2[,i]
    taxname <- colnames(df_tot2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = df_tot2)
    res <- summary(model2)
    confint_model2 <- confint(model2)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model2[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model2[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(taxname, group1 = 0, group2 = 24, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
}

statres <- as.data.frame(statres)
statres <- statres |> 
    arrange(pval, group1) |> 
    mutate(qval = p.adjust(pval, method = "fdr"),
           group1 = as.numeric(group1),
           group2 = as.numeric(group2))
head(statres, n = 20)
maxsig <- statres %>% filter(pval <= 0.05) %>% filter(!duplicated(taxname))

df_tot <- left_join(sgrel_nolog, clindf, by = "SampleID")
plist <- list()
for(i in 1:nrow(maxsig)){
    taxname2 <- maxsig$tax[i]
    pval <- maxsig$pval[i]
    taxnumber <- maxsig$taxname[i]
    df_tot$tax <- df_tot[,maxsig$taxname[i]]
    df_means <- df_tot %>% group_by(Group, FUtime) %>% 
        summarise(mean = mean(tax), sd = sd(tax), n = length(tax))
    res_lmm <- statres %>% filter(taxname == taxname2) |> 
        dplyr::select(-taxname) |> filter(sig != "")
    asvmax <- log10(max(df_tot$tax*1.1))
    pl2 <- ggplot() +
        geom_line(data = df_tot, aes(x = FUtime, y = tax,
                                    color = Group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_tot, aes(x = FUtime, y = tax,
                                     color = Group, group = Group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = FUtime, y = mean, 
                                          color = Group, group = Group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = FUtime, y = mean, 
                                           color = Group, group = Group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = FUtime,
                          color = Group), width=0.1) +
        stat_pvalue_manual(res_lmm, y.position = asvmax*0.7, label = "{sig}", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_manual(values = c(pal_nejm()(8)[c(8,6)])) + 
        scale_y_log10() +
        scale_x_continuous(breaks = c(0,6,12,24)) +
        theme_Publication() +
        labs(x = "Months", y = "log10(relative abundance)", title = taxname2,
             color = "") 
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:2],
          nrow = 1, ncol = 2))
ggsave(plots, filename = "results/microbiome/lmer/lmer_plots_hmstrains.pdf", width = 8, height = 4)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.png", width = 14, height = 8)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.svg", width = 14, height = 8)
write.csv2(statres, "results/microbiome/lmer/lmm_results_hmstrains.csv")
