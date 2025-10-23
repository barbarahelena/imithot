# Linear mixed for metabolomics IMITHOT
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
met <- readRDS('data/metabolomics_clean.RDS')
head(met)[,1:8]
tail(met)[,1:8]
met <- as.data.frame(met)
metno <- ncol(met)
met$ID
met <- met |> filter(Group != "Healthy")

statres <- c()
for(i in c(3:(metno-4))) {
    met1 <- met %>% filter(FUtime %in% c(0, 6))
    met1$metabolite <- met1[[i]]
    metname <- colnames(met1)[i]
    print(i); print(metname)
    model1 <- lmer(metabolite ~ Group*FUtime + (1|ID), data = met1)
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
    statres_line <- cbind(metname, group1 = 0, group2 = 6, pval, 
                          sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
    
    met2 <- met %>% filter(FUtime %in% c(0, 12))
    met2$tax <- met2[,i]
    metname <- colnames(met2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = met2)
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
    statres_line <- cbind(metname, group1 = 0, group2 = 12, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
  
  met2 <- met %>% filter(FUtime %in% c(0, 24))
    met2$tax <- met2[,i]
    metname <- colnames(met2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = met2)
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
    statres_line <- cbind(metname, group1 = 0, group2 = 24, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
}

statres <- as.data.frame(statres)
statres <- statres |> 
    # filter(group2 == "24") |> 
    arrange(pval, group1) |> 
    mutate(qval = p.adjust(pval, method = "fdr"),
           group1 = as.numeric(group1),
           group2 = as.numeric(group2))
head(statres, n = 20)
maxsig <- statres %>% filter(pval <= 0.05) %>% filter(!duplicated(metname))
maxsig6 <- statres |> filter(group2 == 6 & pval <= 0.05)
maxsig12 <- statres |> filter(group2 == 12 & pval <= 0.05)
maxsig24 <- statres |> filter(group2 == 24 & pval <= 0.05)

# There is more overlap 1 and 2 years than 6 months with 1/2 years
maxsig6 |> filter(metname %in% maxsig12$metname) # 2 met
maxsig6 |> filter(metname %in% maxsig24$metname) # 2 met
maxsig12 |> filter(metname %in% maxsig24$metname) # 21 met

plist <- list()
for(i in 1:nrow(maxsig6)){
    metname2 <- maxsig6$metname[i]
    pval <- maxsig6$pval[i]
    met$metabolite <- met[,maxsig6$metname[i]]
    df_means <- met %>% group_by(Group, FUtime) %>% 
        summarise(mean = mean(metabolite), sd = sd(metabolite), n = length(metabolite))
    res_lmm <- statres |> 
          filter(metname == metname2) |> 
          dplyr::select(-metname) |> 
          filter(sig != "") |> 
          arrange(group2)
    metmax <- max(met$metabolite*1.1)
    posy <- if (nrow(res_lmm) == 1) {
                0.7 * metmax
              } else if (nrow(res_lmm) == 2) {
                c(0.7, 0.6) * metmax
              } else if (nrow(res_lmm) == 3) {
                c(0.7, 0.67, 0.64) * metmax
              }
        metname2 <- maxsig6$metname[i]
    metname2 <- gsub("\\\\n", "\n", str_wrap(metname2, width = 25))
  print(metname2)
    pl2 <- ggplot() +
        geom_line(data = met, aes(x = FUtime, y = metabolite,
                                    color = Group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = met, aes(x = FUtime, y = metabolite,
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
        scale_x_continuous(breaks = c(0,6,12,24)) +
        theme_Publication() +
        labs(x = "Months", y = "log10(relative abundance)", title = metname2,
             color = "") 
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist[1:21], common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:21],
          nrow = 6, ncol = 5))
ggsave(plots, filename = "results/metabolomics/lmer/lmer_plots_6m.pdf", width = 16, height = 21)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.png", width = 12, height = 8)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.svg", width = 12, height = 8)
write.csv2(statres, "results/metabolomics/lmer/lmm_results.csv")


plist <- list()
for(i in 1:nrow(maxsig12)){
    metname2 <- maxsig12$metname[i]
    pval <- maxsig12$pval[i]
    met$metabolite <- met[,maxsig12$metname[i]]
    df_means <- met %>% group_by(Group, FUtime) %>% 
        summarise(mean = mean(metabolite), sd = sd(metabolite), n = length(metabolite))
    res_lmm <- statres |> 
          filter(metname == metname2) |> 
          dplyr::select(-metname) |> 
          filter(sig != "") |> 
          arrange(group2)
    metmax <- max(met$metabolite*1.1)
    posy <- if (nrow(res_lmm) == 1) {
                0.7 * metmax
              } else if (nrow(res_lmm) == 2) {
                c(0.7, 0.55) * metmax
              } else if (nrow(res_lmm) == 3) {
                c(0.77, 0.65, 0.52) * metmax
              }
    metname2 <- maxsig12$metname[i]
    metname2 <- gsub("\\\\n", "\n", str_wrap(metname2, width = 25))
  print(metname2)
    pl2 <- ggplot() +
        geom_line(data = met, aes(x = FUtime, y = metabolite,
                                    color = Group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = met, aes(x = FUtime, y = metabolite,
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
        scale_x_continuous(breaks = c(0,6,12,24)) +
        theme_Publication() +
        labs(x = "Months", y = "log10(relative abundance)", title = metname2,
             color = "") 
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist[1:43], common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:43],
          nrow = 9, ncol = 5))
ggsave(plots, filename = "results/metabolomics/lmer/lmer_plots_12m.pdf", width = 16, height = 32)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.png", width = 12, height = 8)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.svg", width = 12, height = 8)

plist <- list()
for(i in 1:nrow(maxsig24)){
    metname2 <- maxsig24$metname[i]
    pval <- maxsig24$pval[i]
    met$metabolite <- met[,maxsig24$metname[i]]
    df_means <- met %>% group_by(Group, FUtime) %>% 
        summarise(mean = mean(metabolite), sd = sd(metabolite), n = length(metabolite))
    res_lmm <- statres |> 
          filter(metname == metname2) |> 
          dplyr::select(-metname) |> 
          filter(sig != "") |> 
          arrange(group2)
    metmax <- max(met$metabolite*1.1)
    posy <- if (nrow(res_lmm) == 1) {
                0.7 * metmax
              } else if (nrow(res_lmm) == 2) {
                c(0.7, 0.55) * metmax
              } else if (nrow(res_lmm) == 3) {
                c(0.77, 0.65, 0.52) * metmax
              }
    metname2 <- maxsig24$metname[i]
    metname2 <- gsub("\\\\n", "\n", str_wrap(metname2, width = 25))
  print(metname2)
    pl2 <- ggplot() +
        geom_line(data = met, aes(x = FUtime, y = metabolite,
                                    color = Group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = met, aes(x = FUtime, y = metabolite,
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
        scale_x_continuous(breaks = c(0,6,12,24)) +
        theme_Publication() +
        labs(x = "Months", y = "log10(relative abundance)", title = metname2,
             color = "") 
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist[1:80], common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:80],
          nrow = 9, ncol = 9))
ggsave(plots, filename = "results/metabolomics/lmer/lmer_plots_24m.pdf", width = 32, height = 32)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.png", width = 12, height = 8)
# ggsave(plots, filename = "results/microbiome/lmer/lmer_plots.svg", width = 12, height = 8)


# Baseline differences 
topranked <- rio::import("hashimoto_metabolites/output_XGB_class_hashimoto_metabolites_2025_10_23__23-29-49/feature_importance.txt")
topranked <- topranked |> arrange(-RelFeatImp) |> slice(1:20)

(featimppl <- ggplot(topranked |> mutate(FeatName = fct_rev(fct_inorder(FeatName)))) +
    geom_bar(aes(x = FeatName, y = RelFeatImp), 
        stat = "identity", fill = pal_nejm()(7)[6]) +
    labs(x = "", y = "Relative importance (%)") +
    coord_flip() +
    theme_Publication())
ggsave("results/metabolomics/featimp_hashimoto.pdf", width = 8, height = 7)

mettop <- met |> select(SampleID, ID, Group, FUtime, any_of(topranked$FeatName))

statres <- c()
for(i in c(5:24)) {
    met1 <- mettop %>% filter(FUtime %in% c(0, 6))
    met1$metabolite <- met1[[i]]
    metname <- colnames(met1)[i]
    print(i); print(metname)
    model1 <- lmer(metabolite ~ Group*FUtime + (1|ID), data = met1)
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
    statres_line <- cbind(metname, group1 = 0, group2 = 6, pval, 
                          sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
    
    met2 <- mettop %>% filter(FUtime %in% c(0, 12))
    met2$tax <- met2[,i]
    metname <- colnames(met2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = met2)
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
    statres_line <- cbind(metname, group1 = 0, group2 = 12, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
  
  met2 <- mettop %>% filter(FUtime %in% c(0, 24))
    met2$tax <- met2[,i]
    metname <- colnames(met2)[i]
    model2 <- lmer(tax ~ Group*FUtime + (1|ID), 
                   data = met2)
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
    statres_line <- cbind(metname, group1 = 0, group2 = 24, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
}

statres <- as.data.frame(statres)
statres <- statres |> 
    # filter(group2 == "24") |> 
    arrange(pval, group1) |> 
    mutate(qval = p.adjust(pval, method = "fdr"),
           group1 = as.numeric(group1),
           group2 = as.numeric(group2))
head(statres, n = 20)
maxsig <- statres %>% filter(pval <= 0.05) %>% filter(!duplicated(metname))

plist <- list()
for(i in 1:nrow(maxsig)){
    metname2 <- maxsig$metname[i]
    pval <- maxsig$pval[i]
    mettop$metabolite <- mettop[,maxsig$metname[i]]
    df_means <- mettop %>% group_by(Group, FUtime) %>% 
        summarise(mean = mean(metabolite), sd = sd(metabolite), n = length(metabolite))
    res_lmm <- statres |> 
          filter(metname == metname2) |> 
          dplyr::select(-metname) |> 
          filter(sig != "") |> 
          arrange(group2)
    metmax <- max(mettop$metabolite*1.1)
    posy <- if (nrow(res_lmm) == 1) {
                0.7 * metmax
              } else if (nrow(res_lmm) == 2) {
                c(0.7, 0.55) * metmax
              } else if (nrow(res_lmm) == 3) {
                c(0.77, 0.65, 0.52) * metmax
              }
    metname2 <- maxsig$metname[i]
    metname2 <- gsub("\\\\n", "\n", str_wrap(metname2, width = 25))
  print(metname2)
    pl2 <- ggplot() +
        geom_line(data = mettop, aes(x = FUtime, y = metabolite,
                                    color = Group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = mettop, aes(x = FUtime, y = metabolite,
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
        scale_x_continuous(breaks = c(0,6,12,24)) +
        theme_Publication() +
        labs(x = "Months", y = "log10(relative abundance)", title = metname2,
             color = "") 
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist[1:3], common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:3],
          nrow = 1, ncol = 3))
ggsave(plots, filename = "results/metabolomics/lmer/lmer_plots_hmtop.pdf", width = 8, height = 4)
