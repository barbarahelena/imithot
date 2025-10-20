# TCAM plots IMITHOT
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

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


plot_loadings <- function(data, title) {
  data$comp <- data[[2]]
  ggplot(data, aes(x = reorder(X, comp), y = comp)) +
    geom_bar(stat = "identity", 
        fill = c(rep("firebrick", 10), rep("royalblue",10))) +
    coord_flip() +
    labs(title = title, x = "", y = "Loading") +
    theme_Publication()
}

## Load data
tcam <- read.csv("results/microbiome/tcam/pythonoutput_plot.csv") |> 
  select(1:17) %>% 
  mutate(Group = as.factor(Group))
head(tcam)
load <- read.csv("results/microbiome/tcam/df_loadings.csv")
head(load)
mb <- readRDS("data/filtered_imputed_microbiome_intervention.RDS")
meta <- readRDS("data/imithot_meta.RDS")

## PERMANOVAs ##
tcam <- tcam %>% filter(FUtime == 0) 
f1n <- names(tcam)[7]
f2n <- names(tcam)[8]

set.seed(1234)
# PERMANOVA for overall Group effect
tcam_dist <- tcam %>% select(all_of(c(f1n, f2n))) %>% dist(method = "euclidean")
(permanova_group <- adonis2(tcam_dist ~ Group, data = tcam, permutations = 999, method = "euclidean"))

annotation_text <- str_c("p = ", permanova_group$`Pr(>F)`[1])

### TCAM plot ###
(tcampl <- ggplot(data = tcam, aes(x = .data[[f1n]], y = .data[[f2n]], 
                        color = Group, fill = Group)) +
        #stat_ellipse(geom = "polygon", alpha = 0.3) +
        geom_point(size = 3) +
        ggtitle('TCAM') +
        scale_color_manual(values = pal_nejm()(6)[c(6,3)]) +
        scale_fill_manual(values = pal_nejm()(6)[c(6,3)]) +
        labs(x=str_c(str_replace(f1n, "[.]", " "), "%"),
            y=str_c(str_replace(f2n, "[.]", " "), "%")) +
        annotate("text", x = Inf, y = Inf, 
                label = annotation_text,
                hjust = 1, vjust = 1, 
                size = 4) +
        theme_Publication() +
        theme(legend.title = element_blank()))
ggsave("results/microbiome/tcam/f1f2_scatter.pdf", width = 6, height = 6) 

## Loading of component plots ##
head(load)
names(load)[1:10]
last <- nrow(load)
min10 <- last - 9
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:10, min10:last),]
f2 <- load %>% select(X, all_of(f2n)) %>% arrange(-.data[[f2n]])
f2 <- f2[c(1:10, min10:last),]

(f1plot <- plot_loadings(f1, "TCAM F1"))
ggsave("results/microbiome/tcam/loading_pc1.pdf", width = 8, height = 7)

plot_loadings(f2, "TCAM F2")
ggsave("results/microbiome/tcam/loading_pc2.pdf", width = 8, height = 7)

## Lineplots ##
f1 <- load %>% select(X, all_of(f1n)) %>% arrange(-.data[[f1n]])
f1 <- f1[c(1:10, (min10:last)),]
pseudocounts <- mb %>%  select(any_of(f1$X)) %>%
  summarise(across(everything(), ~ min(.x[.x > 0], na.rm = TRUE) / 2))
mbsel <- mb %>% ungroup(.) %>% select(any_of(f1$X), ID, Group, FUtime) %>%
  mutate( across(f1$X, ~log10(.x + pseudocounts[[cur_column()]]))
          ) %>%
  rename_at(c(f1$X), ~str_remove(.x, " \\(.*\\)$"))
head(mbsel)
means_ses <- mbsel %>% # Calculate means and standard deviations per Genotype
  group_by(Group, FUtime) %>%
  summarise(across(1:20, list(mean = ~mean(.x, na.rm = TRUE), 
                        se = ~(sd(.x, na.rm = TRUE) / sqrt(n()))),
                      .names = "{.col}-{.fn}"))
means_ses

means_ses_long <- means_ses %>% # to long format for plotting
  pivot_longer(cols = c(-Group, -FUtime), 
                        names_to = c("microbe", ".value"), names_sep = "-")
head(means_ses_long)

mbsel_long <- mbsel %>%
  pivot_longer(cols = c(-ID, -Group, -FUtime), 
                        names_to = "microbe", values_to = "value")

(lineplots <- ggplot(means_ses_long, aes(x = FUtime, y = mean, 
                                    group = Group, color = Group)) +
                      geom_line() +
                      geom_point(size = 0.5) +
                      geom_jitter(data = mbsel_long, aes(x = FUtime, y = value, color = Group),
                                    width = 0.1, alpha = 0.7) +
                      scale_color_manual(values = pal_nejm()(8)[c(8,6)]) +
                      geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
                      labs(title = "Abundance over time - Autologous vs Allogenic",
                          x = "Follow-up time (months)",
                          y = "log10(relative abundance)") +
                      scale_x_continuous(breaks = c(0,3,6,12,24)) +
                      facet_wrap(~ microbe, scales = "free", nrow = 4) +
                      theme_Publication() +
                      theme(strip.text = element_text(size = 8)))
ggsave("results/microbiome/tcam/lineplots_groups_log10.pdf", width = 20, height = 17)

