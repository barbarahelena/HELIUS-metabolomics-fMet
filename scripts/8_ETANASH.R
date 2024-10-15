## fMet ELISA plots

library(dplyr)
library(ggsci)
library(tidyverse)
library(ComplexHeatmap)
library(Hmisc)

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
                axis.text.x = element_text(angle = 0, vjust=1), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.title = element_text(size = rel(0.7)),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

## Data
df <- rio::import("data/FMET_ETANASH.xlsx")

## Cleaning
meta <- df %>% select(Subject = patient.ID, fMet = Abundance, AB = ABX, 
                      MMT, NASH) %>% 
    filter(NASH == "NASH") %>% 
    mutate(ID = as.factor(Subject)) %>% 
    pivot_wider(., id_cols = Subject, names_from = c(AB, MMT),
                values_from = fMet, names_sep = " - ")

long <- meta %>% pivot_longer(., 2:5, names_to = c("ABX", "MMT"), names_sep = " - ") %>% 
    mutate(ABX = fct_relevel(ABX, "No ABX")) %>% 
    filter(MMT == "Fasted" & !is.na(value)) %>% 
    filter(Subject != 3)


# Plot ELISA data
mycolors <- colorRampPalette(pal_jco()(9))(18)

ggplot(long, aes(x = ABX, y = value, color = as.factor(Subject),
                  group = Subject)) +
    geom_point(size = 2.5, alpha = 0.75) +
    geom_line(size = 1, alpha = 0.75) +
    scale_y_continuous(limits = c(0, 1.5)) +
    scale_color_jco() +
    theme_Publication() +
    #facet_wrap(~platform, scales = "free")+
    labs(title = "fMet change after antibiotics",
         color = "Subject")
ggsave("Metabolon_ABchange.pdf", width = 5, height = 5)

(pl2 <- ggplot(data = long, aes(x = ABX, y = value)) +
        geom_line(aes(group = Subject), 
                  alpha = 0.5, color = "grey40") +
        geom_jitter(aes(x = ABX, y = value, group = Subject, 
                    color = ABX), alpha = 0.5, width = 0) +
        gghalves::geom_half_violin(aes(x = ABX, y = value, fill = ABX),
                                   side = c("l", "r"), nudge = 0.05) +
        gghalves::geom_half_boxplot(data = long %>% filter(ABX == "No ABX"), 
                                    aes(x = ABX, y = value), 
                                    nudge = 0.05, side = "l", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        gghalves::geom_half_boxplot(data = long %>% filter(ABX == "ABX"), 
                                    aes(x = ABX, y = value), 
                                    nudge = 0.05, side = "r", fill = "white", color = "black", 
                                    width = 0.25, outlier.shape = NA, errorbar.draw = FALSE) +
        stat_compare_means(method = 'wilcox.test', paired = TRUE) +
        scale_y_continuous(limits = c(0, 1.5)) +
        scale_color_jama(guide = "none") +
        scale_fill_jama(guide = "none") +
        theme_Publication() +
        labs(x = "Treatment", y = "Relative concentration", title = "fMet change",
             color = ""))
ggsave("results/fMetchange.pdf", width = 4, height = 6)
ggsave("results/fMetchange.svg", width = 4, height = 6)
