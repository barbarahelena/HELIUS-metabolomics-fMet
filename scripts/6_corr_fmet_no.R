## Correlations with NO metabolites and fMet
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
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
                strip.text = element_text(face="bold")
        ))
    
} 

## Data 
setwd("~/Documents/VUmc/fMet/fmet_analyses/XGB")
helius <- readRDS("data/heliusdf.RDS")
helius$ID <- as.character(helius$Heliusnr)
mat2 <- readRDS('data/Met_HELIUS_noXenobiotics.RDS')
infomet <- readRDS('data/infomet.RDS')
nomet <- infomet$met[which(infomet$sub=='Urea cycle, Arginine and Proline Metabolism')]
cyst <- infomet$met[which(str_detect(infomet$met,'cysteine'))]
shortlist <- c(nomet, "N-formylmethionine", cyst)
shortlist <- shortlist[-22]
mat3 <- mat2[which(rownames(mat2) %in% shortlist),]

mets <- as.data.frame(t(as.matrix(mat3)))
mets <- mets[1:nrow(mets)-1,]
mets$ID <- rownames(mets)
mets <- mets %>% mutate_at(c(1:24), as.numeric)
metshelius <- left_join(mets, helius, by = "ID")

list2 <- list()
for(i in shortlist){
    df <- metshelius
    df$met <- df[,i]
    pl <- ggplot(df, aes(x=log10(`N-formylmethionine`), y=log10(met+0.0001))) +
        geom_point(alpha = 0.5, color = pal_lancet()(1)) +
        geom_smooth(formula = y~x, method = "glm", se = FALSE, color = pal_lancet()(2)[2])+
        stat_cor(method = "spearman", alternative = "two.sided") +
        labs(title = i, y = i, x = "Formylmethionine") +
        theme_Publication() +
        theme(legend.position = "none")
    list2[[i]] <- pl
}

ggarrange(plotlist = list2,
          nrow = 4, ncol = 6, common.legend = TRUE,
          labels = c(LETTERS[1:24]))
ggsave("results/corr_no_fMet_1.pdf", width = 10, height = 15)
ggsave("results/corr_no_fMet_1.svg", width = 10, height = 15)

ggarrange(list2[[13]], list2[[14]],list2[[15]], list2[[16]],
          list2[[17]],
          nrow = 2, ncol = 3, common.legend = TRUE,
          labels = c(LETTERS[13:18]))
ggsave("results/corr_no_fMet_2.pdf", width = 10, height = 8)
ggsave("results/corr_no_fMet_2.svg", width = 10, height = 8)


(pl <- ggplot(metshelius, aes(x=HT, y=log10(`N-formylmethionine`))) +
    geom_violin(aes(fill = HT)) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ggsci::scale_fill_jco(guide = "none") +
    stat_compare_means()+
    # geom_point(alpha = 0.5, color = pal_lancet()(1)) +
    # geom_smooth(formula = y~x, method = "glm", se = FALSE, color = pal_lancet()(2)[2])+
    # stat_cor(method = "spearman", alternative = "two.sided") +
    theme_Publication())

(pl <- ggplot(metshelius, aes(x=Ethnicity, y=log10(`N-formylmethionine`))) +
        geom_violin(aes(fill = Ethnicity)) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        ggsci::scale_fill_jco(guide = "none") +
        stat_compare_means()+
        # geom_point(alpha = 0.5, color = pal_lancet()(1)) +
        # geom_smooth(formula = y~x, method = "glm", se = FALSE, color = pal_lancet()(2)[2])+
        # stat_cor(method = "spearman", alternative = "two.sided") +
        theme_Publication())

(pl <- ggplot(metshelius, aes(x=Sex, y=log10(`N-formylmethionine`))) +
        geom_violin(aes(fill = Sex)) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        ggsci::scale_fill_jco(guide = "none") +
        stat_compare_means()+
        # geom_point(alpha = 0.5, color = pal_lancet()(1)) +
        # geom_smooth(formula = y~x, method = "glm", se = FALSE, color = pal_lancet()(2)[2])+
        # stat_cor(method = "spearman", alternative = "two.sided") +
        theme_Publication())

(pl <- ggplot(metshelius, aes(x=, y=log10(`N-formylmethionine`))) +
        geom_violin(aes(fill = Ethnicity)) +
        geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        ggsci::scale_fill_jco(guide = "none") +
        stat_compare_means()+
        # geom_point(alpha = 0.5, color = pal_lancet()(1)) +
        # geom_smooth(formula = y~x, method = "glm", se = FALSE, color = pal_lancet()(2)[2])+
        # stat_cor(method = "spearman", alternative = "two.sided") +
        theme_Publication())


(pl <- ggplot(metshelius, aes(x=Age, y=log10(`N-formylmethionine`))) +
        # geom_violin(aes(fill = Ethnicity)) +
        # geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
        # ggsci::scale_fill_jco(guide = "none") +
        # stat_compare_means()+
        geom_point(alpha = 0.5, color = pal_lancet()(1)) +
        geom_smooth(formula = y~x, method = "glm", se = FALSE, color = pal_lancet()(2)[2])+
        stat_cor(method = "pearson", alternative = "two.sided") +
        theme_Publication())
