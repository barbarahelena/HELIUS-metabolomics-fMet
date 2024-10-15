## Explained variances plot
library(tidyverse)
library(ggsci)
library(ggpubr)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_blank(),
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

ev_list <- list()
groups <- c("Fem", "Male", "Old", "Young", "DM", "All")
for(g in groups){
    li <- list.files(path = g)
    for (i in li) {
        li2 <- list.files(path = str_c(g, i, sep = '/'))
        a <- str_detect(li2, i)
        b <- !(str_detect(li2, 'PERMUTED'))
        ev_list[[str_c(g, i, sep = '_')]] <- rio::import(file.path(g, i, li2[which(a&b)], 'aggregated_metrics_regression.txt'))
    }
}

df <- data.frame()
for (i in c(1:24)) {
    group <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,1]
    bp <- str_split(names(ev_list[i]), pattern = '_', 2, simplify = T)[,2]
    ev <- ev_list[[i]]$`Median Explained Variance`
    row <- cbind(group, bp, ev)
    df <- rbind(df, row)
}

df<- df %>% 
    mutate(
        ev = as.numeric(ev),
        group = as.factor(group),
        group2 = factor(group, levels = c("All", "DM", "Young", "Old", "Fem", "Male"), 
                        labels = c("All", "No DM", "Young", "Old", "Female", "Male")),
        bp = as.factor(bp),
        bp = fct_rev(fct_relevel(bp, c("SBP", "SBPres", "DBP", "DBPres"))),
        bp_group = ifelse(bp=="DBP"|bp=="DBPres", "DBP", "SBP"),
        bp_group = factor(bp_group, levels = c("SBP", "DBP"), labels = c("Systolic BP", "Diastolic BP")),
        res_group = ifelse(bp=="SBP"|bp=="DBP", "Unadjusted values", "Residuals"),
        res_group = as.factor(res_group)                   
    ) %>% 
    arrange(group, bp)

head(df)
levels(df$res_group)

write.table(df, "expvar.csv", sep=",")

pl <- ggplot(df, aes(x=res_group, y=(ev*100), color = group2)) +
    geom_point(size=2) + 
    geom_segment(aes(x=res_group, xend=res_group, y=0, yend=(ev*100), linetype=res_group, size = res_group)) +
    scale_linetype_manual(values = c("dashed", "solid"), guide = guide_legend(reverse = TRUE)) +
    scale_size_manual(values = c(0.7, 0.9), guide = "none") +
    labs(x = '', y = 'Explained variance in %', linetype = '', size = '') +
    coord_flip() +
    facet_grid(group2~bp_group, switch = 'y') +
    scale_x_discrete(position = "top") +
    scale_color_lancet(guide = "none") +
    theme_Publication() +
    theme(
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        strip.text.y.left = element_text(angle = 0)
    )
pl
ggsave('results/expl_var_allmodels.pdf', device = 'pdf', width = 7, height = 5)
ggsave('results/expl_var_allmodels.svg', device = 'svg', width = 7, height = 5)
