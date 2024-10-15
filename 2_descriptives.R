## Descriptives
## b.j.verhaar@amsterdamumc.nl

library(rio)
library(haven)
library(tidyverse)
library(tableone)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
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

## Data
helius <- readRDS('data/helius.RDS')
metabolomics <- readRDS('data/HELIUS_plasma_metabolites.RDS')
infomet <- readRDS('data/infomet.RDS')

############### Table 1 #######################
table1 <- helius %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, CKDEPI, ACR_KDIGO, HbA1C, TC, LDL, HDL, Trig, FramRisk) %>% 
    CreateTableOne(data=., strata = 'DM', addOverall = TRUE, test = FALSE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))

table2 <- helius %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, CKDEPI, ACR_KDIGO, HbA1C, TC, LDL, HDL, Trig, FramRisk) %>% 
    CreateTableOne(data=., strata = 'Age_cat', test = TRUE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))

table3 <- helius %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, CKDEPI, ACR_KDIGO, HbA1C, TC, LDL, HDL, Trig, FramRisk) %>% 
    CreateTableOne(data=., strata = 'Sex', test = TRUE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))

# table4 <- helius %>% 
#     select(Age, Age_cat, Sex, EthnAfr, BMI, CurrSmoking, CVD, DM, HT,  
#            AntiHT, SBP, DBP, CKDEPI, ACR_KDIGO, HbA1C, TC, LDL, HDL, Trig, FramRisk) %>% 
#     CreateTableOne(data=., strata = 'EthnAfr', test = TRUE) %>% 
#     print(nonnormal=c("FramRisk", "Trig"))

openxlsx::write.xlsx(cbind(table1,table2,table3), "results/table1.xlsx")
cbind(table1,table2,table3)

########################### PCA metabolites ############################
pca <- prcomp(t(mat), center = T, scale. = T)
df <- as.data.frame(pca$x[, 1:2])
head(df)
df$ID <- rownames(df)
df$CKD <- helius$CKD_group[match(df$ID, helius$ID)]
df$CKDT2DM <- helius$group[match(df$ID, helius$ID)]
df$Sex <- helius$Sex[match(df$ID, helius$ID)]
df$HT <- helius$HT[match(df$ID, helius$ID)]
df$DM <- helius$DM[match(df$ID, helius$ID)]
df$Age_cat <- helius$Age_cat[match(df$ID, helius$ID)]
df$EthnAfr <- helius$EthnAfr[match(df$ID, helius$ID)]

(pl <- ggplot(df, aes(x=PC1, y=PC2, color=HT)) +
    geom_point()+
    ggtitle('PCA metabolites - Hypertension') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank()))

(plC <- ggplot(df, aes(x=PC1, y=PC2, color=Sex)) +
    geom_point() +
    ggtitle('PCA plasma metabolites') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme_Publication() +
    theme(legend.title = element_blank()))

(plD <- ggplot(df, aes(x=PC1, y=PC2, color=EthnAfr)) +
    geom_point() +
    ggtitle('Ethnicity') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank()))

(plB <- ggplot(df, aes(x=PC1, y=PC2, color=Age_cat)) +
    geom_point()+
    ggtitle('Age') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank()))

(plA <- ggplot(df, aes(x=PC1, y=PC2, color=DM)) +
    geom_point()+
    ggtitle('Diabetes') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank()))

pdf("PCAs_subgroups.pdf", width = 10, height = 7)
grid.arrange(plA, plB, plC, plD, nrow=2, top=textGrob("PCA metabolites for subgroups", gp=gpar(font=2)))
dev.off()
svg("PCAs_subgroups.svg", width = 10, height = 7)
grid.arrange(plA, plB, plC, plD, nrow=2, top=textGrob("PCA metabolites for subgroups", gp=gpar(font=2)))
dev.off()

(plE <- ggplot(df, aes(x=PC1, y=PC2, color=CKDT2DM)) +
    geom_point()+
    scale_color_lancet() +
    ggtitle('PCA metabolites - CKD groups') +
    theme_Publication() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank()))

(plF <- ggplot(df, aes(x=PC1, y=PC2, color=CKD)) +
    geom_point()+
    ggtitle('PCA metabolites - CKD') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x='PC1 7.7%', y='PC2 5.4%') +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank()))

pdf("PCAs_ckdgroups.pdf", width = 12, height = 6)
grid.arrange(plE, plF, nrow=1, top=textGrob("PCA metabolites for CKD groups", gp=gpar(font=2)))
dev.off()
svg("PCAs_ckdgroups.svg", width = 12, height = 6)
grid.arrange(plE, plF, nrow=1, top=textGrob("PCA metabolites for CKD groups", gp=gpar(font=2)))
dev.off()


######################## Metabolite info ########################
summary(infomet$sup)

colpal <- c("Lipid"="#00468BFF", "Amino Acid"="#ED0000FF", "Nucleotide"="#42B540FF",
            "Cofactors and Vitamins"="#0099B4FF", "Carbohydrate"="#925E9FFF", "Peptide"="#FDAF91FF", "Partially Characterized Molecules"="#AD002AFF", "Energy"="#ADB6B6FF", "Xenobiotics"="#1B1919FF")

infomet %>% 
    group_by(sup) %>% 
    summarize(count=n()) %>% 
    arrange(desc(count)) %>% 
    ggplot(.) +
    geom_histogram(aes(x=reorder(sup, count), y=count, fill = sup), stat = "identity") +
    scale_fill_manual(values = colpal) +
    labs(title="Number of metabolites per group", x = element_blank(), y='Count') +
    coord_flip() +
    theme_Publication() +
    theme(legend.position = "none")
ggsave("metabolites.pdf", device = "pdf", width = 8, height = 6)
ggsave("metabolites.svg", device = "svg", width = 8, height = 6)

(pl_noxeno <- infomet %>% 
        group_by(sup) %>% 
        summarize(count=n()) %>% 
        arrange(desc(count)) %>% 
        filter(sup!="Xenobiotics") %>% 
        ggplot(.) +
        geom_histogram(aes(x=reorder(sup, count), y=count, fill = sup), stat = "identity") +
        scale_fill_manual(values = colpal) +
        labs(title="Groups of metabolites", x = element_blank(), y='Count') +
        coord_flip() +
        theme_Publication() +
        theme(legend.position = "none", axis.line.y = element_blank(),
              axis.ticks.y = element_blank()))
ggsave("metabolites_noxeno.pdf", device = "pdf", width = 8, height = 6)
ggsave("metabolites_noxeno.svg", device = "svg", width = 8, height = 6)

