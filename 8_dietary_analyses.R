# Source of fMet: dietary analyses HELIUS
# Barbara Verhaar

# Library
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

# Data
df <- haven::read_sav("data/dietary_data/221115_HELIUS data Barbara Verhaar.sav")
met <- readRDS("data/HELIUS_plasma_metabolites.RDS")
met$ID <- as.integer(rownames(met))
met <- met %>% select(ID, fMet = `N-formylmethionine`)

df2 <- df %>% select(
    ID = Heliusnr, Sex = H1_geslacht, Age = H1_lft, Ethn = H1_etniciteit,
    BMI = H1_LO_BMI, Hypertension = H1_HT_SelfBPMed, BPlow = H1_Antihypertensiva,
    BPlow_gen = H1_AntihypertensivaC02, Diuretics = H1_AntihypertensivaC03,
    Betablock = H1_AntihypertensivaC07, Calc_ant = H1_AntihypertensivaC08,
    RAASinh = H1_AntihypertensivaC09, eGFR = H1_CKDEPI_eGFR, 
    MicroAlb = H1_Microalbuminurie, ACR_KDIGO = H1_ACR_KDIGO,
    Red_meat:Alcohol_Sum) %>% 
    mutate(ID = as.integer(ID),
           across(c("Sex", "Ethn", "Hypertension", "BPlow", "BPlow_gen", 
                    "Diuretics", "Betablock", "Calc_ant", "RAASinh", "MicroAlb", 
                    "ACR_KDIGO"), ~as_factor(.x, levels = c("labels"), ordered = FALSE)),
           across(where(is.numeric), as.numeric),
           Sex = fct_recode(Sex, "Male" = "man", 
                            "Female" = "vrouw"),
           across(c("Hypertension", "BPlow", "BPlow_gen", 
                    "Diuretics", "Betablock", "Calc_ant", "RAASinh", "MicroAlb"),
                  ~fct_recode(.x, "Yes" = "Ja", "No" = "Nee"))
    ) %>% 
    naniar::replace_with_na_all(condition = ~(.x == "Missing")) %>% 
    droplevels(.)

df_total <- left_join(df2, met, by = "ID")

# Table 1
table1 <- df %>%
    select(Age, Sex, BMI, Ethn,
           eGFR, MicroAlb) %>% 
    CreateTableOne(data=.) %>% 
    print(nonnormal = "Alcohol")
write.csv2(table1, "results/table1.csv")

## Open dietary data
diet <- df_total %>% 
    select(`Total calories` = ENKcal_Sum, `Carbohydrates` = Carbo_Sum, 
           `Dietary fibres` = Fibre_Sum, `Polysaccharides` = Polysac_Sum, 
           `Mono- and disaccharides` = Modisac_Sum, `Proteins` = Prot_Sum, 
           `Protein (animal)` = Prot_ani_Sum, `Protein (veg)` = Prot_veg_Sum, 
           `Total fatty acids` = FattyAcidsTot_Sum, `Saturated fatty acids` = SFA_Sum,
           `Alcohol intake` = Alcohol_Sum, `Red meat` = Red_meat, `Fatty fish` = Fatty_fish,
           `Lean fish` = Lean_fish, Eggs, `High fat cheese` = HF_cheese, 
           `Low fast cheese` = LF_cheese, `Nuts and seeds` = Nuts_seeds,
           Legumes, `High fibre grains` = HF_grains, `Low fibre grains` = LF_grains, fMet)

# Correlations macronutrient groups and fMet
macro <- colnames(diet)[1:(ncol(diet)-1)]

list2 <- list()
for(i in macro){
    df <- diet
    df$diet <- df[[i]]
    print(i)
    print(names(df))
    pl <- ggplot(df, aes(x=diet, y=fMet)) +
        geom_point(alpha = 0.5, color = pal_lancet()(1)) +
        geom_smooth(formula = y~x, method = "lm", se = FALSE, color = pal_lancet()(2)[2])+
        stat_cor(method = "spearman", alternative = "two.sided") +
        labs(title = i, x = i, y = "Formylmethionine") +
        theme_Publication() +
        theme(legend.position = "none")
    list2[[i]] <- pl
}

ggarrange(list2[[1]], list2[[2]],list2[[3]], list2[[4]],
          list2[[5]], list2[[6]],list2[[7]], list2[[8]],
          list2[[9]], list2[[10]],list2[[11]],list2[[12]],
          nrow = 4, ncol = 3, common.legend = TRUE,
          labels = c(LETTERS[1:12]))
ggsave("results/corr_diet_fMet_1.pdf", width = 10, height = 15)
ggsave("results/corr_diet_fMet_1.svg", width = 10, height = 15)

ggarrange(list2[[13]], list2[[14]],list2[[15]], list2[[16]],
          list2[[17]], list2[[18]],list2[[19]], list2[[20]],
          list2[[21]],
          nrow = 3, ncol = 3, common.legend = TRUE,
          labels = c(LETTERS[13:21]))
ggsave("results/corr_diet_fMet_2.pdf", width = 10, height = 11)
ggsave("results/corr_diet_fMet_2.svg", width = 10, height = 11)

