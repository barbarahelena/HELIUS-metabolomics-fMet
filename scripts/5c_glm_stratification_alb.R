## Linear models stratified for ethnicity
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "haven", "dplyr", "tidyverse", "corrplot", "Hmisc", "igraph", 
              "ggraph", "network", "naniar", "ggpubr", "ggsci", "gridExtra",
              "RColorBrewer", "broom", "patchwork")
pacman::p_load(packages, character.only = TRUE)

## Functions
lm_prep <- function(df, metabolites, helius){
    df <- df %>% select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% slice(1:20)
    df_met <- metabolites %>% inner_join(., df, by='Metabolite') %>% 
        arrange(-RelFeatImp)
    rownames(df_met) <- df_met$Metabolite
    df_met$Metabolite <- NULL
    df_met$RelFeatImp <- NULL
    df_met <- as.data.frame(t(as.matrix(df_met)))
    head(df_met)
    df_met$ID <- as.integer(rownames(df_met))
    helius$ID <- helius$Heliusnr
    helius_sub <- left_join(helius, df_met, by='ID')
    return(helius_sub)
}

afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
afronden3 <- function(x) return(as.numeric(format(round(x, 3),3)))

## Opening HELIUS file, metabolite data and best predictor files
helius <- readRDS("data/heliusdf.RDS")
helius$Ethnicity <- case_when(helius$Ethnicity == "NL" ~ "Dutch",
                              helius$Ethnicity == "Hind" ~ "South-Asian Surinamese",
                              helius$Ethnicity == "Creools" ~ "African Surinamese",
                              helius$Ethnicity == "Ghanees" ~ "Ghanaian")
summary <- rio::import('data/Info_plasma_metabolites_b.xlsx')
infomet <- summary %>% select(metabolite=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
metabolites <- rio::import('data/HELIUS_EDTA_plasma_metabolites.xlsx')
best_sbp <- rio::import('All/SBP/output_XGB_reg_SBP_100_iterations_y_scaled_2020_11_06__16-41-19/feature_importance.txt')
best_dbp <- rio::import('All/DBP/output_XGB_reg_DBP_100_iterations_y_scaled_2020_11_06__16-45-06/feature_importance.txt')

## Preparation for models
helius_sbp <- lm_prep(best_sbp, metabolites, helius)
helius_dbp <- lm_prep(best_dbp, metabolites, helius)

## Regression models
### SBP models stratified for ethnic group
res_sbp <- c()
for (i in c(38:57)) {
    helius_sbp$met <- NULL    
    helius_sbp$met <- helius_sbp[,i]
    # Other covariates from other models: age, sex, ckdepi, smoking
    m0 <- lm(SBP ~ scale(met), data = helius_sbp %>% filter(group == "Controls"))
    m1 <- lm(SBP ~ scale(met), data = helius_sbp %>% filter(group != "Controls"))
    
    metname <- colnames(helius_sbp)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    
    resRow <- cbind(metname, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value)
    colnames(resRow) <- c("Metabolite", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p")
    res_sbp <- rbind(res_sbp, resRow)
    
}

ressyss <- as.data.frame(res_sbp)
ressyss2 <- ressyss %>% 
    mutate_at(c(2:9), as.character) %>% 
    mutate_at(c(2:9), as.numeric) %>% 
    mutate_at(c(2:4, 6:8), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr'),
    ) %>% 
    mutate_at(c(5,9), afronden3)
openxlsx::write.xlsx(ressys2, "results/lm_sbp_alb.xlsx")

ressys3 <- ressyss2 %>% 
    pivot_longer(c(2:11), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1"), 
                          labels = c("Controls", "Albuminuria")),
           Metabolite = factor(Metabolite, levels = colnames(sbp_met)[1:20]),
           Metabolite = fct_rev(Metabolite),
           sigp = case_when(p < 0.05 ~ "p<0.05", 
                            .default = "not sig"),
           sigp = as.factor(sigp))

(plsbp <- ggplot(ressys3, aes(x=Metabolite, y=est, color=model, shape=sigp)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.5)) +
        scale_shape_manual(values = c(21,19))+
        scale_y_continuous(breaks = c(-2:12))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.5)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Linear regression: top predictors for SBP",
             x = "", y = "Difference in SBP (mmHg) per SD increase", shape = "", color = "") +
        scale_color_lancet() +
        coord_flip())

ggsave("results/lm_sbp_albuminuria.pdf", width = 9, height = 10)
ggsave("results/lm_sbp_albuminuria.svg", width = 9, height = 10)

### DBP regression stratified for ethnic group
res_dbp <- c()
for (i in c(38:57)) {
    helius_dbp$met <- NULL    
    helius_dbp$met <- helius_dbp[,i]
    # Other covariates from other models: age, sex, ckdepi, smoking
    m0 <- lm(DBP ~ scale(met), data = helius_dbp %>% filter(group == "Controls"))
    m1 <- lm(DBP ~ scale(met), data = helius_dbp %>% filter(group != "Controls"))
    
    metname <- colnames(helius_dbp)[i]
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    resRow <- cbind(metname, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value)
    colnames(resRow) <- c("Metabolite", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p")
    res_dbp <- rbind(res_dbp, resRow)
    
}

resdia <- as.data.frame(res_dbp)
resdia2 <- resdia %>% 
    mutate_at(c(2:9), as.character) %>% 
    mutate_at(c(2:9), as.numeric) %>% 
    mutate_at(c(2:4, 6:8), afronden2) %>% 
    mutate(
        `m0-q` = p.adjust(`m0-p`, 'fdr'),
        `m1-q` = p.adjust(`m1-p`, 'fdr')
    ) %>% 
    mutate_at(c(5,9), afronden3)
openxlsx::write.xlsx(resdia2, "results/lm_dbp_albuminuria.xlsx")

resdia3 <- resdia2 %>% 
    pivot_longer(c(2:11), names_to=c("model", "cat"), 
                 names_prefix="m", 
                 names_sep='-',
                 values_to="value") %>% 
    pivot_wider(names_from = cat, values_from = value) %>% 
    mutate(model = factor(model, levels = c("0", "1"), 
                          labels = c("Controls", "Albuminuria")),
           Metabolite = factor(Metabolite, levels = colnames(dbp_met)[1:20]),
           Metabolite = fct_rev(Metabolite),
           sigp = case_when(p < 0.05 ~ "p<0.05", .default = "not sig"),
           sigp = as.factor(sigp))

(pldbp <- ggplot(resdia3, aes(x=Metabolite, y=est, color=model, shape=sigp)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.5)) +
        scale_shape_manual(values = c(21,19))+
        scale_y_continuous(breaks = c(-5:5))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.5)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Linear regression: top predictors for DBP",
             x = "", y = "Difference in DBP (mmHg) per SD increase", shape = "",
             color = "") +
        scale_color_lancet() +
        coord_flip())

ggsave("results/lm_dbp_albuminuria.pdf", width = 9, height = 10)
ggsave("results/lm_dbp_albuminuria.svg", width = 9, height = 10)

ggarrange(plsbp, pldbp, nrow = 2, labels = LETTERS[1:2], common.legend = TRUE, legend = "bottom")
ggsave("results/lm_tot_albuminuria.svg", width = 9, height = 16)
ggsave("results/lm_tot_albuminuria.pdf", width = 9, height = 16)

## Only for formylmethionine
restot <- rbind(resdia3 %>% mutate(Outcome = "Diastolic BP"), 
                ressys3 %>% mutate(Outcome = "Systolic BP")) %>% 
    filter(Metabolite == "N-formylmethionine")
(pltot <- ggplot(restot, aes(x=Outcome, y=est, color=model, shape=sigp)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.5)) +
        scale_shape_manual(values = c(21,19))+
        scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14))+
        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.5)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Associations fMet with BP: stratified for ethnicity",
             x = "", y = "Difference in BP (mmHg) per SD increase in fMet", shape = "",
             color = "") +
        scale_color_lancet() +
        coord_flip())