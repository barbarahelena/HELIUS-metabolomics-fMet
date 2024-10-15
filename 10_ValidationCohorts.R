## Linear models in EPIC-Norfolk
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "haven", "dplyr", "tidyverse", "corrplot", "Hmisc", "igraph", 
              "ggraph", "network", "naniar", "ggpubr", "ggsci", "gridExtra",
              "RColorBrewer", "broom", "patchwork")
pacman::p_load(packages, character.only = TRUE)

## Opening file with results Epic-Norfolk
epic <- import("data/validation_epicnorfolk.xlsx") %>% 
    filter(effect == "M02829") %>% 
    rename("model" = "model*")

## Calculate 95%CI
epic <- epic %>% 
    mutate(
        margin = qt(0.975,df=10142-1)*SE,
        lower = Estimate - margin,
        upper = Estimate + margin,
        outcome = str_extract(model, "(?<=_)(diastol|systol)(?=_)"),
        model = str_extract(model, "[0-9]")
        ) %>% 
    filter(model %in% c(1:3)) %>% 
    mutate(model = forcats::fct_recode(model, "Unadjusted" = "1",
                                    "Age + sex" = "2",
                                    "+ eGFR" = "3",
                                    # "+ BMI" = "4",
                                    # "+ fasting time" = "5"
                                    ),
            outcome = forcats::fct_recode(outcome, "systolic BP" = "systol",
                                                    "diastolic BP" = "diastol"),
           sigp = case_when(`p-value` < 0.05 ~ "p<0.05", .default = "not sig")
    )

(plepic <- ggplot(epic, aes(x=outcome, y=Estimate, color=model, shape=sigp)) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.5)) +
        scale_shape_manual(values = c(19))+
        scale_y_continuous(breaks = c(-2:3))+
        geom_errorbar(aes(ymin=lower,ymax=upper,width=.3), position=position_dodge(-0.5)) +
        theme_Publication()+
        theme(legend.position = "bottom")+
        labs(title = "Linear regression: fMet and BP in EPIC-Norfolk cohort",
             x = "", y = "Difference in BP (mmHg) per SD increase", shape = "", color = "") +
        scale_color_lancet() +
        coord_flip())

ggsave("results/lm_epicnorfolk.pdf", width = 6, height = 3)
ggsave("results/lm_epicnorfolk.svg", width = 6, height = 3)

str(epic)
