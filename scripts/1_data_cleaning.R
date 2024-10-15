## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl


##### Data loading #####
heliusdf <- read_sav("data/201104_HELIUS data Barbara Verhaar.sav")
metabolites <- rio::import('data/HELIUS_EDTA_plasma_metabolites.xlsx')
infomet <- rio::import('data/Info_plasma_metabolites.xlsx')

##### Data wrangling #####
helius <- heliusdf %>% 
    select(ID=Heliusnr, Age=H1_lft, Sex=H1_geslacht, Ethnicity=H1_EtnTotaal, 
           Smoking=H1_Roken, BMI=H1_LO_BMI, WHR=H1_LO_WHR, SBP=H1_LO_GemBPSysZit, 
           DBP=H1_LO_GemBPDiaZit, HT=H1_HT_SelfBPMed, DM=H1_Diabetes_SelfGlucMed,
           DMMed=H1_Diabetesmiddelen, Claudicatio=H1_CI_Rose, PossInf=H1_possINF_Rose, 
           AP=H1_AP_Rose, CVD=H1_CVD_Rose, AntiHT=H1_Antihypertensiva, 
           MDRD=H1_MDRD_eGFR, CKDEPI=H1_CKDEPI_eGFR, CKDStage=H1_CKDEPI_stage, 
           MetSyn=H1_MetSyn_MetabolicSyndrome, TC = H1_Lab_UitslagCHOL, 
           LDL=H1_Lab_uitslagRLDL, HDL=H1_Lab_UitslagHDLS, Trig=H1_Lab_UitslagTRIG,
           FramRisk=H1_Fram_CVD, SCORENL=H1_SCORE_CVDmort_NL,
           ACR_KDIGO=H1_ACR_KDIGO, Microalb=H1_Microalbuminurie, 
           HbA1C=H1_Lab_UitslagIH1C, Kreat=H1_Lab_UitslagKREA_HP, Glucose = H1_Lab_UitslagGLUC) %>% 
    filter(ID %in% sampleList) %>%
    mutate(
        group = case_when(
            ACR_KDIGO==1 ~ 'Controls',
            DM==1 & ACR_KDIGO==2 ~ 'CKD+T2DM',
            DM==0 & ACR_KDIGO==2 ~ 'CKD-T2DM'
        ),
        group = factor(group, levels = c('Controls', 'CKD+T2DM', 'CKD-T2DM')),
        CKD_group = case_when(
            group =='Controls' ~ 'Controls',
            group == 'CKD+T2DM' ~ 'CKD',
            group == 'CKD-T2DM' ~ 'CKD'
        ),
        CKD_group = factor(CKD_group, levels = c('Controls', 'CKD')),
        Sex = as_factor(Sex, levels = c("labels"), ordered = FALSE),
        Sex = fct_recode(Sex, "Male" = "man", "Female" = "vrouw"),
        SexBin = case_when(Sex == "Male" ~ 1,
                           Sex == "Female" ~ 0),
        Ethnicity = as_factor(Ethnicity, levels = c("labels"), ordered = FALSE),
        Smoking = as_factor(Smoking, levels = c("labels"), ordered = FALSE),
        CurrSmoking = fct_recode(Smoking, "Yes" = "Ja", "No" = "Nee, ik heb nooit gerookt", "No" = "Nee, maar vroeger wel"),
        CurrSmoking = fct_infreq(CurrSmoking),
        Glucose = as.numeric(Glucose),
        HT = as_factor(HT, levels = c("labels"), ordered = FALSE),
        DM = as_factor(DM, levels = c("labels"), ordered = FALSE),
        DM = fct_recode(DM, "Diabetes" = "Ja", "No diabetes" = "Nee"),
        DMMed = as_factor(DMMed, levels = c("labels"), ordered = FALSE),
        Claudicatio = as_factor(Claudicatio, levels = c("labels"), ordered = FALSE),
        PossInf = as_factor(PossInf, levels = c("labels"), ordered = FALSE),
        AP = as_factor(AP, levels = c("labels"), ordered = FALSE),
        CVD = as_factor(CVD, levels = c("labels"), ordered = FALSE),
        AntiHT = as_factor(AntiHT, levels = c("labels"), ordered = FALSE),
        CKDStage = as_factor(CKDStage, levels = c("labels"), ordered = FALSE),
        MetSyn=as_factor(MetSyn, levels = c("labels"), ordered = FALSE),
        ACR_KDIGO = as_factor(ACR_KDIGO, levels = c("labels"), ordered = FALSE),
        Microalb = as_factor(Microalb, levels = c("labels"), ordered = FALSE),
        EthnAfr = ifelse(Ethnicity == "Hind" | Ethnicity == "Dutch", "Non-African", "African"),
        EthnBin = case_when(
            EthnAfr == "Non-African" ~ paste("non-black"),
            EthnAfr == "African" ~ paste("black")
        ),
        Age_cat = ifelse(Age<=50, "Young", "Old")
    ) %>% 
    mutate_if(is.numeric, as.numeric) %>% 
    droplevels()

saveRDS(helius, 'data/helius.RDS')

#### Metabolomics data ####
dim(metabolites)
metabolites[1:10, 1]
mat <- metabolites[,-1]
rownames(mat) <- metabolites[,1]

#### Infofile ####
noxeno <- infomet$metabolite[which(infomet$sup!='Xenobiotics')]
mat2 <- mat[which(rownames(mat) %in% noxeno),]
mat2$metabolite <- rownames(mat2)
saveRDS(mat2, 'Met_HELIUS_noXenobiotics.RDS')
