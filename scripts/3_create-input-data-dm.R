# create XGB input files for all datasets

library(dplyr)

rm(list=ls())
setwd('XGB')

# writes input data files for XGB models as tab-delimited 
# subject ids and feature ids are written as separate tab-delimited files
# write X data / predictors
write_data <- function(x, data_path){
    dir.create(data_path)
    x <- as.matrix(x)
    if(any(is.na(x))){
        cat('There are missing values in the input data!\n')
    }
    write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write y / predicted outcome
write_y <- function(x, name_y, data_path){
    dir.create(data_path)
    if(missing(name_y)){
        cat('\n\nYou need to provide a name for the y data file!\n')
    }
    if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
        cat('\nThe file name is not compatible with XGBeast!\n' )
    }
    if(any(is.na(x))){
        cat('\nThere are missing values in the outcome data!\n')
    }
    write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}

## Open HELIUS BP dataframe with BP and residuals data
d <- readRDS('data/helius_bpdm.RDS')
head(d)
d$ID<-paste0('S', d$ID)

## Open RDS file with metabolite data (HELIUS IDs variables; met rownames)
dd <- readxl::read_xlsx('data/HELIUS_EDTA_plasma_metabolites.xlsx')
met <- dd$Metabolite
dd$Metabolite <- NULL
dd2 <- as.data.frame(t(as.matrix(dd)))
dim(dd)
colnames(dd2)
colnames(dd2) <- met
rownames(dd2) <- paste0('S', rownames(dd2))

infomet <- readxl::read_xlsx('data/Info_plasma_metabolites_b.xlsx')
colnames(infomet)
infomet <- infomet %>% select(metabolite = `BIOCHEMICAL`, sup = `SUPER PATHWAY`, sub = `SUB PATHWAY`)
infomet$sup <- as.factor(infomet$sup)
summary(infomet$sup)

noxeno <- infomet$metabolite[which(infomet$sup!='Xenobiotics')]
noxeno[1:5]
colnames(dd2)
ddd <- dd2[,which(colnames(dd2) %in% noxeno)]
rownames(ddd)
ddd <- ddd[which(rownames(ddd) %in% d$ID),]
dim(ddd)

# Put d and dd in same sequence of IDs
d <- d[match(rownames(ddd), d$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(d$ID == rownames(ddd)) # TRUE
d$ID
rownames(ddd)

# make data for machine learning XGB classification models
# make input data for SBP
path <- 'SBP_dm'
dir.create(path)
write_data(ddd, file.path(path, 'input_data'))
y <- as.data.frame(scale(d$SBP))
y_unscaled <- as.data.frame(d$SBP)
y
y_unscaled
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))
write_y(y_unscaled, name_y = 'y_unscaled.txt', file.path(path, 'input_data'))

# make input data for DBP
path <- 'DBP_dm'
dir.create(path)
write_data(ddd, file.path(path, 'input_data'))
y <- as.data.frame(scale(d$DBP))
y_unscaled <- as.data.frame(d$DBP)
y
y_unscaled
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))
write_y(y_unscaled, name_y = 'y_unscaled.txt', file.path(path, 'input_data'))

# make input data for resSBP
path <- 'SBPres_dm'
dir.create(path)
write_data(ddd, file.path(path, 'input_data'))
y <- as.data.frame(d$SBPres)
y
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))


# make input data for resDBP
path <- 'DBPres_dm'
dir.create(path)
write_data(ddd, file.path(path, 'input_data'))
y <- as.data.frame(d$DBPres)
y
write_y(y, name_y = 'y_reg.txt', file.path(path, 'input_data'))

