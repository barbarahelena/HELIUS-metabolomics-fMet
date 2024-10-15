## Functions regression

compared_to_permuted_reg <- function(path_true_outcome, path_permuted_outcome){
    theme_Publication <- function(base_size=14, base_family="sans") {
        library(grid)
        library(ggthemes)
        library(stringr)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.2), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    
    # load permuted results
    path <- file.path(path_permuted_outcome,'model_results_per_iteration.txt')
    df <- rio::import(path)
    df$`Explained Variance` <- df$`Explained Variance` * 100
    n_iter <- nrow(df)
    path_true <- file.path(path_true_outcome,'model_results_per_iteration.txt')
    df_true <- rio::import(path_true)
    df_true$`Explained Variance` <- df_true$`Explained Variance` * 100
    true_result_explained_variance <- median(df_true$`Explained Variance`)
    df$Model <- 'Permuted'
    df_true$Model <- 'True'
    dff <- rbind(df, df_true)
    
    pl <- ggdensity(df$`Explained Variance`, 
                    main = paste("Density plot of", n_iter,'permuted iterations'),
                    xlab = "Explained Variance [%]", fill=pal_lancet()(1), alpha = 0.7)+
        geom_vline(xintercept = 0, linetype=1, size=1, color="black")
    #print(pl)
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_permuted_iterations.pdf', width = 6, height = 4)
    
    # test if the permuted iterations is significantly smaller than the mean of the true model
    t_test <- t.test(df$`Explained Variance`, mu = true_result_explained_variance, alternative = "less")
    p <- format(x=t_test$p.value, scientific=T, digits = 3)
    pl <- ggplot(df, aes(x=`Explained Variance`))+
        geom_density(fill=pal_lancet()(1), color='black', alpha = 0.7)+
        theme_Publication()+
        labs(x='Explained variance (%)', title=paste0("True model compared to permuted outcome\np-value = ",p))+
        geom_vline(xintercept = true_result_explained_variance, linetype=2, size=1, color='black')+
        geom_vline(xintercept = median(df$`Explained Variance`), linetype=2, size=1, color='black')+
        xlim(c(min(df$`Explained Variance`), true_result_explained_variance+1))
    #print(pl)
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted_with_t_test.pdf', width = 6, height = 6)
    
    pl <- ggplot(dff, aes(x=`Explained Variance`, fill=Model))+
        scale_fill_manual(values=pal_lancet()(2))+
        geom_density(color='black', alpha = 0.7)+
        theme_Publication()+
        labs(x='Explained variance (%)', title="True model compared to permuted outcome")+
        geom_vline(xintercept = 0, linetype=1, size=1, color='black')+
        geom_vline(xintercept = median(df_true$`Explained Variance`), linetype=2, size=1, color='black')+
        geom_vline(xintercept = median(df$`Explained Variance`), linetype=2, size=1, color='black')
    #print(pl)
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted_with_all_true_iterations.pdf', width = 6, height = 6)
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted_with_all_true_iterations.svg', width = 6, height = 6)
    
    # test if  true iterations are significantly larger than zero
    t_test <- t.test(df_true$`Explained Variance`, mu = 0, alternative = "greater")
    p <- format(x=t_test$p.value, scientific=T, digits = 3)
    pl <- ggplot(df_true, aes(x=`Explained Variance`))+
        geom_density(fill=pal_lancet()(2)[2], color='black', alpha = 0.7)+
        theme_Publication()+
        xlab('Explained variance [%]')+
        geom_vline(xintercept = true_result_explained_variance, linetype=2, size=1, color='black')+
        geom_vline(xintercept = 0, linetype=1, size=1, color='black')+
        ggtitle(paste0("True model compared to zero\np-value = ",p))
    #print(pl)
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_zero_with_t_test.pdf', width = 6, height = 6)
} 

plot_features_tests_reg <- function(input_path, output_path, top_n=10, outcome_name='Outcome', x_lab='Levels', residuals=FALSE){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    yfile <- ifelse(residuals==TRUE, 'y_reg.txt', 'y_unscaled.txt')
    
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\n
        Showing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, yfile))
    dd$y <- as.numeric(as.character(y$V1))
    
    for(j in 1:(ncol(dd)-1)){
        df <- dd[, c(j, ncol(dd))]
        names(df)[1] <- 'Feature'
        cc <- cor.test(df$Feature, df$y, method = 'spearman')
        pl <- ggplot(df, aes(x=Feature, y=y))+
            geom_point(color='black', fill=pal_lancet()(2)[1], shape=21)+
            geom_smooth(method='lm', color=pal_lancet()(2)[2])+
            scale_x_log10() +
            theme_Publication() +
            theme(legend.position = 'right') +
            labs(x=x_lab, y=outcome_name, title=paste(paste(strwrap(r$FeatName[j], width = 45), collapse = "\n"), "\nrho =", round(cc$estimate, digits = 3), ' p =', format(cc$p.value, digits = 2, scientific = T)))+
            #scale_fill_gradient2_tableau(trans='reverse')+
            theme(legend.key.height = unit(1, "cm"))+
            labs(fill=outcome_name)
        fname <- r$FeatName[j]
        cat(j, fname, '\n')
        fname <- str_replace_all(fname, "[*\";,:/\\\\ ]","_")
        ggsave(pl, path = plot_path, filename = paste0(j, '_',fname, '.pdf'), device = 'pdf', width = 5, height = 5)
        ggsave(pl, path = plot_path, filename = paste0(j, '_',fname, '.svg'), device = 'svg', width = 5, height = 5)
    }
}

plot_features_top_n_reg <- function(input_path, output_path, top_n=20, nrow=4, outcome_name='Outcome', x_lab='Levels', residuals=FALSE){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold", family = 'Helvetica',
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    yfile <- ifelse(residuals==TRUE, 'y_reg.txt', 'y_unscaled.txt')
    
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, yfile))
    dd$y <- y$V1
    dd$y <- as.numeric(as.character(dd$y))
    df <- dd %>% pivot_longer(-y, names_to = 'features', values_to = 'values')
    df$features <- factor(df$features, r$FeatName)
    pl <- ggplot(df, aes(x=values, y=y))+
        geom_point(fill=pal_lancet()(1), color='black', shape=21)+
        geom_smooth(method='lm', color=pal_lancet()(2)[2])+
        scale_x_log10()+
        theme_minimal()+
        theme_Publication()+
        theme(legend.position = 'right')+
        labs(x=x_lab, y=outcome_name)+
        #scale_fill_gradient2_tableau()+
        theme(legend.key.height = unit(1, "cm"))+
        facet_wrap(~ features, nrow=nrow, scales = 'free')
    pl
    ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features_scatterplots.pdf'), 
           device = 'pdf', width = 22, height = 16)
}

plot_feature_importance_microbiome <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    cols <- list(low = '#ECE7F2',
                 mid = '#0570B0',
                 high = '#034E7B')
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    mp <- mean(r$RelFeatImp)
    r$Tax <- factor(r$FeatName, levels = rev(r$FeatName))
    r <- r[1:top_n, ]
    mp <- mean(r$RelFeatImp)
    pl <- ggplot(data=r, aes(y=RelFeatImp, x=Tax, fill=RelFeatImp)) + 
        theme_Publication()+
        scale_fill_gradient2(low=cols$low, mid = cols$mid, high=cols$high, space='Lab', 
                             name="", midpoint = 65)+
        geom_bar(stat="identity")+
        coord_flip() +
        labs(y="Relative Importance %", x='', title = 'Metabolite relative importance') +
        theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10))+
        theme(legend.key.size= unit(0.5, "cm"))+
        theme(legend.position = 'none')
    ggsave(pl, path = path_true, filename = str_c('plot_Feature_Importance_Metabolites.pdf'), device = 'pdf', width = 8, height = 7)
    ggsave(pl, path = path_true, filename = str_c('plot_Feature_Importance_Metabolites.svg'), device = 'svg', width = 8, height = 7)
}

plot_feature_importance_class_microbiome <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
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
    
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    a <- rio::import('../data/Info_plasma_metabolites_b.xlsx')
    a <- a %>% select(FeatName = BIOCHEMICAL, sup = `SUPER PATHWAY`) %>% 
        mutate(sup = as.factor(sup))
    summ <- a %>% 
        group_by(sup) %>% 
        dplyr::summarise(count=n()) %>% 
        filter(sup!='Xenobiotics') %>%
        arrange(desc(count)) %>% 
        mutate(sup = reorder(sup, count)) %>% 
        droplevels()
    ra <- left_join(r, a, by='FeatName')
    mp <- mean(r$RelFeatImp)
    ra$Tax <- factor(ra$FeatName, levels = rev(ra$FeatName))
    ra$sup <- factor(ra$sup, rev(levels(summ$sup)))
    ra <- ra[1:top_n, ]
    mp <- mean(ra$RelFeatImp)
    colpal <- c("Lipid"="#00468BFF", "Amino Acid"="#ED0000FF", "Nucleotide"="#42B540FF",
                "Cofactors and Vitamins"="#0099B4FF", "Carbohydrate"="#925E9FFF", "Peptide"="#FDAF91FF", "Partially Characterized Molecules"="#AD002AFF", "Energy"="#ADB6B6FF")
    pl <- ggplot(data=ra, aes(y=RelFeatImp, x=Tax, fill=sup)) + 
        theme_Publication() +
        geom_bar(stat="identity", alpha=0.8) +
        scale_fill_manual(values = colpal)+
        coord_flip() +
        labs(y = 'Relative Importance %', x='', 
             title = 'Metabolite relative importance',
             fill = '') +
        theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10)) +
        theme(legend.key.size= unit(0.5, "cm")) +
        theme(legend.position = 'bottom', legend.justification = 'center')
    pl
    ggsave(pl, path = path_true, filename = 'plot_Feature_Importance_Metabolites_color.pdf', 
            device = 'pdf', width = 10, height = 7)
    ggsave(pl, path = path_true, filename = 'plot_Feature_Importance_Metabolites_color.svg', 
            device = 'svg', width = 10, height = 7)
}


## Functions classification

compared_to_permuted_class <- function(path_true_outcome, path_permuted_outcome){
    theme_Publication <- function(base_size=14, base_family="sans") {
        library(grid)
        library(ggthemes)
        library(stringr)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.2), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    
        # load permuted results
    path <- file.path(path_permuted_outcome,'model_results_per_iteration.txt')
    df_permuted <- rio::import(path)
    n_iter <- nrow(df_permuted)
    path_true <- file.path(path_true_outcome,'model_results_per_iteration.txt')
    df_true <- rio::import(path_true)
    true_result_AUC <- median(df_true$ROC_AUC_scores)
    df_permuted$Model <- 'Permuted'
    df_true$Model <- 'True'
    dff <- rbind(df_permuted, df_true)
    
    # test if the permuted iterations is significantly smaller than the mean of the true model
    mwu_test <- wilcox.test(df_permuted$ROC_AUC_scores, mu = true_result_AUC, alternative = "less", paired = F)
    p <- format(x = mwu_test$p.value, scientific=T, digits = 3)
    
    # plot all iterations 
    pl <- ggplot(df_permuted, aes(x=ROC_AUC_scores))+
        geom_density(fill=pal_lancet()(2)[2], alpha = 0.5)+
        theme_Publication()+
        xlab('AUC')+
        geom_vline(xintercept = true_result_AUC, linetype=1, size=1, color=pal_lancet()(1))+
        geom_vline(xintercept = median(df_permuted$ROC_AUC_scores), linetype=1, size=1, color='black')+
        ggtitle(paste0("True model median compared to \nall permuted outcome iterations\np-value = ", p))+
        xlim(c(min(df_permuted$ROC_AUC_scores), 1))
    pl
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted_with_t_test.pdf', width = 6, height = 6)
    
    # plot all iterations of true model and all iterations of permuted model
    pl <- ggplot(dff, aes(x=ROC_AUC_scores, fill=Model))+
        geom_density(alpha=0.5)+
        scale_fill_manual(values = c(pal_lancet()(2)[2], pal_lancet()(1)))+
        theme_Publication()+
        xlab('AUC')+
        geom_vline(xintercept = median(df_true$ROC_AUC_scores), linetype=1, size=1, color='black')+
        geom_vline(xintercept = median(df_permuted$ROC_AUC_scores), linetype=1, size=1, color='black')+
        ggtitle("True model compared to permuted outcome")+
        geom_vline(xintercept = 1, linetype=1)+
        geom_vline(xintercept = 0.50, linetype=2, color='black')
    pl
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_permuted.pdf', width = 6, height = 6)
    
    # test if  true iterations are significantly larger than zero
    mwu_test <- wilcox.test(df_true$ROC_AUC_scores, mu = 0, alternative = "greater", paired = F)
    p <- format(x=mwu_test$p.value, scientific=T, digits = 3)
    pl <- ggplot(df_true, aes(x=ROC_AUC_scores))+
        geom_density(fill=pal_lancet()(2)[1], alpha = 0.5)+
        theme_Publication()+
        xlab('AUC')+
        geom_vline(xintercept = true_result_AUC, linetype=1, size=1, color='black')+
        geom_vline(xintercept = 0.50, linetype=2, size=1, color='black')+
        ggtitle(paste0("True model compared to zero\np-value = ",p))
    pl
    ggsave(pl, path = path_true_outcome, filename = 'density_plot_true_versus_coin_flip.pdf', width = 6, height = 6)
} 

plot_feature_importance_class <- function(path_true, top_n){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5,
                                              family = 'Helvetica'),
                    text = element_text(family = 'Helvetica'),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
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
    cols <- list(low = '#ECE7F2',
                 mid = '#0570B0',
                 high = '#034E7B')
    r <- rio::import(file.path(path_true, 'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    mp <- mean(r$RelFeatImp)
    r$Tax <- factor(r$FeatName, levels = rev(r$FeatName))
    r <- r[1:top_n, ]
    mp <- mean(r$RelFeatImp)
    pl <- ggplot(data=r, aes(y=RelFeatImp, x=Tax, fill=RelFeatImp)) + 
        theme_Publication()+
        scale_fill_gradient2(low=cols$low, mid = cols$mid, high=cols$high, space='Lab', name="",
                             midpoint = 50)+
        geom_bar(stat="identity")+
        coord_flip() +
        ylab("Metabolite Relative Importance (%)")+
        xlab("") +
        theme(axis.text.x = element_text(size=12)) + 
        theme(axis.text.y = element_text(size=10))+
        theme(legend.key.size= unit(0.5, "cm"))+
        theme(legend.position = 'right')
    ggsave(path = path_true, filename = 'plot_Feature_Importance.pdf', device = 'pdf', width = 8, height = 7)
}

plot_features_tests_class <- function(input_path, output_path, top_n=10){
    theme_Publication <- function(base_size=11, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    ## edit depending on outcome
    dd$y <- factor(ifelse(dd$y==0, 'Female','Male'))
    comps <- list(c('Female','Male'))
    cl <-  c(pal_lancet()(1), pal_lancet()(2)[2])
    for(j in 1:(ncol(dd)-1)){
        df <- dd[, c(j, ncol(dd))]
        names(df)[1] <- 'Feature'
        df$y <- factor(df$y)
        library(ggpubr)
        pl <- ggplot(df, aes(x=y, y=Feature))+
            geom_violin(aes(fill=y)) +
            theme_Publication()+
            theme(legend.position = 'none')+
            xlab('Group')+
            ylab('Metabolite concentration')+
            ggtitle(r$FeatName[j])+
            stat_compare_means(comparisons = comps, paired = F)+
            scale_color_manual(values = cl)+
            scale_fill_manual(values = cl)
        fname <- r$FeatName[j]
        cat(j, fname, '\n')
        fname <- str_replace_all(fname, "[*\";,:/\\\\ ]","_")
        
        ggsave(pl, path = plot_path, filename = paste0(j, '_',fname, '.pdf'), device = 'pdf', width = 5, height = 5)
    }
}

plot_features_tests_top <- function(input_path, output_path, top_n=20, nrow=4){
    theme_Publication <- function(base_size=11, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(1)),
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
    plot_path <- file.path(output_path, 'plots')
    dir.create(plot_path)
    r <- rio::import(file.path(output_path,'feature_importance.txt'))
    r <- r %>% arrange(-RelFeatImp)
    input_data <- rio::import(file.path(input_path, 'X_data.txt'))
    feature_names <- read.csv(file.path(input_path, 'feat_ids.txt'), sep = '\t', header = F)
    names(input_data) <- feature_names$V1
    if(top_n > ncol(input_data)){
        cat('\n\nRequested no. of features is higher than total number of features in model.\nShowing all features in model.\n\n')
        top_n <- ncol(input_data)
    }
    features_tk <- r$FeatName[1:top_n]
    features_tk <- features_tk[! features_tk %in% c('random_variable1', 'random_variable2')]
    dd <- input_data %>% dplyr::select(any_of(features_tk))
    y <- rio::import(file.path(input_path, 'y_binary.txt'))
    dd$y <- y$V1
    ## edit depending on outcome
    ## edit depending on outcome
    dd$y <- factor(ifelse(dd$y==0, 'Female','Male'))
    comps <- list(c('Female','Male'))
    cl <-  c(pal_lancet()(2)[2], pal_lancet()(2)[1])
    df <- dd %>% pivot_longer(-y, names_to = 'features', values_to = 'values')
    df$y <- factor(df$y)
    library(ggpubr)
    pl <- ggplot(df, aes(x=y, y=values))+
        geom_violin(aes(fill=y), alpha = 0.8)+
        scale_fill_manual(values = cl) +
        theme_Publication()+
        theme(legend.position = 'none')+
        xlab('Group')+
        ylab('Metabolite concentration')+
        stat_compare_means(comparisons = comps, paired = F)+
        facet_wrap(~ features, nrow=nrow, scales = 'free')
    ggsave(pl, path = plot_path, filename = paste0('top_',top_n,'_features_boxplots.pdf'), device = 'pdf', width = 22, height = 16)
}

