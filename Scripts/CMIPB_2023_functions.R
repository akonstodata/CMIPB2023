# CMIPB 2023 functions

# Load datasets, combine different experiments from the training data
combine_exps<-function(df_source, output_request){
  
  metaDf <- data.frame(df_source[["subject_specimen"]])
  metaDf["age_at_boost"] <- as.numeric(round(difftime(metaDf$date_of_boost, metaDf$year_of_birth,units="weeks")/52, 2))
  metaDf <- metaDf[, output_request$META_COLS] %>%
    data.frame()
  
  abtiterDf <- df_source[["abtiter_wide"]]$normalized_data %>%
    t() %>%
    data.frame() 
  
  abtiterDf$specimen_id <- as.numeric(rownames(abtiterDf))
  abtiterDf <- data.frame(abtiterDf[, c("specimen_id", output_request$ABTITER_COLS)])
  
  
  rnaDf <- df_source[["pbmc_gene_expression"]]$raw_data %>%
    t() %>%
    data.frame() 
  
  rnaDf$specimen_id <- as.numeric(rownames(rnaDf))
  tasks_seq <- c('ENSG00000277632')
  for (i in 1:length(tasks_seq)){
    rnaDf <- data.frame(rnaDf %>% rename_at(vars(starts_with(tasks_seq[i])), ~output_request$RNA_COLS[i]))
  }
  rnaDf <- data.frame(rnaDf[, c("specimen_id", output_request$RNA_COLS)])
  
  
  cellDf <- df_source[["pbmc_cell_frequency"]]$normalized_data %>%
    t() %>%
    data.frame()  
  
  cellDf$specimen_id <- as.numeric(rownames(cellDf))
  cellDf <- data.frame(cellDf[, c("specimen_id", output_request$CELL_COLS)])
  
  list_df <- list(metaDf, cellDf, abtiterDf, rnaDf)
  df_merge <- list_df %>% reduce(full_join, by="specimen_id")
  df_merge <- df_merge[df_merge$timepoint %in% output_request$TIMEPOINTS, ]
  
  df_pivot <- df_merge[, names(df_merge)!="specimen_id"] %>%
    pivot_wider(id_cols=c("subject_id", "dataset", "biological_sex",
                          "infancy_vac", "age_at_boost"),
                names_from = timepoint,
                values_from = all_of(c(output_request$CELL_COLS, output_request$RNA_COLS, output_request$ABTITER_COLS)),
                names_sep = "_D")
  
  df_pivot <- df_pivot[df_pivot$dataset %in% output_request$DATASET, ]
  df_pivot <- data.frame(df_pivot %>%
                           mutate(across(everything(),  ~ case_when(.x >=0 ~ .x))))
}

## Dataset fold change and ranking
fold_change <- function(df_pivot){
  targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
  targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")
  
  fc_cols <- paste(targetY, "FC", sep="_")
  
  ranked_cols <- paste(c("age_at_boost", targetX, targetY, fc_cols), "Rank", sep="_")
  
  df <- df_pivot[, c("subject_id", "dataset", DEMOGRAPHY_COLS, targetX, targetY)]
  
  rankingFunction <- function(x) {
    as.numeric(rank(-x, ties.method = "min", na.last = "keep"))
  }
  
  targetX <- c("Monocytes_D0", "CCL3_D0", "IgG_PT_D0")
  targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")
  
  df[,"Monocytes_D1_FC"] <- df[, "Monocytes_D1"] / df[, "Monocytes_D0"]
  df[,"CCL3_D3_FC"] <- df[, "CCL3_D3"] / df[, "CCL3_D0"]
  df[,"IgG_PT_D14_FC"] <- df[, "IgG_PT_D14"] / df[, "IgG_PT_D0"]
  
  df[, ranked_cols] <- apply(df[, c("age_at_boost", targetX, targetY, fc_cols)],
                             2, rankingFunction)
  df <- data.frame(df)
  
  
}

## Load datasets, combine different experiments from test data
combine_exps_test<-function(df_source){
  DATASET <- c("2022_dataset")
  TIMEPOINTS <- c(0)
  
  META_COLS <- c("specimen_id", "subject_id", "timepoint", "dataset", 
                 "biological_sex", "infancy_vac", "age_at_boost")
  ABTITER_COLS <- c("IgG_PT")
  RNA_COLS <- c("CCL3")
  CELL_COLS <- c("Monocytes")
  
  DEMOGRAPHY_COLS <- c("age_at_boost", "biological_sex", "infancy_vac")
  
  BASE_COLS <- c("Monocytes_D0", "IgG_PT_D0", "CCL3_D0")
  
  
  metaDf <- data.frame(df_source[["subject_specimen"]])
  metaDf["age_at_boost"] <- as.numeric(round(difftime(metaDf$date_of_boost, metaDf$year_of_birth,units="weeks")/52, 2))
  metaDf <- metaDf[, META_COLS] %>%
    data.frame()
  
  abtiterDf <- df_source[["abtiter"]]$processed_similar_to_training %>%
    t() %>%
    data.frame() 
  
  abtiterDf$specimen_id <- as.numeric(rownames(abtiterDf))
  abtiterDf <- data.frame(abtiterDf[, c("specimen_id", ABTITER_COLS)])
  
  
  rnaDf <- df_source[["pbmc_gene_expression"]]$processed_similar_to_training %>%
    t() %>%
    data.frame() 
  
  rnaDf$specimen_id <- as.numeric(rownames(rnaDf))
  tasks_seq <- c('ENSG00000277632')
  for (i in 1:length(tasks_seq)){
    rnaDf <- data.frame(rnaDf %>% rename_at(vars(starts_with(tasks_seq[i])), ~RNA_COLS[i]))
  }
  rnaDf <- data.frame(rnaDf[, c("specimen_id", RNA_COLS)])
  
  
  cellDf <- df_source[["pbmc_cell_frequency"]]$processed_similar_to_training %>%
    t() %>%
    data.frame()  
  
  cellDf$specimen_id <- as.numeric(rownames(cellDf))
  cellDf <- data.frame(cellDf[, c("specimen_id", CELL_COLS)])
  
  list_df <- list(metaDf, cellDf, abtiterDf, rnaDf)
  df_merge <- list_df %>% reduce(full_join, by="specimen_id")
  df_merge <- df_merge[df_merge$timepoint %in% TIMEPOINTS, ]
  
  df_pivot <- df_merge[, names(df_merge)!="specimen_id"] %>%
    pivot_wider(id_cols=c("subject_id", "dataset", "biological_sex",
                          "infancy_vac", "age_at_boost"),
                names_from = timepoint,
                values_from = all_of(c(CELL_COLS, RNA_COLS, ABTITER_COLS)),
                names_sep = "_D")
  
  df_pivot <- df_pivot[df_pivot$dataset %in% DATASET, ]
  df_pivot <- data.frame(df_pivot %>%
                           mutate(across(everything(),  ~ case_when(.x >=0 ~ .x))))
}

# Analyze correlation between predicted and true values in leave-one-out cross-validation assay

pred_cor_calc <- function(gs_filt_use, train_tasks){
  pred_cor<-data.frame(matrix(nrow=ncol(train_tasks),ncol=1))
  rownames(pred_cor)<-colnames(train_tasks)
  colnames(pred_cor)<-c('cor.pred.true')
  # Loop through all tasks
  for (i in 1:ncol(train_tasks)){
    print(i)
    all_preds<-c()
    all_true<-c()
    set.seed(1)
    
    x_out <-data.frame(train_tasks[,i])
    rownames(x_out)<-rownames(train_tasks)
    
    x_out$temp<-'temp'
    names(x_out)<-c('Y','temp')
    x_out_r<-x_out
    #x_out_r<-na.omit(x_out)
    #print(x_out_r)
    
    #row_int<-intersect(rownames(x_out_r),rownames(gs_full_filt2))
    x_out_r<-x_out_r[rownames(gs_full_filt2),]
    #gs_full_filt<-gs_full[rownames(x_out_r),]
    
    train = 1:nrow(gs_full_filt2)
    
    # train all leave-one-out models
    for (j in 1:nrow(gs_full_filt2)){
      train = 1:nrow(gs_full_filt2)
      train = train[-c(j)]
      
      # create lasso model
      cvfit_out<-cv.glmnet(x=as.matrix(gs_full_filt2[train,]), x_out_r[train,'Y'], family='gaussian',
                           alpha=1, nfolds=nrow(gs_full_filt2[train,]-1)) 
      preds<-predict(cvfit_out,newx=as.matrix(data.frame(gs_full_filt2[-train,])),s='lambda.min')
      all_preds<-c(all_preds,preds)
      all_true<-c(all_true,x_out_r[-train,'Y']) 
    }
    
    pred_cor[i,'cor.pred.true']<-cor.test(all_preds,all_true,method='spearman')$estimate[[1]] # Can also use spearman correlation
  }
  print(pred_cor)
}

# Returns top features for a given input dataset and prediction tasks
top_feat_model <- function(gs_filt_use, train_tasks){
  
  all_models_coef<-vector(mode='list',length=ncol(train_tasks))
  all_models_names<-vector(mode='list',length=ncol(train_tasks))
  all_models<-vector(mode='list',length=ncol(train_tasks))
  
  pred_cor<-data.frame(matrix(nrow=ncol(train_tasks),ncol=1))
  rownames(pred_cor)<-colnames(train_tasks)
  colnames(pred_cor)<-c('cor.pred.true')
  # Loop through all tasks
  for (i in 1:ncol(train_tasks)){
    set.seed(1)
    all_preds<-c()
    all_true<-c()
    set.seed(1)
    
    x_out <-data.frame(train_tasks[,i])
    rownames(x_out)<-rownames(train_tasks)
    
    x_out$temp<-'temp'
    names(x_out)<-c('Y','temp')
    x_out_r<-x_out
    x_out_r<-x_out_r[rownames(gs_full_filt2),]
    
    # create lasso model
    cvfit_out<-cv.glmnet(x=as.matrix(gs_full_filt2), x_out_r[,'Y'], family='gaussian',
                         alpha=1, nfolds=nrow(gs_full_filt2-1)) 
    all_models_coef[i]=list(coef(cvfit_out, s = 'lambda.min')[coef(cvfit_out, s = 'lambda.min')[,1]!= 0])
    all_models_names[i]=list(rownames(coef(cvfit_out, s = 'lambda.min'))[coef(cvfit_out, s = 'lambda.min')[,1]!= 0])
  }
  
  names(all_models_coef)<-colnames(train_tasks)
  names(all_models_names)<-colnames(train_tasks)
  
  for (i in 1:ncol(train_tasks)){
    all_models[[i]] = data.frame(cbind(all_models_names[[i]],all_models_coef[[i]]))
    colnames(all_models[[i]])<-c("Variable","Coefficient")
    all_models[[i]]$Coefficient<-as.numeric(all_models[[i]]$Coefficient)
    all_models[[i]]$Coefficient=round(all_models[[i]]$Coefficient,4)
    all_models[[i]]<-all_models[[i]] %>% arrange(desc(abs(Coefficient)))
  }
  names(all_models)<-colnames(train_tasks)
  
  
  return(all_models)
  
}


# Prediction for test data
pred_cor_calc_test <- function(task,model_train,new_data){
  val_predict<-data.frame(matrix(nrow = nrow(new_data), ncol=1))
  colnames(val_predict)<-c(task)
  all_preds<-c()
  set.seed(1)
  x_out <-data.frame(train_tasks[,task])
  rownames(x_out)<-rownames(train_tasks)
  x_out$temp<-'temp'
  names(x_out)<-c('Y','temp')
  x_out_r<-x_out
  x_out_r<-x_out_r[rownames(model_train),]
  
  # create lasso model using training data
  cvfit_out<-cv.glmnet(x=as.matrix(model_train), x_out_r$Y, family='gaussian',
                       alpha=1,nfolds=nrow(model_train),type.measure="mse") 
  
  # make predictions on new data
  preds<-data.frame(predict(cvfit_out,newx=as.matrix(new_data),s='lambda.min'))
  val_predict[,1]<-preds
}

# Generate rankings
rankingFunction <- function(x) {
  as.numeric(rank(-x, ties.method = "min", na.last = "keep"))
}

targetX <- c("Monocytes_D0", "CCL3_D0")
targetY <- c("Monocytes_D1","CCL3_D3", "IgG_PT_D14")
#ranked_cols <- paste(c(targetX,targetY, fc_cols), "Rank", sep="_")
ranked_cols_base<-paste(targetX,"Rank",sep="_")

#df[, ranked_cols] <- apply(df[, c("age_at_boost", targetX, targetY, fc_cols)],
#                           2, rankingFunction)