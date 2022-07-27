#' Runs time-to-event analysis for the subgroup of non-CKD patients, including time to new onset CKD
#' @param siteid currSiteId (site ID of the site as obtained earlier)
#' @param aki_index_recovery aki_index_recovery table from running the run_time_to_event_analysis() function
#' @param aki_index_death aki_index_death table from running the run_time_to_event_analysis() function
#' @param lists Respective lists from running the run_time_to_event_analysis() function
#' @param obfuscation is_obfuscated
#' @param obfuscation_level obfuscation_value
#' @param restrict_model_corr Whether to restrict models
#' @param factor_threshold Threshold to remove a covariate from analysis if the total number of events in a level is less than this number. Default = 5
#' @param use_custom_output Specifies if a custom output directory is to be used
#' @param use_custom_output_dir Custom output directory name
#' 
run_time_to_event_analysis_nonckd <- function(siteid,
                                              aki_index_recovery,aki_index_death,
                                              med_recovery_list, comorbid_recovery_list, demog_recovery_list,earliest_cr_recovery_list,
                                              med_death_list, comorbid_death_list, demog_death_list,earliest_cr_death_list,
                                              obfuscation,obfuscation_level,
                                              restrict_model_corr, factor_threshold = 5,
                                              use_custom_output = FALSE,use_custom_output_dir = "/4ceData/Output") {
  currSiteId <- siteid
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  restrict_models <- restrict_model_corr
  factor_cutoff <- factor_threshold
  if(isTRUE(use_custom_output)) {
    dir.output <- use_custom_output_dir
  } else {
    dir.output <- getProjectOutputDirectory()
  }
  restrict_list <- ""
  if(isTRUE(restrict_models)) {
    cat("\nWe notice that you are keen to restrict the models to certain variables.")
    cat("\nWe are now going to read in the file CustomModelVariables.txt...")
    restrict_list <- scan("Input/CustomModelVariables.txt",what="")
    message(paste("Variables to restrict analyses to :",restrict_list,collapse=" "))
  }
  
  # CKD/Non-CKD subgroup
  model_ckd_subgroup_1 <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2 <- c("age_group","sex","severe","bronchiectasis","copd","rheum","vte")
  model_ckd_subgroup_3 <- c("age_group","sex","severe","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4 <- c("age_group","sex","severe","aki_kdigo_final","acei_arb_preexposure")
  model_ckd_subgroup_1a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_ckd_subgroup_3a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_3b <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral")
  model_ckd_subgroup_4b <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_1c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_ckd_subgroup_3c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_ckd <- list(model_ckd_subgroup_1,model_ckd_subgroup_2,model_ckd_subgroup_3,model_ckd_subgroup_4,model_ckd_subgroup_1a,model_ckd_subgroup_2a,model_ckd_subgroup_3a,model_ckd_subgroup_4a,model_ckd_subgroup_3b,model_ckd_subgroup_4b,model_ckd_subgroup_1c,model_ckd_subgroup_2c,model_ckd_subgroup_3c,model_ckd_subgroup_4c)
  models_ckd_labels <- c("Model_NonCKD_Subgroup_1","Model_NonCKD_Subgroup_2","Model_NonCKD_Subgroup_3","Model_NonCKD_Subgroup_4","Model_NonCKD_Subgroup_1A","Model_NonCKD_Subgroup_2A","Model_NonCKD_Subgroup_3A","Model_NonCKD_Subgroup_4A","Model_NonCKD_Subgroup_3B","Model_NonCKD_Subgroup_4B","Model_NonCKD_Subgroup_1C","Model_NonCKD_Subgroup_2C","Model_NonCKD_Subgroup_3C","Model_NonCKD_Subgroup_4C")
  
  # Generate tables intended to assess new onset of new CKD
  aki_index_nonckd <- aki_index_death[aki_index_death$ckd == 0,]
  aki_index_nonckd_akionly <- aki_index_recovery[aki_index_recovery$ckd == 0,]
  
  # Same workflow
  # 1) Filter for variables with enough levels
  # 2) Filter variables with enough counts on both sides of levels
  # 3) Run models on the filtered list of covariates
  
  # First do this for the new onset CKD tables
  new_ckd_tmp <- aki_index_nonckd[,c("patient_id","new_ckd",demog_death_list,comorbid_death_list,med_death_list,earliest_cr_death_list)] %>% as.data.frame()
  comorbid_new_ckd_list <- comorbid_death_list
  demog_new_ckd_list <- demog_death_list
  med_new_ckd_list <- med_death_list
  earliest_cr_new_ckd_list <- earliest_cr_death_list
  
  if(length(comorbid_new_ckd_list) > 0) {
    comorbid_new_ckd_list_tmp <- vector(mode="list",length=length(comorbid_new_ckd_list))
    for(i in 1:length(comorbid_new_ckd_list)) {
      new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id",comorbid_new_ckd_list[i],"new_ckd")]
      new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(get(comorbid_new_ckd_list[i]),new_ckd)
      new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(new_ckd == 1)
      if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
        message(paste0(c("Including ",comorbid_new_ckd_list[i]," into the comorbid_new_ckd list...")))
        comorbid_new_ckd_list_tmp[i] <- comorbid_new_ckd_list[i]
      }
      rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
    }
    comorbid_new_ckd_list <- unlist(comorbid_new_ckd_list_tmp[lengths(comorbid_new_ckd_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(demog_new_ckd_list) > 0) {
    demog_new_ckd_list_tmp <- vector(mode="list",length=length(demog_new_ckd_list))
    for(i in 1:length(demog_new_ckd_list)) {
      new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id",demog_new_ckd_list[i],"new_ckd")]
      new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(get(demog_new_ckd_list[i]),new_ckd)
      new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(new_ckd == 1)
      if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
        message(paste0(c("Including ",demog_new_ckd_list[i]," into the demog_new_ckd list...")))
        demog_new_ckd_list_tmp[i] <- demog_new_ckd_list[i]
      }
      rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
    }
    demog_new_ckd_list <- unlist(demog_new_ckd_list_tmp[lengths(demog_new_ckd_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(med_new_ckd_list) > 0) {
    med_new_ckd_list_tmp <- vector(mode="list",length=length(med_new_ckd_list))
    if(length(med_new_ckd_list)>0) {
      for(i in 1:length(med_new_ckd_list)) {
        new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id",med_new_ckd_list[i],"new_ckd")]
        new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(get(med_new_ckd_list[i]),new_ckd)
        new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(new_ckd == 1)
        # if(min(new_ckd_tmp3$n) >= factor_cutoff) {
        #     comorbid_new_ckd_list_tmp[i] <- comorbid_new_ckd_list[i]
        # }
        if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
          message(paste0(c("Including ",med_new_ckd_list[i]," into the med_new_ckd list...")))
          med_new_ckd_list_tmp[i] <- med_new_ckd_list[i]
        }
        rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
      }
    }
    med_new_ckd_list <- unlist(med_new_ckd_list_tmp[lengths(med_new_ckd_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(!is.null(earliest_cr_new_ckd_list)) {
    new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id","preadmit_cr_period","new_ckd")]
    new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(preadmit_cr_period,new_ckd)
    new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(new_ckd == 1)
    if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
      message(paste0(c("Including preadmit_cr_period into the earliest_cr_new_ckd list...")))
    } else {
      earliest_cr_new_ckd_list <- NULL # make this null or else you will get this causing errors
    }
    rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
    invisible(gc())
  }
  
  cat("\nFinal factor list for new_ckd (before user customisation): ",paste(c(demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list),collapse=" "))
  
  if(isTRUE(restrict_models)) {
    demog_new_ckd_list <- demog_new_ckd_list[demog_new_ckd_list %in% restrict_list]
    comorbid_new_ckd_list <- comorbid_new_ckd_list[comorbid_new_ckd_list %in% restrict_list]
    med_new_ckd_list <- med_new_ckd_list[med_new_ckd_list %in% restrict_list]
    if(earliest_cr_new_ckd_list %in% restrict_list) {
      earliest_cr_new_ckd_list <- "preadmit_cr_period"
    } else {
      earliest_cr_new_ckd_list <- NULL
    }
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_new_ckd_list,"\nComorbidities:",comorbid_new_ckd_list,"\nMedications:",med_new_ckd_list,"\nPreadmit Cr Period:",earliest_cr_new_ckd_list,sep = " "))
  }
  variable_list_output <- paste(c("Final New Onset CKD variable list:",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list),collapse=" ")
  var_list_new_ckd_nonckd <- c(demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list)
  
  # ===========
  
  
  
  # Run analysis for New Onset CKD
  cat("\nNow proceeding with time-to-event analysis for New Onset CKD (All Non-CKD patients)...")
  cat("\nGenerating Kaplan-Meier curves (time to new onset CKD, all patients)...")
  try({
    cat("\n(a) by AKI occurrence")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ is_aki")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIvsNonAKI_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIvsNonAKI.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(b) by KDIGO AKI severity (no AKI = 0)")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ aki_kdigo_final")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_KDIGO.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity (ALL patients)
  
  try({
    cat("\n(c) by COVID-19 severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ severe")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_Severe.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (Time to new_ckd, all patients)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list), function(x) as.formula(paste('survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_nonckd)}),
    error = function(c) {
      cat("\nError running univariate formulae (univ_models).")
      return(NULL)
    }
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_new_ckd_all <- do.call("rbind",univ_results)
    write.csv(univ_results_new_ckd_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for New CKD onset on Non-CKD patients\n==============\n")
  for(i in 1:14) {
    cat(paste0("\nGenerating ", models_ckd_labels[i], " (Time to New Onset CKD, All Non-CKD patients)..."))
    try({
      new_ckd_model <- c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list)
      new_ckd_model <- new_ckd_model[new_ckd_model %in% models_ckd[[i]]]
      newckdCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ ",paste(new_ckd_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ ",paste(new_ckd_model,collapse="+")))
      coxph_new_ckd <- survival::coxph(newckdCoxPHFormula, data=aki_index_nonckd)
      coxph_new_ckd_summ <- summary(coxph_new_ckd) 
      print(coxph_new_ckd_summ)
      coxph_new_ckd_hr <- cbind(coxph_new_ckd_summ$coefficients,coxph_new_ckd_summ$conf.int)[,-c(6,7)]
      coxph_new_ckd_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_new_ckd_summ$logtest,coxph_new_ckd_summ$sctest,coxph_new_ckd_summ$waldtest))
      coxph_new_ckd_stats2 <- rbind(data.table::as.data.table(coxph_new_ckd_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_new_ckd_summ$rsq,keep.rownames = T))
      write.csv(coxph_new_ckd_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_new_ckd_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_new_ckd_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_new_ckd_plot <- survminer::ggforest(coxph_new_ckd,data=aki_index_nonckd)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_new_ckd_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  # =================================
  # To repeat the above for mortality
  
  new_ckd_tmp <- aki_index_nonckd[,c("patient_id","deceased",demog_death_list,comorbid_death_list,med_death_list,earliest_cr_death_list)] %>% as.data.frame()
  comorbid_new_ckd_list <- comorbid_death_list
  demog_new_ckd_list <- demog_death_list
  med_new_ckd_list <- med_death_list
  earliest_cr_new_ckd_list <- earliest_cr_death_list
  
  if(length(comorbid_new_ckd_list) > 0) {
    comorbid_new_ckd_list_tmp <- vector(mode="list",length=length(comorbid_new_ckd_list))
    for(i in 1:length(comorbid_new_ckd_list)) {
      new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id",comorbid_new_ckd_list[i],"deceased")]
      new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(get(comorbid_new_ckd_list[i]),deceased)
      new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(deceased == 1)
      if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
        message(paste0(c("Including ",comorbid_new_ckd_list[i]," into the comorbid_new_ckd list...")))
        comorbid_new_ckd_list_tmp[i] <- comorbid_new_ckd_list[i]
      }
      rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
    }
    comorbid_new_ckd_list <- unlist(comorbid_new_ckd_list_tmp[lengths(comorbid_new_ckd_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(demog_new_ckd_list) > 0) {
    demog_new_ckd_list_tmp <- vector(mode="list",length=length(demog_new_ckd_list))
    for(i in 1:length(demog_new_ckd_list)) {
      new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id",demog_new_ckd_list[i],"deceased")]
      new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(get(demog_new_ckd_list[i]),deceased)
      new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(deceased == 1)
      if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
        message(paste0(c("Including ",demog_new_ckd_list[i]," into the demog_new_ckd list...")))
        demog_new_ckd_list_tmp[i] <- demog_new_ckd_list[i]
      }
      rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
    }
    demog_new_ckd_list <- unlist(demog_new_ckd_list_tmp[lengths(demog_new_ckd_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(med_new_ckd_list) > 0) {
    med_new_ckd_list_tmp <- vector(mode="list",length=length(med_new_ckd_list))
    if(length(med_new_ckd_list)>0) {
      for(i in 1:length(med_new_ckd_list)) {
        new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id",med_new_ckd_list[i],"deceased")]
        new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(get(med_new_ckd_list[i]),deceased)
        new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(deceased == 1)
        # if(min(new_ckd_tmp3$n) >= factor_cutoff) {
        #     comorbid_new_ckd_list_tmp[i] <- comorbid_new_ckd_list[i]
        # }
        if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
          message(paste0(c("Including ",med_new_ckd_list[i]," into the med_new_ckd list...")))
          med_new_ckd_list_tmp[i] <- med_new_ckd_list[i]
        }
        rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
      }
    }
    med_new_ckd_list <- unlist(med_new_ckd_list_tmp[lengths(med_new_ckd_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(!is.null(earliest_cr_new_ckd_list)) {
    new_ckd_tmp1 <- new_ckd_tmp[,c("patient_id","preadmit_cr_period","deceased")]
    new_ckd_tmp2 <- new_ckd_tmp1 %>% dplyr::count(preadmit_cr_period,deceased)
    new_ckd_tmp3 <- new_ckd_tmp2 %>% dplyr::filter(deceased == 1)
    if(min(new_ckd_tmp3$n) >= factor_cutoff & nrow(new_ckd_tmp3) > 1) {
      message(paste0(c("Including preadmit_cr_period into the earliest_cr_new_ckd list...")))
    } else {
      earliest_cr_new_ckd_list <- NULL # make this null or else you will get this causing errors
    }
    rm(new_ckd_tmp1,new_ckd_tmp2,new_ckd_tmp3)
    invisible(gc())
  }
  
  cat("\nFinal factor list for new_ckd (before user customisation): ",paste(c(demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list),collapse=" "))
  
  if(isTRUE(restrict_models)) {
    demog_new_ckd_list <- demog_new_ckd_list[demog_new_ckd_list %in% restrict_list]
    comorbid_new_ckd_list <- comorbid_new_ckd_list[comorbid_new_ckd_list %in% restrict_list]
    med_new_ckd_list <- med_new_ckd_list[med_new_ckd_list %in% restrict_list]
    if(earliest_cr_new_ckd_list %in% restrict_list) {
      earliest_cr_new_ckd_list <- "preadmit_cr_period"
    } else {
      earliest_cr_new_ckd_list <- NULL
    }
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_new_ckd_list,"\nComorbidities:",comorbid_new_ckd_list,"\nMedications:",med_new_ckd_list,"\nPreadmit Cr Period:",earliest_cr_new_ckd_list,sep = " "))
  }
  variable_list_output <- paste(c("Final Mortality In Non-CKD variable list:",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list),collapse=" ")
  var_list_death_nonckd <- c(demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list)
  
  # Run analysis for Mortality
  cat("\nNow proceeding with time-to-event analysis for Mortality (All Non-CKD patients)...")
  cat("\nGenerating Kaplan-Meier curves (Time to Death, all patients)...")
  try({
    cat("\n(a) by AKI occurrence")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ is_aki")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIvsNonAKI_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIvsNonAKI.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(b) by KDIGO AKI severity (no AKI = 0)")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ aki_kdigo_final")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_KDIGO.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity (ALL patients)
  
  try({
    cat("\n(c) by COVID-19 severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_Severe.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (Time to Death, Non-CKD patients)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_nonckd)}),
    error = function(c) {
      cat("\nError running univariate formulae (univ_models).")
      return(NULL)
    }
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_new_ckd_all <- do.call("rbind",univ_results)
    write.csv(univ_results_new_ckd_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for New CKD onset on Non-CKD patients\n==============\n")
  for(i in 1:14) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to Death, Non-CKD patients)..."))
    try({
      new_ckd_model <- c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list)
      new_ckd_model <- new_ckd_model[new_ckd_model %in% models_ckd[[i]]]
      newckdCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(new_ckd_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(new_ckd_model,collapse="+")))
      coxph_new_ckd <- survival::coxph(newckdCoxPHFormula, data=aki_index_nonckd)
      coxph_new_ckd_summ <- summary(coxph_new_ckd) 
      print(coxph_new_ckd_summ)
      coxph_new_ckd_hr <- cbind(coxph_new_ckd_summ$coefficients,coxph_new_ckd_summ$conf.int)[,-c(6,7)]
      coxph_new_ckd_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_new_ckd_summ$logtest,coxph_new_ckd_summ$sctest,coxph_new_ckd_summ$waldtest))
      coxph_new_ckd_stats2 <- rbind(data.table::as.data.table(coxph_new_ckd_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_new_ckd_summ$rsq,keep.rownames = T))
      write.csv(coxph_new_ckd_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_new_ckd_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_new_ckd_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_new_ckd_plot <- survminer::ggforest(coxph_new_ckd,data=aki_index_nonckd)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_new_ckd_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  # ===================
  
  # Now do New Onset CKD in Non-CKD but AKI patients only
  # Also add on mortality and AKI recovery outcomes (just in case)
  
  cat("\n=====================\nNow doing factor filtering for Non-CKD patients with AKI only (new CKD outcome)...\n")
  
  new_ckd_akionly_tmp <- aki_index_nonckd_akionly[,c("patient_id","new_ckd",demog_recovery_list,comorbid_recovery_list,med_recovery_list,earliest_cr_recovery_list)] %>% as.data.frame()
  comorbid_new_ckd_akionly_list <- comorbid_recovery_list[comorbid_recovery_list %in% colnames(new_ckd_akionly_tmp)]
  demog_new_ckd_akionly_list <- demog_recovery_list
  med_new_ckd_akionly_list <- med_recovery_list
  earliest_cr_new_ckd_akionly_list <- earliest_cr_recovery_list
  
  if(length(comorbid_new_ckd_akionly_list) > 0) {
    comorbid_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(comorbid_new_ckd_akionly_list))
    for(i in 1:length(comorbid_new_ckd_akionly_list)) {
      new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",comorbid_new_ckd_akionly_list[i],"new_ckd")]
      new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(comorbid_new_ckd_akionly_list[i]),new_ckd)
      new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(new_ckd == 1)
      if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
        message(paste0(c("Including ",comorbid_new_ckd_akionly_list[i]," into the comorbid_new_ckd_akionly list...")))
        comorbid_new_ckd_akionly_list_tmp[i] <- comorbid_new_ckd_akionly_list[i]
      }
      rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    }
    comorbid_new_ckd_akionly_list <- unlist(comorbid_new_ckd_akionly_list_tmp[lengths(comorbid_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(demog_new_ckd_akionly_list) > 0) {
    demog_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(demog_new_ckd_akionly_list))
    for(i in 1:length(demog_new_ckd_akionly_list)) {
      new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",demog_new_ckd_akionly_list[i],"new_ckd")]
      new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(demog_new_ckd_akionly_list[i]),new_ckd)
      new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(new_ckd == 1)
      if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
        message(paste0(c("Including ",demog_new_ckd_akionly_list[i]," into the demog_new_ckd_akionly list...")))
        demog_new_ckd_akionly_list_tmp[i] <- demog_new_ckd_akionly_list[i]
      }
      rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    }
    demog_new_ckd_akionly_list <- unlist(demog_new_ckd_akionly_list_tmp[lengths(demog_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(med_new_ckd_akionly_list) > 0) {
    med_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(med_new_ckd_akionly_list))
    if(length(med_new_ckd_akionly_list)>0) {
      for(i in 1:length(med_new_ckd_akionly_list)) {
        new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",med_new_ckd_akionly_list[i],"new_ckd")]
        new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(med_new_ckd_akionly_list[i]),new_ckd)
        new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(new_ckd == 1)
        # if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff) {
        #     comorbid_new_ckd_akionly_list_tmp[i] <- comorbid_new_ckd_akionly_list[i]
        # }
        if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
          message(paste0(c("Including ",med_new_ckd_akionly_list[i]," into the med_new_ckd_akionly list...")))
          med_new_ckd_akionly_list_tmp[i] <- med_new_ckd_akionly_list[i]
        }
        rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
      }
    }
    med_new_ckd_akionly_list <- unlist(med_new_ckd_akionly_list_tmp[lengths(med_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(!is.null(earliest_cr_new_ckd_akionly_list)) {
    new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id","preadmit_cr_period","new_ckd")]
    new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(preadmit_cr_period,new_ckd)
    new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(new_ckd == 1)
    if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
      message(paste0(c("Including preadmit_cr_period into the earliest_cr_new_ckd_akionly list...")))
    } else {
      earliest_cr_new_ckd_akionly_list <- NULL # make this null or else you will get this causing errors
    }
    rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    invisible(gc())
  }
  
  cat("\nFinal factor list for new_ckd_akionly (before user customisation): ",paste(c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list),collapse=" "))
  var_list_new_ckd_nonckd_akionly <- c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list)
  
  if(isTRUE(restrict_models)) {
    demog_new_ckd_akionly_list <- demog_new_ckd_akionly_list[demog_new_ckd_akionly_list %in% restrict_list]
    comorbid_new_ckd_akionly_list <- comorbid_new_ckd_akionly_list[comorbid_new_ckd_akionly_list %in% restrict_list]
    med_new_ckd_akionly_list <- med_new_ckd_akionly_list[med_new_ckd_akionly_list %in% restrict_list]
    if(earliest_cr_new_ckd_akionly_list %in% restrict_list) {
      earliest_cr_new_ckd_akionly_list <- "preadmit_cr_period"
    } else {
      earliest_cr_new_ckd_akionly_list <- NULL
    }
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_new_ckd_akionly_list,"\nComorbidities:",comorbid_new_ckd_akionly_list,"\nMedications:",med_new_ckd_akionly_list,"\nPreadmit Cr Period:",earliest_cr_new_ckd_akionly_list,sep = " "))
  }
  variable_list_output <- paste(c("Final New Onset CKD (AKI only) variable list:",demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list),collapse=" ")
  var_list_new_ckd_nonckd_akionly <- c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list)
  
  # Run analysis for new_ckd
  cat("\nNow proceeding with time-to-event analysis for New Onset CKD (All Non-CKD patients)...")
  cat("\nGenerating Kaplan-Meier curves (time to new onset CKD, all patients)...")
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(a) by KDIGO AKI severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ aki_kdigo_final")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd_akionly)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd_akionly)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_KDIGO.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity (ALL patients)
  
  try({
    cat("\n(b) by COVID-19 severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ severe")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd_akionly)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd_akionly)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_Severe.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (Time to new_ckd, all patients)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list), function(x) as.formula(paste('survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_nonckd_akionly)}),
    error = function(c) {
      cat("\nError running univariate formulae (univ_models).")
      return(NULL)
    }
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_new_ckd_all <- do.call("rbind",univ_results)
    write.csv(univ_results_new_ckd_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for New CKD onset on Non-CKD patients with AKI only\n==============\n")
  for(i in 1:14) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to New Onset CKD, Non-CKD patients with AKI only)..."))
    try({
      new_ckd_model <- c("severe","aki_kdigo_final",demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list)
      new_ckd_model <- new_ckd_model[new_ckd_model %in% models_ckd[[i]]]
      newckdCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ ",paste(new_ckd_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_new_ckd,event=new_ckd) ~ ",paste(new_ckd_model,collapse="+")))
      coxph_new_ckd <- survival::coxph(newckdCoxPHFormula, data=aki_index_nonckd_akionly)
      coxph_new_ckd_summ <- summary(coxph_new_ckd) 
      print(coxph_new_ckd_summ)
      coxph_new_ckd_hr <- cbind(coxph_new_ckd_summ$coefficients,coxph_new_ckd_summ$conf.int)[,-c(6,7)]
      coxph_new_ckd_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_new_ckd_summ$logtest,coxph_new_ckd_summ$sctest,coxph_new_ckd_summ$waldtest))
      coxph_new_ckd_stats2 <- rbind(data.table::as.data.table(coxph_new_ckd_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_new_ckd_summ$rsq,keep.rownames = T))
      write.csv(coxph_new_ckd_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_new_ckd_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_new_ckd_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_new_ckd_plot <- survminer::ggforest(coxph_new_ckd,data=aki_index_nonckd_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_new_ckd_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  # Run analysis for Mortality
  cat("\nNow proceeding with time-to-event analysis for Mortality (Non-CKD patients, AKI Only)...")
  
  cat("\nFirst, refiltering the variables for mortality analysis\n")
  
  new_ckd_akionly_tmp <- aki_index_nonckd_akionly[,c("patient_id","deceased",demog_recovery_list,comorbid_recovery_list,med_recovery_list,earliest_cr_recovery_list)] %>% as.data.frame()
  comorbid_new_ckd_akionly_list <- comorbid_recovery_list
  demog_new_ckd_akionly_list <- demog_recovery_list
  med_new_ckd_akionly_list <- med_recovery_list
  earliest_cr_new_ckd_akionly_list <- earliest_cr_recovery_list
  
  if(length(comorbid_new_ckd_akionly_list) > 0) {
    comorbid_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(comorbid_new_ckd_akionly_list))
    for(i in 1:length(comorbid_new_ckd_akionly_list)) {
      new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",comorbid_new_ckd_akionly_list[i],"deceased")]
      new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(comorbid_new_ckd_akionly_list[i]),deceased)
      new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(deceased == 1)
      if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
        message(paste0(c("Including ",comorbid_new_ckd_akionly_list[i]," into the comorbid_new_ckd_akionly list...")))
        comorbid_new_ckd_akionly_list_tmp[i] <- comorbid_new_ckd_akionly_list[i]
      }
      rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    }
    comorbid_new_ckd_akionly_list <- unlist(comorbid_new_ckd_akionly_list_tmp[lengths(comorbid_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(demog_new_ckd_akionly_list) > 0) {
    demog_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(demog_new_ckd_akionly_list))
    for(i in 1:length(demog_new_ckd_akionly_list)) {
      new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",demog_new_ckd_akionly_list[i],"deceased")]
      new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(demog_new_ckd_akionly_list[i]),deceased)
      new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(deceased == 1)
      if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
        message(paste0(c("Including ",demog_new_ckd_akionly_list[i]," into the demog_new_ckd_akionly list...")))
        demog_new_ckd_akionly_list_tmp[i] <- demog_new_ckd_akionly_list[i]
      }
      rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    }
    demog_new_ckd_akionly_list <- unlist(demog_new_ckd_akionly_list_tmp[lengths(demog_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(med_new_ckd_akionly_list) > 0) {
    med_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(med_new_ckd_akionly_list))
    if(length(med_new_ckd_akionly_list)>0) {
      for(i in 1:length(med_new_ckd_akionly_list)) {
        new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",med_new_ckd_akionly_list[i],"deceased")]
        new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(med_new_ckd_akionly_list[i]),deceased)
        new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(deceased == 1)
        # if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff) {
        #     comorbid_new_ckd_akionly_list_tmp[i] <- comorbid_new_ckd_akionly_list[i]
        # }
        if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
          message(paste0(c("Including ",med_new_ckd_akionly_list[i]," into the med_new_ckd_akionly list...")))
          med_new_ckd_akionly_list_tmp[i] <- med_new_ckd_akionly_list[i]
        }
        rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
      }
    }
    med_new_ckd_akionly_list <- unlist(med_new_ckd_akionly_list_tmp[lengths(med_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(!is.null(earliest_cr_new_ckd_akionly_list)) {
    new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id","preadmit_cr_period","deceased")]
    new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(preadmit_cr_period,deceased)
    new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(deceased == 1)
    if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
      message(paste0(c("Including preadmit_cr_period into the earliest_cr_new_ckd_akionly list...")))
    } else {
      earliest_cr_new_ckd_akionly_list <- NULL # make this null or else you will get this causing errors
    }
    rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    invisible(gc())
  }
  
  cat("\nFinal factor list for new_ckd_akionly (before user customisation): ",paste(c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list),collapse=" "))
  
  if(isTRUE(restrict_models)) {
    demog_new_ckd_akionly_list <- demog_new_ckd_akionly_list[demog_new_ckd_akionly_list %in% restrict_list]
    comorbid_new_ckd_akionly_list <- comorbid_new_ckd_akionly_list[comorbid_new_ckd_akionly_list %in% restrict_list]
    med_new_ckd_akionly_list <- med_new_ckd_akionly_list[med_new_ckd_akionly_list %in% restrict_list]
    if(earliest_cr_new_ckd_akionly_list %in% restrict_list) {
      earliest_cr_new_ckd_akionly_list <- "preadmit_cr_period"
    } else {
      earliest_cr_new_ckd_akionly_list <- NULL
    }
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_new_ckd_akionly_list,"\nComorbidities:",comorbid_new_ckd_akionly_list,"\nMedications:",med_new_ckd_akionly_list,"\nPreadmit Cr Period:",earliest_cr_new_ckd_akionly_list,sep = " "))
  }
  variable_list_output <- paste(c("Final Mortality (Non-CKD patients, AKI only) variable list:",demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list),collapse=" ")
  var_list_death_nonckd_akionly <- c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list)
  
  cat("\nGenerating Kaplan-Meier curves (Non-CKD patients, AKI Only)...")
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(a) by KDIGO AKI severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ aki_kdigo_final")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd_akionly)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd_akionly)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_KDIGO.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity
  try({
    cat("\n(b) by COVID-19 severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd_akionly)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd_akionly)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_Severe.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (Non-CKD patients, AKI Only)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_nonckd_akionly)}),
    error = function(c) {
      cat("\nError running univariate formulae (univ_models).")
      return(NULL)
    }
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_new_ckd_all <- do.call("rbind",univ_results)
    write.csv(univ_results_new_ckd_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for Mortality (Non-CKD patients, AKI Only)\n==============\n")
  for(i in 1:14) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to Death, Non-CKD patients, AKI Only)..."))
    try({
      new_ckd_model <- c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list)
      new_ckd_model <- new_ckd_model[new_ckd_model %in% models_ckd[[i]]]
      newckdCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(new_ckd_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(new_ckd_model,collapse="+")))
      coxph_new_ckd <- survival::coxph(newckdCoxPHFormula, data=aki_index_nonckd_akionly)
      coxph_new_ckd_summ <- summary(coxph_new_ckd) 
      print(coxph_new_ckd_summ)
      coxph_new_ckd_hr <- cbind(coxph_new_ckd_summ$coefficients,coxph_new_ckd_summ$conf.int)[,-c(6,7)]
      coxph_new_ckd_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_new_ckd_summ$logtest,coxph_new_ckd_summ$sctest,coxph_new_ckd_summ$waldtest))
      coxph_new_ckd_stats2 <- rbind(data.table::as.data.table(coxph_new_ckd_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_new_ckd_summ$rsq,keep.rownames = T))
      write.csv(coxph_new_ckd_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_new_ckd_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_new_ckd_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_new_ckd_plot <- survminer::ggforest(coxph_new_ckd,data=aki_index_nonckd_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_new_ckd_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  # --------------------
  # Run analysis for AKI Recovery
  cat("\nNow proceeding with time-to-event analysis for AKI Recovery (Non-CKD patients, AKI Only)...")
  
  cat("\nFirst, refiltering the variables for AKI Recovery analysis\n")
  
  new_ckd_akionly_tmp <- aki_index_nonckd_akionly[,c("patient_id","recover_1.25x",demog_recovery_list,comorbid_recovery_list,med_recovery_list,earliest_cr_recovery_list)] %>% as.data.frame()
  comorbid_new_ckd_akionly_list <- comorbid_recovery_list
  demog_new_ckd_akionly_list <- demog_recovery_list
  med_new_ckd_akionly_list <- med_recovery_list
  earliest_cr_new_ckd_akionly_list <- earliest_cr_recovery_list
  
  if(length(comorbid_new_ckd_akionly_list) > 0) {
    comorbid_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(comorbid_new_ckd_akionly_list))
    for(i in 1:length(comorbid_new_ckd_akionly_list)) {
      new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",comorbid_new_ckd_akionly_list[i],"recover_1.25x")]
      new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(comorbid_new_ckd_akionly_list[i]),recover_1.25x)
      new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(recover_1.25x == 1)
      if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
        message(paste0(c("Including ",comorbid_new_ckd_akionly_list[i]," into the comorbid_new_ckd_akionly list...")))
        comorbid_new_ckd_akionly_list_tmp[i] <- comorbid_new_ckd_akionly_list[i]
      }
      rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    }
    comorbid_new_ckd_akionly_list <- unlist(comorbid_new_ckd_akionly_list_tmp[lengths(comorbid_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(demog_new_ckd_akionly_list) > 0) {
    demog_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(demog_new_ckd_akionly_list))
    for(i in 1:length(demog_new_ckd_akionly_list)) {
      new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",demog_new_ckd_akionly_list[i],"recover_1.25x")]
      new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(demog_new_ckd_akionly_list[i]),recover_1.25x)
      new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(recover_1.25x == 1)
      if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
        message(paste0(c("Including ",demog_new_ckd_akionly_list[i]," into the demog_new_ckd_akionly list...")))
        demog_new_ckd_akionly_list_tmp[i] <- demog_new_ckd_akionly_list[i]
      }
      rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    }
    demog_new_ckd_akionly_list <- unlist(demog_new_ckd_akionly_list_tmp[lengths(demog_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(length(med_new_ckd_akionly_list) > 0) {
    med_new_ckd_akionly_list_tmp <- vector(mode="list",length=length(med_new_ckd_akionly_list))
    if(length(med_new_ckd_akionly_list)>0) {
      for(i in 1:length(med_new_ckd_akionly_list)) {
        new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id",med_new_ckd_akionly_list[i],"recover_1.25x")]
        new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(get(med_new_ckd_akionly_list[i]),recover_1.25x)
        new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(recover_1.25x == 1)
        # if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff) {
        #     comorbid_new_ckd_akionly_list_tmp[i] <- comorbid_new_ckd_akionly_list[i]
        # }
        if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
          message(paste0(c("Including ",med_new_ckd_akionly_list[i]," into the med_new_ckd_akionly list...")))
          med_new_ckd_akionly_list_tmp[i] <- med_new_ckd_akionly_list[i]
        }
        rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
      }
    }
    med_new_ckd_akionly_list <- unlist(med_new_ckd_akionly_list_tmp[lengths(med_new_ckd_akionly_list_tmp) > 0L])
    invisible(gc())
  }
  
  if(!is.null(earliest_cr_new_ckd_akionly_list)) {
    new_ckd_akionly_tmp1 <- new_ckd_akionly_tmp[,c("patient_id","preadmit_cr_period","recover_1.25x")]
    new_ckd_akionly_tmp2 <- new_ckd_akionly_tmp1 %>% dplyr::count(preadmit_cr_period,recover_1.25x)
    new_ckd_akionly_tmp3 <- new_ckd_akionly_tmp2 %>% dplyr::filter(recover_1.25x == 1)
    if(min(new_ckd_akionly_tmp3$n) >= factor_cutoff & nrow(new_ckd_akionly_tmp3) > 1) {
      message(paste0(c("Including preadmit_cr_period into the earliest_cr_new_ckd_akionly list...")))
    } else {
      earliest_cr_new_ckd_akionly_list <- NULL # make this null or else you will get this causing errors
    }
    rm(new_ckd_akionly_tmp1,new_ckd_akionly_tmp2,new_ckd_akionly_tmp3)
    invisible(gc())
  }
  
  cat("\nFinal factor list for new_ckd_akionly (before user customisation): ",paste(c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list),collapse=" "))
  
  if(isTRUE(restrict_models)) {
    demog_new_ckd_akionly_list <- demog_new_ckd_akionly_list[demog_new_ckd_akionly_list %in% restrict_list]
    comorbid_new_ckd_akionly_list <- comorbid_new_ckd_akionly_list[comorbid_new_ckd_akionly_list %in% restrict_list]
    med_new_ckd_akionly_list <- med_new_ckd_akionly_list[med_new_ckd_akionly_list %in% restrict_list]
    if(earliest_cr_new_ckd_akionly_list %in% restrict_list) {
      earliest_cr_new_ckd_akionly_list <- "preadmit_cr_period"
    } else {
      earliest_cr_new_ckd_akionly_list <- NULL
    }
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_new_ckd_akionly_list,"\nComorbidities:",comorbid_new_ckd_akionly_list,"\nMedications:",med_new_ckd_akionly_list,"\nPreadmit Cr Period:",earliest_cr_new_ckd_akionly_list,sep = " "))
  }
  variable_list_output <- paste(c("Final AKI Recovery (Non-CKD patients, AKI only) variable list:",demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list),collapse=" ")
  var_list_recovery_nonckd_akionly <- c(demog_new_ckd_akionly_list,comorbid_new_ckd_akionly_list,med_new_ckd_akionly_list,earliest_cr_new_ckd_akionly_list)
  
  cat("\nGenerating Kaplan-Meier curves (AKI Recovery, Non-CKD patients, AKI Only)...")
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(a) by KDIGO AKI severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ aki_kdigo_final")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd_akionly)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd_akionly)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_KDIGO.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity
  try({
    cat("\n(b) by COVID-19 severity")
    newckdPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe")
    fit_newckd <- survminer::surv_fit(newckdPlotFormula, data=aki_index_nonckd_akionly)
    plot_newckd <- survminer::ggsurvplot(fit_newckd,data=aki_index_nonckd_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_newckd_summ <- survminer::surv_summary(fit_newckd,data=aki_index_nonckd_akionly)
    plot_newckd_summ_table <- plot_newckd$data.survtable
    write.csv(plot_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_Severe.png")),plot=print(plot_newckd,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (AKI Recovery, Non-CKD patients, AKI Only)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list), function(x) as.formula(paste('survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_nonckd_akionly)}),
    error = function(c) {
      cat("\nError running univariate formulae (univ_models).")
      return(NULL)
    }
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_new_ckd_all <- do.call("rbind",univ_results)
    write.csv(univ_results_new_ckd_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for AKI Recovery onset (Non-CKD patients, AKI Only)\n==============\n")
  for(i in 1:14) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to AKI Recovery, Non-CKD patients, AKI Only)..."))
    try({
      new_ckd_model <- c("severe","aki_kdigo_final",demog_new_ckd_list,comorbid_new_ckd_list,med_new_ckd_list,earliest_cr_new_ckd_list)
      new_ckd_model <- new_ckd_model[new_ckd_model %in% models_ckd[[i]]]
      newckdCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(new_ckd_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(new_ckd_model,collapse="+")))
      coxph_new_ckd <- survival::coxph(newckdCoxPHFormula, data=aki_index_nonckd_akionly)
      coxph_new_ckd_summ <- summary(coxph_new_ckd) 
      print(coxph_new_ckd_summ)
      coxph_new_ckd_hr <- cbind(coxph_new_ckd_summ$coefficients,coxph_new_ckd_summ$conf.int)[,-c(6,7)]
      coxph_new_ckd_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_new_ckd_summ$logtest,coxph_new_ckd_summ$sctest,coxph_new_ckd_summ$waldtest))
      coxph_new_ckd_stats2 <- rbind(data.table::as.data.table(coxph_new_ckd_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_new_ckd_summ$rsq,keep.rownames = T))
      write.csv(coxph_new_ckd_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_new_ckd_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_new_ckd_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_new_ckd_plot <- survminer::ggforest(coxph_new_ckd,data=aki_index_nonckd_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_new_ckd_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  return(list("aki_index_nonckd" = aki_index_nonckd, "aki_index_nonckd_akionly" = aki_index_nonckd_akionly,"var_list_new_ckd_nonckd_akionly" = var_list_new_ckd_nonckd_akionly, "var_list_recovery_nonckd_akionly" = var_list_recovery_nonckd_akionly, "var_list_death_nonckd_akionly" = var_list_death_nonckd_akionly, "var_list_death_nonckd" = var_list_death_nonckd, "var_list_new_ckd_nonckd" = var_list_new_ckd_nonckd))
}
