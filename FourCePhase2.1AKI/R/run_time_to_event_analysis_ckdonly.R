#' Runs time-to-event analysis for the subgroup of CKD patients
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
run_time_to_event_analysis_ckdonly <- function(siteid,
                                              aki_index_recovery,aki_index_death,
                                              var_list_recovery,var_list_death,
                                              obfuscation,obfuscation_level,
                                              restrict_model_corr, factor_threshold = 5,
                                              use_custom_output = FALSE,use_custom_output_dir = "/4ceData/Output") {
  currSiteId <- siteid
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  restrict_models <- restrict_model_corr
  factor_cutoff <- factor_threshold
  var_list_death_ckdonly <- var_list_death
  var_list_recovery_all <- var_list_recovery
  
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
  
  # CKD/CKD subgroup
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
  
  # D: Include serum creatinine cutoffs (CKD diagnosis coding may be defined using urine criteria as well)
  model_ckd_subgroup_1d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_ckd_subgroup_3d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_ckd <- list(model_ckd_subgroup_1,model_ckd_subgroup_2,model_ckd_subgroup_3,model_ckd_subgroup_4,model_ckd_subgroup_1a,model_ckd_subgroup_2a,model_ckd_subgroup_3a,model_ckd_subgroup_4a,model_ckd_subgroup_3b,model_ckd_subgroup_4b,model_ckd_subgroup_1c,model_ckd_subgroup_2c,model_ckd_subgroup_3c,model_ckd_subgroup_4c,model_ckd_subgroup_1d,model_ckd_subgroup_2d,model_ckd_subgroup_3d,model_ckd_subgroup_4d)
  models_ckd_labels <- c("Model_CKD_Subgroup_1","Model_CKD_Subgroup_2","Model_CKD_Subgroup_3","Model_CKD_Subgroup_4","Model_CKD_Subgroup_1A","Model_CKD_Subgroup_2A","Model_CKD_Subgroup_3A","Model_CKD_Subgroup_4A","Model_CKD_Subgroup_3B","Model_CKD_Subgroup_4B","Model_CKD_Subgroup_1C","Model_CKD_Subgroup_2C","Model_CKD_Subgroup_3C","Model_CKD_Subgroup_4C","Model_CKD_Subgroup_1D","Model_CKD_Subgroup_2D","Model_CKD_Subgroup_3D","Model_CKD_Subgroup_4D")
  
  # Generate tables intended to assess new onset of new CKD
  aki_index_ckdonly <- aki_index_death[aki_index_death$ckd == 1,]
  aki_index_ckdonly_akionly <- aki_index_recovery[aki_index_recovery$ckd == 1,]
  
  # Same workflow
  # 1) Filter for variables with enough levels
  # 2) Filter variables with enough counts on both sides of levels
  # 3) Run models on the filtered list of covariates
  
  ckdonly_tmp <- aki_index_ckdonly[,c("patient_id","deceased",var_list_death_ckdonly)] %>% as.data.frame()
  
  if(length(var_list_death_ckdonly) > 0) {
    var_list_death_ckdonly_tmp <- vector(mode="list",length=length(var_list_death_ckdonly))
    for(i in 1:length(var_list_death_ckdonly)) {
      ckdonly_tmp1 <- ckdonly_tmp[,c("patient_id",var_list_death_ckdonly[i],"deceased")]
      ckdonly_tmp2 <- ckdonly_tmp1 %>% dplyr::count(get(var_list_death_ckdonly[i]),deceased)
      ckdonly_tmp3 <- ckdonly_tmp2 %>% dplyr::filter(deceased == 1)
      if(min(ckdonly_tmp3$n) >= factor_cutoff & nrow(ckdonly_tmp3) > 1) {
        message(paste0(c("Including ",var_list_death_ckdonly[i]," into the comorbid_ckdonly list...")))
        var_list_death_ckdonly_tmp[i] <- var_list_death_ckdonly[i]
      }
      rm(ckdonly_tmp1,ckdonly_tmp2,ckdonly_tmp3)
    }
    var_list_death_ckdonly <- unlist(var_list_death_ckdonly_tmp[lengths(var_list_death_ckdonly_tmp) > 0L])
    invisible(gc())
  }
  cat("\nFinal factor list for ckdonly (before user customisation): ",paste(var_list_death_ckdonly,collapse=" "))
  
  if(isTRUE(restrict_models)) {
    var_list_death_ckdonly <- var_list_death_ckdonly[var_list_death_ckdonly %in% restrict_list]
    message(paste(c("\nAfter filtering for custom-specified variables, we have the following:",var_list_death_ckdonly),collapse = " "))
  }
  variable_list_output <- paste(c("Final Mortality In CKD variable list:",var_list_death_ckdonly),collapse=" ")

  # Run analysis for Mortality
  cat("\nNow proceeding with time-to-event analysis for Mortality (All CKD patients)...")
  cat("\nGenerating Kaplan-Meier curves (Time to Death, all patients)...")
  try({
    cat("\n(a) by AKI occurrence")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ is_aki")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIvsNonAKI_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIvsNonAKI.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(b) by KDIGO AKI severity (no AKI = 0)")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ aki_kdigo_final")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_KDIGO.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity (ALL patients)
  
  try({
    cat("\n(c) by COVID-19 severity")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_Severe.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (Time to Death, CKD patients)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",var_list_death_ckdonly), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_ckdonly)}),
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
    univ_results_ckdonly_all <- do.call("rbind",univ_results)
    write.csv(univ_results_ckdonly_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for Mortality on CKD patients\n==============\n")
  for(i in 1:14) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to Death, CKD patients)..."))
    try({
      ckdonly_model <- c("severe","aki_kdigo_final",var_list_death_ckdonly)
      ckdonly_model <- ckdonly_model[ckdonly_model %in% models_ckd[[i]]]
      ckdonlyCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
      coxph_ckdonly <- survival::coxph(ckdonlyCoxPHFormula, data=aki_index_ckdonly)
      coxph_ckdonly_summ <- summary(coxph_ckdonly) 
      print(coxph_ckdonly_summ)
      coxph_ckdonly_hr <- cbind(coxph_ckdonly_summ$coefficients,coxph_ckdonly_summ$conf.int)[,-c(6,7)]
      coxph_ckdonly_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_ckdonly_summ$logtest,coxph_ckdonly_summ$sctest,coxph_ckdonly_summ$waldtest))
      coxph_ckdonly_stats2 <- rbind(data.table::as.data.table(coxph_ckdonly_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_ckdonly_summ$rsq,keep.rownames = T))
      write.csv(coxph_ckdonly_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_ckdonly_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_ckdonly_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  # ===================
  
  # Run analysis for Mortality, CKD with AKI only
  cat("\nNow proceeding with time-to-event analysis for Mortality (CKD patients, AKI Only)...")
  
  cat("\nFirst, refiltering the variables for mortality analysis\n")
  # Reset the var list
  var_list_death_ckdonly_akionly <- var_list_recovery_all
  ckdonly_akionly_tmp <- aki_index_ckdonly_akionly[,c("patient_id","deceased",var_list_death_ckdonly_akionly)] %>% as.data.frame()
  
  if(length(var_list_death_ckdonly_akionly) > 0) {
    var_list_death_ckdonly_akionly_tmp <- vector(mode="list",length=length(var_list_death_ckdonly_akionly))
    for(i in 1:length(var_list_death_ckdonly_akionly)) {
      ckdonly_akionly_tmp1 <- ckdonly_akionly_tmp[,c("patient_id",var_list_death_ckdonly_akionly[i],"deceased")]
      ckdonly_akionly_tmp2 <- ckdonly_akionly_tmp1 %>% dplyr::count(get(var_list_death_ckdonly_akionly[i]),deceased)
      ckdonly_akionly_tmp3 <- ckdonly_akionly_tmp2 %>% dplyr::filter(deceased == 1)
      if(min(ckdonly_akionly_tmp3$n) >= factor_cutoff & nrow(ckdonly_akionly_tmp3) > 1) {
        message(paste0(c("Including ",var_list_death_ckdonly_akionly[i]," into the comorbid_ckdonly_akionly list...")))
        var_list_death_ckdonly_akionly[i] <- var_list_death_ckdonly_akionly[i]
      }
      rm(ckdonly_akionly_tmp1,ckdonly_akionly_tmp2,ckdonly_akionly_tmp3)
    }
    var_list_death_ckdonly_akionly <- unlist(var_list_death_ckdonly_akionly_tmp[lengths(var_list_death_ckdonly_akionly_tmp) > 0L])
    invisible(gc())
  }
  
  cat("\nFinal factor list for ckdonly_akionly (before user customisation): ",paste(var_list_death_ckdonly_akionly,collapse=" "),"\n")
  
  if(isTRUE(restrict_models)) {
    var_list_death_ckdonly_akionly <- var_list_death_ckdonly_akionly[var_list_death_ckdonly_akionly %in% restrict_list]
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\n",var_list_death_ckdonly_akionly,sep = " "))
  }
  variable_list_output <- paste(c("Final Mortality (CKD patients, AKI only) variable list:",var_list_death_ckdonly_akionly),collapse=" ")
  
  cat("\nGenerating Kaplan-Meier curves (CKD patients, AKI Only)...")
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(a) by KDIGO AKI severity")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ aki_kdigo_final")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly_akionly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly_akionly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_KDIGO.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity
  try({
    cat("\n(b) by COVID-19 severity")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly_akionly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly_akionly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_Severe.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (CKD patients, AKI Only)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",var_list_death_ckdonly_akionly), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_ckdonly_akionly)}),
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
    univ_results_ckdonly_all <- do.call("rbind",univ_results)
    write.csv(univ_results_ckdonly_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for Mortality (CKD patients, AKI Only)\n==============\n")
  for(i in 1:length(models_ckd_labels)) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to Death, CKD patients, AKI Only)..."))
    try({
      ckdonly_model <- c("severe","aki_kdigo_final",var_list_death_ckdonly_akionly)
      ckdonly_model <- ckdonly_model[ckdonly_model %in% models_ckd[[i]]]
      ckdonlyCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
      coxph_ckdonly <- survival::coxph(ckdonlyCoxPHFormula, data=aki_index_ckdonly_akionly)
      coxph_ckdonly_summ <- summary(coxph_ckdonly) 
      print(coxph_ckdonly_summ)
      coxph_ckdonly_hr <- cbind(coxph_ckdonly_summ$coefficients,coxph_ckdonly_summ$conf.int)[,-c(6,7)]
      coxph_ckdonly_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_ckdonly_summ$logtest,coxph_ckdonly_summ$sctest,coxph_ckdonly_summ$waldtest))
      coxph_ckdonly_stats2 <- rbind(data.table::as.data.table(coxph_ckdonly_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_ckdonly_summ$rsq,keep.rownames = T))
      write.csv(coxph_ckdonly_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_ckdonly_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_ckdonly_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  # --------------------
  # Run analysis for AKI Recovery
  cat("\nNow proceeding with time-to-event analysis for AKI Recovery (CKD patients, AKI Only)...")
  
  cat("\nFirst, refiltering the variables for AKI Recovery analysis\n")
  var_list_recovery_ckdonly_akionly <- var_list_recovery_all
  ckdonly_akionly_tmp <- aki_index_ckdonly_akionly[,c("patient_id","recover_1.25x",var_list_recovery_ckdonly_akionly)] %>% as.data.frame()
  
  if(length(var_list_recovery_ckdonly_akionly) > 0) {
    var_list_recovery_ckdonly_akionly_tmp <- vector(mode="list",length=length(var_list_recovery_ckdonly_akionly))
    for(i in 1:length(var_list_recovery_ckdonly_akionly)) {
      ckdonly_akionly_tmp1 <- ckdonly_akionly_tmp[,c("patient_id",var_list_recovery_ckdonly_akionly[i],"recover_1.25x")]
      ckdonly_akionly_tmp2 <- ckdonly_akionly_tmp1 %>% dplyr::count(get(var_list_recovery_ckdonly_akionly[i]),recover_1.25x)
      ckdonly_akionly_tmp3 <- ckdonly_akionly_tmp2 %>% dplyr::filter(recover_1.25x == 1)
      if(min(ckdonly_akionly_tmp3$n) >= factor_cutoff & nrow(ckdonly_akionly_tmp3) > 1) {
        message(paste0(c("Including ",var_list_recovery_ckdonly_akionly[i]," into the comorbid_ckdonly_akionly list...")))
        var_list_recovery_ckdonly_akionly_tmp[i] <- var_list_recovery_ckdonly_akionly[i]
      }
      rm(ckdonly_akionly_tmp1,ckdonly_akionly_tmp2,ckdonly_akionly_tmp3)
    }
    var_list_recovery_ckdonly_akionly <- unlist(var_list_recovery_ckdonly_akionly_tmp[lengths(var_list_recovery_ckdonly_akionly_tmp) > 0L])
    invisible(gc())
  }
  
  cat("\nFinal factor list for ckdonly_akionly (before user customisation): ",paste(var_list_recovery_ckdonly_akionly,collapse=" "))
  
  
  if(isTRUE(restrict_models)) {
    var_list_recovery_ckdonly_akionly <- var_list_recovery_ckdonly_akionly[var_list_recovery_ckdonly_akionly %in% restrict_list]
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\n",var_list_recovery_ckdonly_akionly,sep = " "))
  }
  variable_list_output <- paste(c("Final AKI Recovery (CKD patients, AKI only) variable list:",var_list_recovery_ckdonly_akionly),collapse=" ")
  
  cat("\nGenerating Kaplan-Meier curves (AKI Recovery, CKD patients, AKI Only)...")
  # Survival curves stratified by KDIGO stage
  try( {
    cat("\n(a) by KDIGO AKI severity")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ aki_kdigo_final")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly_akionly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly_akionly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_KDIGO_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_KDIGO.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  # Survival curves stratified by COVID-19 severity
  try({
    cat("\n(b) by COVID-19 severity")
    ckdonlyPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe")
    fit_ckdonly <- survminer::surv_fit(ckdonlyPlotFormula, data=aki_index_ckdonly_akionly)
    plot_ckdonly <- survminer::ggsurvplot(fit_ckdonly,data=aki_index_ckdonly_akionly,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_ckdonly_summ <- survminer::surv_summary(fit_ckdonly,data=aki_index_ckdonly_akionly)
    plot_ckdonly_summ_table <- plot_ckdonly$data.survtable
    write.csv(plot_ckdonly_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    plot.new()
    ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_Severe.png")),plot=print(plot_ckdonly,newpage=F),width=12,height=12,units="cm")
  })
  
  cat("\nNow proceeding with Cox PH model time-to-event analysis...")
  cat("\nGenerating univariate Cox PH models (AKI Recovery, CKD patients, AKI Only)...")
  univ_formulas <- tryCatch({ 
    sapply(c("severe","aki_kdigo_final",var_list_recovery_ckdonly_akionly), function(x) as.formula(paste('survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ', x)))
  }, error = function(c) {
    cat("\nError running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_ckdonly_akionly)}),
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
    univ_results_ckdonly_all <- do.call("rbind",univ_results)
    write.csv(univ_results_ckdonly_all,file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  cat("\n=================\nRunning Models (minus CKD), Both Original and Supp, for AKI Recovery onset (CKD patients, AKI Only)\n==============\n")
  for(i in 1:length(models_ckd_labels)) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to AKI Recovery, CKD patients, AKI Only)..."))
    try({
      ckdonly_model <- c("severe","aki_kdigo_final",var_list_recovery_ckdonly_akionly)
      ckdonly_model <- ckdonly_model[ckdonly_model %in% models_ckd[[i]]]
      ckdonlyCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(ckdonly_model,collapse="+")))
      message(paste("Formula for ", models_ckd_labels[i],": survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(ckdonly_model,collapse="+")))
      coxph_ckdonly <- survival::coxph(ckdonlyCoxPHFormula, data=aki_index_ckdonly_akionly)
      coxph_ckdonly_summ <- summary(coxph_ckdonly) 
      print(coxph_ckdonly_summ)
      coxph_ckdonly_hr <- cbind(coxph_ckdonly_summ$coefficients,coxph_ckdonly_summ$conf.int)[,-c(6,7)]
      coxph_ckdonly_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_ckdonly_summ$logtest,coxph_ckdonly_summ$sctest,coxph_ckdonly_summ$waldtest))
      coxph_ckdonly_stats2 <- rbind(data.table::as.data.table(coxph_ckdonly_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_ckdonly_summ$rsq,keep.rownames = T))
      write.csv(coxph_ckdonly_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],".csv")),row.names=TRUE)
      write.csv(coxph_ckdonly_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
      write.csv(coxph_ckdonly_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
      coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  return(list("aki_index_ckdonly" = aki_index_ckdonly, "aki_index_ckdonly_akionly" = aki_index_ckdonly_akionly,"var_list_recovery_ckdonly_akionly" = var_list_recovery_ckdonly_akionly, "var_list_death_ckdonly_akionly" = var_list_death_ckdonly_akionly, "var_list_death_ckdonly" = var_list_death_ckdonly))
}
