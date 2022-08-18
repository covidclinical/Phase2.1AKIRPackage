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
                                              use_custom_output = FALSE,use_custom_output_dir = "/4ceData/Output",input) {
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
    restrict_list <- scan(file.path(input,"CustomModelVariables.txt"),what="")
    message(paste("Variables to restrict analyses to :",restrict_list,collapse=" "))
  }
  
  # CKD/CKD subgroup
  model_ckd_subgroup_1 <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2 <- c("age_group","sex","severe","copd","rheum","vte")
  model_ckd_subgroup_3 <- c("age_group","sex","severe","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4 <- c("age_group","sex","severe","aki_kdigo_final","acei_arb_preexposure")
  model_ckd_subgroup_1a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte")
  model_ckd_subgroup_3a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_3b <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral")
  model_ckd_subgroup_4b <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_1c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte")
  model_ckd_subgroup_3c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # D: Include serum creatinine cutoffs (CKD diagnosis coding may be defined using urine criteria as well)
  model_ckd_subgroup_1d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld")
  model_ckd_subgroup_2d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte")
  model_ckd_subgroup_3d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_ckd_subgroup_4d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # E: Use Models 4A,B,C,D but all the age groups as strata
  model_ckd_subgroup_1e <- c("age_group_original","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2e <-  c("age_group_original","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3e <- c("age_group_original","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4e <- c("age_group_original","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # F: Use Models 4A,B,C,D but split age to three strata (< 50, 50-69, >=70)
  model_ckd_subgroup_1f <- c("age_group_2","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2f <-  c("age_group_2","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3f <- c("age_group_2","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4f <- c("age_group_2","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # G: Use Models 4A,B,C,D but split age by strata (<50,50-69,70-79,>=80)
  model_ckd_subgroup_1g <- c("age_group_3","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2g <-  c("age_group_3","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3g <- c("age_group_3","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4g <- c("age_group_3","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # H: Use Models 4A,B,C,D but split age by >=50 instead of >= 70
  model_ckd_subgroup_1h <- c("age_group_4","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2h <-  c("age_group_4","sex","severe","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3h <- c("age_group_4","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4h <- c("age_group_4","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_ckd <- list(model_ckd_subgroup_1,model_ckd_subgroup_2,model_ckd_subgroup_3,model_ckd_subgroup_4,model_ckd_subgroup_1a,model_ckd_subgroup_2a,model_ckd_subgroup_3a,model_ckd_subgroup_4a,model_ckd_subgroup_3b,model_ckd_subgroup_4b,model_ckd_subgroup_1c,model_ckd_subgroup_2c,model_ckd_subgroup_3c,model_ckd_subgroup_4c,model_ckd_subgroup_1d,model_ckd_subgroup_2d,model_ckd_subgroup_3d,model_ckd_subgroup_4d)
  models_ckd_labels <- c("Model_CKD_Subgroup_1","Model_CKD_Subgroup_2","Model_CKD_Subgroup_3","Model_CKD_Subgroup_4","Model_CKD_Subgroup_1A","Model_CKD_Subgroup_2A","Model_CKD_Subgroup_3A","Model_CKD_Subgroup_4A","Model_CKD_Subgroup_3B","Model_CKD_Subgroup_4B","Model_CKD_Subgroup_1C","Model_CKD_Subgroup_2C","Model_CKD_Subgroup_3C","Model_CKD_Subgroup_4C","Model_CKD_Subgroup_1D","Model_CKD_Subgroup_2D","Model_CKD_Subgroup_3D","Model_CKD_Subgroup_4D")
  
  models_ckd_age <- list(model_ckd_subgroup_1e,model_ckd_subgroup_2e,model_ckd_subgroup_3e,model_ckd_subgroup_4e,
                         model_ckd_subgroup_1f,model_ckd_subgroup_2f,model_ckd_subgroup_3f,model_ckd_subgroup_4f,
                         model_ckd_subgroup_1g,model_ckd_subgroup_2g,model_ckd_subgroup_3g,model_ckd_subgroup_4g,
                         model_ckd_subgroup_1h,model_ckd_subgroup_2h,model_ckd_subgroup_3h,model_ckd_subgroup_4h)
  
  models_ckd_age_labels <- c("Model_CKD_Subgroup_1E","Model_CKD_Subgroup_2E","Model_CKD_Subgroup_3E","Model_CKD_Subgroup_4E",
                             "Model_CKD_Subgroup_1F","Model_CKD_Subgroup_2F","Model_CKD_Subgroup_3F","Model_CKD_Subgroup_4F",
                             "Model_CKD_Subgroup_1G","Model_CKD_Subgroup_2G","Model_CKD_Subgroup_3G","Model_CKD_Subgroup_4G",
                             "Model_CKD_Subgroup_1H","Model_CKD_Subgroup_2H","Model_CKD_Subgroup_3H","Model_CKD_Subgroup_4H")
  
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
  
  # Now fix the factoring for age groups
  if("age_group_original" %in% var_list_death_ckdonly) {
    aki_index_ckdonly$age_group_original <- factor(aki_index_ckdonly$age_group_original,levels=c("18to25", "26to49", "50to69", "70to79", "80plus"))
  }
  if("age_group_2" %in% var_list_death_ckdonly) {
    aki_index_ckdonly$age_group_2 <- factor(aki_index_ckdonly$age_group_2,levels=c("below_50", "50to69", "70_and_above"))
  }
  if("age_group_3" %in% var_list_death_ckdonly) {
    aki_index_ckdonly$age_group_3 <- factor(aki_index_ckdonly$age_group_3,levels=c("below_50", "50to69", "70to79", "80plus"))
  }
  if("age_group_4" %in% var_list_death_ckdonly) {
    aki_index_ckdonly$age_group_4 <- factor(aki_index_ckdonly$age_group_4,levels=c("below_50", "50_and_above"))
  }
  

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
  for(i in 1:length(models_ckd_labels)) {
    cat(paste0("\nGenerating", models_ckd_labels[i], "(Time to Death, CKD patients)..."))
    try({
      ckdonly_model <- c("severe","aki_kdigo_final",var_list_death_ckdonly)
      if(isTRUE(restrict_models)) {
        ckdonly_model <- ckdonly_model[ckdonly_model %in% restrict_list]
      }
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
      coxph_ckdonly_vcov <- vcov(coxph_ckdonly)
      write.csv(coxph_ckdonly_vcov,file=file.path(dir.output,paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_labels[i],"_vcov.csv")))
      coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  cat("\n=================\nRunning Models (minus CKD), Supp (Age), for Mortality on CKD patients\n==============\n")
  for(i in 1:length(models_ckd_age_labels)) {
    cat(paste0("\nGenerating ", models_ckd_age_labels[i], " (Time to Death, CKD Patients Only)..."))
    if((i %in% c(1:4) & "age_group_original" %in% var_list_death_ckdonly) |
       (i %in% c(5:8) & "age_group_2" %in% var_list_death_ckdonly) |
       (i %in% c(9:12) & "age_group_3" %in% var_list_death_ckdonly) |
       (i %in% c(13:16) & "age_group_4" %in% var_list_death_ckdonly)) {
      try({
        ckdonly_model <- c("severe","aki_kdigo_final",var_list_death_ckdonly)
        if(isTRUE(restrict_models)) {
          ckdonly_model <- ckdonly_model[ckdonly_model %in% restrict_list]
        }
        ckdonly_model <- ckdonly_model[ckdonly_model %in% models_ckd_age[[i]]]
        ckdonlyCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
        message(paste("Formula for ", models_ckd_age_labels[i],": survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
        coxph_ckdonly <- survival::coxph(ckdonlyCoxPHFormula, data=aki_index_ckdonly)
        coxph_ckdonly_summ <- summary(coxph_ckdonly) 
        print(coxph_ckdonly_summ)
        coxph_ckdonly_hr <- cbind(coxph_ckdonly_summ$coefficients,coxph_ckdonly_summ$conf.int)[,-c(6,7)]
        coxph_ckdonly_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_ckdonly_summ$logtest,coxph_ckdonly_summ$sctest,coxph_ckdonly_summ$waldtest))
        coxph_ckdonly_stats2 <- rbind(data.table::as.data.table(coxph_ckdonly_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_ckdonly_summ$rsq,keep.rownames = T))
        write.csv(coxph_ckdonly_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_age_labels[i],".csv")),row.names=TRUE)
        write.csv(coxph_ckdonly_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_age_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
        write.csv(coxph_ckdonly_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_age_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
        coxph_ckdonly_vcov <- vcov(coxph_ckdonly)
        write.csv(coxph_ckdonly_vcov,file=file.path(dir.output,paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_age_labels[i],"_vcov.csv")))
        coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly)
        ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_CoxPH_",models_ckd_age_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_ckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
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
  
  # Now fix the factoring for age groups
  if("age_group_original" %in% var_list_death_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_original <- factor(aki_index_ckdonly_akionly$age_group_original,levels=c("18to25", "26to49", "50to69", "70to79", "80plus"))
  }
  if("age_group_2" %in% var_list_death_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_2 <- factor(aki_index_ckdonly_akionly$age_group_2,levels=c("below_50", "50to69", "70_and_above"))
  }
  if("age_group_3" %in% var_list_death_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_3 <- factor(aki_index_ckdonly_akionly$age_group_3,levels=c("below_50", "50to69", "70to79", "80plus"))
  }
  if("age_group_4" %in% var_list_death_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_4 <- factor(aki_index_ckdonly_akionly$age_group_4,levels=c("below_50", "50_and_above"))
  }
  
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
      if(isTRUE(restrict_models)) {
        ckdonly_model <- ckdonly_model[ckdonly_model %in% restrict_list]
      }
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
      coxph_ckdonly_vcov <- vcov(coxph_ckdonly)
      write.csv(coxph_ckdonly_vcov,file=file.path(dir.output,paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],"_vcov.csv")))
      coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  cat("\n=================\nRunning Models (minus CKD), Supp (Age), for Mortality on CKD and AKI patients\n==============\n")
  for(i in 1:length(models_ckd_age_labels)) {
    cat(paste0("\nGenerating ", models_ckd_age_labels[i], " (Time to Mortality, CKD, AKI Patients Only)..."))
    if((i %in% c(1:4) & "age_group_original" %in% var_list_death_ckdonly_akionly) |
       (i %in% c(5:8) & "age_group_2" %in% var_list_death_ckdonly_akionly) |
       (i %in% c(9:12) & "age_group_3" %in% var_list_death_ckdonly_akionly) |
       (i %in% c(13:16) & "age_group_4" %in% var_list_death_ckdonly_akionly)) {
      try({
        ckdonly_model <- c("severe","aki_kdigo_final",var_list_death_ckdonly_akionly)
        if(isTRUE(restrict_models)) {
          ckdonly_model <- ckdonly_model[ckdonly_model %in% restrict_list]
        }
        ckdonly_model <- ckdonly_model[ckdonly_model %in% models_ckd_age[[i]]]
        ckdonlyCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
        message(paste("Formula for ", models_ckd_age_labels[i],": survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(ckdonly_model,collapse="+")))
        coxph_ckdonly <- survival::coxph(ckdonlyCoxPHFormula, data=aki_index_ckdonly_akionly)
        coxph_ckdonly_summ <- summary(coxph_ckdonly) 
        print(coxph_ckdonly_summ)
        coxph_ckdonly_hr <- cbind(coxph_ckdonly_summ$coefficients,coxph_ckdonly_summ$conf.int)[,-c(6,7)]
        coxph_ckdonly_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_ckdonly_summ$logtest,coxph_ckdonly_summ$sctest,coxph_ckdonly_summ$waldtest))
        coxph_ckdonly_stats2 <- rbind(data.table::as.data.table(coxph_ckdonly_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_ckdonly_summ$rsq,keep.rownames = T))
        write.csv(coxph_ckdonly_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],".csv")),row.names=TRUE)
        write.csv(coxph_ckdonly_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
        write.csv(coxph_ckdonly_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
        coxph_ckdonly_vcov <- vcov(coxph_ckdonly)
        write.csv(coxph_ckdonly_vcov,file=file.path(dir.output,paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],"_vcov.csv")))
        coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly_akionly)
        ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Death_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_ckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
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
  
  # Now fix the factoring for age groups
  if("age_group_original" %in% var_list_recovery_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_original <- factor(aki_index_ckdonly_akionly$age_group_original,levels=c("18to25", "26to49", "50to69", "70to79", "80plus"))
  }
  if("age_group_2" %in% var_list_recovery_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_2 <- factor(aki_index_ckdonly_akionly$age_group_2,levels=c("below_50", "50to69", "70_and_above"))
  }
  if("age_group_3" %in% var_list_recovery_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_3 <- factor(aki_index_ckdonly_akionly$age_group_3,levels=c("below_50", "50to69", "70to79", "80plus"))
  }
  if("age_group_4" %in% var_list_recovery_ckdonly_akionly) {
    aki_index_ckdonly_akionly$age_group_4 <- factor(aki_index_ckdonly_akionly$age_group_4,levels=c("below_50", "50_and_above"))
  }
  
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
      if(isTRUE(restrict_models)) {
        ckdonly_model <- ckdonly_model[ckdonly_model %in% restrict_list]
      }
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
      coxph_ckdonly_vcov <- vcov(coxph_ckdonly)
      write.csv(coxph_ckdonly_vcov,file=file.path(dir.output,paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],"_vcov.csv")))
      coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly_akionly)
      ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
      invisible(gc())
    })
  }
  
  cat("\n=================\nRunning Models (minus CKD), Supp (Age), for AKI Recovery on CKD patients\n==============\n")
  for(i in 1:length(models_ckd_age_labels)) {
    cat(paste0("\nGenerating ", models_ckd_age_labels[i], " (Time to AKI Recovery, CKD Patients Only)..."))
    if((i %in% c(1:4) & "age_group_original" %in% var_list_recovery_ckdonly_akionly) |
       (i %in% c(5:8) & "age_group_2" %in% var_list_recovery_ckdonly_akionly) |
       (i %in% c(9:12) & "age_group_3" %in% var_list_recovery_ckdonly_akionly) |
       (i %in% c(13:16) & "age_group_4" %in% var_list_recovery_ckdonly_akionly)) {
      try({
        ckdonly_model <- c("severe","aki_kdigo_final",var_list_recovery_ckdonly_akionly)
        if(isTRUE(restrict_models)) {
          ckdonly_model <- ckdonly_model[ckdonly_model %in% restrict_list]
        }
        ckdonly_model <- ckdonly_model[ckdonly_model %in% models_ckd_age[[i]]]
        ckdonlyCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(ckdonly_model,collapse="+")))
        message(paste("Formula for ", models_ckd_age_labels[i],": survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(ckdonly_model,collapse="+")))
        coxph_ckdonly <- survival::coxph(ckdonlyCoxPHFormula, data=aki_index_ckdonly_akionly)
        coxph_ckdonly_summ <- summary(coxph_ckdonly) 
        print(coxph_ckdonly_summ)
        coxph_ckdonly_hr <- cbind(coxph_ckdonly_summ$coefficients,coxph_ckdonly_summ$conf.int)[,-c(6,7)]
        coxph_ckdonly_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_ckdonly_summ$logtest,coxph_ckdonly_summ$sctest,coxph_ckdonly_summ$waldtest))
        coxph_ckdonly_stats2 <- rbind(data.table::as.data.table(coxph_ckdonly_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_ckdonly_summ$rsq,keep.rownames = T))
        write.csv(coxph_ckdonly_hr,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],".csv")),row.names=TRUE)
        write.csv(coxph_ckdonly_stats1,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],"_teststats.csv")),row.names=FALSE,col.names = FALSE)
        write.csv(coxph_ckdonly_stats2,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],"_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
        coxph_ckdonly_vcov <- vcov(coxph_ckdonly)
        write.csv(coxph_ckdonly_vcov,file=file.path(dir.output,paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],"_vcov.csv")))
        coxph_ckdonly_plot <- survminer::ggforest(coxph_ckdonly,data=aki_index_ckdonly_akionly)
        ggplot2::ggsave(filename=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKDOnly_AKIOnly_CoxPH_",models_ckd_age_labels[i],".png")),plot=print(coxph_ckdonly_plot),width=20,height=20,units="cm")
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_ckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
  }
  
  return(list("aki_index_ckdonly" = aki_index_ckdonly, "aki_index_ckdonly_akionly" = aki_index_ckdonly_akionly,"var_list_recovery_ckdonly_akionly" = var_list_recovery_ckdonly_akionly, "var_list_death_ckdonly_akionly" = var_list_death_ckdonly_akionly, "var_list_death_ckdonly" = var_list_death_ckdonly))
}
