#' Generated cumulative incidence curves and Fine-Gray models for competing risks of either recovery / new onset CKD (event of interest) against mortality
#' 
run_cic_analysis <- function(currSiteId,aki_index_recovery,aki_index_nonckd,aki_index_nonckd_akionly,aki_index_ckdonly_akionly,
                             var_list_recovery_all,
                             var_list_new_ckd_nonckd_akionly,
                             var_list_recovery_nonckd_akionly,
                             var_list_new_ckd_nonckd,
                             var_list_recovery_ckdonly_akionly,use_custom_output = FALSE,use_custom_output_dir = '/4ceData/Output') {
  if(isTRUE(use_custom_output)) {
    dir.output <- use_custom_output_dir
  } else {
    dir.output <- getProjectOutputDirectory()
  }
  
  cat("\nFirst listing models\n")
  model1 <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld")
  model2 <- c("age_group","sex","severe","bronchiectasis","copd","rheum","vte")
  model3 <- c("age_group","sex","severe","COAGA","COAGB","covid_rx")
  model4 <- c("age_group","sex","severe","aki_kdigo_final","ckd","acei_arb_preexposure")

  # Additional models as requested by Lancet reviewers
  # A: include all variables in stepwise manner
  model_2a <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_3a <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_4a <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # B: Split by antiviral type
  model_3b <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral")
  model_4b <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  
  # C: Adjust for timing of pre-admission Cr
  model_1c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld")
  model_2c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_3c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_4c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  model_1d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld")
  model_2d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_3d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_4d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models <- list(model1,model2,model3,model4,model_2a,model_3a,model_4a,model_3b,model_4b,model_1c,model_2c,model_3c,model_4c,model_1d,model_2d,model_3d,model_4d)
  models_labels <- c("Model1","Model2","Model3","Model4","Model_2A","Model_3A","Model_4A","Model_3B","Model_4B","Model_1C","Model_2C","Model_3C","Model_4C","Model_1D","Model_2D","Model_3D","Model_4D")
  
  # Now create a new list where the models are modified to not include ckd
  model1 <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld")
  model2 <- c("age_group","sex","severe","bronchiectasis","copd","rheum","vte")
  model3 <- c("age_group","sex","severe","COAGA","COAGB","covid_rx")
  model4 <- c("age_group","sex","severe","aki_kdigo_final","acei_arb_preexposure")
  model_2a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_3a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_4a <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_3b <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral")
  model_4b <- c("age_group","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_1c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld")
  model_2c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte")
  model_3c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx")
  model_4c <- c("age_group","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_nockd <- list(model1,model2,model3,model4,model_2a,model_3a,model_4a,model_3b,model_4b,model_1c,model_2c,model_3c,model_4c)
  models_nockd_labels <- c("Model1","Model2","Model3","Model4","Model_2A","Model_3A","Model_4A","Model_3B","Model_4B","Model_1C","Model_2C","Model_3C","Model_4C")
  
  var_list_recovery_all_aki <- var_list_recovery_all
  
  cat("\nProceeding to transforming tables\n")
  recovery_all_cic <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
    recover_1.25x == 1 ~ 1,
    recover_1.25x == 0 & deceased == 1 ~ 2,
    TRUE ~ 0
  )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
    event == 0 ~ time_to_death_km,
    event == 1 ~ time_to_ratio1.25,
    event == 2 ~ time_to_death_km
  )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T)
  recovery_all_cic$event <- factor(recovery_all_cic$event,levels=c(0,1,2),labels = c("Censored","Recovery","Death"))
  
  recovery_nonckd_cic <- aki_index_nonckd_akionly %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
    recover_1.25x == 1 ~ 1,
    recover_1.25x == 0 & deceased == 1 ~ 2,
    TRUE ~ 0
  )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
    event == 0 ~ time_to_death_km,
    event == 1 ~ time_to_ratio1.25,
    event == 2 ~ time_to_death_km
  )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T)
  recovery_nonckd_cic$event <- factor(recovery_nonckd_cic$event,levels=c(0,1,2),labels = c("Censored","Recovery","Death"))
  
  recovery_ckd_cic <- aki_index_ckdonly_akionly %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
    recover_1.25x == 1 ~ 1,
    recover_1.25x == 0 & deceased == 1 ~ 2,
    TRUE ~ 0
  )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
    event == 0 ~ time_to_death_km,
    event == 1 ~ time_to_ratio1.25,
    event == 2 ~ time_to_death_km
  )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T)
  recovery_ckd_cic$event <- factor(recovery_ckd_cic$event,levels=c(0,1,2),labels = c("Censored","Recovery","Death"))
  
  new_ckd_nonckd_cic <- aki_index_nonckd %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
    new_ckd == 1 ~ 1,
    new_ckd == 0 & deceased == 1 ~ 2,
    TRUE ~ 0
  )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
    event == 0 ~ time_to_death_km,
    event == 1 ~ time_to_new_ckd,
    event == 2 ~ time_to_death_km
  )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T)
  new_ckd_nonckd_cic$event <- factor(new_ckd_nonckd_cic$event,levels=c(0,1,2),labels = c("Censored","CKD Onset","Death"))
  
  new_ckd_nonckd_akionly_cic <- aki_index_nonckd_akionly %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
    new_ckd == 1 ~ 1,
    new_ckd == 0 & deceased == 1 ~ 2,
    TRUE ~ 0
  )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
    event == 0 ~ time_to_death_km,
    event == 1 ~ time_to_new_ckd,
    event == 2 ~ time_to_death_km
  )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T)
  new_ckd_nonckd_akionly_cic$event <- factor(new_ckd_nonckd_akionly_cic$event,levels=c(0,1,2),labels = c("Censored","CKD Onset","Death"))
  
  # ==========================
  # Cumulative Incidence Plots
  # ==========================
  cat("\nNow preparing plots\n")
  # ===============
  recovery_formula <- as.formula("survival::Surv(time_to_ratio1.25,event) ~ aki_kdigo_final")
  cuminc_recovery_all <- tidycmprsk::cuminc(recovery_formula,data=recovery_all_cic)
  # cuminc_recovery_all %>% tidycmprsk::autoplot(conf.int=T)
  cuminc_recovery_nonckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_nonckd_cic)
  cuminc_recovery_ckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_ckd_cic)
  
  cuminc_recovery_all_tidy <- tidycmprsk::tidy(cuminc_recovery_all)
  cuminc_recovery_all_stats <- tidycmprsk::glance(cuminc_recovery_all) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_all_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_KDIGO_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_all_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_KDIGO_stats.csv")),row.names=F)
  
  cuminc_recovery_nonckd_tidy <- tidycmprsk::tidy(cuminc_recovery_nonckd)
  cuminc_recovery_nonckd_stats <- tidycmprsk::glance(cuminc_recovery_nonckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_KDIGO_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_nonckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_KDIGO_Stats.csv")),row.names=F)
  
  cuminc_recovery_ckd_tidy <- tidycmprsk::tidy(cuminc_recovery_ckd)
  cuminc_recovery_ckd_stats <- tidycmprsk::glance(cuminc_recovery_ckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_ckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_KDIGO_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_ckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_KDIGO_Stats.csv")),row.names=F)
  
  # ===============
  recovery_formula <- as.formula("survival::Surv(time_to_ratio1.25,event) ~ severe")
  cuminc_recovery_all <- tidycmprsk::cuminc(recovery_formula,data=recovery_all_cic)
  cuminc_recovery_nonckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_nonckd_cic)
  cuminc_recovery_ckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_ckd_cic)
  
  cuminc_recovery_all_tidy <- tidycmprsk::tidy(cuminc_recovery_all)
  cuminc_recovery_all_stats <- tidycmprsk::glance(cuminc_recovery_all) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_all_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_Severe_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_all_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_Severe_stats.csv")),row.names=F)
  
  cuminc_recovery_nonckd_tidy <- tidycmprsk::tidy(cuminc_recovery_nonckd)
  cuminc_recovery_nonckd_stats <- tidycmprsk::glance(cuminc_recovery_nonckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_Severe_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_nonckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_Severe_Stats.csv")),row.names=F)
  
  cuminc_recovery_ckd_tidy <- tidycmprsk::tidy(cuminc_recovery_ckd)
  cuminc_recovery_ckd_stats <- tidycmprsk::glance(cuminc_recovery_ckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_ckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_Severe_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_ckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_Severe_Stats.csv")),row.names=F)
  
  # ===============
  recovery_formula <- as.formula("survival::Surv(time_to_ratio1.25,event) ~ ckd")
  cuminc_recovery_all <- tidycmprsk::cuminc(recovery_formula,data=recovery_all_cic)
  
  cuminc_recovery_all_tidy <- tidycmprsk::tidy(cuminc_recovery_all)
  cuminc_recovery_all_stats <- tidycmprsk::glance(cuminc_recovery_all) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
  write.csv(cuminc_recovery_all_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_CKDvsNonCKD_Plot.csv")),row.names=F)
  write.csv(cuminc_recovery_all_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_CKDvsNonCKD_stats.csv")),row.names=F)
  
  # ===============
  
  newckd_formula <- as.formula("survival::Surv(time_to_new_ckd,event) ~ aki_kdigo_final")
  cuminc_new_ckd_nonckd <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_cic)
  cuminc_new_ckd_nonckd_akionly <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_akionly_cic)
  
  cuminc_new_ckd_nonckd_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd)
  write.csv(cuminc_new_ckd_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_All_NonCKD_KDIGO_Plot.csv")),row.names=F)
  cuminc_new_ckd_nonckd_akionly_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd_akionly)
  write.csv(cuminc_new_ckd_nonckd_akionly_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_NonCKD_AKIOnly_KDIGO_Plot.csv")),row.names=F)
  
  
  newckd_formula <- as.formula("survival::Surv(time_to_new_ckd,event) ~ severe")
  cuminc_new_ckd_nonckd <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_cic)
  cuminc_new_ckd_nonckd_akionly <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_akionly_cic)
  
  cuminc_new_ckd_nonckd_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd)
  write.csv(cuminc_new_ckd_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_All_NonCKD_Severe_Plot.csv")),row.names=F)
  cuminc_new_ckd_nonckd_akionly_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd_akionly)
  write.csv(cuminc_new_ckd_nonckd_akionly_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_NonCKD_AKIOnly_Severe_Plot.csv")),row.names=F)
  
  # ==================
  
  for(i in 1:length(models_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_labels[i], ") (time to recovery, Non-CKD + CKD, AKI patients only)..."))
    try({
      recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_all_aki)
      recovery_model <- recovery_model[recovery_model %in% models[[i]]]
      recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      message(paste("Formula for ", models_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_all_cic,failcode = "Recovery")
      
      FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
      FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_labels[i],".csv")),row.names=F)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_labels[i],"_stats.csv")),row.names=F)
      invisible(gc())
    })
  }
  
  for(i in 1:length(models_nockd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_nockd_labels[i], ") (time to recovery, Non-CKD Only, AKI patients only)..."))
    try({
      recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_nonckd_akionly)
      recovery_model <- recovery_model[recovery_model %in% models_nockd[[i]]]
      recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      message(paste("Formula for ", models_nockd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_nonckd_cic,failcode = "Recovery")
      FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
      FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nockd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nockd_labels[i],"_stats.csv")),row.names=F)
      invisible(gc())
    })
  }
  
  for(i in 1:length(models_nockd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_nockd_labels[i], ") (time to recovery, CKD Only, AKI patients only)..."))
    try({
      recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_ckdonly_akionly)
      recovery_model <- recovery_model[recovery_model %in% models_nockd[[i]]]
      recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      message(paste("Formula for ", models_nockd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_ckd_cic,failcode = "Recovery")
      FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
      FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_nockd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_nockd_labels[i],"_stats.csv")),row.names=F)
      invisible(gc())
    })
  }
  
  for(i in 1:length(models_nockd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models_nockd (", models_nockd_labels[i], ") (time to New Onset CKD, Non-CKD Only, All)..."))
    try({
      newckd_model <- c("severe","aki_kdigo_final",var_list_new_ckd_nonckd_akionly)
      newckd_model <- newckd_model[newckd_model %in% models_nockd[[i]]]
      newckdFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      message(paste("Formula for ", models_nockd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      FineGray_newckd <- tidycmprsk::crr(newckdFineGrayFormula, data=new_ckd_nonckd_cic,failcode = "CKD Onset")
      FineGray_newckd_summ <- tidycmprsk::tidy(FineGray_newckd,exponentiate=T,conf.int=T) 
      FineGray_newckd_stats <- tidycmprsk::glance(FineGray_newckd)
      write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_FineGray_",models_nockd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_FineGray_",models_nockd_labels[i],"_stats.csv")),row.names=F)
      invisible(gc())
    })
  }
  
  for(i in 1:length(models_nockd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models_nockd (", models_nockd_labels[i], ") (time to New Onset CKD, Non-CKD Only, AKI patients only)..."))
    try({
      newckd_model <- c("severe","aki_kdigo_final",var_list_new_ckd_nonckd_akionly)
      newckd_model <- newckd_model[newckd_model %in% models_nockd[[i]]]
      newckdFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      message(paste("Formula for ", models_nockd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      FineGray_newckd <- tidycmprsk::crr(newckdFineGrayFormula, data=new_ckd_nonckd_akionly_cic,failcode = "CKD Onset")
      FineGray_newckd_summ <- tidycmprsk::tidy(FineGray_newckd,exponentiate=T,conf.int=T) 
      FineGray_newckd_stats <- tidycmprsk::glance(FineGray_newckd)
      write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_AKIOnly_FineGray_",models_nockd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_AKIOnly_FineGray_",models_nockd_labels[i],"_stats.csv")),row.names=F)
      invisible(gc())
    })
  }
  
}
