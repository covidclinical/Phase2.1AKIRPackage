#' Generated cumulative incidence curves and Fine-Gray models for competing risks of either recovery / new onset CKD (event of interest) against mortality
#' 
run_cic_analysis <- function(currSiteId,aki_index_recovery,aki_index_nonckd,aki_index_nonckd_akionly,aki_index_ckdonly_akionly,
                             var_list_recovery_all,
                             var_list_new_ckd_nonckd_akionly,
                             var_list_recovery_nonckd_akionly,
                             var_list_new_ckd_nonckd,
                             var_list_recovery_ckdonly_akionly,use_custom_output = FALSE,use_custom_output_dir = '/4ceData/Output',ckd_present=TRUE) {
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
  model_4d <- c("age_group","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_nonckd <- list(model1,model2,model3,model4,model_2a,model_3a,model_4a,model_3b,model_4b,model_1c,model_2c,model_3c,model_4c,model_4d)
  models_nonckd_labels <- c("Model1","Model2","Model3","Model4","Model_2A","Model_3A","Model_4A","Model_3B","Model_4B","Model_1C","Model_2C","Model_3C","Model_4C","Model_4D")
  
  # Create age group supplementary models
  # E: Use Models 4A,B,C,D but all the age groups as strata
  model_1e <- c("age_group_original","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_2e <-  c("age_group_original","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_3e <- c("age_group_original","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_4e <- c("age_group_original","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # F: Use Models 4A,B,C,D but split age to three strata (< 50, 50-69, >=70)
  model_1f <- c("age_group_2","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_2f <-  c("age_group_2","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_3f <- c("age_group_2","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_4f <- c("age_group_2","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # G: Use Models 4A,B,C,D but split age by strata (<50,50-69,70-79,>=80)
  model_1g <- c("age_group_3","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_2g <-  c("age_group_3","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_3g <- c("age_group_3","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_4g <- c("age_group_3","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # H: Use Models 4A,B,C,D but split age by >=50 instead of >= 70
  model_1h <- c("age_group_4","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_2h <-  c("age_group_4","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_3h <- c("age_group_4","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_4h <- c("age_group_4","sex","severe","ckd_stage","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_age <- list(model_1e,model_2e,model_3e,model_4e,model_1f,model_2f,model_3f,model_4f,model_1g,model_2g,model_3g,model_4g,model_1h,model_2h,model_3h,model_4h)
  models_age_labels <- c("Model_1E","Model_2E","Model_3E","Model_4E",
                         "Model_1F","Model_2F","Model_3F","Model_4F",
                         "Model_1G","Model_2G","Model_3G","Model_4G",
                         "Model_1H","Model_2H","Model_3H","Model_4H")
  
  # E: Use Models 4A,B,C,D but all the age groups as strata
  model_ckd_subgroup_1e <- c("age_group_original","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2e <-  c("age_group_original","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3e <- c("age_group_original","sex","severe","preadmit_cr_period","aki_kdigo_final","ckd","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4e <- c("age_group_original","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # F: Use Models 4A,B,C,D but split age to three strata (< 50, 50-69, >=70)
  model_ckd_subgroup_1f <- c("age_group_2","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2f <-  c("age_group_2","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3f <- c("age_group_2","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4f <- c("age_group_2","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # G: Use Models 4A,B,C,D but split age by strata (<50,50-69,70-79,>=80)
  model_ckd_subgroup_1g <- c("age_group_3","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2g <-  c("age_group_3","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3g <- c("age_group_3","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4g <- c("age_group_3","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  # H: Use Models 4A,B,C,D but split age by >=50 instead of >= 70
  model_ckd_subgroup_1h <- c("age_group_4","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_2h <-  c("age_group_4","sex","severe","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","remdesivir","covidviral","acei_arb_preexposure")
  model_ckd_subgroup_3h <- c("age_group_4","sex","severe","preadmit_cr_period","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  model_ckd_subgroup_4h <- c("age_group_4","sex","severe","ckd_stage","aki_kdigo_final","htn","ihd","cld","bronchiectasis","copd","rheum","vte","COAGA","COAGB","covid_rx","acei_arb_preexposure")
  
  models_ckd_age <- list(model_ckd_subgroup_1e,model_ckd_subgroup_2e,model_ckd_subgroup_3e,model_ckd_subgroup_4e,
                         model_ckd_subgroup_1f,model_ckd_subgroup_2f,model_ckd_subgroup_3f,model_ckd_subgroup_4f,
                         model_ckd_subgroup_1g,model_ckd_subgroup_2g,model_ckd_subgroup_3g,model_ckd_subgroup_4g,
                         model_ckd_subgroup_1h,model_ckd_subgroup_2h,model_ckd_subgroup_3h,model_ckd_subgroup_4h)
  
  models_nonckd_age <- list(model_ckd_subgroup_1e,model_ckd_subgroup_2e,model_ckd_subgroup_3e,
                         model_ckd_subgroup_1f,model_ckd_subgroup_2f,model_ckd_subgroup_3f,
                         model_ckd_subgroup_1g,model_ckd_subgroup_2g,model_ckd_subgroup_3g,
                         model_ckd_subgroup_1h,model_ckd_subgroup_2h,model_ckd_subgroup_3h)
  
  models_ckd_age_labels <- c("Model_1E","Model_2E","Model_3E","Model_4E",
                             "Model_1F","Model_2F","Model_3F","Model_4F",
                             "Model_1G","Model_2G","Model_3G","Model_4G",
                             "Model_1H","Model_2H","Model_3H","Model_4H")
  
  models_nonckd_age_labels <- c("Model_1E","Model_2E","Model_3E",
                             "Model_1F","Model_2F","Model_3F",
                             "Model_1G","Model_2G","Model_3G",
                             "Model_1H","Model_2H","Model_3H")
  
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
  
  if(isTRUE(ckd_present)) {
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
    
  } else { # the analysis is equivalent to a non-CKD subgroup already
    
    recovery_nonckd_cic <- recovery_all_cic
    recovery_ckd_cic <- NULL
    
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
    
    new_ckd_nonckd_akionly_cic <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
      new_ckd == 1 ~ 1,
      new_ckd == 0 & deceased == 1 ~ 2,
      TRUE ~ 0
    )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
      event == 0 ~ time_to_death_km,
      event == 1 ~ time_to_new_ckd,
      event == 2 ~ time_to_death_km
    )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T)
    new_ckd_nonckd_akionly_cic$event <- factor(new_ckd_nonckd_akionly_cic$event,levels=c(0,1,2),labels = c("Censored","CKD Onset","Death"))
  }
  
  
  # ==========================
  # Cumulative Incidence Plots
  # ==========================
  cat("\nNow preparing plots\n")
  # ===============
  recovery_formula <- as.formula("survival::Surv(time_to_ratio1.25,event) ~ aki_kdigo_final")
  try({
    cat("\nCreating CumInc curve for AKI Recovery (All) - KDIGO\n")
    cuminc_recovery_all <- tidycmprsk::cuminc(recovery_formula,data=recovery_all_cic)
    # cuminc_recovery_all %>% tidycmprsk::autoplot(conf.int=T)
    cuminc_recovery_all_tidy <- tidycmprsk::tidy(cuminc_recovery_all)
    cuminc_recovery_all_stats <- tidycmprsk::glance(cuminc_recovery_all) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
    write.csv(cuminc_recovery_all_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_KDIGO_Plot.csv")),row.names=F)
    write.csv(cuminc_recovery_all_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_KDIGO_stats.csv")),row.names=F)
  })
  
  try({
    cat("\nCreating CumInc curve for AKI Recovery (Non-CKD only) - KDIGO\n")
    cuminc_recovery_nonckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_nonckd_cic)
    cuminc_recovery_nonckd_tidy <- tidycmprsk::tidy(cuminc_recovery_nonckd)
    cuminc_recovery_nonckd_stats <- tidycmprsk::glance(cuminc_recovery_nonckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
    write.csv(cuminc_recovery_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_KDIGO_Plot.csv")),row.names=F)
    write.csv(cuminc_recovery_nonckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_KDIGO_Stats.csv")),row.names=F)
  })
  if(isTRUE(ckd_present)) {
    try({
      cat("\nCreating CumInc curve for AKI Recovery (CKD only) - KDIGO\n")
      cuminc_recovery_ckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_ckd_cic)
      cuminc_recovery_ckd_tidy <- tidycmprsk::tidy(cuminc_recovery_ckd)
      cuminc_recovery_ckd_stats <- tidycmprsk::glance(cuminc_recovery_ckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
      write.csv(cuminc_recovery_ckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_KDIGO_Plot.csv")),row.names=F)
      write.csv(cuminc_recovery_ckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_KDIGO_Stats.csv")),row.names=F)
    })
  }
  
  # ===============
  recovery_formula <- as.formula("survival::Surv(time_to_ratio1.25,event) ~ severe")
  
  try({
    cat("\nCreating CumInc curve for AKI Recovery (All) - COVID-19 Severity\n")
    cuminc_recovery_all <- tidycmprsk::cuminc(recovery_formula,data=recovery_all_cic)
    cuminc_recovery_all_tidy <- tidycmprsk::tidy(cuminc_recovery_all)
    cuminc_recovery_all_stats <- tidycmprsk::glance(cuminc_recovery_all) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
    write.csv(cuminc_recovery_all_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_Severe_Plot.csv")),row.names=F)
    write.csv(cuminc_recovery_all_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_Severe_stats.csv")),row.names=F)
  })
  try({
    cat("\nCreating CumInc curve for AKI Recovery (Non-CKD Only) - COVID-19 Severity\n")
    cuminc_recovery_nonckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_nonckd_cic)
    cuminc_recovery_nonckd_tidy <- tidycmprsk::tidy(cuminc_recovery_nonckd)
    cuminc_recovery_nonckd_stats <- tidycmprsk::glance(cuminc_recovery_nonckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
    write.csv(cuminc_recovery_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_Severe_Plot.csv")),row.names=F)
    write.csv(cuminc_recovery_nonckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_NonCKD_Severe_Stats.csv")),row.names=F)
  })
  if(isTRUE(ckd_present)) {
    try({
      cat("\nCreating CumInc curve for AKI Recovery (CKD Only) - COVID-19 Severity\n")
      cuminc_recovery_ckd <- tidycmprsk::cuminc(recovery_formula,data=recovery_ckd_cic)
      cuminc_recovery_ckd_tidy <- tidycmprsk::tidy(cuminc_recovery_ckd)
      cuminc_recovery_ckd_stats <- tidycmprsk::glance(cuminc_recovery_ckd) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
      write.csv(cuminc_recovery_ckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_Severe_Plot.csv")),row.names=F)
      write.csv(cuminc_recovery_ckd_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_CKD_Severe_Stats.csv")),row.names=F)
    })
  }
  
  # ===============
  if(isTRUE(ckd_present)) {
    recovery_formula <- as.formula("survival::Surv(time_to_ratio1.25,event) ~ ckd")
    try({
      cat("\nCreating CumInc curve for AKI Recovery (All) - CKD\n")
      cuminc_recovery_all <- tidycmprsk::cuminc(recovery_formula,data=recovery_all_cic)
      cuminc_recovery_all_tidy <- tidycmprsk::tidy(cuminc_recovery_all)
      cuminc_recovery_all_stats <- tidycmprsk::glance(cuminc_recovery_all) %>% tidyr::pivot_longer(tidyr::everything(),names_to = c(".value", "outcome_id"),names_pattern = "(.*)_(.*)")
      write.csv(cuminc_recovery_all_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_CKDvsNonCKD_Plot.csv")),row.names=F)
      write.csv(cuminc_recovery_all_stats,file=file.path(dir.output, paste0(currSiteId, "_CumInc_Recovery_All_AKI_CKDvsNonCKD_stats.csv")),row.names=F)
    })
  }
  
  # ===============
  
  newckd_formula <- as.formula("survival::Surv(time_to_new_ckd,event) ~ aki_kdigo_final")
  
  try({
    cat("\nCreating CumInc curve for Time to New CKD Onset (Non-CKD only) - KDIGO\n")
    cuminc_new_ckd_nonckd <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_cic)
    cuminc_new_ckd_nonckd_akionly <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_akionly_cic)
    cuminc_new_ckd_nonckd_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd)
    write.csv(cuminc_new_ckd_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_All_NonCKD_KDIGO_Plot.csv")),row.names=F)
    cuminc_new_ckd_nonckd_akionly_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd_akionly)
    write.csv(cuminc_new_ckd_nonckd_akionly_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_NonCKD_AKIOnly_KDIGO_Plot.csv")),row.names=F)
  })
  
  
  try({
    cat("\nCreating CumInc curve for Time to New CKD Onset (Non-CKD only) - COVID-19 Severity\n")
    newckd_formula <- as.formula("survival::Surv(time_to_new_ckd,event) ~ severe")
    cuminc_new_ckd_nonckd <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_cic)
    cuminc_new_ckd_nonckd_akionly <- tidycmprsk::cuminc(newckd_formula,data=new_ckd_nonckd_akionly_cic)
    
    cuminc_new_ckd_nonckd_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd)
    write.csv(cuminc_new_ckd_nonckd_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_All_NonCKD_Severe_Plot.csv")),row.names=F)
    cuminc_new_ckd_nonckd_akionly_tidy <- tidycmprsk::tidy(cuminc_new_ckd_nonckd_akionly)
    write.csv(cuminc_new_ckd_nonckd_akionly_tidy,file=file.path(dir.output, paste0(currSiteId, "_CumInc_NewCKD_NonCKD_AKIOnly_Severe_Plot.csv")),row.names=F)
  })
  
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
      write.csv(FineGray_recovery_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_labels[i],"_stats.csv")),row.names=F)
      
      FineGray_recovery_vcov <- vcov(FineGray_recovery)
      write.csv(FineGray_recovery_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_labels[i],"_vcov.csv")))
      
      invisible(gc())
    })
  }
  
  cat("\n============================\nFine-Gray, Age Supp Models, AKI Recovery\n============================\n")
  for(i in 1:length(models_age_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_age_labels[i], ") (time to recovery, Non-CKD + CKD, AKI patients only)..."))
    if((i %in% c(1:4) & "age_group_original" %in% var_list_recovery_all_aki) |
       (i %in% c(5:8) & "age_group_2" %in% var_list_recovery_all_aki) |
       (i %in% c(9:12) & "age_group_3" %in% var_list_recovery_all_aki) |
       (i %in% c(13:16) & "age_group_4" %in% var_list_recovery_all_aki)) {
      try({
        recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_all_aki)
        recovery_model <- recovery_model[recovery_model %in% models_age[[i]]]
        recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
        message(paste("Formula for ", models_age_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
        FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_all_cic,failcode = "Recovery")
        FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
        FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
        
        write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_age_labels[i],".csv")),row.names=F)
        write.csv(FineGray_recovery_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_age_labels[i],"_stats.csv")),row.names=F)
        
        FineGray_recovery_vcov <- vcov(FineGray_recovery)
        write.csv(FineGray_recovery_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_All_AKI_FineGray_",models_age_labels[i],"_vcov.csv")))
        
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
  }
  
  for(i in 1:length(models_nonckd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_nonckd_labels[i], ") (time to recovery, Non-CKD Only, AKI patients only)..."))
    try({
      recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_nonckd_akionly)
      recovery_model <- recovery_model[recovery_model %in% models_nonckd[[i]]]
      recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      message(paste("Formula for ", models_nonckd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
      FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_nonckd_cic,failcode = "Recovery")
      FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
      FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
      write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nonckd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_recovery_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nonckd_labels[i],"_stats.csv")),row.names=F)
      FineGray_recovery_vcov <- vcov(FineGray_recovery)
      write.csv(FineGray_recovery_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nonckd_labels[i],"_vcov.csv")))
      invisible(gc())
    })
  }
  
  cat("\n============================\nFine-Gray, Age Supp Models, AKI Recovery, Non-CKD only\n============================\n")
  for(i in 1:length(models_nonckd_age_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_nonckd_age_labels[i], ") (time to recovery, Non-CKD, AKI patients only)..."))
    if((i %in% c(1:3) & "age_group_original" %in% var_list_recovery_nonckd_akionly) |
       (i %in% c(4:6) & "age_group_2" %in% var_list_recovery_nonckd_akionly) |
       (i %in% c(7:9) & "age_group_3" %in% var_list_recovery_nonckd_akionly) |
       (i %in% c(10:12) & "age_group_4" %in% var_list_recovery_nonckd_akionly)) {
      try({
        recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_nonckd_akionly)
        recovery_model <- recovery_model[recovery_model %in% models_nonckd_age[[i]]]
        recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
        message(paste("Formula for ", models_nonckd_age_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
        FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_nonckd_cic,failcode = "Recovery")
        FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
        FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
        write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nonckd_age_labels[i],".csv")),row.names=F)
        write.csv(FineGray_recovery_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nonckd_age_labels[i],"_stats.csv")),row.names=F)
        FineGray_recovery_vcov <- vcov(FineGray_recovery)
        write.csv(FineGray_recovery_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_NonCKD_AKI_FineGray_",models_nonckd_age_labels[i],"_vcov.csv")))
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_nonckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
  }
  
  
  if(isTRUE(ckd_present)) {
    for(i in 1:length(models_nonckd_labels)) {
      cat(paste0("\nGenerating Fine-Gray models (", models_nonckd_labels[i], ") (time to recovery, CKD Only, AKI patients only)..."))
      try({
        recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_ckdonly_akionly)
        recovery_model <- recovery_model[recovery_model %in% models_nonckd[[i]]]
        recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
        message(paste("Formula for ", models_nonckd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
        FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_ckd_cic,failcode = "Recovery")
        FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
        FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
        write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_nonckd_labels[i],".csv")),row.names=F)
        write.csv(FineGray_recovery_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_nonckd_labels[i],"_stats.csv")),row.names=F)
        FineGray_recovery_vcov <- vcov(FineGray_recovery)
        write.csv(FineGray_recovery_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_nonckd_labels[i],"_vcov.csv")))
        invisible(gc())
      })
    }
    cat("\n============================\nFine-Gray, Age Supp Models, AKI Recovery, CKD only\n============================\n")
    for(i in 1:length(models_ckd_age_labels)) {
      cat(paste0("\nGenerating Fine-Gray models (", models_ckd_age_labels[i], ") (time to recovery, CKD, AKI patients only)..."))
      if((i %in% c(1:4) & "age_group_original" %in% var_list_recovery_ckdonly_akionly) |
         (i %in% c(5:8) & "age_group_2" %in% var_list_recovery_ckdonly_akionly) |
         (i %in% c(9:12) & "age_group_3" %in% var_list_recovery_ckdonly_akionly) |
         (i %in% c(13:16) & "age_group_4" %in% var_list_recovery_ckdonly_akionly)) {
        try({
          recovery_model <- c("severe","aki_kdigo_final",var_list_recovery_ckdonly_akionly)
          recovery_model <- recovery_model[recovery_model %in% models_ckd_age[[i]]]
          recoveryFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
          message(paste("Formula for ", models_ckd_age_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(recovery_model,collapse="+")))
          FineGray_recovery <- tidycmprsk::crr(recoveryFineGrayFormula, data=recovery_ckd_cic,failcode = "Recovery")
          FineGray_recovery_summ <- tidycmprsk::tidy(FineGray_recovery,exponentiate=T,conf.int=T) 
          FineGray_recovery_stats <- tidycmprsk::glance(FineGray_recovery)
          write.csv(FineGray_recovery_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_ckd_age_labels[i],".csv")),row.names=F)
          write.csv(FineGray_recovery_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_ckd_age_labels[i],"_stats.csv")),row.names=F)
          FineGray_recovery_vcov <- vcov(FineGray_recovery)
          write.csv(FineGray_recovery_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_Recovery_CKD_FineGray_",models_ckd_age_labels[i],"_vcov.csv")))
          invisible(gc())
        })
      } else {
        cat("\nThe age group variable for Model ",models_ckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
      }
    }
  }
  
  for(i in 1:length(models_nonckd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_nonckd_labels[i], ") (time to New Onset CKD, Non-CKD Only, All)..."))
    try({
      newckd_model <- c("severe","aki_kdigo_final",var_list_new_ckd_nonckd)
      newckd_model <- newckd_model[newckd_model %in% models_nonckd[[i]]]
      newckdFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      message(paste("Formula for ", models_nonckd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      FineGray_newckd <- tidycmprsk::crr(newckdFineGrayFormula, data=new_ckd_nonckd_cic,failcode = "CKD Onset")
      FineGray_newckd_summ <- tidycmprsk::tidy(FineGray_newckd,exponentiate=T,conf.int=T) 
      FineGray_newckd_stats <- tidycmprsk::glance(FineGray_newckd)
      write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_FineGray_",models_nonckd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_newckd_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_FineGray_",models_nonckd_labels[i],"_stats.csv")),row.names=F)
      FineGray_newckd_vcov <- vcov(FineGray_newckd)
      write.csv(FineGray_newckd_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_NonCKD_FineGray_",models_nonckd_labels[i],"_vcov.csv")))
      invisible(gc())
    })
  }
  
  cat("\n=================\nRunning Models (minus CKD), Supp (Age), for New Onset CKD, Non-CKD patients\n==============\n")
  for(i in 1:length(models_nonckd_age_labels)) {
    cat(paste0("\nGenerating ", models_nonckd_age_labels[i], " (Time to New Onset CKD, Non-CKD Only, All)..."))
    if((i %in% c(1:3) & "age_group_original" %in% var_list_new_ckd_nonckd) |
       (i %in% c(4:6) & "age_group_2" %in% var_list_new_ckd_nonckd) |
       (i %in% c(7:9) & "age_group_3" %in% var_list_new_ckd_nonckd) |
       (i %in% c(10:12) & "age_group_4" %in% var_list_new_ckd_nonckd)) {
      try({
        newckd_model <- c("severe","aki_kdigo_final",var_list_new_ckd_nonckd)
        newckd_model <- newckd_model[newckd_model %in% models_nonckd_age[[i]]]
        newckdFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
        message(paste("Formula for ", models_nonckd_age_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
        FineGray_newckd <- tidycmprsk::crr(newckdFineGrayFormula, data=new_ckd_nonckd_cic,failcode = "CKD Onset")
        FineGray_newckd_summ <- tidycmprsk::tidy(FineGray_newckd,exponentiate=T,conf.int=T) 
        FineGray_newckd_stats <- tidycmprsk::glance(FineGray_newckd)
        write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_FineGray_",models_nonckd_age_labels[i],".csv")),row.names=F)
        write.csv(FineGray_newckd_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_FineGray_",models_nonckd_age_labels[i],"_stats.csv")),row.names=F)
        FineGray_newckd_vcov <- vcov(FineGray_newckd)
        write.csv(FineGray_newckd_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_NonCKD_FineGray_",models_nonckd_age_labels[i],"_vcov.csv")))
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_nonckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
  }
  
  for(i in 1:length(models_nonckd_labels)) {
    cat(paste0("\nGenerating Fine-Gray models (", models_nonckd_labels[i], ") (time to New Onset CKD, Non-CKD Only, AKI patients only)..."))
    try({
      newckd_model <- c("severe","aki_kdigo_final",var_list_new_ckd_nonckd_akionly)
      newckd_model <- newckd_model[newckd_model %in% models_nonckd[[i]]]
      newckdFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      message(paste("Formula for ", models_nonckd_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
      FineGray_newckd <- tidycmprsk::crr(newckdFineGrayFormula, data=new_ckd_nonckd_akionly_cic,failcode = "CKD Onset")
      FineGray_newckd_summ <- tidycmprsk::tidy(FineGray_newckd,exponentiate=T,conf.int=T) 
      FineGray_newckd_stats <- tidycmprsk::glance(FineGray_newckd)
      write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_AKIOnly_FineGray_",models_nonckd_labels[i],".csv")),row.names=F)
      write.csv(FineGray_newckd_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_AKIOnly_FineGray_",models_nonckd_labels[i],"_stats.csv")),row.names=F)
      FineGray_newckd_vcov <- vcov(FineGray_newckd)
      write.csv(FineGray_newckd_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_NonCKD_AKIOnly_FineGray_",models_nonckd_labels[i],"_vcov.csv")))
      invisible(gc())
    })
  }
  
  cat("\n=================\nRunning Models (minus CKD), Supp (Age), for New Onset CKD, Non-CKD, AKI patients only\n==============\n")
  for(i in 1:length(models_nonckd_age_labels)) {
    cat(paste0("\nGenerating ", models_nonckd_age_labels[i], " (Time to New Onset CKD, Non-CKD Only, AKI Patients Only)..."))
    if((i %in% c(1:3) & "age_group_original" %in% var_list_new_ckd_nonckd_akionly) |
       (i %in% c(4:6) & "age_group_2" %in% var_list_new_ckd_nonckd_akionly) |
       (i %in% c(7:9) & "age_group_3" %in% var_list_new_ckd_nonckd_akionly) |
       (i %in% c(10:12) & "age_group_4" %in% var_list_new_ckd_nonckd_akionly)) {
      try({
        newckd_model <- c("severe","aki_kdigo_final",var_list_new_ckd_nonckd_akionly)
        newckd_model <- newckd_model[newckd_model %in% models_nonckd_age[[i]]]
        newckdFineGrayFormula <- as.formula(paste("survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
        message(paste("Formula for ", models_nonckd_age_labels[i],"survival::Surv(time_to_event,event) ~ ",paste(newckd_model,collapse="+")))
        FineGray_newckd <- tidycmprsk::crr(newckdFineGrayFormula, data=new_ckd_nonckd_akionly_cic,failcode = "CKD Onset")
        FineGray_newckd_summ <- tidycmprsk::tidy(FineGray_newckd,exponentiate=T,conf.int=T) 
        FineGray_newckd_stats <- tidycmprsk::glance(FineGray_newckd)
        write.csv(FineGray_newckd_summ,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_AKIOnly_FineGray_",models_nonckd_age_labels[i],".csv")),row.names=F)
        write.csv(FineGray_newckd_stats,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_New_CKD_NonCKD_AKIOnly_FineGray_",models_nonckd_age_labels[i],"_stats.csv")),row.names=F)
        FineGray_newckd_vcov <- vcov(FineGray_newckd)
        write.csv(FineGray_newckd_vcov,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_NewCKD_NonCKD_AKIOnly_FineGray_",models_nonckd_age_labels[i],"_vcov.csv")))
        invisible(gc())
      })
    } else {
      cat("\nThe age group variable for Model ",models_nonckd_age_labels[i]," is not available. Check if the lowest quantity in all levels meets the required threshhold of events.\n")
    }
  }
  
}
