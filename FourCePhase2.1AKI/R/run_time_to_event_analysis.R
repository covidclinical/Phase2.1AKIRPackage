#' Runs time-to-event analysis for the 4 pre-defined models - Kaplan-Meier curves, Cox PH models.
#' If specified, can also filter the output to pre-specified patients with pre-admission sCr values.
#' @param siteid currSiteId (site ID of the site as obtained earlier)
#' @param base_table peak_trend
#' @param aki_episodes aki_index
#' @param demog_table The demographics_filt equivalent
#' @param demog_toe Specify the appropriate demog_time_to_event table
#' @param demog_var_list Specify the appropriate demog_list list
#' @param comorbid_table The comorbid equivalent
#' @param comorbid_header The comorbid_list equivalent
#' @param kdigo_grade_table kdigo_grade
#' @param ckd_valid ckd_present (Whether at least one patient has CKD as a comorbidity)
#' @param coaga_valid coaga_present
#' @param coagb_valid coagb_present
#' @param covid_valid covid19antiviral_present (Novel COVID-19 antiviral use)
#' @param remdesivir_valid remdesivir_present (Remdesivir use)
#' @param acei_valid acei_present (ACE-inhibitor use)
#' @param arb_valid arb_present (ARB use)
#' @param med_coaga med_coaga_new (List of patients with COAGA)
#' @param med_coagb med_coagb_new (List of patients with COAGB)
#' @param med_covid19 med_covid19_new (List of patients with Novel COVID-19 antiviral/remdesivir use)
#' @param med_acearb med_acearb_chronic (List of patients previously on ACE-inhibitors/ARBs)
#' @param preadmit_cr_list Pre-specified patient list with pre-admission sCr labs. Default value is NULL
#' @param preadmit_only_analysis Specifies whether demographics are to be narrowed down to list specified by preadmit_cr_list. Default is FALSE.
#' @param obfuscation is_obfuscated
#' @param obfuscation_level obfuscation_value

run_time_to_event_analysis <- function(siteid, base_table, aki_episodes,
                                       demog_table, demog_toe = NULL, demog_var_list = NULL,
                                       comorbid_table = NULL,comorbid_header = NULL,
                                       kdigo_grade_table = NULL,
                                       ckd_valid,
                                       coaga_valid,coagb_valid,
                                       covid_valid,remdesivir_valid,
                                       acei_valid,arb_valid,
                                       med_coaga = NULL,med_coagb = NULL,med_covid19 = NULL,med_acearb = NULL,
                                       preadmit_cr_list = NULL,preadmit_only_analysis = FALSE,
                                       obfuscation,obfuscation_level,
                                       restrict_model_corr) {
  currSiteId <- siteid
  peak_trend <- base_table
  aki_index <- aki_episodes
  demographics <- demog_table
  comorbid <- comorbid_table
  comorbid_list <- comorbid_header
  kdigo_grade <- kdigo_grade_table
  ckd_present <- ckd_valid
  coaga_present <- coaga_valid
  coagb_present <- coaga_valid
  covid19antiviral_present <- covid_valid
  remdesivir_present <- remdesivir_valid
  acei_present <- acei_valid
  arb_present <- arb_valid
  med_coaga_new <- med_coaga
  med_coagb_new <- med_coagb
  med_covid19_new <- med_covid19
  med_acearb_chronic <- med_acearb
  
  patients_with_preadmit_cr <- preadmit_cr_list
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  
  restrict_models <- restrict_models_corr
  
  file_prefix <- ""
  if(isTRUE(preadmit_only_analysis)) {
    file_prefix <- "_PreAdmitCrOnly_"
    currSiteId <- paste0(currSiteId,file_prefix)
    peak_trend <- peak_trend[peak_trend$patient_id %in% patients_with_preadmit_cr,]
    aki_index <- aki_index[aki_index$patient_id %in% patients_with_preadmit_cr,]
    demographics <- demographics[demographics$patient_id %in% patients_with_preadmit_cr,]
    try({comorbid <- comorbid[comorbid$patient_id %in% patients_with_preadmit_cr,]})
    try({kdigo_grade <- kdigo_grade[kdigo_grade$patient_id %in% patients_with_preadmit_cr,]})
    try({med_coaga_new <- med_coaga_new[med_coaga_new$patient_id %in% patients_with_preadmit_cr]})
    try({med_coagb_new <- med_coagb_new[med_coagb_new$patient_id %in% patients_with_preadmit_cr]})
    try({med_covid19_new <- med_covid19_new[med_covid19_new$patient_id %in% patients_with_preadmit_cr]})
    try({med_acearb_chronic <- med_acearb_chronic[med_acearb_chronic$patient_id %in% patients_with_preadmit_cr]})
    
    # Set the ratio used for analysis to the ratio calculated from pre-admission Cr only
    peak_trend$ratio <- peak_trend$ratio_prioronly
  }
  
  # If user wishes to customize the Cox PH equations used for recovery and death analysis, we will read in
  # custom files specifying the factors to restrict analyses to.
  # This may be helpful in cases where there may not be enough events for certain factors, causing model
  # convergence to fail.
  # These should be listed as space-separated names in a file "such as the example shown below:
  # age sex ckd cld htn hld ihd
  
  restrict_list <- ""
  model1 <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","ihd","cld")
  model2 <- c("age_group","sex","severe","bronchiectasis","copd","rheum","vte")
  model3 <- c("age_group","sex","severe","COAGA","COAGB","covid_rx")
  model4 <- c("age_group","sex","severe","aki_kdigo_final","ckd","acei_arb_preexposure")
  
  if(restrict_models == TRUE) {
    message("\nWe notice that you are keen to restrict the models to certain variables.")
    message("We are now going to read in the file CustomModelVariables.txt...")
    restrict_list <- scan("Input/CustomModelVariables.txt",what="")
    message(paste("Variables to restrict analyses to :",restrict_list,collapse=" "))
    
  }
  
  message("\n============================\nPart 1: Time to Recovery\n============================")
  # First generate table where we find the time taken to achieve 1.25x baseline ratio
  labs_cr_recovery <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak >= 0) %>% tidyr::fill(severe)
  labs_cr_recovery_tmp <- labs_cr_recovery %>% dplyr::group_by(patient_id) %>% tidyr::complete(time_from_peak = tidyr::full_seq(time_from_peak,1)) %>% dplyr::mutate(ratio = zoo::na.fill(ratio,Inf))
  time_to_ratio1.25 <- labs_cr_recovery_tmp %>% split(.$patient_id) %>% purrr::map(~get_day(.$ratio,.$time_from_peak,target=1.25)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
  colnames(time_to_ratio1.25)[2] <- "time_to_ratio1.25"
  
  labs_aki_summ_index <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission))
  # Get index AKI grade
  index_aki_grade <- labs_aki_summ_index %>% dplyr::select(patient_id,aki_kdigo_final) %>% dplyr::group_by(patient_id) %>% dplyr::filter(aki_kdigo_final == max(aki_kdigo_final)) %>% dplyr::ungroup()
  message("Filtering for AKI patients only...")
  # Filter for AKI cases only
  aki_index_recovery <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::filter(severe %in% c(2,4,5)) %>% dplyr::mutate(severe=ifelse(severe==2,0,1))
  aki_index_recovery <- merge(aki_index_recovery,time_to_ratio1.25,by="patient_id",all.x=TRUE) 
  aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(recover_1.25x = ifelse(is.na(time_to_ratio1.25),0,1))
  message("Computing death and recovery times...")
  # Get death times/censor times
  discharge_day <- demographics %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death_km = dplyr::if_else(deceased==0,as.integer(days_since_admission),as.integer(as.Date(death_date) - as.Date(admission_date)))) %>% dplyr::select(patient_id,deceased,time_to_death_km)
  aki_index_recovery <- merge(aki_index_recovery,discharge_day,by="patient_id",all.x=TRUE)
  #aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_ratio1.25 = dplyr::if_else(recover_1.25x == 0,as.integer(time_to_death_km),as.integer(time_to_ratio1.25)))
  aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_ratio1.25 = dplyr::if_else(recover_1.25x == 0,as.integer(time_to_death_km)-as.integer(peak_cr_time),as.integer(time_to_ratio1.25)))
  
  # Correct the death times for peak Cr
  aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death_km = as.integer(time_to_death_km) - as.integer(peak_cr_time),time_adm_to_death = as.integer(time_to_death_km))
  
  aki_index_recovery <- merge(aki_index_recovery,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE)
  
  message("\nDoing initial filter for medications with more than one factor level.")
  med_recovery_list <- c("COAGA","COAGB","covid_rx","acei_arb_preexposure")
  med_recovery_list <- med_recovery_list[c(coaga_present,coagb_present,dplyr::if_else((covid19antiviral_present == TRUE | remdesivir_present == TRUE),TRUE,FALSE))]
  message("\nAvailable medications: ",paste(med_recovery_list,collapse=" "))
  # First create a temporary table where we filter out the medications with only one factor level
  
  med_recovery_tmp <- aki_index_recovery
  if(isTRUE(coaga_present)) {
    med_recovery_tmp <- merge(med_recovery_tmp,med_coaga_new,by="patient_id",all.x=TRUE)
    med_recovery_tmp <- med_recovery_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(COAGA = dplyr::if_else(is.na(COAGA),0,COAGA)) %>% dplyr::ungroup()
  }
  if(isTRUE(coagb_present)) {
    med_recovery_tmp <- merge(med_recovery_tmp,med_coagb_new,by="patient_id",all.x=TRUE)
    med_recovery_tmp <- med_recovery_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(COAGB=dplyr::if_else(is.na(COAGB),0,COAGB)) %>% dplyr::ungroup()
  }
  if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
    med_recovery_tmp <- merge(med_recovery_tmp,med_covid19_new,by="patient_id",all.x=TRUE)
    med_recovery_tmp <- med_recovery_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = dplyr::if_else(is.na(covid_rx),0,covid_rx)) %>% dplyr::ungroup()
  }
  if(isTRUE(acei_present) | isTRUE(arb_present)) {
    med_recovery_tmp <- merge(med_recovery_tmp,med_acearb_chronic,by="patient_id",all.x=TRUE)
    med_recovery_tmp <- med_recovery_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = dplyr::if_else(is.na(acei_arb_preexposure),0,acei_arb_preexposure)) %>% dplyr::ungroup()
  }
  med_recovery_tmp <- med_recovery_tmp[med_recovery_list]
  med_recovery_tmp <- data.table::as.data.table(lapply(med_recovery_tmp,factor))
  med_recovery_tmp <- data.table::as.data.table(med_recovery_tmp)[,sapply(med_recovery_tmp,function(col) nlevels(col) > 1),with=FALSE]
  med_recovery_list <- colnames(med_recovery_tmp)
  
  # aki_index_recovery <- merge(aki_index_recovery,comorbid,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  # aki_index_recovery[is.na(aki_index_recovery)] <- 0
  # aki_index_recovery[c("severe","aki_kdigo_final",comorbid_list)] <- lapply(aki_index_recovery[c("severe","aki_kdigo_final",comorbid_list)],factor)
  
  message("\nDoing initial filter for comorbids with more than one factor level.")
  # First create a temporary table where we filter out the comorbids with only one factor level
  comorbid_recovery_tmp <- merge(aki_index_recovery,comorbid,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  comorbid_recovery_tmp[is.na(comorbid_recovery_tmp)] <- 0
  comorbid_recovery_tmp <- comorbid_recovery_tmp[comorbid_list]
  comorbid_recovery_tmp <- data.table::as.data.table(lapply(comorbid_recovery_tmp,factor))
  comorbid_recovery_tmp <- data.table::as.data.table(comorbid_recovery_tmp)[,sapply(comorbid_recovery_tmp,function(col) nlevels(col) > 1),with=FALSE] 
  comorbid_recovery_list <- colnames(comorbid_recovery_tmp)
  
  # Then run the actual merging
  aki_index_recovery <- merge(aki_index_recovery,comorbid[c("patient_id",comorbid_recovery_list)],by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  aki_index_recovery <- merge(aki_index_recovery,demog_time_to_event,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  if(isTRUE(coaga_present)) {
    aki_index_recovery <- merge(aki_index_recovery,med_coaga_new,by="patient_id",all.x=TRUE)
  }
  if(isTRUE(coagb_present)) {
    aki_index_recovery <- merge(aki_index_recovery,med_coagb_new,by="patient_id",all.x=TRUE)
  }
  if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
    aki_index_recovery <- merge(aki_index_recovery,med_covid19_new,by="patient_id",all.x=TRUE)
  }
  if(isTRUE(acei_present) | isTRUE(arb_present)) {
    aki_index_recovery <- merge(aki_index_recovery,med_acearb_chronic,by="patient_id",all.x=TRUE)
  }
  aki_index_recovery[is.na(aki_index_recovery)] <- 0
  aki_index_recovery[c("severe","aki_kdigo_final",demog_list,comorbid_recovery_list,med_recovery_list)] <- lapply(aki_index_recovery[c("severe","aki_kdigo_final",demog_list,comorbid_recovery_list,med_recovery_list)],factor)
  
  message("Current factor list for recovery: ",paste(c(demog_list,comorbid_recovery_list,med_recovery_list),collapse=", "))
  
  message("Filtering factor list down further for CoxPH models...")
  # This portion of code deals with the issue of Cox PH models generating large coefficients and/or overfitting
  # We are going to select for the variables where there are at least 5 occurrences of an event for each factor level
  # We will then modify comorbid_recovery_list to only include variable names where this criteria is fulfilled
  # This does NOT require the aki_index_recovery table to be modified
  recovery_tmp <- aki_index_recovery[,c("patient_id","recover_1.25x",demog_list,comorbid_recovery_list,med_recovery_list)] %>% as.data.frame()
  
  comorbid_recovery_list_tmp <- vector(mode="list",length=length(comorbid_recovery_list))
  for(i in 1:length(comorbid_recovery_list)) {
    recovery_tmp1 <- recovery_tmp[,c("patient_id",comorbid_recovery_list[i],"recover_1.25x")]
    recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(comorbid_recovery_list[i]),recover_1.25x)
    recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.25x == 1)
    # if(min(recovery_tmp3$n) >= factor_cutoff) {
    #     comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
    # }
    if(min(recovery_tmp3$n) >= factor_cutoff & nrow(recovery_tmp3) > 1) {
      message(paste0(c("Including ",comorbid_recovery_list[i]," into the comorbid_recovery list...")))
      comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
    }
  }
  comorbid_recovery_list <- unlist(comorbid_recovery_list_tmp[lengths(comorbid_recovery_list_tmp) > 0L])
  
  demog_recovery_list_tmp <- vector(mode="list",length=length(demog_list))
  for(i in 1:length(demog_list)) {
    recovery_tmp1 <- recovery_tmp[,c("patient_id",demog_list[i],"recover_1.25x")]
    recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(demog_list[i]),recover_1.25x)
    recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.25x == 1)
    # if(min(recovery_tmp3$n) >= factor_cutoff) {
    #     comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
    # }
    if(min(recovery_tmp3$n) >= factor_cutoff & nrow(recovery_tmp3) > 1) {
      message(paste0(c("Including ",demog_list[i]," into the demog_recovery list...")))
      demog_recovery_list_tmp[i] <- demog_list[i]
    }
  }
  demog_recovery_list <- unlist(demog_recovery_list_tmp[lengths(demog_recovery_list_tmp) > 0L])
  
  med_recovery_list_tmp <- vector(mode="list",length=length(med_recovery_list))
  if(length(med_recovery_list)>0) {
    for(i in 1:length(med_recovery_list)) {
      recovery_tmp1 <- recovery_tmp[,c("patient_id",med_recovery_list[i],"recover_1.25x")]
      recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(med_recovery_list[i]),recover_1.25x)
      recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.25x == 1)
      # if(min(recovery_tmp3$n) >= factor_cutoff) {
      #     comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
      # }
      if(min(recovery_tmp3$n) >= factor_cutoff & nrow(recovery_tmp3) > 1) {
        message(paste0(c("Including ",med_recovery_list[i]," into the med_recovery list...")))
        med_recovery_list_tmp[i] <- med_recovery_list[i]
      }
    }
  }
  med_recovery_list <- unlist(med_recovery_list_tmp[lengths(med_recovery_list_tmp) > 0L])
  
  message("\nFinal factor list for recovery (before user customisation): ",paste(c(demog_recovery_list,comorbid_recovery_list,med_recovery_list),collapse=" "))
  
  if(restrict_models == TRUE) {
    demog_recovery_list <- demog_recovery_list[demog_recovery_list %in% restrict_list]
    comorbid_recovery_list <- comorbid_recovery_list[comorbid_recovery_list %in% restrict_list]
    med_recovery_list <- med_recovery_list[med_recovery_list %in% restrict_list]
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_recovery_list,"\nComorbidities:",comorbid_recovery_list,"\nMedications:",med_recovery_list,sep = " "))
  }
  variable_list_output <- paste(c("Final Recovery variable list:",demog_recovery_list,comorbid_recovery_list,med_recovery_list),collapse=" ")
  readr::write_lines(variable_list_output,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_custom_equation.txt")),append=F)
  
  message("Now proceeding to time-to-Cr recovery analysis...")
  # Now run the actual time-to-event analysis
  
  # Kaplan Meier plot for COVID-19 severity
  message("Generating Kaplan-Meier plots...")
  recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe")
  fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
  plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
  plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
  plot_recover_summ_table <- plot_recover$data.survtable
  write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_PlotSummStats.csv")),row.names=TRUE)
  write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Plot.csv")),row.names=FALSE)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe.png")),plot=print(plot_recover),width=12,height=12,units="cm")
  
  # Kaplan Meier plot for KDIGO grades
  recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ aki_kdigo_final")
  fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
  plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
  plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
  write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO_PlotSummStats.csv")),row.names=TRUE)
  write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO_Plot.csv")),row.names=FALSE)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO.png")),plot=print(plot_recover),width=12,height=12,units="cm")
  
  # # Kaplan Meier plot for KDIGO grades - collapsing KDIGO2/3 into one group
  # # This also generates the same collapsed table for Model 1B later
  # aki_index_recovery_collapse <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki_kdigo_final = dplyr::if_else(as.numeric(aki_kdigo_final) >=2,2,1)) %>% dplyr::ungroup()
  # aki_index_recovery_collapse$aki_kdigo_final <- as.factor(aki_index_recovery_collapse$aki_kdigo_final)
  # fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery_collapse)
  # plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery_collapse,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
  # plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery_collapse)
  # write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3_PlotSummStats.csv")),row.names=TRUE)
  # write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3_Plot.csv")),row.names=FALSE)
  # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3.png")),plot=print(plot_recover),width=12,height=12,units="cm")
  # 
  
  if(isTRUE(ckd_present)) {
    try({
      # Kaplan Meier plot for CKD
      recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ckd")
      fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
      plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
      plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
      write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CKD_PlotSummStats.csv")),row.names=TRUE)
      write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CKD_Plot.csv")),row.names=FALSE)
      ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CKD.png")),plot=print(plot_recover),width=12,height=12,units="cm")
    })
  }
  
  # CoxPH model
  # Generate univariate analyses first
  message("Generating univariate Cox PH models (time to recovery, AKI patients only)...")
  univ_formulas <- sapply(c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list), function(x) as.formula(paste('survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ', x)))
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_recovery)}),
    error = function(c) "Problem generating univariate models."
  )
  # Extract data 
  
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_recover <- do.call("rbind",univ_results)
    write.csv(univ_results_recover,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  message("\nGenerating Model 1 (time to recovery, AKI patients only)...")
  try({
    recovery_model1 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model1 <- recovery_model1[recovery_model1 %in% model1]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model1,collapse="+")))
    message(paste("Formula for Model 1: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model1,collapse="+")))
    coxph_recover1 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover1_summ <- summary(coxph_recover1) 
    print(coxph_recover1_summ)
    coxph_recover1_hr <- cbind(coxph_recover1_summ$coefficients,coxph_recover1_summ$conf.int)[,-c(6,7)]
    coxph_recover1_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover1_summ$logtest,coxph_recover1_summ$sctest,coxph_recover1_summ$waldtest))
    coxph_recover1_stats2 <- rbind(data.table::as.data.table(coxph_recover1_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover1_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover1_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model1.csv")),row.names=TRUE)
    write.csv(coxph_recover1_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover1_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover1_plot <- survminer::ggforest(coxph_recover1,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model1.png")),plot=print(coxph_recover1_plot),width=20,height=20,units="cm")
  })
  
  # message("\nGenerating Model 1B - collapsing KDIGO 2/3 to single group (time to recovery, AKI patients only)...")
  # try({
  #   recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model1,collapse="+")))
  #   message(paste("Formula for Model 1: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model1,collapse="+")))
  #   coxph_recovery1b <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery_collapse)
  #   coxph_recovery1b_summ <- summary(coxph_recovery1b) 
  #   print(coxph_recovery1b_summ)
  #   coxph_recovery1b_hr <- cbind(coxph_recovery1b_summ$coefficients,coxph_recovery1b_summ$conf.int)[,-c(6,7)]
  #   coxph_recovery1b_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recovery1b_summ$logtest,coxph_recovery1b_summ$sctest,coxph_recovery1b_summ$waldtest))
  #   coxph_recovery1b_stats2 <- rbind(data.table::as.data.table(coxph_recovery1b_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recovery1b_summ$rsq,keep.rownames = T))
  #   write.csv(coxph_recovery1b_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3_CoxPH_Model1.csv")),row.names=TRUE)
  #   write.csv(coxph_recovery1b_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
  #   write.csv(coxph_recovery1b_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
  #   coxph_recovery1b_plot <- survminer::ggforest(coxph_recovery1b,data=aki_index_recovery_collapse)
  #   ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_KDIGO1vs2+3_CoxPH_Model1.png")),plot=print(coxph_recovery1b_plot),width=20,height=20,units="cm")
  # })
  
  
  message("\nGenerating Model 2 (time to recovery, AKI patients only)...")
  try({
    recovery_model2 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model2 <- recovery_model2[recovery_model2 %in% model2]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model2,collapse="+")))
    message(paste("Formula for Model 2: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model2,collapse="+")))
    coxph_recover2 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover2_summ <- summary(coxph_recover2) 
    print(coxph_recover2_summ)
    coxph_recover2_hr <- cbind(coxph_recover2_summ$coefficients,coxph_recover2_summ$conf.int)[,-c(6,7)]
    coxph_recover2_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover2_summ$logtest,coxph_recover2_summ$sctest,coxph_recover2_summ$waldtest))
    coxph_recover2_stats2 <- rbind(data.table::as.data.table(coxph_recover2_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover2_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover2_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model2.csv")),row.names=TRUE)
    write.csv(coxph_recover2_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model2_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover2_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model2_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover2_plot <- survminer::ggforest(coxph_recover2,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model2.png")),plot=print(coxph_recover2_plot),width=20,height=20,units="cm")
  })
  
  message("\nGenerating Model 3 (time to recovery, AKI patients only)...")
  try({
    recovery_model3 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model3 <- recovery_model3[recovery_model3 %in% model3]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
    message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
    # if(("covid_rx" %in% recovery_model3) == TRUE) {
    #     recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model3,"severe:covid_rx"),collapse="+")))
    #     message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model3,"severe:covid_rx"),collapse="+")))
    # } else {
    #     recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
    #     message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
    # }
    message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
    coxph_recover3 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover3_summ <- summary(coxph_recover3) 
    print(coxph_recover3_summ)
    coxph_recover3_hr <- cbind(coxph_recover3_summ$coefficients,coxph_recover3_summ$conf.int)[,-c(6,7)]
    coxph_recover3_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover3_summ$logtest,coxph_recover3_summ$sctest,coxph_recover3_summ$waldtest))
    coxph_recover3_stats2 <- rbind(data.table::as.data.table(coxph_recover3_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover3_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover3_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model3.csv")),row.names=TRUE)
    write.csv(coxph_recover3_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model3_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover3_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model3_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover3_plot <- survminer::ggforest(coxph_recover3,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model3.png")),plot=print(coxph_recover3_plot),width=20,height=20,units="cm")
  })
  message("\nGenerating Model 4 with ACE-i/ARBs (time to recovery, AKI patients only)...")
  try({
    recovery_model4 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model4 <- recovery_model4[recovery_model4 %in% model4]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model4,collapse="+")))
    message(paste("Formula for Model 4: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model4,collapse="+")))
    coxph_recover4 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover4_summ <- summary(coxph_recover4) 
    print(coxph_recover4_summ)
    coxph_recover4_hr <- cbind(coxph_recover4_summ$coefficients,coxph_recover4_summ$conf.int)[,-c(6,7)]
    coxph_recover4_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover4_summ$logtest,coxph_recover4_summ$sctest,coxph_recover4_summ$waldtest))
    coxph_recover4_stats2 <- rbind(data.table::as.data.table(coxph_recover4_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover4_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover4_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model4.csv")),row.names=TRUE)
    write.csv(coxph_recover4_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model4_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover4_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model4_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover4_plot <- survminer::ggforest(coxph_recover4,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH_Model4.png")),plot=print(coxph_recover4_plot),width=20,height=20,units="cm")
  })
  
  message("If you are getting any errors with model generation - do note that it may actually be normal to get errors\nif your site numbers are low (especially for model 3). Please check your data to see if the appropriate\nnumber of events occur for each factor level.")
  
  message("Now proceeding to time-to-death analysis for AKI patients only...")
  
  deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
  fit_death_aki_only <- survminer::surv_fit(deathPlotFormula, data=aki_index_recovery)
  plot_death_aki_only <- survminer::ggsurvplot(fit_death_aki_only,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
  plot_death_aki_only_summ <- survminer::surv_summary(fit_death_aki_only,data=aki_index_recovery)
  plot_death_aki_only_summ_table <- plot_death_aki_only$data.survtable
  write.csv(fit_death_aki_only$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_PlotSummStats.csv")),row.names=TRUE)
  write.csv(plot_death_aki_only_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
  write.csv(plot_death_aki_only_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Table.csv")),row.names=FALSE)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe.png")),plot=print(plot_death_aki_only),width=12,height=12,units="cm")
  
  deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ aki_kdigo_final")
  fit_death_aki_only <- survminer::surv_fit(deathPlotFormula, data=aki_index_recovery)
  plot_death_aki_only <- survminer::ggsurvplot(fit_death_aki_only,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
  plot_death_aki_only_summ <- survminer::surv_summary(fit_death_aki_only,data=aki_index_recovery)
  write.csv(fit_death_aki_only$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO_PlotSummStats.csv")),row.names=TRUE)
  write.csv(plot_death_aki_only_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO_Plot.csv")),row.names=FALSE)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO.png")),plot=print(plot_death_aki_only),width=12,height=12,units="cm")
  
  if(isTRUE(ckd_present)) {
    try({
      deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ ckd")
      fit_death_aki_only <- survminer::surv_fit(deathPlotFormula, data=aki_index_recovery)
      plot_death_aki_only <- survminer::ggsurvplot(fit_death_aki_only,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
      plot_death_aki_only_summ <- survminer::surv_summary(fit_death_aki_only,data=aki_index_recovery)
      write.csv(fit_death_aki_only$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CKD_PlotSummStats.csv")),row.names=TRUE)
      write.csv(plot_death_aki_only_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CKD_Plot.csv")),row.names=FALSE)
      ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CKD.png")),plot=print(plot_death_aki_only),width=12,height=12,units="cm")
      
    })
  }
  
  message("Generating univariate Cox PH models (Time to death, AKI patients only)...")
  univ_formulas <- sapply(c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_recovery)}),
    error = "Problem generating univariate models."
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_death_akionly <- do.call("rbind",univ_results)
    write.csv(univ_results_death_akionly,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  message("Generating Model 1 (Time to death, AKI patients only)...")
  try({
    death_aki_only_model1 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    death_aki_only_model1 <- death_aki_only_model1[death_aki_only_model1 %in% model1]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model1,collapse="+")))
    message("Formula for Model 1: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model1,collapse="+")))
    coxph_death_akionly1 <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery)
    coxph_death_akionly1_summ <- summary(coxph_death_akionly1)
    print(coxph_death_akionly1_summ)
    coxph_death_akionly1_hr <- cbind(coxph_death_akionly1_summ$coefficients,coxph_death_akionly1_summ$conf.int)[,-c(6,7)]
    coxph_death_akionly1_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_akionly1_summ$logtest,coxph_death_akionly1_summ$sctest,coxph_death_akionly1_summ$waldtest))
    coxph_death_akionly1_stats2 <- rbind(data.table::as.data.table(coxph_death_akionly1_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_akionly1_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_akionly1_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model1.csv")),row.names=TRUE)
    write.csv(coxph_death_akionly1_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_akionly1_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_akionly1_plot <- survminer::ggforest(coxph_death_akionly1,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model1.png")),plot=print(coxph_death_akionly1_plot),width=20,height=20,units="cm")
  })
  # message("Generating Model 1b - collapsing KDIGO2/3 into single group (Time to death, AKI patients only)...")
  # try({
  #   coxph_death_akionly1b <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery_collapse)
  #   coxph_death_akionly1b_summ <- summary(coxph_death_akionly1b)
  #   print(coxph_death_akionly1b_summ)
  #   coxph_death_akionly1b_hr <- cbind(coxph_death_akionly1b_summ$coefficients,coxph_death_akionly1b_summ$conf.int)[,-c(6,7)]
  #   coxph_death_akionly1b_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_akionly1b_summ$logtest,coxph_death_akionly1b_summ$sctest,coxph_death_akionly1b_summ$waldtest))
  #   coxph_death_akionly1b_stats2 <- rbind(data.table::as.data.table(coxph_death_akionly1b_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_akionly1b_summ$rsq,keep.rownames = T))
  #   write.csv(coxph_death_akionly1b_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO1vs2+3_CoxPH_Model1.csv")),row.names=TRUE)
  #   write.csv(coxph_death_akionly1b_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO1vs2+3_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
  #   write.csv(coxph_death_akionly1b_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO1vs2+3_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
  #   coxph_death_akionly1b_plot <- survminer::ggforest(coxph_death_akionly1b,data=aki_index_recovery)
  #   ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_KDIGO1vs2+3_CoxPH_Model1.png")),plot=print(coxph_death_akionly1b_plot),width=20,height=20,units="cm")
  # })
  
  message("Generating Model 2 (Time to death, AKI patients only)...")
  try({
    death_aki_only_model2 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    death_aki_only_model2 <- death_aki_only_model2[death_aki_only_model2 %in% model2]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model2,collapse="+")))
    message("Formula for Model 2: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model2,collapse="+")))
    coxph_death_akionly2 <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery)
    coxph_death_akionly2_summ <- summary(coxph_death_akionly2) 
    print(coxph_death_akionly2_summ)
    coxph_death_akionly2_hr <- cbind(coxph_death_akionly2_summ$coefficients,coxph_death_akionly2_summ$conf.int)[,-c(6,7)]
    coxph_death_akionly2_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_akionly2_summ$logtest,coxph_death_akionly2_summ$sctest,coxph_death_akionly2_summ$waldtest))
    coxph_death_akionly2_stats2 <- rbind(data.table::as.data.table(coxph_death_akionly2_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_akionly2_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_akionly2_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model2.csv")),row.names=TRUE)
    write.csv(coxph_death_akionly2_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model2_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_akionly2_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model2_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_akionly2_plot <- survminer::ggforest(coxph_death_akionly2,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model2.png")),plot=print(coxph_death_akionly2_plot),width=20,height=20,units="cm")
  })
  
  message("Generating Model 3 (Time to death, AKI patients only)...")
  try({
    death_aki_only_model3 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    death_aki_only_model3 <- death_aki_only_model3[death_aki_only_model3 %in% model3]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model3,collapse="+")))
    message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model3,collapse="+")))
    # if(("covid_rx" %in% death_aki_only_model3) == TRUE) {
    #     deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_aki_only_model3,"severe * covid_rx"),collapse="+")))
    #     message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_aki_only_model3,"severe * covid_rx"),death_aki_only_model3,collapse="+")))
    # } else {
    #     deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model3,collapse="+")))
    #     message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model3,collapse="+")))
    # }
    coxph_death_akionly3 <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery)
    coxph_death_akionly3_summ <- summary(coxph_death_akionly3) 
    print(coxph_death_akionly3_summ)
    coxph_death_akionly3_hr <- cbind(coxph_death_akionly3_summ$coefficients,coxph_death_akionly3_summ$conf.int)[,-c(6,7)]
    coxph_death_akionly3_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_akionly3_summ$logtest,coxph_death_akionly3_summ$sctest,coxph_death_akionly3_summ$waldtest))
    coxph_death_akionly3_stats2 <- rbind(data.table::as.data.table(coxph_death_akionly3_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_akionly3_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_akionly3_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model3.csv")),row.names=TRUE)
    write.csv(coxph_death_akionly3_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model3_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_akionly3_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model3_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_akionly3_plot <- survminer::ggforest(coxph_death_akionly3,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model3.png")),plot=print(coxph_death_akionly3_plot),width=20,height=20,units="cm")
  })
  
  message("Generating Model 4 (Time to death, AKI patients only)...")
  try({
    death_aki_only_model4 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    death_aki_only_model4 <- death_aki_only_model4[death_aki_only_model4 %in% model4]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model4,collapse="+")))
    message("Formula for Model 4: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_aki_only_model4,collapse="+")))
    coxph_death_akionly4 <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery)
    coxph_death_akionly4_summ <- summary(coxph_death_akionly4)
    print(coxph_death_akionly4_summ)
    coxph_death_akionly4_hr <- cbind(coxph_death_akionly4_summ$coefficients,coxph_death_akionly4_summ$conf.int)[,-c(6,7)]
    coxph_death_akionly4_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_akionly4_summ$logtest,coxph_death_akionly4_summ$sctest,coxph_death_akionly4_summ$waldtest))
    coxph_death_akionly4_stats2 <- rbind(data.table::as.data.table(coxph_death_akionly4_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_akionly4_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_akionly4_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model4.csv")),row.names=TRUE)
    write.csv(coxph_death_akionly4_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model4_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_akionly4_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model4_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_akionly4_plot <- survminer::ggforest(coxph_death_akionly4,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH_Model4.png")),plot=print(coxph_death_akionly4_plot),width=20,height=20,units="cm")
  })
  
  message("If you are getting any errors with model generation - do note that it may actually be normal to get errors\nif your site numbers are low (especially for model 3). Please check your data to see if the appropriate\nnumber of events occur for each factor level.")
  
  # We now do the same to the time to death analyses:
  # 1) Create a temporary comorbid table and obtain the valid comorbids with more than 1 factor level
  
  message("Part 2: Time to Death analysis...")
  aki_index_death <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(is_aki=ifelse(severe %in% c(2,4,5),1,0)) %>% dplyr::mutate(severe=ifelse(severe %in% c(3,4,5),1,0))
  aki_index_death <- merge(aki_index_death,discharge_day,by="patient_id",all.x=TRUE) # merge in time_to_death_km
  # Correct death times for peak sCr time
  aki_index_death <- aki_index_death %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death_km = as.integer(time_to_death_km) - as.integer(peak_cr_time), time_adm_to_death = as.integer(time_to_death_km)) %>% dplyr::ungroup()
  aki_index_death <- merge(aki_index_death,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE) # AKI KDIGO grade
  
  message("\nDoing initial filter for medications with more than one factor level.")
  med_death_list <- c("COAGA","COAGB","covid_rx")
  med_death_list <- med_death_list[c(coaga_present,coagb_present,dplyr::if_else((covid19antiviral_present == TRUE | remdesivir_present == TRUE),TRUE,FALSE))]
  message("\nAvailable medications: ",paste(med_death_list,collapse=" "))
  # First create a temporary table where we filter out the medications with only one factor level
  
  med_death_tmp <- aki_index_death
  if(coaga_present == TRUE) {
    med_death_tmp <- merge(med_death_tmp,med_coaga_new,by="patient_id",all.x=TRUE)
    med_death_tmp <- med_death_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(COAGA = dplyr::if_else(is.na(COAGA),0,COAGA)) %>% dplyr::ungroup()
  }
  if(coagb_present == TRUE) {
    med_death_tmp <- merge(med_death_tmp,med_coagb_new,by="patient_id",all.x=TRUE)
    med_death_tmp <- med_death_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(COAGB=dplyr::if_else(is.na(COAGB),0,COAGB)) %>% dplyr::ungroup()
  }
  if(covid19antiviral_present == TRUE | remdesivir_present == TRUE) {
    med_death_tmp <- merge(med_death_tmp,med_covid19_new,by="patient_id",all.x=TRUE)
    med_death_tmp <- med_death_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = dplyr::if_else(is.na(covid_rx),0,covid_rx)) %>% dplyr::ungroup()
  }
  
  if(isTRUE(acei_present) | isTRUE(arb_present)) {
    med_death_tmp <- merge(med_death_tmp,med_acearb_chronic,by="patient_id",all.x=TRUE)
    med_death_tmp <- med_death_tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = dplyr::if_else(is.na(acei_arb_preexposure),0,acei_arb_preexposure)) %>% dplyr::ungroup()
  }
  
  med_death_tmp <- med_death_tmp[med_death_list]
  med_death_tmp <- data.table::as.data.table(lapply(med_death_tmp,factor))
  med_death_tmp <- data.table::as.data.table(med_death_tmp)[,sapply(med_death_tmp,function(col) nlevels(col) > 1),with=FALSE]
  med_death_list <- colnames(med_death_tmp)
  
  comorbid_death_tmp <- merge(aki_index_death,comorbid,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  comorbid_death_tmp[is.na(comorbid_death_tmp)] <- 0
  comorbid_death_tmp <- comorbid_death_tmp[comorbid_list]
  comorbid_death_tmp <- data.table::as.data.table(lapply(comorbid_death_tmp,factor))
  comorbid_death_tmp <- data.table::as.data.table(comorbid_death_tmp)[,sapply(comorbid_death_tmp,function(col) nlevels(col) > 1),with=FALSE] 
  comorbid_death_list <- colnames(comorbid_death_tmp)
  
  message("Factor list for Death Analysis before filtering for CoxPH: ",paste(c(demog_list,comorbid_death_list,med_death_list),collapse = " "))
  # 2) Create a new table with the cleaned up comorbids
  aki_index_death <- merge(aki_index_death,comorbid[c("patient_id",comorbid_death_list)],by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  aki_index_death <- merge(aki_index_death,demog_time_to_event,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
  if(coaga_present == TRUE) {
    aki_index_death <- merge(aki_index_death,med_coaga_new,by="patient_id",all.x=TRUE)
  }
  if(coagb_present == TRUE) {
    aki_index_death <- merge(aki_index_death,med_coagb_new,by="patient_id",all.x=TRUE)
  }
  if(covid19antiviral_present == TRUE | remdesivir_present == TRUE) {
    aki_index_death <- merge(aki_index_death,med_covid19_new,by="patient_id",all.x=TRUE)
  }
  if(isTRUE(acei_present) | isTRUE(arb_present)) {
    aki_index_death <- merge(aki_index_death,med_acearb_chronic,by="patient_id",all.x=TRUE)
  }
  
  aki_index_death[is.na(aki_index_death)] <- 0
  aki_index_death <- aki_index_death %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = dplyr::if_else(!is.na(severe_to_aki),as.integer(min(severe_to_aki)),NA_integer_)) %>% dplyr::distinct()
  aki_index_death[c("severe","aki_kdigo_final","is_aki",demog_list,comorbid_death_list,med_death_list)] <- lapply(aki_index_death[c("severe","aki_kdigo_final","is_aki",demog_list,comorbid_death_list,med_death_list)],factor)
  
  # 3) This portion of code deals with the issue of Cox PH models generating large coefficients and/or overfitting
  # We are going to select for the variables where there are at least 5 occurrences of an event for each factor level
  # We will then modify comorbid_death_list to only include variable names where this criteria is fulfilled
  # This does NOT require the aki_index_death table to be modified
  death_tmp <- aki_index_death[,c("patient_id","deceased",demog_list,comorbid_death_list,med_death_list)] %>% as.data.frame()
  comorbid_death_list_tmp <- vector(mode="list",length=length(comorbid_death_list))
  for(i in 1:length(comorbid_death_list)) {
    death_tmp1 <- death_tmp[,c("patient_id",comorbid_death_list[i],"deceased")]
    death_tmp2 <- death_tmp1 %>% dplyr::count(get(comorbid_death_list[i]),deceased)
    death_tmp3 <- death_tmp2 %>% dplyr::filter(deceased == 1)
    if(min(death_tmp3$n) >= factor_cutoff & nrow(death_tmp3) > 1) {
      message(paste0(c("Including ",comorbid_death_list[i]," into the comorbid_death list...")))
      comorbid_death_list_tmp[i] <- comorbid_death_list[i]
    }
  }
  comorbid_death_list <- unlist(comorbid_death_list_tmp[lengths(comorbid_death_list_tmp) > 0L])
  
  demog_death_list_tmp <- vector(mode="list",length=length(demog_list))
  for(i in 1:length(demog_list)) {
    death_tmp1 <- death_tmp[,c("patient_id",demog_list[i],"deceased")]
    death_tmp2 <- death_tmp1 %>% dplyr::count(get(demog_list[i]),deceased)
    death_tmp3 <- death_tmp2 %>% dplyr::filter(deceased == 1)
    # if(min(death_tmp3$n) >= factor_cutoff) {
    #     comorbid_death_list_tmp[i] <- comorbid_death_list[i]
    # }
    if(min(death_tmp3$n) >= factor_cutoff & nrow(death_tmp3) > 1) {
      message(paste0(c("Including ",demog_list[i]," into the demog_death list...")))
      demog_death_list_tmp[i] <- demog_list[i]
    }
  }
  demog_death_list <- unlist(demog_death_list_tmp[lengths(demog_death_list_tmp) > 0L])
  
  med_death_list_tmp <- vector(mode="list",length=length(med_death_list))
  if(length(med_death_list)>0) {
    for(i in 1:length(med_death_list)) {
      death_tmp1 <- death_tmp[,c("patient_id",med_death_list[i],"deceased")]
      death_tmp2 <- death_tmp1 %>% dplyr::count(get(med_death_list[i]),deceased)
      death_tmp3 <- death_tmp2 %>% dplyr::filter(deceased == 1)
      # if(min(death_tmp3$n) >= factor_cutoff) {
      #     comorbid_death_list_tmp[i] <- comorbid_death_list[i]
      # }
      if(min(death_tmp3$n) >= factor_cutoff & nrow(death_tmp3) > 1) {
        message(paste0(c("Including ",med_death_list[i]," into the med_death list...")))
        med_death_list_tmp[i] <- med_death_list[i]
      }
    }
  }
  med_death_list <- unlist(med_death_list_tmp[lengths(med_death_list_tmp) > 0L])
  
  message("\nFinal factor list for death (before user customisation):",paste(c(demog_death_list,comorbid_death_list,med_death_list),collapse=" "))
  
  if(restrict_models == TRUE) {
    demog_death_list <- demog_death_list[demog_death_list %in% restrict_list]
    comorbid_death_list <- comorbid_death_list[comorbid_death_list %in% restrict_list]
    med_death_list <- med_death_list[med_death_list %in% restrict_list]
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_death_list,"\nComorbidities:",comorbid_death_list,"\nMedications:",med_death_list,sep = " "))
  }
  variable_list_death <- paste(c("Final Death variable list: ",demog_death_list,comorbid_death_list,med_death_list),collapse=" ")
  # readr::write_lines(variable_list_death,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_custom_equation.txt")),append=T)
  
  # 4) Run analysis
  message("Now proceeding with time-to-event analysis...")
  message("Generating Kaplan-Meier curves (time to death, all patients)...")
  deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ is_aki")
  fit_death <- survminer::surv_fit(deathPlotFormula, data=aki_index_death)
  plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
  plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death)
  plot_death_summ_table <- plot_death$data.survtable
  write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_Plot.csv")),row.names=FALSE)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI.png")),plot=print(plot_death),width=12,height=12,units="cm")
  
  # Survival curves stratified by KDIGO stage
  deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ aki_kdigo_final")
  fit_death <- survminer::surv_fit(deathPlotFormula, data=aki_index_death)
  plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
  plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death)
  plot_death_summ_table <- plot_death$data.survtable
  write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_KDIGO_Plot.csv")),row.names=FALSE)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_KDIGO.png")),plot=print(plot_death),width=12,height=12,units="cm")
  
  # # Collapse KDIGO2/3 into single group
  # aki_index_death_collapse <- aki_index_death %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki_kdigo_final = dplyr::if_else(as.numeric(aki_kdigo_final) >=2,2,as.numeric(aki_kdigo_final))) %>% dplyr::ungroup()
  # aki_index_death_collapse$aki_kdigo_final <- as.factor(aki_index_death_collapse$aki_kdigo_final)                                                                           
  # fit_death <- survminer::surv_fit(deathPlotFormula, data=aki_index_death_collapse)
  # plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death_collapse,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
  # plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death_collapse)
  # plot_death_summ_table <- plot_death$data.survtable
  # write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_KDIGO1vs2+3_Plot.csv")),row.names=FALSE)
  # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_KDIGO1vs2+3.png")),plot=print(plot_death),width=12,height=12,units="cm")
  # 
  
  if(isTRUE(ckd_present)) {
    try({
      deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ ckd")
      fit_death <- survminer::surv_fit(deathPlotFormula, data=aki_index_death)
      plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
      plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death)
      plot_death_summ_table <- plot_death$data.survtable
      write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CKD_Plot.csv")),row.names=FALSE)
      ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CKD.png")),plot=print(plot_death),width=12,height=12,units="cm")
      
    })
  }
  
  message("Generating univariate Cox PH models (Time to death, all patients)...")
  univ_formulas <- sapply(c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_death)}),
    error = "Problem generating univariate models."
  )
  # Extract data 
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_death_all <- do.call("rbind",univ_results)
    write.csv(univ_results_death_all,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Univariate.csv")),row.names=TRUE)
    # save(univ_results_death_all,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Univariate.rdata")))
  })
  
  message("Generating Model 1 (Time to death, all patients)...")
  try({
    death_model1 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
    death_model1 <- death_model1[death_model1 %in% model1]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model1,collapse="+")))
    message("Formula for Model 1: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model1,collapse="+")))
    coxph_death_all1 <- survival::coxph(deathCoxPHFormula, data=aki_index_death)
    coxph_death_all1_summ <- summary(coxph_death_all1)
    print(coxph_death_all1_summ)
    coxph_death_all1_hr <- cbind(coxph_death_all1_summ$coefficients,coxph_death_all1_summ$conf.int)[,-c(6,7)]
    coxph_death_all1_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_all1_summ$logtest,coxph_death_all1_summ$sctest,coxph_death_all1_summ$waldtest))
    coxph_death_all1_stats2 <- rbind(data.table::as.data.table(coxph_death_all1_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_all1_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_all1_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model1.csv")),row.names=TRUE)
    write.csv(coxph_death_all1_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_all1_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_all1_plot <- survminer::ggforest(coxph_death_all1,data=aki_index_death)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model1.png")),plot=print(coxph_death_all1_plot),width=20,height=20,units="cm")
  })
  # 
  # message("Generating Model 1b - collapsing KDIGO2/3 into single group (Time to death, all patients)...")
  # try({
  #   death_model1 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
  #   death_model1 <- death_model1[death_model1 %in% model1]
  #   deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model1,collapse="+")))
  #   message("Formula for Model 1: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model1,collapse="+")))
  #   coxph_death_all1b <- survival::coxph(deathCoxPHFormula, data=aki_index_death_collapse)
  #   coxph_death_all1b_summ <- summary(coxph_death_all1b)
  #   print(coxph_death_all1b_summ)
  #   coxph_death_all1b_hr <- cbind(coxph_death_all1b_summ$coefficients,coxph_death_all1b_summ$conf.int)[,-c(6,7)]
  #   coxph_death_all1b_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_all1b_summ$logtest,coxph_death_all1b_summ$sctest,coxph_death_all1b_summ$waldtest))
  #   coxph_death_all1b_stats2 <- rbind(data.table::as.data.table(coxph_death_all1b_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_all1b_summ$rsq,keep.rownames = T))
  #   write.csv(coxph_death_all1b_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_KDIGO1vs2+3_CoxPH_Model1.csv")),row.names=TRUE)
  #   write.csv(coxph_death_all1b_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_KDIGO1vs2+3_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
  #   write.csv(coxph_death_all1b_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_KDIGO1vs2+3_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
  #   coxph_death_all1b_plot <- survminer::ggforest(coxph_death_all1b,data=aki_index_death)
  #   ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_KDIGO1vs2+3_CoxPH_Model1.png")),plot=print(coxph_death_all1b_plot),width=20,height=20,units="cm")
  # })
  
  message("Generating Model 2 (Time to death, all patients)...")
  try({
    death_model2 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
    death_model2 <- death_model2[death_model2 %in% model2]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model2,collapse="+")))
    message("Formula for Model 2: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model2,collapse="+")))
    coxph_death_all2 <- survival::coxph(deathCoxPHFormula, data=aki_index_death)
    coxph_death_all2_summ <- summary(coxph_death_all2) 
    print(coxph_death_all2_summ)
    coxph_death_all2_hr <- cbind(coxph_death_all2_summ$coefficients,coxph_death_all2_summ$conf.int)[,-c(6,7)]
    coxph_death_all2_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_all2_summ$logtest,coxph_death_all2_summ$sctest,coxph_death_all2_summ$waldtest))
    coxph_death_all2_stats2 <- rbind(data.table::as.data.table(coxph_death_all2_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_all2_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_all2_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model2.csv")),row.names=TRUE)
    write.csv(coxph_death_all2_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model2_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_all2_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model2_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_all2_plot <- survminer::ggforest(coxph_death_all2,data=aki_index_death)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model2.png")),plot=print(coxph_death_all2_plot),width=20,height=20,units="cm")
  })
  
  message("Generating Model 3 (Time to death, all patients)...")
  try({
    death_model3 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
    death_model3 <- death_model3[death_model3 %in% model3]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model3,collapse="+")))
    message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model3,collapse="+")))
    coxph_death_all3 <- survival::coxph(deathCoxPHFormula, data=aki_index_death)
    coxph_death_all3_summ <- summary(coxph_death_all3) 
    print(coxph_death_all3_summ)
    coxph_death_all3_hr <- cbind(coxph_death_all3_summ$coefficients,coxph_death_all3_summ$conf.int)[,-c(6,7)]
    coxph_death_all3_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_all3_summ$logtest,coxph_death_all3_summ$sctest,coxph_death_all3_summ$waldtest))
    coxph_death_all3_stats2 <- rbind(data.table::as.data.table(coxph_death_all3_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_all3_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_all3_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model3.csv")),row.names=TRUE)
    write.csv(coxph_death_all3_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model3_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_all3_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model3_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_all3_plot <- survminer::ggforest(coxph_death_all3,data=aki_index_death)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model3.png")),plot=print(coxph_death_all3_plot),width=20,height=20,units="cm")
  })
  
  message("Generating Model 4 (Time to death, all patients)...")
  try({
    death_model4 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
    death_model4 <- death_model4[death_model4 %in% model4]
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model4,collapse="+")))
    message("Formula for Model 4: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model4,collapse="+")))
    coxph_death_all4 <- survival::coxph(deathCoxPHFormula, data=aki_index_death)
    coxph_death_all4_summ <- summary(coxph_death_all4)
    print(coxph_death_all4_summ)
    coxph_death_all4_hr <- cbind(coxph_death_all4_summ$coefficients,coxph_death_all4_summ$conf.int)[,-c(6,7)]
    coxph_death_all4_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_death_all4_summ$logtest,coxph_death_all4_summ$sctest,coxph_death_all4_summ$waldtest))
    coxph_death_all4_stats2 <- rbind(data.table::as.data.table(coxph_death_all4_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_death_all4_summ$rsq,keep.rownames = T))
    write.csv(coxph_death_all4_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model4.csv")),row.names=TRUE)
    write.csv(coxph_death_all4_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model4_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_death_all4_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model4_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_death_all4_plot <- survminer::ggforest(coxph_death_all4,data=aki_index_death)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model4.png")),plot=print(coxph_death_all4_plot),width=20,height=20,units="cm")
  })
  
  message("If you are getting any errors with model generation - do note that it may actually be normal to get errors\nif your site numbers are low (especially for model 3). Please check your data to see if the appropriate\nnumber of events occur for each factor level.")
}