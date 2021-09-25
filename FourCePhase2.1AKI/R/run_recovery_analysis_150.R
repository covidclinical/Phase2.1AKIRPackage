#' Runs time-to-event analysis for the 4 pre-defined models - Kaplan-Meier curves, Cox PH models.
#' This version uses the 150% baseline sCr as the definition for renal recovery instead, and only runs recovery analysis.
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

run_recovery_analysis_150 <- function(siteid, base_table, aki_episodes,aki_labs,
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
                                       restrict_model_corr, factor_threshold = 5) {
  currSiteId <- siteid
  peak_trend <- base_table
  aki_index <- aki_episodes
  labs_aki_summ <- aki_labs
  demographics <- demog_table
  demog_time_to_event <- demog_toe
  demog_list <- demog_var_list
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
  factor_cutoff <- factor_threshold
  
  patients_with_preadmit_cr <- preadmit_cr_list
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  
  restrict_models <- restrict_model_corr
  
  file_prefix <- ""
  message("\n============================\nSet-up for Time-to-Event Analysis\n============================")
  if(isTRUE(preadmit_only_analysis)) {
    file_prefix <- "_PreAdmitCrOnly"
    currSiteId <- paste0(currSiteId,file_prefix)
    peak_trend <- peak_trend[peak_trend$patient_id %in% patients_with_preadmit_cr,]
    aki_index <- aki_index[aki_index$patient_id %in% patients_with_preadmit_cr,]
    demographics <- demographics[demographics$patient_id %in% patients_with_preadmit_cr,]
    try({comorbid <- comorbid[comorbid$patient_id %in% patients_with_preadmit_cr,]})
    try({kdigo_grade <- kdigo_grade[kdigo_grade$patient_id %in% patients_with_preadmit_cr,]})
    try({med_coaga_new <- med_coaga_new[med_coaga_new$patient_id %in% patients_with_preadmit_cr,]})
    try({med_coagb_new <- med_coagb_new[med_coagb_new$patient_id %in% patients_with_preadmit_cr,]})
    try({med_covid19_new <- med_covid19_new[med_covid19_new$patient_id %in% patients_with_preadmit_cr,]})
    try({med_acearb_chronic <- med_acearb_chronic[med_acearb_chronic$patient_id %in% patients_with_preadmit_cr,]})
    
    # Set the ratio used for analysis to the ratio calculated from pre-admission Cr only
    peak_trend$ratio <- peak_trend$ratio_prioronly
    
    message("Performing analysis for patients with pre-admission Cr only.")
    message("Recovery threshold: 150%")
  } else {
    message("Performing analysis for ALL patients.")
    message("Recovery threshold: 150%")
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
  # First generate table where we find the time taken to achieve 1.50x baseline ratio
  labs_cr_recovery <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak >= 0) %>% tidyr::fill(severe)
  labs_cr_recovery_tmp <- labs_cr_recovery %>% dplyr::group_by(patient_id) %>% tidyr::complete(time_from_peak = tidyr::full_seq(time_from_peak,1)) %>% dplyr::mutate(ratio = zoo::na.fill(ratio,Inf))
  time_to_ratio1.50 <- labs_cr_recovery_tmp %>% split(.$patient_id) %>% purrr::map(~get_day(.$ratio,.$time_from_peak,target=1.50)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
  colnames(time_to_ratio1.50)[2] <- "time_to_ratio1.50"
  
  labs_aki_summ_index <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission))
  # Get index AKI grade
  index_aki_grade <- labs_aki_summ_index %>% dplyr::select(patient_id,aki_kdigo_final) %>% dplyr::group_by(patient_id) %>% dplyr::filter(aki_kdigo_final == max(aki_kdigo_final)) %>% dplyr::ungroup()
  message("Filtering for AKI patients only...")
  # Filter for AKI cases only
  aki_index_recovery <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::filter(severe %in% c(2,4,5)) %>% dplyr::mutate(severe=ifelse(severe==2,0,1))
  aki_index_recovery <- merge(aki_index_recovery,time_to_ratio1.50,by="patient_id",all.x=TRUE) 
  aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(recover_1.50x = ifelse(is.na(time_to_ratio1.50),0,1))
  message("Computing death and recovery times...")
  # Get death times/censor times
  discharge_day <- demographics %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death_km = dplyr::if_else(deceased==0,as.integer(days_since_admission),as.integer(as.Date(death_date) - as.Date(admission_date)))) %>% dplyr::select(patient_id,deceased,time_to_death_km)
  aki_index_recovery <- merge(aki_index_recovery,discharge_day,by="patient_id",all.x=TRUE)
  #aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_ratio1.50 = dplyr::if_else(recover_1.50x == 0,as.integer(time_to_death_km),as.integer(time_to_ratio1.50)))
  aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_ratio1.50 = dplyr::if_else(recover_1.50x == 0,as.integer(time_to_death_km)-as.integer(peak_cr_time),as.integer(time_to_ratio1.50)))
  
  # Correct the death times for peak Cr
  aki_index_recovery <- aki_index_recovery %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death_km = as.integer(time_to_death_km) - as.integer(peak_cr_time),time_adm_to_death = as.integer(time_to_death_km))
  
  aki_index_recovery <- merge(aki_index_recovery,labs_aki_summ_index[,c("patient_id","aki_kdigo_final")],by="patient_id",all.x=TRUE)
  
  message("\nDoing initial filter for medications with more than one factor level.")
  med_recovery_list <- c("COAGA","COAGB","covid_rx","acei_arb_preexposure")
  med_recovery_list <- med_recovery_list[c(coaga_present,coagb_present,dplyr::if_else((isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)),TRUE,FALSE),dplyr::if_else((isTRUE(acei_present) | isTRUE(arb_present)),TRUE,FALSE))]
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
  med_recovery_tmp <- med_recovery_tmp[,med_recovery_list]
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
  recovery_tmp <- aki_index_recovery[,c("patient_id","recover_1.50x",demog_list,comorbid_recovery_list,med_recovery_list)] %>% as.data.frame()
  
  if(length(comorbid_recovery_list) > 0) {
    comorbid_recovery_list_tmp <- vector(mode="list",length=length(comorbid_recovery_list))
    for(i in 1:length(comorbid_recovery_list)) {
      recovery_tmp1 <- recovery_tmp[,c("patient_id",comorbid_recovery_list[i],"recover_1.50x")]
      recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(comorbid_recovery_list[i]),recover_1.50x)
      recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.50x == 1)
      # if(min(recovery_tmp3$n) >= factor_cutoff) {
      #     comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
      # }
      if(min(recovery_tmp3$n) >= factor_cutoff & nrow(recovery_tmp3) > 1) {
        message(paste0(c("Including ",comorbid_recovery_list[i]," into the comorbid_recovery list...")))
        comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
      }
    }
    comorbid_recovery_list <- unlist(comorbid_recovery_list_tmp[lengths(comorbid_recovery_list_tmp) > 0L])
  }
  
  if(length(demog_list) > 0) {
    demog_recovery_list_tmp <- vector(mode="list",length=length(demog_list))
    for(i in 1:length(demog_list)) {
      recovery_tmp1 <- recovery_tmp[,c("patient_id",demog_list[i],"recover_1.50x")]
      recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(demog_list[i]),recover_1.50x)
      recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.50x == 1)
      # if(min(recovery_tmp3$n) >= factor_cutoff) {
      #     comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
      # }
      if(min(recovery_tmp3$n) >= factor_cutoff & nrow(recovery_tmp3) > 1) {
        message(paste0(c("Including ",demog_list[i]," into the demog_recovery list...")))
        demog_recovery_list_tmp[i] <- demog_list[i]
      }
    }
    demog_recovery_list <- unlist(demog_recovery_list_tmp[lengths(demog_recovery_list_tmp) > 0L])
  }
  
  if(length(med_recovery_list) > 0) {
    med_recovery_list_tmp <- vector(mode="list",length=length(med_recovery_list))
    if(length(med_recovery_list)>0) {
      for(i in 1:length(med_recovery_list)) {
        recovery_tmp1 <- recovery_tmp[,c("patient_id",med_recovery_list[i],"recover_1.50x")]
        recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(med_recovery_list[i]),recover_1.50x)
        recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.50x == 1)
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
  }
  
  message("\nFinal factor list for recovery (before user customisation): ",paste(c(demog_recovery_list,comorbid_recovery_list,med_recovery_list),collapse=" "))
  
  if(restrict_models == TRUE) {
    demog_recovery_list <- demog_recovery_list[demog_recovery_list %in% restrict_list]
    comorbid_recovery_list <- comorbid_recovery_list[comorbid_recovery_list %in% restrict_list]
    med_recovery_list <- med_recovery_list[med_recovery_list %in% restrict_list]
    message(paste("\nAfter filtering for custom-specified variables, we have the following:\nDemographics: ",demog_recovery_list,"\nComorbidities:",comorbid_recovery_list,"\nMedications:",med_recovery_list,sep = " "))
  }
  variable_list_output <- paste(c("Final Recovery variable list:",demog_recovery_list,comorbid_recovery_list,med_recovery_list),collapse=" ")
  readr::write_lines(variable_list_output,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_custom_equation_150.txt")),append=F)
  
  message("Now proceeding to time-to-Cr recovery analysis...")
  # Now run the actual time-to-event analysis
  
  # Kaplan Meier plot for COVID-19 severity
  message("Generating Kaplan-Meier plots...")
  try({
    message("(a) COVID-19 severity")
    recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ severe")
    fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
    plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
    plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
    plot_recover_summ_table <- plot_recover$data.survtable
    write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_Severe_PlotSummStats.csv")),row.names=TRUE)
    write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_Severe_Plot.csv")),row.names=FALSE)
    # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_Severe.png")),plot=plot_recover$plot,width=12,height=12,units="cm")
    plot.new()
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_Severe.png")),plot=print(plot_recover,newpage=F),width=12,height=12,units="cm")
  })
 
  # Kaplan Meier plot for KDIGO grades
  try({
    message("(b) AKI KDIGO stage")
    recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ aki_kdigo_final")
    fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
    plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
    plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
    write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_KDIGO_PlotSummStats.csv")),row.names=TRUE)
    write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_KDIGO_Plot.csv")),row.names=FALSE)
    # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_KDIGO.png")),plot=plot_recover,width=12,height=12,units="cm")
    plot.new()
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_KDIGO.png")),plot=print(plot_recover,newpage=F),width=12,height=12,units="cm")
  })
 
  if(isTRUE(ckd_present)) {
    try({
      # Kaplan Meier plot for CKD
      message("(c) CKD comorbidity")
      recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ckd")
      fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
      plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
      plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
      write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CKD_PlotSummStats.csv")),row.names=TRUE)
      write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CKD_Plot.csv")),row.names=FALSE)
      # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CKD.png")),plot=plot_recover,width=12,height=12,units="cm")
      plot.new()
      ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CKD.png")),plot=print(plot_recover,newpage=F),width=12,height=12,units="cm")
    })
  }
  
  # CoxPH model
  # Generate univariate analyses first
  message("Generating univariate Cox PH models (time to recovery, AKI patients only)...")
  univ_formulas <- tryCatch({
    sapply(c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list), function(x) as.formula(paste('survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ', x)))
  },error = function(c) {
    message("Error running univariate formulae (univ_formulas).")
    return(NULL)
  })
  univ_models <- tryCatch(
    lapply(univ_formulas, function(x){survival::coxph(x, data = aki_index_recovery)}),
    error = function(c) {
      message("Error running univariate formulae (univ_models).")
      return(NULL)
    }
  )
  # Extract data 
  
  try({
    univ_results <- lapply(univ_models,function(x){
      x <- summary(x)
      return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
    })
    univ_results_recover <- do.call("rbind",univ_results)
    write.csv(univ_results_recover,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Univariate.csv")),row.names=TRUE)
  })
  
  message("\nGenerating Model 1 (time to recovery, AKI patients only)...")
  try({
    recovery_model1 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model1 <- recovery_model1[recovery_model1 %in% model1]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model1,collapse="+")))
    message(paste("Formula for Model 1: survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model1,collapse="+")))
    coxph_recover1 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover1_summ <- summary(coxph_recover1) 
    print(coxph_recover1_summ)
    coxph_recover1_hr <- cbind(coxph_recover1_summ$coefficients,coxph_recover1_summ$conf.int)[,-c(6,7)]
    coxph_recover1_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover1_summ$logtest,coxph_recover1_summ$sctest,coxph_recover1_summ$waldtest))
    coxph_recover1_stats2 <- rbind(data.table::as.data.table(coxph_recover1_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover1_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover1_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model1.csv")),row.names=TRUE)
    write.csv(coxph_recover1_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model1_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover1_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model1_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover1_plot <- survminer::ggforest(coxph_recover1,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model1.png")),plot=print(coxph_recover1_plot),width=20,height=20,units="cm")
  })

  message("\nGenerating Model 2 (time to recovery, AKI patients only)...")
  try({
    recovery_model2 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model2 <- recovery_model2[recovery_model2 %in% model2]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model2,collapse="+")))
    message(paste("Formula for Model 2: survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model2,collapse="+")))
    coxph_recover2 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover2_summ <- summary(coxph_recover2) 
    print(coxph_recover2_summ)
    coxph_recover2_hr <- cbind(coxph_recover2_summ$coefficients,coxph_recover2_summ$conf.int)[,-c(6,7)]
    coxph_recover2_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover2_summ$logtest,coxph_recover2_summ$sctest,coxph_recover2_summ$waldtest))
    coxph_recover2_stats2 <- rbind(data.table::as.data.table(coxph_recover2_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover2_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover2_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model2.csv")),row.names=TRUE)
    write.csv(coxph_recover2_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model2_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover2_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model2_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover2_plot <- survminer::ggforest(coxph_recover2,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model2.png")),plot=print(coxph_recover2_plot),width=20,height=20,units="cm")
  })
  
  message("\nGenerating Model 3 (time to recovery, AKI patients only)...")
  try({
    recovery_model3 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model3 <- recovery_model3[recovery_model3 %in% model3]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model3,collapse="+")))
    message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model3,collapse="+")))
    coxph_recover3 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover3_summ <- summary(coxph_recover3) 
    print(coxph_recover3_summ)
    coxph_recover3_hr <- cbind(coxph_recover3_summ$coefficients,coxph_recover3_summ$conf.int)[,-c(6,7)]
    coxph_recover3_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover3_summ$logtest,coxph_recover3_summ$sctest,coxph_recover3_summ$waldtest))
    coxph_recover3_stats2 <- rbind(data.table::as.data.table(coxph_recover3_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover3_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover3_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model3.csv")),row.names=TRUE)
    write.csv(coxph_recover3_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model3_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover3_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model3_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover3_plot <- survminer::ggforest(coxph_recover3,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model3.png")),plot=print(coxph_recover3_plot),width=20,height=20,units="cm")
  })
  message("\nGenerating Model 4 with ACE-i/ARBs (time to recovery, AKI patients only)...")
  try({
    recovery_model4 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
    recovery_model4 <- recovery_model4[recovery_model4 %in% model4]
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model4,collapse="+")))
    message(paste("Formula for Model 4: survival::Surv(time=time_to_ratio1.50,event=recover_1.50x) ~ ",paste(recovery_model4,collapse="+")))
    coxph_recover4 <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover4_summ <- summary(coxph_recover4) 
    print(coxph_recover4_summ)
    coxph_recover4_hr <- cbind(coxph_recover4_summ$coefficients,coxph_recover4_summ$conf.int)[,-c(6,7)]
    coxph_recover4_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_recover4_summ$logtest,coxph_recover4_summ$sctest,coxph_recover4_summ$waldtest))
    coxph_recover4_stats2 <- rbind(data.table::as.data.table(coxph_recover4_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_recover4_summ$rsq,keep.rownames = T))
    write.csv(coxph_recover4_hr,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model4.csv")),row.names=TRUE)
    write.csv(coxph_recover4_stats1,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model4_teststats.csv")),row.names=FALSE,col.names = FALSE)
    write.csv(coxph_recover4_stats2,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model4_concord_rsq.csv")),row.names=FALSE,col.names = FALSE)
    coxph_recover4_plot <- survminer::ggforest(coxph_recover4,data=aki_index_recovery)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover150_CoxPH_Model4.png")),plot=print(coxph_recover4_plot),width=20,height=20,units="cm")
  })
  
  message("If you are getting any errors with model generation - do note that it may actually be normal to get errors\nif your site numbers are low (especially for model 3). Please check your data to see if the appropriate\nnumber of events occur for each factor level.")
  aki_index_recovery
}

