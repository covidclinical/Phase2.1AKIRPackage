#' Creates (non)-obfuscated demographic tables and also pre-processes demographics for use in time-to-event analyses.
#' If specified, can also filter the output to pre-specified patients with pre-admission sCr values.
#' @return A list: demog_toe gives a table with the filtered demographics ready for merging for time-to-event 
#' analyses, while demog_var_list provides the headers of relevant demographics to be included during the automated filtering 
#' process for model construction.
#' @param siteid currSiteId
#' @param demog_table The demographics_filt equivalent
#' @param comorbid_table The comorbid equivalent
#' @param comorbid_header The comorbid_list equivalent
#' @param kdigo_grade_table kdigo_grade
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
#' @param cirrhosis_valid cirrhosis_present (Whether cirrhotic patients are in the cohort)
#' @param obfuscation is_obfuscated
#' @param obfuscation_level obfuscation_value
#' 

generate_demog_files <- function(siteid, demog_table,aki_labs,
                                 comorbid_table,comorbid_header,
                                 kdigo_grade_table,
                                 coaga_valid,coagb_valid,
                                 covid_valid,remdesivir_valid,
                                 acei_valid,arb_valid,
                                 med_coaga = NULL,med_coagb = NULL,med_covid19 = NULL,med_acearb = NULL,
                                 cirrhosis_valid,
                                 obfuscation,obfuscation_level,aki_first_visit_list,
                                 use_custom_output = FALSE,use_custom_output_dir = "/4ceData/Output",ckd_stage_table,earliest_cr_table) {
  currSiteId <- siteid
  demographics_filt <- demog_table
  labs_aki_summ <- aki_labs
  comorbid <- comorbid_table
  comorbid_list <- comorbid_header
  kdigo_grade <- kdigo_grade_table
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
  cirrhosis_present <- cirrhosis_valid
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  aki_list_first_admit <- aki_first_visit_list
  ckd_staging <- ckd_stage_table
  earliest_cr <- earliest_cr_table
  
  if(isTRUE(use_custom_output)) {
    dir.output <- use_custom_output_dir
  } else {
    dir.output <- getProjectOutputDirectory()
  }
  
  file_suffix <- ".csv"
  
  demog_summ <- demographics_filt %>% dplyr::select(patient_id,sex,age_group,race,severe,deceased,time_to_severe,time_to_death) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
  demog_summ <- merge(demog_summ,comorbid,by="patient_id",all.x=TRUE)
  demog_summ <- merge(demog_summ,kdigo_grade,by="patient_id",all.x=TRUE)
  demog_summ <- merge(demog_summ,ckd_staging,by="patient_id",all.x=TRUE)
  demog_summ <- merge(demog_summ,earliest_cr,by="patient_id",all.x=TRUE)
  
  demog_summ$ckd_stage[is.na(demog_summ$ckd_stage)] <- "egfr_90_and_above"
  demog_summ$preadmit_cr_period[is.na(demog_summ$preadmit_cr_period)] <- "zero_to_90_days"
  
  # Add in COAGA and COAGB information
  if(isTRUE(coaga_present)) {
    demog_summ <- merge(demog_summ,med_coaga_new,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$COAGA <- factor(demog_summ$COAGA,levels=c(0,1),labels=c("No Antiplatelets","Antiplatelets"))
  }
  if(isTRUE(coagb_present)) {
    demog_summ <- merge(demog_summ,med_coagb_new,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$COAGB <- factor(demog_summ$COAGB,levels=c(0,1),labels=c("No Anticoagulation","Anticoagulation"))
  }
  # Add COVID-19 antiviral information
  if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
    demog_summ <- merge(demog_summ,med_covid19_new,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$covid_rx <- factor(demog_summ$covid_rx,levels=c(0,1),labels=c("No COVID-19 Therapy","COVID-19 Therapy"))
    if(isTRUE(covid19antiviral_present)) {
      demog_summ$covidviral <- factor(demog_summ$covidviral,levels=c(0,1),labels=c("No Antiviral","Antivirals"))
    }
    if(isTRUE(remdesivir_present)) {
      demog_summ$remdesivir <- factor(demog_summ$remdesivir,levels=c(0,1),labels=c("No Remdesivir","Remdesivir"))
    }
  }
  # Add in RAAS blockade information
  if(isTRUE(acei_present) | isTRUE(arb_present)) {
    demog_summ <- merge(demog_summ,med_acearb_chronic,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$acei_arb_preexposure <- factor(demog_summ$acei_arb_preexposure,levels=c(0,1),labels=c("No ACE-i/ARB","ACE-i and/or ARB"))
  }
  
  demog_summ[is.na(demog_summ)] <- 0
  demog_summ$aki <- 0
  demog_summ$aki[demog_summ$patient_id %in% aki_list_first_admit] <- 1
  
  demog_summ$severe <- factor(demog_summ$severe,levels=c(0,1),labels=c("Non-severe","Severe"))
  demog_summ$deceased <- factor(demog_summ$deceased,levels=c(0,1),labels=c("Alive","Deceased"))
  demog_summ$aki <- factor(demog_summ$aki,levels=c(0,1),labels=c("No AKI","AKI"))
  demog_summ$aki_kdigo_grade <- factor(demog_summ$aki_kdigo_grade,levels=c(0,1,2,3),labels=c("No AKI","Stage 1","Stage 2","Stage 3"))
  demog_summ$ckd_stage <- factor(demog_summ$ckd_stage,levels=c("egfr_90_and_above","egfr_60_to_90","egfr_30_to_60"))
  demog_summ$preadmit_cr_period <- factor(demog_summ$preadmit_cr_period,levels=c("zero_to_90_days","91_to_180_days","181_to_365_days"))
  demog_summ[comorbid_list] <- lapply(demog_summ[comorbid_list],factor)
  demog_summ <- demog_summ %>% dplyr::distinct()
  
  cat(paste0(c("\nObfuscation cutoff: ",obfuscation_value,"\n")))
  # Obfuscation requirements by certain sites
  if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
    comorbid_demog_summ_tmp <- vector(mode="list",length=length(comorbid_list))
    for(i in 1:length(comorbid_list)) {
      demog_summ_tmp1 <- demog_summ[,c("patient_id","aki",comorbid_list[i])]
      demog_summ_tmp2 <- demog_summ_tmp1 %>% dplyr::group_by(aki) %>% dplyr::count(get(comorbid_list[i]))
      cat(paste0(c("\nPerforming demographics table comorbid filtering for: ",comorbid_list[i]," with lowest count ",min(demog_summ_tmp2$n)," and obfuscation cutoff ",obfuscation_value,"\n")))
      if(min(demog_summ_tmp2$n) >= obfuscation_value) {
        comorbid_demog_summ_tmp[i] <- comorbid_list[i]
      }
    }
    comorbid_demog_summ <- unlist(comorbid_demog_summ_tmp[lengths(comorbid_demog_summ_tmp) > 0L])
  } else {
    comorbid_demog_summ <- comorbid_list
  }
  
  med_summ <- NULL
  if(isTRUE(coaga_present) & isTRUE(coagb_present)) {
    med_summ <- c("COAGA","COAGB")
  } else if(isTRUE(coaga_present)) {
    med_summ <- "COAGA"
  } else if(isTRUE(coagb_present)) {
    med_summ <- "COAGB"
  }
  
  if(isTRUE(acei_present) | isTRUE(arb_present)) {
    med_summ <- c(med_summ,"acei_arb_preexposure")
  }
  
  if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
    med_summ <- c(med_summ,"covid_rx")
    if(isTRUE(covid19antiviral_present)) {
      med_summ <- c(med_summ,"covidviral")
    }
    if(isTRUE(remdesivir_present)) {
      med_summ <- c(med_summ,"remdesivir")
    }
  }
  
  table_one_vars <- c("sex","age_group","race","severe","deceased","aki_kdigo_grade","ckd_stage","preadmit_cr_period",comorbid_demog_summ,med_summ)
  
  # Create obfuscated table one for sites which require it
  if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
    # First print table stratified by AKI status
    demog_obf <- demog_summ %>% dplyr::group_by(aki) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "aki",values_from = "n")
    demog_obf$category <- "n"
    demog_obf <- demog_obf[,c(3,1,2)]
    demog_obf$p_val <- NA
    colnames(demog_obf) <- c("category","No_AKI","AKI","p_val")
    for(i in 1:length(table_one_vars)) {
      try({
        tmp <- demog_summ %>% dplyr::group_by(aki) %>% dplyr::count(get(table_one_vars[i])) %>% tidyr::pivot_wider(names_from="aki",values_from = "n")
        tmp[is.na(tmp)] <- 0
        colnames(tmp) <- c("category","No_AKI","AKI")
        tmp <- tmp %>% dplyr::mutate(category = paste0(table_one_vars[i],"_",category))
        
        # tryCatch statements attempt to catch instances where the Fisher's test may fail due to insufficient convergent cycles
        p_value <- tryCatch({
          cat(paste0(c("Attempting Fisher's test for ",table_one_vars[i],"\n")))
          fisher.test(data.frame(tmp[-1],row.names=tmp$category))$p.value
        }, error = function(e) {
          message("Failed running Fisher's exact test for ",table_one_vars[i],", proceeding with Monte Carlo simulation.")
          tryCatch({
            fisher.test(data.frame(tmp[-1],row.names=tmp$category),simulate.p.value=TRUE)$p.value
          },error=function(e){
            message("Failed running Fisher's exact test even with Monte Carlo simulation. A default P-value of NA will be printed - please inspect data.")
            return(NA)
          })
        })
        tmp$p_val = p_value
        colnames(tmp) <- c("category","No_AKI","AKI","p_val")
        demog_obf <- rbind(demog_obf,tmp)
        rm(tmp)
        rm(p_value)
      })
    }
    demog_obf$No_AKI[demog_obf$No_AKI < obfuscation_value] <- 0
    demog_obf$AKI[demog_obf$AKI < obfuscation_value] <- 0
    demog_obf <- demog_obf %>% dplyr::group_by(category) %>% dplyr::mutate(total = No_AKI + AKI) %>% dplyr::ungroup()
    no_aki_total <- demog_obf$No_AKI[1]
    aki_total <- demog_obf$AKI[1]
    total_pop <- demog_obf$total[1]
    demog_obf <- demog_obf %>% dplyr::group_by(category) %>% dplyr::mutate(No_AKI_perc = No_AKI / no_aki_total * 100, AKI_perc = AKI/aki_total * 100, total_perc = total/total_pop * 100) %>% dplyr::ungroup()
    write.csv(demog_obf,file=file.path(dir.output, paste0(currSiteId, "_TableOne_obfuscated",file_suffix)),row.names=F,na="NA")
    
    # Then print table stratified by CKD staging status
    cat("\nAttempting to generate demographics table stratified by CKD stage and obfuscated\n")
    try({
      demog_ckd_obf <- demog_summ %>% dplyr::group_by(ckd_stage) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "ckd_stage",values_from = "n")
      ckd_stage_ref <- c("egfr_90_and_above","egfr_60_to_90","egfr_30_to_60")
      ckd_stage_ref <- ckd_stage_ref[ckd_stage_ref %in% colnames(demog_ckd_obf)]
      demog_ckd_obf$category <- "n"
      demog_ckd_obf <- demog_ckd_obf[,c("category",ckd_stage_ref)]
      demog_ckd_obf$p_val <- NA
      for(i in 1:length(table_one_vars)) {
        try({
          tmp <- demog_summ %>% dplyr::group_by(ckd_stage) %>% dplyr::count(get(table_one_vars[i])) %>% tidyr::pivot_wider(names_from="ckd_stage",values_from = "n")
          tmp[is.na(tmp)] <- 0
          colnames(tmp) <- c("category",ckd_stage_ref)
          tmp <- tmp %>% dplyr::mutate(category = paste0(table_one_vars[i],"_",category))
          
          # tryCatch statements attempt to catch instances where the Fisher's test may fail due to insufficient convergent cycles
          p_value <- tryCatch({
            cat(paste0(c("Attempting Fisher's test for ",table_one_vars[i],"\n")))
            fisher.test(data.frame(tmp[-1],row.names=tmp$category))$p.value
          }, error = function(e) {
            message("Failed running Fisher's exact test for ",table_one_vars[i],", proceeding with Monte Carlo simulation.")
            tryCatch({
              fisher.test(data.frame(tmp[-1],row.names=tmp$category),simulate.p.value=TRUE)$p.value
            },error=function(e){
              message("Failed running Fisher's exact test even with Monte Carlo simulation. A default P-value of NA will be printed - please inspect data.")
              return(NA)
            })
          })
          tmp$p_val = p_value
          colnames(tmp) <- c("category",ckd_stage_ref,"p_val")
          demog_ckd_obf <- rbind(demog_ckd_obf,tmp)
          rm(tmp)
          rm(p_value)
        })
      }
      
      demog_ckd_obf$total <- 0
      demog_ckd_obf$total[demog_ckd_obf$category == "n"] <- sum(demog_ckd_obf[demog_ckd_obf$category == "n",ckd_stage_ref],na.rm=T)
      total_pop <- demog_ckd_obf$total[1]
      if("egfr_90_and_above" %in% ckd_stage_ref) {
        demog_ckd_obf$egfr_90_and_above[demog_ckd_obf$egfr_90_and_above < obfuscation_value] <- 0
        egfr_90_and_above_total <- demog_ckd_obf$egfr_90_and_above[1]
        demog_ckd_obf$total <- demog_ckd_obf$total + demog_ckd_obf$egfr_90_and_above
        demog_ckd_obf <- demog_ckd_obf %>% dplyr::group_by(category) %>% dplyr::mutate(egfr_90_and_above_perc = egfr_90_and_above / egfr_90_and_above_total * 100) %>% dplyr::ungroup()
      }
      if("egfr_60_to_90" %in% ckd_stage_ref) {
        demog_ckd_obf$egfr_60_to_90[demog_ckd_obf$egfr_60_to_90 < obfuscation_value] <- 0
        egfr_60_to_90_total <- demog_ckd_obf$egfr_60_to_90[1]
        demog_ckd_obf$total <- demog_ckd_obf$total + demog_ckd_obf$egfr_60_to_90
        demog_ckd_obf <- demog_ckd_obf %>% dplyr::group_by(category) %>% dplyr::mutate(egfr_60_to_90_perc = egfr_60_to_90 / egfr_60_to_90_total * 100) %>% dplyr::ungroup()
      }
      if("egfr_30_to_60" %in% ckd_stage_ref) {
        demog_ckd_obf$egfr_30_to_60[demog_ckd_obf$egfr_30_to_60 < obfuscation_value] <- 0
        egfr_30_to_60_total <- demog_ckd_obf$egfr_30_to_60[1]
        demog_ckd_obf$total <- demog_ckd_obf$total + demog_ckd_obf$egfr_30_to_60
        demog_ckd_obf <- demog_ckd_obf %>% dplyr::group_by(category) %>% dplyr::mutate(egfr_30_to_60_perc = egfr_30_to_60 / egfr_30_to_60_total * 100) %>% dplyr::ungroup()
      }
      demog_ckd_obf <- demog_ckd_obf %>% dplyr::group_by(category) %>% dplyr::mutate(total_perc = total / total_pop * 100) %>% dplyr::ungroup()
      write.csv(demog_ckd_obf,file=file.path(dir.output, paste0(currSiteId, "_TableOne_obfuscated_CKDStage",file_suffix)),row.names=F,na="NA")
      
    })
  }
  
  if(obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
    table_one <- tableone::CreateTableOne(data=demog_summ,vars=table_one_vars,strata="aki",addOverall = T)
    export_table_one <- print(table_one,showAllLevels=TRUE,formatOptions=list(big.mark=","))
    write.csv(export_table_one,file=file.path(dir.output, paste0(currSiteId, "_TableOne",file_suffix)))
    
    try({
      table_one_ckd <- tableone::CreateTableOne(data=demog_summ,vars=table_one_vars,strata="ckd_stage",addOverall = T)
      export_table_one_ckd <- print(table_one_ckd,showAllLevels=TRUE,formatOptions=list(big.mark=","))
      write.csv(export_table_one_ckd,file=file.path(dir.output, paste0(currSiteId, "_TableOne_CKD",file_suffix)))
    })
  }
  cat("\nTableOne with patient demographics should have been generated in CSV files at this point. Check for any errors.\n")
  
  demog_time_to_event_tmp <- demog_summ[,c("sex","age_group","race")] 
  demog_time_to_event_tmp <- data.table::as.data.table(lapply(demog_time_to_event_tmp,factor))
  demog_time_to_event_tmp <- data.table::as.data.table(demog_time_to_event_tmp)[,sapply(demog_time_to_event_tmp,function(col) nlevels(col) > 1),with=FALSE] 
  demog_list <- colnames(demog_time_to_event_tmp)
  demog_time_to_event <- demog_summ[,c("patient_id",demog_list)]
  try({demog_time_to_event <- demog_time_to_event %>% dplyr::group_by(patient_id) %>% dplyr::mutate(age_group = dplyr::if_else(age_group == "70to79" | age_group == "80plus","70_and_above","below_70")) %>% dplyr::ungroup()})
  #in case there is only one age_group level
  
  # Deals with the special case where there is a sex category of others (e.g. KUMC)
  try({demog_time_to_event <- demog_time_to_event %>% dplyr::group_by(patient_id) %>% dplyr::mutate(sex = dplyr::if_else(sex == "male","male","female")) %>% dplyr::ungroup()})
  # in case there is only one sex level
  demog_output <- list("demog_toe" = demog_time_to_event, "demog_var_list" = demog_list)
  return(demog_output)
}
