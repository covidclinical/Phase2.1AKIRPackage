#' Generates risk tables for combining into a single aggregated life table for computing cumulative incidence curves for AKI recovery and New Onset CKD (competing risk is mortality)
#' @param aki_index_recovery The aki_index_recovery table generated from run_time_to_event_analysis()
#' @param aki_index_nonckd The aki_index_nonckd table generated from run_time_to_event_analysis_nonckd()
#' @param aki_index_nonckd_akionly The aki_index_nonckd_akionly table generated from run_time_to_event_analysis_nonckd()
#' @param aki_index_ckdonly_akionly The aki_index_ckdonly_akionly table generated from run_time_to_event_analysis_nonckd()

generate_life_table_competing_risk <- function(currSiteId,aki_index_recovery,aki_index_nonckd,aki_index_nonckd_akionly,aki_index_ckdonly_akionly, use_custom_output = FALSE,use_custom_output_dir = '/4ceData/Output') {
  cat("\nGenerating preliminary tables for calculating the cause-specific cumulatic incidence curves for
      \n1) AKI recovery
      \n2) Time to New CKD")
  
  if(isTRUE(use_custom_output)) {
    dir.output <- use_custom_output_dir
  } else {
    dir.output <- getProjectOutputDirectory()
  }
  
  variables <- c("patient_id","aki_kdigo_final","severe","ckd","time_to_death_km","time_to_ratio1.25","time_to_new_ckd","deceased","recover_1.25x","new_ckd")
  table_names <- c("AKI_Recovery_All_AKI","New_CKD_All_NonCKD","New_CKD_NonCKD_and_AKI_Only","AKI_Recovery_NonCKD_Only","AKI_Recovery_CKD_Only")
  data_list <- list(aki_index_recovery,aki_index_nonckd,aki_index_nonckd_akionly,aki_index_ckdonly_akionly)
  for(i in 1:4) {
    message(i)
    tmp <- data_list[[i]][,variables[variables %in% colnames(data_list[[i]])]]
    if(i == 1) { # aki_index_recovery - recovery outcome
      tmp <- tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
        recover_1.25x == 1 ~ 1,
        recover_1.25x == 0 & deceased == 1 ~ 2,
        TRUE ~ 0
      )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
        event == 0 ~ time_to_death_km,
        event == 1 ~ time_to_ratio1.25,
        event == 2 ~ time_to_death_km
      )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T) %>% dplyr::select(patient_id,aki_kdigo_final,severe,ckd,time_to_event,event)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,aki_kdigo_final) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_KDIGO.csv")),row.names=FALSE)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,severe) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_Severe.csv")),row.names=FALSE)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,ckd) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_CKD.csv")),row.names=FALSE)
    } else if (i == 2) { #Non-CKD, All - New CKD outcome
      tmp <- tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
        new_ckd == 1 ~ 1,
        new_ckd == 0 & deceased == 1 ~ 2,
        TRUE ~ 0
      )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
        event == 0 ~ time_to_death_km,
        event == 1 ~ time_to_new_ckd,
        event == 2 ~ time_to_death_km
      )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T) %>% dplyr::select(patient_id,aki_kdigo_final,severe,ckd,time_to_event,event)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,aki_kdigo_final) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_KDIGO.csv")),row.names=FALSE)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,severe) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_Severe.csv")),row.names=FALSE)
    } else if(i == 3) { # Non-CKD, AKI Only - two outcomes, New CKD onset + AKI Recovery
      tmp_1 <- tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
        new_ckd == 1 ~ 1,
        new_ckd == 0 & deceased == 1 ~ 2,
        TRUE ~ 0
      )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
        event == 0 ~ time_to_death_km,
        event == 1 ~ time_to_new_ckd,
        event == 2 ~ time_to_death_km
      )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T) %>% dplyr::select(patient_id,aki_kdigo_final,severe,time_to_event,event)
      tmp_table <- tmp_1 %>% dplyr::group_by(time_to_event,aki_kdigo_final) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_KDIGO.csv")),row.names=FALSE)
      tmp_table <- tmp_1 %>% dplyr::group_by(time_to_event,severe) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i],"_Severe.csv")),row.names=FALSE)
      
      tmp_2 <- tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
        new_ckd == 1 ~ 1,
        new_ckd == 0 & deceased == 1 ~ 2,
        TRUE ~ 0
      )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
        event == 0 ~ time_to_death_km,
        event == 1 ~ time_to_new_ckd,
        event == 2 ~ time_to_death_km
      )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T) %>% dplyr::select(patient_id,aki_kdigo_final,severe,time_to_event,event)
      tmp_table <- tmp_2 %>% dplyr::group_by(time_to_event,aki_kdigo_final) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i+1],"_KDIGO.csv")),row.names=FALSE)
      tmp_table <- tmp_2 %>% dplyr::group_by(time_to_event,severe) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i+1],"_Severe.csv")),row.names=FALSE)
    } else {
      tmp <- tmp %>% dplyr::group_by(patient_id) %>% dplyr::mutate(event = dplyr::case_when(
        recover_1.25x == 1 ~ 1,
        recover_1.25x == 0 & deceased == 1 ~ 2,
        TRUE ~ 0
      )) %>% dplyr::mutate(time_to_event = dplyr::case_when(
        event == 0 ~ time_to_death_km,
        event == 1 ~ time_to_ratio1.25,
        event == 2 ~ time_to_death_km
      )) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,.keep_all = T) %>% dplyr::select(patient_id,aki_kdigo_final,severe,time_to_event,event)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,aki_kdigo_final) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i+1],"_KDIGO.csv")),row.names=FALSE)
      tmp_table <- tmp %>% dplyr::group_by(time_to_event,severe) %>% dplyr::summarise(recover = length(event[event == 1]),deceased = length(event[event == 2]),censor = length(event[event == 0])) %>% dplyr::ungroup()
      write.csv(tmp_table,file=file.path(dir.output, paste0(currSiteId, "_TimeToEvent_CIC_",table_names[i+1],"_Severe.csv")),row.names=FALSE)
    }
  }
  cat("Done!")
}
