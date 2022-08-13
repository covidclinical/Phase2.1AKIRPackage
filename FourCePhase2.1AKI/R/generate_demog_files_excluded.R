#' Creates (non)-obfuscated demographic tables for excluded patients and also pre-processes demographics for use in time-to-event analyses.
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

generate_demog_files_excluded <- function(siteid,
                                          demog_table,
                                          comorbid_table,
                                          comorbid_header,
                                          coaga_valid,
                                          coagb_valid,
                                          covid_valid,
                                          remdesivir_valid,
                                          acei_valid,
                                          arb_valid,
                                          med_coaga = NULL,
                                          med_coagb = NULL,
                                          med_covid19 = NULL,
                                          med_acearb = NULL,
                                          obfuscation,
                                          obfuscation_level,
                                          use_custom_output = FALSE,
                                          use_custom_output_dir = "/4ceData/Output") {
  currSiteId <- siteid
  demographics_filt <- demog_table
  comorbid <- comorbid_table
  comorbid_list <- comorbid_header
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
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  
  if (isTRUE(use_custom_output)) {
    dir.output <- use_custom_output_dir
  } else {
    dir.output <- getProjectOutputDirectory()
  }
  
  file_suffix <- ".csv"
  
  demog_summ <-
    demographics_filt %>% dplyr::select(patient_id,
                                        sex,
                                        age_group,
                                        race,
                                        severe,
                                        deceased,
                                        time_to_severe,
                                        time_to_death,excluded_type) %>% dplyr::distinct(patient_id, .keep_all = TRUE)
  demog_summ <-
    merge(demog_summ, comorbid, by = "patient_id", all.x = TRUE)
  
  # Add in COAGA and COAGB information
  if (isTRUE(coaga_present)) {
    demog_summ <-
      merge(demog_summ,
            med_coaga_new,
            by = "patient_id",
            all.x = TRUE) %>% dplyr::distinct(patient_id, .keep_all = TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$COAGA <-
      factor(
        demog_summ$COAGA,
        levels = c(0, 1),
        labels = c("No Antiplatelets", "Antiplatelets")
      )
  }
  if (isTRUE(coagb_present)) {
    demog_summ <-
      merge(demog_summ,
            med_coagb_new,
            by = "patient_id",
            all.x = TRUE) %>% dplyr::distinct(patient_id, .keep_all = TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$COAGB <-
      factor(
        demog_summ$COAGB,
        levels = c(0, 1),
        labels = c("No Anticoagulation", "Anticoagulation")
      )
  }
  # Add COVID-19 antiviral information
  if (isTRUE(covid19antiviral_present) |
      isTRUE(remdesivir_present)) {
    demog_summ <-
      merge(demog_summ,
            med_covid19_new,
            by = "patient_id",
            all.x = TRUE) %>% dplyr::distinct(patient_id, .keep_all = TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$covid_rx <-
      factor(
        demog_summ$covid_rx,
        levels = c(0, 1),
        labels = c("No COVID-19 Therapy", "COVID-19 Therapy")
      )
    if (isTRUE(covid19antiviral_present)) {
      demog_summ$covidviral <-
        factor(
          demog_summ$covidviral,
          levels = c(0, 1),
          labels = c("No Antiviral", "Antivirals")
        )
    }
    if (isTRUE(remdesivir_present)) {
      demog_summ$remdesivir <-
        factor(
          demog_summ$remdesivir,
          levels = c(0, 1),
          labels = c("No Remdesivir", "Remdesivir")
        )
    }
  }
  # Add in RAAS blockade information
  if (isTRUE(acei_present) | isTRUE(arb_present)) {
    demog_summ <-
      merge(demog_summ,
            med_acearb_chronic,
            by = "patient_id",
            all.x = TRUE) %>% dplyr::distinct(patient_id, .keep_all = TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$acei_arb_preexposure <-
      factor(
        demog_summ$acei_arb_preexposure,
        levels = c(0, 1),
        labels = c("No ACE-i/ARB", "ACE-i and/or ARB")
      )
  }
  
  demog_summ[is.na(demog_summ)] <- 0
  demog_summ$excluded_type <- factor(demog_summ$excluded_type)
  demog_summ$severe <-
    factor(
      demog_summ$severe,
      levels = c(0, 1),
      labels = c("Non-severe", "Severe")
    )
  demog_summ$deceased <-
    factor(
      demog_summ$deceased,
      levels = c(0, 1),
      labels = c("Alive", "Deceased")
    )
  demog_summ[comorbid_list] <-
    lapply(demog_summ[comorbid_list], factor)
  demog_summ <- demog_summ %>% dplyr::distinct()
  
  cat(paste0(c(
    "\nObfuscation cutoff: ", obfuscation_value, "\n"
  )))
  # Obfuscation requirements by certain sites
  if (isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
    comorbid_demog_summ_tmp <-
      vector(mode = "list", length = length(comorbid_list))
    for (i in 1:length(comorbid_list)) {
      demog_summ_tmp1 <- demog_summ[, c("patient_id", comorbid_list[i])]
      demog_summ_tmp2 <-
        demog_summ_tmp1 %>% dplyr::count(get(comorbid_list[i]))
      cat(paste0(
        c(
          "\nPerforming demographics table comorbid filtering for: ",
          comorbid_list[i],
          " with lowest count ",
          min(demog_summ_tmp2$n),
          " and obfuscation cutoff ",
          obfuscation_value,
          "\n"
        )
      ))
      if (min(demog_summ_tmp2$n) >= obfuscation_value) {
        comorbid_demog_summ_tmp[i] <- comorbid_list[i]
      }
    }
    comorbid_demog_summ <-
      unlist(comorbid_demog_summ_tmp[lengths(comorbid_demog_summ_tmp) > 0L])
  } else {
    comorbid_demog_summ <- comorbid_list
  }
  
  med_summ <- NULL
  if (isTRUE(coaga_present) & isTRUE(coagb_present)) {
    med_summ <- c("COAGA", "COAGB")
  } else if (isTRUE(coaga_present)) {
    med_summ <- "COAGA"
  } else if (isTRUE(coagb_present)) {
    med_summ <- "COAGB"
  }
  
  if (isTRUE(acei_present) | isTRUE(arb_present)) {
    med_summ <- c(med_summ, "acei_arb_preexposure")
  }
  
  if (isTRUE(covid19antiviral_present) |
      isTRUE(remdesivir_present)) {
    med_summ <- c(med_summ, "covid_rx")
    if (isTRUE(covid19antiviral_present)) {
      med_summ <- c(med_summ, "covidviral")
    }
    if (isTRUE(remdesivir_present)) {
      med_summ <- c(med_summ, "remdesivir")
    }
  }
  
  table_one_vars <-
    c("sex",
      "age_group",
      "race",
      "severe",
      "deceased","excluded_type",
      comorbid_demog_summ,
      med_summ)
  
  # Create obfuscated table one for sites which require it
  if (isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
    # First print table stratified by AKI status
    excluded_type_list <- c("rrt_surrogate","rrt_procedure_code","baseline_more_than_2.25","insufficient_cr")
    excluded_type_list <- excluded_type_list[excluded_type_list %in% unique(demog_summ$excluded_type)]
    demog_obf <- demog_summ %>% dplyr::group_by(excluded_type) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "excluded_type",values_from = "n")
    demog_obf$category <- "n"
    demog_obf <- demog_obf[,c("category",excluded_type_list)]
    demog_obf <- demog_obf %>% dplyr::mutate(total= sum(dplyr::c_across(tidyselect::vars_select_helpers$where(is.numeric))))
    demog_obf$p_val <- NA_real_
    for (i in 1:length(table_one_vars)) {
      try({
        tmp <- demog_summ %>% dplyr::group_by(excluded_type) %>% dplyr::count(get(table_one_vars[i])) %>% tidyr::pivot_wider(names_from = "excluded_type",values_from = "n")
        colnames(tmp)[1] <- "category"
        tmp[is.na(tmp)] <- 0
        tmp <- tmp %>% dplyr::mutate(category = paste0(table_one_vars[i], "_", category))
        tmp <- tmp[,c("category",excluded_type_list)]
        
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
        tmp <- tmp %>% dplyr::rowwise() %>% dplyr::mutate(total= sum(dplyr::c_across(tidyselect::vars_select_helpers$where(is.numeric))))
        tmp <- tmp %>% dplyr::mutate_if(is.numeric, ~ . * (. > obfuscation_value))
        tmp$p_val = p_value
        colnames(tmp) <- c("category", excluded_type_list,"total","p_val")
        demog_obf <- rbind(demog_obf, tmp)
        rm(tmp)
        rm(p_value)
      })
    }
    demog_obf$total[demog_obf$total < obfuscation_value] <- 0
    write.csv(
      demog_obf,
      file = file.path(
        dir.output,
        paste0(
          currSiteId,
          "_TableOne_Excluded_Pts_obfuscated",
          file_suffix
        )
      ),
      row.names = F,
      na = "NA"
    )
  }
  
  if (obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
    table_one <-
      tableone::CreateTableOne(data = demog_summ, vars = table_one_vars,strata="excluded_type",addOverall = T)
    export_table_one <-
      print(table_one,
            showAllLevels = TRUE,
            formatOptions = list(big.mark = ","))
    write.csv(export_table_one, file = file.path(
      dir.output,
      paste0(currSiteId, "_TableOne_Excluded_Pts", file_suffix)
    ))

  }
  cat(
    "\nTableOne with patient demographics should have been generated in CSV files at this point. Check for any errors.\n"
  )

}
