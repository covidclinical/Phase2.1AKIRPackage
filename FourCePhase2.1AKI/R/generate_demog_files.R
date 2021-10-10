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
#' @param meld_valid meld_analysis_valid
#' @param labs_meld_admit labs_meld_admission
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
#' @param preadmit_cr_list Pre-specified patient list with pre-admission sCr labs. Default value is NULL
#' @param preadmit_only_analysis Specifies whether demographics are to be narrowed down to list specified by preadmit_cr_list. Default is FALSE.
#' @param obfuscation is_obfuscated
#' @param obfuscation_level obfuscation_value
#' 

generate_demog_files <- function(siteid, demog_table,aki_labs,
                                 comorbid_table,comorbid_header,
                                 kdigo_grade_table,
                                 meld_valid = FALSE,labs_meld_admit = NULL,
                                 coaga_valid,coagb_valid,
                                 covid_valid,remdesivir_valid,
                                 acei_valid,arb_valid,
                                 med_coaga = NULL,med_coagb = NULL,med_covid19 = NULL,med_acearb = NULL,
                                 cirrhosis_valid,
                                 preadmit_cr_list = NULL,preadmit_only_analysis = FALSE,
                                 obfuscation,obfuscation_level,aki_first_visit_list) {
  currSiteId <- siteid
  demographics_filt <- demog_table
  labs_aki_summ <- aki_labs
  comorbid <- comorbid_table
  comorbid_list <- comorbid_header
  kdigo_grade <- kdigo_grade_table
  meld_analysis_valid <- meld_valid
  labs_meld_admission <- labs_meld_admit
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
  patients_with_preadmit_cr <- preadmit_cr_list
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  aki_list_first_admit <- aki_first_visit_list
  
  file_suffix <- ".csv"
  if(isTRUE(preadmit_only_analysis)) {
    file_suffix <- "_preadmit_cr_only.csv"
    demographics_filt <- demographics_filt[demographics_filt$patient_id %in% patients_with_preadmit_cr,]
    try({comorbid <- comorbid[comorbid$patient_id %in% patients_with_preadmit_cr,]})
    try({kdigo_grade <- kdigo_grade[kdigo_grade$patient_id %in% patients_with_preadmit_cr,]})
    try({med_coaga_new <- med_coaga_new[med_coaga_new$patient_id %in% patients_with_preadmit_cr,]})
    try({med_coagb_new <- med_coagb_new[med_coagb_new$patient_id %in% patients_with_preadmit_cr,]})
    try({med_covid19_new <- med_covid19_new[med_covid19_new$patient_id %in% patients_with_preadmit_cr,]})
    try({med_acearb_chronic <- med_acearb_chronic[med_acearb_chronic$patient_id %in% patients_with_preadmit_cr,]})
    try({labs_meld_admission <- labs_meld_admission[labs_meld_admission$patient_id %in% patients_with_preadmit_cr,]})
  }
  
  demog_summ <- demographics_filt %>% dplyr::select(patient_id,sex,age_group,race,severe,deceased,time_to_severe,time_to_death) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
  demog_summ <- merge(demog_summ,comorbid,by="patient_id",all.x=TRUE)
  demog_summ <- merge(demog_summ,kdigo_grade,by="patient_id",all.x=TRUE)
  
  if(isTRUE(meld_analysis_valid)) {
    demog_summ <- merge(demog_summ,labs_meld_admission[,-2],by="patient_id",all.x=TRUE)
    demog_summ$meld_admit_severe[is.na(demog_summ$meld_admit_severe)] <- 0
    demog_summ$meld_admit_severe <- factor(demog_summ$meld_admit_severe,levels=c(0,1),labels = c("MELD < 20","MELD >= 20"))
  }
  
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
    demog_summ$covid_rx <- factor(demog_summ$covid_rx,levels=c(0,1),labels=c("No Novel Antiviral","Novel Antiviral"))
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
  
  demog_summ$preadmit_cr <- 0
  demog_summ$preadmit_cr[demog_summ$patient_id %in% patients_with_preadmit_cr] <- 1
  
  demog_summ$severe <- factor(demog_summ$severe,levels=c(0,1),labels=c("Non-severe","Severe"))
  demog_summ$deceased <- factor(demog_summ$deceased,levels=c(0,1),labels=c("Alive","Deceased"))
  demog_summ$aki <- factor(demog_summ$aki,levels=c(0,1),labels=c("No AKI","AKI"))
  demog_summ$aki_kdigo_grade <- factor(demog_summ$aki_kdigo_grade,levels=c(0,1,2,3),labels=c("No AKI","Stage 1","Stage 2","Stage 3"))
  demog_summ$preadmit_cr <- factor(demog_summ$preadmit_cr,levels=c(0,1),labels=c("No Pre-admission sCr","Pre-admission sCr Present"))
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
  }
  
  table_one_vars <- c("sex","age_group","race","severe","deceased","aki_kdigo_grade","preadmit_cr",comorbid_demog_summ,med_summ)
  if(isTRUE(meld_analysis_valid)) {
    table_one_meld_vars <- c("sex","age_group","race","aki","meld_admit_severe","deceased","aki_kdigo_grade","preadmit_cr",comorbid_demog_summ,med_summ)
  }
  #capture.output(summary(table_one),file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne_Missingness.txt")))
  demog_meld_files <- NULL
  if(isTRUE(cirrhosis_present)) {
    demog_cld_summ <- demog_summ %>% dplyr::filter(cld == 1)
  }
  
  # Create obfuscated table one for sites which require it
  if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
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
    demog_obf <- demog_obf %>% dplyr::group_by(category) %>% dplyr::mutate(No_AKI_perc = No_AKI / no_aki_total * 100, AKI_perc = AKI/aki_total * 100, total_perc = total/total_pop) %>% dplyr::ungroup()
    write.csv(demog_obf,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne_obfuscated",file_suffix)),row.names=F,na="NA")
    
    tryCatch({
      if(isTRUE(cirrhosis_present)) {
        # First generate the subgroup table stratified by AKI for cirrhotic patients only
        cat("\nCreating temporary demographics table for cirrhotic patients\n")
        demog_cld_obf <- demog_cld_summ %>% dplyr::group_by(aki) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "aki",values_from = "n")
        demog_cld_obf$category <- "n"
        demog_cld_obf <- demog_cld_obf[,c(3,1,2)]
        demog_cld_obf$p_val <- NA
        colnames(demog_cld_obf) <- c("category","No_AKI","AKI","p_val")
        cat("\nFiltering variables for cirrhotic patient demographics table...\n")
        for(i in 1:length(table_one_vars)) {
          try({
            tmp <- demog_cld_summ %>% dplyr::group_by(aki) %>% dplyr::count(get(table_one_vars[i])) %>% tidyr::pivot_wider(names_from="aki",values_from = "n")
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
            demog_cld_obf <- rbind(demog_cld_obf,tmp)
            rm(tmp)
            rm(p_value)
          })
        }
        cat("\nComputing counts and percentages\n")
        demog_cld_obf$No_AKI[demog_cld_obf$No_AKI < obfuscation_value] <- 0
        demog_cld_obf$AKI[demog_cld_obf$AKI < obfuscation_value] <- 0
        demog_cld_obf <- demog_cld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(total = No_AKI + AKI) %>% dplyr::ungroup()
        no_aki_total <- demog_cld_obf$No_AKI[1]
        aki_total <- demog_cld_obf$AKI[1]
        total_pop <- demog_cld_obf$total[1]
        demog_cld_obf <- demog_cld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(No_AKI_perc = No_AKI / no_aki_total * 100, AKI_perc = AKI/aki_total * 100, total_perc = total/total_pop) %>% dplyr::ungroup()
        demog_meld_files <- "demog_cld_obf"
        cat("\nNo issues with creating demog_cld_obf\n")
      }
    },error = function(e) {
      cat("\nHaving issues generating demographic table for cirrhosis patients only. Check for any error messages that appear.")
      cat("\nOriginal error message:\n",e,"\n")
      message("Having issues generating demographic table for cirrhosis patients only. Check for any error messages that appear.")
      message("Original error message:\n",e)
    })
    try({
      if(isTRUE(meld_analysis_valid)) {
        # If possible to split by MELD, then generate the demographics table split by MELD score cutoff 20
        cat("\nMELD analysis possible. Generating demographics table for cirrhotics stratified by MELD.\n")
        demog_meld_summ <- demog_cld_summ[demog_cld_summ$meld_labs_valid == 1,]
        demog_meld_obf <- demog_meld_summ %>% dplyr::group_by(meld_admit_severe) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "meld_admit_severe",values_from = "n")
        demog_meld_obf$category <- "n"
        demog_meld_obf <- demog_meld_obf[,c(3,1,2)]
        demog_meld_obf$p_val <- NA
        colnames(demog_meld_obf) <- c("category","MELD_less20","MELD_20ormore","p_val")
        cat("\nPerforming filtering of variables for demographics table for cirrhotics by MELD\n")
        for(i in 1:length(table_one_meld_vars)) {
          try({
            tmp <- demog_meld_summ %>% dplyr::group_by(meld_admit_severe) %>% dplyr::count(get(table_one_meld_vars[i])) %>% tidyr::pivot_wider(names_from="meld_admit_severe",values_from = "n")
            tmp[is.na(tmp)] <- 0
            colnames(tmp) <- c("category","MELD_less20","MELD_20ormore")
            tmp <- tmp %>% dplyr::mutate(category = paste0(table_one_meld_vars[i],"|",category))
            
            # tryCatch statements attempt to catch instances where the Fisher's test may fail due to insufficient convergent cycles
            p_value <- tryCatch({
              cat(paste0(c("Attempting Fisher's test for ",table_one_meld_vars[i],"\n")))
              fisher.test(data.frame(tmp[-1],row.names=tmp$category))$p.value
            }, error = function(e) {
              message("Failed running Fisher's exact test for ",table_one_meld_vars[i],", proceeding with Monte Carlo simulation.")
              tryCatch({
                fisher.test(data.frame(tmp[-1],row.names=tmp$category),simulate.p.value=TRUE)$p.value
              },error=function(e){
                message("Failed running Fisher's exact test even with Monte Carlo simulation. A default P-value of NA will be printed - please inspect data.")
                return(NA)
              })
            })
            tmp$p_val = p_value
            colnames(tmp) <- c("category","MELD_less20","MELD_20ormore","p_val")
            demog_meld_obf <- rbind(demog_meld_obf,tmp)
            rm(tmp)
            rm(p_value)
          })
        }
        cat("\nComputing admission labs averages\n")
        #Calculate stats of admission labs
        lab_meld_list <- NULL
        lab_anova <- NULL
        lab_meld_stats <- NULL
        lab_anova_raw <- NULL
        tryCatch({
          try({
            cat("Two-way ANOVA for AST (MELD >= 20, AKI as variables)...\n")
            ast_stats <- demog_meld_summ %>% dplyr::select(patient_id,first_ast,meld_admit_severe,aki)
            model <- lm(first_ast ~ meld_admit_severe * aki,data=ast_stats)
            ast_anova_raw <- car::Anova(model,type=3)
            ast_anova <- broom::tidy(ast_anova_raw)
            lab_anova_raw <- c(lab_anova_raw,"ast_anova_raw")
            lab_meld_list <- c(lab_meld_list,"first_ast")
            lab_anova <- c(lab_anova,"ast_anova")
            ast_stats <- ast_stats %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
              mean_admit_ast = mean(first_ast,na.rm=TRUE),
              sd_admit_ast = sd(first_ast,na.rm=TRUE),
              n_admit_ast = sum(!is.na(first_ast))
            ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
            lab_meld_stats <- c(lab_meld_stats,"ast_stats")
          })
          try({
            cat("Two-way ANOVA for ALT (MELD >= 20, AKI as variables)...\n")
            alt_stats <- demog_meld_summ %>% dplyr::select(patient_id,first_alt,meld_admit_severe,aki)
            model <- lm(first_alt ~ meld_admit_severe * aki,data=alt_stats)
            alt_anova_raw <- car::Anova(model,type=3)
            alt_anova <- broom::tidy(alt_anova_raw)
            lab_anova_raw <- c(lab_anova_raw,"alt_anova_raw")
            lab_meld_list <- c(lab_meld_list,"first_alt")
            lab_anova <- c(lab_anova,"alt_anova")
            alt_stats <- alt_stats %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
              mean_admit_alt = mean(first_alt,na.rm=TRUE),
              sd_admit_alt = sd(first_alt,na.rm=TRUE),
              n_admit_alt = sum(!is.na(first_alt))
            ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
            lab_meld_stats <- c(lab_meld_stats,"alt_stats")
          })
          try({
            cat("Two-way ANOVA for bilirubin (MELD >= 20, AKI as variables)...\n")
            bil_stats <- demog_meld_summ %>% dplyr::select(patient_id,max_bil,meld_admit_severe,aki)
            model <- lm(max_bil ~ meld_admit_severe * aki,data=bil_stats)
            bil_anova_raw <- car::Anova(model,type=3)
            bil_anova <- broom::tidy(bil_anova_raw)
            lab_anova_raw <- c(lab_anova_raw,"bil_anova_raw")
            lab_meld_list <- c(lab_meld_list,"max_bil")
            lab_anova <- c(lab_anova,"bil_anova")
            bil_stats <- bil_stats %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
              mean_admit_bil = mean(max_bil,na.rm=TRUE),
              sd_admit_bil = sd(max_bil,na.rm=TRUE),
              n_admit_bil = sum(!is.na(max_bil))
            ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
            lab_meld_stats <- c(lab_meld_stats,"bil_stats")
          })
          try({
            cat("Two-way ANOVA for INR (MELD >= 20, AKI as variables)...\n")
            inr_stats <- demog_meld_summ %>% dplyr::select(patient_id,max_inr,meld_admit_severe,aki)
            model <- lm(max_inr ~ meld_admit_severe * aki,data=inr_stats)
            inr_anova_raw <- car::Anova(model,type=3)
            inr_anova <- broom::tidy(inr_anova_raw)
            lab_anova_raw <- c(lab_anova_raw,"inr_anova_raw")
            lab_meld_list <- c(lab_meld_list,"max_inr")
            lab_anova <- c(lab_anova,"inr_anova")
            inr_stats <- inr_stats %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
              mean_admit_inr = mean(max_inr,na.rm=TRUE),
              sd_admit_inr = sd(max_inr,na.rm=TRUE),
              n_admit_inr = sum(!is.na(max_inr))
            ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
            lab_meld_stats <- c(lab_meld_stats,"inr_stats")
          })
          try({
            cat("Two-way ANOVA for alb (MELD >= 20, AKI as variables)...\n")
            alb_stats <- demog_meld_summ %>% dplyr::select(patient_id,first_alb,meld_admit_severe,aki)
            model <- lm(first_alb ~ meld_admit_severe * aki,data=alb_stats)
            alb_anova_raw <- car::Anova(model,type=3)
            alb_anova <- broom::tidy(alb_anova_raw)
            lab_anova_raw <- c(lab_anova_raw,"alb_anova_raw")
            lab_meld_list <- c(lab_meld_list,"first_alb")
            lab_anova <- c(lab_anova,"alb_anova")
            alb_stats <- alb_stats %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
              mean_admit_alb = mean(first_alb,na.rm=TRUE),
              sd_admit_alb = sd(first_alb,na.rm=TRUE),
              n_admit_alb = sum(!is.na(first_alb))
            ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
            lab_meld_stats <- c(lab_meld_stats,"alb_stats")
          })
          try({
            cat("Two-way ANOVA for Na (MELD >= 20, AKI as variables)...\n")
            na_stats <- demog_meld_summ %>% dplyr::select(patient_id,min_na,meld_admit_severe,aki)
            model <- lm(min_na ~ meld_admit_severe * aki,data=ast_stats)
            na_anova_raw <- car::Anova(model,type=3)
            na_anova <- broom::tidy(na_anova_raw)
            lab_anova_raw <- c(lab_anova_raw,"na_anova_raw")
            lab_meld_list <- c(lab_meld_list,"min_na")
            lab_anova <- c(lab_anova,"na_anova")
            na_stats <- na_stats %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
              mean_admit_na = mean(first_na,na.rm=TRUE),
              sd_admit_na = sd(first_na,na.rm=TRUE),
              n_admit_na = sum(!is.na(first_na))
            ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
            lab_meld_stats <- c(lab_meld_stats,"na_stats")
          })
          lab_meld_stats <- mget(lab_meld_stats) %>% purrr::reduce(dplyr::full_join,by=c("meld_admit_severe","aki"))
        }, error = function(e) {
          message("\nError in processing labs. Check error messages.\n")
          message("\nOriginal error:\n",e)
          cat("\nError in processing labs. Check error messages.\n")
          cat("\nOriginal error:\n",e,"\n")
        })
        try({
          demog_meld_obf$MELD_less20[demog_meld_obf$MELD_less20 < obfuscation_value] <- 0
          demog_meld_obf$MELD_20ormore[demog_meld_obf$MELD_20ormore < obfuscation_value] <- 0
          demog_meld_obf <- demog_meld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(total = MELD_less20 + MELD_20ormore) %>% dplyr::ungroup()
          meld_less20_total <- demog_meld_obf$MELD_less20[1]
          meld_20ormore_total <- demog_meld_obf$MELD_20ormore[1]
          total_pop <- demog_meld_obf$total[1]
          demog_meld_obf <- demog_meld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(MELD_less20_perc = MELD_less20 / meld_less20_total * 100, MELD_20ormore_perc = MELD_20ormore/meld_20ormore_total * 100, total_perc = total/total_pop * 100) %>% dplyr::ungroup()
          cat("\nComputed final demog_meld_obf table.\n")
        })
        
        if(obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
          table_one_meld <- tableone::CreateTableOne(data=demog_meld_summ,vars=table_one_meld_vars,strata="meld_admit_severe")
          export_table_one_meld <- print(table_one_meld,showAllLevels=TRUE,formatOptions=list(big.mark=","))
          if(exists("demog_meld_obf")) {
            demog_meld_files <- c(demog_meld_files,"demog_meld_obf","table_one_meld","export_table_one_meld","lab_meld_stats",lab_anova,lab_anova_raw)
          } else {
            demog_meld_files <- c(demog_meld_files,"table_one_meld","export_table_one_meld","lab_meld_stats",lab_anova,lab_anova_raw)
          }
        } else {
          demog_meld_files <- c(demog_meld_files,"demog_meld_obf","lab_meld_stats",lab_anova,lab_anova_raw)
        }
      }
    })
    
  }
  
  if(obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
    table_one <- tableone::CreateTableOne(data=demog_summ,vars=table_one_vars,strata="aki")
    export_table_one <- print(table_one,showAllLevels=TRUE,formatOptions=list(big.mark=","))
    write.csv(export_table_one,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne",file_suffix)))
    try({
      if(isTRUE(cirrhosis_present)) {
        table_one_cld <- tableone::CreateTableOne(data=demog_cld_summ,vars=table_one_vars,strata="aki")
        export_table_one_cld <- print(table_one_cld,showAllLevels=TRUE,formatOptions=list(big.mark=","))
      }
    })
  }
  cat("\nAttempting to save demographics tables for cirrhotic patients.\n")
  try({
    if(isTRUE(cirrhosis_present)) {
      if(obfuscation_value == 0 & isTRUE(is_obfuscated)) {
        demog_meld_files <- c(demog_meld_files,"table_one_cld","export_table_one_cld","demog_cld_obf")
      } else if (obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
        demog_meld_files <- c(demog_meld_files,"table_one_cld","export_table_one_cld")
      } else if (isTRUE(is_obfuscated)) {
        demog_meld_files <- c(demog_meld_files,"demog_cld_obf")
      }
    }
  })
  if(!is.null(demog_meld_files)) {
    file_rda_suffix <- ".rdata"
    if(isTRUE(preadmit_only_analysis)) {
      file_rda_suffix <- "_preadmit_cr_only.rdata"
    }
    try({
      save(list=demog_meld_files,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_MELD_Cirrhosis",file_rda_suffix)),compress="bzip2")
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