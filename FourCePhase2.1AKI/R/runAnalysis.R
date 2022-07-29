
#' Runs the analytic workflow for the AKI project
#'
#' @keywords 4CE
#' @export

runAnalysis <- function(is_obfuscated=TRUE,factor_cutoff = 5, ckd_cutoff = 2.25, restrict_models = FALSE, docker = TRUE, input = "/4ceData/Input", siteid_nodocker = "", skip_qc = FALSE, offline = FALSE,allow_no_prior_cr = FALSE, use_rrt_surrogate = TRUE,print_rrt_surrogate = FALSE,modified_severity = FALSE,debug_on=FALSE,date_cutoff = "2020-09-10",lab_date_cutoff = "2021-09-10",custom_output = FALSE, custom_output_dir = "/4ceData/Output") {
    
    if(isFALSE(offline)) {
        ## make sure this instance has the latest version of the quality control and data wrangling code available
        devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)
    }
    
    ## ========================================
    ## PART 1: Read in Data Tables
    ## ========================================
    cat("\nPlease ensure that your working directory is set to /4ceData")
    cat("\nReading in Input/LocalPatientSummary.csv and Input/LocalPatientObservations.csv...\n")
    
    if(isTRUE(custom_output)) {
        dir.output <- custom_output_dir
    } else {
        dir.output <- getProjectOutputDirectory()
    }
    
    if(isTRUE(docker)) {
        ## get the site identifier associated with the files stored in the /4ceData/Input directory that 
        ## is mounted to the container
        currSiteId = toupper(FourCePhase2.1Data::getSiteId())
        ## run the quality control
        FourCePhase2.1Data::runQC(currSiteId)
        demographics <- FourCePhase2.1Data::getLocalPatientSummary(currSiteId)
        observations <- FourCePhase2.1Data::getLocalPatientObservations(currSiteId)
        course <- FourCePhase2.1Data::getLocalPatientClinicalCourse(currSiteId)
    } else {
        currSiteId = siteid_nodocker
        if(isFALSE(skip_qc)) {
            FourCePhase2.1Data::runQC_nodocker(currSiteId,input)
        }
        demographics <- FourCePhase2.1Data::getLocalPatientSummary_nodocker(currSiteId,input)
        observations <- FourCePhase2.1Data::getLocalPatientObservations_nodocker(currSiteId,input)
        course <- FourCePhase2.1Data::getLocalPatientClinicalCourse_nodocker(currSiteId,input)
    }
    
    error_log <- file(file.path(dir.output,paste0("/",currSiteId,"_error.log")))
    sink(error_log,append=TRUE,split=TRUE)
    if(isTRUE(debug_on)) {
        cat("\n===================\nYou have enabled debugging for this session. Warnings and messages will be redirected to the error log.\nDue to the nature of error handling in R, this output will not be displayed on the console.\nOutput via the cat() or print() functions will still be visible.\n===================\n")
        sink(error_log,append=TRUE,type="message")
    }
    
    cat("\n4CE AKI Analysis\n======================\nVersion",paste0(packageVersion("FourCePhase2.1AKI")),"\n\n")
    
    # Detects if there are issues with site data extraction being inadequate for study.
    try({
        data_date_offset <- demographics %>% dplyr::group_by(patient_num) %>% dplyr::mutate(offset = as.Date(admission_date) + days_since_admission)
        data_date_offset <- unlist(min(data_date_offset$offset))
        lab_date_cutoff <- as.Date(lab_date_cutoff)
        if(data_date_offset > lab_date_cutoff) {
            gap <- as.numeric(data_date_offset - lab_date_cutoff)
            cat("Data exceeds required date of ",as.character(lab_date_cutoff)," by ",gap," days.\n")
            demographics <- demographics %>% dplyr::group_by(patient_num) %>% dplyr::mutate(days_since_admission = days_since_admission - gap)
            date_offset <- demographics[,c("patient_num","days_since_admission")]
            colnames(date_offset)[2] <- "max_days"
            observations <- merge(observations,date_offset,by="patient_num",all.x=TRUE) %>% dplyr::group_by(patient_num) %>% dplyr::filter(days_since_admission <= max_days) %>% dplyr::ungroup()
            course <- course %>% dplyr::filter(calendar_date <= lab_date_cutoff)
            invisible(gc())
        } else if(data_date_offset == lab_date_cutoff) {
            cat("Date of extraction matches ",lab_date_cutoff, " exactly.\n")
        } else {
            message("Warning: please see below (INCOMPLETE DATA)\n")
            cat("Warning: time-to-event data may be inaccurate as last lab values were before ",lab_date_cutoff," as extraction date was ",data_date_offset,".\nPlease ensure that the correct duration of patient data was extracted - modify and re-run any custom SQL scripts to fit the date criteria.\n")
        }
    })
    
    obfuscation_value = as.numeric(FourCePhase2.1Data::getObfuscation(currSiteId))
    cat(paste0(c("\nObfuscation level set to ",obfuscation_value)))

    cat("\nTransforming the Summary and Observations tables to generate tables for demographics, diagnoses, procedures.")
    # first generate a unique ID for each patient
    demographics <- demographics %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    observations <- observations %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    course <- course %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    
    # From this point on, we will be using our custom-generated patient_id as a unique patient identifier
    # Reorder the columns in each table to bring patient_id to the first column and remove patient_num
    demographics <- demographics %>% dplyr::select(patient_id,siteid,admission_date,days_since_admission,last_discharge_date,still_in_hospital,severe_date,severe,death_date,deceased,sex,age_group,race,race_collected)
    observations <- observations %>% dplyr::select(patient_id,siteid,days_since_admission,concept_type,concept_code,value)
    demographics$death_date <- as.Date(demographics$death_date,format="%Y-%m-%d")
    demographics$admission_date <- as.Date(demographics$admission_date,format="%Y-%m-%d")
    
    # Identify patients without any age or < 18 years and remove them from analysis
    young_patients <- unlist(demographics$patient_id[demographics$age_group %in% c("other", "00to02", "03to05", "06to11", "12to17")])
    if(length(young_patients) > 0) {
        demographics <- demographics[!(demographics$patient_id %in% young_patients),]
        observations <- observations[!(observations$patient_id %in% young_patients),]
        course <- course[!(course$patient_id %in% young_patients),]
    }
    
    # Plot histogram of admission dates for site, binned by month
    cat("\nCreating histogram of admission dates per month.\n")
    admission_date_histogram <- ggplot2::ggplot(demographics,ggplot2::aes(x=admission_date)) + ggplot2::stat_bin(binwidth = 30,position = "identity") + ggplot2::scale_x_date(date_breaks = "months", date_labels = "%b-%Y", guide=ggplot2::guide_axis(angle=60)) + ggplot2::theme_bw() + ggplot2::xlab("Admission Month") + ggplot2::ylab("No. of Admissions")
    ggplot2::ggsave(filename=file.path(dir.output,paste0(currSiteId, "_AdmissionDates_Histogram_Waves.png")),plot=print(admission_date_histogram))
    
    # Remove patients who were admitted after Dec 31, 2020 (or other pre-specific date cutoff)
    later_patients <- unlist(demographics$patient_id[demographics$admission_date > as.Date(date_cutoff)])
    cat("\nNo. of patients with admission dates past ",date_cutoff,": ",length(later_patients),"\n")
    if(length(later_patients) > 0) {
        demographics <- demographics[!(demographics$patient_id %in% later_patients),]
        observations <- observations[!(observations$patient_id %in% later_patients),]
        course <- course[!(course$patient_id %in% later_patients),]
    }
    
    # Generate a diagnosis table
    diagnosis <- observations[observations$concept_type %in% c("DIAG-ICD9","DIAG-ICD10"),-6]
    colnames(diagnosis) <- c("patient_id","siteid","days_since_admission","concept_type","icd_code")
    diag_icd9 <- diagnosis[diagnosis$concept_type == "DIAG-ICD9",]
    diag_icd10 <- diagnosis[diagnosis$concept_type == "DIAG-ICD10",]
    
    # Generate a procedures table
    procedures <- observations[observations$concept_type %in% c("PROC-ICD9", "PROC-ICD10"),-6]
    colnames(procedures) <- c("patient_id","siteid","days_since_admission","icd_version","procedure_code")
    
    demographics_filt <- demographics %>% dplyr::mutate(time_to_severe = ifelse(severe == 1, as.numeric(as.Date(severe_date) - as.Date(admission_date)),NA))
    demographics_filt <- demographics_filt %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_to_death = ifelse(deceased == 1, as.numeric(death_date - admission_date),NA))
    demographics_filt <- demographics_filt %>% dplyr::mutate(length_stay = ifelse(still_in_hospital==1,days_since_admission,as.numeric(as.Date(last_discharge_date) - as.Date(admission_date))))
    
    # Reorder the columns to be more readable
    demographics_filt <- demographics_filt %>% dplyr::select(patient_id,siteid,sex,age_group,race,length_stay,severe,time_to_severe,deceased,time_to_death)
    
    # Final headers for demographics_filt
    # patient_id  siteid sex age_group race  length_stay severe  time_to_severe  deceased  time_to_death
    
    # Since AKIs can occur in later admissions but not in the index admission (especially when expanding)
    # the window for the baseline Cr, we will need to first extract the first discharge date/day for each
    # patient, and then use this as a basis to filter patients
    # We cannot naively cut off Cr values as their true baseline may be achieved even after first discharge
    
    first_discharge <- course %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(in_hospital == 0) %>% dplyr::filter(days_since_admission == min(days_since_admission,na.rm=TRUE))
    first_discharge <- first_discharge %>% dplyr::select(patient_id,days_since_admission)
    colnames(first_discharge) <- c("patient_id","first_discharge_day")
    
    
    # Now we need to fix a bug in the labeling of severe patients
    # The problem with the original approach of labeling a patient as severe (0,1) as per Griffin's SQL scripts is that
    # the presence of any of the markers beyond the first admission would automatically confer severe status to a patient
    # even if the corresponding markers appear at a significant period of time post-COVID-19 infection (e.g. patient 
    # admitted for appendicitis requiring ICU admission a year after a clinically mild COVID-19 infection).
    # The following lines attempt to correct the labeling using the severe dates instead as a marker
    demographics_filt <- merge(demographics_filt,first_discharge,by="patient_id",all.x=T)
    demographics_filt <- demographics_filt %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_not_on_first_admit = dplyr::if_else( !is.na(time_to_severe)| (time_to_severe > first_discharge_day),1,0)) %>% dplyr::ungroup()
    cat("No. of patients with severe markers past index admission (demographics_filt):",sum(demographics_filt$severe_not_on_first_admit,na.rm=T),"\n")
    demographics_filt <- demographics_filt %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = dplyr::if_else(severe == 0,0,dplyr::if_else(time_to_severe > first_discharge_day,0,1))) %>% dplyr::ungroup()
    demographics_filt <- demographics_filt %>% dplyr::select(patient_id,siteid,sex,age_group,race,length_stay,severe,time_to_severe,deceased,time_to_death)
    
    demographics <- merge(demographics[,c("patient_id","siteid","admission_date","days_since_admission","last_discharge_date","still_in_hospital","severe_date","death_date","deceased","sex","age_group","race","race_collected")],demographics_filt[,c("patient_id","severe")],by="patient_id",all.x=T)
    demographics <- demographics %>% dplyr::select(patient_id,siteid,admission_date,days_since_admission,last_discharge_date,still_in_hospital,severe_date,severe,death_date,deceased,sex,age_group,race,race_collected)
    
    invisible(gc())
    # Time to RRT
    # ==================
    cat("\nCreating table for RRT/Kidney transplant...")
    rrt_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","3E1.M39Z","54.98","39.95","55.61","55.69","0TY00Z0","0TY00Z1","0TY00Z2","0TY10Z0","0TY10Z1","0TY10Z2")
    rrt_diag_icd10_code <- c("T82","Y60","Y61","Y62","Y84","Z49","Z99")
    rrt_diag_icd9_code <- c("458","792","V45","V56")
    #hd_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","399.5")
    #pd_code <- c("3E1.M39Z","549.8")
    
    rrt <- procedures[procedures$procedure_code %in% rrt_code,-c(2,4)]
    rrt <- rrt[,c(1,2)]
    # rrt <- rrt[order(rrt$patient_id,rrt$days_since_admission),]
    # rrt <- rrt[!duplicated(rrt$patient_id),]
    
    rrt_icd9 <- diag_icd9[diag_icd9$icd_code %in% rrt_diag_icd9_code,c("patient_id","days_since_admission")]
    rrt_icd10 <- diag_icd10[diag_icd10$icd_code %in% rrt_diag_icd9_code,c("patient_id","days_since_admission")]
    rrt_diag <- dplyr::bind_rows(rrt_icd9,rrt_icd10)
    rrt <- dplyr::bind_rows(rrt,rrt_diag) %>% dplyr::arrange(patient_id,days_since_admission) %>% dplyr::distinct(patient_id,.keep_all = TRUE)
    
    # Generate list of patients already on RRT prior to admission
    # This list can be used to exclude ESRF patients in subsequent analyses
    rrt_old <- unlist(rrt$patient_id[rrt$days_since_admission < 0])
    
    # Generate list of patients who were only initiated on RRT for the first time ever during admission for COVID-19
    rrt_new <- unlist(rrt$patient_id[!(rrt$patient_id %in% rrt_old)])
    
    rrt_index_admit <- merge(rrt[!(rrt$patient_id %in% rrt_old),],first_discharge,by="patient_id",all.x=TRUE)
    rrt_index_admit <- rrt_index_admit[!is.na(rrt_index_admit$first_discharge_day),] %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= first_discharge_day) %>% dplyr::ungroup()
    rrt_index_admit <- unlist(rrt_index_admit$patient_id)
    
    rrt_future_admit <- setdiff(rrt_new,rrt_index_admit)
    
    # For debugging purposes
    cat(paste("\nNumber of patients on RRT in total: ",length(rrt$patient_id),"\nRRT previously: ",length(rrt_old),"\nFirst RRT during and/or after first admission: ",length(rrt_new),"\nFirst RRT during index admission: ",length(rrt_index_admit),"\nFirst RRT after index admission:",length(rrt_future_admit)))
    
    # =============
    # Medications
    # =============
    cat("\nNow generating the medications tables:")
    medications <- observations[observations$concept_type == "MED-CLASS",]
    medications <- medications[,-c(4,6)]
    medications <- medications %>% dplyr::arrange(patient_id,days_since_admission,concept_code)
    # Use 15 days as the cutoff for chronic medications
    cat("\nGenerating table for chronic medications...")
    med_chronic <- medications[medications$days_since_admission < -15,]
    med_new <- medications[medications$days_since_admission >= -15,]
    
    # Re-code chronic medications into wide format
    med_chronic <- med_chronic[!duplicated(med_chronic[,c(1,2,4)]),]
    med_chronic <- med_chronic[,-c(2,3)]
    med_chronic$concept_code <- paste("old",med_chronic$concept_code,sep="_")
    med_chronic$present <- 1
    med_chronic <- med_chronic %>% tidyr::spread(concept_code,present)
    med_chronic[is.na(med_chronic)] <- 0
    
    # Create subtable for ACE-i/ARB pre-exposure
    cat("\nGenerating table for prior ACE-inhibitor (ACEI) or angiotensin receptor blocker (ARB) use...")
    acei_present = ("old_ACEI" %in% colnames(med_chronic))
    arb_present = ("old_ARB" %in% colnames(med_chronic))
    med_acearb_chronic <- NULL
    if(acei_present == TRUE && arb_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ACEI,old_ARB)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ACEI + old_ARB > 0,1,0)) %>% dplyr::ungroup()
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    } else if (acei_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ACEI)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ACEI > 0,1,0)) %>% dplyr::ungroup()
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    } else if (arb_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ARB)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ARB > 0,1,0)) %>% dplyr::ungroup()
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    }
    
    # med_acearb_chronic
    # Headers: patient_id   acei_arb_preexposure
    
    # For simplicity of initial analysis, we will use the earliest date where each new medication class is used.
    cat("\nGenerating table for new medications started during the admission...")
    med_new <- med_new[!duplicated(med_new[,c(1,2,4)]),]
    med_new <- med_new[,c(1,4,3)]
    colnames(med_new)[3] <- "start_day"
    
    # Generate another table with the start date of the new medications
    med_new <- med_new %>% tidyr::spread(concept_code,start_day)
    #med_new[is.na(med_new)] <- -999
    
    # Generate simplified table for determining who were started on COAGA near admission
    cat("\nGenerating table for COAGA use during the admission...")
    med_coaga_new <- NULL
    coaga_present <- tryCatch({
        med_coaga_new <- med_new %>% dplyr::select(patient_id,COAGA) %>% dplyr::group_by(patient_id) %>% dplyr::filter(COAGA == min(COAGA,na.rm=TRUE)) %>% dplyr::ungroup() %>% dplyr::distinct()
        med_coaga_new$COAGA[med_coaga_new$COAGA < -15] <- 0
        med_coaga_new$COAGA[med_coaga_new$COAGA >= -15] <- 1
        TRUE
    },error= function(c) FALSE)
    
    # Generate simplified table for determining who were started on COAGB near admission
    cat("\nGenerating table for COAGB use during the admission...")
    med_coagb_new <- NULL
    coagb_present <- tryCatch({
        med_coagb_new <- med_new %>% dplyr::select(patient_id,COAGB) %>% dplyr::group_by(patient_id) %>% dplyr::filter(COAGB == min(COAGB,na.rm=TRUE)) %>% dplyr::ungroup() %>% dplyr::distinct()
        med_coagb_new$COAGB[med_coagb_new$COAGB < -15] <- 0
        med_coagb_new$COAGB[med_coagb_new$COAGB >= -15] <- 1
        TRUE
    },error=function(c) FALSE)
    
    # Generate simplified table for determining who were started on novel antivirals
    covid19antiviral_present = ("COVIDVIRAL" %in% colnames(med_new))
    remdesivir_present = ("REMDESIVIR" %in% colnames(med_new))
    med_covid19_new <- NULL
    if(isTRUE(covid19antiviral_present) && isTRUE(remdesivir_present)){
        cat("\nGenerating table for experimental COVID-19 treatment... (both COVIDVIRAL and REMDESIVIR)")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,COVIDVIRAL,REMDESIVIR) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covidviral=ifelse(is.na(COVIDVIRAL),0,COVIDVIRAL),remdesivir=ifelse(is.na(REMDESIVIR),0,REMDESIVIR)) %>% dplyr::ungroup()
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(covidviral > 0 | remdesivir > 0,1,0)) %>% dplyr::ungroup()
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx,covidviral,remdesivir)
    } else if(isTRUE(covid19antiviral_present)){
        cat("\nGenerating table for experimental COVID-19 treatment... (COVIDVIRAL only)")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,COVIDVIRAL) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covidviral=ifelse(is.na(COVIDVIRAL),0,COVIDVIRAL)) %>% dplyr::ungroup()
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(covidviral > 0,1,0)) %>% dplyr::ungroup()
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx,covidviral)
    } else if(isTRUE(remdesivir_present)){
        cat("\nGenerating table for experimental COVID-19 treatment... (REMDESIVIR only)")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,REMDESIVIR) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(remdesivir=ifelse(is.na(REMDESIVIR),0,REMDESIVIR)) %>% dplyr::ungroup()
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(remdesivir > 0,1,0)) %>% dplyr::ungroup()
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx,remdesivir)
    }
    
    invisible(gc())
    
    # ===============================================================
    # Lancet Reviewer edit: Create new definition of severe COVID-19 excluding ABGs
    # ===============================================================
    if(isTRUE(modified_severity)) {
      cat("\nGenerating the modified severity definition as requested by Lancet reviewers...\nProcedures...\n")
      procedures_severe <- procedures[procedures$procedure_code %in% c("0BH17EZ","0BH","5A0","5A093","5A0935","5A09357","5A09358","5A09359","5A0935A","5A0935B","5A0935Z","5A094","5A0945","5A09457","5A09458","5A09459","5A0945A","5A0945B","5A0945Z","5A095","5A0955","5A09557","5A09558","5A09559","5A0955A","5A0955B","5A0955Z","96","096","09604","9604","96.04","096.04","96.70","9670","096.70","09670","96.71","9671","096.71","09671","96.72","9672","096.72","09672"),]
      procedures_severe <- merge(procedures_severe,first_discharge,by="patient_id",all.x=T)
      procedures_severe <- procedures_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= first_discharge_day & days_since_admission >= 0) %>% dplyr::mutate(procedure_severe_status = 1) %>% dplyr::ungroup() %>% dplyr::select(patient_id,procedure_severe_status) %>% dplyr::distinct()
      if(nrow(procedures_severe) == 0) {
        procedures_severe <- NULL
        cat("\nEmpty Severe Procedures List.")
      }
      
      cat("\nMedications...")
      severe_meds_list <- c("SIANES","SICARDIAC")
      med_new_list <- colnames(med_new)
      severe_meds_present_list <- med_new_list[med_new_list %in% severe_meds_list]
      if(length(severe_meds_present_list) > 0) {
        severe_meds_new <- merge(med_new[,c("patient_id",severe_meds_present_list)],first_discharge,by="patient_id",all.x=T)
        if(("SIANES" %in% colnames(severe_meds_new)) & ("SICARDIAC" %in% colnames(severe_meds_new))) {
          severe_meds_new <- severe_meds_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(sianes_severe = ifelse(is.na(SIANES) | SIANES > first_discharge_day,0,1),sicardiac_severe = ifelse(is.na(SICARDIAC) | SICARDIAC > first_discharge_day,0,1)) %>% dplyr::mutate(meds_severe_status = ifelse(sianes_severe + sicardiac_severe > 0,1,0)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,meds_severe_status)
        } else if("SIANES" %in% colnames(severe_meds_new)) {
          severe_meds_new <- severe_meds_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meds_severe_status = ifelse(is.na(SIANES) | SIANES > first_discharge_day,0,1)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,meds_severe_status)
        } else if("SICARDIAC" %in% colnames(severe_meds_new)) {
          severe_meds_new <- severe_meds_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meds_severe_status = ifelse(is.na(SICARDIAC) | SICARDIAC > first_discharge_day,0,1)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,meds_severe_status)
        } else {
          severe_meds_new <- NULL
          cat("\nEmpty severe medications list.")
        }
      } else {
        severe_meds_new <- NULL
        cat("\nEmpty severe medications list.")
      }
      
      cat("\nDiagnoses...")
      diagnosis_severe <- diagnosis[diagnosis$icd_code %in% c("J80","518","J95","997","51882","518.82","J95.851","J95851","997.31","99731"),c("patient_id","days_since_admission")]
      if(nrow(diagnosis_severe) > 0) {
        diagnosis_severe <- merge(diagnosis_severe,first_discharge,by="patient_id",all.x=T)
        diagnosis_severe <- diagnosis_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0 & days_since_admission <= first_discharge_day) %>% dplyr::mutate(diagnosis_severe_status = 1) %>% dplyr::ungroup() %>% dplyr::select(patient_id,diagnosis_severe_status) %>% dplyr::distinct()
        if(nrow(diagnosis_severe) == 0) {
          diagnosis_severe <- NULL
          cat("\nEmpty Severe Procedures List.")
        }
      } else {
        diagnosis_severe <- NULL
        cat("\nEmpty Severe Procedures List.")
      }
      
      cat("\nMerging tables...")
      severe_lancet_table <- NULL
      if(is.null(procedures_severe) & is.null(severe_meds_new) & is.null(diagnosis_severe)) {
        cat("\nAll lists are empty. Will proceed with rest of analysis as usual, with the ORIGINAL severity definition. This means severity is predominated by ABGs.")
      } else {
        if(is.null(procedures_severe) & is.null(severe_meds_new)) {
          severe_lancet_table <- diagnosis_severe
          colnames(severe_lancet_table)[2] <- "severe"
        } else if(is.null(procedures_severe) & is.null(diagnosis_severe)) {
          severe_lancet_table <- severe_meds_new
          colnames(severe_lancet_table)[2] <- "severe"
        } else if(is.null(severe_meds_new) & is.null(diagnosis_severe)) {
          severe_lancet_table <- procedures_severe
          colnames(severe_lancet_table)[2] <- "severe"
        } else if(is.null(procedures_severe)) {
          severe_lancet_table <- merge(severe_meds_new,diagnosis_severe,by="patient_id",all=T)
          severe_lancet_table[is.na(severe_lancet_table)] <- 0
          severe_lancet_table <- severe_lancet_table %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(meds_severe_status + diagnosis_severe_status > 0, 1,0)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,severe) %>% dplyr::distinct()
        } else if(is.null(severe_meds_new)) {
          severe_lancet_table <- merge(procedures_severe,diagnosis_severe,by="patient_id",all=T)
          severe_lancet_table[is.na(severe_lancet_table)] <- 0
          severe_lancet_table <- severe_lancet_table %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(procedures_severe_status + diagnosis_severe_status > 0, 1,0)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,severe) %>% dplyr::distinct()
        } else if(is.null(diagnosis_severe)) {
          severe_lancet_table <- merge(procedures_severe,severe_meds_new,by="patient_id",all=T)
          severe_lancet_table[is.na(severe_lancet_table)] <- 0
          severe_lancet_table <- severe_lancet_table %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(procedures_severe_status + meds_severe_status > 0, 1,0)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,severe) %>% dplyr::distinct()
        } else {
          severe_lancet_table <- merge(procedures_severe,severe_meds_new,by="patient_id",all=T)
          severe_lancet_table <- merge(severe_lancet_table,diagnosis_severe,by="patient_id",all=T)
          severe_lancet_table[is.na(severe_lancet_table)] <- 0
          severe_lancet_table <- severe_lancet_table %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(procedures_severe_status + meds_severe_status + diagnosis_severe_status > 0, 1,0)) %>% dplyr::ungroup() %>% dplyr::select(patient_id,severe) %>% dplyr::distinct()
        }
        
        demographics_filt <- demographics_filt %>% dplyr::select(patient_id,siteid,sex,age_group,race,length_stay,time_to_severe,deceased,time_to_death)
        demographics_filt <- merge(demographics_filt,severe_lancet_table,by="patient_id",all.x=T)
        demographics_filt$severe[is.na(demographics_filt$severe)] <- 0
        demographics_filt <- demographics_filt %>% dplyr::select(patient_id,siteid,sex,age_group,race,length_stay,time_to_severe,deceased,time_to_death)
      }
      invisible(gc())
      
    }
    
    # ====================
    # PART 2: AKI Detection Code
    # ====================
    # For the purposes of generating a table of AKI events, we will have to create a new data.table
    # with the serum creatinine levels.
    # We will then use this table, labs_cr_aki, to generate a summary table containing details of 
    # each AKI event and the post-AKI recovery labs
    cat("\nNow proceeding to AKI detection code:")
    
    if(!is.null(rrt_old)) {
        cat("\nFirst removing patients who were previously on RRT (based on procedure codes)...")
        demographics_filt <- demographics_filt[!(demographics_filt$patient_id %in% rrt_old),]
        observations <- observations[!(observations$patient_id %in% rrt_old),]
        course <- course[!(course$patient_id %in% rrt_old),]
    }
    
    cat("\nExtracting serum creatinine values...")
    # Extract serum Cr levels
    labs_cr_aki <- observations[observations$concept_code == '2160-0',] #LOINC code for Cr 2160-0
    # Remove unnecessary columns
    labs_cr_aki <- labs_cr_aki[,-c(4,5)]
    # Filter for labs >= -365 days
    labs_cr_aki <- labs_cr_aki %>% dplyr::filter(days_since_admission >= -365)
    cat("\nInitial extraction - total number of patients: ",length(unique(labs_cr_aki$patient_id)))
    
    # Generate separate demographics table for patients who do not have any sCr values fulfilling 
    # the above (e.g. all the labs are before -365 days or patient has no sCr value in index admission)
    cat("\nRemoving patients who do not have any serum creatinine values during admission...")
    pts_valid_cr <- labs_cr_aki %>% dplyr::filter(days_since_admission >= 0)
    pts_valid_cr <- merge(pts_valid_cr,first_discharge,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    pts_valid_cr <- pts_valid_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= first_discharge_day) %>% dplyr::ungroup()
    pts_valid_cr <- unique(pts_valid_cr$patient_id)
    cat("\nPatients with at least 1 sCr value in first admission (before subsetting tables): ",length(pts_valid_cr))
    cat("\nCurrent total number of patients in demographics_filt prior to filtering: ",length(demographics_filt$patient_id))
    cat("\nNo. of patients without any serum creatinine during first admission: ",length(demographics_filt$patient_id)-length(pts_valid_cr))
    demog_no_cr <- demographics_filt[!(demographics_filt$patient_id %in% pts_valid_cr),]
    demographics_filt <- demographics_filt[demographics_filt$patient_id %in% pts_valid_cr,]
    labs_cr_aki <- labs_cr_aki[labs_cr_aki$patient_id %in% pts_valid_cr,]
    cat("\nPatients with at least 1 sCr value during first admission: ",length(unique(labs_cr_aki$patient_id)))
    
    labs_cr_aki <- labs_cr_aki %>% dplyr::arrange(patient_id,days_since_admission) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(rec = 1) %>% dplyr::mutate(count_prior_cr = cumsum(rec)) %>% dplyr::ungroup()
    labs_cr_aki <- labs_cr_aki[,-5]
    
    # Find patients who have only 1 serum creatinine reading by the end of first admission and remove these patients, i.e. this single serum creatinine value is only occuring during the admission itself and not prior to the admission
    # Lancet reviewer edit: to provide a flag that allows this removal to be disabled
    if(isFALSE(allow_no_prior_cr)) {
      pts_insufficient_cr <- merge(labs_cr_aki,first_discharge,by="patient_id",all.x=TRUE) %>% dplyr::distinct() %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= first_discharge_day) %>% dplyr::filter(max(count_prior_cr) <= 1) %>% dplyr::ungroup()
      pts_insufficient_cr <- unlist(unique(pts_insufficient_cr$patient_id))
      cat("\nNo. of patients without sCr prior to admission: ",length(pts_insufficient_cr))
      labs_cr_aki <- labs_cr_aki[!(labs_cr_aki$patient_id %in% pts_insufficient_cr),]
      cat("\nFinal number of patients after filtering: ",length(unique(labs_cr_aki$patient_id)))
    } else {
      cat("\nPatients without pre-admission sCr allowed into analysis.\n")
    }
    
    
    # There are two possible scenarios which we have to consider when detecting each AKI event:
    # (1) AKI occurs after admission
    #   - easy to detect with the formal KDIGO definition as we only need to use older data points
    #   - we can use an approach similar to what is used in the MIMIC-III data-set: use a rolling
    #     time frame and detect the lowest Cr value to use as our baseline
    # (2) Patient presents with an AKI at the time of admission with no pre-admission labs
    #   - this makes it harder for us to determine what is the true baseline especially with limited
    #     longitudinal data
    #
    # Main issue now is determining baseline sCr for patients without any pre-admission labs.
    #
    # Methods used in literature:
    # (a) First admission sCr
    #     - Most stringent
    #     - problem is that this will underestimate AKI incidence especially as most patients will have 
    #       community-acquired AKI
    # (b) Back-deriving from MDRD with assumption of eGFR 75mL/min/1.73m2 (most commonly used method in literature)
    #     - over-estimates AKI in high proportion of CKD patients
    #     - under-estimates kidney function recovery
    # (c) Lowest sCr recorded in first admission
    #     - may be useful for patients who have no CKD and/or are younger (since more likely to recover 
    #       in-patient to baseline sCr)
    #     - may under-estimate in cases where AKIs are severe enough to cause irreversible kidney function decline
    #     - in cases of prolonged hospitalization and/or severely fluid overloaded, may (1) falsely over-estimate AKI
    #       diagnoses and misclassify AKI staging, (2) falsely increase rate of kidney function recovery
    #
    # The performance of these methods of estimating baseline sCr has been evaluated in limited settings.
    # 
    # Bagshaw et al (Nephrol Dial Transplant (2009) 24: 2739-2744, https://doi.org/10.1093/ndt/gfp159) used data from
    # the BEST Kidney study (prospective observational study from 54 ICUs in 23 countries) and compared observed premorbid
    # and estimated baseline sCr values (using method (b)). They found that using the estimated baseline sCr misclassified
    # 18.8% of patients as having AKI at the time of ICU admission and 11.7% of patients at the time of study enrolment. 
    # These numbers, when excluding CKD patients, improved to 6.6% and 4.0% respectively.
    #
    # More recently, Cooper et al (Kidney Int Rep (2021) 6, 645-656; https://doi.org/10.1016/j.ekir.2020.12.020).
    # evaluated methods (b), (c), back-calculation with MDRD assuming eGFR of 100mL/min/1.73m2, and assigning age and 
    # sex-reference GFR values. The latter three performed similarly (AUCROC 0.85, 0.85 and 0.87 respectively), whereas (b)
    # had an AUCROC of 0.72 and underestimated AKI incidence by 50%. The study population mostly younger patients with
    # low CKD premorbidity.
    #
    # Considerations when using sCr values post-AKI to estimate baseline
    # 1) Too short a window may fail to capture true baseline sCr especially in cases of prolonged time to kidney function
    #    recovery
    # 2) Too long a window, especially in lack of RRT information and fluctuating sCr in CKD patients, may falsely capture a lower 
    #    baseline sCr than is physiological (e.g. if a patient progressed to needing RRT by the 2nd admission in 3 months, but 
    #    procedure codes do not capture RRT, then we will be falsely using a lower sCr)
    # 3) Considerations of time-frame of long-term outcomes to avoid a "cyclical" argument (which does not really make sense but
    #    some reviewers still view it that way that you cannot use retrospective AKI diagnoses??)
    # 4) The method used should ideally have been evaluated in the literature before (in terms of sensitivity and specificity)
    #
    # In the initial manuscript submission, about 23% of the AKI cohort had a diagnosis code for CKD, and the actual number
    # (if we were to calculate eGFR) is likely higher. Majority of the patients (64654 out of 77927, 82.9%) were also older than 
    # 50 years old which meant a higher likelihood of CKD.
    #
    # Main considerations for minimizing errors would thus be
    # 1) Keep the estimates of baseline sCr as close to the index AKI episode as possible to prevent confounding with
    #    post-AKI recovery (i.e. the estimate must reflect the baseline sCr as close as possible to the AKI episode)
    # 2) Use pre-admission labs wherever possible (unless subsequent values show patient to recover to even lower values)
    #
    # One must also be CONSISTENT in the baseline serum creatinine definition when looking at longitudinal outcomes at different
    # time points or this will be drawn into question as to why different baselines are used and whether these outcomes can be
    # interpreted in the same context.
    #
    # Hence, main algorithm for determining baseline sCr
    # 1) If patient has pre-admission labs (or at least 2 sCr values prior to AKI onset), use these preferentially
    # 2) If patient does not have sufficient labs prior to AKI onset, then use the lowest sCr in the INDEX admission as reference
    #
    # Alternatively, these can be re-expressed as such:
    # - Use the lowest sCr value from -365 days all the way to the last day of the INDEX admission
    
    # Generate the lowest in-patient sCr from -365 days to last day of the first admission - one measure of baseline Cr
    baseline_cr_index_admit <- merge(labs_cr_aki,first_discharge,by="patient_id",all.x=TRUE) %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= ifelse(first_discharge_day < 7,7,first_discharge_day) & days_since_admission >= -365) %>% dplyr::summarise(min_cr_index_admit = min(value)) %>% dplyr::ungroup()
    
    # Scenario (1): AKI occurs during admission
    # Find minimum Cr level in a rolling 365 day timeframe
    
    # For some reason, when rlang is updated to 1.0.0 (in this case, required for tidycmprsk), the data.table statement for filling up the missing timepoints messes up when run indirectly from the package
    # If you run the data.table statement one by itself (i.e. run this code line-by-line), there is no issues in evaluation
    # So use the tidyverse-equivalent for now until we can figure a way around this
    
    cat("\nGenerating minimum Cr level in the past 365 days\n")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    attach(labs_cr_aki)
    # labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission)-365,max(days_since_admission)))]
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission=(min(days_since_admission)-365):(max(days_since_admission)+365))
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_365d = RcppRoll::roll_min(value,366,fill=NA,na.rm=TRUE,align="right")) %>% dplyr::filter(!is.na(value))
    # Find minimum Cr level in a rolling 2 day timeframe (48h)
    cat("\nGenerating minimum Cr level in the past 48h\n")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    # labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission)-2,max(days_since_admission)))]
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission=(min(days_since_admission)-2):max(days_since_admission))
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_48h = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,align="right")) %>% dplyr::filter(!is.na(value))
    
    # Scenario (2): Patient presents with an AKI already on board
    # Find minimum Cr level in a rolling 365 day timeframe
    cat("\nGenerating minimum Cr level 365 days in the future\n")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    # labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+365))]
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission=min(days_since_admission):(max(days_since_admission)+365))
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_retro_365d = RcppRoll::roll_min(value,366,fill=NA,na.rm=TRUE,align="left")) %>% dplyr::filter(!is.na(value))
    # Find minimum Cr level in a rolling 2 day timeframe (48h)
    cat("\nGenerating minimum Cr level 48h in the future\n")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    # labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+2))]
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission=min(days_since_admission):(max(days_since_admission)+2))
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_48h_retro = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,align="left")) %>% dplyr::filter(!is.na(value))
    
    
    # Another outcome we are interested in is to look at acute kidney disease, AKD (in between AKI and CKD)
    # We will use the definitions proposed for AKD as described by Chawla et. al. 2017 (ref (1))
    # We are interested in renal recovery at the 7-day and 90-day timepoint
    # We will use a cutoff of recovery to 1.25x baseline Cr as recovery, as used in ref (2)
    # References:
    # 1. Chawla, L., Bellomo, R., Bihorac, A. et al. Acute kidney disease and renal recovery: consensus report
    #    of the Acute Disease Quality Initiative (ADQI) 16 Workgroup. Nat Rev Nephrol 13, 241-257 (2017). 
    #    https://doi.org/10.1038/nrneph.2017.2
    # 2. Pannu, N., James, M., Hemmelgarn, B. & Klarenbach, S. Association between AKI, Recovery of Renal 
    #    Function, and Long-Term Outcomes after Hospital Discharge. Clinical Journal of the American Society of
    #    Nephrology 8, 194-202 (2013). https://doi.org/10.2215/CJN.06480612
    
    # Generate sCr levels at +90d (cr_90d) and +180d (cr_180d) timepoints (for determining post-AKI recovery, AKD)
    # labs_cr_aki <- data.table::setDT(labs_cr_aki)[,':='(cr_180d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+180,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][]
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    cat("\nNow finding future creatinine values at 90 days, 180 days and 365 days...\n")
    cr_90d_range <- 10 # means that the 90-day creatinine will be in the interval [80,100]
    cr_180d_range <- 30 # means that the 180-day creatinine will be in the interval [150,210]
    cr_365d_range <- 35 # means that the 365-day creatinine will be in the interval [330,400]
    
    invisible(gc())
    
    message("Warning: the package may appear to freeze at this stage. Do NOT stop the package!")
    # First find the nearest retrospective value to 90,180 and 365 days
    labs_cr_aki <- eval(data.table::setDT(labs_cr_aki)[,':='(cr_180d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+180,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1),cr_365d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+365,incbounds = TRUE)],1),cr_180d_day = tail(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+180,incbounds = TRUE)],1),cr_90d_day = tail(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1),cr_365d_day = tail(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+365,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][],envir = globalenv())
    
    invisible(gc())
    
    # Now find the nearest future value
    cat("\nNow finding the nearest future value closest to the 90,180 and 365-day mark...\n")
    message("Warning again: the package may appear to freeze at this stage. Do NOT stop the package!")
    labs_cr_aki <- eval(data.table::setDT(labs_cr_aki)[,':='(cr_180d_ul = head(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission+180,days_since_admission+180+cr_180d_range,incbounds = TRUE)],1),cr_90d_ul = head(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission+90,days_since_admission+90+cr_90d_range,incbounds = TRUE)],1),cr_365d_ul = head(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission+365,days_since_admission+365+cr_365d_range,incbounds = TRUE)],1),cr_180d_ul_day = head(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission+180,days_since_admission+180+cr_180d_range,incbounds = TRUE)],1),cr_90d_ul_day = head(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission+90,days_since_admission+90+cr_90d_range,incbounds = TRUE)],1),cr_365d_ul_day = head(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission+365,days_since_admission+365+cr_365d_range,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][],envir = globalenv())
    
    invisible(gc())
    
    detach(labs_cr_aki)
    cat("\nNow getting the closest possible value to each of the time points...")
    # 
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(cr_90d_closer = ifelse((90 - cr_90d_day) <= (cr_90d_ul_day - 90) | is.na(cr_90d_ul_day),TRUE,FALSE),cr_180d_closer = ifelse(((180 - cr_180d_day) <= (cr_180d_ul_day - 180) | is.na(cr_180d_ul_day)),TRUE,FALSE),cr_365d_closer = ifelse(((365 - cr_365d_day) <= (cr_365d_ul_day - 365) | is.na(cr_365d_ul_day)),TRUE,FALSE)) %>% dplyr::mutate(cr_90d = ifelse(isTRUE(cr_90d_closer),cr_90d,cr_90d_ul),cr_180d = ifelse(isTRUE(cr_180d_closer),cr_180d,cr_180d_ul),cr_365d = ifelse(isTRUE(cr_365d_closer),cr_365d,cr_365d_ul)) %>% dplyr::ungroup() 
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(cr_90d_day = ifelse(isTRUE(cr_90d_closer),cr_90d_day,cr_90d_ul_day),cr_180d_day = ifelse(isTRUE(cr_180d_closer),cr_180d_day,cr_180d_ul_day),cr_365d_day = ifelse(isTRUE(cr_365d_closer),cr_365d_day,cr_365d_ul_day)) %>% dplyr::ungroup()
    labs_cr_aki <- labs_cr_aki %>% dplyr::select(patient_id,days_since_admission,siteid,value,min_cr_365d,min_cr_48h,min_cr_retro_365d,min_cr_48h_retro,cr_180d,cr_90d,cr_365d,cr_180d_day,cr_90d_day,cr_365d_day,count_prior_cr)
    
    # Points to consider - 
    # 1) Should this be an average instead (in case this ends up being in the middle of another AKI)? What is the window
    #    for qualifying a Cr value as 180 day/90 day? (given that patients may not be followed up exactly at 90 days or 180 days)
    # 2) Should we use this to generate the 1-year outcome? Will need more code to determine if there were sufficient prior
    #    Cr values to determine baseline before admission
    
    # At this point, our table has these headers:
    # patient_id  siteid  days_since_admission  value min_cr_365d min_cr_48h  min_cr_retro_365d min_cr_48h_retro  cr_180d cr_90d
    
    # Now we have to start grading AKI severity at each time point
    # This approach is similar to how the MIMIC-III dataset generates AKI severity
    # Generate two columns using both the formal KDIGO AKI definition and the modified retrospective AKI definition
    cat("\nGenerating KDIGO severity grades for each serum Cr value")
    labs_cr_aki$aki_kdigo <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade)
    labs_cr_aki$aki_kdigo_retro <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade_retro)
    # labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(aki_kdigo_final = max(aki_kdigo))
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(aki_kdigo_final = max(aki_kdigo,aki_kdigo_retro))

    # Merge in baseline_cr_index_admit and do grading of AKIs based on this
    labs_cr_aki <- merge(labs_cr_aki,baseline_cr_index_admit,by="patient_id",all.x=TRUE)
    labs_cr_aki$aki_kdigo_index_baseline <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade_index_baseline)
    
    # Now we are going to generate the start days of each AKI
    labs_cr_aki_tmp <- labs_cr_aki
    labs_cr_aki_tmp$valid = 1
    
    # Find the day of the minimum Cr used for grading AKIs (taken as baseline)
    cat("\nNow finding the day at which the minimum serum Cr is achieved")
    labs_cr_aki_tmp <- labs_cr_aki_tmp %>% dplyr::group_by(patient_id) %>% tidyr::complete(days_since_admission = tidyr::full_seq(days_since_admission,1)) %>% dplyr::mutate(value = zoo::na.fill(value,Inf))
    labs_cr_aki_tmp2 <- labs_cr_aki_tmp
    labs_cr_aki_tmp3 <- labs_cr_aki_tmp
    labs_cr_aki_tmp2 <- labs_cr_aki_tmp2 %>% split(.$patient_id) %>% purrr::map(~pos_min(.$value,.$days_since_admission)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
    colnames(labs_cr_aki_tmp2)[2] <- "day_min"
    labs_cr_aki_tmp3 <- labs_cr_aki_tmp3 %>% split(.$patient_id) %>% purrr::map(~pos_min(.$value,.$days_since_admission,lag=FALSE)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
    colnames(labs_cr_aki_tmp3)[2] <- "day_min_retro"
    labs_cr_aki_tmp4 <- cbind(labs_cr_aki_tmp,"day_min" = labs_cr_aki_tmp2$day_min,"day_min_retro" = labs_cr_aki_tmp3$day_min_retro)
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[!is.na(labs_cr_aki_tmp4$valid),]
    
    # Generate delta_cr
    cat("\nNow identifying maxima points of serum Cr")
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(min_cr_365d_final = min(min_cr_365d,min_cr_retro_365d,na.rm=TRUE)) %>% dplyr::mutate(delta_cr = value - min_cr_365d_final)
    
    # Use the largest delta_cr to find the peak of each AKI
    labs_cr_aki_delta_maxima <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr %in% delta_cr[which.peaks(delta_cr,decreasing=FALSE)])
    labs_cr_aki_delta_maxima$delta_is_max = 1
    labs_cr_aki_delta_maxima <- labs_cr_aki_delta_maxima %>% dplyr::rename(delta_maxima = delta_cr) %>% dplyr::select(patient_id,days_since_admission,delta_maxima,delta_is_max)
    labs_cr_aki_tmp4 <- merge(labs_cr_aki_tmp4,labs_cr_aki_delta_maxima,by=c("patient_id","days_since_admission"),all.x=TRUE)
    
    # Filter for KDIGO grades > 0
    cat("\nNow generating tables of all AKI events")
    labs_cr_aki_tmp5 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final > 0,]
    labs_cr_aki_tmp5[is.na(labs_cr_aki_tmp5)] <- 0
    # Filter for maxima of delta_cr (which should give us the peaks)
    labs_cr_aki_tmp5 <- labs_cr_aki_tmp5[labs_cr_aki_tmp5$delta_is_max > 0,]
    
    # Filter and reorder columns to generate our final table of all AKI events
    labs_aki_summ <- labs_cr_aki_tmp5 %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_365d,min_cr_48h,min_cr_retro_365d,min_cr_48h_retro,min_cr_365d_final,min_cr_index_admit,cr_180d,cr_90d,cr_365d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,aki_kdigo_index_baseline,cr_180d_day,cr_90d_day,cr_365d_day,count_prior_cr)
    # labs_aki_summ <- labs_cr_aki_tmp5 %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_365d,min_cr_48h,min_cr_retro_365d,min_cr_48h_retro,min_cr_365d_final,min_cr_index_admit,cr_180d,cr_90d,cr_365d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,aki_kdigo_index_baseline,akd_180d,akd_90d,cr_180d_day,cr_90d_day,cr_365d_day)
    
    labs_aki_summ <- labs_aki_summ %>% dplyr::distinct(patient_id,days_since_admission,.keep_all=TRUE)
    labs_aki_summ <- labs_aki_summ %>% dplyr::filter(days_since_admission >= 0)
    
    # Merge in the first discharge day and label the AKIs whether they occurred during the first admission
    labs_aki_summ <- merge(labs_aki_summ,first_discharge,by="patient_id",all.x=TRUE)
    labs_aki_summ <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::mutate(first_admit = dplyr::if_else(days_since_admission >= 0 & days_since_admission <= first_discharge_day,1,0)) %>% dplyr::ungroup() %>% dplyr::distinct()
    
    aki_list_first_admit <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::summarise(first_admit_aki = max(first_admit)) %>% dplyr::filter(first_admit_aki == 1) %>% dplyr::ungroup()
    aki_list_first_admit <- unlist(unique(aki_list_first_admit$patient_id))
    # Final headers for labs_aki_summ:
    # patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_365d,min_cr_48h,min_cr_retro_365d,min_cr_48h_retro,min_cr_365d_final,cr_180d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d,cr_180d_day,cr_90d_day,first_discharge_day,first_admit
    # days_since_admission - time at which peak Cr is achieved
    # day_min - time at which Cr begins to rise
    
    # Generate a separate table (for reference) of all creatinine peaks not fulfilling KDIGO AKI criteria
    
    #Create table of sCr values for those without AKIs in index admission
    cat("\nNow generating a separate table for serum Cr peaks which do not reach AKI definitions")
    labs_cr_nonaki <- labs_cr_aki_tmp4[!(labs_cr_aki_tmp4$patient_id %in% aki_list_first_admit),]
    # labs_cr_nonaki <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final == 0,]
    labs_cr_nonaki[is.na(labs_cr_nonaki)] <- 0
    # labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$delta_is_max > 0,]
    labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_365d,min_cr_48h,min_cr_retro_365d,min_cr_48h_retro,min_cr_365d_final,min_cr_index_admit,cr_180d,cr_90d,cr_365d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,aki_kdigo_index_baseline,cr_180d_day,cr_90d_day,cr_365d_day,count_prior_cr)
    # labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_365d,min_cr_48h,min_cr_retro_365d,min_cr_48h_retro,min_cr_365d_final,min_cr_index_admit,cr_180d,cr_90d,cr_365d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,aki_kdigo_index_baseline,akd_180d,akd_90d,cr_180d_day,cr_90d_day,cr_365d_day)
    labs_cr_nonaki_tmp <- labs_cr_nonaki %>% dplyr::filter(days_since_admission >= 0)
    # Generate the highest Cr peak for non-AKI peaks detected in the first admission
    labs_cr_nonaki_tmp <- merge(labs_cr_nonaki_tmp,first_discharge,by="patient_id",all.x=TRUE) %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= first_discharge_day) %>% dplyr::ungroup()
    
    labs_nonaki_summ <- labs_cr_nonaki_tmp %>% dplyr::group_by(patient_id) %>% dplyr::slice(which.max(delta_cr))
    # Merge in the first discharge day and label the AKIs whether they occurred during the first admission
    # labs_nonaki_summ <- merge(labs_nonaki_summ,first_discharge,by="patient_id",all.x=TRUE)
    labs_nonaki_summ <- labs_nonaki_summ %>% dplyr::group_by(patient_id) %>% dplyr::mutate(first_admit = dplyr::if_else(days_since_admission >= 0 & days_since_admission <= first_discharge_day,1,0)) %>% dplyr::ungroup()
    
    # We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset 
    # (2) how long before/after disease severity
    # These tables will help in segregating the populations for analysis later
    cat("\nGenerating intermediate tables to determine if AKIs occured before or after severe COVID-19 onset")
    severe_time <- demographics_filt %>% dplyr::select(patient_id,severe,time_to_severe)
    labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = ifelse(!is.na(time_to_severe), time_to_severe - day_min,NA))
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0)) %>% dplyr::ungroup()
    labs_aki_severe <- labs_aki_severe %>% dplyr::distinct()
    
    labs_nonaki_severe <- merge(labs_nonaki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_nonaki_severe$severe_to_aki <- NA 
    labs_nonaki_severe$severe_before_aki <- 1 
    
    ## Save the generated AKI tables for future reference / debugging (note: these will NOT be uploaded!!)
    #write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)
    #write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)
    
    invisible(gc())
    
    ## ==================================================================================
    ## PART 2: Serum Creatinine Trends - Plots against Time from Peak Serum Creatinine
    ## ==================================================================================
    
    # Severity indices
    # (1) non-severe, no AKI
    # (2) non-severe, AKI
    # (3) severe, no AKI
    # (4) severe, AKI
    
    # Novel antivirals indices
    # (1) non-severe, no novel anti-COVID-19 agents (i.e. severe == 1 or 2 && covidrx == 0)
    # (2) non-severe, with novel anti-COVID-19 agents (i.e. severe == 1 or 2 && covidrx == 1)
    # (3) severe, no novel anti-COVID-19 agents (i.e. severe == 3, 4 && covidrx == 0)
    # (4) severe, with novel anti-COVID-19 agents (i.e. severe == 3, 4 && covidrx == 1)
    
    # This analysis uses the index AKI episode
    # Feel free to modify the code to look at other types of AKI episodes, e.g. you can dplyr::filter for the most severe AKI episode
    
    cat("\nNow generating tables for use for plotting normalised serum Cr values against time")
    cat("\nFirst creating table for AKI patients only...")
    # First, dplyr::filter the labs_aki_severe table to show only the index AKI episodes
    # aki_only_index <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% tidyr::fill(severe) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission,na.rm=TRUE)) %>% dplyr::filter(delta_cr == max(delta_cr)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE)
    # patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_365d	min_cr_48h	min_cr_retro_365d	min_cr_48h_retro	min_cr_365d_final	cr_180d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	severe  time_to_severe	severe_to_aki	severe_before_aki
    
    aki_only_index <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% tidyr::fill(severe) %>% dplyr::filter(first_admit == 1) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(count_prior_cr >= 2) %>% dplyr::filter(days_since_admission == min(days_since_admission,na.rm=TRUE)) %>% dplyr::filter(delta_cr == max(delta_cr,na.rm=TRUE)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE)
    
    patients_no_prior_cr <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(first_admit == 1) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(count_prior_cr < 2) %>% dplyr::filter(days_since_admission == min(days_since_admission,na.rm=TRUE)) %>% dplyr::filter(delta_cr == max(delta_cr,na.rm=TRUE)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE) %>% dplyr::ungroup()
    patients_no_prior_cr <- unlist(unique(patients_no_prior_cr$patient_id))
    patients_no_prior_cr <- setdiff(patients_no_prior_cr,unlist(unique(aki_only_index$patient_id)))
    
    
    # Generate the patient list including (1) severity indices from this dplyr::filtered table (2) day of peak Cr
    # severe - 2 = never severe, 4 = severe, AKI
    # aki_only_index <- aki_only_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = dplyr::if_else(!is.na(time_to_severe),4,2)) %>% dplyr::ungroup()
    cat("\nAssigning new severity labels...")
    aki_only_index <- aki_only_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = 2 * severe + 2) %>% dplyr::ungroup() # if severe = 0 initially, will be recoded as 2; if severe = 1 initially, will be recoded as 4
    
    # create the change in baseline index table
    cat("\nCreating baseline shift tables for only AKI patients...")
    aki_only_index_baseline_shift <- aki_only_index %>% dplyr::select(patient_id,severe,cr_365d,cr_180d,cr_90d,cr_365d_day,cr_180d_day,cr_90d_day)
    aki_only_index <- aki_only_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki,count_prior_cr)
    colnames(aki_only_index)[2] <- "peak_cr_time"
    colnames(aki_only_index)[4] <- "aki_start"
    # Headers of aki_only_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki
    cat("\n\nNow creating list of non-AKI patients...")
    no_aki_list <- demographics_filt %>% dplyr::select(patient_id,severe)
    no_aki_list <- no_aki_list[!(no_aki_list$patient_id %in% aki_only_index$patient_id),]
    no_aki_list <- no_aki_list[(no_aki_list$patient_id %in% pts_valid_cr),]
    no_aki_list <- no_aki_list %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(is.na(severe),1,2 * severe + 1)) %>% dplyr::ungroup() # if severe = 0, will be recoded as 1; if severe = 1, will be recoded as 3
    
    labs_nonaki_summ <- labs_nonaki_summ[labs_nonaki_summ$patient_id %in% no_aki_list$patient_id,]
    labs_nonaki_severe <- labs_nonaki_severe[labs_nonaki_severe$patient_id %in% no_aki_list$patient_id,]
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$patient_id %in% no_aki_list$patient_id,]
    
    cat("\nCreating table for non-AKI patients...")
    # Create a non-AKI equivalent for aki_only_index - except that this takes the largest delta_cr (and the earliest occurrence of such a delta_cr)
    no_aki_index <- labs_nonaki_severe %>% dplyr::group_by(patient_id) %>% dplyr::arrange(severe,.by_group = TRUE) %>% tidyr::fill(severe) %>% dplyr::filter(delta_cr == max(delta_cr, na.rm = TRUE)) %>% dplyr::filter(days_since_admission == min(days_since_admission, na.rm = TRUE)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE) %>% dplyr::ungroup()
    no_aki_index <- no_aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(is.na(severe),1,2 * severe + 1))
    
    # create the change in baseline index table
    no_aki_index_baseline_shift <- no_aki_index %>% dplyr::select(patient_id,severe,cr_365d,cr_180d,cr_90d,cr_365d_day,cr_180d_day,cr_90d_day)
    no_aki_index <- no_aki_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki,count_prior_cr)
    colnames(no_aki_index)[2] <- "peak_cr_time"
    colnames(no_aki_index)[4] <- "aki_start"
    
    #no_aki_index$severe_to_aki <- -999
    cat("\nCreating final table of peak Cr for all patients...")
    aki_index <- dplyr::bind_rows(aki_only_index,no_aki_index)
    
    # =================================
    # Lancet reviewer edit: Generate table categorizing occurrence of baseline sCr value
    # =================================
    cat("\nNow attempting to categorize the day of baseline sCr...\n")
    
    earliest_cr <- aki_index[,c("patient_id","aki_start")] %>% dplyr::group_by(patient_id) %>% dplyr::filter(aki_start == min(aki_start,na.rm=T)) %>% dplyr::ungroup() %>% dplyr::distinct()
    # earliest_cr <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::summarise(earliest_cr_day = min(days_since_admission,na.rm=T)) %>% dplyr::ungroup()
    earliest_cr <- earliest_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(preadmit_cr_period = dplyr::case_when(
      aki_start < -180 ~ "181_to_365_days",
      aki_start < -90 ~ "91_to_180_days",
      TRUE ~ "zero_to_90_days"
    )) %>% dplyr::ungroup()
    
    # =================================
    # End edit
    # =================================
    
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        cat("\nAdding on temporal information for experimental COVID-19 antivirals...")
        aki_index <- merge(aki_index,med_covid19_new,by="patient_id",all.x=TRUE)
        aki_index$covid_rx[is.na(aki_index$covid_rx)] <- 0
        aki_index <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covidrx_grp = dplyr::if_else(severe <= 2, dplyr::if_else(covid_rx == 0,1,2),dplyr::if_else(covid_rx == 0,3,4))) %>% dplyr::ungroup()
        cat(paste0(c("Column names for aki_index",colnames(aki_index))))
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,covidrx_grp,count_prior_cr)
    } else {
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,count_prior_cr)
    }
    # Headers of aki_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki  covidrx_grp
    aki_index <- aki_index %>% dplyr::arrange(patient_id,peak_cr_time,desc(severe))%>% dplyr::distinct(patient_id,peak_cr_time,.keep_all = TRUE)
    
    cat("\nFinal table aki_index created.")
    
    # Uncomment the following line to remove patients who were previously on RRT prior to admission
    # aki_index <- aki_index[!(aki_index$patient_id %in% rrt_new),]
    cat("\nNow generating table of normalised serum Cr...")
    # Create a common labs_cr_all table containing the serum Cr values, the severity groupings and anti-viral groupings
    labs_cr_aki_tmp <- labs_cr_aki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_365d,min_cr_retro_365d)
    labs_cr_nonaki_tmp <- labs_cr_nonaki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_365d,min_cr_retro_365d)
    labs_cr_all <- dplyr::bind_rows(labs_cr_aki_tmp,labs_cr_nonaki_tmp)
    labs_cr_all <- merge(labs_cr_all,aki_index,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    labs_cr_all <- labs_cr_all[labs_cr_all$patient_id %in% aki_index$patient_id,]
    # Now, generate a table containing lab values with timepoints calculated from time of peak cr
    peak_trend <- labs_cr_all
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_365d,min_cr_retro_365d)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_365d,min_cr_retro_365d)
    }
    # patient_id  severe  covidrx_grp  days_since_admission  peak_cr_time  value min_cr_365d min_cr_retro_365d
    
    # Calculate the day from peak Cr
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_peak = days_since_admission - peak_cr_time) %>% dplyr::ungroup()
    ## dplyr::filter this table for Cr values that fall within the 7 day window
    # peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(between(time_from_peak,0,7)) %>% dplyr::ungroup()
    # Normalise to baseline values used for AKI calculation
    # peak_trend <- peak_trend %>% dplyr::group_by(patient_id,time_from_peak) %>% dplyr::mutate(baseline_cr = min(min_cr_365d)) %>% dplyr::ungroup()
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id,time_from_peak) %>% dplyr::mutate(baseline_cr = min(min_cr_365d,min_cr_retro_365d)) %>% dplyr::ungroup()
    
    # In the event longitudinal data becomes very long, we will create a column where the very first baseline Cr for the index AKI is generated for each patient
    first_baseline <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak == 0) %>% dplyr::filter(baseline_cr == min(baseline_cr,na.rm=TRUE)) %>% dplyr::select(patient_id,baseline_cr) %>% dplyr::distinct()
    colnames(first_baseline)[2] <- "first_baseline_cr"
    
    first_baseline_prioronly <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak == 0) %>% dplyr::filter(min_cr_365d == min(min_cr_365d,na.rm=TRUE)) %>% dplyr::select(patient_id,min_cr_365d) %>% dplyr::distinct()
    colnames(first_baseline_prioronly)[2] <- "first_baseline_cr_prioronly"
    
    # Merge the first_baseline table back in
    peak_trend <- merge(peak_trend,first_baseline,by="patient_id",all.x=TRUE)
    
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_365d,min_cr_retro_365d,time_from_peak,baseline_cr,first_baseline_cr)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_365d,min_cr_retro_365d,time_from_peak,baseline_cr,first_baseline_cr)
    }
    
    # TESTING: add in the (-365,last day of admission) definition of baseline Cr
    peak_trend <- merge(peak_trend,baseline_cr_index_admit,by="patient_id",all.x=TRUE)
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::arrange(severe,first_baseline_cr) %>% tidyr::fill(severe,first_baseline_cr) %>% dplyr::mutate(ratio = value/first_baseline_cr, ratio_prioronly = value/min_cr_365d, ratio_baseline_index = value / min_cr_index_admit) %>% dplyr::ungroup() %>% dplyr::distinct()
    # peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::arrange(severe,first_baseline_cr) %>% tidyr::fill(severe,first_baseline_cr) %>% dplyr::mutate(ratio = value/first_baseline_cr, ratio_prioronly = value/min_cr_365d) %>% dplyr::ungroup() %>% dplyr::distinct()
    cat("\nFinal table of peak Cr for all patients - peak_trend - created.")
    # peak_trend will now be a common table to plot from the selected AKI peak
    
    invisible(gc())
    
    
    
    # =================================
    # Lancet edit: Add crude classification of CKD staging using sCr cutoffs (2009 CKD-EPI Equation, Male 60 year old, Non-black)
    # eGFR 90 - 0.9206857472204mg/dL
    # eGFR 60 - 1.2875430484355mg/dL
    # =================================
    ckd_staging <- first_baseline %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ckd_stage = dplyr::case_when(
      first_baseline_cr <= 0.92 ~ 0,
      first_baseline_cr <= 1.29 ~ 1,
      TRUE ~ 2
    )) %>% dplyr::ungroup()
    
    # =======================================================================================
    # Filter the demographics, peak_trend, comorbid, meds tables to remove those patients who 
    # do not have at least 1 measurement prior to peak Cr
    # =======================================================================================
    cat("\nFiltering to exclude patients without at least 2 prior Cr values to peak / only one sCr value in first admission.\n")
    cat("No. of patients without at least 2 Cr values prior to peak: ",length(patients_no_prior_cr),"\n")
    cat("No. of patients with only one sCr in first admission: ",length(pts_insufficient_cr),"\n")
    patients_no_prior_cr <- unique(c(patients_no_prior_cr,pts_insufficient_cr))
    cat("Total number of patients to be excluded for insufficient sCr: ",length(patients_no_prior_cr),"\n")
    
    if(isFALSE(allow_no_prior_cr)) {
      demographics_filt <- demographics_filt[!(demographics_filt$patient_id %in% patients_no_prior_cr),]
      observations <- observations[!(observations$patient_id %in% patients_no_prior_cr),]
      course <- course[!(course$patient_id %in% patients_no_prior_cr),]
      first_discharge <- first_discharge[!(first_discharge$patient_id %in% patients_no_prior_cr),]
      diagnosis <- diagnosis[!(diagnosis$patient_id %in% patients_no_prior_cr),]
      procedures <- procedures[!(procedures$patient_id %in% patients_no_prior_cr),]
      
      labs_cr_aki <- labs_cr_aki[!(labs_cr_aki$patient_id %in% patients_no_prior_cr),]
      labs_cr_nonaki <- labs_cr_nonaki[!(labs_cr_nonaki$patient_id %in% patients_no_prior_cr),]
      labs_aki_summ <- labs_aki_summ[!(labs_aki_summ$patient_id %in% patients_no_prior_cr),]
      labs_aki_severe <- labs_aki_severe[!(labs_aki_severe$patient_id %in% patients_no_prior_cr),]
      labs_nonaki_summ <- labs_nonaki_summ[!(labs_nonaki_summ$patient_id %in% patients_no_prior_cr),]
      labs_nonaki_severe <- labs_nonaki_severe[!(labs_nonaki_severe$patient_id %in% patients_no_prior_cr),]
      labs_cr_all <- labs_cr_all[!(labs_cr_all$patient_id %in% patients_no_prior_cr),]
      
      aki_only_index <- aki_only_index[!(aki_only_index$patient_id %in% patients_no_prior_cr),]
      no_aki_index <- no_aki_index[!(no_aki_index$patient_id %in% patients_no_prior_cr),]
      aki_index <- aki_index[!(aki_index$patient_id %in% patients_no_prior_cr),]
      
      aki_only_index_baseline_shift <- aki_only_index_baseline_shift[!(aki_only_index_baseline_shift$patient_id %in% patients_no_prior_cr),]
      no_aki_index_baseline_shift <- no_aki_index_baseline_shift[!(no_aki_index_baseline_shift$patient_id %in% patients_no_prior_cr),]
      peak_trend <- peak_trend[!(peak_trend$patient_id %in% patients_no_prior_cr),]
      
      medications <- medications[!(medications$patient_id %in% patients_no_prior_cr),]
      med_chronic <- med_chronic[!(med_chronic$patient_id %in% patients_no_prior_cr),]
      med_new <- med_new[!(med_new$patient_id %in% patients_no_prior_cr),]
      if(isTRUE(coaga_present)) {
        med_coaga_new <- med_coaga_new[!(med_coaga_new$patient_id %in% patients_no_prior_cr),]
      }
      if(isTRUE(coagb_present)) {
        med_coagb_new <- med_coagb_new[!(med_coagb_new$patient_id %in% patients_no_prior_cr),]
      }
      if(isTRUE(remdesivir_present) | isTRUE(covid19antiviral_present)) {
        med_covid19_new <- med_covid19_new[!(med_covid19_new$patient_id %in% patients_no_prior_cr),]
      }
    } else {
      cat("\nNo exclusion done for no prior sCr.")
    }
    
    
    # ================================
    # Identifying patients to exclude
    # ================================
    
    # We now need to filter for patients with CKD4/5 using surrogate cutoffs (within the constraints of Phase 2.1)
    # This is where the ckd_cutoff value comes in useful
    # Default: 2.25mg/dL
    cat("\nExcluding patients with baseline Cr >= ",ckd_cutoff,"mg/dL...\n")
    esrf_list <- unlist(first_baseline$patient_id[first_baseline$first_baseline_cr >= ckd_cutoff])
    cat("\nNo. of patients with baseline Cr >= 2.25mg/dL: ",length(esrf_list))
    # We will now use the surrogate detection for RRTs
    # We are using this method as there are lack of RRT procedure codes in some sites and Phase 1.1/2.1 does not define these clearly
    # We also do not have access to urea values
    # Briefly, if the peak Cr >= 3mg/dL and there is a drop of >= 25% in sCr in 24h or less, this is highly likely to indicate RRT
    if(isTRUE(print_rrt_surrogate) | isTRUE(use_rrt_surrogate)) {
        cat("\n\n==========================\nRRT Surrogate Detection\n==========================\n")
        cat("\nPerforming RRT Surrogate detection with Cr cutoff = 3 and ratio = 0.75...")
        rrt_detection <- detect_rrt_drop(labs_cr_aki,cr_abs_cutoff = 3,ratio = 0.75)
        if(isTRUE(print_rrt_surrogate)) {
            cat("\nPrinting RRT Surrogate Detection Tables. Ensure this is NOT uploaded as it contains patient-level data!")
            cat("\nLook in ~/FourCePhase2.1AKI for the following files:")
            cat("\na) DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv")
            cat("\nb) DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv")
            cat("\nThese files contain patient_id and days corresponding to suspected RRT episodes. Use these for manual chart review.")
            write.csv(rrt_detection,file=file.path(dir.output, "DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv"),row.names=FALSE)
            write.csv(rrt_detection[!(rrt_detection$patient_id %in% esrf_list),],file=file.path(dir.output, "DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv"),row.names=FALSE)
        }
        
    }

    if(isTRUE(use_rrt_surrogate)) {
        cat("\nYou have opted to use the RRT surrogate results as exclusion criteria as well. Excluding these patients.")
        rrt_detection_prior_list <- unlist(rrt_detection$patient_id[rrt_detection$days_since_admission < 0])
        cat("\nFinal breakdown of omission by RRT/RRT surrogate criteria:
            \nExceeding baseline cutoff of ",ckd_cutoff,"mg/dL: ",length(esrf_list),
            "\nRRT surrogate criteria: ",length(rrt_detection_prior_list),
            "\nRRT codes: ",length(rrt_old))
        exclusion_criteria <- data.frame(c("RRT","No_Prior_sCr","Baseline_2.25","RRT_Surrogate","Total"),c(length(rrt_old),length(patients_no_prior_cr),length(esrf_list),length(rrt_detection_prior_list),length(unique(c(patients_no_prior_cr,esrf_list,rrt_detection_prior_list,rrt_old)))))
        colnames(exclusion_criteria) <- c("Exclusion Criteria","Counts")
        exclusion_criteria$siteid <- currSiteId
        write.csv(exclusion_criteria,file=file.path(dir.output, paste0(currSiteId,"_ExclusionCriteria_Counts.csv")),row.names=FALSE)
        esrf_list <- unique(c(esrf_list,rrt_detection_prior_list,rrt_old))
    } else {
      cat("\nFinal breakdown of omission by RRT/RRT surrogate criteria:
            \nExceeding baseline cutoff of ",ckd_cutoff,"mg/dL: ",length(esrf_list),
          "\nRRT codes: ",length(rrt_old))
      exclusion_criteria <- data.frame(c("RRT","No_Prior_sCr","Baseline_2.25","RRT_Surrogate","Total"),c(length(rrt_old),length(patients_no_prior_cr),length(esrf_list),0,length(unique(c(patients_no_prior_cr,esrf_list,rrt_old)))))
      colnames(exclusion_criteria) <- c("Exclusion Criteria","Counts")
      exclusion_criteria$siteid <- currSiteId
      write.csv(exclusion_criteria,file=file.path(dir.output, paste0(currSiteId,"_ExclusionCriteria_Counts.csv")),row.names=FALSE)
      esrf_list <- unique(c(esrf_list,rrt_old))
    }
    cat("Total excluded for either baseline sCr/RRT code/RRT surrogate criteria: ",length(esrf_list))
    
    # =======================================================================================
    # Filter the demographics, peak_trend, comorbid, meds tables to remove those with CKD4/5
    # =======================================================================================
    cat("\n\nFiltering to exclude patients with CKD4/5 and/or RRT/kidney transplant procedure/diagnoses codes.")
    demographics_filt <- demographics_filt[!(demographics_filt$patient_id %in% esrf_list),]
    observations <- observations[!(observations$patient_id %in% esrf_list),]
    course <- course[!(course$patient_id %in% esrf_list),]
    first_discharge <- first_discharge[!(first_discharge$patient_id %in% esrf_list),]
    diagnosis <- diagnosis[!(diagnosis$patient_id %in% esrf_list),]
    procedures <- procedures[!(procedures$patient_id %in% esrf_list),]
    
    labs_cr_aki <- labs_cr_aki[!(labs_cr_aki$patient_id %in% esrf_list),]
    labs_cr_nonaki <- labs_cr_nonaki[!(labs_cr_nonaki$patient_id %in% esrf_list),]
    labs_aki_summ <- labs_aki_summ[!(labs_aki_summ$patient_id %in% esrf_list),]
    labs_aki_severe <- labs_aki_severe[!(labs_aki_severe$patient_id %in% esrf_list),]
    labs_nonaki_summ <- labs_nonaki_summ[!(labs_nonaki_summ$patient_id %in% esrf_list),]
    labs_nonaki_severe <- labs_nonaki_severe[!(labs_nonaki_severe$patient_id %in% esrf_list),]
    labs_cr_all <- labs_cr_all[!(labs_cr_all$patient_id %in% esrf_list),]
    
    aki_only_index <- aki_only_index[!(aki_only_index$patient_id %in% esrf_list),]
    no_aki_index <- no_aki_index[!(no_aki_index$patient_id %in% esrf_list),]
    aki_index <- aki_index[!(aki_index$patient_id %in% esrf_list),]
    
    aki_only_index_baseline_shift <- aki_only_index_baseline_shift[!(aki_only_index_baseline_shift$patient_id %in% esrf_list),]
    no_aki_index_baseline_shift <- no_aki_index_baseline_shift[!(no_aki_index_baseline_shift$patient_id %in% esrf_list),]
    peak_trend <- peak_trend[!(peak_trend$patient_id %in% esrf_list),]
    
    medications <- medications[!(medications$patient_id %in% esrf_list),]
    med_chronic <- med_chronic[!(med_chronic$patient_id %in% esrf_list),]
    med_new <- med_new[!(med_new$patient_id %in% esrf_list),]
    if(isTRUE(coaga_present)) {
        med_coaga_new <- med_coaga_new[!(med_coaga_new$patient_id %in% esrf_list),]
    }
    if(isTRUE(coagb_present)) {
        med_coagb_new <- med_coagb_new[!(med_coagb_new$patient_id %in% esrf_list),]
    }
    if(isTRUE(remdesivir_present) | isTRUE(covid19antiviral_present)) {
        med_covid19_new <- med_covid19_new[!(med_covid19_new$patient_id %in% esrf_list),]
    }

    # ============================================
    # PREPARING OTHER PRE-REQUISITE TABLES
    # 1) Comorbidities
    # 2) Intubation (not used in current analyses)
    # 
    # NOTE: In versions prior to 0.1.5, this segment was placed before the AKI detection code.
    # However, in view of implementing filtering for CKD4/5 patients, these segments have been shifted here to facilitate
    # filtering and avoid having to add additional unnecessary lines of code (and make this more readable)
    # ============================================
    
    # Comorbidities & Prothrombotic Events
    # ======================================
    cat("\nNow creating tables for comorbidities and admission diagnoses...")
    cat("\nCreating Comorbidities table...")
    # Filter comorbids from all diagnoses
    comorbid_icd9 <- diag_icd9[diag_icd9$icd_code %in% comorbid_icd9_ref$icd_code,]
    comorbid_icd10 <- diag_icd10[diag_icd10$icd_code %in% comorbid_icd10_ref$icd_code,]
    # Filter comorbids to be restricted to diagnosis codes made -15days
    comorbid_icd9 <- comorbid_icd9[comorbid_icd9$days_since_admission < -15,-c(2:4)]
    comorbid_icd10 <- comorbid_icd10[comorbid_icd10$days_since_admission < -15,-c(2:4)]
    # Map the comorbid codes
    comorbid_icd9 <- merge(comorbid_icd9,comorbid_icd9_ref,by="icd_code",all.x=TRUE)
    comorbid_icd10 <- merge(comorbid_icd10,comorbid_icd10_ref,by="icd_code",all.x=TRUE)
    comorbid <- rbind(comorbid_icd9,comorbid_icd10)
    comorbid <- comorbid %>% dplyr::select(patient_id,comorbid_type)
    comorbid$present <- 1
    comorbid <- comorbid %>% dplyr::distinct()
    comorbid <- comorbid %>% tidyr::spread(comorbid_type,present)
    comorbid[is.na(comorbid)] <- 0
    
    comorbid_list <- colnames(comorbid)[-1]
    cat("\nComorbids in data set: ",comorbid_list)
    # Final headers for comorbid table
    # Note: order of columns may depend on the overall characteristics of your patient population
    # patient_id	ihd,htn,dm,asthma,copd,bronchiectasis,ild,ckd,pe,dvt,cancer
    # Values stored are in binary, 1 = present, 0 = absent
    
    # Time to Intubation
    # ====================
    cat("\nCreating table for intubation...")
    # (1) Determine intubation from procedure code
    intubation_code <- c("0BH13EZ","0BH17EZ","0BH18EZ","0B21XEZ","5A09357","5A09358","5A09359","5A0935B","5A0935Z","5A09457","5A09458","5A09459","5A0945B","5A0945Z","5A09557","5A09558","5A09559","5A0955B","5A0955Z","96.7","96.04","96.70","96.71","96.72")
    intubation <- procedures[procedures$procedure_code %in% intubation_code,-c(2,4)]
    intubation <- intubation[,c(1,2)]
    intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission),]
    intubation <- intubation[!duplicated(intubation$patient_id),]
    # (2) In some cases intubation may not be coded as a procedure. Hence surrogate way is to determine if 
    #     patient had been diagnosed with ARDS and/or VAP
    vap_ards_codes <- c("J80","J95.851","518.82","997.31","518","997","J95")
    vap_ards_diag <- diagnosis[diagnosis$icd_code %in% vap_ards_codes,]
    intubation <- rbind(vap_ards_diag[,c(1,3)],intubation)
    intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission),]
    intubation <- intubation[!duplicated(intubation$patient_id),]
    
    # Final table headers:
    # patient_id	days_since_admission
    # days_since_admission = time to first intubation event
    
    
    # ==============================
    # Baseline Shift Tables
    # ==============================
    cat("\nCreating baseline shift counts (125% sCr definition)...")
    baseline_shift <- dplyr::bind_rows(aki_only_index_baseline_shift,no_aki_index_baseline_shift) %>% dplyr::distinct()
    baseline_shift <- merge(baseline_shift,first_baseline,by="patient_id",all.x=T)
    
    # Filter shifts based on omission criteria for Cr at 90d (80-100), 180d (150-210), 365d (330-390)
    baseline_shift <- baseline_shift %>% dplyr::group_by(patient_id) %>% dplyr::mutate(cr_90d = ifelse(cr_90d_day >= 80 & cr_90d_day <= 100,cr_90d,NA_real_),cr_180d = ifelse(cr_180d_day >= 150 & cr_180d_day <= 210,cr_180d,NA_real_),cr_365d = ifelse(cr_365d_day >= 330 & cr_365d_day <= 390,cr_365d,NA_real_)) %>% dplyr::filter(!is.na(cr_365d) | !is.na(cr_180d) | !is.na(cr_90d)) %>% dplyr::ungroup()
    baseline_shift <- baseline_shift %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio_365d = cr_365d/first_baseline_cr, ratio_180d = cr_180d/first_baseline_cr,ratio_90d = cr_90d/first_baseline_cr) %>% dplyr::select(patient_id,severe,ratio_365d,ratio_180d,ratio_90d)
    baseline_shift_summ <- baseline_shift %>% dplyr::group_by(severe) %>% dplyr::summarise(n_all=dplyr::n(),n_all_365d = sum(ratio_365d > 0, na.rm=T), n_all_180d = sum(ratio_180d > 0, na.rm = T), n_all_90d = sum(ratio_90d > 0, na.rm = T), n_shift_365d = sum(ratio_365d >= 1.25, na.rm=T),n_shift_180d = sum(ratio_180d >= 1.25, na.rm=T),n_shift_90d = sum(ratio_90d >= 1.25, na.rm=T)) %>% dplyr::ungroup()
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        baseline_shift_summ$n_all[baseline_shift_summ$n_all < obfuscation_value] <- NA
        baseline_shift_summ$n_all_365d[baseline_shift_summ$n_all_365d < obfuscation_value] <- NA
        baseline_shift_summ$n_all_180d[baseline_shift_summ$n_all_180d < obfuscation_value] <- NA
        baseline_shift_summ$n_all_90d[baseline_shift_summ$n_all_90d < obfuscation_value] <- NA
        baseline_shift_summ$n_shift_365d[baseline_shift_summ$n_shift_365d < obfuscation_value] <- NA
        baseline_shift_summ$n_shift_180d[baseline_shift_summ$n_shift_180d < obfuscation_value] <- NA
        baseline_shift_summ$n_shift_90d[baseline_shift_summ$n_shift_90d < obfuscation_value] <- NA
    }
    write.csv(baseline_shift_summ,file=file.path(dir.output, paste0(currSiteId, "_BaselineShift_Counts.csv")),row.names=FALSE)
    cat("\nCreating baseline shift counts (150% sCr definition)...")
    baseline_shift150 <- dplyr::bind_rows(aki_only_index_baseline_shift,no_aki_index_baseline_shift) %>% dplyr::distinct()
    baseline_shift150 <- merge(baseline_shift150,first_baseline,by="patient_id",all.x=T)
    baseline_shift150 <- baseline_shift150 %>% dplyr::group_by(patient_id) %>% dplyr::mutate(cr_90d = ifelse(cr_90d_day >= 80 & cr_90d_day <= 100,cr_90d,NA_real_),cr_180d = ifelse(cr_180d_day >= 150 & cr_180d_day <= 210,cr_180d,NA_real_),cr_365d = ifelse(cr_365d_day >= 330 & cr_365d_day <= 390,cr_365d,NA_real_)) %>% dplyr::filter(!is.na(cr_365d) | !is.na(cr_180d) | !is.na(cr_90d)) %>% dplyr::ungroup()
    baseline_shift150 <- baseline_shift150 %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio_365d = cr_365d/first_baseline_cr, ratio_180d = cr_180d/first_baseline_cr,ratio_90d = cr_90d/first_baseline_cr) %>% dplyr::select(patient_id,severe,ratio_365d,ratio_180d,ratio_90d)
    baseline_shift150_summ <- baseline_shift150 %>% dplyr::group_by(severe) %>% dplyr::summarise(n_all=dplyr::n(),n_all_365d = sum(ratio_365d > 0, na.rm=T), n_all_180d = sum(ratio_180d > 0, na.rm = T), n_all_90d = sum(ratio_90d > 0, na.rm = T), n_shift_365d = sum(ratio_365d >= 1.5, na.rm=T),n_shift_180d = sum(ratio_180d >= 1.5, na.rm=T),n_shift_90d = sum(ratio_90d >= 1.5, na.rm=T)) %>% dplyr::ungroup()
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        baseline_shift150_summ$n_all[baseline_shift150_summ$n_all < obfuscation_value] <- NA
        baseline_shift150_summ$n_all_365d[baseline_shift150_summ$n_all_365d < obfuscation_value] <- NA
        baseline_shift150_summ$n_all_180d[baseline_shift150_summ$n_all_180d < obfuscation_value] <- NA
        baseline_shift150_summ$n_all_90d[baseline_shift150_summ$n_all_90d < obfuscation_value] <- NA
        baseline_shift150_summ$n_shift_365d[baseline_shift150_summ$n_shift_365d < obfuscation_value] <- NA
        baseline_shift150_summ$n_shift_180d[baseline_shift150_summ$n_shift_180d < obfuscation_value] <- NA
        baseline_shift150_summ$n_shift_90d[baseline_shift150_summ$n_shift_90d < obfuscation_value] <- NA
    }
    write.csv(baseline_shift150_summ,file=file.path(dir.output, paste0(currSiteId, "_BaselineShift150_Counts.csv")),row.names=FALSE)
    
    cat("\nCreating baseline shift counts (125% sCr definition) using baseline sCr defined as lowest sCr from (-365) days to last day of index admission...")
    baseline_shift_index_admit_lowest <- dplyr::bind_rows(aki_only_index_baseline_shift,no_aki_index_baseline_shift) %>% dplyr::distinct()
    baseline_shift_index_admit_lowest <- merge(baseline_shift_index_admit_lowest,baseline_cr_index_admit,by="patient_id",all.x=T)
    baseline_shift_index_admit_lowest <- baseline_shift_index_admit_lowest %>% dplyr::group_by(patient_id) %>% dplyr::mutate(cr_90d = ifelse(cr_90d_day >= 80 & cr_90d_day <= 100,cr_90d,NA_real_),cr_180d = ifelse(cr_180d_day >= 150 & cr_180d_day <= 210,cr_180d,NA_real_),cr_365d = ifelse(cr_365d_day >= 330 & cr_365d_day <= 390,cr_365d,NA_real_)) %>% dplyr::filter(!is.na(cr_365d) | !is.na(cr_180d) | !is.na(cr_90d)) %>% dplyr::ungroup()
    baseline_shift_index_admit_lowest <- baseline_shift_index_admit_lowest %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio_365d = cr_365d/min_cr_index_admit, ratio_180d = cr_180d/min_cr_index_admit,ratio_90d = cr_90d/min_cr_index_admit) %>% dplyr::select(patient_id,severe,ratio_365d,ratio_180d,ratio_90d)
    baseline_shift_index_admit_lowest_summ <- baseline_shift_index_admit_lowest %>% dplyr::group_by(severe) %>% dplyr::summarise(n_all=dplyr::n(),n_all_365d = sum(ratio_365d > 0, na.rm=T), n_all_180d = sum(ratio_180d > 0, na.rm = T), n_all_90d = sum(ratio_90d > 0, na.rm = T), n_shift_365d = sum(ratio_365d >= 1.25, na.rm=T),n_shift_180d = sum(ratio_180d >= 1.25, na.rm=T),n_shift_90d = sum(ratio_90d >= 1.25, na.rm=T)) %>% dplyr::ungroup()
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        baseline_shift_index_admit_lowest_summ$n_all[baseline_shift_index_admit_lowest_summ$n_all < obfuscation_value] <- NA
        baseline_shift_index_admit_lowest_summ$n_all_365d[baseline_shift_index_admit_lowest_summ$n_all_365d < obfuscation_value] <- NA
        baseline_shift_index_admit_lowest_summ$n_all_180d[baseline_shift_index_admit_lowest_summ$n_all_180d < obfuscation_value] <- NA
        baseline_shift_index_admit_lowest_summ$n_all_90d[baseline_shift_index_admit_lowest_summ$n_all_90d < obfuscation_value] <- NA
        baseline_shift_index_admit_lowest_summ$n_shift_365d[baseline_shift_index_admit_lowest_summ$n_shift_365d < obfuscation_value] <- NA
        baseline_shift_index_admit_lowest_summ$n_shift_180d[baseline_shift_index_admit_lowest_summ$n_shift_180d < obfuscation_value] <- NA
        baseline_shift_index_admit_lowest_summ$n_shift_90d[baseline_shift_index_admit_lowest_summ$n_shift_90d < obfuscation_value] <- NA
    }
    write.csv(baseline_shift_index_admit_lowest_summ,file=file.path(dir.output, paste0(currSiteId, "_BaselineShift_Counts_IndexAdmitLowestCr.csv")),row.names=FALSE)
    cat("\nCreating baseline shift counts (150% sCr definition) using baseline sCr defined as lowest sCr from (-365) days to last day of index admission...")
    baseline_shift150_index_admit_lowest <- dplyr::bind_rows(aki_only_index_baseline_shift,no_aki_index_baseline_shift) %>% dplyr::distinct()
    baseline_shift150_index_admit_lowest  <- merge(baseline_shift150_index_admit_lowest,baseline_cr_index_admit,by="patient_id",all.x=T)
    baseline_shift150_index_admit_lowest <- baseline_shift150_index_admit_lowest %>% dplyr::group_by(patient_id) %>% dplyr::mutate(cr_90d = ifelse(cr_90d_day >= 80 & cr_90d_day <= 100,cr_90d,NA_real_),cr_180d = ifelse(cr_180d_day >= 150 & cr_180d_day <= 210,cr_180d,NA_real_),cr_365d = ifelse(cr_365d_day >= 330 & cr_365d_day <= 390,cr_365d,NA_real_)) %>% dplyr::filter(!is.na(cr_365d) | !is.na(cr_180d) | !is.na(cr_90d)) %>% dplyr::ungroup()
    baseline_shift150_index_admit_lowest  <- baseline_shift150_index_admit_lowest %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio_365d = cr_365d/min_cr_index_admit, ratio_180d = cr_180d/min_cr_index_admit,ratio_90d = cr_90d/min_cr_index_admit) %>% dplyr::select(patient_id,severe,ratio_365d,ratio_180d,ratio_90d)
    baseline_shift150_index_admit_lowest_summ <- baseline_shift150_index_admit_lowest %>% dplyr::group_by(severe) %>% dplyr::summarise(n_all=dplyr::n(),n_all_365d = sum(ratio_365d > 0, na.rm=T), n_all_180d = sum(ratio_180d > 0, na.rm = T), n_all_90d = sum(ratio_90d > 0, na.rm = T), n_shift_365d = sum(ratio_365d >= 1.5, na.rm=T),n_shift_180d = sum(ratio_180d >= 1.5, na.rm=T),n_shift_90d = sum(ratio_90d >= 1.5, na.rm=T)) %>% dplyr::ungroup()
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        baseline_shift150_index_admit_lowest_summ$n_all[baseline_shift150_index_admit_lowest_summ$n_all < obfuscation_value] <- NA
        baseline_shift150_index_admit_lowest_summ$n_all_365d[baseline_shift150_index_admit_lowest_summ$n_all_365d < obfuscation_value] <- NA
        baseline_shift150_index_admit_lowest_summ$n_all_180d[baseline_shift150_index_admit_lowest_summ$n_all_180d < obfuscation_value] <- NA
        baseline_shift150_index_admit_lowest_summ$n_all_90d[baseline_shift150_index_admit_lowest_summ$n_all_90d < obfuscation_value] <- NA
        
        baseline_shift150_index_admit_lowest_summ$n_shift_365d[baseline_shift150_index_admit_lowest_summ$n_shift_365d < obfuscation_value] <- NA
        baseline_shift150_index_admit_lowest_summ$n_shift_180d[baseline_shift150_index_admit_lowest_summ$n_shift_180d < obfuscation_value] <- NA
        baseline_shift150_index_admit_lowest_summ$n_shift_90d[baseline_shift150_index_admit_lowest_summ$n_shift_90d < obfuscation_value] <- NA
    }
    write.csv(baseline_shift150_index_admit_lowest_summ,file=file.path(dir.output, paste0(currSiteId, "_BaselineShift150_Counts_IndexAdmitLowestCr.csv")),row.names=FALSE)
    
    invisible(gc())
    
    # =======================================================================================
    # Figure 1(a): Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
    # =======================================================================================

    # Create KDIGO Stage table for demographics table
    # kdigo_grade <- peak_aki_vs_non_aki %>% dplyr::filter(time_from_peak == 0) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki_kdigo_stage = ifelse(aki == 0,0,ifelse(ratio < 2,1,ifelse(ratio<3,2,3)))) %>% dplyr::ungroup()
    # kdigo_grade <- kdigo_grade %>% dplyr::select(patient_id,aki_kdigo_stage)
    cat("\nExtracting KDIGO stages for index AKI episodes in the index admission only (does NOT extract for AKIs in future admissions)...")
    # kdigo_grade <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission,na.rm=TRUE)) %>% dplyr::filter(aki_kdigo_final == max(aki_kdigo_final,na.rm=TRUE)) %>% dplyr::select(patient_id,aki_kdigo_final) %>% dplyr::distinct(patient_id,.keep_all = TRUE) %>% dplyr::ungroup()
    kdigo_grade <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(first_admit == 1) %>% dplyr::filter(days_since_admission == min(days_since_admission,na.rm=TRUE)) %>% dplyr::filter(aki_kdigo_final == max(aki_kdigo_final,na.rm=TRUE)) %>% dplyr::select(patient_id,aki_kdigo_final) %>% dplyr::distinct(patient_id,.keep_all = TRUE) %>% dplyr::ungroup()
    colnames(kdigo_grade)[2] <- "aki_kdigo_grade"
    kdigo_grade <- merge(kdigo_grade,demographics_filt[,c("patient_id","siteid")],by="patient_id",all=TRUE)
    kdigo_grade$aki_kdigo_grade[is.na(kdigo_grade$aki_kdigo_grade)] <- 0 # the NAs are likely due to the AKIs occuring after first admission, so assign them to the KDIGO 0 group
    kdigo_grade <- kdigo_grade %>% dplyr::select(patient_id,aki_kdigo_grade)
    
    # In the rare circumstance where a non-AKI patient gets mistakenly detected as a KDIGO 3 AKI patient (when the sCr gradually increases to >= 4mg/dL but does not fulfill AKI criteria), this line will help reset the non-AKI patients to KDIGO 0
    kdigo_grade$aki_kdigo_grade[(kdigo_grade$patient_id %in% unique(no_aki_index$patient_id))] <- 0
    
    
    # If patient received RRT in the index admission, the AKI is likely to be KDIGO 3 grade but this would not necessarily
    # be captured in the serum creatinine trends
    # Hence will need to check back against the rrt_index_admit list and replace the grading with 3
    tryCatch({
        if(!is.null(rrt_index_admit) | length(rrt_index_admit) > 0) {
            kdigo_grade$aki_kdigo_grade[(kdigo_grade$patient_id %in% rrt_index_admit) & (kdigo_grade$patient_id %in% aki_only_index$patient_id)] <- 3
        }
    }, error = function(e) {
        cat("\nIssue with re-assigning KDIGO grade for AKI patients who underwent RRT.\nError message:\n",e,"\n")
    })
    
    
    # Create the CKD table
    cat("\nNow checking if CKD is a comorbid present in the cohort...")
    ckd_list <- NULL
    ckd_present <- ("ckd" %in% comorbid_list)
    if(isTRUE(ckd_present)) {
        ckd_list <- comorbid %>% dplyr::select(patient_id,ckd)
    }
    cat("\nCKD present: ",ckd_present)
        
    cat("\nGenerating sCr plots for ALL patients...")
    cr_graphs <- generate_cr_graphs(currSiteId,peak_trend,is_obfuscated,obfuscation_value,kdigo_grade,ckd_present,ckd_list,labs_cr_all,use_index_baseline = FALSE,baseline_cr_index_admit,custom_output,custom_output_dir)
    cat("\nDone!")
    
    cat("\nNow doing the same for ALL patients - but using (-365,last_day_of_index_admit) as timeframe for baseline sCr...")
    cr_graphs_baselineindex <- generate_cr_graphs(currSiteId,peak_trend,is_obfuscated,obfuscation_value,kdigo_grade,ckd_present,ckd_list,labs_cr_all,use_index_baseline = TRUE,baseline_cr_index_admit,custom_output,custom_output_dir)
    cat("\nDone!")
    
    peak_trend_severe <- cr_graphs$peak_trend_severe
    adm_to_aki_cr <- cr_graphs$adm_to_aki_cr
    aki_from_start <- cr_graphs$aki_from_start
    
    ## ===============================================================================
    ## Initializing analysis for cirrhotic patients + plotting serum creatinine graphs
    ## ===============================================================================
    
    cirrhosis_present <- ("cld" %in% comorbid_list)
    inr_loinc <- c("6301-6","34714-6","38875-1","46418-0","52129-4","61189-7","72281-9","92891-1")
    sodium_loinc <- c("2947-0","32717-1","39792-7","41657-8","39791-9","2951-2","77139-4")
    meld_analysis_valid <- FALSE
    labs_meld_admission <- NULL
    meld_labs_day_cutoff <- 3
    tryCatch({
        if(isTRUE(cirrhosis_present)) {
            cat("\n=====================================\n")
            cat("\nProcessing serum creatinine graphs for cirrhotic patients.\nFirst checking if admission MELD labs are available.")
            dir.create(file.path(dir.output,paste0(currSiteId,"_Cirrhosis")))
            cirrhosis_list <- comorbid %>% dplyr::select(patient_id,cld) %>% dplyr::filter(cld == 1)
            labs_cirrhosis <- observations[observations$patient_id %in% cirrhosis_list$patient_id,] %>% dplyr::filter(concept_type == "LAB-LOINC")
            labs_meld_first72h <- labs_cirrhosis %>% dplyr::filter(days_since_admission >= 0 & days_since_admission <= meld_labs_day_cutoff)
            labs_list <- unlist(unique(labs_meld_first72h$concept_code))
            inr_present <- (length(intersect(inr_loinc,labs_list)) > 0)
            bil_present <- ('1975-2' %in% labs_list)
            cr_present <- ('2160-0' %in% labs_list)
            sodium_present <- (length(intersect(sodium_loinc,labs_list)) > 0)
            if(isTRUE(sodium_present)){
                cat("\nSodium values present for MELD score correction.")
            }
            
            if(isTRUE(inr_present) & isTRUE(bil_present) & isTRUE(cr_present)) {
                cat("\nFound INR + bilirubin + creatinine values in first 72h in the Observations table. Will proceed with sub-group analysis for hepatorenal syndrome.")
                meld_analysis_valid <- TRUE
                cat("\nExtracting and binning INR")
                labs_meld_list <- NULL
                try({
                    labs_inr <- labs_meld_first72h[labs_meld_first72h$concept_code %in% inr_loinc,-c(4,5)]
                    labs_inr <- labs_inr %>% dplyr::filter(days_since_admission >= 0)
                    labs_inr <- labs_inr %>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_inr = mean(na.omit(value)),min_inr = min(na.omit(value)),max_inr = max(na.omit(value)),first_inr = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_inr")
                })
                cat("\nExtracting and binning bilirubin")
                try({
                    labs_bil <- labs_meld_first72h[labs_meld_first72h$concept_code == '1975-2',-c(4,5)]
                    labs_bil <- labs_bil %>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_bil = mean(na.omit(value)),min_bil = min(na.omit(value)),max_bil = max(na.omit(value)),first_bil = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_bil")
                })
                cat("\nExtracting and binning Cr")
                try({
                    labs_cr <- labs_meld_first72h[labs_meld_first72h$concept_code == '2160-0',-c(4,5)]
                    labs_cr <- labs_cr%>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_cr = mean(na.omit(value)),min_cr = min(na.omit(value)),max_cr = max(na.omit(value)),first_cr = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_cr")
                })
                cat("\nExtracting and binning AST")
                try({
                    labs_ast <- labs_meld_first72h[labs_meld_first72h$concept_code == '1920-8',-c(4,5)]
                    labs_ast <- labs_ast %>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_ast = mean(na.omit(value)),min_ast = min(na.omit(value)),max_ast = max(na.omit(value)),first_ast = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_ast")
                })
                cat("\nExtracting and binning ALT")
                try({
                    labs_alt <- labs_meld_first72h[labs_meld_first72h$concept_code == '1742-6',-c(4,5)]
                    labs_alt <- labs_alt %>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_alt = mean(na.omit(value)),min_alt = min(na.omit(value)),max_alt = max(na.omit(value)),first_alt = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_alt")
                })
                cat("\nExtracting and binning albumin")
                try({
                    labs_alb <- labs_meld_first72h[labs_meld_first72h$concept_code == '1751-7',-c(4,5)]
                    labs_alb <- labs_alb %>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_alb = mean(na.omit(value)),min_alb = min(na.omit(value)),max_alb = max(na.omit(value)),first_alb = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_alb")
                })
                
                if(isTRUE(sodium_present)) {
                    cat("\nAdding in sodium data")
                    labs_na <- labs_meld_first72h[labs_meld_first72h$concept_code %in% sodium_loinc,-c(4,5)]
                    labs_na <- labs_na %>% dplyr::group_by(patient_id) %>% dplyr::summarise(mean_na = mean(na.omit(value)),min_na = min(na.omit(value)),max_na = max(na.omit(value)),first_na = dplyr::first(na.omit(value)))
                    labs_meld_list <- c(labs_meld_list,"labs_na")
                }
                cat("Valid variables: ",labs_meld_list)
                
                cat("\nMerging all tables with binned data\n")
                
                labs_meld <- mget(labs_meld_list) %>% purrr::reduce(dplyr::full_join,by="patient_id") %>% dplyr::distinct()

                severe_label <- data.table::data.table(c(1,2,3,4),c("Non-severe, no AKI","Non-severe, AKI","Severe, no AKI","Severe, AKI"))
                colnames(severe_label) <- c("severe","severe_label")
                
                cat("\nFiltering for patients with complete MELD labs on admission.\n")
                labs_meld <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld_valid = dplyr::if_else((is.na(first_inr) | is.na(first_bil) | is.na(first_cr)), 0,1)) %>% dplyr::ungroup() %>% dplyr::filter(meld_valid == 1)
                
                if(length(unique(labs_meld$patient_id)) >= obfuscation_value) {
                    cat("\nCalculating MELD score...")
                    if(isTRUE(sodium_present)) {
                        labs_meld <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld = FourCePhase2.1AKI:::meld_score(bil = first_bil,inr = first_inr,sCr = first_cr,Na = min_na)) %>% dplyr::ungroup()
                    } else {
                        labs_meld <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld = FourCePhase2.1AKI:::meld_score(bil = first_bil,inr = first_inr,sCr = first_cr)) %>% dplyr::ungroup()
                    }
                    cat("\nExtracting admission MELD score...")
                    if(isTRUE(sodium_present)) {
                        labs_meld_admission <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld_admit_severe = dplyr::if_else(meld >= 20,1,0),sodium_valid = dplyr::if_else(is.na(min_na),0,1)) %>% dplyr::ungroup()
                    } else {
                        labs_meld_admission <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld_admit_severe = dplyr::if_else(meld >= 20,1,0)) %>% dplyr::ungroup()
                    }
                    labs_meld_admission <- labs_meld_admission[!is.na(labs_meld_admission$meld_admit_severe),]
                    meld_severe_list <- labs_meld_admission %>% dplyr::select(patient_id,meld,meld_admit_severe) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
                    
                    # Final headers
                    # labs_meld:
                    # patient_id, day_bin, mean/min/max/first of labs, meld
                    #
                    # labs_meld_admission:
                    # patient_id, day_bin, mean/min/max/first of labs, meld (integer score), meld_admit_severe (0/1)
                    # Possible that some patient_ids may not be inside this if there are no labs in days 0-3 (i.e. no day_bin == "[0,3]")
                    cat("\nNow creating creatinine graphs with MELD scores")
                    peak_trend_meld <- peak_trend[peak_trend$patient_id %in% meld_severe_list$patient_id,]
                    peak_trend_meld <- merge(peak_trend_meld,meld_severe_list,by="patient_id",all.x=TRUE)
                    
                    peak_trend_meld <- peak_trend_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = dplyr::if_else(severe == 1 | severe == 3,0,1)) %>% dplyr::mutate(meld_admit_severe = dplyr::if_else(aki == 1,dplyr::if_else(meld_admit_severe == 1,4,2),dplyr::if_else(meld_admit_severe == 1,3,1)))
                    peak_cr_meld_summ <- peak_trend_meld %>% dplyr::group_by(meld_admit_severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
                    
                    
                    meld_label <- data.table::data.table(c(1,2,3,4),c("MELD < 20, no AKI","MELD < 20, AKI","MELD >= 20, no AKI","MELD >= 20, AKI"))
                    colnames(meld_label) <- c("meld_admit_severe","meld_severe_label")
                    peak_cr_meld_summ <- merge(peak_cr_meld_summ,meld_label,by="meld_admit_severe",all.x=TRUE)
                    if(isTRUE(is_obfuscated)) {
                        cat("\nObfuscating the MELD AKI graphs...")
                        peak_cr_meld_summ <- peak_cr_meld_summ[peak_cr_meld_summ$n >= obfuscation_value,]
                    }
                    peak_cr_meld_timeplot <- ggplot2::ggplot(peak_cr_meld_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=meld_severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(meld_severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("MELD < 20, AKI"="#bc3c29","MELD < 20, no AKI"="#0072b5","MELD >= 20, AKI" = "#e18727","MELD >= 20, no AKI"="#20854e")) + ggplot2::theme_minimal()
                    ggplot2::ggsave(filename=file.path(dir.output,paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_PeakCr_MELD_AKI.png")),plot=peak_cr_meld_timeplot,width=12,height=9,units="cm")
                    
                    # Plot from start of admission to 30 days post-peak AKI (if no AKI, then from peak Cr)
                    adm_meld_cr <- labs_cr_all[labs_cr_all$patient_id %in% meld_severe_list$patient_id,]
                    adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(peak_cr_time == min(peak_cr_time)) %>% dplyr::distinct() %>% dplyr::ungroup()
                    adm_meld_cr <- adm_meld_cr[order(adm_meld_cr$patient_id,adm_meld_cr$days_since_admission),]
                    adm_meld_cr <- merge(adm_meld_cr,meld_severe_list,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
                    adm_meld_cr$meld_admit_severe[is.na(adm_meld_cr$meld_admit_severe)] <- 0
                    adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_365d,min_cr_retro_365d)) %>% dplyr::ungroup()
                    adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
                    adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = dplyr::if_else(severe == 1 | severe == 3,0,1)) %>% dplyr::mutate(meld_admit_severe = dplyr::if_else(aki == 1,dplyr::if_else(meld_admit_severe == 1,4,2),dplyr::if_else(meld_admit_severe == 1,3,1)))
                    adm_meld_summ <- adm_meld_cr %>% dplyr::group_by(meld_admit_severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
                    adm_meld_summ <- merge(adm_meld_summ,meld_label,by="meld_admit_severe",all.x=TRUE)
                    if(isTRUE(is_obfuscated)) {
                        # adm_to_aki_summ <- adm_to_aki_summ %>% dplyr::filter(n >= obfuscation_value)
                        cat("\nObfuscating the admission to AKI graphs...")
                        adm_meld_summ <- adm_meld_summ[adm_meld_summ$n >= obfuscation_value,]
                    }
                    adm_meld_timeplot <- ggplot2::ggplot(adm_meld_summ[which(adm_meld_summ$days_since_admission <= 30 & adm_meld_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=meld_severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(meld_severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(0,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("MELD < 20, AKI"="#bc3c29","MELD < 20, no AKI"="#0072b5","MELD >= 20, AKI" = "#e18727","MELD >= 20, no AKI"="#20854e")) + ggplot2::theme_minimal()
                    ggplot2::ggsave(filename=file.path(dir.output,paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_CrfromAdm_MELD_AKI.png")),plot=adm_meld_timeplot,width=12,height=9,units="cm")
                    
                    # Plot from start of AKI to 30 days later 
                    aki_start_meld <- labs_cr_all[labs_cr_all$patient_id %in% meld_severe_list$patient_id,]
                    # Uncomment the following line to restrict analysis to AKI patients only
                    #aki_start_meld <- aki_start_meld[aki_start_meld$severe %in% c(2,4,5),]
                    aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_start = days_since_admission - aki_start) %>% dplyr::ungroup()
                    aki_start_meld <- aki_start_meld[order(aki_start_meld$patient_id,aki_start_meld$days_since_admission),] %>% dplyr::distinct()
                    aki_start_meld <- merge(aki_start_meld,meld_severe_list,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
                    aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_365d,min_cr_retro_365d)) %>% dplyr::ungroup()
                    aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
                    #aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
                    aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = dplyr::if_else(severe == 1 | severe == 3,0,1)) %>% dplyr::mutate(meld_admit_severe = dplyr::if_else(aki == 1,dplyr::if_else(meld_admit_severe == 1,4,2),dplyr::if_else(meld_admit_severe == 1,3,1)))
                    aki_start_meld_summ <- aki_start_meld %>% dplyr::group_by(meld_admit_severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
                    aki_start_meld_summ <- merge(aki_start_meld_summ,meld_label,by="meld_admit_severe",all.x=TRUE)
                    if(isTRUE(is_obfuscated)) {
                        cat("\nObfuscating the start of AKI graphs...")
                        # aki_start_meld_summ <- aki_start_meld_summ %>% dplyr::filter(n >= obfuscation_value)
                        aki_start_meld_summ <- aki_start_meld_summ[aki_start_meld_summ$n >= obfuscation_value,]
                    }
                    aki_start_meld_timeplot <- ggplot2::ggplot(aki_start_meld_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=meld_severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(meld_severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("MELD < 20, AKI"="#bc3c29","MELD < 20, no AKI"="#0072b5","MELD >= 20, AKI" = "#e18727","MELD >= 20, no AKI"="#20854e")) + ggplot2::theme_minimal()
                    ggplot2::ggsave(file=file.path(dir.output,paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_CrFromStart_MELD_AKI.png")),plot=aki_start_meld_timeplot,width=12,height=9,units="cm")
                    
                } else {
                    meld_analysis_valid = FALSE
                }
                
                cat("\nCreating tables for graphs sorted by COVID-19 severity...")
                peak_trend_severe_cld <- peak_trend_severe[peak_trend_severe$patient_id %in% cirrhosis_list$patient_id,]
                peak_cr_cld_summ <- peak_trend_severe_cld %>% dplyr::group_by(severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
                peak_cr_cld_summ <- merge(peak_cr_cld_summ,severe_label,by="severe",all.x=TRUE)
                if(isTRUE(is_obfuscated)) {
                    # peak_cr_summ <- peak_cr_summ %>% dplyr::filter(n >= obfuscation_value)
                    cat("\nObfuscating the AKI with severity graphs...")
                    peak_cr_cld_summ <- peak_cr_cld_summ[peak_cr_cld_summ$n >= obfuscation_value,]
                }
                peak_cr_cld_timeplot <- ggplot2::ggplot(peak_cr_cld_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
                ggplot2::ggsave(filename=file.path(dir.output,paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_PeakCr_AKI_CLD_only.png")),plot=peak_cr_cld_timeplot,width=12,height=9,units="cm")
                cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of cirrhotic patients (with severity) should have been generated.")
                
                adm_to_aki_cld_summ <- adm_to_aki_cr[adm_to_aki_cr$patient_id %in% cirrhosis_list$patient_id,]
                adm_to_aki_cld_summ <- adm_to_aki_cld_summ %>% dplyr::group_by(severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
                adm_to_aki_cld_summ <- merge(adm_to_aki_cld_summ,severe_label,by="severe",all.x=TRUE)
                if(isTRUE(is_obfuscated)) {
                    # adm_to_aki_cld_summ <- adm_to_aki_cld_summ %>% dplyr::filter(n >= obfuscation_value)
                    cat("\nObfuscating the admission to AKI graphs...")
                    adm_to_aki_cld_summ <- adm_to_aki_cld_summ[adm_to_aki_cld_summ$n >= obfuscation_value,]
                }
                adm_to_aki_cld_timeplot <- ggplot2::ggplot(adm_to_aki_cld_summ[which(adm_to_aki_cld_summ$days_since_admission <= 30 & adm_to_aki_cld_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(0,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
                ggplot2::ggsave(filename=file.path(dir.output,paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_CrfromAdm_AKI_CLD_only.png")),plot=adm_to_aki_cld_timeplot,width=12,height=9,units="cm")
                
                cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of cirrhotic patients, plotted from first day of admission, should have been generated.")
                
                aki_start_cld_summ <- aki_from_start[aki_from_start$patient_id %in% cirrhosis_list$patient_id,]
                aki_start_cld_summ <- aki_start_cld_summ %>% dplyr::group_by(severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
                aki_start_cld_summ <- merge(aki_start_cld_summ,severe_label,by="severe",all.x=TRUE)
                if(isTRUE(is_obfuscated)) {
                    cat("\nObfuscating the start of AKI graphs...")
                    # aki_start_cld_summ <- aki_start_cld_summ %>% dplyr::filter(n >= obfuscation_value)
                    aki_start_cld_summ <- aki_start_cld_summ[aki_start_cld_summ$n >= obfuscation_value,]
                }
                aki_start_cld_timeplot <- ggplot2::ggplot(aki_start_cld_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
                ggplot2::ggsave(file=file.path(dir.output,paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_CrFromStart_AKI_CLD_only.png")),plot=aki_start_cld_timeplot,width=12,height=9,units="cm")
                cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from start of AKI/creatinine increase, should have been generated.")
            }
        }
        
        if(isTRUE(meld_analysis_valid)) {
            try({save(peak_cr_cld_summ,peak_cr_meld_summ,adm_to_aki_cld_summ,adm_meld_summ,aki_start_cld_summ,aki_start_meld_summ,file=file.path(dir.output, paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_MELD_CLD_graphs.rda")),compress="bzip2")})
        } else {
            try({save(peak_cr_cld_summ,adm_to_aki_cld_summ,aki_start_cld_summ,file=file.path(dir.output, paste0(currSiteId,"_Cirrhosis"), paste0(currSiteId,"_CLD_graphs.rda")),compress="bzip2")})
        }
    }, error = function(e) {
        message("\nEncountered error while processing cirrhosis graphs.\nSpecific message:\n",e)
    })
    
    
    # =====================
    # Demographics Table
    # =====================
    cat("\nNow generating the demographics table...")
    
    demog_full_analysis <- generate_demog_files(currSiteId,demographics_filt,labs_aki_summ,comorbid,comorbid_list,kdigo_grade,coaga_present,coagb_present,covid19antiviral_present,remdesivir_present,acei_present,arb_present,med_coaga_new,med_coagb_new,med_covid19_new,med_acearb_chronic,cirrhosis_valid = isTRUE(cirrhosis_present),is_obfuscated,obfuscation_value,aki_list_first_admit,custom_output,custom_output_dir,ckd_staging,earliest_cr)

    demog_time_to_event <- demog_full_analysis$demog_toe
    demog_list <- demog_full_analysis$demog_var_list

    ## ====================================
    ## PART 4: Time To Event Analysis
    ## ====================================
    
    cat("\nNow proceeding to time-to-event analysis. Ensure that the survival and survminer packages are installed in your R environment.")
    
    ckd_staging <- ckd_staging[,c("patient_id","ckd_stage")]
    
    main_analysis <- run_time_to_event_analysis(currSiteId,peak_trend,aki_index,labs_aki_summ,demographics,demog_time_to_event,demog_list,comorbid,comorbid_list,kdigo_grade,ckd_present,coaga_present,coagb_present,covid19antiviral_present,remdesivir_present,acei_present,arb_present,med_coaga_new,med_coagb_new,med_covid19_new,med_acearb_chronic,earliest_cr,ckd_staging,is_obfuscated,obfuscation_value,restrict_models,factor_cutoff,custom_output,custom_output_dir)
    
    aki_index_recovery <- main_analysis$aki_index_recovery
    aki_index_death <- main_analysis$aki_index_death
    med_recovery_list <- main_analysis$med_recovery_list
    comorbid_recovery_list <- main_analysis$comorbid_recovery_list
    demog_recovery_list <- main_analysis$demog_recovery_list
    earliest_cr_recovery_list <- main_analysis$earliest_cr_recovery_list
    med_death_list <- main_analysis$med_death_list
    comorbid_death_list <- main_analysis$comorbid_death_list
    demog_death_list <- main_analysis$demog_death_list
    earliest_cr_death_list <- main_analysis$earliest_cr_death_list
    
    var_list_recovery_all <- main_analysis$var_list_recovery_all
    var_list_death_all <- main_analysis$var_list_death_all
    
    cat("\nProceeding to non-CKD subgroup analysis...\n")
    
    nonckd_analysis <- run_time_to_event_analysis_nonckd(currSiteId,
                                                         aki_index_recovery,aki_index_death,
                                                         med_recovery_list, comorbid_recovery_list, demog_recovery_list,earliest_cr_recovery_list,
                                                         med_death_list, comorbid_death_list, demog_death_list,earliest_cr_death_list,
                                                         is_obfuscated,obfuscation_value,
                                                         restrict_models, factor_cutoff,
                                                         custom_output,custom_output_dir,ckd_present)
    
    aki_index_nonckd <- nonckd_analysis$aki_index_nonckd
    aki_index_nonckd_akionly <- nonckd_analysis$aki_index_nonckd_akionly
    var_list_new_ckd_nonckd_akionly <- nonckd_analysis$var_list_new_ckd_nonckd_akionly
    var_list_recovery_nonckd_akionly <- nonckd_analysis$var_list_recovery_nonckd_akionly
    var_list_death_nonckd_akionly <- nonckd_analysis$var_list_death_nonckd_akionly
    var_list_death_nonckd <- nonckd_analysis$var_list_death_nonckd
    var_list_new_ckd_nonckd <- nonckd_analysis$var_list_new_ckd_nonckd
    
    if(isTRUE(ckd_present)) {
      cat("\nCKD present: ",ckd_present,"\nProceeding to CKD subgroup analysis")
      
      ckd_analysis <- run_time_to_event_analysis_ckdonly(currSiteId,aki_index_recovery,aki_index_death,
                                                         var_list_recovery_all, var_list_death_all,
                                                         is_obfuscated,obfuscation_value,
                                                         restrict_models, factor_cutoff,
                                                         custom_output,custom_output_dir)
      aki_index_ckdonly <- ckd_analysis$aki_index_ckdonly
      aki_index_ckdonly_akionly <- ckd_analysis$aki_index_ckdonly_akionly
      var_list_recovery_ckdonly_akionly <- ckd_analysis$var_list_recovery_ckdonly_akionly
      var_list_death_ckdonly_akionly <- ckd_analysis$var_list_death_ckdonly_akionly
      var_list_death_ckdonly <- ckd_analysis$var_list_death_ckdonly
    } else {
      cat("\nNoted that no CKD comorbidities are present in the cohort. Will proceed to skip subgroup analyses since all patients are non-CKD.\n")
      aki_index_ckdonly <- NULL
      aki_index_ckdonly_akionly <- NULL
      var_list_recovery_ckdonly_akionly <- NULL
      var_list_death_ckdonly_akionly <- NULL
      var_list_death_ckdonly <- NULL
    }
    
    cat("\nNow generating the life tables for competing risks...\n")
    generate_life_table_competing_risk(currSiteId,aki_index_recovery,aki_index_nonckd,aki_index_nonckd_akionly,aki_index_ckdonly_akionly,custom_output,custom_output_dir,ckd_present)
    cat("\nNow running the competing risk analyses...\n")
    run_cic_analysis(currSiteId,aki_index_recovery,aki_index_nonckd,aki_index_nonckd_akionly,aki_index_ckdonly_akionly,
                     var_list_recovery_all,
                     var_list_new_ckd_nonckd_akionly,
                     var_list_recovery_nonckd_akionly,
                     var_list_new_ckd_nonckd,
                     var_list_recovery_ckdonly_akionly,
                     custom_output,custom_output_dir,ckd_present)
    # aki_index_recovery <- main_analysis$aki_index_recovery
    # med_recovery_list <- main_analysis$med_recovery_list
    # comorbid_recovery_list <- main_analysis$comorbid_recovery_list
    # aki_index_death <- main_analysis$aki_index_death
    
    # aki_index_recovery150 <- run_recovery_analysis_150(currSiteId,peak_trend,aki_index,labs_aki_summ,demographics,demog_time_to_event,demog_list,comorbid,comorbid_list,kdigo_grade,ckd_present,coaga_present,coagb_present,covid19antiviral_present,remdesivir_present,acei_present,arb_present,med_coaga_new,med_coagb_new,med_covid19_new,med_acearb_chronic,is_obfuscated,obfuscation_value,restrict_models,factor_cutoff,custom_output,custom_output_dir)


    cat("\n\n========================================\nAnalysis complete.\n")
    if(isTRUE(print_rrt_surrogate)) {
        message("Final reminder:\nYou have opted to print patient-level data files for the purposes of manual chart review for RRT procedures.")
        message("The files are located at ",dir.output, " and are named as:")
        message("a) DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv")
        message("b) DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv")
        message("\nEnsure that these files are removed before uploading to GitHub!")
    }
    sink()
    if(isTRUE(debug_on)) {
        sink(type="message")
        cat("\n\nPlease check the error log in the project output folder for details.\n")
    }
    closeAllConnections()
    invisible(gc())
}

