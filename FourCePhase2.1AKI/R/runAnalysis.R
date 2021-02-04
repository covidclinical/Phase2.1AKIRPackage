
#' Runs the analytic workflow for the AKI project
#'
#' @keywords 4CE
#' @export

runAnalysis <- function(is_obfuscated=TRUE,obfuscation_value=3,factor_cutoff = 5) {

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

    ## get the site identifier associated with the files stored in the /4ceData/Input directory that 
    ## is mounted to the container
    currSiteId = FourCePhase2.1Data::getSiteId()

    ## run the quality control
    FourCePhase2.1Data::runQC(currSiteId)

    ## DO NOT CHANGE ANYTHING ABOVE THIS LINE

    ## To Do: implement analytic workflow, saving results to a site-specific 
    ## file to be sent to the coordinating site later via submitAnalysis()

    ## ========================================
    ## PART 1: Read in Data Tables
    ## ========================================
    message("Reading in Input/LocalPatientSummary.csv and Input/LocalPatientObservations.csv...")
    demographics <- read.csv("Input/LocalPatientSummary.csv",na.strings = '1/1/1900')
    observations <- read.csv("Input/LocalPatientObservations.csv",na.strings = '-999')
    message("Using data() to load internal tables (may result in errors)...")
    data(thromb_ref)
    data(comorbid_ref)
    data(thromb_icd9_ref)
    data(comorbid_icd9_ref)
    data(thromb_icd10_ref)
    data(comorbid_icd10_ref)
    
    message("Transforming the Summary and Observations tables to generate tables for demographics, diagnoses, procedures")
    # first generate a unique ID for each patient
    demographics <- demographics %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    # course <- course %>% mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    observations <- observations %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
    
    # From this point on, we will be using our custom-generated patient_id as a unique patient identifier
    # Reorder the columns in each table to bring patient_id to the first column and remove patient_num
    demographics <- demographics %>% dplyr::select(patient_id,siteid,admission_date,days_since_admission,last_discharge_date,still_in_hospital,severe_date,severe,death_date,deceased,sex,age_group,race,race_collected)
    # course <- course[,c(8,1,3:7)]
    observations <- observations %>% dplyr::select(patient_id,siteid,days_since_admission,concept_type,concept_code,value)
    
    # Generate a diagnosis table
    diagnosis <- observations[observations$concept_type %in% c("DIAG-ICD9","DIAG-ICD10"),-6]
    colnames(diagnosis) <- c("patient_id","siteid","days_since_admission","concept_type","icd_code")
    
    # Generate a procedures table
    procedures <- observations[observations$concept_type %in% c("PROC-ICD9", "PROC-ICD10"),-6]
    colnames(procedures) <- c("patient_id","siteid","days_since_admission","icd_version","procedure_code")
    
    demographics_filt <- demographics %>% dplyr::mutate(time_to_severe = ifelse(severe == 1, as.numeric(as.Date(severe_date) - as.Date(admission_date)),NA))
    demographics_filt <- demographics_filt %>% dplyr::mutate(time_to_death = ifelse(deceased == 1, as.numeric(as.Date(death_date) - as.Date(admission_date)),NA))
    demographics_filt <- demographics_filt %>% dplyr::mutate(length_stay = ifelse(still_in_hospital==1,days_since_admission,as.numeric(as.Date(last_discharge_date) - as.Date(admission_date))))
    
    # Reorder the columns to be more readable
    demographics_filt <- demographics_filt %>% dplyr::select(patient_id,siteid,sex,age_group,race,length_stay,severe,time_to_severe,deceased,time_to_death)
    
    # Final headers for demographics_filt
    # patient_id  siteid sex age_group race  length_stay severe  time_to_severe  deceased  time_to_death
    
    # Comorbidities & Prothrombotic Events
    # ======================================
    message("Now creating tables for comorbidities and admission diagnoses...")
    diag_icd9 <- diagnosis[diagnosis$concept_type == "DIAG-ICD9",]
    diag_icd10 <- diagnosis[diagnosis$concept_type == "DIAG-ICD10",]
    
    message("Creating Comorbidities table...")
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
    # 
    # # Filter prothrombotic events from all diagnoses
    # thromb_icd9 <- diag_icd9[diag_icd9$icd_code %in% thromb_icd9_ref$icd_code,]
    # thromb_icd10 <- diag_icd10[diag_icd10$icd_code %in% thromb_icd10_ref$icd_code,]
    # # Filter prothrombotic diagnoses to be restricted to diagnosis codes made after -15days
    # thromb_icd9 <- thromb_icd9[thromb_icd9$days_since_admission >= -15,-c(2,4)]
    # thromb_icd10 <- thromb_icd10[thromb_icd10$days_since_admission >= -15,-c(2,4)]
    # # Map the prothrombotic codes - store day diagnosed
    # thromb_icd9 <- merge(thromb_icd9,thromb_icd9_ref,by="icd_code",all.x=TRUE)
    # thromb_icd10 <- merge(thromb_icd10,thromb_icd10_ref,by="icd_code",all.x=TRUE)
    # thromb_diag <- rbind(thromb_icd9,thromb_icd10)
    # thromb_diag <- thromb_diag %>% dplyr::select(patient_id,type) %>% dplyr::mutate(present = 1) %>% dplyr::distinct()
    # thromb_diag <- thromb_diag %>% tidyr::spread(type,present)
    # thromb_diag[is.na(thromb_diag)] <- NA
    
    # Final headers for comorbid table
    # Note: order of columns may depend on the overall characteristics of your patient population
    # patient_id	ihd,htn,dm,asthma,copd,bronchiectasis,ild,ckd,pe,dvt,cancer
    # Values stored are in binary, 1 = present, 0 = absent
    
    # Final headers for thromb_diag table
    # Note: order of columns may depend on the overall characteristics of your patient population
    # patient_id	dvt,vt,pe,mi
    # Values stored are in binary, 1 = present, 0 = absent
    
    # Time to Intubation
    # ====================
    message("Creating table for intubation...")
    # (1) Determine intubation from procedure code
    intubation_code <- c("0BH13EZ","0BH17EZ","0BH18EZ","0B21XEZ","5A09357","5A09358","5A09359","5A0935B","5A0935Z","5A09457","5A09458","5A09459","5A0945B","5A0945Z","5A09557","5A09558","5A09559","5A0955B","5A0955Z","96.7","96.04","96.70","96.71","96.72")
    intubation <- procedures[procedures$procedure_code %in% intubation_code,-c(2,4)]
    intubation <- intubation[,c(1,2)]
    intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission),]
    intubation <- intubation[!duplicated(intubation$patient_id),]
    # (2) In some cases intubation may not be coded as a procedure. Hence surrogate way is to determine if 
    #     patient had been diagnosed with ARDS and/or VAP
    vap_ards_codes <- c("J80","J95.851","518.82","997.31")
    vap_ards_diag <- diagnosis[diagnosis$icd_code %in% vap_ards_codes,]
    intubation <- rbind(vap_ards_diag[,c(1,3)],intubation)
    intubation <- intubation[order(intubation$patient_id,intubation$days_since_admission),]
    intubation <- intubation[!duplicated(intubation$patient_id),]
    
    # Final table headers:
    # patient_id	days_since_admission
    # days_since_admission = time to first intubation event
    
    # Time to RRT
    # ==================
    message("Creating table for RRT...")
    rrt_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","3E1.M39Z","549.8","399.5")
    #hd_code <- c("5A1D70Z","5A1D80Z","5A1D90Z","399.5")
    #pd_code <- c("3E1.M39Z","549.8")
    rrt <- procedures[procedures$procedure_code %in% rrt_code,-c(2,4)]
    rrt <- rrt[,c(1,2)]
    rrt <- rrt[order(rrt$patient_id,rrt$days_since_admission),]
    rrt <- rrt[!duplicated(rrt$patient_id),]
    
    # Generate list of patients already on RRT prior to admission
    # This list can be used to exclude ESRF patients in subsequent analyses
    patients_already_rrt <- rrt$patient_id[rrt$days_since_admission < 0]
    
    # Generate list of patients who were only initiated on RRT for the first time ever
    # during admission for COVID-19
    rrt_new <- rrt[!(rrt$patient_id %in% patients_already_rrt),]
    
    # Final table headers for rrt and rrt_new:
    # patient_id	days_since_admission
    # days_since_admission = time to first RRT event
    
    # ====================
    # PART 2: AKI Detection Code
    # ====================
    # For the purposes of generating a table of AKI events, we will have to create a new data.table
    # with the serum creatinine levels.
    # We will then use this table, labs_cr_aki, to generate a summary table containing details of 
    # each AKI event and the post-AKI recovery labs
    message("Now proceeding to AKI detection code:")
    message("Extracting serum creatinine values...")
    # Extract serum Cr levels
    labs_cr_aki <- observations[observations$concept_code == '2160-0',] #LOINC code for Cr 2160-0
    # Remove unnecessary columns
    labs_cr_aki <- labs_cr_aki[,-c(4,5)]
    # Filter for labs >= -90 days
    labs_cr_aki <- labs_cr_aki %>% dplyr::filter(days_since_admission >= -60)
    
    # Generate separate demographics table for patients who do not have any sCr values fulfilling 
    # the above (e.g. all the labs are before t= -90days or patient has no sCr value)
    message("Removing patients who do not have any serum creatinine values...")
    pts_valid_cr <- unique(labs_cr_aki$patient_id)
    demog_no_cr <- demographics_filt[!(demographics_filt$patient_id %in% pts_valid_cr),]
    demographics_filt <- demographics_filt[demographics_filt$patient_id %in% pts_valid_cr,]
    # There are two possible scenarios which we have to consider when detecting each AKI event:
    # (1) AKI occurs after admission
    #   - easy to detect with the formal KDIGO definition as we only need to use older data points
    #   - we can use an approach similar to what is used in the MIMIC-III data-set: use a rolling
    #     time frame and detect the lowest Cr value to use as our baseline
    # (2) Patient presents with an AKI at the time of admission 
    #   - this makes it harder for us to determine what is the true baseline especially with limited
    #     longitudinal data
    #   - Hence one way to get around this is to generate a "retrospective" baseline (i.e. look into
    #     future Cr values and find the minimum in a rolling timeframe) and use this as a surrogate
    #     baseline
    #   - Such a workaround is the only feasible way of dealing with missing retrospective data though
    #     it is likely that we may miss quite a number of true AKIs using this method
    
    # Scenario (1): AKI occurs during admission
    # Find minimum Cr level in a rolling 90 day timeframe
    message("Generating minimum Cr level in the past 90 days")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission)-90,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_90d = RcppRoll::roll_min(value,91,fill=NA,na.rm=TRUE,align="right")) %>% dplyr::filter(!is.na(value))
    # Find minimum Cr level in a rolling 2 day timeframe (48h)
    message("Generating minimum Cr level in the past 48h")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission)-2,max(days_since_admission)))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_48h = RcppRoll::roll_min(value,3,fill=NA,na.rm=TRUE,align="right")) %>% dplyr::filter(!is.na(value))
    
    # Scenario (2): Patient presents with an AKI already on board
    # Find minimum Cr level in a rolling 7 day timeframe
    message("Generating minimum Cr level 7 days in the future")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+7))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(min_cr_retro_7day = RcppRoll::roll_min(value,8,fill=NA,na.rm=TRUE,align="left")) %>% dplyr::filter(!is.na(value))
    # Find minimum Cr level in a rolling 2 day timeframe (48h)
    message("Generating minimum Cr level 48h in the future")
    labs_cr_aki <- data.table::data.table(labs_cr_aki,key=c("patient_id","days_since_admission"))
    labs_cr_aki <- labs_cr_aki[data.table::CJ(unique(patient_id),seq(min(days_since_admission),max(days_since_admission)+2))] # temporarily fill in null rows for missing days - the minimum rolling code only works for consecutive data
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
    
    # Generate sCr levels at +7d (cr_7d) and +90d (cr_90d) timepoints (for determining post-AKI recovery, AKD)
    labs_cr_aki <- data.table::setDT(labs_cr_aki)[,':='(cr_7d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+7,incbounds = TRUE)],1),cr_90d = tail(labs_cr_aki$value[labs_cr_aki$patient_id==patient_id][data.table::between(labs_cr_aki$days_since_admission[labs_cr_aki$patient_id==patient_id],days_since_admission,days_since_admission+90,incbounds = TRUE)],1)),by=c('patient_id','days_since_admission')][]
    
    # At this point, our table has these headers:
    # patient_id  siteid  days_since_admission  value min_cr_90d min_cr_48h  min_cr_retro_7day min_cr_48h_retro  cr_7d cr_90d
    
    # Now we have to start grading AKI severity at each time point
    # This approach is similar to how the MIMIC-III dataset generates AKI severity
    # Generate two columns using both the formal KDIGO AKI definition and the modified retrospective AKI definition
    message("Generating KDIGO severity grades for each serum Cr value")
    labs_cr_aki$aki_kdigo <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade)
    labs_cr_aki$aki_kdigo_retro <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::aki_kdigo_grade_retro)
    labs_cr_aki <- labs_cr_aki %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(aki_kdigo_final = max(aki_kdigo,aki_kdigo_retro))
    
    # Generate two columns grading AKD severity at 7d and 90d (grade 0B/C is coded as 0.5)
    labs_cr_aki$akd_7d <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::akd_grade_7d)
    labs_cr_aki$akd_90d <- apply(labs_cr_aki,1,FourCePhase2.1AKI:::akd_grade_90d)
    
    # Now we are going to generate the start days of each AKI
    labs_cr_aki_tmp <- labs_cr_aki
    labs_cr_aki_tmp$valid = 1
    
    # Find the day of the minimum Cr used for grading AKIs (taken as baseline)
    message("Now finding the day at which the minimum serum Cr is achieved")
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
    message("Now identifying maxima points of serum Cr")
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id,days_since_admission) %>% dplyr::mutate(min_cr_7d_final = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::mutate(delta_cr = value - min_cr_7d_final)
    
    # Use the largest delta_cr to find the peak of each AKI
    labs_cr_aki_delta_maxima <- labs_cr_aki_tmp4 %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr %in% delta_cr[which.peaks(delta_cr,decreasing=FALSE)])
    labs_cr_aki_delta_maxima$delta_is_max = 1
    labs_cr_aki_delta_maxima <- labs_cr_aki_delta_maxima %>% dplyr::rename(delta_maxima = delta_cr) %>% dplyr::select(patient_id,days_since_admission,delta_maxima,delta_is_max)
    labs_cr_aki_tmp4 <- merge(labs_cr_aki_tmp4,labs_cr_aki_delta_maxima,by=c("patient_id","days_since_admission"),all.x=TRUE)
    
    # Generate a separate table (for reference) of all creatinine peaks not fulfilling KDIGO AKI criteria
    message("Now generating a separate table for serum Cr peaks which do not reach AKI definitions")
    labs_cr_nonaki <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final == 0,]
    labs_cr_nonaki[is.na(labs_cr_nonaki)] <- 0
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$delta_is_max > 0,]
    labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
    
    # Filter for KDIGO grades > 0
    message("Now generating tables of all AKI events")
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final > 0,]
    labs_cr_aki_tmp4[is.na(labs_cr_aki_tmp4)] <- 0
    # Filter for maxima of delta_cr (which should give us the peaks)
    labs_cr_aki_tmp4 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$delta_is_max > 0,]
    
    # Filter and reorder columns to generate our final table of all AKI events
    labs_aki_summ <- labs_cr_aki_tmp4 %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
    
    labs_aki_summ <- labs_aki_summ %>% dplyr::distinct(patient_id,days_since_admission,.keep_all=TRUE)
    
    # Final headers for labs_aki_summ:
    # patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
    # days_since_admission - time at which peak Cr is achieved
    # day_min - time at which Cr begins to rise
    
    # Generate the highest Cr peak for non-AKI peaks detected
    labs_nonaki_summ <- labs_cr_nonaki %>% dplyr::group_by(patient_id) %>% dplyr::slice(which.max(delta_cr))
    
    # We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset 
    # (2) how long before/after disease severity
    # These tables will help in segregating the populations for analysis later
    message("Generating intermediate tables to determine if AKIs occured before or after severe COVID-19 onset")
    severe_time <- demographics_filt %>% dplyr::select(patient_id,severe,time_to_severe)
    labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = ifelse(!is.na(time_to_severe), time_to_severe - day_min,NA))
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0))
    # Final headers for labs_aki_severe:
    # patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d  severe time_to_severe  severe_to_aki severe_before_aki
    
    labs_nonaki_severe <- merge(labs_nonaki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_nonaki_severe$severe_to_aki <- NA 
    labs_nonaki_severe$severe_before_aki <- 1 
    
    ## Save the generated AKI tables for future reference / debugging (note: these will NOT be uploaded!!)
    #write.csv(labs_aki_summ,"PatientAKIEvents.csv",row.names=FALSE)
    #write.csv(labs_aki_severe,"PatientAKIEvents_Severe.csv",row.names=FALSE)
    
    # =============
    # Medications
    # =============
    message("Now generating the medications tables:")
    medications <- observations[observations$concept_type == "MED-CLASS",]
    medications <- medications[,-c(4,6)]
    medications <- medications %>% dplyr::arrange(patient_id,days_since_admission,concept_code)
    # Use 15 days as the cutoff for chronic medications
    message("Generating table for chronic medications...")
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
    message("Generating table for prior ACE-inhibitor (ACEI) or angiotensin receptor blocker (ARB) use...")
    acei_present = ("old_ACEI" %in% colnames(med_chronic))
    arb_present = ("old_ARB" %in% colnames(med_chronic))
    if(acei_present == TRUE && arb_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ACEI,old_ARB)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ACEI + old_ARB > 0,1,0))
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    } else if (acei_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ACEI)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ACEI > 0,1,0))
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    } else if (arb_present == TRUE) {
        med_acearb_chronic <- med_chronic %>% dplyr::select(patient_id,old_ARB)
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::group_by(patient_id) %>% dplyr::mutate(acei_arb_preexposure = ifelse(old_ARB > 0,1,0))
        med_acearb_chronic <- med_acearb_chronic %>% dplyr::select(patient_id,acei_arb_preexposure)
    }
    
    # For simplicity of initial analysis, we will use the earliest date where each new medication class is
    # used. However, if we are to incorporate a recurrent neural network model to account for temporal changes
    # in medications, this approach cannot be used.
    message("Generating table for new medications started during the admission...")
    med_new <- med_new[!duplicated(med_new[,c(1,2,4)]),]
    med_new <- med_new[,c(1,4,3)]
    colnames(med_new)[3] <- "start_day"
    # Compute the new medications as a time difference from the start of AKI
    # If the value is < 0 (and presumably > -365) then the medication was initiated before the start of AKI
    # Otherwise it means the medication had been started after the peak Cr had been achieved
    # This will give insight into which medications may be potentially nephrotoxic
    aki_start_time <- labs_aki_summ[,c(1,5)]
    med_new_aki <- merge(med_new,aki_start_time,by="patient_id",all.x=TRUE)
    med_new_aki <- med_new_aki[!is.na(med_new_aki$day_min),]
    med_new_aki <- med_new_aki %>% dplyr::distinct()
    med_new_aki <- med_new_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(offset_aki = start_day - day_min) %>% dplyr::filter(offset_aki == min(offset_aki))
    # Re-code whether medication was given before AKI - 1 = yes, 0 = no
    med_new_aki <- med_new_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(med_before_aki = ifelse(offset_aki <=0,1,0))
    med_new_aki <- med_new_aki[,c(1,2,6)]
    med_new_aki <- med_new_aki %>% tidyr::spread(concept_code,med_before_aki)
    #med_new_aki[is.na(med_new_aki)] <- -999
    
    # Generate another table with the start date of the new medications
    med_new <- med_new %>% tidyr::spread(concept_code,start_day)
    #med_new[is.na(med_new)] <- -999
    
    # Generate simplified table for determining who were started on COAGB near admission
    message("Generating table for COAGB use during the admission...")
    med_coagb_new <- med_new %>% dplyr::select(patient_id,COAGB)
    med_coagb_new$COAGB[med_coagb_new$COAGB < -15] <- 0
    med_coagb_new$COAGB[med_coagb_new$COAGB >= -15] <- 1
    
    # Generate simplified table for determining who were started on novel antivirals
    covid19antiviral_present = ("COVIDVIRAL" %in% colnames(med_new))
    if(covid19antiviral_present == TRUE){
        message("Generating table for experimental COVID-19 treatment...")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,COVIDVIRAL)
        # Uncomment following line to select for patients who were started on novel antivirals less than 72h from admission
        #med_covid19_new <- med_covid19_new[med_covid19_new$COVIDVIRAL <= 3,]
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(COVIDVIRAL >= 0,1,0))
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx)
    }
    
    # =====================
    # Demographics Table
    # =====================
    message("Now generating the demographics table...")
    demog_summ <- demographics_filt %>% dplyr::select(patient_id,sex,age_group,race,severe,deceased,time_to_severe,time_to_death)
    demog_summ <- merge(demog_summ,comorbid,by="patient_id",all.x=TRUE)
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$aki <- 0
    demog_summ$aki[demog_summ$patient_id %in% labs_aki_summ$patient_id] <- 1
    demog_summ$severe <- factor(demog_summ$severe,levels=c(0,1),labels=c("Non-severe","Severe"))
    demog_summ$deceased <- factor(demog_summ$deceased,levels=c(0,1),labels=c("Alive","Deceased"))
    demog_summ$aki <- factor(demog_summ$aki,levels=c(0,1),labels=c("No AKI","AKI"))
    demog_summ[comorbid_list] <- lapply(demog_summ[comorbid_list],factor)
    table_one_vars <- c("sex","age_group","race","severe","deceased","time_to_severe","time_to_death",comorbid_list)
    table_one <- tableone::CreateTableOne(data=demog_summ,vars=table_one_vars,strata="aki")
    export_table_one <- print(table_one,showAllLevels=TRUE,formatOptions=list(big.mark=","))
    write.csv(export_table_one,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne.csv")))
    capture.output(summary(table_one),file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne_Missingness.txt")))
    message("TableOne with patient demographics should have been generated in CSV files at this point. Check for any errors.")
    
    demog_time_to_event_tmp <- demog_summ[,c("sex","age_group","race")]
    demog_time_to_event_tmp <- data.table::as.data.table(lapply(demog_time_to_event_tmp,factor))
    demog_time_to_event_tmp <- data.table::as.data.table(demog_time_to_event_tmp)[,sapply(demog_time_to_event_tmp,function(col) nlevels(col) > 1),with=FALSE] 
    demog_list <- colnames(demog_time_to_event_tmp)
    demog_time_to_event <- demog_summ[,c("patient_id",demog_list)]
    
    ## ==================================================================================
    ## PART 3: Serum Creatinine Trends - Plots against Time from Peak Serum Creatinine
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
    
    message("Now generating tables for use for plotting normalised serum Cr values against time")
    message("First creating table for AKI patients only...")
    # First, dplyr::filter the labs_aki_severe table to show only the index AKI episodes
    aki_only_index <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::ungroup()
    # patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_90d	min_cr_48h	min_cr_retro_7day	min_cr_48h_retro	min_cr_7d_final	cr_7d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	severe  time_to_severe	severe_to_aki	severe_before_aki
    
    # Generate the patient list including (1) severity indices from this dplyr::filtered table (2) day of peak Cr
    # severe - 2 = never severe, 4 = severe, AKI
    aki_only_index <- aki_only_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = dplyr::if_else(!is.na(time_to_severe),4,2)) %>% dplyr::ungroup()
    aki_only_index <- aki_only_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
    colnames(aki_only_index)[2] <- "peak_cr_time"
    colnames(aki_only_index)[4] <- "aki_start"
    # Headers of aki_only_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki
    message("Now creating list of non-AKI patients...")
    no_aki_list <- demographics_filt %>% dplyr::select(patient_id,severe)
    no_aki_list <- no_aki_list[!(no_aki_list$patient_id %in% aki_only_index$patient_id),]
    no_aki_list <- no_aki_list %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = dplyr::if_else(severe == 1,3,1))
    
    labs_nonaki_summ <- labs_nonaki_summ[labs_nonaki_summ$patient_id %in% no_aki_list$patient_id,]
    labs_nonaki_severe <- labs_nonaki_severe[labs_nonaki_severe$patient_id %in% no_aki_list$patient_id,]
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$patient_id %in% no_aki_list$patient_id,]
    
    message("Creating table for non-AKI patients...")
    # Create a non-AKI equivalent for aki_only_index - except that this takes the largest delta_cr (and the earliest occurence of such a delta_cr)
    no_aki_index <- labs_nonaki_severe %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr == max(delta_cr)) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::ungroup()
    no_aki_index <- no_aki_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
    no_aki_index <- no_aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(severe == 1,3,1))
    colnames(no_aki_index)[2] <- "peak_cr_time"
    colnames(no_aki_index)[4] <- "aki_start"
    #no_aki_index$severe_to_aki <- -999
    message("Creating final table of peak Cr for all patients...")
    aki_index <- dplyr::bind_rows(aki_only_index,no_aki_index)
    if(covid19antiviral_present == TRUE) {
        message("Adding on temporal information for experimental COVID-19 antivirals...")
        aki_index <- merge(aki_index,med_covid19_new,by="patient_id",all.x=TRUE)
        aki_index$covid_rx[is.na(aki_index$covid_rx)] <- 0
        aki_index <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covidrx_grp = ifelse(severe <= 2, ifelse(covid_rx == 0,1,2),ifelse(covid_rx == 0,3,4))) %>% dplyr::ungroup()
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,covidrx_grp)
    } else {
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki)
    }
    # Headers of aki_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki  covidrx_grp
    aki_index <- aki_index %>% dplyr::arrange(patient_id,peak_cr_time,desc(severe))%>% dplyr::distinct(patient_id,peak_cr_time,.keep_all = TRUE)
    message("Final table aki_index created.")
    # Uncomment the following line to remove patients who were previously on RRT prior to admission
    # aki_index <- aki_index[!(aki_index$patient_id %in% patients_already_rrt),]
    message("Now generating table of normalised serum Cr...")
    # Create a common labs_cr_all table containing the serum Cr values, the severity groupings and anti-viral groupings
    labs_cr_aki_tmp <- labs_cr_aki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_90d,min_cr_retro_7day)
    labs_cr_nonaki_tmp <- labs_cr_nonaki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_90d,min_cr_retro_7day)
    labs_cr_all <- dplyr::bind_rows(labs_cr_aki_tmp,labs_cr_nonaki_tmp)
    labs_cr_all <- merge(labs_cr_all,aki_index,by="patient_id",all.x=TRUE)
    
    # Now, generate a table containing lab values with timepoints calculated from time of peak cr
    peak_trend <- labs_cr_all
    if(covid19antiviral_present == TRUE) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day)
    }
    # patient_id  severe  covidrx_grp  days_since_admission  peak_cr_time  value min_cr_90d min_cr_retro_7day
    
    # Calculate the day from peak Cr
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_peak = days_since_admission - peak_cr_time) %>% dplyr::ungroup()
    ## dplyr::filter this table for Cr values that fall within the 7 day window
    # peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(between(time_from_peak,0,7)) %>% dplyr::ungroup()
    # Normalise to baseline values used for AKI calculation
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    
    # In the event longitudinal data becomes very long, we will create a column where the very first baseline Cr for the index AKI is generated for each patient
    first_baseline <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak == 0) %>% dplyr::select(patient_id,baseline_cr) %>% dplyr::distinct()
    colnames(first_baseline)[2] <- "first_baseline_cr"
    peak_trend <- merge(peak_trend,first_baseline,by="patient_id",all.x=TRUE)
    if(covid19antiviral_present == TRUE) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day,time_from_peak,baseline_cr,first_baseline_cr)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day,time_from_peak,baseline_cr,first_baseline_cr)
    }
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/first_baseline_cr) %>% dplyr::ungroup()
    #peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    message("Final table of peak Cr for all patients - peak_trend - created.")
    # peak_trend will now be a common table to plot from the dplyr::selected AKI peak
    
    # =======================================================================================
    # Figure 1(a): Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
    # =======================================================================================
    
    # First create a plot of the creatinine trends of AKI vs non-AKI patients from the day of the first AKI peak (or the highest Cr peak for non-AKI patients)
    peak_aki_vs_non_aki <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio)
    colnames(peak_aki_vs_non_aki) <- c("patient_id","aki","time_from_peak","ratio")
    peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = ifelse((aki == 2 | aki == 4 | aki == 5),1,0))
    peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki %>% dplyr::group_by(aki,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    aki_label <- data.table::data.table(c(0,1),c("Non-AKI","AKI"))
    colnames(aki_label) <- c("aki","aki_label")
    peak_aki_vs_non_aki_summ <- merge(peak_aki_vs_non_aki_summ,aki_label,by="aki",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki_summ %>% dplyr::group_by(aki,time_from_peak) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(peak_aki_vs_non_aki_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrFromPeak_AKI_vs_NonAKI.csv")),row.names=FALSE)
    peak_aki_vs_non_aki_timeplot <- ggplot2::ggplot(peak_aki_vs_non_aki_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=aki_label))+ggplot2::geom_line(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio, color=factor(aki_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "AKI Group") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("AKI"="#bc3c29","Non-AKI"="#0072b5")) + ggplot2::theme_minimal()
    print(peak_aki_vs_non_aki_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_AKI_vs_NonAKI.png")),plot=peak_aki_vs_non_aki_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients should have been generated.")

    # Now, derive our first table peak_trend_severe to compare across the different severity groups
    peak_trend_severe <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio)
    # Headers: patient_id  severe  time_from_peak  ratio
    peak_trend_severe <- peak_trend_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    
    # Calculate mean and SD each for severe and non-severe groups
    peak_cr_summ <- peak_trend_severe %>% dplyr::group_by(severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    severe_label <- data.table::data.table(c(1,2,3,4),c("Non-severe, no AKI","Non-severe, AKI","Severe, no AKI","Severe, AKI"))
    colnames(severe_label) <- c("severe","severe_label")
    peak_cr_summ <- merge(peak_cr_summ,severe_label,by="severe",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        peak_cr_summ <- peak_cr_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(peak_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrfromPeak_Severe_AKI.csv")),row.names=FALSE)
    # Plot the graphs
    peak_cr_timeplot <- ggplot2::ggplot(peak_cr_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
    print(peak_cr_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromPeak_Severe_AKI.png")),plot=peak_cr_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients (with severity) should have been generated.")
    
    # Plot from start of admission to 30 days post-peak AKI (if no AKI, then from peak Cr)
    adm_to_aki_cr <- labs_cr_all
    # adm_to_aki_cr$peak_cr_time[is.na(adm_to_aki_cr$peak_cr_time)] <- 0
    #adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(dplyr::between(days_since_admission,0,peak_cr_time+30)) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(peak_cr_time == min(peak_cr_time)) %>% dplyr::distinct() %>% dplyr::filter(days_since_admission >= 0 & days_since_admission <= peak_cr_time+30) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    adm_to_aki_summ <- adm_to_aki_cr %>% dplyr::group_by(severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    adm_to_aki_summ <- merge(adm_to_aki_summ,severe_label,by="severe",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        adm_to_aki_summ <- adm_to_aki_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(peak_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.csv")),row.names=FALSE)
    adm_to_aki_timeplot <- ggplot2::ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 30),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(0,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
    print(adm_to_aki_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.png")),plot=adm_to_aki_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from first day of admission, should have been generated.")
    
    # Plot from start of AKI to 30 days later 
    
    aki_30d_cr <- labs_cr_all
    # Uncomment the following line to restrict analysis to AKI patients only
    #aki_30d_cr <- aki_30d_cr[aki_30d_cr$severe %in% c(2,4,5),]
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_start = days_since_admission - aki_start) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr[order(aki_30d_cr$patient_id,aki_30d_cr$days_since_admission),]
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    aki_30d_cr_summ <- aki_30d_cr %>% dplyr::group_by(severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    aki_30d_cr_summ <- merge(aki_30d_cr_summ,severe_label,by="severe",all.x=TRUE)
    if(is_obfuscated==TRUE) {
        aki_30d_cr_summ <- aki_30d_cr_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
    }
    write.csv(aki_30d_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.csv")),row.names=FALSE)
    aki_30d_cr_timeplot <- ggplot2::ggplot(aki_30d_cr_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
    print(aki_30d_cr_timeplot)
    ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.png")),plot=aki_30d_cr_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from start of AKI/creatinine increase, should have been generated.")
    
    # ================================================================================================================================
    # Figure 1(b) Comparing serum creatinine trends of severe and non-severe patients, with or without remdesivir/lopinavir+ritonavir
    # ================================================================================================================================
    if(covid19antiviral_present == TRUE) {
        message("Now preparing plots for serum Cr trends in experimental COVID-19 treatment...")
        # Plotting from peak
        peak_trend_covidviral <- peak_trend %>% dplyr::select(patient_id,covidrx_grp,time_from_peak,ratio)
        peak_cr_covidviral_summ <- peak_trend_covidviral %>% dplyr::group_by(covidrx_grp,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
        covidrx_grp_label <- data.table::data.table(c(1:4),c("Non-severe, no novel COVID-19 treatment","Non-severe, with novel COVID-19 treatment","Severe, no novel COVID-19 treatment","Severe, with novel COVID-19 treatment"))
        colnames(covidrx_grp_label) <- c("covidrx_grp","covidrx_label")
        peak_cr_covidviral_summ <- merge(peak_cr_covidviral_summ,covidrx_grp_label,by="covidrx_grp",all.x=TRUE)
        if(is_obfuscated==TRUE) {
            peak_cr_covidviral_summ  <- peak_cr_covidviral_summ  %>% dplyr::group_by(covidrx_grp) %>% dplyr::filter(n >= obfuscation_value)
        }
        write.csv(peak_cr_covidviral_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromPeak_CovidViral.csv")),row.names=FALSE)
        peak_cr_covidviral_timeplot <- ggplot2::ggplot(peak_cr_covidviral_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=covidrx_label))+ggplot2::geom_line(ggplot2::aes(color = factor(covidrx_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(covidrx_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color=factor(covidrx_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity + COVID-19 Treatment") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, no novel COVID-19 treatment"="#bc3c29","Non-severe, with novel COVID-19 treatment"="#0072b5","Severe, no novel COVID-19 treatment" = "#e18727","Severe, with novel COVID-19 treatment"="#20854e")) + ggplot2::theme_minimal()
        print(peak_cr_covidviral_timeplot)
        ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromPeak_CovidViral.png")),plot=peak_cr_covidviral_timeplot,width=20,height=9,units="cm")
        
        # Plotting from initiation of novel anti-virals
        cr_from_covidrx_trend <- merge(peak_trend,med_covid19_new_date,by="patient_id",all.x=TRUE)
        # dplyr::filter out patients who have never received any of the novel antivirals
        #cr_from_covidrx_trend$covid_rx_start[is.na(cr_from_covidrx_trend$covid_rx_start)] <- -999
        cr_from_covidrx_trend$covid_rx[is.na(cr_from_covidrx_trend$covid_rx)] <- 0
        cr_from_covidrx_trend <- cr_from_covidrx_trend[cr_from_covidrx_trend$covid_rx == 1,]
        
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,covid_rx_start,peak_cr_time,value,min_cr_90d,min_cr_retro_7day)
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(severe <= 2,0,1)) %>% dplyr::mutate(covidrx_grp = ifelse((covidrx_grp == 2 || covidrx_grp == 4),1,0)) %>% dplyr::ungroup()
        # Calculate the day from initiation of novel anti-virals
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_covidrx = days_since_admission - ifelse(covidrx_grp == 1,covid_rx_start,0)) %>% dplyr::ungroup()
        # Filter this table for Cr values that fall within the 7 day window
        # cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(between(time_from_covidrx,0,30)) %>% dplyr::ungroup()
        # Normalise to baseline values used for AKI calculation
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
        cr_from_covidrx_trend_severe <- cr_from_covidrx_trend %>% dplyr::select(patient_id,severe,time_from_covidrx,ratio)
        # Headers: patient_id  severe (coded as 0/1)  time_from_covidrx  ratio
        cr_from_covidrx_summ <- cr_from_covidrx_trend_severe %>% dplyr::group_by(severe,time_from_covidrx) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
        covidrx_severe_label <- data.table::data.table(c(0,1),c("Non-severe","Severe"))
        colnames(covidrx_severe_label) <- c("severe","severe_label")
        cr_from_covidrx_summ <- merge(cr_from_covidrx_summ,covidrx_severe_label,by="severe",all.x=TRUE)
        if(is_obfuscated==TRUE) {
            cr_from_covidrx_summ  <- cr_from_covidrx_summ %>% dplyr::group_by(severe) %>% dplyr::filter(n >= obfuscation_value)
        }
        
        write.csv(peak_cr_covidviral_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromCovidRx_Severe.csv")),row.names=FALSE)
        cr_from_covidrx_timeplot <- ggplot2::ggplot(cr_from_covidrx_summ,ggplot2::aes(x=time_from_covidrx,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from COVIDVIRAL Start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe"="#bc3c29","Severe"="#0072b5")) + ggplot2::theme_minimal()
        print(cr_from_covidrx_timeplot)
        ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromCovidRx_Severe.png")),plot=cr_from_covidrx_timeplot,width=20,height=9,units="cm")
        message("If no errors at this point, plots for COVID-19 experimental treatment should be done.")
    }
    
    ## ====================================
    ## PART 4: Time To Event Analysis
    ## ====================================
    
    message("Now proceeding to time-to-event analysis. Ensure that the survival and survminer packages are installed in your R environment.")
    message("Part 1: Time to Recovery:")
    # First generate table where we find the time taken to achieve 1.25x baseline ratio
    labs_cr_recovery <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::filter(time_from_peak >= 0)
    labs_cr_recovery_tmp <- labs_cr_recovery %>% dplyr::group_by(patient_id) %>% tidyr::complete(time_from_peak = tidyr::full_seq(time_from_peak,1)) %>% dplyr::mutate(ratio = zoo::na.fill(ratio,Inf))
    time_to_ratio1.25 <- labs_cr_recovery_tmp %>% split(.$patient_id) %>% purrr::map(~get_day(.$ratio,.$time_from_peak,target=1.25)) %>% purrr::map_df(~dplyr::data_frame(.x),.id='patient_id')
    colnames(time_to_ratio1.25)[2] <- "time_to_ratio1.25"

    labs_aki_summ_index <- labs_aki_summ %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission))
    # Get index AKI grade
    index_aki_grade <- labs_aki_summ_index %>% dplyr::select(patient_id,aki_kdigo_final)
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
    aki_index_recovery <- merge(aki_index_recovery,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE)
    
    # aki_index_recovery <- merge(aki_index_recovery,comorbid,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    # aki_index_recovery[is.na(aki_index_recovery)] <- 0
    # aki_index_recovery[c("severe","aki_kdigo_final",comorbid_list)] <- lapply(aki_index_recovery[c("severe","aki_kdigo_final",comorbid_list)],factor)
    
    message("Doing initial filter for comorbids with more than one factor level.")
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
    aki_index_recovery[is.na(aki_index_recovery)] <- 0
    aki_index_recovery[c("severe","aki_kdigo_final",demog_list,comorbid_recovery_list)] <- lapply(aki_index_recovery[c("severe","aki_kdigo_final",demog_list,comorbid_recovery_list)],factor)
    
    message("Current factor list for recovery: ",paste(comorbid_recovery_list,demog_list,sep=", "))
    message("Filtering factor list down further for CoxPH models...")
    # This portion of code deals with the issue of Cox PH models generating large coefficients and/or overfitting
    # We are going to select for the variables where there are at least 5 occurrences of an event for each factor level
    # We will then modify comorbid_recovery_list to only include variable names where this criteria is fulfilled
    # This does NOT require the aki_index_recovery table to be modified
    recovery_tmp <- aki_index_recovery[,c("patient_id","recover_1.25x",demog_list,comorbid_recovery_list)] %>% as.data.frame()
   
    comorbid_recovery_list_tmp <- vector(mode="list",length=length(comorbid_recovery_list))
    for(i in 1:length(comorbid_recovery_list)) {
        recovery_tmp1 <- recovery_tmp[,c("patient_id",comorbid_recovery_list[i],"recover_1.25x")]
        recovery_tmp2 <- recovery_tmp1 %>% dplyr::count(get(comorbid_recovery_list[i]),recover_1.25x)
        recovery_tmp3 <- recovery_tmp2 %>% dplyr::filter(recover_1.25x == 1)
        # if(min(recovery_tmp3$n) >= factor_cutoff) {
        #     comorbid_recovery_list_tmp[i] <- comorbid_recovery_list[i]
        # }
        if(min(recovery_tmp3$n) >= factor_cutoff & nrow(recovery_tmp3) > 1) {
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
            demog_recovery_list_tmp[i] <- demog_list[i]
        }
    }
    demog_recovery_list <- unlist(demog_recovery_list_tmp[lengths(demog_recovery_list_tmp) > 0L])
    
    message("Final factor list for recovery: ",paste(comorbid_recovery_list,demog_recovery_list,sep=", "))
    
    message("Now proceeding to time-to-Cr recovery analysis...")
    # Now run the actual time-to-event analysis
    recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe")
    recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list),collapse="+")))
    #surv_recover <- survival::Surv(time=aki_index_recovery$time_to_ratio1.25,event=aki_index_recovery$recover_1.25x)
    fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
    plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
    plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
    plot_recover_summ_table <- plot_recover$data.survtable
    write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_recover_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe.png")),plot=print(plot_recover),width=12,height=12,units="cm")
    coxph_recover <- survival::coxph(recoverCoxPHFormula, data=aki_index_recovery)
    coxph_recover_plot <- survminer::ggforest(coxph_recover,data=aki_index_recovery)
    coxph_recover_summ <- summary(coxph_recover) 
    write.csv(coxph_recover_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH.csv")),row.names=TRUE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_CoxPH.png")),plot=print(coxph_recover_plot),width=20,height=20,units="cm")
    
    message("Now proceeding to time-to-death analysis for AKI patients only...")
    deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~",paste(c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list),collapse="+")))
    #surv_death_aki_only<- survival::Surv(time=aki_index_recovery$time_to_death_km,event=aki_index_recovery$deceased)
    fit_death_aki_only <- survminer::surv_fit(deathPlotFormula, data=aki_index_recovery)
    plot_death_aki_only <- survminer::ggsurvplot(fit_death_aki_only,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
    plot_death_aki_only_summ <- survminer::surv_summary(fit_death_aki_only,data=aki_index_recovery)
    plot_death_aki_only_summ_table <- plot_death_aki_only$data.survtable
    write.csv(plot_death_aki_only_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_death_aki_only_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe.png")),plot=print(plot_death_aki_only),width=12,height=12,units="cm")
    coxph_death_aki_only <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery)
    coxph_death_aki_only_plot <- survminer::ggforest(coxph_death_aki_only,data=aki_index_recovery)
    coxph_death_aki_only_summ <- summary(coxph_death_aki_only) 
    write.csv(coxph_death_aki_only_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH.csv")),row.names=TRUE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH.png")),plot=print(coxph_death_aki_only_plot),width=20,height=20,units="cm")
    
    # We now do the same to the time to death analyses:
    # 1) Create a temporary comorbid table and obtain the valid comorbids with more than 1 factor level
    
    message("Part 2: Time to Death analysis...")
    aki_index_death <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(is_aki=ifelse(severe %in% c(2,4,5),1,0)) %>% dplyr::mutate(severe=ifelse(severe %in% c(3,4,5),1,0))
    aki_index_death <- merge(aki_index_death,discharge_day,by="patient_id",all.x=TRUE) # merge in time_to_death_km
    aki_index_death <- merge(aki_index_death,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE) # AKI KDIGO grade
    
    comorbid_death_tmp <- merge(aki_index_death,comorbid,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    comorbid_death_tmp[is.na(comorbid_death_tmp)] <- 0
    comorbid_death_tmp <- comorbid_death_tmp[comorbid_list]
    comorbid_death_tmp <- data.table::as.data.table(lapply(comorbid_death_tmp,factor))
    comorbid_death_tmp <- data.table::as.data.table(comorbid_death_tmp)[,sapply(comorbid_death_tmp,function(col) nlevels(col) > 1),with=FALSE] 
    comorbid_death_list <- colnames(comorbid_death_tmp)
    message("Factor list for Death Analysis before filtering for CoxPH: ",paste(comorbid_death_list,demog_list,sep = ","))
    # 2) Create a new table with the cleaned up comorbids
    aki_index_death <- merge(aki_index_death,comorbid[c("patient_id",comorbid_death_list)],by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    aki_index_death <- merge(aki_index_death,demog_time_to_event,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    aki_index_death[is.na(aki_index_death)] <- 0
    aki_index_death <- aki_index_death %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = dplyr::if_else(!is.na(severe_to_aki),as.integer(min(severe_to_aki)),NA_integer_)) %>% dplyr::distinct()
    aki_index_death[c("severe","aki_kdigo_final","is_aki",demog_list,comorbid_death_list)] <- lapply(aki_index_death[c("severe","aki_kdigo_final","is_aki",demog_list,comorbid_death_list)],factor)

    # 3) This portion of code deals with the issue of Cox PH models generating large coefficients and/or overfitting
    # We are going to select for the variables where there are at least 5 occurrences of an event for each factor level
    # We will then modify comorbid_death_list to only include variable names where this criteria is fulfilled
    # This does NOT require the aki_index_death table to be modified
    death_tmp <- aki_index_death[,c("patient_id","deceased",demog_list,comorbid_death_list)] %>% as.data.frame()
    comorbid_death_list_tmp <- vector(mode="list",length=length(comorbid_death_list))
    for(i in 1:length(comorbid_death_list)) {
        death_tmp1 <- death_tmp[,c("patient_id",comorbid_death_list[i],"deceased")]
        death_tmp2 <- death_tmp1 %>% dplyr::count(get(comorbid_death_list[i]),deceased)
        death_tmp3 <- death_tmp2 %>% dplyr::filter(deceased == 1)
        if(min(death_tmp3$n) >= factor_cutoff & nrow(death_tmp3) > 1) {
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
            demog_death_list_tmp[i] <- demog_list[i]
        }
    }
    demog_death_list <- unlist(demog_death_list_tmp[lengths(demog_death_list_tmp) > 0L])
    
    message("Final factor list for death: ",paste(comorbid_death_list,demog_death_list,sep=", "))

    # 4) Run analysis
    deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ is_aki")
    deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~",paste(c("is_aki","severe",demog_death_list,comorbid_death_list),collapse="+")))
    
    #surv_death <- survival::Surv(time=aki_index_death$time_to_death_km,event=aki_index_death$deceased)
    fit_death <- survminer::surv_fit(deathPlotFormula, data=aki_index_death)
    plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death)
    plot_death_summ_table <- plot_death$data.survtable
    write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_Plot.csv")),row.names=FALSE)
    write.csv(plot_death_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI.png")),plot=print(plot_death),width=12,height=12,units="cm")
    coxph_death <- survival::coxph(deathCoxPHFormula, data=aki_index_death)
    coxph_death_plot <- survminer::ggforest(coxph_death,data=aki_index_death)
    coxph_death_summ <- summary(coxph_death) 
    write.csv(coxph_death_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_CoxPH.csv")),row.names=TRUE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CoxPH.png")),plot=print(coxph_death_plot),width=20,height=20,units="cm")
    
    # # ================================================
    # # Part 3: Extending analyses to Thrombotic Events
    # # ================================================
    # # Questions
    # # 1) Does incidence of AKI in COVID-19 correlate with thrombotic events occurring during COVID-19 illness?
    # # 2) Are there differences in serum Cr trends between patients with thrombotic events and those without?
    # 
    # aki_index_thromb <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(is_aki=ifelse(severe %in% c(3,4,5),1,0)) %>% dplyr::mutate(severe=ifelse(severe %in% c(2,4,5),1,0))
    # aki_index_thromb <- merge(aki_index_thromb,discharge_day,by="patient_id",all.x=TRUE) # merge in time_to_death_km
    # aki_index_thromb <- merge(aki_index_thromb,labs_aki_summ_index[,c(1,17)],by="patient_id",all.x=TRUE) # aki_kdigo_final
    # aki_index_thromb <- merge(aki_index_thromb,time_to_ratio1.25,by="patient_id",all.x=TRUE)
    # aki_index_thromb <- aki_index_thromb %>% dplyr::group_by(patient_id) %>% dplyr::mutate(recover_1.25x = ifelse(is.na(time_to_ratio1.25),0,1))
    # 
    # aki_index_thromb <- merge(aki_index_thromb,comorbid[c("patient_id",thromb_list)],by="patient_id",all.x=TRUE) %>% dplyr::distinct()
    # aki_index_thromb[is.na(aki_index_thromb)] <- 0
    # aki_index_thromb <- aki_index_thromb %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = dplyr::if_else(!is.na(severe_to_aki),as.integer(min(severe_to_aki)),NA_integer_)) %>% dplyr::distinct()
    # aki_index_thromb[c("severe","aki_kdigo_final","is_aki",thromb_list)] <- lapply(aki_index_thromb[c("severe","aki_kdigo_final","is_aki",thromb_list)],factor)
    # 
    # # This portion of code deals with the issue of Cox PH models generating large coefficients and/or overfitting
    # # We are going to select for the variables where there are at least 5 occurrences of an event for each factor level
    # # We will then modify thromb_list to only include variable names where this criteria is fulfilled
    # # This does NOT require the aki_index_thromb table to be modified
    # message("thromb_list for Thromb Analysis before filtering for regression: ",paste(thromb_list,sep = ","))
    # thromb_tmp <- aki_index_thromb[,c("patient_id","is_aki",thromb_list)]
    # thromb_list_tmp <- vector(mode="list",length=length(thromb_list))
    # for(i in 1:length(thromb_list)) {
    #     thromb_tmp1 <- thromb_tmp[,c("patient_id",thromb_list[i],"is_aki")]
    #     thromb_tmp2 <- thromb_tmp1 %>% dplyr::count(get(thromb_list[i]),is_aki)
    #     thromb_tmp3 <- thromb_tmp2 %>% dplyr::filter(is_aki == 1)
    #     if(min(thromb_tmp3$n) >= factor_cutoff) {
    #         thromb_list_tmp[i] <- thromb_list[i]
    #     }
    # }
    # thromb_list <- unlist(thromb_list_tmp[lengths(thromb_list_tmp) > 0L])
    # message("thromb_list for Thromb Analysis after filtering for CoxPH: ",paste(thromb_list,sep = ","))
    # 
    # aki_thromb_formula <- as.formula(paste("is_aki ~ severe",thromb_list),sep="+")
    # aki_thromb_logit <- glm(aki_thromb_formula,family = binomial,data=aki_index_thromb)
    # writeLines(capture.output(summary(aki_thromb_logit)),con=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_ThrombGLMSummary.txt"))
    # aki_thromb_logit_tidy <- aki_thromb_logit %>% broom::tidy(exponentiate=T,conf.int=T) %>% knitr::kable(align="l")
    # 
}

