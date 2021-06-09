
#' Runs the analytic workflow for the AKI project
#'
#' @keywords 4CE
#' @export

runAnalysis <- function(is_obfuscated=TRUE,factor_cutoff = 5,restrict_models = FALSE) {
    
    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)
    
    ## get the site identifier associated with the files stored in the /4ceData/Input directory that 
    ## is mounted to the container
    currSiteId = toupper(FourCePhase2.1Data::getSiteId())
    
    ## run the quality control
    FourCePhase2.1Data::runQC(currSiteId)
    
    ## DO NOT CHANGE ANYTHING ABOVE THIS LINE
    
    ## To Do: implement analytic workflow, saving results to a site-specific 
    ## file to be sent to the coordinating site later via submitAnalysis()
    obfuscation_value = as.numeric(FourCePhase2.1Data::getObfuscation(currSiteId))
    message(paste0(c("Obfuscation level set to ",obfuscation_value)))
    ## ========================================
    ## PART 1: Read in Data Tables
    ## ========================================
    message("\nPlease ensure that your working directory is set to /4ceData")
    message("Reading in Input/LocalPatientSummary.csv and Input/LocalPatientObservations.csv...")
    #demographics <- read.csv(paste0(FourCePhase2.1Data::getInputDataDirectoryName(),"/LocalPatientSummary.csv"),na.strings = '1/1/1900')
    demographics <- read.csv("Input/LocalPatientSummary.csv",na.strings = '1/1/1900')
    #observations <- read.csv(paste0(FourCePhase2.1Data::getInputDataDirectoryName(),"/LocalPatientObservations.csv"),na.strings = '-999')
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
    # thromb_diag[is.na(thromb_diag)] <- 0
    
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
    rrt_old <- rrt$patient_id[rrt$days_since_admission < 0]
    
    # Generate list of patients who were only initiated on RRT for the first time ever
    # during admission for COVID-19
    rrt_new <- rrt[!(rrt$patient_id %in% rrt_old),]
    
    # For debugging purposes
    message(paste("Number of patients on RRT in total: ",nrow(rrt),"\nRRT previously: ",nrow(rrt_old),"\nRRT during admission: ",nrow(rrt_new)))
    
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
    message("Removing patients who do not have any serum creatinine values during admission...")
    pts_valid_cr <- labs_cr_aki %>% dplyr::filter(days_since_admission >= 0)
    pts_valid_cr <- unique(pts_valid_cr$patient_id)
    demog_no_cr <- demographics_filt[!(demographics_filt$patient_id %in% pts_valid_cr),]
    demographics_filt <- demographics_filt[demographics_filt$patient_id %in% pts_valid_cr,]
    labs_cr_aki <- labs_cr_aki[labs_cr_aki$patient_id %in% pts_valid_cr,]
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
    
    # Filter for KDIGO grades > 0
    message("Now generating tables of all AKI events")
    labs_cr_aki_tmp5 <- labs_cr_aki_tmp4[labs_cr_aki_tmp4$aki_kdigo_final > 0,]
    labs_cr_aki_tmp5[is.na(labs_cr_aki_tmp5)] <- 0
    # Filter for maxima of delta_cr (which should give us the peaks)
    labs_cr_aki_tmp5 <- labs_cr_aki_tmp5[labs_cr_aki_tmp5$delta_is_max > 0,]
    
    # Filter and reorder columns to generate our final table of all AKI events
    labs_aki_summ <- labs_cr_aki_tmp5 %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
    
    labs_aki_summ <- labs_aki_summ %>% dplyr::distinct(patient_id,days_since_admission,.keep_all=TRUE)
    labs_aki_summ <- labs_aki_summ %>% dplyr::filter(days_since_admission >= 0)
    
    # Final headers for labs_aki_summ:
    # patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d
    # days_since_admission - time at which peak Cr is achieved
    # day_min - time at which Cr begins to rise
    
    # Generate a separate table (for reference) of all creatinine peaks not fulfilling KDIGO AKI criteria
    message("Now generating a separate table for serum Cr peaks which do not reach AKI definitions")
    labs_cr_nonaki <- labs_cr_aki_tmp4[!(labs_cr_aki_tmp4$patient_id %in% labs_aki_summ$patient_id),]
    labs_cr_nonaki[is.na(labs_cr_nonaki)] <- 0
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$delta_is_max > 0,]
    labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::select(patient_id,siteid,days_since_admission,value,day_min,day_min_retro,min_cr_90d,min_cr_48h,min_cr_retro_7day,min_cr_48h_retro,min_cr_7d_final,cr_7d,cr_90d,delta_cr,aki_kdigo,aki_kdigo_retro,aki_kdigo_final,akd_7d,akd_90d)
    labs_cr_nonaki <- labs_cr_nonaki %>% dplyr::filter(days_since_admission >= 0)
    # Generate the highest Cr peak for non-AKI peaks detected
    labs_nonaki_summ <- labs_cr_nonaki %>% dplyr::group_by(patient_id) %>% dplyr::slice(which.max(delta_cr))
    
    # We also want to generate tables to determine (1) if AKIs occurred before/after severe disease onset 
    # (2) how long before/after disease severity
    # These tables will help in segregating the populations for analysis later
    message("Generating intermediate tables to determine if AKIs occured before or after severe COVID-19 onset")
    severe_time <- demographics_filt %>% dplyr::select(patient_id,severe,time_to_severe)
    labs_aki_severe <- merge(labs_aki_summ,severe_time,by="patient_id",all.x=TRUE)
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_to_aki = ifelse(!is.na(time_to_severe), time_to_severe - day_min,NA))
    labs_aki_severe <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe_before_aki = ifelse(severe_to_aki < 0,1,0)) %>% dplyr::ungroup()
    labs_aki_severe <- labs_aki_severe %>% dplyr::distinct()
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
    
    # Generate simplified table for determining who were started on COAGA near admission
    message("Generating table for COAGA use during the admission...")
    coaga_present <- tryCatch({
        med_coaga_new <- med_new %>% dplyr::select(patient_id,COAGA) %>% dplyr::filter(COAGA == min(COAGA)) %>% dplyr::distinct()
        med_coaga_new$COAGA[med_coaga_new$COAGA < -15] <- 0
        med_coaga_new$COAGA[med_coaga_new$COAGA >= -15] <- 1
        TRUE
    },error= function(c) FALSE)
    
    # Generate simplified table for determining who were started on COAGB near admission
    message("Generating table for COAGB use during the admission...")
    coagb_present <- tryCatch({
        med_coagb_new <- med_new %>% dplyr::select(patient_id,COAGB) %>% dplyr::filter(COAGB == min(COAGB)) %>% dplyr::distinct()
        med_coagb_new$COAGB[med_coagb_new$COAGB < -15] <- 0
        med_coagb_new$COAGB[med_coagb_new$COAGB >= -15] <- 1
        TRUE
    },error=function(c) FALSE)
    
    # Generate simplified table for determining who were using ACEI near admission
    message("Generating table for ACEI use during the admission...")
    acei_present <- tryCatch({
        med_acei_new <- med_new %>% dplyr::select(patient_id,ACEI) %>% dplyr::filter(ACEI == min(ACEI)) %>% dplyr::distinct()
        med_acei_new$ACEI[med_acei_new$ACEI < -15] <- 0
        med_acei_new$ACEI[med_acei_new$ACEI >= -15] <- 1
        TRUE
    },error=function(c) FALSE)
    arb_present <- tryCatch({
        med_arb_new <- med_new %>% dplyr::select(patient_id,ARB) %>% dplyr::filter(ARB == min(ARB)) %>% dplyr::distinct()
        med_arb_new$ARB[med_arb_new$ARB < -15] <- 0
        med_arb_new$ARB[med_arb_new$ARB >= -15] <- 1
        TRUE
    },error=function(c) FALSE)
    
    if(isTRUE(acei_present) & isTRUE(arb_present)){
        med_raas_new <- med_acei_new
        colnames(med_raas_new)[2] <- "raas_new"
        tmp <- med_arb_new
        colnames(tmp)[2] <- "raas_new"
        med_raas_new <- rbind(med_raas_new,tmp) %>% dplyr::distinct() %>% dplyr::group_by(patient_id) %>% dplyr::filter(raas_new == max(raas_new)) %>% dplyr::ungroup()
    }
    
    # Generate simplified table for determining who were started on novel antivirals
    covid19antiviral_present = ("COVIDVIRAL" %in% colnames(med_new))
    remdesivir_present = ("REMDESIVIR" %in% colnames(med_new))
    if(isTRUE(covid19antiviral_present) && isTRUE(remdesivir_present)){
        message("Generating table for experimental COVID-19 treatment... (both COVIDVIRAL and REMDESIVIR)")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,COVIDVIRAL,REMDESIVIR) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(COVIDVIRAL=ifelse(is.na(COVIDVIRAL),0,COVIDVIRAL),REMDESIVIR=ifelse(is.na(REMDESIVIR),0,REMDESIVIR)) %>% dplyr::ungroup()
        # Uncomment following line to select for patients who were started on novel antivirals less than 72h from admission
        #med_covid19_new <- med_covid19_new[med_covid19_new$COVIDVIRAL <= 3,]
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(COVIDVIRAL >= 0 | REMDESIVIR >= 0,1,0)) %>% dplyr::ungroup()
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx)
    } else if(isTRUE(covid19antiviral_present)){
        message("Generating table for experimental COVID-19 treatment... (COVIDVIRAL only)")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,COVIDVIRAL) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(COVIDVIRAL=ifelse(is.na(COVIDVIRAL),0,COVIDVIRAL)) %>% dplyr::ungroup()
        # Uncomment following line to select for patients who were started on novel antivirals less than 72h from admission
        #med_covid19_new <- med_covid19_new[med_covid19_new$COVIDVIRAL <= 3,]
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(COVIDVIRAL >= 0,1,0)) %>% dplyr::ungroup()
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx)
    } else if(isTRUE(remdesivir_present)){
        message("Generating table for experimental COVID-19 treatment... (REMDESIVIR only)")
        med_covid19_new <- med_new %>% dplyr::select(patient_id,REMDESIVIR) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(REMDESIVIR=ifelse(is.na(REMDESIVIR),0,REMDESIVIR)) %>% dplyr::ungroup()
        # Uncomment following line to select for patients who were started on novel antivirals less than 72h from admission
        #med_covid19_new <- med_covid19_new[med_covid19_new$COVIDVIRAL <= 3,]
        med_covid19_new <- med_covid19_new %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covid_rx = ifelse(REMDESIVIR >= 0,1,0)) %>% dplyr::ungroup()
        med_covid19_new_date <- med_covid19_new
        colnames(med_covid19_new_date)[2] <- "covid_rx_start"
        med_covid19_new <- med_covid19_new %>% dplyr::select(patient_id,covid_rx)
    }
    
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
    aki_only_index <- labs_aki_severe %>% dplyr::group_by(patient_id) %>% tidyr::fill(severe) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::filter(delta_cr == max(delta_cr)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE)
    # patient_id	site_id	days_since_admission	value	day_min	day_min_retro	min_cr_90d	min_cr_48h	min_cr_retro_7day	min_cr_48h_retro	min_cr_7d_final	cr_7d	cr_90d	delta_cr	aki_kdigo	aki_kdigo_retro	aki_kdigo_final	akd_7d	akd_90d	severe  time_to_severe	severe_to_aki	severe_before_aki
    
    # Generate the patient list including (1) severity indices from this dplyr::filtered table (2) day of peak Cr
    # severe - 2 = never severe, 4 = severe, AKI
    #aki_only_index <- aki_only_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = dplyr::if_else(!is.na(time_to_severe),4,2)) %>% dplyr::ungroup()
    aki_only_index <- aki_only_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = 2 * severe + 2) %>% dplyr::ungroup() # if severe = 0 initially, will be recoded as 2; if severe = 1 initially, will be recoded as 4
    
    # create the change in baseline index table
    aki_only_index_baseline_shift <- aki_only_index %>% dplyr::select(patient_id,severe,cr_7d,cr_90d)
    
    aki_only_index <- aki_only_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
    colnames(aki_only_index)[2] <- "peak_cr_time"
    colnames(aki_only_index)[4] <- "aki_start"
    # Headers of aki_only_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki
    message("Now creating list of non-AKI patients...")
    no_aki_list <- demographics_filt %>% dplyr::select(patient_id,severe)
    no_aki_list <- no_aki_list[!(no_aki_list$patient_id %in% aki_only_index$patient_id),]
    no_aki_list <- no_aki_list %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(is.na(severe),1,2 * severe + 1)) %>% dplyr::ungroup() # if severe = 0, will be recoded as 1; if severe = 1, will be recoded as 3
    
    labs_nonaki_summ <- labs_nonaki_summ[labs_nonaki_summ$patient_id %in% no_aki_list$patient_id,]
    labs_nonaki_severe <- labs_nonaki_severe[labs_nonaki_severe$patient_id %in% no_aki_list$patient_id,]
    labs_cr_nonaki <- labs_cr_nonaki[labs_cr_nonaki$patient_id %in% no_aki_list$patient_id,]
    
    message("Creating table for non-AKI patients...")
    # Create a non-AKI equivalent for aki_only_index - except that this takes the largest delta_cr (and the earliest occurence of such a delta_cr)
    no_aki_index <- labs_nonaki_severe %>% dplyr::group_by(patient_id) %>% dplyr::arrange(severe,.by_group = TRUE) %>% tidyr::fill(severe) %>% dplyr::filter(delta_cr == max(delta_cr, na.rm = TRUE)) %>% dplyr::filter(days_since_admission == min(days_since_admission, na.rm = TRUE)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE) %>% dplyr::ungroup()
    # no_aki_index <- labs_nonaki_severe %>% dplyr::group_by(patient_id) %>% dplyr::arrange(dplyr::desc(delta_cr),.by_group=TRUE) %>% tidyr::fill(severe,delta_cr) %>% dplyr::ungroup()
    # try({no_aki_index$delta_cr[is.na(no_aki_index$delta_cr)] <- 0})
    # no_aki_index <- no_aki_index %>% dplyr::group_by(patient_id) %>% dplyr::filter(delta_cr == max(delta_cr)) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::distinct(days_since_admission,.keep_all = TRUE) %>% dplyr::ungroup()
    no_aki_index <- no_aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse(is.na(severe),1,2 * severe + 1))
    # create the change in baseline index table
    no_aki_index_baseline_shift <- no_aki_index %>% dplyr::select(patient_id,severe,cr_7d,cr_90d)
    
    no_aki_index <- no_aki_index %>% dplyr::select(patient_id,days_since_admission,severe,day_min,severe_to_aki)
    colnames(no_aki_index)[2] <- "peak_cr_time"
    colnames(no_aki_index)[4] <- "aki_start"
    
    #no_aki_index$severe_to_aki <- -999
    message("Creating final table of peak Cr for all patients...")
    aki_index <- dplyr::bind_rows(aki_only_index,no_aki_index)
    
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        message("Adding on temporal information for experimental COVID-19 antivirals...")
        aki_index <- merge(aki_index,med_covid19_new,by="patient_id",all.x=TRUE)
        aki_index$covid_rx[is.na(aki_index$covid_rx)] <- 0
        aki_index <- aki_index %>% dplyr::group_by(patient_id) %>% dplyr::mutate(covidrx_grp = dplyr::if_else(severe <= 2, dplyr::if_else(covid_rx == 0,1,2),dplyr::if_else(covid_rx == 0,3,4))) %>% dplyr::ungroup()
        message(paste0(c("Column names for aki_index",colnames(aki_index))))
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki,covidrx_grp)
    } else {
        aki_index <- aki_index %>% dplyr::select(patient_id,peak_cr_time,severe,aki_start,severe_to_aki)
    }
    # Headers of aki_index: patient_id  peak_cr_time  severe  aki_start  severe_to_aki  covidrx_grp
    aki_index <- aki_index %>% dplyr::arrange(patient_id,peak_cr_time,desc(severe))%>% dplyr::distinct(patient_id,peak_cr_time,.keep_all = TRUE)
    message("Final table aki_index created.")
    
    # Uncomment the following line to remove patients who were previously on RRT prior to admission
    # aki_index <- aki_index[!(aki_index$patient_id %in% rrt_new),]
    message("Now generating table of normalised serum Cr...")
    # Create a common labs_cr_all table containing the serum Cr values, the severity groupings and anti-viral groupings
    labs_cr_aki_tmp <- labs_cr_aki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_90d,min_cr_retro_7day)
    labs_cr_nonaki_tmp <- labs_cr_nonaki %>% dplyr::select(patient_id,days_since_admission,value,min_cr_90d,min_cr_retro_7day)
    labs_cr_all <- dplyr::bind_rows(labs_cr_aki_tmp,labs_cr_nonaki_tmp)
    labs_cr_all <- merge(labs_cr_all,aki_index,by="patient_id",all.x=TRUE)
    
    # Now, generate a table containing lab values with timepoints calculated from time of peak cr
    peak_trend <- labs_cr_all
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
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
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day,time_from_peak,baseline_cr,first_baseline_cr)
    } else {
        peak_trend <- peak_trend %>% dplyr::select(patient_id,severe,days_since_admission,peak_cr_time,value,min_cr_90d,min_cr_retro_7day,time_from_peak,baseline_cr,first_baseline_cr)
    }
    peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::arrange(severe,first_baseline_cr) %>% tidyr::fill(severe,first_baseline_cr) %>% dplyr::mutate(ratio = value/first_baseline_cr) %>% dplyr::ungroup() %>% dplyr::distinct()
    #peak_trend <- peak_trend %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    message("Final table of peak Cr for all patients - peak_trend - created.")
    # peak_trend will now be a common table to plot from the dplyr::selected AKI peak
    
    baseline_shift <- dplyr::bind_rows(aki_only_index_baseline_shift,no_aki_index_baseline_shift) %>% dplyr::distinct()
    baseline_shift <- merge(baseline_shift,first_baseline,by="patient_id",all.x=T)
    baseline_shift <- baseline_shift %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio_7d = cr_7d/first_baseline_cr,ratio_90d = cr_90d/first_baseline_cr) %>% dplyr::select(patient_id,severe,ratio_7d,ratio_90d)
    baseline_shift <- baseline_shift %>% dplyr::group_by(severe) %>% dplyr::summarise(n_all=dplyr::n(),n_shift_7d = sum(ratio_7d >= 1.25, na.rm=T),n_shift_90d = sum(ratio_90d >= 1.25, na.rm=T)) %>% dplyr::ungroup()
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        baseline_shift$n_all[baseline_shift$n_all < obfuscation_value] <- NA
        baseline_shift$n_shift_7d[baseline_shift$n_shift_7d < obfuscation_value] <- NA
        baseline_shift$n_shift_90d[baseline_shift$n_shift_90d < obfuscation_value] <- NA
    }
    write.csv(baseline_shift,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_BaselineShift_Counts.csv")),row.names=FALSE)
    
    # =======================================================================================
    # Figure 1(a): Cr trends from start of AKI / after peak Cr for severe vs non-severe groups
    # =======================================================================================
    
    # First create a plot of the creatinine trends of AKI vs non-AKI patients from the day of the first AKI peak (or the highest Cr peak for non-AKI patients)
    peak_aki_vs_non_aki <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio) %>% dplyr::arrange(patient_id,severe,time_from_peak,ratio)
    peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% dplyr::group_by(patient_id,time_from_peak) %>% tidyr::fill(ratio,severe) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,time_from_peak,.keep_all=TRUE)
    colnames(peak_aki_vs_non_aki) <- c("patient_id","aki","time_from_peak","ratio")
    peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = ifelse((aki == 2 | aki == 4),1,0))
    peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki %>% dplyr::group_by(aki,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    aki_label <- data.table::data.table(c(0,1),c("Non-AKI","AKI"))
    colnames(aki_label) <- c("aki","aki_label")
    peak_aki_vs_non_aki_summ <- merge(peak_aki_vs_non_aki_summ,aki_label,by="aki",all.x=TRUE)
    if(isTRUE(is_obfuscated)) {
        # peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki_summ %>% dplyr::filter(n >= obfuscation_value) %>% dplyr::arrange(aki,time_from_peak)
        message("Obfuscating the AKI vs non-AKI graphs...")
        peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki_summ[peak_aki_vs_non_aki_summ$n >= obfuscation_value,]
    }
    write.csv(peak_aki_vs_non_aki_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrFromPeak_AKI_vs_NonAKI.csv")),row.names=FALSE)
    peak_aki_vs_non_aki_timeplot <- ggplot2::ggplot(peak_aki_vs_non_aki_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=aki_label))+ggplot2::geom_line(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio, color=factor(aki_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "AKI Group") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("AKI"="#bc3c29","Non-AKI"="#0072b5")) + ggplot2::theme_minimal()
    print(peak_aki_vs_non_aki_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_AKI_vs_NonAKI.png")),plot=peak_aki_vs_non_aki_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients should have been generated.")
    
    # Create KDIGO Stage table for demographics table
    kdigo_grade <- peak_aki_vs_non_aki %>% dplyr::filter(time_from_peak == 0) %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki_kdigo_stage = ifelse(aki == 0,0,ifelse(ratio < 2,1,ifelse(ratio<3,2,3)))) %>% dplyr::ungroup()
    kdigo_grade <- kdigo_grade %>% dplyr::select(patient_id,aki_kdigo_stage)
    
    # Now, derive our first table peak_trend_severe to compare across the different severity groups
    peak_trend_severe <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio) %>% dplyr::arrange(patient_id,severe,time_from_peak,ratio)
    peak_trend_severe <- peak_trend_severe %>% dplyr::group_by(patient_id,time_from_peak) %>% tidyr::fill(ratio,severe) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,time_from_peak,.keep_all=TRUE)
    # Headers: patient_id  severe  time_from_peak  ratio
    # peak_trend_severe <- peak_trend_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    
    # Calculate mean and SD each for severe and non-severe groups
    peak_cr_summ <- peak_trend_severe %>% dplyr::group_by(severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    severe_label <- data.table::data.table(c(1,2,3,4),c("Non-severe, no AKI","Non-severe, AKI","Severe, no AKI","Severe, AKI"))
    colnames(severe_label) <- c("severe","severe_label")
    peak_cr_summ <- merge(peak_cr_summ,severe_label,by="severe",all.x=TRUE)
    if(isTRUE(is_obfuscated)) {
        # peak_cr_summ <- peak_cr_summ %>% dplyr::filter(n >= obfuscation_value)
        message("Obfuscating the AKI with severity graphs...")
        peak_cr_summ <- peak_cr_summ[peak_cr_summ$n >= obfuscation_value,]
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
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(peak_cr_time == min(peak_cr_time)) %>% dplyr::distinct() %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    #adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    adm_to_aki_summ <- adm_to_aki_cr %>% dplyr::group_by(severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    adm_to_aki_summ <- merge(adm_to_aki_summ,severe_label,by="severe",all.x=TRUE)
    if(isTRUE(is_obfuscated)) {
        # adm_to_aki_summ <- adm_to_aki_summ %>% dplyr::filter(n >= obfuscation_value)
        message("Obfuscating the admission to AKI graphs...")
        adm_to_aki_summ <- adm_to_aki_summ[adm_to_aki_summ$n >= obfuscation_value,]
    }
    write.csv(adm_to_aki_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.csv")),row.names=FALSE)
    adm_to_aki_timeplot <- ggplot2::ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 30 & adm_to_aki_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(0,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
    print(adm_to_aki_timeplot)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.png")),plot=adm_to_aki_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from first day of admission, should have been generated.")
    
    # Plot from start of AKI to 30 days later 
    
    aki_30d_cr <- labs_cr_all
    # Uncomment the following line to restrict analysis to AKI patients only
    #aki_30d_cr <- aki_30d_cr[aki_30d_cr$severe %in% c(2,4,5),]
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_start = days_since_admission - aki_start) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr[order(aki_30d_cr$patient_id,aki_30d_cr$days_since_admission),] %>% dplyr::distinct()
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
    aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
    #aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
    aki_30d_cr_summ <- aki_30d_cr %>% dplyr::group_by(severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    aki_30d_cr_summ <- merge(aki_30d_cr_summ,severe_label,by="severe",all.x=TRUE)
    if(isTRUE(is_obfuscated)) {
        message("Obfuscating the start of AKI graphs...")
        # aki_30d_cr_summ <- aki_30d_cr_summ %>% dplyr::filter(n >= obfuscation_value)
        aki_30d_cr_summ <- aki_30d_cr_summ[aki_30d_cr_summ$n >= obfuscation_value,]
    }
    write.csv(aki_30d_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.csv")),row.names=FALSE)
    aki_30d_cr_timeplot <- ggplot2::ggplot(aki_30d_cr_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
    print(aki_30d_cr_timeplot)
    ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.png")),plot=aki_30d_cr_timeplot,width=12,height=9,units="cm")
    message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from start of AKI/creatinine increase, should have been generated.")
    
    ## ===============================================================================
    ## Initializing analysis for cirrhotic patients + plotting serum creatinine graphs
    ## ===============================================================================
    
    cirrhosis_present <- ("cld" %in% comorbid_list)
    inr_loinc <- c("6301-6","34714-6","38875-1","46418-0","52129-4","61189-7","72281-9","92891-1")
    meld_analysis_valid <- FALSE
    if(isTRUE(cirrhosis_present)) {
        cirrhosis_list <- comorbid %>% dplyr::select(patient_id,cld) %>% dplyr::filter(cld == 1)
        labs_cirrhosis <- observations[observations$patient_id %in% cirrhosis_list$patient_id,] %>% dplyr::filter(concept_type == "LAB-LOINC")
        labs_list <- unique(labs_cirrhosis$concept_code)
        inr_present <- FALSE
        if(length(intersect(inr_loinc,labs_list)) > 0) {
            inr_present <- TRUE
        }
        if(isTRUE(inr_present)) {
            message("=====================================")
            message("Found cirrhotic patients and INR values in the Observations table. Will proceed with sub-group analysis for hepatorenal syndrome.")
            meld_analysis_valid <- TRUE
            # platelet_loinc <- c("13056-7","26515-7","49497-1","74464-9","777-3","778-1")
            
            sodium_loinc <- c("2947-0","32717-1","39792-7","41657-8","39791-9","2951-2","77139-4")
            # platelet_present <- FALSE
            # if(length(intersect(platelet_loinc,labs_list)) > 0) {
            #     platelet_present <- TRUE
            # }
            sodium_present <- FALSE
            if(length(intersect(sodium_loinc,labs_list)) > 0) {
                sodium_present <- TRUE
                message("Sodium values present for MELD score correction.")
            }
            
            message("Extracting first discharge dates...")
            admissions <- read.csv("Input/LocalPatientClinicalCourse.csv")
            admissions <- admissions %>% dplyr::mutate(patient_id=paste(currSiteId,patient_num,sep="_"))
            first_discharge <- admissions %>% dplyr::group_by(patient_id) %>% dplyr::filter(in_hospital == 0) %>% dplyr::filter(days_since_admission >= 0) %>% dplyr::filter(days_since_admission == min(days_since_admission)) %>% dplyr::ungroup()
            first_discharge <- first_discharge %>% dplyr::select(patient_id,days_since_admission)
            colnames(first_discharge)[2] <- "first_discharge_day"
            
            message("Restricting labs to first admission only")
            labs_cirrhosis_firstdischarge <- merge(labs_cirrhosis,first_discharge,by="patient_id",all.x=TRUE) %>% dplyr::group_by(patient_id) %>% dplyr::filter(days_since_admission <= first_discharge_day & days_since_admission >= 0) %>% dplyr::ungroup()
            message("Extracting and binning INR")
            try({
                labs_inr <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code %in% c("6301-6","34714-6","38875-1","46418-0","52129-4","61189-7","72281-9","92891-1"),]
                labs_inr <- labs_inr[,-c(4,5)]
                labs_inr <- labs_inr %>% dplyr::filter(days_since_admission >= 0)
                labs_inr <- labs_inr %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_inr = mean(na.omit(value)),min_inr = min(na.omit(value)),max_inr = max(na.omit(value)),first_inr = dplyr::first(na.omit(value)))
            })
            message("Extracting and binning bilirubin")
            try({
                labs_bil <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code == '1975-2',]
                labs_bil <- labs_bil[,-c(4,5)]
                labs_bil <- labs_bil %>% dplyr::filter(days_since_admission >= 0)
                labs_bil <- labs_bil %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_bil = mean(na.omit(value)),min_bil = min(na.omit(value)),max_bil = max(na.omit(value)),first_bil = dplyr::first(na.omit(value)))
            })
            message("Extracting and binning Cr")
            try({
                labs_cr <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code == '2160-0',]
                labs_cr <- labs_cr[,-c(4,5)]
                labs_cr <- labs_cr %>% dplyr::filter(days_since_admission >= 0)
                labs_cr <- labs_cr %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_cr = mean(na.omit(value)),min_cr = min(na.omit(value)),max_cr = max(na.omit(value)),first_cr = dplyr::first(na.omit(value)))
            })
            message("Extracting and binning AST")
            try({
                labs_ast <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code == '1920-8',]
                labs_ast <- labs_ast[,-c(4,5)]
                labs_ast <- labs_ast %>% dplyr::filter(days_since_admission >= 0)
                labs_ast <- labs_ast %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_ast = mean(na.omit(value)),min_ast = min(na.omit(value)),max_ast = max(na.omit(value)),first_ast = dplyr::first(na.omit(value)))
            })
            message("Extracting and binning ALT")
            try({
                labs_alt <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code == '1742-6',]
                labs_alt <- labs_alt[,-c(4,5)]
                labs_alt <- labs_alt %>% dplyr::filter(days_since_admission >= 0)
                labs_alt <- labs_alt %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_alt = mean(na.omit(value)),min_alt = min(na.omit(value)),max_alt = max(na.omit(value)),first_alt = dplyr::first(na.omit(value)))
            })
            message("Extracting and binning albumin")
            try({
                labs_alb <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code == '1751-7',]
                labs_alb <- labs_alb[,-c(4,5)]
                labs_alb <- labs_alb %>% dplyr::filter(days_since_admission >= 0)
                labs_alb <- labs_alb %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_alb = mean(na.omit(value)),min_alb = min(na.omit(value)),max_alb = max(na.omit(value)),first_alb = dplyr::first(na.omit(value)))
            })
            message("Merging all tables with binned data")
            labs_meld <- merge(labs_inr,labs_bil,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            labs_meld <- merge(labs_meld,labs_alb,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            labs_meld <- merge(labs_meld,labs_ast,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            labs_meld <- merge(labs_meld,labs_alt,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            labs_meld <- merge(labs_meld,labs_cr,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            
            if(isTRUE(sodium_present)) {
                message("Adding in sodium data")
                labs_na <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code %in% sodium_loinc,]
                labs_na <- labs_na[,-c(4,5)]
                labs_na <- labs_na %>% dplyr::filter(days_since_admission >= 0)
                labs_na <- labs_na %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_na = mean(na.omit(value)),min_na = min(na.omit(value)),max_na = max(na.omit(value)),first_na = dplyr::first(na.omit(value)))
                labs_meld <- merge(labs_meld,labs_na,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            }
            
            # if(isTRUE(platelet_present)) {
            #     message("Adding in platelet data")
            #     labs_plt <- labs_cirrhosis_firstdischarge[labs_cirrhosis_firstdischarge$concept_code %in% platelet_loinc,]
            #     labs_plt <- labs_plt[,-c(4,5)]
            #     labs_plt <- labs_plt %>% dplyr::filter(days_since_admission >= 0)
            #     labs_plt <- labs_plt %>% dplyr::mutate(day_bin = ggplot2::cut_width(days_since_admission,width=3,boundary=0)) %>% dplyr::group_by(patient_id,day_bin) %>% dplyr::summarise(mean_plt = mean(na.omit(value)),min_plt = min(na.omit(value)),max_plt = max(na.omit(value)),first_plt = dplyr::first(na.omit(value)))
            #     labs_meld <- merge(labs_meld,labs_plt,by=c("patient_id","day_bin"),all=T) %>% dplyr::distinct()
            # }
            
            message("Imputing empty fields prior to MELD score calculation")
            labs_meld <- labs_meld %>% dplyr::group_by(patient_id,day_bin) %>% tidyr::fill(dplyr::everything()) %>% dplyr::distinct()
            message("Calculating MELD score...")
            if(isTRUE(sodium_present)) {
                labs_meld <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld = FourCePhase2.1AKI:::meld_score(bil = max_bil,inr = max_inr,sCr = max_cr,Na = min_na)) %>% dplyr::ungroup()
            } else {
                labs_meld <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(meld = FourCePhase2.1AKI:::meld_score(bil = max_bil,inr = max_inr,sCr = max_cr)) %>% dplyr::ungroup()
            }
            message("Extracting admission MELD score...")
            labs_meld_admission <- labs_meld %>% dplyr::group_by(patient_id) %>% dplyr::filter(day_bin == "[0,3]") %>% dplyr::mutate(meld_admit_severe = dplyr::if_else(meld >= 20,1,0)) %>% dplyr::ungroup()
            labs_meld_admission$meld_admit_severe[is.na(labs_meld_admission$meld_admit_severe)] <- 0
            meld_severe_list <- labs_meld_admission %>% dplyr::select(patient_id,meld,meld_admit_severe) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
            
            # Final headers
            # labs_meld:
            # patient_id, day_bin, mean/min/max/first of labs, meld
            #
            # labs_meld_admission:
            # patient_id, day_bin, mean/min/max/first of labs, meld (integer score), meld_admit_severe (0/1)
            # Possible that some patient_ids may not be inside this if there are no labs in days 0-3 (i.e. no day_bin == "[0,3]")
            peak_trend_meld <- peak_trend[peak_trend$patient_id %in% cirrhosis_list$patient_id,]
            peak_trend_meld <- merge(peak_trend_meld,meld_severe_list,by="patient_id",all.x=TRUE)
            peak_trend_meld$meld_admit_severe[peak_trend_meld$meld_admit_severe] <- 0
            peak_cr_meld_summ <- peak_trend_meld %>% dplyr::group_by(meld_admit_severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
            meld_label <- data.table::data.table(c(1,2,3,4),c("MELD < 20, no AKI","MELD < 20, AKI","MELD ≥ 20, no AKI","MELD ≥ 20, AKI"))
            colnames(meld_label) <- c("meld_admit_severe","meld_severe_label")
            peak_cr_meld_summ <- merge(peak_cr_meld_summ,meld_label,by="meld_admit_severe",all.x=TRUE)
            if(isTRUE(is_obfuscated)) {
                message("Obfuscating the MELD AKI graphs...")
                peak_cr_meld_summ <- peak_cr_meld_summ[peak_cr_meld_summ$n >= obfuscation_value,]
            }
            peak_cr_meld_timeplot <- ggplot2::ggplot(peak_cr_meld_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=meld_severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(meld_severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("MELD < 20, AKI"="#bc3c29","MELD < 20, no AKI"="#0072b5","MELD ≥ 20, AKI" = "#e18727","MELD ≥ 20, no AKI"="#20854e")) + ggplot2::theme_minimal()
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromPeak_MELD_AKI.png")),plot=peak_cr_meld_timeplot,width=12,height=9,units="cm")
            
            peak_trend_severe_cld <- peak_trend_severe[peak_trend_severe$patient_id %in% cirrhosis_list$patient_id,]
            peak_cr_cld_summ <- peak_trend_severe_cld %>% dplyr::group_by(severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
            peak_cr_cld_summ <- merge(peak_cr_cld_summ,severe_label,by="severe",all.x=TRUE)
            if(isTRUE(is_obfuscated)) {
                # peak_cr_summ <- peak_cr_summ %>% dplyr::filter(n >= obfuscation_value)
                message("Obfuscating the AKI with severity graphs...")
                peak_cr_cld_summ <- peak_cr_cld_summ[peak_cr_cld_summ$n >= obfuscation_value,]
            }
            peak_cr_cld_timeplot <- ggplot2::ggplot(peak_cr_cld_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromPeak_AKI_CLD_only.png")),plot=peak_cr_timeplot,width=12,height=9,units="cm")
            message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of cirrhotic patients (with severity) should have been generated.")
            
            # Plot from start of admission to 30 days post-peak AKI (if no AKI, then from peak Cr)
            adm_meld_cr <- labs_cr_all[labs_cr_all$patient_id %in% cirrhosis_list$patient_id,]
            adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(peak_cr_time == min(peak_cr_time)) %>% dplyr::distinct() %>% dplyr::ungroup()
            adm_meld_cr <- adm_meld_cr[order(adm_meld_cr$patient_id,adm_meld_cr$days_since_admission),]
            adm_meld_cr <- merge(adm_meld_cr,meld_severe_list,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
            adm_meld_cr$meld_admit_severe[is.na(adm_meld_cr$meld_admit_severe)] <- 0
            adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
            adm_meld_cr <- adm_meld_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
            adm_meld_summ <- adm_meld_cr %>% dplyr::group_by(meld_admit_severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
            adm_meld_summ <- merge(adm_meld_summ,meld_label,by="meld_admit_severe",all.x=TRUE)
            if(isTRUE(is_obfuscated)) {
                # adm_to_aki_summ <- adm_to_aki_summ %>% dplyr::filter(n >= obfuscation_value)
                message("Obfuscating the admission to AKI graphs...")
                adm_meld_summ <- adm_meld_summ[adm_meld_summ$n >= obfuscation_value,]
            }
            adm_meld_timeplot <- ggplot2::ggplot(adm_meld_summ[which(adm_meld_summ$days_since_admission <= 30 & adm_meld_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=meld_severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(meld_severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(0,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("MELD < 20, AKI"="#bc3c29","MELD < 20, no AKI"="#0072b5","MELD ≥ 20, AKI" = "#e18727","MELD ≥ 20, no AKI"="#20854e")) + ggplot2::theme_minimal()
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdm_MELD_AKI.png")),plot=adm_meld_timeplot,width=12,height=9,units="cm")
            
            adm_to_aki_cld_summ <- adm_to_aki_cr[adm_to_aki_cr$patient_id %in% cirrhosis_list$patient_id,]
            adm_to_aki_cld_summ <- adm_to_aki_cld_summ %>% dplyr::group_by(severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
            adm_to_aki_cld_summ <- merge(adm_to_aki_cld_summ,severe_label,by="severe",all.x=TRUE)
            if(isTRUE(is_obfuscated)) {
                # adm_to_aki_cld_summ <- adm_to_aki_cld_summ %>% dplyr::filter(n >= obfuscation_value)
                message("Obfuscating the admission to AKI graphs...")
                adm_to_aki_cld_summ <- adm_to_aki_cld_summ[adm_to_aki_cld_summ$n >= obfuscation_value,]
            }
            adm_to_aki_cld_timeplot <- ggplot2::ggplot(adm_to_aki_cld_summ[which(adm_to_aki_cld_summ$days_since_admission <= 30 & adm_to_aki_cld_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(0,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdm_AKI_CLD_only.png")),plot=adm_to_aki_cld_timeplot,width=12,height=9,units="cm")
            
            message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of cirrhotic patients, plotted from first day of admission, should have been generated.")
            
            # Plot from start of AKI to 30 days later 
            
            aki_start_meld <- labs_cr_all[labs_cr_all$patient_id %in% cirrhosis_list$patient_id,]
            # Uncomment the following line to restrict analysis to AKI patients only
            #aki_start_meld <- aki_start_meld[aki_start_meld$severe %in% c(2,4,5),]
            aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_start = days_since_admission - aki_start) %>% dplyr::ungroup()
            aki_start_meld <- aki_start_meld[order(aki_start_meld$patient_id,aki_start_meld$days_since_admission),] %>% dplyr::distinct()
            aki_start_meld <- merge(aki_start_meld,meld_severe_list,by="patient_id",all.x=TRUE) %>% dplyr::distinct()
            aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_90d,min_cr_retro_7day)) %>% dplyr::ungroup()
            aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
            #aki_start_meld <- aki_start_meld %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
            aki_start_meld_summ <- aki_start_meld %>% dplyr::group_by(meld_admit_severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
            aki_start_meld_summ <- merge(aki_start_meld_summ,meld_label,by="meld_admit_severe",all.x=TRUE)
            if(isTRUE(is_obfuscated)) {
                message("Obfuscating the start of AKI graphs...")
                # aki_start_meld_summ <- aki_start_meld_summ %>% dplyr::filter(n >= obfuscation_value)
                aki_start_meld_summ <- aki_start_meld_summ[aki_start_meld_summ$n >= obfuscation_value,]
            }
            aki_start_meld_timeplot <- ggplot2::ggplot(aki_start_meld_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=meld_severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(meld_severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(meld_severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("MELD < 20, AKI"="#bc3c29","MELD < 20, no AKI"="#0072b5","MELD ≥ 20, AKI" = "#e18727","MELD ≥ 20, no AKI"="#20854e")) + ggplot2::theme_minimal()
            ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_MELD_AKI.png")),plot=aki_start_meld_timeplot,width=12,height=9,units="cm")
            
            aki_start_cld_summ <- aki_30d_cr[aki_30d_cr$patient_id %in% cirrhosis_list$patient_id,]
            aki_start_cld_summ <- aki_start_cld_summ %>% dplyr::group_by(severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
            aki_start_cld_summ <- merge(aki_start_cld_summ,severe_label,by="severe",all.x=TRUE)
            if(isTRUE(is_obfuscated)) {
                message("Obfuscating the start of AKI graphs...")
                # aki_start_cld_summ <- aki_start_cld_summ %>% dplyr::filter(n >= obfuscation_value)
                aki_start_cld_summ <- aki_start_cld_summ[aki_start_cld_summ$n >= obfuscation_value,]
            }
            aki_start_cld_timeplot <- ggplot2::ggplot(aki_start_cld_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
            ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_AKI_CLD_only.png")),plot=aki_start_cld_timeplot,width=12,height=9,units="cm")
            message("At this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from start of AKI/creatinine increase, should have been generated.")
        }
    }
    
    if(isTRUE(meld_analysis_valid)) {
        save(peak_cr_cld_summ,peak_cr_meld_summ,peak_cr_cld_timeplot,peak_cr_meld_timeplot,adm_to_aki_cld_summ,adm_meld_summ,adm_to_aki_cld_timeplot,adm_meld_timeplot,aki_start_cld_summ,aki_start_meld_summ,aki_start_cld_timeplot,aki_start_meld_timeplot,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_MELD_CLD_graphs.rda")))
    }
    
    # =====================
    # Demographics Table
    # =====================
    message("Now generating the demographics table...")
    demog_summ <- demographics_filt %>% dplyr::select(patient_id,sex,age_group,race,severe,deceased,time_to_severe,time_to_death) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
    demog_summ <- merge(demog_summ,comorbid,by="patient_id",all.x=TRUE)
    demog_summ <- merge(demog_summ,kdigo_grade,by="patient_id",all.x=TRUE)
    
    if(isTRUE(meld_analysis_valid)) {
        demog_summ <- merge(demog_summ,labs_meld_admission[,-2],by="patient_id",all.x=TRUE)
        demog_summ$meld_admit_severe[is.na(demog_summ$meld_admit_severe)] <- 0
        demog_summ$meld_admit_severe <- factor(demog_summ$meld_admit_severe,levels=c(0,1),labels = c("MELD < 20","MELD ≥ 20"))
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
    if(isTRUE(acei_present) & isTRUE(arb_present)) {
        demog_summ <- merge(demog_summ,med_raas_new,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
        demog_summ[is.na(demog_summ)] <- 0
        demog_summ$raas_new <- factor(demog_summ$raas_new,levels=c(0,1),labels=c("No ACE-i/ARB","ACE-i and/or ARB"))
    } else if(isTRUE(acei_present)) {
        demog_summ <- merge(demog_summ,med_acei_new,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
        demog_summ[is.na(demog_summ)] <- 0
        demog_summ$ACEI <- factor(demog_summ$ACEI,levels=c(0,1),labels=c("No ACE-i","ACE-i"))
    } else if(isTRUE(arb_present)) {
        demog_summ <- merge(demog_summ,med_arb_new,by="patient_id",all.x=TRUE) %>% dplyr::distinct(patient_id,.keep_all=TRUE)
        demog_summ[is.na(demog_summ)] <- 0
        demog_summ$ARB <- factor(demog_summ$ARB,levels=c(0,1),labels=c("No ARB","ARB"))
    }
    
    demog_summ[is.na(demog_summ)] <- 0
    demog_summ$aki <- 0
    demog_summ$aki[demog_summ$patient_id %in% labs_aki_summ$patient_id] <- 1
    demog_summ$severe <- factor(demog_summ$severe,levels=c(0,1),labels=c("Non-severe","Severe"))
    demog_summ$deceased <- factor(demog_summ$deceased,levels=c(0,1),labels=c("Alive","Deceased"))
    demog_summ$aki <- factor(demog_summ$aki,levels=c(0,1),labels=c("No AKI","AKI"))
    demog_summ$aki_kdigo_stage <- factor(demog_summ$aki_kdigo_stage,levels=c(0,1,2,3),labels=c("No AKI","Stage 1","Stage 2","Stage 3"))
    demog_summ[comorbid_list] <- lapply(demog_summ[comorbid_list],factor)
    demog_summ <- demog_summ %>% dplyr::distinct()
    
    message(paste0(c("Obfuscation cutoff: ",obfuscation_value)))
    # Obfuscation requirements by certain sites
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        comorbid_demog_summ_tmp <- vector(mode="list",length=length(comorbid_list))
        for(i in 1:length(comorbid_list)) {
            demog_summ_tmp1 <- demog_summ[,c("patient_id","aki",comorbid_list[i])]
            demog_summ_tmp2 <- demog_summ_tmp1 %>% dplyr::group_by(aki) %>% dplyr::count(get(comorbid_list[i]))
            message(paste0(c("Performing demographics table comorbid filtering for: ",comorbid_list[i]," with lowest count ",min(demog_summ_tmp2$n)," and obfuscation cutoff ",obfuscation_value)))
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
    
    if(isTRUE(acei_present) & isTRUE(arb_present)) {
        med_summ <- c(med_summ,"raas_new")
    } else if(isTRUE(acei_present)) {
        med_summ <- c(med_summ,"ACEI")
    } else if(isTRUE(arb_present)) {
        med_summ <- c(med_summ,"ARB")
    }
    
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        med_summ <- c(med_summ,"covid_rx")
    }
    
    table_one_vars <- c("sex","age_group","race","severe","deceased","aki_kdigo_stage",comorbid_demog_summ,med_summ)
    if(isTRUE(meld_analysis_valid)) {
        table_one_meld_vars <- c("sex","age_group","race","aki","meld_admit_severe","deceased","aki_kdigo_stage",comorbid_demog_summ,med_summ)
    }
    #capture.output(summary(table_one),file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne_Missingness.txt")))
    demog_meld_files <- NULL
    if(isTRUE(cirrhosis_present)) {
        demog_cld_summ <- demog_summ %>% dplyr::filter(cld == 1)
    }
    
    # Create obfuscated table one for sites which require it
    if(isTRUE(is_obfuscated) & !is.null(obfuscation_value)) {
        # obfuscated_table <- data.frame(export_table_one)
        # var_names <- rownames(obfuscated_table)
        # obfuscated_table <- within(obfuscated_table,No.AKI <- data.frame(do.call('rbind',strsplit(as.character(No.AKI),' (',fixed=T))))
        # obfuscated_table <- within(obfuscated_table,AKI <- data.frame(do.call('rbind',strsplit(as.character(AKI),' (',fixed=T))))
        # obfuscated_table <- as.data.frame(gsub("[),]","",as.matrix(obfuscated_table)))
        # obfuscated_table[c(2:6)] <- lapply(obfuscated_table[c(2:6)],as.numeric)
        # obfuscated_table <- obfuscated_table[,-7]
        # colnames(obfuscated_table) <- c("level","NoAKI_n","NoAKI_perc","AKI_n","AKI_perc","p")
        # obfuscated_table$names <- var_names
        # obfuscated_table <- lapply(obfuscated_table,function(x) { replace(x,grep("[X]",x),NA)})
        # obfuscated_table$names <- zoo::na.locf(obfuscated_table$names)
        # obfuscated_table <- as.data.frame(obfuscated_table)
        # obfuscated_table <- obfuscated_table %>% dplyr::select(names,level,NoAKI_n,NoAKI_perc,AKI_n,AKI_perc,p)
        # obfuscated_table$NoAKI_n[obfuscated_table$NoAKI_n < obfuscation_value] <- NA
        # obfuscated_table$AKI_n[obfuscated_table$AKI_n < obfuscation_value] <- NA
        # # obfuscated_table <- obfuscated_table %>% dplyr::mutate(remove = ifelse(NoAKI_n < obfuscation_value | AKI_n < obfuscation_value,1,0))
        # # obfuscated_table <- obfuscated_table %>% dplyr::filter(remove == 0)
        # obfuscated_table$names <- stringr::str_remove(obfuscated_table$names,stringr::fixed("...."))
        # # obfuscated_table <- obfuscated_table[,-8]
        
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
                    message(paste0(c("Attempting Fisher's test for ",table_one_vars[i])))
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
        write.csv(demog_obf,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne_obfuscated.csv")),row.names=F,na="NA")
        
        tryCatch({
            if(isTRUE(cirrhosis_present)) {
                # First generate the subgroup table stratified by AKI for cirrhotic patients only
                message("Creating temporary demographics table for cirrhotic patients")
                demog_cld_obf <- demog_cld_summ %>% dplyr::group_by(aki) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "aki",values_from = "n")
                demog_cld_obf$category <- "n"
                demog_cld_obf <- demog_cld_obf[,c(3,1,2)]
                demog_cld_obf$p_val <- NA
                colnames(demog_cld_obf) <- c("category","No_AKI","AKI","p_val")
                message("Filtering variables for cirrhotic patient demographics table...")
                for(i in 1:length(table_one_vars)) {
                    try({
                        tmp <- demog_cld_summ %>% dplyr::group_by(aki) %>% dplyr::count(get(table_one_vars[i])) %>% tidyr::pivot_wider(names_from="aki",values_from = "n")
                        tmp[is.na(tmp)] <- 0
                        colnames(tmp) <- c("category","No_AKI","AKI")
                        tmp <- tmp %>% dplyr::mutate(category = paste0(table_one_vars[i],"_",category))
                        
                        # tryCatch statements attempt to catch instances where the Fisher's test may fail due to insufficient convergent cycles
                        p_value <- tryCatch({
                            message(paste0(c("Attempting Fisher's test for ",table_one_vars[i])))
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
                message("Computing counts and percentages")
                demog_cld_obf$No_AKI[demog_cld_obf$No_AKI < obfuscation_value] <- 0
                demog_cld_obf$AKI[demog_cld_obf$AKI < obfuscation_value] <- 0
                demog_cld_obf <- demog_cld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(total = No_AKI + AKI) %>% dplyr::ungroup()
                no_aki_total <- demog_cld_obf$No_AKI[1]
                aki_total <- demog_cld_obf$AKI[1]
                total_pop <- demog_cld_obf$total[1]
                demog_cld_obf <- demog_cld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(No_AKI_perc = No_AKI / no_aki_total * 100, AKI_perc = AKI/aki_total * 100, total_perc = total/total_pop) %>% dplyr::ungroup()
                demog_meld_files <- "demog_cld_obf"
                message("No issues with creating demog_cld_obf")
            }
        },error = function(e) {
            message("Having issues generating demographic table for cirrhosis patients only. Check for any error messages that appear.")
        })
        try({
            if(isTRUE(meld_analysis_valid)) {
                # If possible to split by MELD, then generate the demographics table split by MELD score cutoff 20
                message("MELD analysis possible. Generating demographics table for cirrhotics stratified by MELD.")
                demog_meld_summ <- demog_cld_summ
                demog_meld_obf <- demog_meld_summ %>% dplyr::group_by(meld_admit_severe) %>% dplyr::count() %>% tidyr::pivot_wider(names_from = "meld_admit_severe",values_from = "n")
                demog_meld_obf$category <- "n"
                demog_meld_obf <- demog_meld_obf[,c(3,1,2)]
                demog_meld_obf$p_val <- NA
                colnames(demog_meld_obf) <- c("category","MELD_less20","MELD_20ormore","p_val")
                message("Performing filtering of variables for demographics table for cirrhotics by MELD")
                for(i in 1:length(table_one_meld_vars)) {
                    try({
                        tmp <- demog_meld_summ %>% dplyr::group_by(meld_admit_severe) %>% dplyr::count(get(table_one_meld_vars[i])) %>% tidyr::pivot_wider(names_from="meld_admit_severe",values_from = "n")
                        tmp[is.na(tmp)] <- 0
                        colnames(tmp) <- c("category","MELD_less20","MELD_20ormore")
                        tmp <- tmp %>% dplyr::mutate(category = paste0(table_one_meld_vars[i],"|",category))
                        
                        # tryCatch statements attempt to catch instances where the Fisher's test may fail due to insufficient convergent cycles
                        p_value <- tryCatch({
                            message(paste0(c("Attempting Fisher's test for ",table_one_meld_vars[i])))
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
                message("Computing admission labs averages")
                #Calculate stats of admission labs
                lab_meld_list <- NULL
                lab_anova <- NULL
                lab_meld_stats <- NULL
                tryCatch({
                    lab_meld_stats <- demog_meld_summ %>% dplyr::group_by(meld_admit_severe,aki) %>% dplyr::summarise(
                        mean_admit_ast = mean(first_ast,na.rm=TRUE),
                        mean_admit_alt = mean(first_alt,na.rm=TRUE),
                        mean_admit_bil = mean(first_bil,na.rm=TRUE),
                        mean_admit_inr = mean(first_inr,na.rm=TRUE),
                        mean_admit_alb = mean(first_alb,na.rm=TRUE),
                        sd_admit_ast = sd(first_ast,na.rm=TRUE),
                        sd_admit_alt = sd(first_alt,na.rm=TRUE),
                        sd_admit_bil = sd(first_bil,na.rm=TRUE),
                        sd_admit_inr = sd(first_inr,na.rm=TRUE),
                        sd_admit_alb = sd(first_alb,na.rm=TRUE),
                        n_admit_ast = sum(!is.na(first_ast)),
                        n_admit_alt = sum(!is.na(first_alt)),
                        n_admit_bil = sum(!is.na(first_bil)),
                        n_admit_inr = sum(!is.na(first_inr)),
                        n_admit_alb = sum(!is.na(first_alb))
                    ) %>% dplyr::ungroup() %>% dplyr::arrange(meld_admit_severe,aki)
                    ast_anova <- stats::aov(first_ast ~ meld_admit_severe * aki,data=demog_meld_summ)
                    alt_anova <- stats::aov(first_alt ~ meld_admit_severe * aki,data=demog_meld_summ)
                    bil_anova <- stats::aov(first_bil ~ meld_admit_severe * aki,data=demog_meld_summ)
                    inr_anova <- stats::aov(first_inr ~ meld_admit_severe * aki,data=demog_meld_summ)
                    alb_anova <- stats::aov(first_alb ~ meld_admit_severe * aki,data=demog_meld_summ)
                    lab_meld_list <- c("first_ast","first_alt","first_bil","first_inr","first_alb")
                    lab_anova <- c("ast_anova","alt_anova",'bil_anova","inr_anova',"alb_anova")
                    lab_meld_stats <- "lab_meld_stats"
                }, error = function(e) {
                    message("Error in processing labs. Check error messages.")
                })
                try({
                    demog_meld_obf$MELD_less20[demog_meld_obf$MELD_less20 < obfuscation_value] <- 0
                    demog_meld_obf$MELD_20ormore[demog_meld_obf$MELD_20ormore < obfuscation_value] <- 0
                    demog_meld_obf <- demog_meld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(total = MELD_less20 + MELD_20ormore) %>% dplyr::ungroup()
                    meld_less20_total <- demog_meld_obf$MELD_less20[1]
                    meld_20ormore_total <- demog_meld_obf$MELD_20ormore[1]
                    total_pop <- demog_meld_obf$total[1]
                    demog_meld_obf <- demog_meld_obf %>% dplyr::group_by(category) %>% dplyr::mutate(MELD_less20_perc = MELD_less20 / meld_less20_total * 100, MELD_20ormore_perc = MELD_20ormore/meld_20ormore_total * 100, total_perc = total/total_pop * 100) %>% dplyr::ungroup()
                    message("Computed final demog_meld_obf table.")
                })
                
                if(obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
                    table_one_meld <- tableone::CreateTableOne(data=demog_meld_summ,vars=table_one_Meld_vars,strata="meld_admit_severe")
                    export_table_one_meld <- print(table_one_meld,showAllLevels=TRUE,formatOptions=list(big.mark=","))
                    if(exists("demog_meld_obf")) {
                        demog_meld_files <- c(demog_meld_files,"demog_meld_obf","table_one_meld","export_table_one_meld",lab_meld_stats,lab_anova)
                    } else {
                        demog_meld_files <- c(demog_meld_files,"table_one_meld","export_table_one_meld",lab_meld_stats,lab_anova)
                    }
                } else {
                    demog_meld_files <- c(demog_meld_files,"demog_meld_obf",lab_meld_stats,lab_anova)
                }
            }
        })
        
    }

    if(obfuscation_value == 0 | isTRUE(!is_obfuscated)) {
        table_one <- tableone::CreateTableOne(data=demog_summ,vars=table_one_vars,strata="aki")
        export_table_one <- print(table_one,showAllLevels=TRUE,formatOptions=list(big.mark=","))
        write.csv(export_table_one,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TableOne.csv")))
        try({
            if(isTRUE(cirrhosis_present)) {
                table_one_cld <- tableone::CreateTableOne(data=demog_cld_summ,vars=table_one_vars,strata="aki")
                export_table_one_cld <- print(table_one_cld,showAllLevels=TRUE,formatOptions=list(big.mark=","))
            }
        })
    }
    message("Attempting to save demographics tables for cirrhotic patients.")
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
        save(list=demog_meld_files,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_MELD_Cirrhosis.rda")),compress="bzip2")
    }
    message("TableOne with patient demographics should have been generated in CSV files at this point. Check for any errors.")
    
    demog_time_to_event_tmp <- demog_summ[,c("sex","age_group","race")] 
    demog_time_to_event_tmp <- data.table::as.data.table(lapply(demog_time_to_event_tmp,factor))
    demog_time_to_event_tmp <- data.table::as.data.table(demog_time_to_event_tmp)[,sapply(demog_time_to_event_tmp,function(col) nlevels(col) > 1),with=FALSE] 
    demog_list <- colnames(demog_time_to_event_tmp)
    demog_time_to_event <- demog_summ[,c("patient_id",demog_list)] %>% dplyr::group_by(patient_id) %>% dplyr::mutate(age_group = dplyr::if_else(age_group == "70to79" | age_group == "80plus","70_and_above","below_70")) %>% dplyr::ungroup()
    
    # Deals with the special case where there is a sex category of others (e.g. KUMC)
    demog_time_to_event <- demog_time_to_event %>% dplyr::group_by(patient_id) %>% dplyr::mutate(sex = dplyr::if_else(sex == "male","male","female")) %>% dplyr::ungroup()
    
    
    # ================================================================================================================================
    # Figure 1(b) Comparing serum creatinine trends of severe and non-severe patients, with or without remdesivir/lopinavir+ritonavir
    # ================================================================================================================================
    if(isTRUE(covid19antiviral_present) | isTRUE(remdesivir_present)) {
        message("Now preparing plots for serum Cr trends in experimental COVID-19 treatment...")
        # Plotting from peak
        peak_trend_covidviral <- peak_trend %>% dplyr::select(patient_id,covidrx_grp,time_from_peak,ratio) %>% dplyr::distinct(patient_id,time_from_peak,.keep_all=TRUE)
        peak_cr_covidviral_summ <- peak_trend_covidviral %>% dplyr::group_by(covidrx_grp,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio),sem_ratio = sd(ratio)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
        covidrx_grp_label <- data.table::data.table(c(1:4),c("Non-severe, no novel COVID-19 treatment","Non-severe, with novel COVID-19 treatment","Severe, no novel COVID-19 treatment","Severe, with novel COVID-19 treatment"))
        colnames(covidrx_grp_label) <- c("covidrx_grp","covidrx_label")
        peak_cr_covidviral_summ <- merge(peak_cr_covidviral_summ,covidrx_grp_label,by="covidrx_grp",all.x=TRUE)
        if(isTRUE(is_obfuscated)) {
            peak_cr_covidviral_summ  <- peak_cr_covidviral_summ[peak_cr_covidviral_summ$n >= obfuscation_value,]
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
        
        cr_from_covidrx_trend <- cr_from_covidrx_trend %>% dplyr::select(patient_id,severe,covidrx_grp,days_since_admission,covid_rx_start,peak_cr_time,value,min_cr_90d,min_cr_retro_7day) %>% dplyr::distinct(patient_id,days_since_admission,.keep_all=TRUE)
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
        if(isTRUE(is_obfuscated)) {
            cr_from_covidrx_summ  <- cr_from_covidrx_summ[cr_from_covidrx_summ$n >= obfuscation_value,]
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
    
    # If user wishes to customize the Cox PH equations used for recovery and death analysis, we will read in
    # custom files specifying the factors to restrict analyses to.
    # This may be helpful in cases where there may not be enough events for certain factors, causing model
    # convergence to fail.
    # These should be listed as space-separated names in a file "such as the example shown below:
    # age sex ckd cld htn hld ihd
    
    restrict_list <- ""
    model1 <- c("age_group","sex","severe","aki_kdigo_final","ckd","htn","hld","ihd","cld")
    model2 <- c("age_group","sex","severe","bronchiectasis","copd","rheum","vte")
    model3 <- c("age_group","sex","severe","COAGA","COAGB","covid_rx")
    model4 <- c("age_group","sex","severe","aki_kdigo_final","ckd")
    
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
    
    message("\nDoing initial filter for medications with more than one factor level.")
    med_recovery_list <- c("COAGA","COAGB","covid_rx")
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
    
    # Kaplan Meier plot
    message("Generating Kaplan-Meier plots...")
    recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe")
    #surv_recover <- survival::Surv(time=aki_index_recovery$time_to_ratio1.25,event=aki_index_recovery$recover_1.25x)
    fit_km_recover <- survminer::surv_fit(recoverPlotFormula, data=aki_index_recovery)
    plot_recover <- survminer::ggsurvplot(fit_km_recover,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
    plot_recover_summ <- survminer::surv_summary(fit_km_recover,data=aki_index_recovery)
    plot_recover_summ_table <- plot_recover$data.survtable
    write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_PlotSummStats.csv")),row.names=TRUE)
    write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_recover_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe.png")),plot=print(plot_recover),width=12,height=12,units="cm")
    
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
        #recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model1,"severe * aki_kdigo_final"),collapse="+")))
        recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model1,collapse="+")))
        message(paste("Formula for Model 1: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model1,"severe * aki_kdigo_final"),collapse="+")))
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
    
    try({
        recovery_model4 <- c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list,med_recovery_list)
        recovery_model4 <- recovery_model4[recovery_model4 %in% model4]
        recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model4,"severe * aki_kdigo_final"),collapse="+")))
        #recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model4,collapse="+")))
        message(paste("Formula for Model 4: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model4,"severe * aki_kdigo_final"),collapse="+")))
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
    # deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c("severe","aki_kdigo_final",demog_recovery_list,comorbid_recovery_list),sep="+")))
    # surv_death_aki_only<- survival::Surv(time=aki_index_recovery$time_to_death_km,event=aki_index_recovery$deceased)
    fit_death_aki_only <- survminer::surv_fit(deathPlotFormula, data=aki_index_recovery)
    plot_death_aki_only <- survminer::ggsurvplot(fit_death_aki_only,data=aki_index_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
    plot_death_aki_only_summ <- survminer::surv_summary(fit_death_aki_only,data=aki_index_recovery)
    plot_death_aki_only_summ_table <- plot_death_aki_only$data.survtable
    write.csv(fit_death_aki_only$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_PlotSummStats.csv")),row.names=TRUE)
    write.csv(plot_death_aki_only_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Plot.csv")),row.names=FALSE)
    write.csv(plot_death_aki_only_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_Severe.png")),plot=print(plot_death_aki_only),width=12,height=12,units="cm")
    # coxph_death_aki_only <- survival::coxph(deathCoxPHFormula, data=aki_index_recovery)
    # coxph_death_aki_only_plot <- survminer::ggforest(coxph_death_aki_only,data=aki_index_recovery)
    # coxph_death_aki_only_summ <- summary(coxph_death_aki_only) 
    # write.csv(coxph_death_aki_only_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH.csv")),row.names=TRUE)
    # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIOnly_CoxPH.png")),plot=print(coxph_death_aki_only_plot),width=20,height=20,units="cm")
    #
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
        # deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_aki_only_model1,"severe * aki_kdigo_final"),collapse="+")))
        # message("Formula for Model 1: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_aki_only_model1,"severe * aki_kdigo_final"),collapse="+")))
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
        deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_aki_only_model4,"severe * aki_kdigo_final"),collapse="+")))
        message("Formula for Model 4: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_aki_only_model4,"severe * aki_kdigo_final"),collapse="+")))
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
    #deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~",paste(c("is_aki","severe",demog_death_list,comorbid_death_list),collapse="+")))
    
    #surv_death <- survival::Surv(time=aki_index_death$time_to_death_km,event=aki_index_death$deceased)
    fit_death <- survminer::surv_fit(deathPlotFormula, data=aki_index_death)
    plot_death <- survminer::ggsurvplot(fit_death,data=aki_index_death,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw())
    plot_death_summ <- survminer::surv_summary(fit_death,data=aki_index_death)
    plot_death_summ_table <- plot_death$data.survtable
    #write.csv(fit_death$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_PlotSummStats.csv")),row.names=TRUE)
    write.csv(plot_death_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_Plot.csv")),row.names=FALSE)
    #write.csv(plot_death_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_Table.csv")),row.names=FALSE)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI.png")),plot=print(plot_death),width=12,height=12,units="cm")
    # coxph_death <- survival::coxph(deathCoxPHFormula, data=aki_index_death)
    # coxph_death_plot <- survminer::ggforest(coxph_death,data=aki_index_death)
    # coxph_death_summ <- summary(coxph_death) 
    # write.csv(coxph_death_summ$coefficients,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_AKIvsNonAKI_CoxPH.csv")),row.names=TRUE)
    # ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CoxPH.png")),plot=print(coxph_death_plot),width=20,height=20,units="cm")
    # 
    
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
        # deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_model1,"severe * aki_kdigo_final"),collapse="+")))
        # message("Formula for Model 1: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_model1,"severe * aki_kdigo_final"),collapse="+")))
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
        # save(death_model1,coxph_death_all1_summ,coxph_death_all1_hr,coxph_death_all1_stats1, coxph_death_all1_stats2,coxph_death_all1_plot,file =file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model1.rdata")))
        ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model1.png")),plot=print(coxph_death_all1_plot),width=20,height=20,units="cm")
    })
    
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
        # save(death_model2,coxph_death_all2_summ,coxph_death_all2_hr,coxph_death_all2_stats1, coxph_death_all2_stats2,coxph_death_all2_plot,file =file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model2.rdata")))
        ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model2.png")),plot=print(coxph_death_all2_plot),width=20,height=20,units="cm")
    })
    
    message("Generating Model 3 (Time to death, all patients)...")
    try({
        death_model3 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
        death_model3 <- death_model3[death_model3 %in% model3]
        deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model3,collapse="+")))
        message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model3,collapse="+")))
        # if(("covid_rx" %in% death_model3) == TRUE) {
        #     deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_model3,"severe:covid_rx"),collapse="+")))
        #     message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_model3,"severe:covid_rx"),collapse="+")))
        # } else {
        #     deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model3,collapse="+")))
        #     message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(death_model3,collapse="+")))
        # }
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
        # save(death_model3,coxph_death_all3_summ,coxph_death_all3_hr,coxph_death_all3_stats1, coxph_death_all3_stats2,coxph_death_all3_plot,file =file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model3.rdata")))
        ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model3.png")),plot=print(coxph_death_all3_plot),width=20,height=20,units="cm")
    })
    
    message("Generating Model 4 (Time to death, all patients)...")
    try({
        death_model4 <- c("severe","aki_kdigo_final",demog_death_list,comorbid_death_list,med_death_list)
        death_model4 <- death_model4[death_model4 %in% model4]
        deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_model4,"severe * aki_kdigo_final"),collapse="+")))
        message("Formula for Model 4: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(c(death_model4,"severe * aki_kdigo_final"),collapse="+")))
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
        #save(death_model4,coxph_death_all4_summ,coxph_death_all4_hr,coxph_death_all4_stats1, coxph_death_all4_stats2,coxph_death_all4_plot,file =file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model4.rdata")))
        ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_All_CoxPH_Model4.png")),plot=print(coxph_death_all4_plot),width=20,height=20,units="cm")
    })
    
    message("If you are getting any errors with model generation - do note that it may actually be normal to get errors\nif your site numbers are low (especially for model 3). Please check your data to see if the appropriate\nnumber of events occur for each factor level.")
    
    # ================================================
    # Part 3: Hepatorenal Syndrome Analyses
    # ================================================
    
    if(isTRUE(cirrhosis_present)) {
        message("Doing time-to-event analyses for cirrhotic patients using MELD scoring")
        cirrhotic_recovery <- aki_index_recovery %>% dplyr::filter(cld == 1)
        if(isTRUE(meld_analysis_valid)) {
            cirrhotic_recovery <- merge(cirrhotic_recovery,meld_severe_list[,c(1,3)],by="patient_id",all.x=TRUE) %>% dplyr::distinct()
            cirrhotic_recovery$meld_admit_severe[is.na(cirrhotic_recovery$meld_admit_severe)] <- 0
            cirrhotic_recovery$meld_admit_severe <- factor(cirrhotic_recovery$meld_admit_severe,levels=c(0,1),labels = c("MELD < 20","MELD ≥ 20"))
        }
        
        message("Filtering factor list down further for CoxPH models...")
        comorbid_recovery_list <- comorbid_recovery_list[comorbid_recovery_list != "cld"]
        recovery_tmp <- cirrhotic_recovery[,c("patient_id","recover_1.25x",demog_list,comorbid_recovery_list,med_recovery_list)] %>% as.data.frame()
        
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
        readr::write_lines(variable_list_output,file.path(getProjectOutputDirectory(), paste0(currSiteId, "_cirrhotic_custom_equation.txt")),append=F)
        
        message("Now proceeding to time-to-Cr recovery analysis...")
        # Now run the actual time-to-event analysis
        
        # Kaplan Meier plot
        message("Generating Kaplan-Meier plots...")
        recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ severe")
        #surv_recover <- survival::Surv(time=cirrhotic_recovery$time_to_ratio1.25,event=cirrhotic_recovery$recover_1.25x)
        fit_km_cirrhotic_recover <- survminer::surv_fit(recoverPlotFormula, data=cirrhotic_recovery)
        fit_km_cirrhotic_recover_table <- fit_km_recover$table
        plot_cirrhotic_recover <- survminer::ggsurvplot(fit_km_recover,data=cirrhotic_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
        plot_cirrhotic_recover_summ <- survminer::surv_summary(fit_km_recover,data=cirrhotic_recovery)
        plot_cirrhotic_recover_summ_table <- plot_recover$data.survtable
        # write.csv(fit_km_recover$table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_PlotSummStats.csv")),row.names=TRUE)
        # write.csv(plot_recover_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Plot.csv")),row.names=FALSE)
        # write.csv(plot_recover_summ_table,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Recover_Severe_Table.csv")),row.names=FALSE)
        ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Recover_Severe.png")),plot=print(plot_cirrhotic_recover),width=12,height=12,units="cm")
        
        cirrhotic_files <- c("fit_km_cirrhotic_recover_table","plot_cirrhotic_recover_summ","plot_cirrhotic_recover_summ_table","plot_cirrhotic_recover")
        if(isTRUE(meld_analysis_valid)) {
            try({
                recoverPlotFormula <- as.formula("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ meld_admit_severe")
                fit_km_meld_recover <- survminer::surv_fit(recoverPlotFormula, data=cirrhotic_recovery)
                fit_km_meld_recover_table <- fit_km_recover$table
                plot_meld_recover <- survminer::ggsurvplot(fit_km_recover,data=cirrhotic_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),fun="event",xlim=c(0,90),break.x.by=30)
                plot_meld_recover_summ <- survminer::surv_summary(fit_km_recover,data=cirrhotic_recovery)
                plot_meld_recover_summ_table <- plot_recover$data.survtable
                ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_MELD_Recover_Severe.png")),plot=print(plot_meld_recover),width=12,height=12,units="cm")
                cirrhotic_files <- c(cirrhotic_files,"fit_km_meld_recover_table","plot_meld_recover_summ","plot_meld_recover_summ_table","plot_meld_recover")
            })
        }
        
        # CoxPH model
        # Generate univariate analyses first
        meld_var <- NULL
        if(isTRUE(meld_analysis_valid)) {
            meld_var <- "meld_admit_severe"
        }
        message("Generating univariate Cox PH models (time to recovery, cirrhotic AKI patients only)...")
        univ_formulas <- sapply(c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list), function(x) as.formula(paste('survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ', x)))
        univ_models <- tryCatch(
            lapply(univ_formulas, function(x){survival::coxph(x, data = cirrhotic_recovery)}),
            error = function(c) "Problem generating univariate models."
        )
        # Extract data 
        
        try({
            univ_results <- lapply(univ_models,function(x){
                x <- summary(x)
                return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
            })
            univ_results_recover_cirrhotic <- do.call("rbind",univ_results)
            cirrhotic_files <- c(cirrhotic_files,"univ_results_recover_cirrhotic")
        })
        
        message("\nGenerating Model 1 (time to recovery, cirrhotic AKI patients only)...")
        try({
            recovery_model1 <- c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list)
            recovery_model1 <- recovery_model1[recovery_model1 %in% model1]
            #recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model1,"severe * aki_kdigo_final"),collapse="+")))
            recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model1,collapse="+")))
            message(paste("Formula for Model 1: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(c(recovery_model1,"severe * aki_kdigo_final"),collapse="+")))
            coxph_cirrhotic_recover1 <- survival::coxph(recoverCoxPHFormula, data=cirrhotic_recovery)
            coxph_cirrhotic_recover1_summ <- summary(coxph_cirrhotic_recover1) 
            print(coxph_cirrhotic_recover1_summ)
            coxph_cirrhotic_recover1_hr <- cbind(coxph_cirrhotic_recover1_summ$coefficients,coxph_cirrhotic_recover1_summ$conf.int)[,-c(6,7)]
            coxph_cirrhotic_recover1_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_cirrhotic_recover1_summ$logtest,coxph_cirrhotic_recover1_summ$sctest,coxph_cirrhotic_recover1_summ$waldtest))
            coxph_cirrhotic_recover1_stats2 <- rbind(data.table::as.data.table(coxph_cirrhotic_recover1_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_cirrhotic_recover1_summ$rsq,keep.rownames = T))
            coxph_cirrhotic_recover1_plot <- survminer::ggforest(coxph_cirrhotic_recover1,data=cirrhotic_recovery)
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Recover_CoxPH_Model1.png")),plot=print(coxph_cirrhotic_recover1_plot),width=20,height=20,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"coxph_cirrhotic_recover1_summ","coxph_cirrhotic_recover1_hr","coxph_cirrhotic_recover1_stats1","coxph_cirrhotic_recover1_stats2","coxph_cirrhotic_recover1_plot")
        })
        
        message("\nGenerating Model 2 (time to recovery, cirrhotic AKI patients only)...")
        try({
            recovery_model2 <- c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list)
            recovery_model2 <- recovery_model2[recovery_model2 %in% model2]
            recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model2,collapse="+")))
            message(paste("Formula for Model 2: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model2,collapse="+")))
            coxph_cirrhotic_recover2 <- survival::coxph(recoverCoxPHFormula, data=cirrhotic_recovery)
            coxph_cirrhotic_recover2_summ <- summary(coxph_cirrhotic_recover2) 
            print(coxph_cirrhotic_recover2_summ)
            coxph_cirrhotic_recover2_hr <- cbind(coxph_cirrhotic_recover2_summ$coefficients,coxph_cirrhotic_recover2_summ$conf.int)[,-c(6,7)]
            coxph_cirrhotic_recover2_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_cirrhotic_recover2_summ$logtest,coxph_cirrhotic_recover2_summ$sctest,coxph_cirrhotic_recover2_summ$waldtest))
            coxph_cirrhotic_recover2_stats2 <- rbind(data.table::as.data.table(coxph_cirrhotic_recover2_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_cirrhotic_recover2_summ$rsq,keep.rownames = T))
            coxph_cirrhotic_recover2_plot <- survminer::ggforest(coxph_cirrhotic_recover2,data=cirrhotic_recovery)
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Recover_CoxPH_Model2.png")),plot=print(coxph_cirrhotic_recover2_plot),width=20,height=20,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"coxph_cirrhotic_recover2_summ","coxph_cirrhotic_recover2_hr","coxph_cirrhotic_recover2_stats1","coxph_cirrhotic_recover2_stats2","coxph_cirrhotic_recover2_plot")
        })
        
        message("\nGenerating Model 3 (time to recovery, cirrhotic AKI patients only)...")
        try({
            recovery_model3 <- c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list)
            recovery_model3 <- recovery_model3[recovery_model3 %in% model3]
            recoverCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
            message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
            message(paste("Formula for Model 3: survival::Surv(time=time_to_ratio1.25,event=recover_1.25x) ~ ",paste(recovery_model3,collapse="+")))
            coxph_cirrhotic_recover3 <- survival::coxph(recoverCoxPHFormula, data=cirrhotic_recovery)
            coxph_cirrhotic_recover3_summ <- summary(coxph_cirrhotic_recover3) 
            print(coxph_cirrhotic_recover3_summ)
            coxph_cirrhotic_recover3_hr <- cbind(coxph_cirrhotic_recover3_summ$coefficients,coxph_cirrhotic_recover3_summ$conf.int)[,-c(6,7)]
            coxph_cirrhotic_recover3_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_cirrhotic_recover3_summ$logtest,coxph_cirrhotic_recover3_summ$sctest,coxph_cirrhotic_recover3_summ$waldtest))
            coxph_cirrhotic_recover3_stats2 <- rbind(data.table::as.data.table(coxph_cirrhotic_recover3_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_cirrhotic_recover3_summ$rsq,keep.rownames = T))
            coxph_cirrhotic_recover3_plot <- survminer::ggforest(coxph_cirrhotic_recover3,data=cirrhotic_recovery)
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Recover_CoxPH_Model3.png")),plot=print(coxph_cirrhotic_recover3_plot),width=20,height=20,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"coxph_cirrhotic_recover3_summ","coxph_cirrhotic_recover3_hr","coxph_cirrhotic_recover3_stats1","coxph_cirrhotic_recover3_stats2","coxph_cirrhotic_recover3_plot")
        })
        
        message("If you are getting any errors with model generation - do note that it may actually be normal to get errors\nif your site numbers are low (especially for model 3). Please check your data to see if the appropriate\nnumber of events occur for each factor level.")
        
        message("Now proceeding to time-to-death analysis for AKI patients only...")
        
        deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ severe")
        fit_death_cirrhotic_aki <- survminer::surv_fit(deathPlotFormula, data=cirrhotic_recovery)
        fit_death_cirrhotic_aki_table <- fit_death_cirrhotic_aki$table
        plot_death_cirrhotic_aki <- survminer::ggsurvplot(fit_death_cirrhotic_aki,data=cirrhotic_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
        plot_death_cirrhotic_aki_summ <- survminer::surv_summary(fit_death_cirrhotic_aki,data=cirrhotic_recovery)
        plot_death_cirrhotic_aki_summ_table <- plot_death_cirrhotic_aki$data.survtable
        ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_CirrhoticAKI_Severe.png")),plot=print(plot_death_cirrhotic_aki),width=12,height=12,units="cm")
        cirrhotic_files <- c(cirrhotic_files,"fit_death_cirrhotic_aki_table","plot_death_cirrhotic_aki","plot_death_cirrhotic_aki_summ","plot_death_cirrhotic_aki_summ_table")
        
        try({
            deathPlotFormula <- as.formula("survival::Surv(time=time_to_death_km,event=deceased) ~ meld_admit_severe")
            fit_death_meld_aki <- survminer::surv_fit(deathPlotFormula, data=cirrhotic_recovery)
            fit_death_meld_aki_table <- fit_death_meld_aki$table
            plot_death_meld_aki <- survminer::ggsurvplot(fit_death_meld_aki,data=cirrhotic_recovery,pval=TRUE,conf.int=TRUE,risk.table=TRUE,risk.table.col = "strata", linetype = "strata",surv.median.line = "hv",ggtheme = ggplot2::theme_bw(),xlim=c(0,365),break.x.by=30)
            plot_death_meld_aki_summ <- survminer::surv_summary(fit_death_meld_aki,data=cirrhotic_recovery)
            plot_death_meld_aki_summ_table <- plot_death_meld_aki$data.survtable
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Death_MELD_Severe.png")),plot=print(plot_death_meld_aki),width=12,height=12,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"fit_death_meld_aki_table","plot_death_meld_aki","plot_death_meld_aki_summ","plot_death_meld_aki_summ_table")
        })
        
        message("Generating univariate Cox PH models (Time to death, Cirrhotic AKI patients only)...")
        
        univ_formulas <- sapply(c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list), function(x) as.formula(paste('survival::Surv(time=time_to_death_km,event=deceased) ~ ', x)))
        univ_models <- tryCatch(
            lapply(univ_formulas, function(x){survival::coxph(x, data = cirrhotic_recovery)}),
            error = "Problem generating univariate models."
        )
        # Extract data 
        try({
            univ_results <- lapply(univ_models,function(x){
                x <- summary(x)
                return(cbind(x$coefficients,x$conf.int)[,-c(6,7)])
            })
            univ_results_death_cirrhotic<- do.call("rbind",univ_results)
            cirrhotic_files <- c(cirrhotic_files,"univ_results_death_cirrhotic")
        })
        
        message("Generating Model 1 (Time to death, cirrhotic AKI patients only)...")
        try({
            cirrhotic_death_aki_model1 <- c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list)
            cirrhotic_death_aki_model1 <- cirrhotic_death_aki_model1[cirrhotic_death_aki_model1 %in% model1]
            deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(cirrhotic_death_aki_model1,collapse="+")))
            message("Formula for Model 1: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(cirrhotic_death_aki_model1,collapse="+")))
            coxph_cirrhotic_death1 <- survival::coxph(deathCoxPHFormula, data=cirrhotic_recovery)
            coxph_cirrhotic_death1_summ <- summary(coxph_cirrhotic_death1)
            print(coxph_cirrhotic_death1_summ)
            coxph_cirrhotic_death1_hr <- cbind(coxph_cirrhotic_death1_summ$coefficients,coxph_cirrhotic_death1_summ$conf.int)[,-c(6,7)]
            coxph_cirrhotic_death1_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_cirrhotic_death1_summ$logtest,coxph_cirrhotic_death1_summ$sctest,coxph_cirrhotic_death1_summ$waldtest))
            coxph_cirrhotic_death1_stats2 <- rbind(data.table::as.data.table(coxph_cirrhotic_death1_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_cirrhotic_death1_summ$rsq,keep.rownames = T))
            coxph_cirrhotic_death1_plot <- survminer::ggforest(coxph_cirrhotic_death1,data=cirrhotic_recovery)
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Death_CoxPH_Model1.png")),plot=print(coxph_cirrhotic_death1_plot),width=20,height=20,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"coxph_cirrhotic_death1_summ","coxph_cirrhotic_death1_hr","coxph_cirrhotic_death1_stats1","coxph_cirrhotic_death1_stats2","coxph_cirrhotic_death1_plot")
        })
        
        message("Generating Model 2 (Time to death, cirrhotic AKI patients only)...")
        try({
            cirrhotic_death_aki_model2 <- c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list)
            cirrhotic_death_aki_model2 <- cirrhotic_death_aki_model2[cirrhotic_death_aki_model2 %in% model2]
            deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(cirrhotic_death_aki_model2,collapse="+")))
            message("Formula for Model 2: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(cirrhotic_death_aki_model2,collapse="+")))
            coxph_cirrhotic_death2 <- survival::coxph(deathCoxPHFormula, data=cirrhotic_recovery)
            coxph_cirrhotic_death2_summ <- summary(coxph_cirrhotic_death2) 
            print(coxph_cirrhotic_death2_summ)
            coxph_cirrhotic_death2_hr <- cbind(coxph_cirrhotic_death2_summ$coefficients,coxph_cirrhotic_death2_summ$conf.int)[,-c(6,7)]
            coxph_cirrhotic_death2_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_cirrhotic_death2_summ$logtest,coxph_cirrhotic_death2_summ$sctest,coxph_cirrhotic_death2_summ$waldtest))
            coxph_cirrhotic_death2_stats2 <- rbind(data.table::as.data.table(coxph_cirrhotic_death2_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_cirrhotic_death2_summ$rsq,keep.rownames = T))
            coxph_cirrhotic_death2_plot <- survminer::ggforest(coxph_cirrhotic_death2,data=cirrhotic_recovery)
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Death_CoxPH_Model2.png")),plot=print(coxph_cirrhotic_death2_plot),width=20,height=20,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"coxph_cirrhotic_death2_summ","coxph_cirrhotic_death2_hr","coxph_cirrhotic_death2_stats1","coxph_cirrhotic_death2_stats2","coxph_cirrhotic_death2_plot")
        })
        
        message("Generating Model 3 (Time to death, cirrhotic AKI patients only)...")
        try({
            cirrhotic_death_aki_model3 <- c("severe","aki_kdigo_final",meld_var,demog_recovery_list,comorbid_recovery_list,med_recovery_list)
            cirrhotic_death_aki_model3 <- cirrhotic_death_aki_model3[cirrhotic_death_aki_model3 %in% model3]
            deathCoxPHFormula <- as.formula(paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(cirrhotic_death_aki_model3,collapse="+")))
            message("Formula for Model 3: ",paste("survival::Surv(time=time_to_death_km,event=deceased) ~ ",paste(cirrhotic_death_aki_model3,collapse="+")))
            coxph_cirrhotic_death3 <- survival::coxph(deathCoxPHFormula, data=cirrhotic_recovery)
            coxph_cirrhotic_death3_summ <- summary(coxph_cirrhotic_death3) 
            print(coxph_cirrhotic_death3_summ)
            coxph_cirrhotic_death3_hr <- cbind(coxph_cirrhotic_death3_summ$coefficients,coxph_cirrhotic_death3_summ$conf.int)[,-c(6,7)]
            coxph_cirrhotic_death3_stats1 <- cbind(c("logtest","sctest","waldtest"),rbind(coxph_cirrhotic_death3_summ$logtest,coxph_cirrhotic_death3_summ$sctest,coxph_cirrhotic_death3_summ$waldtest))
            coxph_cirrhotic_death3_stats2 <- rbind(data.table::as.data.table(coxph_cirrhotic_death3_summ$concordance,keep.rownames = T),data.table::as.data.table(coxph_cirrhotic_death3_summ$rsq,keep.rownames = T))
            coxph_cirrhotic_death3_plot <- survminer::ggforest(coxph_cirrhotic_death3,data=cirrhotic_recovery)
            ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_TimeToEvent_Cirrhotic_Death_CoxPH_Model3.png")),plot=print(coxph_cirrhotic_death3_plot),width=20,height=20,units="cm")
            cirrhotic_files <- c(cirrhotic_files,"coxph_cirrhotic_death3_summ","coxph_cirrhotic_death3_hr","coxph_cirrhotic_death3_stats1","coxph_cirrhotic_death3_stats2","coxph_cirrhotic_death3_plot")
        })
    }
    if(isTRUE(exists("cirrhotic_files"))) {
        if(!is.null(cirrhotic_files)) {
            save(list=cirrhotic_files,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_MELD_CLD_TimeToEvent.rda")),compress="bzip2")
        }
    }
    
    
    
    # # ================================================
    # # Part 4: Extending analyses to Thrombotic Events
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