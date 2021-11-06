#' Generates serum creatinine graphs for patients.
#' @param base_table peak_trend
#' @param obfuscation is_obfuscated
#' @param obfuscation_level obfuscation_value
#' @param kdigo_grade_table kdigo_grade
#' @param ckd_valid ckd_present
#' @param ckd_table ckd_list
#' @param labs_cr labs_cr_all 
#' @param use_index_baseline Indicates if the new definition of baseline sCr (-365 days to last day of index admission) is to be used.
#' @param index_admit_baselines baseline_cr_index_admit

generate_cr_graphs <- function(siteid,base_table,obfuscation,obfuscation_level,kdigo_grade_table,ckd_valid,ckd_table,labs_cr,use_index_baseline = FALSE, index_admit_baselines) {
  currSiteId <- siteid
  peak_trend <- base_table
  is_obfuscated <- obfuscation
  obfuscation_value <- obfuscation_level
  kdigo_grade <- kdigo_grade_table
  ckd_present <- ckd_valid
  ckd_list <- ckd_table
  labs_cr_all <- labs_cr
  baseline_cr_index_admit <- index_admit_baselines
  
  if(isTRUE(use_index_baseline)) {
    currSiteId <- paste0(currSiteId,"_IndexBaselineCr")
    peak_trend$ratio <- peak_trend$ratio_baseline_index
  }
  
  # First create a plot of the creatinine trends of AKI vs non-AKI patients from the day of the first AKI peak (or the highest Cr peak for non-AKI patients)
  peak_aki_vs_non_aki <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio,value) %>% dplyr::arrange(patient_id,severe,time_from_peak,ratio)
  peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% dplyr::group_by(patient_id,time_from_peak) %>% tidyr::fill(ratio,value,severe) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,time_from_peak,.keep_all=TRUE)
  colnames(peak_aki_vs_non_aki) <- c("patient_id","aki","time_from_peak","ratio","value")
  peak_aki_vs_non_aki <- peak_aki_vs_non_aki %>% dplyr::group_by(patient_id) %>% dplyr::mutate(aki = ifelse((aki == 2 | aki == 4),1,0)) %>% dplyr::ungroup()
  peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki %>% dplyr::group_by(aki,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  aki_label <- data.table::data.table(c(0,1),c("Non-AKI","AKI"))
  colnames(aki_label) <- c("aki","aki_label")
  peak_aki_vs_non_aki_summ <- merge(peak_aki_vs_non_aki_summ,aki_label,by="aki",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki_summ %>% dplyr::filter(n >= obfuscation_value) %>% dplyr::arrange(aki,time_from_peak)
    cat("\nObfuscating the AKI vs non-AKI graphs...")
    peak_aki_vs_non_aki_summ <- peak_aki_vs_non_aki_summ[peak_aki_vs_non_aki_summ$n >= obfuscation_value,]
  }
  write.csv(peak_aki_vs_non_aki_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_PeakCr_AKI_vs_NonAKI.csv")),row.names=FALSE)
  peak_aki_vs_non_aki_timeplot <- ggplot2::ggplot(peak_aki_vs_non_aki_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=aki_label))+ggplot2::geom_line(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio, color=factor(aki_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "AKI Group") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,6) + ggplot2::scale_color_manual(values=c("AKI"="#bc3c29","Non-AKI"="#0072b5")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_AKI_vs_NonAKI.png")),plot=peak_aki_vs_non_aki_timeplot,width=12,height=9,units="cm")
  peak_aki_vs_non_aki_timeplot_raw <- ggplot2::ggplot(peak_aki_vs_non_aki_summ,ggplot2::aes(x=time_from_peak,y=mean_value,group=aki_label))+ggplot2::geom_line(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(aki_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value, color=factor(aki_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr (mg/dL)", color = "AKI Group") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,6) + ggplot2::scale_color_manual(values=c("AKI"="#bc3c29","Non-AKI"="#0072b5")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_AKI_vs_NonAKI_RawCr.png")),plot=peak_aki_vs_non_aki_timeplot_raw,width=12,height=9,units="cm")
  cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients should have been generated.")
  
  
  # Now, derive our first table peak_trend_severe to compare across the different severity groups
  peak_trend_severe <- peak_trend %>% dplyr::select(patient_id,severe,time_from_peak,ratio,value) %>% dplyr::arrange(patient_id,severe,time_from_peak,ratio)
  peak_trend_severe <- peak_trend_severe %>% dplyr::group_by(patient_id,time_from_peak) %>% tidyr::fill(ratio,value,severe) %>% dplyr::ungroup() %>% dplyr::distinct(patient_id,time_from_peak,.keep_all=TRUE)
  # Headers: patient_id  severe  time_from_peak  ratio    value
  # peak_trend_severe <- peak_trend_severe %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
  # Merge in KDIGO stage
  peak_trend_severe <- merge(peak_trend_severe,kdigo_grade,by="patient_id",all.x=TRUE)
  if(isTRUE(ckd_present)) {
    peak_trend_severe <- merge(peak_trend_severe,ckd_list,by="patient_id",all.x=TRUE)
    peak_trend_severe$ckd[is.na(peak_trend_severe$ckd)] <- 0
  }
  # Calculate mean and SD each for severe and non-severe groups
  peak_cr_summ <- peak_trend_severe %>% dplyr::group_by(severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  severe_label <- data.table::data.table(c(1,2,3,4),c("Non-severe, no AKI","Non-severe, AKI","Severe, no AKI","Severe, AKI"))
  colnames(severe_label) <- c("severe","severe_label")
  peak_cr_summ <- merge(peak_cr_summ,severe_label,by="severe",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # peak_cr_summ <- peak_cr_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the AKI with severity graphs...")
    peak_cr_summ <- peak_cr_summ[peak_cr_summ$n >= obfuscation_value,]
  }
  write.csv(peak_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_PeakCr_Severe_AKI.csv")),row.names=FALSE)
  # Plot the graphs
  peak_cr_timeplot <- ggplot2::ggplot(peak_cr_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,6) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_Severe_AKI.png")),plot=peak_cr_timeplot,width=12,height=9,units="cm")
  peak_cr_timeplot_raw <- ggplot2::ggplot(peak_cr_summ,ggplot2::aes(x=time_from_peak,y=mean_value,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_Severe_AKI_RawCr.png")),plot=peak_cr_timeplot_raw,width=12,height=9,units="cm")
  
  # Calculate mean and SD each for each KDIGO stage
  peak_cr_kdigo_summ <- peak_trend_severe %>% dplyr::group_by(aki_kdigo_grade,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  kdigo_label <- data.table::data.table(c(0,1,2,3),c("No AKI","KDIGO Stage 1","KDIGO Stage 2","KDIGO Stage 3"))
  colnames(kdigo_label) <- c("aki_kdigo_grade","kdigo_label")
  peak_cr_kdigo_summ <- merge(peak_cr_kdigo_summ,kdigo_label,by="aki_kdigo_grade",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # peak_cr_kdigo_summ <- peak_cr_kdigo_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the AKI with severity graphs...")
    peak_cr_kdigo_summ <- peak_cr_kdigo_summ[peak_cr_kdigo_summ$n >= obfuscation_value,]
  }
  write.csv(peak_cr_kdigo_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_PeakCr_KDIGOStage_AKI.csv")),row.names=FALSE)
  # Plot the graphs
  peak_cr_kdigo_timeplot <- ggplot2::ggplot(peak_cr_kdigo_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,6) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_KDIGOStage_AKI.png")),plot=peak_cr_kdigo_timeplot,width=12,height=9,units="cm")
  peak_cr_kdigo_timeplot_raw <- ggplot2::ggplot(peak_cr_kdigo_summ,ggplot2::aes(x=time_from_peak,y=mean_value,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_KDIGOStage_AKI_RawCr.png")),plot=peak_cr_kdigo_timeplot_raw,width=12,height=9,units="cm")
  
  severe_simplified_label <- data.table::data.table(c(0,1),c("Non-severe","Severe"))
  colnames(severe_simplified_label) <- c("severe","severe_label")
  
  try({
    # Stratify by KDIGO stage and COVID-19 severity
    peak_cr_kdigo_and_severe_summ <- peak_trend_severe %>% dplyr::group_by(aki_kdigo_grade,severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
    peak_cr_kdigo_and_severe_summ <- merge(peak_cr_kdigo_and_severe_summ,kdigo_label,by="aki_kdigo_grade",all.x=TRUE) %>% dplyr::mutate(severe = ifelse(severe > 2,1,0))
    peak_cr_kdigo_and_severe_summ <- merge(peak_cr_kdigo_and_severe_summ,severe_simplified_label,by="severe",all.x=TRUE)
    if(isTRUE(is_obfuscated)) {
      # peak_cr_kdigo_summ <- peak_cr_kdigo_summ %>% dplyr::filter(n >= obfuscation_value)
      cat("\nObfuscating the AKI with severity graphs...")
      peak_cr_kdigo_and_severe_summ <- peak_cr_kdigo_and_severe_summ[peak_cr_kdigo_and_severe_summ$n >= obfuscation_value,]
    }
    write.csv(peak_cr_kdigo_and_severe_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_PeakCr_KDIGOStage_CovidSevere_AKI.csv")),row.names=FALSE)
    # Plot the graphs
    peak_cr_kdigo_and_severe_timeplot <- ggplot2::ggplot(peak_cr_kdigo_and_severe_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,6) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~severe_label)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_KDIGOStage_CovidSevere_AKI.png")),plot=peak_cr_kdigo_and_severe_timeplot,width=12,height=8,units="cm")
    peak_cr_kdigo_and_severe_timeplot_raw <- ggplot2::ggplot(peak_cr_kdigo_and_severe_summ,ggplot2::aes(x=time_from_peak,y=mean_value,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~severe_label)
    ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_KDIGOStage_CovidSevere_AKI_RawCr.png")),plot=peak_cr_kdigo_and_severe_timeplot_raw,width=12,height=8,units="cm")
    
    if(isTRUE(ckd_present)) {
      peak_cr_ckd_and_severe_summ <- peak_trend_severe %>% dplyr::group_by(ckd,severe,time_from_peak) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
      ckd_label <- data.table::data.table(c(0,1),c("Non-CKD","CKD"))
      colnames(ckd_label) <- c("ckd","ckd_label")
      peak_cr_ckd_and_severe_summ <- merge(peak_cr_ckd_and_severe_summ,ckd_label,by="ckd",all.x=TRUE)
      peak_cr_ckd_and_severe_summ <- merge(peak_cr_ckd_and_severe_summ,severe_label,by="severe",all.x=TRUE)
      # peak_cr_ckd_and_severe_summ <- merge(peak_cr_ckd_and_severe_summ,ckd_label,by="ckd",all.x=TRUE) %>% dplyr::mutate(severe = ifelse(severe < 2,0,1))
      # severe_simplified_label <- data.table::data.table(c(0,1),c("Non-severe","Severe"))
      # colnames(severe_simplified_label) <- c("severe","severe_label")
      # peak_cr_ckd_and_severe_summ <- merge(peak_cr_ckd_and_severe_summ,severe_simplified_label,by="severe",all.x=TRUE)
      if(isTRUE(is_obfuscated)) {
        # peak_cr_ckd_summ <- peak_cr_ckd_summ %>% dplyr::filter(n >= obfuscation_value)
        cat("\nObfuscating the AKI with severity graphs...")
        peak_cr_ckd_and_severe_summ <- peak_cr_ckd_and_severe_summ[peak_cr_ckd_and_severe_summ$n >= obfuscation_value,]
      }
      write.csv(peak_cr_ckd_and_severe_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_PeakCr_CKD_CovidSevere_AKI.csv")),row.names=FALSE)
      # Plot the graphs
      peak_cr_ckd_and_severe_timeplot <- ggplot2::ggplot(peak_cr_ckd_and_severe_summ,ggplot2::aes(x=time_from_peak,y=mean_ratio,group=ckd_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,6) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~ckd_label)
      ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_CKD_CovidSevere_AKI.png")),plot=peak_cr_ckd_and_severe_timeplot,width=12,height=8,units="cm")
      peak_cr_ckd_and_severe_timeplot_raw <- ggplot2::ggplot(peak_cr_ckd_and_severe_summ,ggplot2::aes(x=time_from_peak,y=mean_value,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI Peak",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~ckd_label)
      ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_PeakCr_CKD_CovidSevere_AKI_RawCr.png")),plot=peak_cr_ckd_and_severe_timeplot_raw,width=12,height=8,units="cm")
      
    }
  })
  
  
  cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients (with severity) should have been generated.")
  
  #-------------------------------------------------------------------------------------
  # Plot from start of admission to 30 days post-peak AKI (if no AKI, then from peak Cr)
  adm_to_aki_cr <- labs_cr_all
  # adm_to_aki_cr$peak_cr_time[is.na(adm_to_aki_cr$peak_cr_time)] <- 0
  #adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(dplyr::between(days_since_admission,0,peak_cr_time+30)) %>% dplyr::ungroup()
  adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::filter(peak_cr_time == min(peak_cr_time,na.rm=TRUE)) %>% dplyr::distinct() %>% dplyr::ungroup()
  adm_to_aki_cr <- adm_to_aki_cr[order(adm_to_aki_cr$patient_id,adm_to_aki_cr$days_since_admission),]
  adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_365d,min_cr_retro_365d,na.rm=TRUE)) %>% dplyr::ungroup()
  adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
  adm_to_aki_cr <- merge(adm_to_aki_cr,kdigo_grade,by="patient_id",all.x = T)
  
  #adm_to_aki_cr <- adm_to_aki_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
  adm_to_aki_summ <- adm_to_aki_cr %>% dplyr::group_by(severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  adm_to_aki_summ <- merge(adm_to_aki_summ,severe_label,by="severe",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # adm_to_aki_summ <- adm_to_aki_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the admission to AKI graphs...")
    adm_to_aki_summ <- adm_to_aki_summ[adm_to_aki_summ$n >= obfuscation_value,]
  }
  write.csv(adm_to_aki_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromAdmToPeak+10D_Severe_AKI.csv")),row.names=FALSE)
  
  adm_to_aki_timeplot <- ggplot2::ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 30 & adm_to_aki_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromAdmToPeak+10D_Severe_AKI.png")),plot=adm_to_aki_timeplot,width=12,height=9,units="cm")
  adm_to_aki_timeplot_raw <- ggplot2::ggplot(adm_to_aki_summ[which(adm_to_aki_summ$days_since_admission <= 30 & adm_to_aki_summ$days_since_admission >= 0),],ggplot2::aes(x=days_since_admission,y=mean_value,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromAdmToPeak+10D_Severe_AKI_RawCr.png")),plot=adm_to_aki_timeplot_raw,width=12,height=9,units="cm")
  cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from first day of admission, should have been generated.")
  
  # Calculate mean and SD each for each KDIGO stage
  adm_to_aki_cr_kdigo_summ <- adm_to_aki_cr %>% dplyr::group_by(aki_kdigo_grade,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  adm_to_aki_cr_kdigo_summ <- merge(adm_to_aki_cr_kdigo_summ,kdigo_label,by="aki_kdigo_grade",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # adm_to_aki_cr_kdigo_summ <- adm_to_aki_cr_kdigo_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the AKI with severity graphs...")
    adm_to_aki_cr_kdigo_summ <- adm_to_aki_cr_kdigo_summ[adm_to_aki_cr_kdigo_summ$n >= obfuscation_value,]
  }
  write.csv(adm_to_aki_cr_kdigo_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrFromAdm_KDIGOStage_AKI.csv")),row.names=FALSE)
  # Plot the graphs
  adm_to_aki_cr_kdigo_timeplot <- ggplot2::ggplot(adm_to_aki_cr_kdigo_summ,ggplot2::aes(x=days_since_admission,y=mean_ratio,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days since admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,6) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromAdm_KDIGOStage_AKI.png")),plot=adm_to_aki_cr_kdigo_timeplot,width=12,height=9,units="cm")
  adm_to_aki_cr_kdigo_timeplot_raw <- ggplot2::ggplot(adm_to_aki_cr_kdigo_summ,ggplot2::aes(x=days_since_admission,y=mean_value,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days since admission",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromAdm_KDIGOStage_AKI_RawCr.png")),plot=adm_to_aki_cr_kdigo_timeplot_raw,width=12,height=9,units="cm")
  
  # Stratify by KDIGO stage and COVID-19 severity
  adm_to_aki_cr_kdigo_and_severe_summ <- adm_to_aki_cr %>% dplyr::group_by(aki_kdigo_grade,severe,days_since_admission) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  adm_to_aki_cr_kdigo_and_severe_summ <- merge(adm_to_aki_cr_kdigo_and_severe_summ,kdigo_label,by="aki_kdigo_grade",all.x=TRUE) %>% dplyr::mutate(severe = ifelse(severe <= 2,0,1))
  adm_to_aki_cr_kdigo_and_severe_summ <- merge(adm_to_aki_cr_kdigo_and_severe_summ,severe_simplified_label,by="severe",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # adm_to_aki_cr_kdigo_summ <- adm_to_aki_cr_kdigo_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the AKI with severity graphs...")
    adm_to_aki_cr_kdigo_and_severe_summ <- adm_to_aki_cr_kdigo_and_severe_summ[adm_to_aki_cr_kdigo_and_severe_summ$n >= obfuscation_value,]
  }
  write.csv(adm_to_aki_cr_kdigo_and_severe_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrFromAdm_KDIGOStage_CovidSevere_AKI.csv")),row.names=FALSE)
  # Plot the graphs
  adm_to_aki_cr_kdigo_and_severe_timeplot <- ggplot2::ggplot(adm_to_aki_cr_kdigo_and_severe_summ,ggplot2::aes(x=days_since_admission,y=mean_ratio,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days since admission",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,6) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~severe_label)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromAdm_KDIGOStage_CovidSevere_AKI.png")),plot=adm_to_aki_cr_kdigo_and_severe_timeplot,width=12,height=8,units="cm")
  adm_to_aki_cr_kdigo_and_severe_timeplot_raw <- ggplot2::ggplot(adm_to_aki_cr_kdigo_and_severe_summ,ggplot2::aes(x=days_since_admission,y=mean_value,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days since admission",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~severe_label)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromAdm_KDIGOStage_CovidSevere_AKI_RawCr.png")),plot=adm_to_aki_cr_kdigo_and_severe_timeplot_raw,width=12,height=8,units="cm")
  
  
  # ----------------------------------------
  # Plot from start of AKI to 30 days later 
  
  aki_30d_cr <- labs_cr_all
  # Uncomment the following line to restrict analysis to AKI patients only
  #aki_30d_cr <- aki_30d_cr[aki_30d_cr$severe %in% c(2,4,5),]
  aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(time_from_start = days_since_admission - aki_start) %>% dplyr::ungroup()
  aki_30d_cr <- aki_30d_cr[order(aki_30d_cr$patient_id,aki_30d_cr$days_since_admission),] %>% dplyr::distinct()
  aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(baseline_cr = min(min_cr_365d,min_cr_retro_365d)) %>% dplyr::ungroup()
  aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(ratio = value/baseline_cr) %>% dplyr::ungroup()
  aki_30d_cr <- merge(aki_30d_cr,kdigo_grade,by="patient_id",all.x = T)
  #aki_30d_cr <- aki_30d_cr %>% dplyr::group_by(patient_id) %>% dplyr::mutate(severe = ifelse((severe == 4 | severe == 5),4,severe))
  aki_30d_cr_summ <- aki_30d_cr %>% dplyr::group_by(severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  aki_30d_cr_summ <- merge(aki_30d_cr_summ,severe_label,by="severe",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    cat("\nObfuscating the start of AKI graphs...")
    # aki_30d_cr_summ <- aki_30d_cr_summ %>% dplyr::filter(n >= obfuscation_value)
    aki_30d_cr_summ <- aki_30d_cr_summ[aki_30d_cr_summ$n >= obfuscation_value,]
  }
  write.csv(aki_30d_cr_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.csv")),row.names=FALSE)
  aki_30d_cr_timeplot <- ggplot2::ggplot(aki_30d_cr_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,3.5) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI.png")),plot=aki_30d_cr_timeplot,width=12,height=9,units="cm")
  aki_30d_cr_timeplot_raw <- ggplot2::ggplot(aki_30d_cr_summ,ggplot2::aes(x=time_from_start,y=mean_value,group=severe_label))+ggplot2::geom_line(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(severe_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(severe_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Days from AKI start",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::scale_color_manual(values=c("Non-severe, AKI"="#bc3c29","Non-severe, no AKI"="#0072b5","Severe, AKI" = "#e18727","Severe, no AKI"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(file=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrFromStart_Severe_AKI_RawCr.png")),plot=aki_30d_cr_timeplot_raw,width=12,height=9,units="cm")
  
  # Calculate mean and SD each for each KDIGO stage
  aki_30d_cr_kdigo_summ <- aki_30d_cr %>% dplyr::group_by(aki_kdigo_grade,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  aki_30d_cr_kdigo_summ <- merge(aki_30d_cr_kdigo_summ,kdigo_label,by="aki_kdigo_grade",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # aki_30d_cr_kdigo_summ <- aki_30d_cr_kdigo_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the AKI with severity graphs...")
    aki_30d_cr_kdigo_summ <- aki_30d_cr_kdigo_summ[aki_30d_cr_kdigo_summ$n >= obfuscation_value,]
  }
  write.csv(aki_30d_cr_kdigo_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrfromAdm_KDIGOStage_AKI.csv")),row.names=FALSE)
  # Plot the graphs
  aki_30d_cr_kdigo_timeplot <- ggplot2::ggplot(aki_30d_cr_kdigo_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Time from AKI Start (days)",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,6) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromStart_KDIGOStage_AKI.png")),plot=aki_30d_cr_kdigo_timeplot,width=12,height=9,units="cm")
  aki_30d_cr_kdigo_timeplot_raw <- ggplot2::ggplot(aki_30d_cr_kdigo_summ,ggplot2::aes(x=time_from_start,y=mean_value,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Time from AKI Start (days)",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal()
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromStart_KDIGOStage_AKI_RawCr.png")),plot=aki_30d_cr_kdigo_timeplot_raw,width=12,height=9,units="cm")
  
  # Stratify by KDIGO stage and COVID-19 severity
  aki_30d_cr_kdigo_and_severe_summ <- aki_30d_cr %>% dplyr::group_by(aki_kdigo_grade,severe,time_from_start) %>% dplyr::summarise(mean_ratio = mean(ratio,na.rm=TRUE),sem_ratio = sd(ratio,na.rm=TRUE)/sqrt(dplyr::n()),mean_value = mean(value,na.rm=TRUE),sem_value = sd(value,na.rm=TRUE)/sqrt(dplyr::n()),n=dplyr::n()) %>% dplyr::ungroup()
  aki_30d_cr_kdigo_and_severe_summ <- merge(aki_30d_cr_kdigo_and_severe_summ,kdigo_label,by="aki_kdigo_grade",all.x=TRUE) %>% dplyr::mutate(severe = ifelse(severe <= 2,0,1))
  aki_30d_cr_kdigo_and_severe_summ <- merge(aki_30d_cr_kdigo_and_severe_summ,severe_simplified_label,by="severe",all.x=TRUE)
  if(isTRUE(is_obfuscated)) {
    # aki_30d_cr_kdigo_summ <- aki_30d_cr_kdigo_summ %>% dplyr::filter(n >= obfuscation_value)
    cat("\nObfuscating the AKI with severity graphs...")
    aki_30d_cr_kdigo_and_severe_summ <- aki_30d_cr_kdigo_and_severe_summ[aki_30d_cr_kdigo_and_severe_summ$n >= obfuscation_value,]
  }
  write.csv(aki_30d_cr_kdigo_and_severe_summ,file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_CrfromAdm_KDIGOStage_CovidSevere_AKI.csv")),row.names=FALSE)
  # Plot the graphs
  aki_30d_cr_kdigo_and_severe_timeplot <- ggplot2::ggplot(aki_30d_cr_kdigo_and_severe_summ,ggplot2::aes(x=time_from_start,y=mean_ratio,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_ratio-sem_ratio,ymax=mean_ratio+sem_ratio,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Time from AKI Start (days)",y = "Serum Cr/Baseline Cr", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(1,6) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~severe_label)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromStart_KDIGOStage_CovidSevere_AKI.png")),plot=aki_30d_cr_kdigo_and_severe_timeplot,width=12,height=8,units="cm")
  aki_30d_cr_kdigo_and_severe_timeplot_raw <- ggplot2::ggplot(aki_30d_cr_kdigo_and_severe_summ,ggplot2::aes(x=time_from_start,y=mean_value,group=kdigo_label))+ggplot2::geom_line(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_point(ggplot2::aes(color = factor(kdigo_label))) + ggplot2::geom_errorbar(ggplot2::aes(ymin=mean_value-sem_value,ymax=mean_value+sem_value,color = factor(kdigo_label)),position=ggplot2::position_dodge(0.05))+ ggplot2::theme(legend.position="right") + ggplot2::labs(x = "Time from AKI Start (days)",y = "Serum Cr (mg/dL)", color = "Severity") + ggplot2::xlim(-30,30) + ggplot2::ylim(0,4) + ggplot2::scale_color_manual(values=c("No AKI"="#bc3c29","KDIGO Stage 1"="#0072b5","KDIGO Stage 2" = "#e18727","KDIGO Stage 3"="#20854e")) + ggplot2::theme_minimal() + ggplot2::facet_wrap(~severe_label)
  ggplot2::ggsave(filename=file.path(getProjectOutputDirectory(), paste0(currSiteId,"_CrfromStart_KDIGOStage_CovidSevere_AKI_RawCr.png")),plot=aki_30d_cr_kdigo_and_severe_timeplot_raw,width=12,height=8,units="cm")
  
  cat("\nAt this point, if there are no errors, graphs and CSV files for normalised creatinine of AKI vs non-AKI patients, plotted from start of AKI/creatinine increase, should have been generated.")

  return(list("peak_aki_vs_non_aki" = peak_aki_vs_non_aki, "peak_trend_severe" = peak_trend_severe, "adm_to_aki_cr" = adm_to_aki_cr, "aki_from_start" = aki_30d_cr))
}
