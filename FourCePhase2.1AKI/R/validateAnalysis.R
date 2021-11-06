
#' Validates the results of the analytic workflow for the AKI project
#'
#' @keywords 4CE
#' @export

validateAnalysis <- function(docker = TRUE, siteid_nodocker = "") {
  if(isTRUE(docker)) {
    currSiteId = toupper(FourCePhase2.1Data::getSiteId())
  } else {
    currSiteId = siteid_nodocker
  }
  
  output_folder = getProjectOutputDirectory()
  demog_obf <- read.csv(file.path(output_folder,paste0(currSiteId,"_TableOne_obfuscated.csv")))
  cr_graph_aki_vs_nonaki <- read.csv(file.path(output_folder,paste0(currSiteId,"_PeakCr_AKI_vs_NonAKI.csv")))
  cr_graph_aki_covidsevere <- read.csv(file.path(output_folder,paste0(currSiteId,"_PeakCr_Severe_AKI.csv")))
  cr_graph_kdigo <- read.csv(file.path(output_folder,paste0(currSiteId,"_PeakCr_KDIGOStage_AKI.csv")))
  cr_graph_kdigo_covidsevere <- read.csv(file.path(output_folder,paste0(currSiteId,"_PeakCr_KDIGOStage_CovidSevere_AKI.csv")))
  
  cat("First checking for consistency within the sCr graphs...\n")
  cr_graph_aki_vs_nonaki <- cr_graph_aki_vs_nonaki %>% dplyr::filter(time_from_peak == 0)
  cr_graph_aki_covidsevere <- cr_graph_aki_covidsevere %>% dplyr::filter(time_from_peak == 0)
  cr_graph_kdigo <- cr_graph_kdigo %>% dplyr::filter(time_from_peak == 0)
  cr_graph_kdigo_covidsevere <- cr_graph_kdigo_covidsevere %>% dplyr::filter(time_from_peak == 0)
  
  cr_graph_aki_total <- as.numeric(cr_graph_aki_vs_nonaki$n[cr_graph_aki_vs_nonaki$aki == 1])
  cr_graph_nonaki_total <- as.numeric(cr_graph_aki_vs_nonaki$n[cr_graph_aki_vs_nonaki$aki == 0])
  
  cr_graph_nonsevere_with_nonaki <- as.numeric(cr_graph_aki_covidsevere$n[cr_graph_aki_covidsevere$severe == 1])
  cr_graph_nonsevere_with_aki <- as.numeric(cr_graph_aki_covidsevere$n[cr_graph_aki_covidsevere$severe == 2])
  cr_graph_severe_with_nonaki <- as.numeric(cr_graph_aki_covidsevere$n[cr_graph_aki_covidsevere$severe == 3])
  cr_graph_severe_with_aki <- as.numeric(cr_graph_aki_covidsevere$n[cr_graph_aki_covidsevere$severe == 4])
  
  cr_graph_nonsevere_total <- cr_graph_nonsevere_with_nonaki + cr_graph_nonsevere_with_aki
  cr_graph_severe_total <- cr_graph_severe_with_nonaki + cr_graph_severe_with_aki
  
  cat("\nComparing AKI vs Non-AKI graph with AKI + COVID-19 severity graph...\n")
  if(cr_graph_aki_total == (cr_graph_nonsevere_with_aki + cr_graph_severe_with_aki)) {
    cat("\nAKI counts match\n")
  } else {
    cat("\nInconsistent AKI counts:\nAKI vs Non-AKI graph: ",cr_graph_aki_total,"\nAKI + COVID-19 severity graph: ",(cr_graph_nonsevere_with_aki + cr_graph_severe_with_aki),"\n")
  }
  if(cr_graph_nonaki_total == (cr_graph_nonsevere_with_nonaki + cr_graph_severe_with_nonaki)) {
    cat("\nNon-AKI counts match\n")
  } else {
    cat("\nInconsistent Non-AKI counts:\nAKI vs Non-AKI graph: ",cr_graph_nonaki_total,"\nAKI + COVID-19 severity graph: ",(cr_graph_nonsevere_with_nonaki + cr_graph_severe_with_nonaki),"\n")
  }
  
  cr_graph_kdigo_akitotal <- as.numeric(sum(cr_graph_kdigo$n[cr_graph_kdigo$aki_kdigo_grade != 0]))
  cr_graph_kdigo_0 <- as.numeric(cr_graph_kdigo$n[cr_graph_kdigo$aki_kdigo_grade == 0])
  cr_graph_kdigo_1 <- as.numeric(cr_graph_kdigo$n[cr_graph_kdigo$aki_kdigo_grade == 1])
  cr_graph_kdigo_2 <- as.numeric(cr_graph_kdigo$n[cr_graph_kdigo$aki_kdigo_grade == 2])
  cr_graph_kdigo_3 <- as.numeric(cr_graph_kdigo$n[cr_graph_kdigo$aki_kdigo_grade == 3])
  
  cat("\nComparing AKI vs Non-AKI graph with KDIGO Stage graph...\n")
  if(cr_graph_aki_total == cr_graph_kdigo_akitotal) {
    cat("\nAKI counts match\n")
  } else {
    cat("\nInconsistent AKI counts:\nAKI vs Non-AKI graph: ",cr_graph_aki_total,"\nKDIGO graph: ",cr_graph_kdigo_akitotal,"\n")
  }
  if(cr_graph_nonaki_total == cr_graph_kdigo_0) {
    cat("\nNon-AKI counts match\n")
  } else {
    cat("\nInconsistent Non-AKI counts:\nAKI vs Non-AKI graph: ",cr_graph_nonaki_total,"\nKDIGO graph: ",cr_graph_kdigo_0,"\n")
  }
  
  cr_graph_kdigo_covid_akitotal <- as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$aki_kdigo_grade != 0]))
  cr_graph_kdigo_covid_nonakitotal <- as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$aki_kdigo_grade == 0]))
  
  cr_graph_kdigo_covid_0 <- as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$aki_kdigo_grade == 0]))
  cr_graph_kdigo_covid_1 <- as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$aki_kdigo_grade == 1]))
  cr_graph_kdigo_covid_2 <- as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$aki_kdigo_grade == 2]))
  cr_graph_kdigo_covid_3 <- as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$aki_kdigo_grade == 3]))
  
  cr_graph_kdigo_covid_severetotal <-  as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$severe == 1]))
  cr_graph_kdigo_covid_nonseveretotal <-  as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$severe == 0]))
  cr_graph_kdigo_covid_nonsevere_nonaki <-  as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$severe == 0 & cr_graph_kdigo_covidsevere$aki_kdigo_grade == 0]))
  cr_graph_kdigo_covid_severe_nonaki <-  as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$severe == 1 & cr_graph_kdigo_covidsevere$aki_kdigo_grade == 0]))
  cr_graph_kdigo_covid_nonsevere_aki <-  as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$severe == 0 & cr_graph_kdigo_covidsevere$aki_kdigo_grade > 0]))
  cr_graph_kdigo_covid_severe_aki <-  as.numeric(sum(cr_graph_kdigo_covidsevere$n[cr_graph_kdigo_covidsevere$severe == 1 & cr_graph_kdigo_covidsevere$aki_kdigo_grade > 0]))
  
  cat("\nComparing AKI vs Non-AKI graph with KDIGO Stage + COVID-19 severity graph...\n")
  if(cr_graph_aki_total == cr_graph_kdigo_covid_akitotal) {
    cat("\nAKI counts match\n")
  } else {
    cat("\nInconsistent AKI counts:\nAKI vs Non-AKI graph: ",cr_graph_aki_total,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_akitotal,"\n")
  }
  if(cr_graph_nonaki_total == cr_graph_kdigo_covid_nonakitotal) {
    cat("\nNon-AKI counts match\n")
  } else {
    cat("\nInconsistent Non-AKI counts:\nAKI vs Non-AKI graph: ",cr_graph_nonaki_total,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_nonakitotal,"\n")
  }
  
  cat("\nComparing AKI + COVID-19 severity graph with KDIGO Stage + COVID-19 severity graph...\n")
  if(cr_graph_nonsevere_with_aki == cr_graph_kdigo_covid_nonsevere_aki) {
    cat("\nNon-severe, AKI counts match\n")
  } else {
    cat("\nInconsistent Non-severe, AKI counts:\nAKI + COVID-19 severity graph: ",cr_graph_nonsevere_with_aki,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_nonsevere_aki,"\n")
  }
  if(cr_graph_nonsevere_with_nonaki == cr_graph_kdigo_covid_nonsevere_nonaki) {
    cat("\nNon-severe, Non-AKI counts match\n")
  } else {
    cat("\nInconsistent Non-severe, Non-AKI counts:\nAKI + COVID-19 severity graph: ",cr_graph_nonsevere_with_nonaki,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_nonsevere_nonaki,"\n")
  }
  if(cr_graph_severe_with_aki == cr_graph_kdigo_covid_severe_aki) {
    cat("\nSevere, AKI counts match\n")
  } else {
    cat("\nInconsistent Severe, AKI counts:\nAKI + COVID-19 severity graph: ",cr_graph_severe_with_aki,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_severe_aki,"\n")
  }
  if(cr_graph_severe_with_nonaki == cr_graph_kdigo_covid_severe_nonaki) {
    cat("\nSevere, Non-AKI counts match\n")
  } else {
    cat("\nInconsistent Severe, Non-AKI counts:\nAKI + COVID-19 severity graph: ",cr_graph_severe_with_nonaki,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_severe_nonaki,"\n")
  }
  
  
  cat("\nComparing KDIGO Stage graph with KDIGO Stage + COVID-19 severity graph...\n")
  if(cr_graph_kdigo_0 == cr_graph_kdigo_covid_0) {
    cat("\nKDIGO Stage 0 match\n")
  } else {
    cat("\nInconsistent KDIGO Stage 0 counts:\nKDIGO Stage graph: ",cr_graph_kdigo_0,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_0,"\n")
  }
  
  if(cr_graph_kdigo_1 == cr_graph_kdigo_covid_1) {
    cat("\nKDIGO Stage 1 match\n")
  } else {
    cat("\nInconsistent KDIGO Stage 1 counts:\nKDIGO Stage graph: ",cr_graph_kdigo_1,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_1,"\n")
  }
  
  if(cr_graph_kdigo_2 == cr_graph_kdigo_covid_2) {
    cat("\nKDIGO Stage 2 match\n")
  } else {
    cat("\nInconsistent KDIGO Stage 2 counts:\nKDIGO Stage graph: ",cr_graph_kdigo_2,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_2,"\n")
  }
  
  if(cr_graph_kdigo_3 == cr_graph_kdigo_covid_3) {
    cat("\nKDIGO Stage 3 match\n")
  } else {
    cat("\nInconsistent KDIGO Stage 3 counts:\nKDIGO Stage graph: ",cr_graph_kdigo_3,"\nKDIGO Stage + COVID-19 severity graph: ",cr_graph_kdigo_covid_3,"\n")
  }

  cat("\n===========================================\n")
  cat("\nNow checking for consistency within the Kaplan-Meier curves.\n\n")
  
  km_death_aki_only_covidsevere <- read.csv(file.path(output_folder,paste0(currSiteId,"_TimeToEvent_Death_AKIOnly_Severe_Plot.csv")))
  km_death_aki_only_kdigo <- read.csv(file.path(output_folder,paste0(currSiteId,"_TimeToEvent_Death_AKIOnly_KDIGO_Plot.csv")))
  km_recover_aki_only_covidsevere <- read.csv(file.path(output_folder,paste0(currSiteId,"_TimeToEvent_Recover_Severe_Plot.csv")))
  km_recover_aki_only_kdigo <- read.csv(file.path(output_folder,paste0(currSiteId,"_TimeToEvent_Recover_KDIGO_Plot.csv")))
  km_death_aki_vs_nonaki <- read.csv(file.path(output_folder,paste0(currSiteId,"_TimeToEvent_Death_AKIvsNonAKI_Plot.csv")))
  km_death_all_kdigo <- read.csv(file.path(output_folder,paste0(currSiteId,"_TimeToEvent_Death_KDIGO_Plot.csv")))
  
  km_death_stats <- km_death_aki_vs_nonaki %>% dplyr::group_by(is_aki) %>% dplyr::summarise(deceased = sum(n.event)) %>% dplyr::ungroup()
  km_death_stats_nonakitotal <- as.numeric(km_death_stats$deceased[km_death_stats$is_aki == 0])
  km_death_stats_akitotal <- as.numeric(km_death_stats$deceased[km_death_stats$is_aki == 1])
  
  km_recover_aki_only_covidsevere <- km_recover_aki_only_covidsevere %>% dplyr::group_by(severe) %>% dplyr::filter(time == min(time)) %>% dplyr::select(n.risk) %>% dplyr::ungroup()
  km_recover_aki_only_kdigo <- km_recover_aki_only_kdigo %>% dplyr::group_by(aki_kdigo_final) %>% dplyr::filter(time == min(time)) %>% dplyr::select(n.risk) %>% dplyr::ungroup()
  km_death_aki_only_covidsevere <- km_death_aki_only_covidsevere %>% dplyr::group_by(severe) %>% dplyr::filter(time == min(time)) %>% dplyr::select(n.risk) %>% dplyr::ungroup()
  km_death_aki_only_kdigo <- km_death_aki_only_kdigo %>% dplyr::group_by(aki_kdigo_final) %>% dplyr::filter(time == min(time)) %>% dplyr::select(n.risk) %>% dplyr::ungroup()
  km_death_aki_vs_nonaki <- km_death_aki_vs_nonaki %>% dplyr::group_by(is_aki) %>% dplyr::filter(time == min(time)) %>% dplyr::select(n.risk) %>% dplyr::ungroup()
  km_death_all_kdigo <- km_death_all_kdigo %>% dplyr::group_by(aki_kdigo_final) %>% dplyr::filter(time == min(time)) %>% dplyr::select(n.risk) %>% dplyr::ungroup()
  
  km_death_aki_only_severetotal <- as.numeric(km_death_aki_only_covidsevere$n.risk[km_death_aki_only_covidsevere$severe == 1])
  km_death_aki_only_nonseveretotal <- as.numeric(km_death_aki_only_covidsevere$n.risk[km_death_aki_only_covidsevere$severe == 0])
  km_death_aki_only_kdigo_1 <- as.numeric(km_death_aki_only_kdigo$n.risk[km_death_aki_only_kdigo$aki_kdigo_final == 1])
  km_death_aki_only_kdigo_2 <- as.numeric(km_death_aki_only_kdigo$n.risk[km_death_aki_only_kdigo$aki_kdigo_final == 2])
  km_death_aki_only_kdigo_3 <- as.numeric(km_death_aki_only_kdigo$n.risk[km_death_aki_only_kdigo$aki_kdigo_final == 3])
  km_death_aki_only_kdigo_total <- km_death_aki_only_kdigo_1 + km_death_aki_only_kdigo_2 + km_death_aki_only_kdigo_3
  
  km_death_aki_vs_nonaki_akitotal <- as.numeric(km_death_aki_vs_nonaki$n.risk[km_death_aki_vs_nonaki$is_aki == 1])
  km_death_aki_vs_nonaki_nonakitotal <- as.numeric(km_death_aki_vs_nonaki$n.risk[km_death_aki_vs_nonaki$is_aki == 0])
  
  km_death_all_kdigo_0 <- as.numeric(km_death_all_kdigo$n.risk[km_death_all_kdigo$aki_kdigo_final == 0])
  km_death_all_kdigo_1 <- as.numeric(km_death_all_kdigo$n.risk[km_death_all_kdigo$aki_kdigo_final == 1])
  km_death_all_kdigo_2 <- as.numeric(km_death_all_kdigo$n.risk[km_death_all_kdigo$aki_kdigo_final == 2])
  km_death_all_kdigo_3 <- as.numeric(km_death_all_kdigo$n.risk[km_death_all_kdigo$aki_kdigo_final == 3])
  km_death_all_akitotal <- km_death_all_kdigo_1 + km_death_all_kdigo_2 + km_death_all_kdigo_3
  
  km_recover_aki_only_severetotal <- as.numeric(km_recover_aki_only_covidsevere$n.risk[km_recover_aki_only_covidsevere$severe == 1])
  km_recover_aki_only_nonseveretotal <- as.numeric(km_recover_aki_only_covidsevere$n.risk[km_recover_aki_only_covidsevere$severe == 0])
  km_recover_aki_only_kdigo_1 <- as.numeric(km_recover_aki_only_kdigo$n.risk[km_recover_aki_only_kdigo$aki_kdigo_final == 1])
  km_recover_aki_only_kdigo_2 <- as.numeric(km_recover_aki_only_kdigo$n.risk[km_recover_aki_only_kdigo$aki_kdigo_final == 2])
  km_recover_aki_only_kdigo_3 <- as.numeric(km_recover_aki_only_kdigo$n.risk[km_recover_aki_only_kdigo$aki_kdigo_final == 3])
  km_recover_aki_only_kdigo_total <- km_recover_aki_only_kdigo_1 + km_recover_aki_only_kdigo_2 + km_recover_aki_only_kdigo_3
  
  cat("Comparing KDIGO patients from Mortality (AKI only) and Mortality (All) graphs...\n")
  if(km_death_aki_only_kdigo_1 == km_death_all_kdigo_1) {
    cat("KDIGO 1 match\n")
  } else {
    cat("KDIGO 1 inconsistency: AKI only = ", km_death_aki_only_kdigo_1,", All = ",km_death_all_kdigo_1,"\n")
  }
  if(km_death_aki_only_kdigo_2 == km_death_all_kdigo_2) {
    cat("KDIGO 2 match\n")
  } else {
    cat("KDIGO 2 inconsistency: AKI only = ", km_death_aki_only_kdigo_2,", All = ",km_death_all_kdigo_2,"\n")
  }
  if(km_death_aki_only_kdigo_3 == km_death_all_kdigo_3) {
    cat("KDIGO 3 match\n")
  } else {
    cat("KDIGO 3 inconsistency: AKI only = ", km_death_aki_only_kdigo_3,", All = ",km_death_all_kdigo_3,"\n")
  }
  if(km_death_aki_only_kdigo_total == km_death_all_akitotal) {
    cat("KDIGO total match\n")
  } else {
    cat("KDIGO total inconsistency: AKI only = ", km_death_aki_only_kdigo_total,", All = ",km_death_all_akitotal,"\n")
  }
  
  cat("Comparing KDIGO patients from Recovery (AKI only) and Mortality (All) graphs...\n")
  if(km_recover_aki_only_kdigo_1 == km_death_all_kdigo_1) {
    cat("KDIGO 1 match\n")
  } else {
    cat("KDIGO 1 inconsistency: Recovery = ", km_recover_aki_only_kdigo_1,", Death = ",km_death_all_kdigo_1,"\n")
  }
  if(km_recover_aki_only_kdigo_2 == km_death_all_kdigo_2) {
    cat("KDIGO 2 match\n")
  } else {
    cat("KDIGO 2 inconsistency: Recovery = ", km_recover_aki_only_kdigo_2,", Death = ",km_death_all_kdigo_2,"\n")
  }
  if(km_recover_aki_only_kdigo_3 == km_death_all_kdigo_3) {
    cat("KDIGO 3 match\n")
  } else {
    cat("KDIGO 3 inconsistency: Recovery = ", km_recover_aki_only_kdigo_3,", Death = ",km_death_all_kdigo_3,"\n")
  }
  if(km_recover_aki_only_kdigo_total == km_death_all_akitotal) {
    cat("KDIGO total match\n")
  } else {
    cat("KDIGO total inconsistency: Recovery = ", km_recover_aki_only_kdigo_total,", Death = ",km_death_all_akitotal,"\n")
  }
  
  cat("Comparing COVID-19 severity patients between Recovery (AKI only) and Mortality (AKI only)...\n")
  if(km_recover_aki_only_severetotal == km_death_aki_only_severetotal) {
    cat("Severe match\n")
  } else {
    cat("Severe inconsistency: Recovery = ", km_recover_aki_only_severetotal,", Mortality = ",km_death_aki_only_severetotal,"\n")
  }
  if(km_recover_aki_only_nonseveretotal == km_death_aki_only_nonseveretotal) {
    cat("Non-severe match\n")
  } else {
    cat("Non-severe inconsistency: Recovery = ", km_recover_aki_only_nonseveretotal,", Mortality = ",km_death_aki_only_nonseveretotal,"\n")
  }
  
  cat("Comparing non-AKI totals between Mortality (AKI vs Non-AKI) and Mortality (KDIGO 0-3) graphs...\n")
  if(km_death_aki_vs_nonaki_nonakitotal == km_death_all_kdigo_0) {
    cat("Non-AKI match\n")
  } else {
    cat("Non-AKI inconsistency: AKI vs Non-AKI = ", km_death_aki_vs_nonaki_nonakitotal,", KDIGO 0-3 = ",km_death_all_kdigo_0,"\n")
  }
  
  cat("\n==============================\n")
  cat("Now cross-comparing across sCr graphs and Kaplan-Meier curves.\n")
  cat("For simplification, we will be doing the following comparisons:\n")
  cat("1) Severe/Non-Severe counts (AKI patients only) - AKI+COVID-19 severity graph, Mortality (AKI only)\n")
  cat("2) KDIGO staging (AKI patients only) - KDIGO Stage graph, Mortality (All, KDIGO staging)\n")
  cat("\n-----------\n")
  cat("1) Severe/Non-Severe counts (AKI patients only) - AKI+COVID-19 severity graph, Mortality (AKI only)\n")
  if(cr_graph_severe_with_aki == km_death_aki_only_severetotal) {
    cat("Severe counts match\n")
  } else {
    cat("Severe counts mismatch: sCr graph = ", cr_graph_severe_with_aki,", KM curve = ", km_death_aki_only_severetotal,"\n")
  }
  if(cr_graph_nonsevere_with_aki == km_death_aki_only_nonseveretotal) {
    cat("Non-severe counts match\n")
  } else {
    cat("Non-severe counts mismatch: sCr graph = ", cr_graph_nonsevere_with_aki,", KM curve = ", km_death_aki_only_nonseveretotal,"\n")
  }
  cat("\n-----------\n")
  cat("2) KDIGO staging (AKI patients only) - KDIGO Stage graph, Mortality (All, KDIGO staging)\n")
  if(cr_graph_kdigo_1 == km_death_all_kdigo_1) {
    cat("KDIGO 1 match\n")
  } else {
    cat("KDIGO 1 inconsistency: sCr graph = ", cr_graph_kdigo_1,", KM curve = ",km_death_all_kdigo_1,"\n")
  }
  if(cr_graph_kdigo_2 == km_death_all_kdigo_2) {
    cat("KDIGO 2 match\n")
  } else {
    cat("KDIGO 2 inconsistency: sCr graph = ", cr_graph_kdigo_2,", KM curve = ",km_death_all_kdigo_2,"\n")
  }
  if(cr_graph_kdigo_3 == km_death_all_kdigo_3) {
    cat("KDIGO 3 match\n")
  } else {
    cat("KDIGO 3 inconsistency: sCr graph = ", cr_graph_kdigo_3,", KM curve = ",km_death_all_kdigo_3,"\n")
  }
  if(cr_graph_kdigo_akitotal == km_death_all_akitotal) {
    cat("KDIGO total match\n")
  } else {
    cat("KDIGO total inconsistency: sCr graph = ", cr_graph_kdigo_akitotal,", KM curve = ",km_death_all_akitotal,"\n")
  }
  
  demog_obf <- demog_obf[demog_obf$category %in% c("n","severe_Non-severe","severe_Severe","aki_kdigo_grade_No AKI","aki_kdigo_grade_Stage 1","aki_kdigo_grade_Stage 2","aki_kdigo_grade_Stage 3","deceased_Alive","deceased_Deceased"),c("category","No_AKI","AKI")]
  cr_demog_comparison_tmp <- data.frame(demog_obf$category,c(cr_graph_nonaki_total,cr_graph_nonsevere_with_nonaki,cr_graph_severe_with_nonaki,(cr_graph_nonaki_total - km_death_stats_nonakitotal),km_death_stats_nonakitotal,cr_graph_kdigo_0,0,0,0),c(cr_graph_aki_total,cr_graph_nonsevere_with_aki,cr_graph_severe_with_aki,(cr_graph_aki_total - km_death_stats_akitotal),km_death_stats_akitotal,0,cr_graph_kdigo_1,cr_graph_kdigo_2,cr_graph_kdigo_3))
  colnames(cr_demog_comparison_tmp) <- c("category","No_AKI_graphs","AKI_graphs")
  demog_obf <- merge(demog_obf,cr_demog_comparison_tmp,by="category",all.x=T)
  demog_obf <- demog_obf %>% dplyr::mutate(diff_nonaki = No_AKI - No_AKI_graphs,diff_aki = AKI - AKI_graphs) %>% dplyr::filter(diff_aki != 0)
  if(length(demog_obf$category) > 0) {
    cat("\nMismatches found in counts between sCr graphs/KM curves and demographics tables.\nRefer to the specific counts below for details.\nIf the errors arise from NAs due to obfuscation, please ignore this message.\nUse this information to aid in debugging.\n")
    print(data.table::as.data.table(demog_obf))
  }
}

