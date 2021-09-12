# 4CE Phase 2.1 - AKI Analysis
R code to run, validate, and submit the analysis for the AKI project.

## Prerequisites
Please ensure that the following is done prior to analysis:
- Laboratory data dating to **365 days prior to admission** (as compared to the default -60 days) have been extracted.
- You are running the R package in either Docker v2.1.0 and above, or in an environment using R 4.0.\* (see below for issues with R 4.1.\*)

We have provided modified SQL scripts for Phase 1.1/2.1 (MSSQL and Oracle) in inst/extdata/phaseX.1_modified_sql, which implement the longer durations of lab value extraction, for reference. These scripts should also extract procedure codes for renal replacement therapy (RRT) and kidney transplants across all time.

## Installation
To install this package in R:
```
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE, force = TRUE)
```
Please ensure you are using 4CE Docker image v2.1.0 to avoid issues with generating forest plots.

## Running Analysis
Please ensure that all input files are stored in /4ceData/Input before running the package!
This code assumes that it is running within the Docker environment and all file paths are made with reference to the Docker environment.
To run the package, run the following command:

```
setwd("/4ceData")
FourCePhase2.1AKI::runAnalysis()
```
## Arguments:
By default, the analysis should not require any additional arguments.
These arguments are only needed if your site is experiencing issues with generating results:
- is_obfuscated - Allows you to toggle obfuscation on and off. Default = TRUE
- factor_cutoff - For Cox proportional hazard modeling, specify the minimum number of events for each factor level. Default = 5
- ckd_cutoff - Sets the baseline serum creatinine cutoff (in mg/dL) beyond which patients are excluded from analysis. Default = 2.25
- restrict_models - **(ADVANCED)** If set to TRUE, you must also include a text file CustomModelVariables.txt which contain the variable names, separated by spaces and all in a single line, to restrict modelling to. Default = FALSE
- docker - **(ADVANCED)** Indicates if running in a Docker environment. Set this to FALSE if you are running on native Windows/Linux/macOS. Default = TRUE
- input - **(ADVANCED)** Specifies the directory where files are stores. **MANDATORY** if docker set to FALSE. Default = '/4ceData/Input'
- siteid_nodocker - **(ADVANCED)** Specifies the siteid to be used if docker is set to FALSE. Default = ''
- skip_qc - **(ADVANCED - DO NOT CHANGE UNLESS NECESSARY)** If docker is set to FALSE, allows the QC step to be bypassed. Mainly for debugging purposes - **DO NOT CHANGE**. Default = FALSE
- use_rrt_surrogate - **(ADVANCED)** If enabled, the algorithm for detecting possible RRT episodes will also be used to exclude patients who may have had RRT based on serum creatinine fluctuations. Default = FALSE
- print_rrt_surrogate - **(ADVANCED - DO NOT CHANGE UNLESS NECESSARY)** If enabled, the algorithm for detecting possible RRT episodes will also print out **patient-level data** indicating the dates of suspected RRT for each patient. The files produced are documented below. **DO NOT CHANGE UNLESS INTENDING TO DO MANUAL CHART REVIEWS!** Default = FALSE

## Detecting Possible RRT Episodes
Starting from version 0.1.5, the AKI R package implements a function where days of possible RRT procedures are detected based on serum creatinine fluctuations.
A suspected RRT episode comprises the following:
- Peak serum creatinine >= 3mg/dL, AND
- Decrease in seurm creatinine by >= 25% in strictly 24h or less


This was implemented as some sites do not have access to procedure codes (even into Phase 1.2/2.2) and thus a surrogate measure besides diagnoses/procedure codes was needed.
The relevant arguments to the runAnalysis() function are as follows:
- use_rrt_surrogate - **(ADVANCED)** If enabled, the algorithm for detecting possible RRT episodes will also be used to exclude patients who may have had RRT based on serum creatinine fluctuations. Default = FALSE
- print_rrt_surrogate - **(ADVANCED - DO NOT CHANGE UNLESS NECESSARY)** If enabled, the algorithm for detecting possible RRT episodes will also print out **patient-level data** indicating the dates of suspected RRT for each patient. The files produced are documented below. **DO NOT CHANGE UNLESS INTENDING TO DO MANUAL CHART REVIEWS!** Default = FALSE


The files produced when **print_rrt_surrogate=TRUE** are as follows:
- DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv - suspected RRT episodes across the WHOLE cohort of patients
- DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv - suspected RRT episodes after excluding patients with baseline serum creatinine greater than or equal to ckd_cutoff

These files contain **patient-level data** and are intended for the purposes of site-level manual chart review only. **These files are not used in analyses, and should be removed prior to submission!**


The files can be found in ~/Phase2.1AKIR (default output folder).

**Before submitting using submitAnalysis(), please ensure that these two files have been removed.**

## Specifying Custom Variables for Cox Models
**(ADVANCED)** This is a feature implemented to address problems that each site may have with Cox proportional hazard model generation despite the semi-automated methods of generating the models. 
If using this feature, you must include a text file CustomModelVariables.txt which contain the variable names, separated by spaces and all in a single line, to restrict modelling to. 

An example of CustomModelVariables.txt is as follows:
```
age_group sex severe aki_kdigo_final ckd htn hld ihd cld bronchiectasis copd rheum vte malig_solid malig_lymph COAGA COAGB covid_rx 
```
This example is also included in the data-raw folder of this package for reference.

Ensure to specify restrict_models = TRUE if you are using this method of customization, and inform the rest of the AKI gorup about any customization you are using with justification.

## Submitting Data
To submit your data to the central repository:

```
FourCePhase2.1AKI::submitAnalysis()
```
If you want to submit your output files in this manner, please ensure you have been granted read/write access to the data repositories. You may contact Byorn Tan (@byorntan) for details.

**Important note**:
If print_rrt_surrogate is set to TRUE, please ensure that these files have been copied to a safe location and **deleted prior to submission**!
- DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv
- DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv

## Known Issues with R 4.1
1) If your site is running a non-Docker based environment running on R 4.1.\*, there is a known issue where the Kaplan-Meier plots produced may be blank
- We are currently in the midst of attempting to fix this
- As a workaround, try to run the package in the older R 4.0.2 environment instead

## Work in Progress
1) Debugging v0.1.5.0 for any issues after extending baseline serum creatinine definitions to prior 365 days only
2) Transitioning package to use Phase 2.2-formatted data when specified
3) Modifying Phase 2.2 SQL scripts to incorporate missing procedure codes for RRT and kidney transplant (e.g. peritoneal dialysis codes not included in most cases)
4) Implementing an age-based serum creatinine cutoff based on the CKD-EPI equation and using an assumed eGFR < 30mL/kg/1.73m2 cutoff (using lower limits of age groups as reference, given the obfuscation of age present in Phase 1.1/2.1)
