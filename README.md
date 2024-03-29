# 4CE Phase 2.1 - AKI Analysis
Last updated: August 13, 2022 (v0.2.0.8)

R code to run, validate, and submit the analysis for the AKI project.

## NEW CHANGES (as of July 27, 2022)
This package has a new dependency on the package `tidycmprsk` which is only available in the CRAN repository but **not** in the MRAN repository referenced in the Docker image. Other new dependencies include the packages `cmprsk`,`broom`, `ellipsis`, `xfun` and `rlang`, whose latest versions (required for the package to run) are **NOT** in MRAN!

Install `tidycmprsk`, `cmprsk`, `broom`, `ellipsis`, `xfun` and `rlang` using the following commands in R prior to installing the package:
```
utils::install.packages(c("tidycmprsk","cmprsk","broom","ellipsis","xfun","rlang"),repos = "https://cloud.r-project.org/")
```
You **MUST** restart your R session before the AKI package can run!

This package has also been tested in the following R versions:
- R 4.2.1 (vanilla)
- R 4.1.3 (vanilla)
- Microsoft R Open 4.0.2 (Docker version)

Only versions R 4.1.* and above may experience issues with generating blank plots (see below).

## Prerequisites
Please ensure that the following is done prior to analysis:
- Only patients admitted from Jan 1, 2020 to **Sep 10, 2020** are extracted **(no earlier than Sep 10, 2020)**
- All lab data/diagnoses/procedures/medications are extracted from Jan 1, **2019** (-365 days) to Sep 10, **2021** (+365 days)
- Procedure codes for RRT (hemodialysis + peritoneal dialysis) and kidney transplants are included as well, across **all** dates (i.e. even those before Jan 1, 2019 and after Sep 10, 2021)
- You are running the R package in either the 4CE Docker environment v2.1.0 and above, or in an environment using R 4.0.\* (see below for issues with R 4.1.\*)

We have provided modified SQL scripts for Phase 1.1/2.1 (MSSQL and Oracle) in `inst/extdata/phaseX.1_modified_sql`, which implement the longer durations of lab value extraction. 

These scripts should also extract procedure codes (ICD9/ICD-10-CM) for renal replacement therapy (RRT) and kidney transplants across all time.

Please modify these scripts if your site uses a different procedure code scheme as per Phase 1.2/2.2 guidelines, e.g. CPT4/CCAM/SNOMED-CT

## Installation
To install this package in R:

### RStudio
```
utils::install.packages(c("tidycmprsk","cmprsk","broom","ellipsis","xfun","rlang"),repos = "https://cloud.r-project.org/")
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE, force = TRUE)
rstudioapi::restartSession()
```
### R in console
If you are running R from a console, you need to first install the dependencies from CRAN followed by the package:
```
utils::install.packages(c("tidycmprsk","cmprsk","broom","ellipsis","xfun","rlang"),repos = "https://cloud.r-project.org/")
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE, force = TRUE)
```

You **MUST** restart your R session prior to running the package to ensure that the updated dependencies are loaded!

An install script, `install_script_v0.2.R`, has been included under the `inst/extdata` folder.

## Running Analysis
Please ensure that all input files are stored in `/4ceData/Input` before running the package!
This code assumes that it is running within the Docker environment and all file paths are made with reference to the Docker environment.
To run the package, run the following command:

```
setwd("/4ceData")
FourCePhase2.1AKI::runAnalysis()
```

To validate results:
```
FourCePhase2.1AKI::validateAnalysis()
```

## Arguments:
By default, the analysis should not require any additional arguments.

These arguments are only needed if your site is experiencing issues with generating results.

Note that the arguments apply for both `runAnalysis()` and `validateAnalysis()` functions.

The `validateAnalysis()` function only has two arguments - `docker` and `siteid_nodocker`.
| Argument              | Default Value    | Description                                                                                                                                                                                                                                                                                                                 |
|-----------------------|------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `is_obfuscated`       | TRUE             | Allows you to toggle obfuscation on and off                                                                                                                                                                                                                                                                                 |
| `factor_cutoff`       | 5                | For Cox proportional hazard modeling, specify the minimum number of events for each  factor level.                                                                                                                                                                                                                          |
| `ckd_cutoff`          | 2.25             | Sets the baseline serum creatinine cutoff (in mg/dL) beyond which patients are excluded  from analysis.                                                                                                                                                                                                                     |
| `restrict_models`     | FALSE            | **(ADVANCED)** If set to TRUE, you must also include a text file CustomModelVariables.txt  which contain the variable names, separated by spaces and all in a single line, to  restrict modelling to.                                                                                                                       |
| `docker`              | TRUE             | **(ADVANCED)** Indicates if running in a Docker environment. Set this to `FALSE` if you  are running on native Windows/Linux/macOS.                                                                                                                                                                                         |
| `input`               | '/4ceData/Input' | **(ADVANCED)** Specifies the directory where files are stores. **MANDATORY** if docker  set to `FALSE`.                                                                                                                                                                                                                     |
| `siteid_nodocker`     | ''               | **(ADVANCED)** Specifies the siteid to be used if docker is set to `FALSE`.                                                                                                                                                                                                                                                 |
| `skip_qc`             | FALSE            | **(ADVANCED - DO NOT CHANGE UNLESS NECESSARY)** If docker is set to FALSE, allows the QC  step to be bypassed. Mainly for debugging purposes - **DO NOT CHANGE**.                                                                                                                                                           |
| `offline`             | FALSE            | Disables the version check for the `FourCePhase2.1Data` package on GitHub. Useful if  running on offline systems.                                                                                                                                                                                                           |
| `use_rrt_surrogate`   | TRUE            | **(ADVANCED)** If enabled, the algorithm for detecting possible RRT episodes will also be  used to exclude patients who may have had RRT based on serum creatinine fluctuations.                                                                                                                                            |
| `print_rrt_surrogate` | FALSE            | **(ADVANCED - DO NOT CHANGE UNLESS NECESSARY)** If enabled, the algorithm for detecting  possible RRT episodes will also print out **patient-level data** indicating the dates of  suspected RRT for each patient. The files produced are documented below. **DO NOT CHANGE  UNLESS INTENDING TO DO MANUAL CHART REVIEWS!** |
| `debug_on`            | FALSE            | Redirects output in `stderr()` (i.e. `message()`, `warning()`, `error()`) to the error log (`~/Phase2.1AKIR/<siteid>_error.log`). Due to limitations in how R handles console output,  warning/error messages will not show up on the console if this option is enabled.                                                    |

## Detecting Possible RRT Episodes
Starting from version 0.1.5, the AKI R package implements a function where days of possible RRT procedures are detected based on serum creatinine fluctuations.
A suspected RRT episode comprises the following:
- Peak serum creatinine >= 3mg/dL, AND
- Decrease in seurm creatinine by >= 25% in strictly 24h or less


This was implemented as some sites do not have access to procedure codes (even into Phase 1.2/2.2) and thus a surrogate measure besides diagnoses/procedure codes was needed.
The relevant arguments to the `runAnalysis()` function are as follows:
- `use_rrt_surrogate` - **(ADVANCED)** If enabled, the algorithm for detecting possible RRT episodes will also be used to exclude patients who may have had RRT based on serum creatinine fluctuations. Default = FALSE
- `print_rrt_surrogate` - **(ADVANCED - DO NOT CHANGE UNLESS NECESSARY)** If enabled, the algorithm for detecting possible RRT episodes will also print out **patient-level data** indicating the dates of suspected RRT for each patient. The files produced are documented below. **DO NOT CHANGE UNLESS INTENDING TO DO MANUAL CHART REVIEWS!** Default = FALSE


The files produced when **`print_rrt_surrogate=TRUE`** are as follows:
- `DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv` - suspected RRT episodes across the WHOLE cohort of patients
- `DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv` - suspected RRT episodes after excluding patients with baseline serum creatinine greater than or equal to ckd_cutoff

These files contain **patient-level data** and are intended for the purposes of site-level manual chart review only. **These files are not used in analyses, and should be removed prior to submission!**


The files can be found in `~/Phase2.1AKIR` (default output folder).

**Before submitting using `submitAnalysis()`, please ensure that these two files have been removed.**

## Specifying Custom Variables for Cox Models
**(ADVANCED)** This is a feature implemented to address problems that each site may have with Cox proportional hazard model generation despite the semi-automated methods of generating the models. 
If using this feature, you must include a text file `CustomModelVariables.txt` which contain the variable names, separated by spaces and all in a single line, to restrict modelling to. 

An example of `CustomModelVariables.txt` is as follows:
```
age_group sex severe aki_kdigo_final ckd htn hld ihd cld bronchiectasis copd rheum vte malig_solid malig_lymph COAGA COAGB covid_rx 
```
This example is also included in the data-raw folder of this package for reference.

Ensure to specify `restrict_models = TRUE` if you are using this method of customization, and inform the rest of the AKI gorup about any customization you are using with justification.

## Submitting Data
To submit your data to the central repository:

```
FourCePhase2.1AKI::submitAnalysis()
```
If you want to submit your output files in this manner, please ensure you have been granted read/write access to the data repositories. You may contact Byorn Tan (@byorntan) for details.

**Important note**:
If `print_rrt_surrogate` is set to TRUE, please ensure that these files have been copied to a safe location and **deleted prior to submission**!
- `DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection.csv`
- `DO_NOT_UPLOAD_PATIENT_LEVEL_DATA_RRT_Surrogate_Detection_ExcludeCrCutoff.csv`

## Known Issues with R 4.1 and above
1) If your site is running a non-Docker based environment running on R 4.1.\*, there is a known issue where the Kaplan-Meier plots produced may be blank
   - We are currently in the midst of attempting to fix this
   - As a workaround, try to run the package in the older R 4.0.2 environment instead
   
2)  If your site is running a non-Docker based environment running R 4.2.\* and above, you may encounter the following error when the package is attempting to parse medications:
```
Error in `tidyr::spread()`:
! Each row of output must be identified by a unique combination of keys.
```
We do not have a feasible workaround for this issue. This issue has not surfaced in our testing on R 4.0.2 and R 4.1.3 so far. Please raise this issue to us if you encounter this issue on versions older than R 4.2 and we will try to rectify this as soon as possible.

## Work in Progress
1) [x] Debugging v0.1.5.0 for any issues after extending baseline serum creatinine definitions to prior 365 days only
2) [ ] Transitioning package to use Phase 2.2-formatted data when specified
3) [ ] Modifying Phase 2.2 SQL scripts to incorporate missing procedure codes for RRT and kidney transplant (e.g. peritoneal dialysis codes not included in most cases)
4) [ ] Implementing an age-based serum creatinine cutoff based on the CKD-EPI equation and using an assumed eGFR < 30mL/kg/1.73m2 cutoff (using lower limits of age groups as reference, given the obfuscation of age present in Phase 1.1/2.1)
