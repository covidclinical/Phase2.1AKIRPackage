# 4CE Phase 2.1 - AKI Analysis
R code to run, validate, and submit the analysis for the AKI project.

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

## Known Issues with R 4.1
1) If your site is running a non-Docker based environment running on R 4.1.\*, there is a known issue where the Kaplan-Meier plots produced may be blank
- We are currently in the midst of attempting to fix this
- As a workaround, try to run the package in the older R 4.0.2 environment instead