# 4CE Phase 2.1 - AKI Analysis
R code to run, validate, and submit the analysis for the AKI project.

## Installation
To install this package in R:
```
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE, force = TRUE)
```
## Issues with 4CE Docker Image v2.0

The latest Docker image as of this writing uses Microsoft R 4.0.2 which uses an older CRAN snapshot in July 2020. THis distribution installs survminer 0.4.7 and broom 0.7.0 by default.
When using survminer 0.4.7 and broom 0.7.0, the package may fail to generate forest plots with the following error:
```
Error in `[.data.frame`(cbind(allTermsDF, coef[inds, ]), , c("var", "level",  :
   undefined columns selected
```
A workaround is to force the Docker installation of Microsoft R to install a newer version of these packages with the following commands:
```
install.packages("survminer", repos = "http://mran.revolutionanalytics.com")
install.packages("broom", repos = "http://mran.revolutionanalytics.com")
```
These commands should be run **AFTER** installing the package, to prevent Microsoft R from overwriting these newer packages with older versions of survminer and broom.

## Running Analysis
Please ensure that all input files are stored in /4ceData/Input before running the package!
This code assumes that it is running within the Docker environment and all file paths are made with reference to the Docker environment.
To run the package, run the following command:

```
setwd("/4ceData")
FourCePhase2.1AKI::runAnalysis(is_obfuscated=TRUE,obfuscation_value=3,factor_cutoff = 5,restrict_models = FALSE)
```
## Arguments:
- is_obfuscated - Specify if obfuscation is applicable to your institution. Default = TRUE
- obfuscation_value - If obfuscation is implemented in your institution, specify the cutoff for obfuscation here. Default = 3
- factor_cutoff - For Cox proportional hazard modeling, specify the minimum number of events for each factor level. Default = 5
- restrict_models - (ADVANCED) This is a new feature implemented to address problems that each site may have with Cox proportional hazard model generation despite the semi-automated methods of generating the models. If set to TRUE, you must also include a text file CustomModelVariables.txt which contain the variable names, separated by spaces and all in a single line, to restrict modelling to. 

## Specifying Custom Variables for Cox Models
(ADVANCED) This is a new feature implemented to address problems that each site may have with Cox proportional hazard model generation despite the semi-automated methods of generating the models. 
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
