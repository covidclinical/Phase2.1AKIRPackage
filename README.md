# 4CE Phase 2.1 - AKI Analysis
R code to run, validate, and submit the analysis for the AKI project.

To install this package in R:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE)
```

Ensure that the data files are placed in your working directory before running the package!

## Note for Diagnoses Codes
Due to the variations in the way ICD codes may be stored, we have adhered to the following conventions:

- icd_code will store the broad category of the diagnosis
- Full ICD diagnosis codes will have the major category and subcategory separated by a period

For example, the ICD-10 diagnosis "Chronic kidney disease, stage 1" with full code N18.1 will have icd_code = "N18" and full_code = "N18.1".
Please ensure that your data tables store diagnoses codes as their broad categories. (e.g. N18.1 will be stored as "N18").
In future versions, this convention may be changed depending on the granularity required for analysis.
