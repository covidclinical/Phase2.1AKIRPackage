# 4CE Phase 2.1 - AKI Analysis
R code to run, validate, and submit the analysis for the AKI project.

To install this package in R:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE, force = TRUE)
```

To run the package, run the following command:

```
FourCePhase2.1AKI::runAnalysis(is_obfuscated = TRUE,obfuscation_value=3)
```
You can specify if obfuscation is applicable to your institution (is_obfuscated) and the cutoff for obfuscation (obfuscation_value).

To submit your data to the central repository:

```
FourCePhase2.1AKI::submitAnalysis()
```

Ensure that the data files are placed in your working directory before running the package!