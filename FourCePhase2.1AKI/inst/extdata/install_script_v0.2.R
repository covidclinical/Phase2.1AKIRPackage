# Install script for FourCePhase2.1AKI package version 0.2.0.0
# =============================================================
# This revised package has a dependency on tidycmprsk to generate cleaned versions of
# tables from cmprsk for competing risk analysis
# However, as this package was released in 2022, this is not available in the MRAN repository (2020)
# Hence, will need to force install this from the CRAN repository
#
# For sites which use an offline system that is not physically connected to the internet, ensure that
# - Minimum of R 4.0
# - You have created an offline repository of the relevant dependencies. If the previous versions had
#   worked correctly, this means that you only need to install tidycmprsk, cmprsk and update the 
#   other dependencies as needed
# Follow the tutorial on https://cran.r-project.org/web/packages/miniCRAN/vignettes/miniCRAN-introduction.html
# for a guide on offline repository creation and installation from an offline repository

message("First updating the Phase 2.1 Data R package...")
devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

message("Now we are going to force installation of the tidycmprsk package from CRAN.")
message("Please note: this will update several other packages which may cause packages to break.")
message("(Case in point: portions of the AKI package had to be rewritten due to changes in rlang 1.0 causing errors when parsing object names in data.table objects within the package environment)")
message("We will also install broom from CRAN as the version in MRAN will cause errors with the survminer functions.")
utils::install.packages(c("tidycmprsk","cmprsk","broom"),repos = "https://cloud.r-project.org/")

message("Finally, install the AKI package itself.")
devtools::install_github("https://github.com/covidclinical/Phase2.1AKIRPackage", subdir="FourCePhase2.1AKI", upgrade=FALSE, force = TRUE)
