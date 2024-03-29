# Flexible age-specific spatio-temporal Bayesian P-splines models

This repository contains the R code to fit in INLA the age-specific spatio-temporal Bayesian P-splines models described in *"Flexible Bayesian P-splines for smoothing age-specific spatio-temporal mortality patterns"* [(Goicoa et al., 2019)](https://doi.org/10.1177/0962280217726802).

## Table of contents

-   [Data](#Data)
-   [R code](#R-code)
-   [References](#References)

# Data {#data}

Simulated female breast cancer mortality data (ICD-10 code 50) in Spanish provinces during the period 1985-2010 by age groups.

-   [**BreastCancer_ESP.txt**](https://github.com/spatialstatisticsupna/Flexible_Psplines_article/blob/master/data/BreastCancer_ESP.txt)

    This .txt file contains a data set with the following variables:

    -   ***Age.group***: numeric vector of age-group identifiers.
    -   ***Age.label***: character vector of age-groups labels.
    -   ***Province***: numeric vector of geographic identifiers (Spanish provinces).
    -   ***Year***: numeric vector of year's identifiers.
    -   ***Obs***: observed number of deaths.
    -   ***Pop***: Population at risk.

    Note that the observed number of deaths has been modified to preserve data confidentiality.

-   [**Carto_ESP.Rdata**](https://github.com/spatialstatisticsupna/Flexible_Psplines_article/blob/master/data/Carto_ESP.Rdata): `sf` object containing the spatial polygons of the Spanish provinces.

-   [**Esp_prov_nb.graph**](https://github.com/spatialstatisticsupna/Flexible_Psplines_article/blob/master/data/Esp_prov_nb.inla): An inla.graph object with the spatial neighbourhood structure of the 50 provinces of Spain.

# R code {#r-code}

R code to fit with INLA (<https://www.r-inla.org/>) the age-specific spatio-temporal Bayesian P-splines models described in Goicoa et al. (2019). All the R files are written by the authors of the paper.

-   [**Psplines_INLA.R**](https://github.com/spatialstatisticsupna/Flexible_Psplines_article/blob/master/R/Psplines_INLA.R)

    Main script including the required functions to fit in INLA the different age-specific spatio-temporal Bayesian P-splines described in the paper. Note that the obtained results will differ from those shown in the original paper, since simulated data are being used here to preserve data confidentiality.

-   [**Figures_and_Tables.R**](https://github.com/spatialstatisticsupna/Flexible_Psplines_article/blob/master/R/Figures_and_Tables.R)

    This R script contains the necessary functions to reproduce some of the figures and tables shown in Goicoa et al. (2019). The final model fitted with INLA (Model9-SI with posterior patterns) can be downloaded from <https://emi-sstcdapp.unavarra.es/Flexible_Psplines_article/Model9_SI.Rdata>.

# Acknowledgements

This work has been supported by the Spanish Ministry of Economy and Competitiveness (project MTM2014-51992-R), and by the Health Department of the Navarre Government (Project 113, Res.2186/2014).

# References {#references}

[Goicoa, T., Adin, A., Etxeberria, J., Militino, A.F., and Ugarte, M.D. (2019). Flexible Bayesian P-splines for smoothing age-specific spatio-temporal mortality patterns. *Statistical Methods in Medical Research*, **28(2)**, 384-403.](https://doi.org/10.1177/0962280217726802)
