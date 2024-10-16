# VDPO

<!-- badges: start -->

[![Build_Status](https://github.com/Pavel-Hernadez-Amaro/VDPO/actions/workflows/build.yml/badge.svg)](https://github.com/Pavel-Hernadez-Amaro/VDPO/actions/workflows/build.yml)

<!-- badges: end -->

The **VDPO** package provides tools for working with and analyzing functional data of varying lengths. This variation in length can occur in two different scenarios: **Variable Domain Data and Partially Observed Data.** This refers to cases where the domain over which the data is observed changes between observations or when the functional data are not fully observed over the entire domain of interest, respectively. For instance, in growth curve analysis, each individual might have measurements starting and ending at different ages, leading to varying observation ranges. Similarly, in environmental studies, different locations might have data collected over distinct time periods, creating domains of different lengths.

## Related Papers

-   Pavel Hernandez-Amaro, Maria Durban, M. Carmen Aguilera-Morillo, Cristobal Esteban Gonzalez, Inmaculada Arostegui. "Modelling physical activity profiles in COPD patients: a fully functional approach to variable domain functional regression models." doi: [10.48550/arXiv.2401.05839](https://doi.org/10.48550/arXiv.2401.05839)

## Installation

To install the package from GitHub, the [remotes](https://cran.r-project.org/package=remotes) package is required:

``` r
# install.packages("remotes")
remotes::install_github("Pavel-Hernadez-Amaro/VDPO")
```

## Web

The web page of the package can be accessed from [this link](https://pavel-hernadez-amaro.github.io/VDPO/).
