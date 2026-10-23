# *POInT*: *P*robabilistic *O*rigin *In*ference using *T*rees

## Table of Contents
- [Author information](#author-information)
- [Overview](#overview)
- [Installation](#installation)
- [Instructions](#instructions)
- [Code organization](#code-organization)

## Author information

Name (ORCID link)| Contact | Affiliation(s)
:-----:|:-----|:-----
[Alexander Lewanski](https://orcid.org/0000-0001-5843-0837) | [allewanski AT gmail DOT com](mailto:allewanski@gmail.com)|Michigan State University, Kellogg Biological Station, University of Michigan


## Overview
The `point` repo holds the R package that I am currently calling POInT: **P**robabilistic **O**rigin **In**ference using **T**rees. This package (and the method it implements) is early in development so check back later when it is more developed and tested!


## Installation
To install the current version of the `point` package, you can use the [`devtools`](https://devtools.r-lib.org) package:

```r
install.packages("devtools")
devtools::install_github("AlexLewanski/point")
library("point")
```

## Instructions
Stay tuned for instructions on using the package when it is ready for use.


## Code organization
Functions in [`R/`](R/) are loosely organized by topic:

- [`em_helpers.R`](R/em_helpers.R): Helper functions for the EM algorithms such as likelihood calculations.
- [`em_location_initialization.R`](R/em_location_initialization.R): Initialization of origin locations for EM estimation.
- [`em_main.R`](R/em_main.R): Core functions for estimating geographic origins with expectation-maximization (EM), including the EM algorithms and wrappers to run and process the EM algorithms multiple times.
- [`em_processing.R`](R/em_processing.R): Processing EM results including selecting top runs and matching groups across runs.
- [`mvn_calcs.R`](R/mvn_calcs.R): Functions for calculating conditional normal distributions from trees and coordinate data.
- [`single_loc_est.R`](R/single_loc_est.R): Supporting functions for estimating a single origin without an EM algorithm.
- [`tree_root_rates.R`](R/tree_root_rates.R): Functions for estimating root locations and Brownian rate matrices.


Function documentation is maintained with the functions using roxygen comments.

