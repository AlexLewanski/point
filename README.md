# **POInT**: **P**robabilistic **O**rigin **In**ference using **T**rees

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
The `point` repo holds the R package that I am currently calling POInT: Probabilistic Origin Inference using Trees. This package (and the method it implements) is early in development so check back later when it is more developed and tested!


## Installation
To install the current version of `POInT`, you can use the `devtools` package:

```r
install.packages("devtools")
devtools::install_github("AlexLewanski/point")
library("point")
```

## Instructions
Stay tuned for instructions on using the package when it is ready for use.


## Code organization
Functions in `R/` are (loosely) organized by topic:

- `em_funcs.R`: The core functions for estimating geographic origin via expectation-maximization (EM). This includes implementations of the EM algorithms and likelihood calculations used in these functions.
- `em_processing.R`: Processing the EM results including extracting the top EM run(s) and matching groups across runs. These functions are used within the EM algorithms (stored in `em_funcs.R`).
- `mvn_calcs.R`: Functionality for calculating the conditional normal distributions from the trees and coordinate information.
- `single_loc_est.R`: Some functions for estimating a single origin (when an EM algorithm isn't required).
- `tree_root_rates.R`: Estimating root locations and Brownian rate matrices.

Function documentation is maintained with the functions using roxygen comments.

