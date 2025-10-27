# coxphHorizon

**coxphHorizon** (R code files) implements a *horizon-based approximation strategy* for Cox proportional hazards (Cox PH) models.  
It allows users to share and validate fitted Cox models across different environments **without transferring the original data**.  
The package exports model coefficients and baseline survival probabilities at user-defined time horizons, enabling reproducible computation of survival and risk predictions from only the exported parameters.

---

## 🔍 Overview

Cox proportional hazards models are widely used for time-to-event analysis but often contain embedded baseline hazard functions that are large and tied to the original data source, making them difficult to transfer or share.  
**coxphHorizon** addresses this challenge by providing a *horizon-based approximation*, exporting only essential model parameters (coefficients and baseline survival estimates at selected time horizons).  
This approach preserves predictive accuracy while supporting collaborative validation and model deployment across datasets or institutions—without requiring access to the original data.

---

## ✨ Key Features

- Extract Cox PH model coefficients and baseline survival at specified time horizons  
- Compute survival probabilities for new data using exported parameters  
- Includes a **demonstration example** with simulated data showing the full workflow  
- Fully compatible with standard R survival modeling (`survival::coxph`)  
- Lightweight, reproducible, and data-agnostic (works with any time-to-event data)

---

## 🧩 Run demo example. 

Place all code files in the same directory :

```r
# simply source the simulate_demo.R file 
source("simulate_demo.R")
