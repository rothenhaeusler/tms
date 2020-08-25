# Targeted selection: parameter-specific model selection

This package provides a function for selecting among a set of estimators, if the goal is to minimize the mean-squared error with respect to a one-dimensional parameter of interest.

## How to install

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using  `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("rothenhaeusler/targeted-selection")`

## Usage

```R
gen_data <- function(s,n) {
  Tr <- rbinom(n,1,.5)
  X <- .5*Tr + rnorm(n)
  Y <- .5*X +  rnorm(n) + s^2*Tr
  
  return(as.data.frame(cbind(Tr,X,Y)))
}
df <- gen_data(.1,100)

surrogate_estimator <- function(df) coef(lm(Y~X,data=dff))[2]*coef(lm(X~Tr,data=dff))[2] 
difference_in_means <- function(df) coef(lm(Y~Tr,data=dff))[2] 

targeted(list(surrogate_estimator,difference_in_means),df)
```
