# Targeted selection: parameter-specific model selection

This package provides a function for selecting among a set of estimators, if the goal is to minimize the mean-squared error with respect to a one-dimensional parameter of interest.

## How to install

1. The [devtools](https://github.com/hadley/devtools) package has to be installed. You can install it using  `install.packages("devtools")`.
2. The latest development version can then be installied using `devtools::install_github("rothenhaeusler/tms")`.

## Usage

```R
n <- 100
Tr <- rbinom(n,1,.5)
X <- .5*Tr + rnorm(n)
Y <- .5*X +  rnorm(n) + .01*Tr
df <-  as.data.frame(cbind(Tr,X,Y))

surrogate_estimator <- function(df) coef(lm(Y~X,data=df))[2]*coef(lm(X~Tr,data=df))[2]
difference_in_means <- function(df) coef(lm(Y~Tr,data=df))[2] 

# The first argument should be an asymptotically unbiased estimator for the parameter of interest; it serves as a benchmark.
tms(difference_in_means,surrogate_estimator,df)
```
