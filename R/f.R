#' Targeted model selection
#' @param f1 An asymptotically unbiased estimator of the parameter of interest
#' @param f2 A potentially biased estimator of the parameter of interest
#' @param df The data
#' @param R Number of bootstrap draws
#' @export
targeted <- function(f1,f2,df,R=200) {
  bstrap <- function(){
    indices <- sample.int(nrow(df),size=nrow(df),replace = TRUE)
    return(simplify2array(lapply(X=list(f1,f2),FUN = function(f) f(df[indices,]))))
  }
  bstrap_reps <- replicate(R,bstrap())
  f1_df <- f1(df)
  f2_df <- f2(df)
  criterion <- function(s) {
    max(( (1-s)*f1_df + s*f2_df  - f1_df )^2 - var((1-s)*bstrap_reps[1,]+s*bstrap_reps[2,] - bstrap_reps[1,]),0) + var((1-s)*bstrap_reps[1,]+s*bstrap_reps[2,])
  }
  criterion_vec <- sapply(seq(0,1,length.out=10),criterion)
  criterion_subset <- function() {
    subset <- sample.int(ncol(bstrap_reps),size = 1,replace=TRUE)
    criterion_b <- function(s){
    max(( (1-s)*mean(bstrap_reps[1,subset]) + s*mean(bstrap_reps[2,subset])  - mean(bstrap_reps[1,subset]) )^2 - var((1-s)*bstrap_reps[1,]+s*bstrap_reps[2,] - bstrap_reps[1,]),0) + var((1-s)*bstrap_reps[1,]+s*bstrap_reps[2,])}
    whichmin <- which.min(sapply(seq(0,1,length.out=10),criterion_b))
    estimator_targeted <- (1-seq(0,1,length.out=10)[whichmin])*mean(bstrap_reps[1,subset]) + (seq(0,1,length.out=10)[whichmin])*mean(bstrap_reps[2,subset])
    return(estimator_targeted)
  }
  bstraps <- replicate(criterion_subset(),n=R)
  
  whichmin <- which.min(criterion_vec)
  estimator_targeted <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  ret_vec <- c(estimator_targeted,quantile(bstraps,.025),quantile(bstraps,.975))
  names(ret_vec) <- c("final estimate", "95% bootstrap CI, lower bound", "95% bootstrap CI, upper bound" )
  
  return(ret_vec)
}
