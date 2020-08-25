targeted <- function(f1,f2,df,R=200) {
  bstrap <- function(){
    indices <- sample.int(nrow(df),size=nrow(df),replace = TRUE)
    return(simplify2array(lapply(X=list(f1,f2),FUN = function(f) f(df[indices,]))))
  }
  bstrap_reps <- replicate(R,bstrap())
  criterion <- function(s) {
    max(( (1-s)*f1 + s*f2  - f1 )^2 - var((1-s)*bstrap_reps[1,]+s*bstrap_reps[2,] - bstrap_reps[1,]),0) + var((1-s)*bstrap_reps[1,]+s*bstrap_reps[2,])
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
  estimator_targeted <- (1-seq(0,1,length.out=10)[whichmin])*f1 + (seq(0,1,length.out=10)[whichmin])*f2
  
  ret_vec <- c(estimator_targeted,quantile(bstraps,.025),quantile(bstraps,.975))
  names(ret_vec) <- c("final estimate", "95% bootstrap CI, lower bound", "95% bootstrap CI, upper bound" )
  
  return(ret_vec)
}
