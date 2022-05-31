setwd("~/Documents/git/targeted_cvx/code_ejsrev/")

iterations <- 100
mc_cores <- 32
mc_cores <- 4
# merge code for d>1 and d=1

# Functions for performing targeted model selection
tms <- function(f1,f2,df,R=200){
  bstrap <- function(){
    indices <- sample.int(nrow(df),size=nrow(df),replace = TRUE)
    return(simplify2array(lapply(X=list(f1,f2),FUN = function(f) f(df[indices,]))))
  }
  bstrap_reps <- replicate(R,bstrap())
  f1_df <- f1(df)
  f2_df <- f2(df)
  criterion <- function(s) {
    max( sum(( (1-s)*f1_df + s*f2_df  - f1_df )^2) - sum(diag(cov((1-s)*t(bstrap_reps[,1,])+s*t(bstrap_reps[,2,]) - t(bstrap_reps[,1,])))),0) + sum(diag(cov((1-s)*t(bstrap_reps[,1,])+s*t(bstrap_reps[,2,]))))
  }
  criterion_vec <- sapply(seq(0,1,length.out=10),criterion)
  
  criterion_subset <- function() {
    subset <- sample.int(length(bstrap_reps[1,1,]),size = 1,replace=TRUE)
    criterion_b <- function(s){
      max( sum(( (1-s)*bstrap_reps[,1,subset] + s*bstrap_reps[,2,subset]  - bstrap_reps[,1,subset] )^2) - sum(diag(cov((1-s)*t(bstrap_reps[,1,])+s*t(bstrap_reps[,2,]) - t(bstrap_reps[,1,])))),0) + sum(diag(cov((1-s)*t(bstrap_reps[,1,])+s*t(bstrap_reps[,2,]))))
      }
    whichmin <- which.min(sapply(seq(0,1,length.out=10),criterion_b))
    estimator_targeted <- (1-seq(0,1,length.out=10)[whichmin])*bstrap_reps[,1,subset] + (seq(0,1,length.out=10)[whichmin])*bstrap_reps[,2,subset]
    return(estimator_targeted)
  }
  bstraps <- replicate(criterion_subset(),n=R)
  
  whichmin <- which.min(criterion_vec)
  estimator_targeted <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  ret_vec <- c(estimator_targeted,apply(bstraps,1,quantile,.025),apply(bstraps,1,quantile,.975))
  indices <- 1:length(estimator_targeted)
  names(ret_vec) <- c(paste("final estimate",indices),paste("95% bootstrap CI, lower bound",indices), paste("95% bootstrap CI, upper bound",indices) )
  
  return(ret_vec)
}

tms_1D <- function(f1,f2,df,R=50) {
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

# Usage of model selection function
# n <- 100
# Tr <- rbinom(n,1,.5)
# X <- .5*Tr + rnorm(n)
# Y <- .5*X +  rnorm(n) + .01*Tr
# df <-  as.data.frame(cbind(Tr,X,Y))
# surrogate_estimator <- function(df) coef(lm(Y~X,data=df))[2]*coef(lm(X~Tr,data=df))[2]
# difference_in_means <- function(df) coef(lm(Y~Tr,data=df))[2] 
# f1 <- function(df) { test <- surrogate_estimator(df); return(c(test,test*2,test*3))}
# f2 <- function(df) { test <- difference_in_means(df); return(c(test,test*2,test*3))}
# R <- 200
# tms(f1,f2,df)

# Functions for performing targeted cross-validation
crossval <- function(list_of_funcs,df,folds=10){
  criterion_vec <- rep(0,folds)
  n <- nrow(df)
  for (j in 0:(folds-1)) {
    minn <- ceiling(quantile(1:n,j/folds))
    maxx  <- floor(quantile(1:n,(j+1)/folds))
    n <- nrow(df)
    samples_in <- setdiff(1:n,minn:maxx)
    samples_out <- minn:maxx
    f1_in <- list_of_funcs[[1]](df[samples_in,])
    f2_in <- list_of_funcs[[2]](df[samples_in,])
    test <- list_of_funcs[[1]](df[samples_out,])
    criterion_vec <- cbind(criterion_vec,colMeans(sapply(seq(0,1,length.out=10), function(s){ ( (1-s)*f1_in + s*f2_in  - test)^2})/folds))
    
  }
  criterion_vec <- rowMeans(criterion_vec,na.rm=TRUE)
  
  f1_df <- list_of_funcs[[1]](df)
  f2_df <- list_of_funcs[[2]](df)
  
  whichmin <- which.min(criterion_vec)
  estimator_crossval <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  return(estimator_crossval)
}

cui <- function(list_of_funcs,df,folds=10){
  criterion_vec <- rep(0,folds)
  n <- nrow(df)
  for (j in 0:(folds-1)) {
    minn <- ceiling(quantile(1:n,j/folds))
    maxx  <- floor(quantile(1:n,(j+1)/folds))
    n <- nrow(df)
    samples_in <- setdiff(1:n,minn:maxx)
    samples_out <- minn:maxx
    f1_in <- list_of_funcs[[1]](df[samples_in,])
    f2_in <- list_of_funcs[[2]](df[samples_in,])
    test <- list_of_funcs[[1]](df[samples_out,])
    criterion_vec <- cbind(criterion_vec,colMeans(sapply(seq(0,1,length.out=10), function(s){ pmax( ((1-s)*f1_in + s*f2_in  - f1_in)^2,( (1-s)*f1_in + s*f2_in  - f2_in)^2)})/folds))
    
  }
  criterion_vec <- rowMeans(criterion_vec,na.rm=TRUE)
  
  f1_df <- list_of_funcs[[1]](df)
  f2_df <- list_of_funcs[[2]](df)
  
  whichmin <- which.min(criterion_vec)
  estimator_crossval <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  return(estimator_crossval)
}

crossval_1D <- function(list_of_funcs,df,folds=5){
  criterion_vec <- rep(0,folds)
  n <- nrow(df)
  # To stabilize cross-validation, we randomly split into train and test 10 times.
  for (j in 0:(folds-1)) {
    minn <- ceiling(quantile(1:n,j/folds))
    maxx  <- floor(quantile(1:n,(j+1)/folds))
    n <- nrow(df)
    samples_in <- setdiff(1:n,minn:maxx)
    samples_out <- minn:maxx
    f1_in <- list_of_funcs[[1]](df[samples_in,])
    f2_in <- list_of_funcs[[2]](df[samples_in,])
    test <- list_of_funcs[[1]](df[samples_out,])
    criterion_vec <- cbind(criterion_vec,sapply(seq(0,1,length.out=10), function(s){ ( (1-s)*f1_in + s*f2_in  - test)^2})/folds)
    
  }
  criterion_vec <- rowMeans(criterion_vec,na.rm=TRUE)
  
  f1_df <- list_of_funcs[[1]](df)
  f2_df <- list_of_funcs[[2]](df)
  
  whichmin <- which.min(criterion_vec)
  estimator_crossval <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  return(estimator_crossval)
}



cui_1D <- function(list_of_funcs,df,folds=5){
  criterion_vec <- rep(0,folds)
  n <- nrow(df)
  # To stabilize cross-validation, we randomly split into train and test 10 times.
  for (j in 0:(folds-1)) {
    minn <- ceiling(quantile(1:n,j/folds))
    maxx  <- floor(quantile(1:n,(j+1)/folds))
    n <- nrow(df)
    samples_in <- setdiff(1:n,minn:maxx)
    samples_out <- minn:maxx
    f1_in <- list_of_funcs[[1]](df[samples_in,])
    f2_in <- list_of_funcs[[2]](df[samples_in,])
    test <- list_of_funcs[[1]](df[samples_out,])
    criterion_vec <- cbind(criterion_vec,sapply(seq(0,1,length.out=10), function(s){ max( ((1-s)*f1_in + s*f2_in  - f1_in)^2,((1-s)*f1_in + s*f2_in  - f2_in)^2) })/folds)
    
  }
  criterion_vec <- rowMeans(criterion_vec,na.rm=TRUE)
  
  f1_df <- list_of_funcs[[1]](df)
  f2_df <- list_of_funcs[[2]](df)
  
  whichmin <- which.min(criterion_vec)
  estimator_crossval <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  return(estimator_crossval)
}


#### Section 4.1 Observational studies: heterogeneous treatment effects

gen_data <- function(s,n) {
  E <- sample(1:3,n,replace=TRUE)
  X <- runif(n) <= .5
  Tr <- ifelse(X,runif(n) <= .5, runif(n) <= .1)
  Y <- X + 3*s^2*Tr*X+  rnorm(n) + Tr*(E==1)*.1 + Tr*(E==2)*.2 + Tr*(E==3)*(-.1)
  
  return(as.data.frame(cbind(E,Tr,X,Y)))
}


ate <- function(dff) {
  Trhat <- fitted.values(lm(Tr ~ X,data=dff))
  fit <- lm(Y ~ X + Tr + X*Tr,data=dff)
  dff_temp <- dff
  dff_temp$Tr <- 1
  Yhat1 <- predict(fit,newdata=dff_temp)
  dff_temp$Tr <- 0
  Yhat0 <- predict(fit,newdata=dff_temp)
  AIPW <-   mean(dff$Y*dff$Tr/Trhat - (dff$Tr - Trhat)/Trhat*Yhat1) - mean(dff$Y*(1-dff$Tr)/(1-Trhat) - (1-dff$Tr - (1-Trhat))/(1-Trhat)*Yhat0)
  return(AIPW)
}
ate_E <- function(dff){
  c(ate(dff[dff$E==1,]),ate(dff[dff$E==2,]),ate(dff[dff$E==3,]) )
}
overlap <- function(dff) {
  Trhat <- fitted.values(lm(Tr ~ X,data=dff))
  fit <- lm(Y ~ X + Tr + X*Tr,data=dff)
  dff_temp <- dff
  dff_temp$Tr <- 1
  Yhat1 <- predict(fit,newdata=dff_temp)
  dff_temp$Tr <- 0
  Yhat0 <- predict(fit,newdata=dff_temp)
  AIPW <-   (mean(dff$Y*dff$Tr*(1-Trhat) - (dff$Tr - Trhat)*(1-Trhat)*Yhat1) - mean(dff$Y*(1-dff$Tr)*(Trhat) - (1-dff$Tr - (1-Trhat))*(Trhat)*Yhat0))/mean(Trhat*(1-Trhat))
  return(AIPW)
}
overlap_E <- function(dff){
  c(overlap(dff[dff$E==1,]),overlap(dff[dff$E==2,]),overlap(dff[dff$E==3,]) )
}

lin_reg_E <- function(dff){
  reg <- function(e) {return(coef(lm(Y~ Tr + X, data = dff[dff$E==e,]))[2])}
  return(c(reg(1),reg(2),reg(3)))
}


wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,1000)
    estimator_targeted <- tms(ate_E,overlap_E,df)
    estimator_crossval <- crossval(list(ate_E,overlap_E),df,folds=10)
    estimator_cui <- cui(list(ate_E,overlap_E),df,folds=10)
    
    truth <- 1.5*s^2 + c(.1,.2,-.1)
    coverage <- mean(estimator_targeted[4:6] <= truth & estimator_targeted[7:9] >= truth)
    ret_vec <- c(  mean((estimator_targeted[1:3] - truth)^2),mean((estimator_crossval[1:3] - truth)^2),mean((ate_E(df) - truth)^2),mean((overlap_E(df) - truth)^2),coverage,mean((estimator_cui - truth)^2))
    return(ret_vec)
  }
  library(parallel)
  # The procedures sometimes returns NA. We exclude these cases.
  mat_results <- t(simplify2array(mclapply(1:iterations,wrapper,mc.cores=mc_cores)))
  ret_vec <- rowMeans(t(mat_results[complete.cases(mat_results),]))
  return(ret_vec)
}

wwrapper(0)
wwrapper(.3)
wwrapper(.6)
wwrapper(1)


mat <- sapply(seq(0,1,length.out=20),wwrapper)
mat

print("Average coverage")
mean(mat[5,])


save(mat,file="mat_hetero")
load(file="mat_hetero")


pdf("hetero.pdf",width = 6,height=5)
matplot(t(mat[c(1,2,3,4,6),]),type="l",x = seq(0,1,length.out=20),ylab="average MSE",xlab="gamma",ylim=c(0,.06),lwd = 2)
legend(.01,.06,c("targeted selection" , "cross-validation", "AIPW ATE", "AIPW overlap","selective ML"), col=1:4,fill=1:5,cex=1,bty="n")
dev.off()

mat[5,]
mean(mat[5,])

# Model selection for subset E==1

lin_reg_subset <- function(dff){
  return(lin_reg_E(dff)[1])
}

overlap_subset <- function(dff){
  return(overlap_E(dff)[1])
}

ate_subset <- function(dff){
  return(ate_E(dff)[1])
}


wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,1000)
    
    estimator_targeted <- tms_1D(ate_subset,overlap_subset,df)
    estimator_crossval <- crossval_1D(list(ate_subset,overlap_subset),df,folds=10)
    estimator_cui <- cui_1D(list(ate_subset,overlap_subset),df)
    
    truth <-  1.5*s^2 + .1
    ret_vec <- c(  mean((estimator_targeted[1] - truth)^2),mean((estimator_crossval - truth)^2),mean((ate_subset(df) - truth)^2),mean((overlap_subset(df) - truth)^2),mean(estimator_targeted[2] <= truth & estimator_targeted[3] >= truth), (estimator_cui - truth)^2)
    return(ret_vec)
  }
  library(parallel)
  # The procedures sometimes returns NA. We exclude these cases.
  mat_results <- t(simplify2array(mclapply(1:iterations,wrapper,mc.cores=mc_cores)))
  ret_vec <- rowMeans(t(mat_results[complete.cases(mat_results),]))
  return(ret_vec)
}


mat <- sapply(seq(0,1,length.out=20),wwrapper)
mat

save(mat,file="mat_hetero_single")
load(file="mat_hetero_single")

print("Average coverage")
mean(mat[5,])
mat

pdf("hetero_single.pdf",width = 6,height=5)
matplot(t(mat[c(1,2,3,4,6),]),type="l",x = seq(0,1,length.out=20),ylab="average MSE",xlab="gamma",ylim=c(0,.06),lwd = 2)
legend(.01,.06,c("targeted selection" , "cross-validation", "AIPW ATE", "AIPW overlap","selective ML"), col=1:4,fill=1:5,cex=1,bty="n")
dev.off()


###### Section 4.2: Instrumental variables and data fusion

gen_data <- function(s,n) {
  E <- sample(1:3,n,replace=TRUE)
  I <- rnorm(n)
  H <- rnorm(n)
  Tr <- .5*I + H + rnorm(n) 
  Y <- Tr - s^2*H + rnorm(n) + Tr*(E==1)*.1 + Tr*(E==2)*.2 + Tr*(E==3)*(-.1)
  I[  sample.int(n,size=round(n*1/2,0))] <- NA
  
  return(as.data.frame(cbind(E,Tr,I,Y)))
}


wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,n)
    f1 <- function(dff) cov(dff$I,dff$Y,use="complete.obs")/cov(dff$I,dff$Tr,use="complete.obs")
    f2 <- function(dff) coef(lm(dff$Y~dff$Tr))[2]
    f1_E <- function(dff) { c(f1(dff[dff$E==1,]),f1(dff[dff$E==2,]),f1(dff[dff$E==3,]) )}
    f2_E <- function(dff) { c(f2(dff[dff$E==1,]),f2(dff[dff$E==2,]),f2(dff[dff$E==3,]) )}
    
    
  
    estimator_targeted <- tms(f1_E,f2_E,df,R=1000)
    estimator_crossval <- crossval(list(f1_E,f2_E),df,folds=10)
    estimator_cui <- cui(list(f1_E,f2_E),df,folds=10)
    
    
    truth <- 1 + c(.1,.2,-.1)
    ret_vec <- c(  mean((estimator_targeted[1:3] - truth)^2),mean((estimator_crossval[1:3] - truth)^2),mean((f1_E(df) - truth)^2),mean((f2_E(df) - truth)^2),mean(estimator_targeted[4:6] <= truth & estimator_targeted[7:9] >= truth), mean((estimator_cui - truth)^2))
    ret_vec
    return(ret_vec)
  }
  library(parallel)
  # The procedures sometimes returns NA. We exclude these cases.
  mat_results <- t(simplify2array(mclapply(1:iterations,wrapper,mc.cores=mc_cores)))
  ret_vec <- rowMeans(t(mat_results[complete.cases(mat_results),]))
  
  return(ret_vec)
}

n <- 2000
mat <- sapply(seq(0,2,length.out=16),wwrapper)
mat

print("Average coverage")
mean(mat[5,])

save(mat,file="mat_IV")
load(file="mat_IV")


pdf("IV.pdf",width = 8,height=6)
matplot(t(mat[c(1,2,3,4,6),1:16]),type="l",x = seq(0,2,length.out=16),ylab="average MSE",xlab="gamma",ylim=c(0,.15),lwd = 2)
legend(0,.15,c("targeted selection" , "cross-validation", "IV", "OLS","selective ML"), col=1:4,fill=1:5,bty="n")
dev.off()




##### Section 4.3: Experiment with proxy outcomes

gen_data <- function(s,n) {
  E <- sample(1:3,n,replace=TRUE)
  Tr <- rnorm(n) <= .5
  X <- rnorm(n) <=  Tr
  Y <- X +  rnorm(n) + s^2*Tr*(E==1)*.1 + s^2*Tr*(E==2)*.2 + s^2*Tr*(E==3)*(-.1)
  
  return(as.data.frame(cbind(E,Tr,X,Y)))
}

wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,2000)
    f1 <- function(dff) coef(lm(Y~Tr,data=dff))[2]
    f2 <- function(dff) coef(lm(Y~X,data=dff))[2]*coef(lm(X~Tr,data=dff))[2]
    f1_E <- function(dff) { c(f1(dff[dff$E==1,]),f1(dff[dff$E==2,]),f1(dff[dff$E==3,]) )}
    f2_E <- function(dff) { c(f2(dff[dff$E==1,]),f2(dff[dff$E==2,]),f2(dff[dff$E==3,]) )}
    

    estimator_targeted <- tms(f1_E,f2_E,df,R=1000)
    estimator_crossval <- crossval(list(f1_E,f2_E),df)
    estimator_cui <- cui(list(f1_E,f2_E),df)
    
    estimator_crossval
    estimator_targeted
    
    # causal effect = c(.1,.2,-.1) + s^2 + pnorm(1) - pnorm(0)
    #truth <- .5*(pnorm(c(.1,.2,-.1))-pnorm(0)) + s^2  
    truth <- s^2*c(.1,.2,-.1) + pnorm(1) - pnorm(0)
    
    ret_vec <- c(  mean((estimator_targeted[1:3] - truth)^2),mean((estimator_crossval[1:3] - truth)^2),mean((f1_E(df) - truth)^2),mean((f2_E(df) - truth)^2),mean(estimator_targeted[4:6] <= truth & estimator_targeted[7:9] >= truth), mean( (estimator_cui-truth)^2))
    ret_vec
    return(ret_vec)
  }
  # Cross-validation is unstable and sometimes returns NA. We exclude these cases.
  mat_results <- t(simplify2array(mclapply(1:iterations,wrapper,mc.cores=mc_cores)))
  ret_vec <- rowMeans(t(mat_results[complete.cases(mat_results),]))
  return(ret_vec)
}

library(parallel)
mat <- sapply(seq(0,1.2,length.out=20),wwrapper)
mat

print("Average coverage")
mean(mat[5,])

save(mat,file="mat_proxy")
load(file="mat_proxy")

pdf("proxy.pdf",width = 8,height=6)
matplot(t(mat[c(1,2,3,4,6),]),type="l",x = seq(0,1.2,length.out=20),ylab="average MSE",xlab="gamma",ylim=c(0,.04),lwd = 2)
legend(.00,.04,c("targeted selection" , "cross-validation", "difference-in-means", "proxy model", "selective ML"), col=1:4,fill=1:5,bty="n")
dev.off()

mat[5,]
mean(mat[5,])

