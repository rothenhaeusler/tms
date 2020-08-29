# Targeted model selection
tms <- function(f1,f2,df,R=200) {
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
n <- 100
Tr <- rbinom(n,1,.5)
X <- .5*Tr + rnorm(n)
Y <- .5*X +  rnorm(n) + .01*Tr
df <-  as.data.frame(cbind(Tr,X,Y))

surrogate_estimator <- function(df) coef(lm(Y~X,data=df))[2]*coef(lm(X~Tr,data=df))[2]
difference_in_means <- function(df) coef(lm(Y~Tr,data=df))[2] 
tms(difference_in_means,surrogate_estimator,df)

# Targeted cross-validation
crossval <- function(list_of_funcs,df,folds=10){
  criterion_vec <- rep(0,10)
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
    
    criterion_vec <- criterion_vec + sapply(seq(0,1,length.out=10), function(s){ ( (1-s)*f1_in + s*f2_in  - test)^2})/10
    
  
  }
  
  f1_df <- list_of_funcs[[1]](df)
  f2_df <- list_of_funcs[[2]](df)
  
  whichmin <- which.min(criterion_vec)
  estimator_crossval <- (1-seq(0,1,length.out=10)[whichmin])*f1_df + (seq(0,1,length.out=10)[whichmin])*f2_df
  
  
  return(estimator_crossval)
}



#### Section 4.1 Observational studies: heterogeneous treatment effects

gen_data <- function(s,n) {
  X <- runif(n) <= .5
  Tr <- ifelse(X,runif(n) <= .7, runif(n) <= .05)
  Y <- X/2 + Tr + 3*s^2*Tr*X+  rnorm(n) 
  
  return(as.data.frame(cbind(Tr,X,Y)))
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

wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,1000)
    estimator_targeted <- tms(ate,overlap,df)
    estimator_crossval <- crossval(list(ate,overlap),df,folds=5)
    
    truth <- 1 + 1.5*s^2 
    
    ret_vec <- c((c(estimator_targeted[1],estimator_crossval,ate(df),overlap(df)) - truth)^2,  estimator_targeted[2]<= truth & estimator_targeted[3] >= truth)
    ret_vec
    return(ret_vec)
  }
  library(parallel)
  ret_vec <- rowMeans(simplify2array(mclapply(1:200,wrapper,mc.cores=4)))
  return(ret_vec)
}

mat <- sapply(seq(0,1,length.out=20),wwrapper)
mat

save(mat,file="mat_hetero")

pdf("hetero.pdf",width = 8,height=6)
matplot(t(mat[c(1,2,3,4),]),type="l",x = seq(0,1,length.out=20),ylab="average MSE",xlab="s",ylim=c(0,.04),lwd = 2)
legend(.01,.04,c("targeted selection" , "cross-validation", "AIPW ATE", "AIPW overlap"), col=1:4,fill=1:5)
dev.off()

mat[5,]
mean(mat[5,])

###### Section 4.2: Instrumental variables and data fusion


gen_data <- function(s,n) {
  I <- rnorm(n)
  H <- rnorm(n)
  Tr <- .5*I + H + rnorm(n)
  Y <- Tr - s^2*H + rnorm(n)
  I[  sample.int(n,size=round(n*1/2,0))] <- NA
  
  return(as.data.frame(cbind(Tr,I,Y)))
}


wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,n)
    f1 <- function(dff) cov(dff$I,dff$Y,use="complete.obs")/cov(dff$I,dff$Tr,use="complete.obs")
    f2 <- function(dff) coef(lm(dff$Y~dff$Tr))[2]
    
  
    estimator_targeted <- tms(f1,f2,df,R=1000)
    estimator_crossval <- crossval(list(f1,f2),df)
    
    truth <- 1
    
    ret_vec <- c((c(estimator_targeted[1],estimator_crossval,f1(df),f2(df)) - truth)^2,  estimator_targeted[2]<= truth & estimator_targeted[3] >= truth)
    ret_vec
    return(ret_vec)
  }
  library(parallel)
  ret_vec <- rowMeans(simplify2array(mclapply(1:200,wrapper,mc.cores=4)))
  return(ret_vec)
}

n <- 1000
mat <- sapply(seq(0,2,length.out=20),wwrapper)
mat

save(mat,file="mat_IV")
load(file="mat_IV")


pdf("IV.pdf",width = 8,height=6)
matplot(t(mat[c(1,2,3,4),1:16]),type="l",x = seq(0,2,length.out=16),ylab="average MSE",xlab="s",ylim=c(0,.065),lwd = 2)
legend(0,.065,c("targeted selection" , "cross-validation", "IV", "OLS"), col=1:4,fill=1:5)
dev.off()

mat[5,]
mean(mat[5,])


##### Section 4.3: Experiment with proxy outcomes

gen_data <- function(s,n) {
  Tr <- rnorm(n) <= .5
  X <- rnorm(n) <= Tr
  Y <- .5*X +  rnorm(n) + s^2*Tr
  
  return(as.data.frame(cbind(Tr,X,Y)))
}

wwrapper <- function(s){
  wrapper <- function(i) {
    set.seed(i)
    
    df <- gen_data(s,200)
    f1 <- function(dff) coef(lm(Y~Tr,data=dff))[2]
    f2 <- function(dff) coef(lm(Y~X,data=dff))[2]*coef(lm(X~Tr,data=dff))[2]

    estimator_targeted <- tms(f1,f2,df,R=1000)
    estimator_crossval <- crossval(list(f1,f2),df)
    estimator_crossval
    
    truth <- .5*(pnorm(1)-pnorm(0)) + s^2
    
    ret_vec <- c((c(estimator_targeted[1],estimator_crossval,f1(df),f2(df)) - truth)^2,  estimator_targeted[2]<= truth & estimator_targeted[3] >= truth)
    ret_vec
    return(ret_vec)
  }
  ret_vec <- rowMeans(sapply(1:200,wrapper))
  return(ret_vec)
}


library(parallel)
mat <- simplify2array(mclapply(seq(0,1.2,length.out=20),wwrapper,mc.cores=4))
mat

save(mat,file="mat_proxy")
load(file="mat_proxy")

pdf("proxy.pdf",width = 8,height=6)
matplot(t(mat[c(1,2,3,4),]),type="l",x = seq(0,1.2,length.out=20),ylab="average MSE",xlab="s",ylim=c(0,.04),lwd = 2)
legend(.00,.04,c("targeted selection" , "cross-validation", "difference-in-means", "proxy model"), col=1:4,fill=1:5)
dev.off()

mat[5,]
mean(mat[5,])
