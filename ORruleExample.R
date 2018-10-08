library(gettheta) 

### function for the weighted average of nonparametric estimators of sensitivity and specificity 
### in the rule based on X alone  
increment1 <- function(thres,data_eval){
  ind_D0 <- which(data_eval[,"D"]==0)
  ind_D1 <- which(data_eval[,"D"]==1)
  f <- function(delta){
    w*mean(ifelse(data_eval[ind_D1,"x"]>delta,1,0))+(1-w)*mean(ifelse(data_eval[ind_D0,"x"]<=delta,1,0))
  } 
  result <- f(thres)
  return(result)
}

### function for the weighted average of nonparametric estimators of sensitivity and specificity 
### in the OR rule based on X and Y
increment2 <- function(thres,data_eval){
  ind_D0 <- which(data_eval[,"D"]==0)
  ind_D1 <- which(data_eval[,"D"]==1)
  f <- function(delta){
    w*mean(ifelse(data_eval[ind_D1,"x"]>delta[1,1]|data_eval[ind_D1,"y"]>delta[1,2],1,0))+
      (1-w)*mean(ifelse(data_eval[ind_D0,"x"]<=delta[1,1]&data_eval[ind_D0,"y"]<=delta[1,2],1,0))
  }
  result <- f(thres)
  return(result)
}

### function for obtaining cross-validated estimators of the incremental value ###
getCV <- function(D,X,Y,num){
  est1 <- est2OR <- rep(NA,num)
  for (i in 1:num){
    oo.case <- (1:length(D))[D==1]
    oo.ctrl <- (1:length(D))[D==0]
    size.case <- floor(sum(D)/num)
    size.ctrl <- floor(sum(1-D)/num)
    if (i < num){
      ind.case <- oo.case[((i-1)*size.case+1):(i*size.case)]
      ind.ctrl <- oo.ctrl[((i-1)*size.ctrl+1):(i*size.ctrl)]
    } else {
      ind.case <- oo.case[((i-1)*size.case+1):sum(D)]
      ind.ctrl <- oo.ctrl[((i-1)*size.ctrl+1):sum(1-D)]
    }
    ind <- c(ind.case,ind.ctrl)
    
    thres1 <-  threshold_1d(D[-ind],X[-ind],w)$thresholds[1]
    thres2OR <- threshold_2d(D[-ind],X[-ind],Y[-ind],w,1)$thresholds[1,]
    
    dat <- data.frame(D=D,x=X,y=Y)
    est1[i] <-  increment1(thres1,dat[ind,])
    est2OR[i] <- increment2(thres2OR,dat[ind,])   
  }   
  return(c(mean(est1),mean(est2OR)))
}

bootsample <- function(x){
  out <- sample(x,length(x),replace=TRUE)
  return(out)
}

getCV_ind <- function(t){
  out <- getCV(D[t],x[t],y[t],fold)
  return(out)
}

N <- 10000    ## population sample size 
n <- 100      ## equal numbers of cases and controls sampled from the population
nsample <- 1000 ## bootstrap replications 
fold <- 10      ## 10-fold cross-validation 
c1 <- qnorm(0.5,0,1)    ## threshold for the established biomarker X
c2 <- Inf               ## threshold for a new biomarker Y
beta0 <- 0.064  ## coefficients in a logic model for the risk of disease outcome D
beta1 <- 0.072 
w <- 0.5    ## weight
mu_x <- 0   ## mean of X
mu_y <- 0   ## mean of Y
sd_x <- 1   ## standard deviation of Y
sd_y <- 1   ## standard deviation of Y


###### Generate data based on a logic model for the risk of D ###### 
X <- rnorm(N,mu_x,sd_x)
Y <- rnorm(N,mu_y,sd_y)
p <- beta0 + beta1*ifelse(X>c1|Y>c2,1,0)
D_pop <- rbinom(N,1,p)

###### Sample equal numbers of cases and controls ######
ind_D0 <- sample(which(D_pop==0),n,replace=FALSE)
ind_D1 <- sample(which(D_pop==1),n,replace=FALSE)
ind <- c(ind_D0,ind_D1)
x <- X[ind]
y <- Y[ind]
D <- D_pop[ind]
  
###### Estimate the incremental value #######
## function threshold_1d is for performance estimator using the rule based on X alone 
## function threshold_2d is for performance estimator using the rule based on X and Y
result1_data <-  threshold_1d(D,x,w)$maxtheta 
result2_data <- threshold_2d(D,x,y,w,1)$maxtheta 
incre_data <- result2_data-result1_data  ## estimator of incremental value

###### Bootstrap without cross-validaiton ######
result1_sample <- result2_sample <- incre_sample <- rep(NA,nsample)
for (i in 1:nsample){
  ind_sample_D0 <- sample(1:n,n,replace=TRUE)
  ind_sample_D1 <- sample((n+1):(2*n),n,replace=TRUE)
  ind_sample <- c(ind_sample_D0, ind_sample_D1)
  result1_sample[i] <-  threshold_1d(D[ind_sample],x[ind_sample],w)$maxtheta
  result2_sample[i] <- threshold_2d(D[ind_sample],x[ind_sample],y[ind_sample],w,1)$maxtheta
  incre_sample[i] <- result2_sample[i]-result1_sample[i] ## estimator of incremental value on 1000 bootstrap samples 
}

###### Estimate the incremental value with cross-validation #######  
outt <- getCV(D,x,y,fold)
result1_data_cv <- outt[1]
result2_data_cv <- outt[2]
incre_data_cv <- result2_data_cv - result1_data_cv  ## estimator of incremental value with cross-validation 

###### Bootstrap with cross-validaiton ######  
ind_sample_D0 <- apply(matrix(rep(1:n,nsample),nrow=n),2,bootsample)
ind_sample_D1 <- apply(matrix(rep((n+1):(2*n),nsample),nrow=n),2,bootsample)
ind_sample <- rbind(ind_sample_D0, ind_sample_D1)
outt.boot <- apply(ind_sample,2,getCV_ind)
result1_sample_cv <- outt.boot[1,]
result2_sample_cv <- outt.boot[2,]
incre_sample_cv <- result2_sample_cv-result1_sample_cv


###### One-sided confidence interval based on Bootstrap ######
lb_1side_em <- 2*incre_data-quantile(incre_sample,0.95) ## Empirical bootstrap 
lb_1side_quan <- quantile(incre_sample,0.05)            ## Percentile bootstrap 

###### Wald-test ######
wald_test <- incre_data_cv/sd(incre_sample_cv)
ifelse(wald_test>=qnorm(0.95,0,1),1,0)  ## reject the null if return 1, fail to reject if return 0

###### Fuzzy p-value ######
sigma <- sqrt(2*n)*sd(result1_sample)
xx_boot <- rnorm(1,0,sigma)
pvalue_onesided_boot <- pnorm(sqrt(2*n)*(-incre_data_cv)+xx_boot,0,sigma)  


          
