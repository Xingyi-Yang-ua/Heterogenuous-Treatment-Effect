
library(dplyr)
##------------------Set parameters-----------------------------##
p = c(10, 100, 1000, 10000)
n = 4000 # training data sample size
nSim <- 100
senario <- c("no_effect", "tree", "linear", "nonlinear")
grid <- expand.grid(n = n, p = p, senario = senario, stringsAsFactors = FALSE)

print(grid)

##--------------- Data generating process---------------------##
datagen <-  function(n, p, senario){
  # generate X from uniform distribution
  X <- matrix(runif(n*p, min = -1, max = 1), nrow = n)
  
  print(head(X))
  
  # generate A with equal probability(p=0.5) and use labels -1 and +1.
  A <- rbinom(n, 1, 0.5)
  
  # generate Y from normal distribution
  
  if(senario == "no effect"){
    pdg <- 0
  }
  else if(senario == "tree"){
    pdg <- 2*(X[,1]<=0)*((X[,3]>=0.6)-0.3)+4
  }
  else if(senario == "linear"){
    pdg <- 3*(2-X[,1]+X[,2])
  }
  else {
    pdg <- 8+X[,1]^3+exp(X[,3]^2-X[,5])+1/(1+exp(-20*(X[,2]+4)))-(X[,4]+X[,6])^2
  }
  Y <- 0.3-0.5*X[,1] + 0.2*X[,2] + X[,3] + pdg*A + rnorm(n)
  data <- list(X=X, A=A, Y=Y, pdg = pdg)
  return(data)
}

##-----------------True optimal ITR----------------------##
trueITR <- function(X){
  pdg <- 2*(X[,1]<=0)*((X[,3]>=0.6)-1)+1
  trueoptITR <- (pdg>0)*2 - 1
  return(trueoptITR)
}


## Causal tree function
set.seed(20210510)
#rm(list = ls())
out_li <- lapply(1:nrow(grid), function(i){
  n <- grid$n[i]
  p <- grid$p[i]
  senario <- grid$senario[i]
  mse <- replicate(nSim,{
    training <- datagen(n=n, p=p, senario = senario)
    testing <- datagen(4000, p = p, senario = senario)
    
    df_li <- data.frame(training$Y, training$A, training$X)
    df_test_li <- data.frame(testing$Y, testing$A, testing$X)
    df_test_li_a1 <- df_test_li
    df_test_li_a0 <- df_test_li
    
    df_test_li_a1[,2] <- 1
    df_test_li_a0[,2] <- 0
    
    names(df_li)[names(df_li) == "training.Y"] <- "Y"
    names(df_li)[names(df_li) == "training.A"] <- "A"
    names(df_test_li_a1)[names(df_test_li_a1) == "testing.Y"] <- "Y"
    names(df_test_li_a0)[names(df_test_li_a0) == "testing.Y"] <- "Y"
    names(df_test_li_a1)[names(df_test_li_a1) == "testing.A"] <- "A"
    names(df_test_li_a0)[names(df_test_li_a0) == "testing.A"] <- "A"
    
    print(head(df_li))
    print(head(df_test_li_a0))
    print(head(df_test_li_a1))
    ##------------training---------------------------------------##
    #neg1Mat <- training1 %>% filter(A == 0) %>% select(-A)
    #pos1Mat <- training1 %>% filter(A == 1) %>% select(-A)
    ols_model <- lm(Y~., data = df_li)
    
    ##---------TESTING-------------------------------------------##
    ## prediction
    
    pred_li_0 <- predict(ols_model, newdata = df_test_li_a0)
    pred_li_1 <- predict(ols_model, newdata = df_test_li_a1)
    
    pred_li <- pred_li_1 - pred_li_0
    
    mse <- mean((testing$pdg - pred_li)^2)
  }
  )
  
  ##--------RETURN OUTPUT------------##
  data.frame(senario = senario, p=p , n = n, C.MSE = mean(mse), C.SD = sd(mse))
})
do.call(rbind, out_li)