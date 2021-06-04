##------------------Set parameters-----------------------------##
library(glmnet)
p = c(10, 100, 1000, 10000)
n = c(4000) # training data sample size
nSim <- 100
senario <- c("no_effect", "tree", "linear", "nonlinear")
grid <- expand.grid(n = n, p = p, senario = senario, stringsAsFactors = FALSE)
grid
nrow(grid)

##--------------- Data generating process---------------------##
datagen <-  function(n, p, senario){
  # generate X from uniform distribution
  X <- matrix(runif(n*p, min = -1, max = 1), nrow = n)
  
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

## LASSO

out_la <- lapply(1:nrow(grid), function(i){
  n <- grid$n[i]
  p <- grid$p[i]
  senario <- grid$senario[i]
  mse <- replicate(nSim,{
    training <- datagen(n=n, p=p, senario = senario)
    testing <- datagen(4000, p = p, senario = senario)
    
    df_lasso <- data.frame(training$Y, training$A, training$X)
    df_test_lasso <- data.frame(testing$Y, testing$A, testing$X)
    df_test_lasso_a1 <- df_test_lasso
    df_test_lasso_a0 <- df_test_lasso
    
    df_test_lasso_a1[,2] <- 1
    df_test_lasso_a0[,2] <- 0
    
    
    
    names(df_lasso)[names(df_lasso) == "training.Y"] <- "Y"
    names(df_lasso)[names(df_lasso) == "training.A"] <- "A"
    names(df_test_lasso_a1)[names(df_test_lasso_a1) == "testing.Y"] <- "Y"
    names(df_test_lasso_a0)[names(df_test_lasso_a0) == "testing.Y"] <- "Y"
    names(df_test_lasso_a1)[names(df_test_lasso_a1) == "testing.A"] <- "A"
    names(df_test_lasso_a0)[names(df_test_lasso_a0) == "testing.A"] <- "A"
    
    print(head(df_lasso))
    print(head(df_test_lasso_a0))
    print(head(df_test_lasso_a1))
    
    Y_lasso <- as.matrix(training$Y, ncol=1)
    colnames(Y_lasso)<- c("Y")
    

    A_lasso <- as.matrix(training$A, ncol=1)
    X_lasso <- as.matrix(training$X)
    
    sumx <- paste(c(colnames(df_lasso)[-c(1,2)]), collapse = "+")
    print(sumx)
    
    linear_het <- as.formula(paste("Y", paste("A * (", sumx, ")", sep = ""), sep = "~"))

    train_A <- model.matrix(object = linear_het, data = df_lasso)[,-1]
    
    test_A0 <- model.matrix(object = linear_het, data = df_test_lasso_a0)[,-1]
    
    test_A1 <- model.matrix(object = linear_het, data = df_test_lasso_a1)[,-1]
  
    lasso_A <- cv.glmnet(x = train_A,
                       y = df_lasso[,1],
                       family = "gaussian")
    print(lasso_A)
    
    pred_lasso_0 <- predict(lasso_A, newx = test_A0 , s=lasso_A$lambda.1se)
    pred_lasso_1 <- predict(lasso_A, newx = test_A1 , s=lasso_A$lambda.1se)
  
    pred_lasso <- pred_lasso_1 - pred_lasso_0
    #mean_la <-  mean(pred_lasso)
    #sd_la <- sd(pred_lasso)
    mse <- mean((testing$pdg -pred_lasso)^2)
  }
)

##--------RETURN OUTPUT------------##
data.frame(senario = senario, p=p , n = n, C.MSE = mean(mse), C.SD = sd(mse))
})
results_la <- do.call(rbind, out_la)
print(results_la)