library(causalTree)

##------------------Set parameters-----------------------------##
p = c(10, 100, 1000, 10000)
n = 4000 # training data sample size
nSim <- 100
senario <- c("no_effect", "tree", "linear", "nonlinear")
grid <- expand.grid(n = n, p = p, senario = senario, stringsAsFactors = FALSE)
grid
nrow(grid)


## --------------------------------------------------Data generating process--------------------------------##
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
    pdg <- 3*(1-X[,1]+X[,2])
  }
  else {
    pdg <- 8+X[,1]^3+exp(X[,3]^2-X[,5])+1/(1+exp(-20*(X[,2]+4)))-(X[,4]+X[,6])^2
  }
  Y <- 0.3-0.5*X[,1] + 0.2*X[,2] + X[,3] + pdg*A + rnorm(n)
  data <- list(X=X, A=A, Y=Y, pdg = pdg)
  return(data)
}

## -----------------Causal tree-------------------------------##

out_ct <- lapply(1:nrow(grid), function(i){
  n <- grid$n[i]
  p <- grid$p[i]
  senario <- grid$senario[i]
  mse <- replicate(nSim,{
    training <- datagen(n=n, p=p, senario = senario)
    testing <- datagen(4000, p = p, senario = senario)
    
    # Diving the data 40%-40%-20% into splitting, estimation and validation samples
    
    df <- data.frame(training$Y, training$A, training$X)
    df_test <- data.frame(testing$Y, testing$A, testing$X)
    
    #print(training$A)
    split_size <- floor(nrow(df) * 0.5)
    split_idx <- sample(nrow(df), replace=FALSE, size=split_size)
    
    # Make the splits
    df_split <- df[split_idx,]
    
    #print(colnames(df_split))
    df_est <- df[-split_idx,]
    
    predictor_names <- colnames(df)[-c(1,2)]
    
    #print(predictor_names)
    
    fmla <- as.formula(paste("training.Y ~ ", paste(predictor_names, collapse= "+")) )
    
    #print(fmla)
    
    ##------------training---------------------------------------##
    ct <- honest.causalTree(
      formula = fmla,            # Define the model
      data = df_split,              # Subset used to create tree structure
      est_data = df_est,            # pass the tree objects here, and replace the leaf estimates with new estimate calculated on the estimation sample.

      treatment = df_split$training.A,       # Splitting sample treatment variable
      est_treatment = df_est$training.A,     # Estimation sample treatment variable
      
      split.Rule = "CT",    # Define the splitting option, the outcome is the difference of treated mean and control mena in the leaf it belongs
      cv.option = "TOT",    # Cross validation options, minimize the unbiased estimate of the validation mean-squared error.
      
      split.Honest=TRUE,    # Use honesty when splitting
      cv.Honest = TRUE,     # Use honesty when performing cross-validation
      
      minsize = 40,        # Min. number of treatment and control cases in each leaf
      HonestSampleSize = nrow(df_est)) # Number obs used in estimation after building the tree
    
    ##---------TESTING-------------------------------------------##
    pred_value <- predict(ct, newdata = df_test)
    #mse <- mean(testing$pdg -pred_value)^2
    mse <- mean(pred_value)
    
  }
  )
  
  ##--------RETURN OUTPUT------------##
  data.frame(senario = senario, p=p , n = n, C.MSE = mean(mse), C.SD=sd(mse))
})
do.call(rbind, out_ct)