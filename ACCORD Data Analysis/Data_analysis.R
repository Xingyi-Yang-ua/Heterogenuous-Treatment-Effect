library(tidyverse)
library(dplyr)
library(causalTree)
library(grf) ## causal forest 
library(glmnet)
library(zoo)
library(psych)
##################################################causal tree##############################################

# df_split: the splitting sample, used to build the tree
# df_est: the estimation sample, used to compute the average treatment effect in each leaf data 40%-40%-20% into splitting, estimation and validation samples
split_size <- floor(nrow(df_train_snp) * 0.5)
split_idx <- sample(nrow(df_train_snp), replace=FALSE, size=split_size)

# Make the splits
df_split_snp <- df_train_snp[split_idx,]
df_est_snp <- df_train_snp[-split_idx,]

##-----------------------------------------------training-------------------------------------------------##
fm <- as.formula(paste("change_hba1c ~ ", paste(predictor_names_snp_var, collapse= "+")))

## we estimate causal effects in the leaves of a given tree on an independent estimation sample rather than the data used to build the tree
ct_unpruned_snp <- honest.causalTree(
  formula = fm,            # Define the model
  data = df_split_snp,              # Subset used to create tree structure
  est_data = df_est_snp,            # pass the tree objects here, and replace the leaf estimates with new estimate calculated on the estimation sample.
  
  treatment = df_split_snp$glycemic_arm,       # Splitting sample treatment variable
  est_treatment = df_est_snp$glycemic_arm,     # Estimation sample treatment variable
  
  split.Rule = "CT",    # Define the splitting option, the outcome is the difference of treated mean and control mena in the leaf it belongs
  cv.option = "CT",    # Cross validation options, minimize the unbiased estimate of the validation mean-squared error.
  
  split.Honest=TRUE,    # Use honesty when splitting
  cv.Honest = TRUE,     # Use honesty when performing cross-validation
  
  minsize = 30,        # Min. number of treatment and control cases in each leaf
  HonestSampleSize = nrow(df_est_snp)
) # Number obs used in estimation after building the tree

## Cross validation
# Table of cross-validated values by tuning parameter.
ct_cptable_snp <- as.data.frame(ct_unpruned_snp$cptable)
ct_cptable_snp

# Obtain optimal complexity parameter to prune tree.
selected_cp_snp <- which.min(ct_cptable_snp$xerror) ## cross validation error
optim_cp_ct_snp <- ct_cptable_snp[selected_cp_snp, "CP"]

# Prune the tree at optimal complexity parameter.
ct_pruned_snp <- prune(tree = ct_unpruned_snp, cp=optim_cp_ct_snp)

rpart.plot(ct_pruned_snp)

tauhat_ct_est_snp <- predict(ct_pruned_snp, newdata=df_est_snp)
table(tauhat_ct_est_snp)

##---------------------------------------------------------Prediction------------------------------------------##
## Prediction on test data set
tauhat_ct_test_snp <- predict(ct_pruned_snp, newdata=df_test_snp)
table(tauhat_ct_test_snp)

mean(tauhat_ct_test_snp)
sd(tauhat_ct_test_snp)


## ----------------------------------------------------------performance---------------------------------------##
p <- mean(df_test_snp$glycemic_arm)
p
y_star <- ((df_test_snp$glycemic_arm-p)/(p*(1-p)))* (df_test_snp$change_hba1c-1)
y_star

#sample_ate <- with(df_test_snp, mean(censor_po[glycemic_arm == 1])- mean(censor_po[glycemic_arm == 0]))
mse_ct <- (y_star - tauhat_ct_test_snp)^2
mean(mse_ct)
sd(mse_ct)
hist(tauhat_ct_test_snp)

##################################################causal forest####################################################
##----------------------------------------------training---------------------------------------------------------##
df_train_snp$change_hba1c[is.na(df_train_snp$change_hba1c)] <- mean(df_train_snp$change_hba1c, na.rm=TRUE)
cf_snp <- causal_forest(
  X = as.matrix(df_train_snp[, predictor_names_snp_var]),
  Y = df_train_snp$change_hba1c,
  W = df_train_snp$glycemic_arm,
  num.trees = 5000,
  honesty = TRUE,
  mtry = ceiling(length(predictor_names_snp_var)/3),
  honesty.fraction = 0.5)

var_imp <- c(variable_importance(cf_snp))
names(var_imp) <- predictor_names_snp_var
sorted_var_imp <- sort(var_imp, decreasing=TRUE)
sorted_var_imp[1:50]

##----------------------------------------------prediction--------------------------------------------------------##
test_pred <- predict(cf_snp, newdata=as.matrix(df_test_snp[predictor_names_snp_var]), estimate.variance=TRUE)
tauhat_cf_test <- test_pred$predictions

hist(tauhat_cf_test)
mean(tauhat_cf_test)
sd(tauhat_cf_test)
##-----------------------------------------------performance----------------------------------------------------##
mse_cf <- (y_star - tauhat_cf_test)^2
mean(mse_cf)
sd(mse_cf)

#################################################S learner###############################################################
##----------------------------------------------training---------------------------------------------------------##
df_test_lasso_a1 <- df_test_snp
df_test_lasso_a0 <- df_test_snp

df_test_lasso_a1$glycemic_arm <- 1
df_test_lasso_a0$glycemic_arm <- 0

sumx <- paste(predictor_names_snp_var, collapse = "+")

linear_het <- as.formula(paste("change_hba1c", paste("glycemic_arm * (", sumx, ")", sep = ""), sep = "~"))

df_train_snp1 <- na.aggregate(df_train_snp)

df_train_snp1 <- data.frame(lapply(df_train_snp1, function(x) as.numeric(as.character(x))))

options(na.action='na.pass')
train_A <- model.matrix(object = linear_het, data = df_train_snp1)[,-1]
dim(train_A)

train_A[is.na(train_A)] <- 0

options(na.action='na.pass')
test_A0 <- model.matrix(object = linear_het, data = df_test_lasso_a0)[,-1]
options(na.action='na.pass')
test_A1 <- model.matrix(object = linear_het, data = df_test_lasso_a1)[,-1]

test_A0[is.na(test_A0)] <- 0
test_A1[is.na(test_A1)] <- 0

lasso_A <- cv.glmnet(x = train_A,
                     y = df_train_snp1$change_hba1c,
                     family = "gaussian")
##----------------------------------------------prediction--------------------------------------------------------##

pred_lasso_0 <- predict(lasso_A, newx = test_A0 , s=lasso_A$lambda.1se)
pred_lasso_1 <- predict(lasso_A, newx = test_A1 , s=lasso_A$lambda.1se)

pred_lasso <- pred_lasso_1 - pred_lasso_0
hist(pred_lasso)
mean(pred_lasso)
sd(pred_lasso)
##-----------------------------------------------performance----------------------------------------------------##
mse_lasso <- (y_star - pred_lasso)^2

mean(mse_lasso)
sd(mse_lasso)

#################################################T learner ###############################################################
##----------------------------------------------training---------------------------------------------------------##

coef <- predict(lasso_A, s = lasso_A$lambda.1se, type = "nonzero")
colnames <- colnames(train_A)
selected_vars <- colnames[unlist(coef)]

ols_formu <- as.formula(paste("change_hba1c", paste(append(selected_vars, "glycemic_arm"),collapse = "+"), sep = "~"))


##------------training---------------------------------------##
ols_model <- lm(formula = ols_formu, data = df_train_snp1)
summary(ols_model)

##----------------------------------------------prediction--------------------------------------------------------##

pred_li_0 <- predict(ols_model, newdata = df_test_lasso_a0)
pred_li_1 <- predict(ols_model, newdata = df_test_lasso_a1)

pred_ols <- pred_li_1 - pred_li_0

hist(pred_ols)

mean(pred_ols, na.rm = TRUE)
sd(pred_ols, na.rm = TRUE)

mse_li <- (y_star -pred_ols)^2
mean(mse_li, na.rm =TRUE)
sd(mse_li, na.rm = TRUE)


##--------------------------------------------Model comparison--------------------------------------------------------##
hist(pred_ols, main = "linear Regression", xlab= "Treatment Effect ", breaks = 6)
hist(pred_lasso, main = "LASSO", xlab= "Treatment Effect ", breaks =3)
hist(tauhat_cf_test, main = "Causal Forest", xlab= "Treatment Effect ", breaks =6)
hist(tauhat_ct_test_snp, main = "Causal Tree", xlab= "Treatment Effect ", breaks =6)

mse <- data.frame(Causal_Tree = mse_ct, Causal_Forest = mse_cf, LASSO = mse_lasso, Linear = mse_li)
mse_summary <- describe(mse)[,c("mean", "se")]
mse_summary



