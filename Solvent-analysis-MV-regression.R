#!/usr/local/bin/Rscript

#---
#title: "MV-regression-solvent-analysis"
#author: "Liliana C. Gallegos"
#email: "lilianac.gallegos@colostate.edu"
#---

## 1. Load libraries
library(dplyr)
library(corrplot)
library(cvq2)     # cross-validation analysis
library(car)     # ANOVA analysis
library(ggplot2) # plotting
library(cowplot) # compile plots
library(MuMIn)  # dredge
library(Metrics)  # rmse and mae


## 2. Load data
# setwd("/Volumes/Working-Yoshi/Research/Solvent/Solvent-updated")
all_labels <-read.csv("Keq-solvents-pca-split.csv", header=TRUE, sep = ",")   # Keq-solvents-updated.csv
data_all_unsc <- subset(all_labels, select = -solvent ) # all data
nododecane_labels <- all_labels[-c(2), ] ## No dodecane data

## Standardize datasets
data_all_numeric <- select_if(all_labels, is.numeric)
data_nododecane_numeric <- select_if(nododecane_labels, is.numeric)

data_allx = data.frame(scale(select(data_all_numeric, select = -c("dG") ), center = TRUE, scale = TRUE))
data_all = cbind(data_all_numeric$dG, data_allx)
data_nododecanex = data.frame(scale(select(data_nododecane_numeric, select = -c("dG") ), center = TRUE, scale = TRUE))
data_nododecane = cbind(data_nododecane_numeric$dG, data_nododecanex)

## Rename y-response variable
names(data_all)[names(data_all) == "data_all_numeric$dG"] <- "dG"
names(data_nododecane)[names(data_nododecane) == "data_nododecane_numeric$dG"] <- "dG"


## 3. (a) Figure S33 - Correlation plot using the Pearson algorithm on All data
cov1 = cov(data_all)
corr <- round(cor(data_all), 2)
res <- cor.mtest(data_all, conf.level = .95)
png(file = "figureS33.png",width = 200, height = 200, pointsize = 18, units='mm', res = 300)
corrplot <- corrplot(corr, p.mat = res$p, sig.level = c(.001, .01, .05), pch.cex = .7,
                     insig = "label_sig", pch.col = "white", type = "upper", 
                     order = "hclust", tl.col = "black", tl.srt = 45, 
                     method="square", cl.align = "l" , tl.cex = .5) 
dev.off()
r2_values = round(corr^2, 2)
print(r2_values)  # shows pairwise R2-correlations to each variable

## (b) Figure S34 - Highest correlations to response: Sig2, beta, and gamma. V-descriptor plotted for comparison. 
Sig2_plot <- ggplot(data_all, aes(data_all$Sig2, data_all$dG)) + geom_point() +
  xlab("Sig2") + ylab(" \u2206G (kcal/mol)")+ ggtitle( paste("\u2206G ~ Sig2  \nR2 =", r2_values[1,14]) ) + stat_smooth(method = "lm", se=FALSE,col= "black") + theme_bw(); Sig2_plot 

beta_plot <- ggplot(data_all, aes(data_all$beta, data_all$dG)) + geom_point() +
  xlab("beta") + ggtitle( paste("\u2206G ~ beta  \nR2 =", r2_values[1,3]) ) + stat_smooth(method = "lm", se=FALSE,col= "black") + theme_bw() + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()); beta_plot 

V_plot <- ggplot(data_all, aes(data_all$V, data_all$dG)) + geom_point() +
  xlab("V") + ylab(" \u2206G (kcal/mol)")+ ggtitle( paste("\u2206G ~ V  \nR2 =", r2_values[1,18]) ) + stat_smooth(method = "lm", se=FALSE,col= "black") + theme_bw(); V_plot

gamma_plot <- ggplot(data_all, aes(data_all$gamma, data_all$dG)) + geom_point() +
  xlab("gamma") + ggtitle( paste("\u2206G ~ gamma  \nR2 =", r2_values[1,4]) ) + stat_smooth(method = "lm", se=FALSE,col= "black") + theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()); gamma_plot

theme_set(theme_cowplot(font_size=12))
FigS34 <- plot_grid(Sig2_plot, beta_plot, V_plot, gamma_plot, align = "h"); FigS34
# save_plot("figureS34.jpeg", FigS34, ncol = 2, nrow = 2)


## 4. Select best model
## (a) Forward Stepwise Model Selection
model_data = data_all  #select from data_all or data_nododecane

null1=lm(dG~1, data=model_data)
full1=lm(dG~ . , data=model_data) # to all parameters
best1=step(null1, scope=list(lower=null1, upper=full1), direction="forward", trace=TRUE);

## (b) Dredge Model Selection ranked by AIC
model_data1 = data_all  # select from data_all or data_nododecane
model_data = subset(model_data1, select = -c(2:7) ) #To run dredge faster, modeled only selected amount of descriptors.

full_model <- lm(dG ~ ., data = model_data);
options(na.action = "na.fail")
all_models <- dredge(full_model, extra = c("R^2", F = function(x)
  summary(x)$fstatistic[[1]]))

top3 = head(all_models, 3); top3 # shows top 3 models 

## (c) Table S18 - Calculates variable importance for ALL possible models
## Weights and and Number of times the descriptors shown in the models
importance(all_models)  

## (d) All Models suggested to consider for ALL samples:
model1 <- lm(dG ~ Sig2 + V + sum + dipole, data=data_all); #summary(model1)
model2 <-lm(dG ~ Sig2 + V + sum , data=data_all); #summary(model2) 
model3 <- lm(dG ~ Sig2 + V, data =data_all); #summary(model3)
model4 <- lm(dG ~ Sig2 + volume, data=data_all); #summary(model4)
model5 <- lm(dG ~ Sig2 + area, data=data_all); #summary(model5)

## All Models suggested to consider for 19-solvent samples excluding dodecane:
model6 <- lm(dG ~ Sig2 + area + V, data=data_nododecane); #summary(model6)
model7 <- lm(dG ~ Sig2 + V, data=data_nododecane); #summary(model7)
model8 <- lm(dG ~ Sig2 + MV_boltz, data=data_nododecane); #summary(model8)

## (e) Table S17 - R2 and RMSE of All models
r2 <- function(model) {
  rr <- summary(model)$r.squared
  return(signif(rr,2)) }

rmse <- function(model) { 
  r <- sqrt(mean(model$residuals^2))
  return(signif(r,2)) }

MAE <- function(model) { 
  m <- mean(abs(model$residuals))
  return(signif(m,2)) }

model_list <- c()
rmse_list <- c()
mae_list <- c()
rr_list <- c()
for (i in 1:8) {
  n <- paste("model",i, sep="")
  nam <- eval(parse(text=n))
  print(nam$call$formula)
  
  # Compile data
  model_list[[i]] <- n
  rmse_list[[i]] <- rmse(nam)
  mae_list[[i]] <- MAE(nam)
  rr_list[[i]] <- r2(nam)
}; 
Table_S17 = tibble(model_list, rr_list, rmse_list, mae_list); Table_S17
# write.csv(Table_S17, "TableS17.csv", row.names = FALSE)

## 5. Select best Train-to-Test split 
## Table S20 - Compare R2 and RMSE values in different PCA clusters
cluster <- c()
formula <- c()
R2_train <- c()
R2_rmse <- c()
Q2_test <- c()
Q2_rmse <- c()

for (i in 5:14) {
  ## Select Train:Test split by 'k_cluster' column
  column <- paste('X',i,'.k_clusters', sep = "")
  
  data_presplit <- split(data_all_unsc, data_all_unsc [column])
  data_tr2 <- data_presplit[["TRAIN"]]
  data_tt2 <- data_presplit[["TEST"]]
  
  data_tr <- select_if(data_tr2, is.numeric) 
  data_tt <- select_if(data_tt2, is.numeric) 
  
  datatr_x <- select(data_tr, select = -c("dG")) # keep only x-variables
  datatt_x <- select(data_tt, select = -c("dG")) # keep only x-variables
  datatr_y <- data_tr["dG"]
  
  ## Standardizes
  trainMean <- apply(datatr_x,2,mean)
  trainSd <- apply(datatr_x,2,sd)
  datatr_scaled <- sweep(sweep(datatr_x, 2L, trainMean), 2, trainSd, "/")
  data <- cbind(datatr_y,datatr_scaled )
  
  datatt_y <- data_tt["dG"]
  datatt_scaled <- sweep(sweep(datatt_x, 2L, trainMean), 2, trainSd, "/") 
  data_tt <- cbind(datatt_y,datatt_scaled ) 
  
  assign(paste0("data_tr", i), data)  
  assign(paste0("data_tt", i), data_tt)
  
  ## Final model 
  final_model <- lm(dG ~ Sig2 + V, data = data); summary(final_model)
  fit = final_model[["call"]][["formula"]] 
  LM_formula = paste("dG = ", round(coefficients(final_model)[1],2), " + ", 
                paste( sprintf("%.2f * %s",coefficients(final_model)[-1],names(coefficients(final_model)[-1])), collapse=" + ") )

  data_train <- data[ c(names(final_model$model)[ ])] # selects the model descriptors for training data
  data_test <- data_tt[ c(names(final_model$model)[ ])] # selects the model descriptors for testing data
  result_q2 <- q2(data_train, data_test, fit) 
  # print(result_q2)
  
  # Collect and compile data
  cluster[[i]] <- i
  formula[[i]] <- LM_formula
  R2_train[[i]] <- signif(result_q2@result[["fit"]][["r2"]], 2)
  R2_rmse[[i]] <- signif(result_q2@result[["fit"]][["rmse"]], 2)
  Q2_test[[i]] <- signif(result_q2@result[["pred"]][["q2"]], 2)
  Q2_rmse[[i]] <- signif(result_q2@result[["pred"]][["rmse"]], 2)
}
Table_S20 = na.omit(tibble(cluster, formula, R2_train, R2_rmse, Q2_test, Q2_rmse)); Table_S20
# write.csv(Table_S20, "TableS20.csv", row.names = FALSE)  # save to csv


## 6. Final analysis based on Best model and Best Train:Test split
# f = sapply(data_all_unsc, levels)  # SHOW the levels for categorical variables
# sample_split = 'X12.k_clusters'   # SELECT sample column to split by
# 
# data_presplit <- split(data_all_unsc, data_all_unsc [sample_split])
# data_tr2 <- data_presplit[["TRAIN"]]
# data_tt2 <- data_presplit[["TEST"]]
# 
# data_tr <- select_if(data_tr2, is.numeric) 
# data_tt <- select_if(data_tt2, is.numeric) 
# 
# datatr_x <- select(data_tr, select = -c("dG")) # keep only x-variables
# datatt_x <- select(data_tt, select = -c("dG")) # keep only x-variables
# datatr_y <- data_tr["dG"]
# 
# ## Standardize Train and Test data sets
# trainMean <- apply(datatr_x,2,mean) # Test dataset based on Train's mean
# trainSd <- apply(datatr_x,2,sd)  # Test dataset based on Train's standard deviation
# datatr_scaled <- sweep(sweep(datatr_x, 2L, trainMean), 2, trainSd, "/")
# data <- cbind(datatr_y,datatr_scaled )
# 
# datatt_y <- data_tt["dG"]
# datatt_scaled <- sweep(sweep(datatt_x, 2L, trainMean), 2, trainSd, "/") 
# data_tt <- cbind(datatt_y,datatt_scaled ) 


## 7. Final Model with Train data set only
FINAL_TR_DATA = data_tr12
FINAL_TT_DATA = data_tt12
final_model <- lm(dG ~ Sig2 + V, data = FINAL_TR_DATA); summary(final_model)


## 8. ANOVA Analysis
Anova(final_model, type=3)


## 9. Cross Validation Analysis
fit = final_model[["call"]][["formula"]] 
## (a) Q2-KFold, RMSE
data_train <- FINAL_TR_DATA[ c(names(final_model$model)[ ])] # selects the model descriptors
k_q2 = cvq2( data_train, fit, nFold = 5, nRun = 5); k_q2

## (b) Q2-Test, RMSE-test
data_train = FINAL_TR_DATA[ c(names(final_model$model)[ ])] # selects the model descriptors for training data
data_test = FINAL_TT_DATA[ c(names(final_model$model)[ ])] # selects the model descriptors for testing data
result_q2 = q2(data_train, data_test, fit  ); result_q2


## 10. R2-Test, k and k' slopes, and R2-ratios
## (a) R2 > 0.6
FINAL_TT_DATA$y_hat = predict(final_model, FINAL_TT_DATA)
r2_test = round( cor(FINAL_TT_DATA$y_hat, FINAL_TT_DATA$dG)^2, 3); cat('\n R2-Test: ',r2_test)

## TableS21 - Train
solvs_train = all_labels[rownames(FINAL_TR_DATA),"solvent"]
solvs_test = all_labels[rownames(FINAL_TT_DATA),"solvent"]
solvs = append(solvs_train, solvs_test); solvs

Exp_dG_train = all_labels[rownames(FINAL_TR_DATA),"dG"]
Exp_dG_test = all_labels[rownames(FINAL_TT_DATA),"dG"]
Exp_dGs = append(Exp_dG_train, Exp_dG_test); Exp_dGs

Calc_dG_train = round(as.vector(final_model$fitted.values), 2)
Calc_dG_test = round(FINAL_TT_DATA$y_hat, 2)
Calc_dGs = append(Calc_dG_train, Calc_dG_test); Calc_dGs

Table_S21_train = tibble(solvs_train, Exp_dG_train, Calc_dG_train); Table_S21_train
Table_S21_test = tibble(solvs_test, Exp_dG_test, Calc_dG_test); Table_S21_test

mae_dG_train = round(mae(Table_S21_train$Exp_dG_train, Table_S21_train$Calc_dG_train),2) ; mae_dG_train
sd_dG_train = round(sd(Table_S21_train$Calc_dG_train ),2); sd_dG_train
mae_dG_test = round(mae(Table_S21_test$Exp_dG_test, Table_S21_test$Calc_dG_test),2) ; mae_dG_test
sd_dG_test= round(sd(Table_S21_test$Calc_dG_test ),2); sd_dG_test

# write.csv(Table_S21_train[ order(Table_S21_train$Exp_dG_train), ], "TableS21-train.csv")
# write.csv(Table_S21_test[ order(Table_S21_test$Exp_dG_test), ], "TableS21-test.csv")


## (b) k/k' slopes between 0.85 and 1.15 
k = round( sum( (FINAL_TT_DATA$dG)*(FINAL_TT_DATA$y_hat) ) / sum( (FINAL_TT_DATA$y_hat)^2 ) , 2); cat('\n k slope: ', k)
k_pred = round( sum( (FINAL_TT_DATA$dG)*(FINAL_TT_DATA$y_hat) ) / sum( (FINAL_TT_DATA$dG)^2 ) , 2); cat('\n k-predicted slope: ', k_pred)

## (c) ratios of r2-difference < 0.1
FINAL_TT_DATA$y0 = k*(FINAL_TT_DATA$y_hat) # calculation for y-thru intercept=0
FINAL_TT_DATA$y0_pred = k_pred*(FINAL_TT_DATA$dG) # calculation for y_predicted thru intercept=0
## R2-0 = r2 actual at intecept=0 
r2.0.actual = 1 - ( sum((FINAL_TT_DATA$y_hat - FINAL_TT_DATA$y0)^2)/sum((FINAL_TT_DATA$y_hat- mean(FINAL_TT_DATA$y_hat))^2) ) 
## R'2-0 = r2 predicted at intecept=0 
r2.0.predicted = 1 - ( sum((FINAL_TT_DATA$dG - FINAL_TT_DATA$y0_pred)^2)/sum((FINAL_TT_DATA$dG- mean(FINAL_TT_DATA$dG))^2) )
cat('\n Ratio of R2 difference: ', round((r2_test - r2.0.actual)/ r2_test, 3) )
cat('\n Ratio of R2-predicted difference :', round((r2_test - r2.0.predicted)/ r2_test, 3) )



