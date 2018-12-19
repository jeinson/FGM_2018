# This script runs through the test data for RIVER

library(RIVER)
library(Biobase)
library(glmnet)
library(pROC)

filename <- system.file("extdata", "simulation_RIVER.gz", 
                        package="RIVER")


dataInput <- getData(filename) # import experimental data
Feat <- t(Biobase::exprs(dataInput)) # genomic features (G)
Out <- as.numeric(unlist(dataInput$Outlier))-1 # outlier status (E)

FeatAll <- t(exprs(dataInput))
# all outlier status (E)
OutAll <- as.numeric(unlist(dataInput$Outlier))-1
# G for training models
FeatTrng <- t(exprs(dataInput[,is.na(dataInput$N2pair)]))
# E for training models
OutTrng <- as.numeric(unlist(dataInput$Outlier
                             [is.na(dataInput$N2pair)]))-1
# G for test
FeatTest <-
  t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
          [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
          exprs(dataInput[,!is.na(dataInput$N2pair)])
          [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))

# E for test (1st and then 2nd individuals from N2 pairs)
OutTest1 <-
  as.numeric(unlist(
    c(dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
      dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1

# E for test (2nd and then 1st individuals from N2 pairs)
OutTest2 <-
  as.numeric(unlist(
    c(dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
      dataInput$Outlier[!is.na(dataInput$N2pair)]
      [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1

## Standardization
meanFeat <- apply(FeatAll, 2, mean)
sdFeat <- apply(FeatAll,2,sd)
FeatAll <- scale(FeatAll, center=meanFeat, scale=sdFeat)
FeatTrng <- scale(FeatTrng, center=meanFeat, scale=sdFeat)

## Search a best lambda from a multivariate logistic regression
##         with outlier status with 10 cross-validation
## GAM (genomeic annotation model)
costs <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)
logisticCV <- cv.glmnet(FeatTrng, as.vector(OutTrng), lambda=costs,
                        family="binomial", alpha=0, nfolds=10)

## Compute a P(FR | G) for all data
postprobTest <-
  predict(logisticCV, FeatTest, s="lambda.min", type="response")

theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)
verbose = T
pseudoc = 50

emModel <- integratedEM(FeatTrng, OutTrng, logisticCV$lambda.min,
                        logisticCV$glmnet.fit, pseudoc,
                        theta_init, costs, verbose)


# ## Generate G data for test data (Revised)
FeatTest <- scale(FeatTest, center=meanFeat, scale=sdFeat)

## Compute P(FR | G, E)
dup.post <- testPosteriors(FeatTest, OutTest1, emModel)

## Check performance of models with N2 pairs
RIVER.roc <- roc(OutTest2, dup.post$posterior[,2]) # RIVER
GAM.roc <- roc(OutTest2, as.numeric(postprobTest)) # GAM

# Naive Logistic Regression Model
LR.model <- glm(OutTrng ~ ., family = "binomial", data = as.data.frame(cbind(OutTrng, FeatTrng)))
LR.posteriors <- predict(LR.model, newdata = as.data.frame(FeatTest), type = "response")
LR.roc <- roc(OutTest2, LR.posteriors)

# BBVI Model
source("Documents/FGM/logistic_BBVI_function.R")
BBVI.model <- logistic_reg_BBVI(FeatTrng, OutTrng, S =3)

# Evaluate on the held out data
# Sample from the parameter's posterior distribution
BBVI.beta <- sapply(1:100, function(s) rnorm(BBVI.model$mu, BBVI.model$mu, BBVI.model$sigma))

logistic.BBVI.prediction <- list()
for(j in 1:ncol(BBVI.beta)){
  logistic.BBVI.prediction[[j]] <- as.vector(sigmoid(t(BBVI.beta[,j]) %*% t(FeatTest)))
}
  
BBVI.roc <- sapply(logistic.BBVI.prediction, function(x) roc(OutTest2, x))
par(mar=c(6.1, 6.1, 4.1, 4.1))
plot(NULL, xlim=c(0,1), ylim=c(0,1), 
     xlab="False positive rate", ylab="True positive rate")
abline(0, 1, col="gray")
apply(BBVI.roc, 2, function(r){
  #print(r$sensitivities)
  lines(1-r$specificities, r$sensitivities,
        type="s")
})
### Plot a ROC Curve for each of the testing methods

pdf("Documents/FGM/figures/Figure2.pdf", width = 10, height = 8)
par(mar=c(6.1, 6.1, 4.1, 4.1))
plot(NULL, xlim=c(0,1), ylim=c(0,1), 
     xlab="False positive rate", ylab="True positive rate")
abline(0, 1, col="gray")
lines(1-RIVER.roc$specificities, RIVER.roc$sensitivities, 
      type="s", col='dodgerblue', lwd=2)
lines(1-GAM.roc$specificities, GAM.roc$sensitivities, 
      type="s", col='mediumpurple', lwd=2)
lines(1-LR.roc$specificities, LR.roc$sensitivities,
      type="s", col='red', lwd=2)
lines(1-BBVI.roc[,1]$specificities, BBVI.roc[,1]$sensitivities, 
      type="s", col='lightgreen', lwd=2)
legend(0.7,0.5,
       c(paste("RIVER: ", round(RIVER.roc$auc, 3)),
         paste("LASSO: ", round(GAM.roc$auc, 3)), 
         paste("LR: ", round(LR.roc$auc, 3)),
         paste("BBVI-LR: ", round(BBVI.roc[,1]$auc, 3)))
       , lty=c(1,1), lwd=c(2,2),
       col=c("dodgerblue","mediumpurple", "red", "lightgreen"), 
       bty="n")
dev.off()

# Look at the Mean auc across manu samples from the posterior
pdf("Documents/FGM/figures/histogram_of_AUCs.pdf", width = 5, height = 5)
hist(unlist(BBVI.roc["auc",]), 
     main = bquote(atop("Histogram of AUCs from many draws of" ~ beta ,"from the posterior" ~ lambda, sigma)),
     xlab = "AUC", 
     col = "cornflowerblue")
box()
dev.off()
