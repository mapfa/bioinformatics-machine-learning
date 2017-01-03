#Authors: Lara Vomfell, Marten Pfannenschmidt

### Load Libraries ###
library(randomForest)
library(ROCR)
library(plyr)
library(caret)
library(paprbag)

##########
# Step 1 #
##########

### Check Data ###
data("trainingData")
nrow(trainingData)
sum(trainingData$Labels == T) # balanced classes

##########
# Step 2 #
##########

### create testing and training set ###
set.seed(111)
trunc.data <- createDataPartition(trainingData$Labels, p = 0.5, list = FALSE)
tr <- trainingData[trunc.data,]
te <- trainingData[-trunc.data,]

### train random forest ###

 #find optimal mtry value
mtry.tune <- tuneRF(tr[,-1], tr$Labels, ntreeTry = 1000, stepFactor = 1.5, improve = 0.5, trace = T, plot = T, do.trace = T)
 #suggests mtry = 27

 #create random forest
tr.rf <- randomForest(Labels~., data = tr, importance = T, mtry = 27, keep.forest=TRUE)
print(tr.rf)

 #apply to testset
te.rf <- predict(tr.rf, newdata=te[,-1])

 #confusion matrix, sensitivity, specificity
table(predicted = te.rf, observed = te$Labels)
sensitivity(data= te.rf, reference = te$Labels)
specificity(data= te.rf, reference = te$Labels)

#The confusion matrix yields a value of ~0.68 for sensitivity and ~0.6 for specificity
#which indicates that this randomforest is slighty balanced in points of false classification.
#The accuracy of ~0.64 is better than guessing but not too high either.


##########
# Step 3 #
##########
library(boot)
library(pROC)

#The ROC curve is a curve of true-positive/false-positive which is sensitivity/(1-specificity)
#Hence, calculate a confidence interval for the roc

 #ROC
  #calculate the predicted probabilities for the test set
tr.rf.pr <- predict(tr.rf, type="prob", newdata = te[,-1])[,2] #subset [,2] because only the positive prediction are relevant

  #generate prediction object
tr.rf.pred <- prediction(tr.rf.pr, te$Labels)
  #create performance object
tr.rf.perf <- performance(tr.rf.pred, "tpr", "fpr")
  #plot
plot(tr.rf.perf, main = "ROC Curve")
abline(a=0,b=1,lwd=1)

 #obtain confidence interval
#create dataframe for boot function
y <- cbind.data.frame(te$Labels, tr.rf.pr) ; names(y) <- c("teLabels", "tr.rf.pr")


#function to return confidence interval of ROC
boots <- function(d, i){
  data <- d[i,]
  roc.val <- roc(response = data$teLabels, predictor = data$tr.rf.pr, ci=TRUE)
  return(roc.val$ci)
}

#now bootstrap this r = number of observations # = Leave-one-out-CV?!??
n <- nrow(te)
roc.boots<-boot(y, boots, R= n)
roc.boots