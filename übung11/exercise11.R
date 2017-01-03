#Step 1
setwd("C:/Users/Lara/Dropbox/Uni/AppliedML/Sheet11/a4-data")

# loading the datasets
clinicalSheet <- read.delim("crc_clinical_sheet.txt", header = TRUE)
dataMat <- read.table("crc_mirna_datamatrix_rpm.tsv", header=TRUE, sep="\t")
genes <- dataMat$Gene
dataMat <- as.data.frame(t(dataMat[,-1]))
colnames(dataMat) <- genes
library(stringr)
# string manipulation to get matching patient IDs
row.names(dataMat) <- str_replace_all(substring(row.names(dataMat), 1, 12), "\\.", "-")

# extracting label-vector from clinical sheet and removing obs from data with missing label
label <- clinicalSheet[match(row.names(dataMat), clinicalSheet$patient), "tumor_site"]
dataMat <- dataMat[!is.na(label),]
label <- label[!is.na(label)]


# feature elimination based on correlation and coefficient of variance
variationCoefficient <- apply(dataMat, 2, var)
dataMat <- dataMat[,variationCoefficient > 0.1]
corMat <- cor(dataMat)
#heatmap(corMat)

library(caret)
highCor <- findCorrelation(corMat, cutoff = 0.85)
#heatmap(corMat[-highCor,-highCor])
smallMatrix <- dataMat[,-highCor]

# after removing obs. with missing values, highly correlated (cutoff 0.85) and low 
# variance (cutoff <= 0.1) features we continue with a smallMatrix of 250 obs and 394 variables.

# Step 2
# Input neurons: Number of features (atm 394, depending on variance and corr. cutoff)
# Output neurons: Number of possible tumor_sites: 4 levels (?)
# Hidden layers: ???

# Step 3
# joining datamatrix and labels
finalFrame <- cbind(label, smallMatrix)
colnames(finalFrame)[1] <- "label"

# splitting in train and testset
set.seed(250)
library(nnet)

trainIndex <- sample(x = 1:dim(finalFrame)[1], size = 0.8 * dim(finalFrame)[1])
trainSet <- finalFrame[trainIndex,]
testSet <- finalFrame[-trainIndex,]

# tuning modelparameters with tune() from e1071 package
# or maybe with train() from caret
library(e1071)
trainSet[,1] <- as.numeric(trainSet[,1])-1
alt <- tune(nnet, label~., data=trainSet, 
            ranges=list(decay=2^(-3:1), size=1:4, maxit= seq(100:1100, 500)), 
            tune.control = tune.control(sampling="cross", cross=5), 
            MaxNWts = 200000,
            trace = F)
# best parameters: decay = 0.125, size = 3

# Step 4
prediction1 <- predict(nnet.train$finalModel, testSet[,-ncol(testSet)])
table(testSet[,ncol(testSet)], prediction1)
accuracy <- sum(testSet[,ncol(testSet)] == prediction1) / length(prediction1)

# Step 5

# number of training-iterations
numbIt <- 20

# matrix for saving error-rates
errorMat <- matrix(ncol = 2, nrow = numbIt)

# training-loop
for(i in 1:numbIt){
  
  #to monitor progress
  cat(i,'\n') 
  flush.console()
  
  loopControl <- trainControl( 
    method = "cv" 
    , number = 5
    , verboseIter = F 
    , returnData = T 
    , savePredictions = T 
  ) 
  
  loopGrid <- expand.grid( 
    size = c(1,5,10), 
    decay = 2^(-3:1) 
  ) 
  
  #to ensure same set of random starting weights are used each time
  set.seed(5)
  
  #temporary nnet model
  tempModel <- train(label ~ ., data = trainSet, method = "nnet", linout = T,
                     maxit = i, trace = F, MaxNWts = 10000,
                     trControl = train.control, tuneGrid = tune.grid)
  
  #training and test accuracy
  trainPrediction <- predict(tempModel, trainSet, type = "raw")
  trainAccuracy <- sum(trainPrediction == trainSet$label) / length(trainPrediction)
  testPrediction <- predict(tempModel, testSet, type = "raw")
  testAccuracy <- sum(testPrediction == testSet$label) / length(testPrediction)
  
  #append error after each iteration
  errorMat[i,] <- c(trainAccuracy, testAccuracy)
  
}