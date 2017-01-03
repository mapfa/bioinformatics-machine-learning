require(Rcpp)
require(caret)
require(randomForest)

#=========
# Step 1
# ========
#----------->set you working direction(path to folder with tumor-data) here<----------
# setwd()
data <- read.csv("tumor/Patient1.txt", header=F)
labels <- read.csv("tumor/Patient1Label.txt", head=F)

# labels as factors and combine
labels      <- as.factor(unlist(labels))
data$labels <- labels



#=========
# Step 2
# ========

# Create loop to create training data partitions of 3 to 100 samples
# We need positive and negative labels in every single on of our trainsets.
# As the probability for a negative label is 2 to 3 times higher than for a
# positive label, it is hard to get mixed-label-trainsets for small setlengths.

positiveSub <- subset(data, data$labels == 1)
negativeSub <- subset(data, data$labels == 0)
tr <- list()
for (i in 1:98){
    tr[[i]] <- as.data.frame(positiveSub[sample(nrow(positiveSub), 1),])
    tr[[i]] <- rbind(tr[[i]], negativeSub[sample(nrow(negativeSub), size = 1),])
    tr[[i]] <- rbind(tr[[i]], data[sample(nrow(data), size = i),])
}

# Loop to subset the corresponding test set
te <- list()
for (j in 1:98){
  rows <- unique(unlist(mapply(function(x, y) 
          sapply(setdiff(x, y), function(d) which(x==d)), data, tr[[j]])))
  te[[j]] <- data[rows,]
}



# creating lists to save results
errorRate   <- list() # the test error
tr.error    <- list() # the training error
predictions <- list() # the predictions


# build 98 forests based on the different training sets
for (i in 1:length(tr)){
  # train model
  forest <- randomForest(labels ~ ., data = tr[[i]])
  
  # get training error from confusion matrix
  tr.error[[i]] <- sum(forest$confusion[,3], na.rm=T)
  
  # predict
  prediction <- predict(forest, newdata = te[[i]][,-119])
  predictions[[i]] <- prediction
  result <- table(prediction, te[[i]][,"labels"])

  # calculating the test errorrate: (false positives + false negatives) / P + N
  errorRate[[i]] <- sum(result[2,1], result[1,2]) / sum(result)
  }

#===============
# Visualization
#===============

plot(unlist(tr.error), unlist(errorRate), 
     xlab="Training Error", ylab = "Test Error", 
     main = "Correlation of Training and Test Error")
abline(lm(unlist(tr.error) ~ unlist(errorRate), data = data.frame(unlist(tr.error), unlist(errorRate))))

cor.test(unlist(tr.error), unlist(errorRate))
# Note that the correlation between test and training error is not significant
# as the p-value is > 0.05 so that the null hypothesis
# that the correlation between training and test error is zero cannot be rejected.


#=========
# Step 3
# ========

# As in Step 2, we pick a sample from each class for the training set
# and subset the corresponding test set
tr.init <- rbind(positiveSub[sample(nrow(positiveSub),1),], negativeSub[sample(nrow(negativeSub), size = 1),])
te.init <- data[setdiff(rownames(data),rownames(tr.init)),]

# some parameters for the loop
maxIterations <- 50
tr.3 <- tr.init
te.3 <- te.init

saved.queries <- vector(mode="list", maxIterations) # oder als ein Dataframe?
tr.error.3    <- list()
errorRate.3   <- list()
predictions.3 <- list()

for (i in 1:maxIterations){
  # train model
  forest <- randomForest(labels~., data=tr.3)
  
  # get training error
  tr.error.3 <- sum(forest$confusion[,3], na.rm=T)
  
  # prediction with prob. for most uncertain query
  # (i.e. with probability closest to 0.5)
  pred.prob <- predict(forest, newdata=te.3, type="prob")
  query <- data[attr(which.min(apply(pred.prob, 1, max) -0.5), 'names'),]
  saved.queries[[i]] <- query

  # repeat prediction to obtain error rate
  pred.response <- predict(forest, newdata=te.3[,!(colnames(te.3) %in% "labels")], type="response")
  predictions.3[[i]] <- pred.response
  result.3 <- table(pred.response, te.3[,"labels"])
  
  # calculating the errorrate: (false positives + false negatives) / P + N
  errorRate.3[[i]] <- sum(result.3[2,1], result.3[1,2]) / sum(result.3)
  
  tr.3 <- rbind(tr.3, query)
  te.3 <- te.3[setdiff(rownames(te.3), rownames(tr.3)),]
  
}

#================
# Visualiziation 
#================

# The development of training and test error. 
err.dev <- ts(data = cbind(unlist(errorRate.3), unlist(tr.error.3)), names = c("Testerror", "Trainerror"),
              start=1, end=maxIterations)
plot(err.dev)

# The spatial positions of the queries

spatPositions <- read.csv("tumor/patient1_spatial_coordinates.csv", head=F, col.names = c("yPos", "xPos"))
for (i in 1:length(saved.queries)){
  spatPositions[row.names(saved.queries[[i]]), "obsTime"] <- i
}
ggplot(data = spatPositions, aes(x = xPos, y = yPos)) +
  geom_point(aes(colour = obsTime)) +
  scale_colour_gradientn(colours = topo.colors(2))



#=========
# Step 4
# ========

# include spatial positions
data.4 <- cbind(data, spatPositions[,1:2])


# as before, subset labels and sample
positiveSub.4 <- subset(data.4, data.4$labels == 1)
negativeSub.4 <- subset(data.4, data.4$labels == 0)

tr.init.4 <- rbind(positiveSub.4[sample(nrow(positiveSub.4),1),], negativeSub.4[sample(nrow(negativeSub.4), size = 1),])
te.init.4 <- data.4[setdiff(rownames(data.4), rownames(tr.init.4)),]


# create matrix of euclidian distances between each point
p <- as.matrix(dist(cbind(spatPositions[,1], spatPositions[,2])))


# parameters for training
maxIterations.4 <- 50 
tr.4            <- tr.init.4
te.4            <- te.init.4
saved.queries.4 <- vector(mode="list", maxIterations.4) 
predictions.4   <- list()
tr.error.4      <- list()
errorRate.4     <- list()

for (i in 1:maxIterations.4){
  # train model
  forest.4 <- randomForest(labels~., data=tr.4)
  
  # get training error
  tr.error.4 <- sum(forest.4$confusion[,3], na.rm=T)
  
  # prediction with prob. for most uncertain query
  pred.prob.4 <- predict(forest.4, newdata=te.4[,!colnames(te.4) %in% "labels"], type="prob")
  
  # We want that the new query is sufficiently far from all the other 
  # values in the training set. To this end,
  # we first subset pred.prob.4 to only include candidates that ARE suff. far.
  # We then scan it for the most uncertain value (closest to 0.5)
  # and set it as the new query

  t <- pred.prob.4[intersect(rownames(pred.prob.4), 
                             rownames(data.4[apply(p[rownames(te.4),rownames(tr.4)], 1,
                                                   function(x) all(x > 2)), ])),]
                                                   # look at the p matrix and look for the remaining 
                                                   # test points which are far from the training data
                                                   
  query.4 <- data.4[attr(which.min(apply(t,1,max)-0.5), 'names'),]
  
  # save final query
  saved.queries.4[[i]] <- query.4
  
  # repeat prediction to obtain error rate
  pred.response.4    <- predict(forest.4, newdata=te.4[,!colnames(te.4) %in% "labels"], type="response")
  predictions.4[[i]] <- pred.response.4
  result.4           <- table(pred.response.4, te.4[,"labels"])

  # calculating the errorrate: (false positives + false negatives) / P + N
  errorRate.4[[i]] <- sum(result.4[2,1], result.4[1,2]) / sum(result.4)
  
  tr.4 <- rbind(tr.4, query.4)
  te.4 <- te.4[setdiff(rownames(te.4), rownames(tr.4)),]
  
}



#============
# Step 5
#============

#Visually represent your results. Comment on the following aspects: Starting
#from which number of data points are you more successful with your active
#learning approach? Starting from which number is the difference negligible?
model.evaluation <- ts(data = cbind(unlist(errorRate.4), unlist(errorRate[1:50])), names = c("Active Learning", "Randomforest"),
              start=1, end=maxIterations)
plot(model.evaluation, col=c("blue","red"), plot.type = "single")
legend(x="topright", legend = c("Active Learning", "Randomforest"), col=c("blue","red"), lty = 1)
