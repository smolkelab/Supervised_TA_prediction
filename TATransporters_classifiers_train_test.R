#Load required libraries
library(dplyr)
library(caret)
library(DMwR)
library(purrr)
library(pROC)
library(ggplot2)
library(lattice)
library(gridExtra)

#Read in RNA seq .csv data
logFPKMs <- read.csv(file="D:/CME 250/aba.matrix.FPKM.vf.082511.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

#Remove null entries and unnecessary header rows
logFPKMs$X <- NULL
colnames(logFPKMs)[3:13] <- unlist(lapply(logFPKMs[1,c(3:13)], as.character), use.names=FALSE)
logFPKMs <- na.omit(logFPKMs[-1, ])

#Convert transcript abundance data (columns 3 to 13) from characters to numerics
logFPKMs[,c(3:13)] <- apply(logFPKMs[,c(3:13)], 2, function(x) as.numeric(x))

#Remove unnecessary columns
logFPKMs <- logFPKMs[ ,c(1:13)]

#Remove 'phantom' genes with zero expression across all tissues
logFPKMs <- logFPKMs[rowSums(abs(logFPKMs[ ,c(3:13)])) > 0, ]

#Remove entries with unknown functional annotation (too short for sequence analysis)
logFPKMs <- logFPKMs[!(grepl("unknown function", logFPKMs$Unified.Functional.Annotation, ignore.case = TRUE)), ]

#Read in table of known TA gene (locus ID and functional annotation)
known_TA_genes <- read.table(file="D:/CME 250/known_TA_genes.txt", header = FALSE, sep = "\t")
colnames(known_TA_genes) <- c("Locus.ID", "Unified.Functional.Annotation")

#Add output column for training: TA, non-TA, or unknown
logFPKMs$pathway <- "unknown"
logFPKMs$pathway[logFPKMs$Locus.ID %in% known_TA_genes$Locus.ID] <- "TA"
logFPKMs$pathway[!(logFPKMs$pathway == "TA") & !(grepl("ase|transporter|enzyme|P450|CYP|transcription factor|DNA-binding", logFPKMs$Unified.Functional.Annotation, ignore.case = TRUE))] <- "nonTA"

#Remove all samples with "unknown" pathway class - these will be used for application testing afterwards.
modeldata <- logFPKMs[!(logFPKMs$pathway == "unknown"), ]

#Split into training and test data 75:25
set.seed(1)
index <- createDataPartition(modeldata$pathway, p = 0.75, list = FALSE)
train_data <- modeldata[index,c(3:14)]
test_data <- modeldata[-index,c(3:14)]

#Set up dataframe to store testing results
results <- data.frame(Model = c(rep("glm",8), rep("ranger",8), rep("nnet",8)), 
                      Metric = rep(c(rep("Accuracy",4), rep("ROC",4)), 3),
                      Resampling = rep(c("none", "up", "smote", "rose"), 6), 
                      BalAcc = 0, 
                      CompTime = 0)

#Construct list to store confusion matrix for each trained and cross-validated model
ConfMatrices <- vector(mode = "list", length = 3*2*4)

#Set results table index counter
i = 0

for (model_method in c("glm","ranger","nnet")) {
  
  #print(model_method)
  
  for (model_metric in c("Accuracy", "ROC")) {
    
    #print(model_metric)
    
    for (sm in c("none","up","smote","rose")) {
      
      start_time <- Sys.time()
      i <- i + 1
      #print(i)
      #print(sm)
      set.seed(1)
      
      #Set training control based on model and summary function
      if (model_metric == "Accuracy") {
        control <- trainControl(method = "cv", number = 10, summaryFunction = defaultSummary, classProbs = TRUE, verboseIter = FALSE)
      }
      else {
        control <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE, verboseIter = FALSE)
      }
      
      #Set re-sampling method
      if (sm != "none") {control$sampling <- sm}
      
      #Train model
      model <- train(pathway ~ ., data = train_data, method = model_method, metric = model_metric, trControl = control)
      
      #Prediction on test set and confusion matrix
      predictions <- predict(model, test_data[,1:11])
      cm <- confusionMatrix(predictions, as.factor(test_data[,12]))
      end_time <- Sys.time()
      
      #Add relevant performance metrics and confusion matrices to results dataframe
      results$BalAcc[i] <- cm$byClass[11]
      results$CompTime[i] <- end_time - start_time
      ConfMatrices[[i]] <- cm$table
    }
  }
}
