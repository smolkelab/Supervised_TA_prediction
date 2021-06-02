#Load required libraries
library(dplyr)
library(caret)
library(DMwR)
library(purrr)
library(pROC)
library(ggplot2)
library(lattice)
library(gridExtra)
library(NeuralNetTools)

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

#Set model and learning parameters for best-performing model: Neural network ("nnet"), ROC metric ("ROC"), random oversampling ("up")
model_method <- "nnet"
model_metric <- "ROC"
sampling_method <- "up"

#Start timer
start_time <- Sys.time()

#Set up training control
set.seed(1)
control <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE, verboseIter = FALSE)

#Set re-sampling method
control$sampling <- sampling_method
      
#Train optimized model
model <- train(pathway ~ ., data = train_data, method = model_method, metric = model_metric, trControl = control)
      
#Prediction on test set and confusion matrix
predictions <- predict(model, test_data[,1:11])
cm <- confusionMatrix(predictions, as.factor(test_data[,12]))

#Stop timer
end_time <- Sys.time()
      
#Print out relevant results
paste("Balanced accuracy = ", cm$byClass[11], sep="")
paste("Computation time = ", end_time - start_time, sep="")
cm$table

#Extract data for real-world prediction: all genes marked as "unknown" from logFPKMs
appldata <- logFPKMs[(logFPKMs$pathway == "unknown"), ]

#Prediction on unknown genes
predictions_unknown <- predict(model, appldata[,3:13])

#Append predictions to unknown dataframe
appldata$prediction <- predictions_unknown

#Print out statistics on new predictions, and output predicted novel TA-related genes
table(appldata$prediction)
new_hits <- appldata[(appldata$prediction == "TA"), ]
write.csv(new_hits,"D:/CME 250/new_TA_hits.csv", row.names = FALSE)

#Visualize architecture of final optimized model (nnet)
model_architecture <- data.frame(Connection = names(coef(model$finalModel)), Weight = coef(model$finalModel))
plotnet(model, cex_val = 0.75, circle_col = "azure3", pos_col = "green", neg_col = "red", bord_col = "black")
