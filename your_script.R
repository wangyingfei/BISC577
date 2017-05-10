#Q3
options(download.file.method="libcurl", url.method="libcurl")
source("https://www.bioconductor.org/biocLite.R")
biocLite("DNAshapeR")
install.packages('caret')
require(DNAshapeR)
require(caret)

#Q4
## Predict DNA shapes
fn_fasta <- "Mad.txt.fa"
fn_fasta <- "Max.txt.fa"
fn_fasta <- "Myc.txt.fa"
pred <- getShape(fn_fasta)
## Encode feature vectors
featureType1 <- c("1-mer")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType1)
featureType2 <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType2)
head(featureVector)
## Build MLR model by using Caret
# Data preparation
fn_exp <- "Mad.txt"
fn_exp <- "Max.txt"
fn_exp <- "Myc.txt"
exp_data <- read.table(fn_exp)
df <- data.frame(affinity=exp_data$V2, featureVector)
# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)
# Prediction without L2-regularized
#model <- train (affinity~ ., data = df, trControl=trainControl, method = "lm", preProcess=NULL)
#summary(model)
# Prediction with L2-regularized
model2 <- train(affinity~., data = df, trControl=trainControl, method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model2

#Q5
## Install and initialize packages
#install.packages("ggplot2")
#install.packages("grid")
library(ggplot2)
library(grid)
## Data preparation
data1 <- c(0.7751104,0.7855506,0.7776267)
data2 <- c(0.8631741,0.8644418,0.8550499)
## Ploting
ggplot() +
  geom_point(aes(x = data1, y = data2), color = "red", size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) 

#Q7
# Extract sample sequences
fn <- "bound_500.fa"
fn <- "unbound_500.fa"
# Predict DNA shapes
pred <- getShape(fn)
# Generate ensemble plots
plotShape(pred$MGW)
heatShape(pred$ProT, 20)
plotShape(pred$Roll)
plotShape(pred$HelT)

#Q8
## Install packages
#install.packages("caret")
install.packages("e1071")
install.packages("ROCR")
biocLite("Biostrings")
## Initialization
library(DNAshapeR)
library(caret)
library(ROCR)
library(Biostrings)
workingPath <- "C:/Users/Wang Yingfei/Desktop/BISC 577/Remo/as2/"
## Generate data for the classifcation (assign Y to bound and N to non-bound)
# bound
boundFasta <- readDNAStringSet(paste0(workingPath, "bound_30.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")
# non-bound
nonboundFasta <- readDNAStringSet(paste0(workingPath, "unbound_30.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")
# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(workingPath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)
## DNAshapeR prediction
pred <- getShape(paste0(workingPath, "ctcf.fa"))
## Encode feature vectors
featureType <- c("1-mer")
#featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)
## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl, method = "glm", family = binomial, metric ="ROC")
summary(model)
## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance)
## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc

