#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

require(methods)
library(caret)
library(glmnet)
library(PRROC)

source('SimilarityRegression_HelperFunctions.R')

CurrentFamily <- args[1]
Cores <- as.integer(args[2])
if(Cores > 1){
  print(paste('!Setting up doMC/Paralell (#cores: ', Cores, ')', sep = ''))
  library(doMC)
  registerDoMC(Cores)
}

dir.create(paste('DNA/ByFamily', CurrentFamily, 'Models', sep = '/'), showWarnings = FALSE, recursive = TRUE)

CurrentAlpha <- 0 #Ridge regression 

print('1) Loading Data')
DataFolder <- paste('DNA/ByFamily', CurrentFamily, 'TrainingData/', sep = '/')
# Read & parse Y Data
Y <- read.csv(paste(DataFolder, 'Y_Sims_PctID.csv.gz', sep =''))

#Possible EScoreOverlap Transforms
#Y$plogis <- plogis(Y$EScoreOverlap)
#Y$qlogis <- qlogis(Y$EScoreOverlap)
#Y$qlogis[Y$qlogis == -Inf] <- qlogis(1/32000)
#Y$qlogis[Y$qlogis == Inf] <- qlogis(1- 1/32000)

#Weight positive samples 1/freq (or 10x whichever is higher) higher 
#than negatives because they're usually at a much lower frequency
weightmultiplier <- length(Y$EClass)/sum(Y$EClass)
if(weightmultiplier < 10){weightmultiplier <- 10}
w <- Y$EClass*weightmultiplier
w[w == 0] <- 1

#Make the EClass a factor for glmnet (Y/N > EScoreThreshold) (Highly Similar)
Y$EClass <- as.factor(Y$EClass)
levels(Y$EClass) <- c('N', 'Y')
#Specify True Negatives (TN) based on E-Score Overlap less than 0.2 (Amb vs. Dissimilar)
Y$TN <- as.numeric(Y$EScoreOverlap >= 0.2)

# Read X (predictor) Data
X_Data = list()
X_Data[['PctID']] <- read.csv(paste(DataFolder, 'X_PctID.csv.gz', sep =''), header = T)
X_Data[['PctID (Smooth3)']] <- read.csv(paste(DataFolder, 'X_PctID_Smooth3.csv.gz', sep =''), header = T)
X_Data[['AvgB62']] <- read.csv(paste(DataFolder, 'X_AvgB62.csv.gz', sep =''), header = T)
X_Data[['AvgB62 (Smooth3)']] <-read.csv(paste(DataFolder, 'X_AvgB62_Smooth3.csv.gz', sep =''), header = T)
X_Data <- lapply(X_Data, convertToMat, 3 ) #Skip first 2 columns (Start at 3rd p1)
#X_nearZeroVar <- lapply(X_Data, nearZeroVar, saveMetrics= TRUE)

#### Scale X (predictor) Data & Output info ####
X_Data <- lapply(X_Data, scaleMat ) #Scale data
count <- 0
for(predictor in names(X_Data)){
  count <- count + 1
  u <- as.numeric(attr(X_Data[[predictor]], 'scaled:center'))
  sd <- as.numeric(attr(X_Data[[predictor]], 'scaled:scale'))
  
  if(count == 1){
    Xscaling <- t(data.frame(c(predictor, 'mean', u)))
    colnames(Xscaling) <- c('X', 'Stat', colnames(X_Data[[1]]))
    row.names(Xscaling) <- NULL
  }else{
    Xscaling <- rbind(Xscaling, c(predictor, 'mean', u))
  }
  Xscaling <- rbind(Xscaling, c(predictor, 'sd', sd))
}
write.csv(Xscaling, paste('DNA/ByFamily', CurrentFamily, 'Models/Xscales.csv', sep = '/'),)

#Parse Training/Testing Folds
TestInds <- readLines(paste(DataFolder, 'CVTestIndicies_i0.txt', sep =''))
TestInds <- sapply(TestInds,strsplit,'\t')
TrainInds <- list()
for(l in TestInds){
  hname <- l[1]
  holdouts <- as.integer(l[2:length(l)])
  holdouts <- holdouts + 1 #because r is in 1-based and python is 0 
  TrainInds[[hname]] <- setdiff(seq.int(1,nrow(Y)), holdouts)
}
######## REGRESSION ############
print('2) Building Regression Models')
print('2a) Specify Caret fit control')
fitControl <- trainControl(
  method = "cv",
  #Use TrainInds (Leave-One-Construct-Out CV)
  index = TrainInds,
  verboseIter = T,
  selectionFunction = 'oneSE',
  savePredictions = 'final'
)
print('2b) Specify Training Function')
caretGLMNET.regression <- function(cx){
  #Find the lambda sequence to CV over & Make Caret CV grid 
  if(length(CurrentAlpha) == 1){
    lambdaseq <- glmnet(x = cx, y = Y$EScoreOverlap, weights= w, alpha = CurrentAlpha, lower.limits = 0, standardize = FALSE)$lambda
    glmnetGrid <- expand.grid(alpha=c(CurrentAlpha), lambda=lambdaseq)
  }else{
    lambdaseq <- glmnet(x = cx, y = Y$EScoreOverlap, weights= w, alpha = min(CurrentAlpha), lower.limits = 0, standardize = FALSE)$lambda
    lambdaseq <- c(lambdaseq, glmnet(x = cx, y = Y$EScoreOverlap, weights= w, alpha = max(CurrentAlpha), lower.limits = 0, standardize = FALSE)$lambda)
    lambdaseq <- c(lambdaseq, glmnet(x = cx, y = Y$EScoreOverlap, weights= w, alpha = median(CurrentAlpha), lower.limits = 0, standardize = FALSE)$lambda)
    subsample <- c()
    for(q in seq.int(1, 0, -0.01)){
      subsample <- c(subsample, as.numeric(quantile(lambdaseq,q)))
    }
    glmnetGrid <- expand.grid(alpha=CurrentAlpha, lambda=subsample)  
  }
  #Train Caret Model
  currentModel <- train(x = cx,
                        y = Y$EScoreOverlap, 
                        method = "glmnet",
                        standardize = FALSE,
                        tuneGrid=glmnetGrid,
                        trControl = fitControl,
                        lower.limits = 0,
                        weights= w
  )
  return(currentModel)
}
print('3c) Run Regression CVs')
cvresults.regression <- lapply(X_Data, caretGLMNET.regression)

######## Classification ############
print('3) Building Classification Models')
print('3a) Specify Caret fit control')
fitControl <- trainControl(## LOCO
  method = "cv",
  index = TrainInds,
  verboseIter = T,
  classProbs = T,
  summaryFunction = mnLogLoss,
  selectionFunction = 'oneSE',
  savePredictions = 'final'
)
print('3b) Specify Training Function')
caretGLMNET.logistic <- function(cx){
  if(length(CurrentAlpha) == 1){
    lambdaseq <- glmnet(x = cx, y = Y$EClass,family = 'binomial', weights= w, alpha = CurrentAlpha, lower.limits = 0, standardize = FALSE)$lambda
    glmnetGrid <- expand.grid(alpha=c(CurrentAlpha), lambda=lambdaseq)
  }else{
    lambdaseq <- glmnet(x = cx, y = Y$EClass,family = 'binomial', weights= w, alpha = min(CurrentAlpha), lower.limits = 0,standardize = FALSE)$lambda
    lambdaseq <- c(lambdaseq, glmnet(x = cx, y = Y$EClass,family = 'binomial', weights= w, alpha = max(CurrentAlpha), lower.limits = 0, standardize = FALSE)$lambda)
    lambdaseq <- c(lambdaseq, glmnet(x = cx, y = Y$EClass,family = 'binomial', weights= w, alpha = median(CurrentAlpha), lower.limits = 0, standardize = FALSE)$lambda)
    subsample <- c()
    for(q in seq.int(1, 0, -0.01)){
      subsample <- c(subsample, as.numeric(quantile(lambdaseq,q)))
    }
    glmnetGrid <- expand.grid(alpha=CurrentAlpha, lambda=subsample)  
  }
  #Train Caret Model
  currentModel <- train(x = cx,
                        y = Y$EClass,
                        metric = 'logLoss',
                        method = "glmnet",
                        family = 'binomial',
                        standardize = FALSE,
                        tuneGrid=glmnetGrid,
                        trControl = fitControl,
                        lower.limits = 0,
                        weights = w
  )
  return(currentModel)
}
print('3c) Run Classifier CVs')
cvresults.logistic <- lapply(X_Data, caretGLMNET.logistic)

######## Performance ############
#Join Held-out predicitons from different methods
print('4) Evaluating Performance')
#Initialize heldoutpreds with reggression results
name <- names(cvresults.regression)[1]
modelname <- paste(name, 'Regression', sep = '.')
methods <- c(modelname)
heldoutpreds <- cvresults.regression[[name]]$pred[,c('obs','pred', 'rowIndex', 'Resample')]
colnames(heldoutpreds)[colnames(heldoutpreds)== 'pred'] <- modelname
colnames(heldoutpreds)[colnames(heldoutpreds)== 'obs'] <- 'EScoreOverlap'
for(name in names(cvresults.regression)[2:length(cvresults.regression)]){
  modelname <- paste(name, 'Regression', sep = '.')
  methods <- c(methods, modelname)
  nextpreds <- cvresults.regression[[name]]$pred[,c('pred', 'rowIndex', 'Resample')]
  colnames(nextpreds)[colnames(nextpreds)== 'pred'] <- modelname
  heldoutpreds <- merge(heldoutpreds, nextpreds, by = c('rowIndex', 'Resample'))
}

#Merge in Class
name <- names(cvresults.logistic)[1]
modelname <- paste(name, 'Logistic', sep = '.')
methods <- c(methods, modelname)
nextpreds <- cvresults.logistic[[name]]$pred[,c('obs','Y', 'rowIndex', 'Resample')]
colnames(nextpreds)[colnames(nextpreds)== 'Y'] <- modelname
colnames(nextpreds)[colnames(nextpreds)== 'obs'] <- 'EClass'
nextpreds$EClass <- as.numeric(nextpreds$EClass)-1
heldoutpreds <- merge(heldoutpreds, nextpreds, by = c('rowIndex', 'Resample'))
for(name in names(cvresults.logistic)[2:length(cvresults.logistic)]){
  modelname <- paste(name, 'Logistic', sep = '.')
  methods <- c(methods, modelname)
  nextpreds <- cvresults.logistic[[name]]$pred[,c('Y', 'rowIndex', 'Resample')]
  colnames(nextpreds)[colnames(nextpreds)== 'Y'] <- modelname
  heldoutpreds <- merge(heldoutpreds, nextpreds, by = c('rowIndex', 'Resample'))
}
#Add some extra info 
heldoutpreds$ArrayLenDifference <- Y$ArrayLenDifference[heldoutpreds$rowIndex]
heldoutpreds$MultiAlnFlag <- Y$MultiAlnFlag[heldoutpreds$rowIndex]
heldoutpreds$TN <- Y$TN[heldoutpreds$rowIndex]
write.csv(heldoutpreds, paste('DNA/ByFamily', CurrentFamily, 'Models/Predictions_TestSet.csv', sep = '/'), row.names = FALSE)

#Plot Regression PR Curve
wpr <-pr.curve(scores.class0 = heldoutpreds[,methods[1]], weights.class0 = heldoutpreds$EClass, curve = TRUE, 
               max.compute = T, min.compute = T, rand.compute = T)
AUPRs <- c(wpr$auc.davis.goadrich)
AUPRs_Labels = c(paste(methods[1], ' (', round(wpr$auc.davis.goadrich, 3), ')', sep = ''))
count = 0
for(PRthresh in seq(0.05,0.95,0.05)){
  count = count + 1
  c2mr <- curve2maxRecall(wpr$curve, PRthresh)
  if(count == 1){
    PRThresholdData <- data.frame(Model= methods[1], Precision_TEST= PRthresh, 
                                  Recall_TEST = c2mr[1], Threshold = c2mr[2])
    PRThresholdData$Model <- as.character(PRThresholdData$Model)
  } else{ PRThresholdData <- rbind(PRThresholdData, c(methods[1], PRthresh, c2mr[1], c2mr[2]))}
}
pdf(paste('DNA/ByFamily', CurrentFamily, 'Models/Regression.PRcurve.pdf', sep = '/'))
plot(wpr, main = CurrentFamily, auc.main=F, max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE,  fill.area = F, color = add.alpha(cscale[1], 1))
for(i in 2:4){
  wpr2 <- pr.curve(scores.class0 =  heldoutpreds[,methods[i]], weights.class0 = heldoutpreds$EClass, curve = TRUE, 
                   max.compute = T, min.compute = T, rand.compute = T)
  AUPRs <- c(AUPRs, wpr2$auc.davis.goadrich)
  AUPRs_Labels = c(AUPRs_Labels, paste(methods[i], ' (', round(wpr2$auc.davis.goadrich, 3), ')', sep = ''))
  plot(wpr2, add = TRUE, color = add.alpha(cscale[i], 1))
  for(PRthresh in seq(0.05,0.95,0.05)){
    c2mr <- curve2maxRecall(wpr2$curve, PRthresh)
    PRThresholdData <- rbind(PRThresholdData, c(methods[i], PRthresh, c2mr[1], c2mr[2]))
  }
}
bs <- c('black', 'grey')
b = 0
for(method in c('PctID_L', 'PctID_S')){
  b <- b + 1
  wpr2 <- pr.curve(scores.class0 =  Y[,method], weights.class0 = as.numeric(Y$EClass)-1, curve = TRUE, 
                   max.compute = T, min.compute = T, rand.compute = T)
  AUPRs <- c(AUPRs, wpr2$auc.davis.goadrich)
  AUPRs_Labels = c(AUPRs_Labels, paste(method, ' (', round(wpr2$auc.davis.goadrich, 3), ')', sep = ''))
  plot(wpr2, add = TRUE, color = bs[b])
  for(PRthresh in seq(0.05,0.95,0.05)){
    c2mr <- curve2maxRecall(wpr2$curve, PRthresh)
    PRThresholdData <- rbind(PRThresholdData, c(method, PRthresh, c2mr[1], c2mr[2]))
  }
}
legend('bottomleft', inset = 0.2, AUPRs_Labels, lty=c(1,1), lwd=c(2,2), col = c(cscale[1:4], bs), cex = 0.6)
dev.off()

#Plot Classifier PR Curve
pdf(paste('DNA/ByFamily', CurrentFamily, 'Models/Classification.PRcurve.pdf', sep = '/'))
for(i in 5:8){
  wpr2 <- pr.curve(scores.class0 =  heldoutpreds[,methods[i]], weights.class0 = heldoutpreds$EClass, curve = TRUE, 
                   max.compute = T, min.compute = T, rand.compute = T)
  AUPRs <- c(AUPRs, wpr2$auc.davis.goadrich)
  if(i == 5){
    AUPRs_Labels = c(paste(methods[i], ' (', round(wpr2$auc.davis.goadrich, 3), ')', sep = ''))
    plot(wpr2, main = CurrentFamily, auc.main=F, max.plot = TRUE, min.plot = TRUE, rand.plot = TRUE,  fill.area = F, color = add.alpha(cscale[i-4], 1))
  }else{
    AUPRs_Labels = c(AUPRs_Labels, paste(methods[i], ' (', round(wpr2$auc.davis.goadrich, 3), ')', sep = ''))
    plot(wpr2, add = TRUE, color = add.alpha(cscale[i-4], 1))
  }
  for(PRthresh in seq(0.05,0.95,0.05)){
    c2mr <- curve2maxRecall(wpr2$curve, PRthresh)
    PRThresholdData <- rbind(PRThresholdData, c(methods[i], PRthresh, c2mr[1], c2mr[2]))
  }
}
b = 0
for(method in c('PctID_L', 'PctID_S')){
  b <- b + 1
  wpr2 <- pr.curve(scores.class0 =  Y[,method], weights.class0 = as.numeric(Y$EClass)-1, curve = TRUE, 
                   max.compute = T, min.compute = T, rand.compute = T)
  AUPRs_Labels = c(AUPRs_Labels, paste(method, ' (', round(wpr2$auc.davis.goadrich, 3), ')', sep = ''))
  plot(wpr2, add = TRUE, color = bs[b])
}
legend('bottomleft', inset = 0.2, AUPRs_Labels, lty=c(1,1), lwd=c(2,2), col = c(cscale[1:4], bs), cex = 0.6)
dev.off()

#### NPV/NegativeRecall Data ####
count <- 0
for(method in methods){
  c.npvcurve <- NPVcurve(heldoutpreds$TN, heldoutpreds[,method])
  for(PRthresh in c(seq(.75,0.89,0.05), seq(0.9,1,0.01))){
    count = count + 1
    #print(paste(method, PRthresh, count))
    npvinfo <- c(method, PRthresh)
    #Get LOOCV metrics
    r <- selectNPVCurveCutoff(c.npvcurve, PRthresh)
    npvinfo <- c(npvinfo, as.numeric(r[c('NPV', 'NegativeRecall', 'Cutoff')]))
    #Get FinalModel Stats
    if(count == 1){
      #Manually specify DF
      NPVThresholdData <- data.frame(Model = method,
                                        NPV_Target = PRthresh,
                                        NPV_TEST = r$NPV,
                                        NegativeRecall_TEST = r$NegativeRecall,
                                        Threshold = r$Cutoff)
      NPVThresholdData$Model <- as.character(NPVThresholdData$Model)
    }else{NPVThresholdData <- rbind(NPVThresholdData, npvinfo)}
  }
}
for(method in c('PctID_L', 'PctID_S')){
  c.npvcurve <- NPVcurve(Y$TN, Y[,method])
  for(PRthresh in c(seq(.75,0.89,0.05), seq(0.9,1,0.01))){
    #print(paste(method, PRthresh, count))
    npvinfo <- c(method, PRthresh)
    #Get LOOCV metrics
    r <- selectNPVCurveCutoff(c.npvcurve, PRthresh)
    npvinfo <- c(npvinfo, as.numeric(r[c('NPV', 'NegativeRecall', 'Cutoff')]))
    #Get FinalModel Stats
    NPVThresholdData <- rbind(NPVThresholdData, npvinfo)
  }
}
######### Final Model ###########
#### Predictions ####
print('5) Outputting ModelCoefficents & Final Model Performance Metrics')
FinalPreds <- as.data.frame(as.numeric(Y$EClass)-1)
colnames(FinalPreds) <- 'Class'
for(predictor in names(X_Data)){
  #Regression 
  modelname <- paste(predictor, 'Regression', sep = '.')
  p <- predict(cvresults.regression[[predictor]], X_Data[[predictor]])
  FinalPreds <- cbind(FinalPreds, p)
  colnames(FinalPreds)[colnames(FinalPreds) == 'p'] <- modelname
  #Classification 
  modelname <- paste(predictor, 'Logistic', sep = '.')
  p <- predict(cvresults.logistic[[predictor]], X_Data[[predictor]],type = 'prob')$Y
  FinalPreds <- cbind(FinalPreds, p)
  colnames(FinalPreds)[colnames(FinalPreds) == 'p'] <- modelname
}
FinalPreds <- cbind(FinalPreds, Y[,c('PctID_L', 'PctID_S', 'ArrayLenDifference', 'MultiAlnFlag', 'TN')])
write.csv(FinalPreds, paste('DNA/ByFamily', CurrentFamily, 'Models/Predictions_FinalModel.csv', sep = '/'), row.names = F)

## Positives ##
PRThresholdData[,'Precison_FINAL'] <- NA
PRThresholdData[,'Recall_FINAL'] <- NA
for(i in 1:nrow(PRThresholdData)){
  row <- PRThresholdData[i,]
  if(is.na(row$Threshold) == FALSE){
    curmodel <- row$Model
    curpreds <- FinalPreds[,c('Class', 'TN', curmodel)]
    curpreds <- curpreds[order(curpreds[,curmodel], decreasing = T),]
    curpreds[,curmodel] <- as.numeric(curpreds[,curmodel] >= row$Threshold)
    x <- table(TRUTH = curpreds[,'Class'], PREDICTED = curpreds[,curmodel])
    tryCatch({
      PRThresholdData[i,'Precison_FINAL'] <- x['1','1']/(x['1','1'] + x['0','1'])
      PRThresholdData[i,'Recall_FINAL'] <- x['1','1']/(x['1','1'] + x['1','0'])
    }, error= function(cond){})
  }
}
write.csv(PRThresholdData, paste('DNA/ByFamily', CurrentFamily, 'Models/PRThresholdData.csv', sep = '/'), row.names = F)

## Negatives ##
NPVThresholdData[,'NPV_FINAL'] <- NA
NPVThresholdData[,'NegativeRecall_FINAL'] <- NA
for(i in 1:nrow(NPVThresholdData)){
  row <- NPVThresholdData[i,]
  if(is.na(row$Threshold) == FALSE){
    curmodel <- row$Model
    curpreds <- FinalPreds[,c('TN', curmodel)]
    curpreds[,curmodel] <- as.numeric(curpreds[,curmodel] >= row$Threshold)
    x <- table(TRUTH = curpreds[,'TN'], PREDICTED = curpreds[,curmodel])
    #print(x)
    tryCatch({
      NPVThresholdData[i,'NPV_FINAL'] <- x['0','0']/(x['0','0'] + x['1','0'])
      NPVThresholdData[i,'NegativeRecall_FINAL'] <- x['0','0']/(x['0','0'] + x['0','1'])
    }, error= function(cond){})
  }
}
write.csv(NPVThresholdData, paste('DNA/ByFamily', CurrentFamily, 'Models/NPVThresholdData.csv', sep = '/'), row.names = F)

#### Coefficents ####
count <- 0
for(predictor in names(X_Data)){
  count <- count + 1
  #Regression 
  modelname_R <- paste(predictor, 'Regression', sep = '.')
  coef_R <- coef(cvresults.regression[[predictor]]$finalModel, cvresults.regression[[predictor]]$bestTune$lambda)
  #Classification 
  modelname_L <- paste(predictor, 'Logistic', sep = '.')
  coef_L <- coef(cvresults.logistic[[predictor]]$finalModel, cvresults.logistic[[predictor]]$bestTune$lambda)
  cp <- cbind(coef_R, coef_L)
  colnames(cp) <- c(modelname_R, modelname_L)
  if(count == 1){
    CoefMat <- cp
  }else{
    CoefMat <- cbind(CoefMat, cp)
  }
}
write.csv(as.matrix(CoefMat), paste('DNA/ByFamily', CurrentFamily, 'Models/ModelCoefficents.csv', sep = '/'))

