#For Scaling Data 
convertToMat <- function(df, colstart = 1){
  as.matrix(df[,colstart:ncol(df)])
}

scaleMat <- function(mat, median = FALSE){
  if(median == TRUE){
    mat <- scale(mat, center = apply(mat, 2, median))
  }else{mat <- scale(mat)}
  mat[is.na(mat)] <- 0
  return(mat)
}

scaleMAD <- function(mat){
  medians <- apply(mat, 2, median)
  mads <- apply(mat, 2, mad)
  scaled <- mat - medians
  scaled <- scaled/mads
  return(scaled)
}

scaleIQR <- function(mat){
  medians <- apply(mat, 2, median)
  IQRs <- apply(mat, 2, IQR)
  scaled <- mat - medians
  scaled <- scaled/IQRs
  return(scaled)
}

#For Precision-Recall
curve2maxRecall <- function(c, TargetPrecision){
  c <- c[order(c[,1], c[,2], decreasing = T),] #order by max recall, then precision
  maxRecall <- NA
  maxRecall_Cutoff <- NA
  tryCatch({c <- c[c[,2] >= TargetPrecision,]
  maxRecall <- c[1,1]
  maxRecall_Cutoff <- c[1,3]
  }, error= function(cond){
  })
  return(c(maxRecall, maxRecall_Cutoff))
}

#Calculate Negative Predictive Value (NPV) and Neagtive Recall for each threshold
NPVcurve <- function(labels, scores){
  #Join & Sort Preds
  preds <- cbind(labels, scores)
  preds <- preds[order(preds[,2]),]
  #Calculate total # of Negatives
  N <- nrow(preds) - sum(preds[,1])
  #Create Output Data
  uScores <- unique(preds[,2])
  p <- matrix(nrow = length(uScores), ncol = 3)
  colnames(p) <- c('Cutoff', 'NPV', 'NegativeRecall')
  #Loop through scores
  i <- 0
  for(c in uScores){
    i <- i + 1
    p[i,'Cutoff'] <- c
    c.Preds <- preds[preds[,2] <= c,]
    if(is.matrix(c.Preds)){
      c.NPV <- (nrow(c.Preds) - sum(c.Preds[,1]))/nrow(c.Preds)
      p[i, 'NPV'] <- c.NPV
      c.NegativeRecall <- (nrow(c.Preds) - sum(c.Preds[,1]))/N
      p[i, 'NegativeRecall'] <- c.NegativeRecall
    }else if(is.numeric(c.Preds)){
      p[i, 'NPV'] <- 1 - c.Preds['labels']
      p[i, 'NegativeRecall'] <- (1 - c.Preds['labels'])/N
    }
  }
  as.data.frame(p)
}

selectNPVCurveCutoff <- function(npvcurve, npvTarget, positiveThresh = NULL){
  # 1) Subset threshold below positiveThresh (if given)
  if(!is.null(positiveThresh)) npvcurve <- npvcurve[npvcurve$Cutoff < positiveThresh,]
  # 2) Find thresholds above NPV threshold (npvTarget)
  npvcurve <- npvcurve[npvcurve$NPV >= npvTarget,]
  # 3) Return the last row (because it's already sorted)
  if(dim(npvcurve)[1] == 0){
   result <- c(NA, NA, NA)
   names(result) <- names(npvcurve)
  }else{
    result <- npvcurve[nrow(npvcurve),]
  }
  return(as.list(result))
}

#For Plotting PR-Curves
library('RColorBrewer')
cscale = brewer.pal(7, "Set1")
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

