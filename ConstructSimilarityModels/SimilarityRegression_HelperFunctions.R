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

