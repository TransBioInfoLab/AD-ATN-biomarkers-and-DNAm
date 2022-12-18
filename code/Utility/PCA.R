###### PCA function ######
### code from https://gist.github.com/tiagochst/d3a7b1639acf603916c315d23b1efb3e

plotPCA <- function (pca, dataset, pheno, group_char, ntop){
  
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  # merge pheno info with PCs
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
  d <- d[as.character(pheno$sample) ,]
  
  test <- identical (as.character(pheno$sample), row.names (d))
  
  if (test == TRUE){
    
    plotData <- merge (d, pheno, by.x = "row.names", by.y = "sample")
    
    
    # add lines for 3 SD from mean
    meanPC1 <- mean (plotData$PC1)
    sdPC1   <- sd (plotData$PC1)
    
    meanPC2 <- mean (plotData$PC2)
    sdPC2   <- sd (plotData$PC2)
    
    # add flag for outlier samples
    plotData$outlier <- ifelse ( abs(plotData$PC1) > meanPC1 + 3*sdPC1 | abs(plotData$PC2) > meanPC2 + 3*sdPC2,
                                 1, 0 )
    
    plotData$outlier_name <- ifelse ( abs(plotData$PC1) > meanPC1 + 3*sdPC1 | abs(plotData$PC2) > meanPC2 + 3*sdPC2,
                                      plotData$Row.names, "" )
    
    title <- paste0("dataset = ", dataset, ", top = ", ntop, " probes ")
    subtitle <- paste0(" x: mean +/- 3*sdPC1 = ", round(meanPC1,1), " +/- 3*", round(sdPC1,1) ,
                       "     y: mean +/- 3*sdPC2 = ", round(meanPC2,1), " +/- 3*", round(sdPC2,1))
    
    p <- ggplot(data= plotData, aes_string(x="PC1", y="PC2", color = group_char)) +
      geom_point(size=1) +
      theme_bw() +
      xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
      ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
      ggtitle(title, subtitle = subtitle) +
      geom_hline (yintercept = meanPC2 + 3*sdPC2, linetype = "dashed") +
      geom_hline (yintercept = meanPC2 - 3*sdPC2, linetype = "dashed") +
      
      geom_vline (xintercept = meanPC1 + 3*sdPC1, linetype = "dashed") +
      geom_vline (xintercept = meanPC1 - 3*sdPC1, linetype = "dashed") +
      geom_text_repel (aes(label = outlier_name), show.legend = FALSE)
    
    print (p)
    
    
    return (plotData)
    
  }
}


## this function orders dataset by most variable features first
## exp_mat = expression matrix, with row names = feature ids, columns = sample ids

## returns:
## - a matrix with features sorted by most variable features on top

OrderDataBySd <- function(exp_mat){
  # compute sds for each row
  sds <- matrixStats::rowSds(exp_mat)
  sdsSorted <- order(sds, decreasing = TRUE)
  
  # order by most variable probes on top
  exp_mat[sdsSorted ,]
}