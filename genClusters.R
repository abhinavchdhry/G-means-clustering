# k <- Number of clusters
# d <- Dimension of resulting data
# nump <- A single number denoting number of points in each cluster,
#         or a vector of values denoting number of points in each cluster
# m1 <- A (k x d) matrix of mean of the points in each cluster for each dimension
# m2 <- A (k x d) matrix of SD of points in each cluster for each dimension

genClusters <- function(k, d, nump, m1, m2) {
  
  ## Do param validation here
  if (length(nump) < k) {
    if (length(nump) == 1) {
      nump <- rep(nump, k)
    }
    else {
      print(paste("Length of vector nump is not consistent with number of clusters (k)"))
    }
  }
  
  matrixInitialized = FALSE
  for (i in 1:d) {  # Generate a 
    x <- c()
    
    for (j in 1:k) {
      x <- c(x, rnorm(nump[j], m1[j, i], m2[j, i]))
    }
    
    if (matrixInitialized == FALSE) {
      m <- matrix(x, length(x), 1)
      matrixInitialized = TRUE
    }
    else {
      m <- cbind(m, x)
    }
  }
  
  randomizer <- order(runif(sum(nump), 0, 100))
  m <- m[randomizer, ]
  return(m)
  
}

