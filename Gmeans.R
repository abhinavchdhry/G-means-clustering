####################################################################
# Implement G-means algorithm                                      #
####################################################################

#install.packages("vegan")
#install.packages("nortest")
#library(vegan)
#library(nortest)

# Find the 1-D projection of a matrix of points X on line V
projection <- function(X, v) {

  if (ncol(X) != length(v)) {
    stop("projection:: Dimensions of points in X is not same as dimension of v")
  }

  v <- matrix(v, ncol=1)
  out.vector <- as.vector(X%*%v/sum(v*v))
  
  return (out.vector)
}

getCluster <- function(i, data, cl) {
  if (nrow(data) != length(cl)) {
    stop(paste("Data and cluster vector are not of same length, data =", nrow(data), " cluster =", length(cl)))
    print(cl)
  }
  
  # Return rows which belong to cluster i
  return(data[cl == i, ])
}

# Compute initial split centers based on PCA, Section 2.1
computeInitialSplitCenters <- function(c, cluster) {
  pr.comps <- prcomp(cluster, scale=TRUE)
  prc1 <- pr.comps$rotation[,1]
  eigenval1 <- vegan::eigenvals(pr.comps)[1]
  
  splitFactor <- sqrt(2*eigenval1 / pi)*prc1
  v <- rbind(as.vector(c + splitFactor), as.vector(c - splitFactor))
  return(v)
}

# X is the input dataset(it is a dataframe contain multiple cols reading from csv file)
Gmeans <- function(X,alpha = 0.0001,k=1){
  
  # Convert the data frame to a matrix for easier handling
  M <- data.matrix(X)
  
  # Run kmeans first to find the initial k centers
  kout <- kmeans(M, k)
  centers <- kout$centers
  clusters <- kout$cluster
  
  GmeansIterationCount <- 0
  newCentersToAddInited = FALSE
  
  # Begin Gmeans iterations
  while (1) {
    ncenters <- dim(centers)[1]

    centersToKeep <- c()   # Indices of centers to keep  after this iteration
    newCentersToAddInited = FALSE

    for (i in 1:ncenters) {
      cl <- getCluster(i, M, clusters)
      
      # Split the center of the current cluster using principal component technique
      splitCenters <- computeInitialSplitCenters(centers[i], cl)

      # XXX BUG: kmeans sometimes throws an Error: empty cluster: try a better set of initial centers
      # In such a case, we do not provide the set of splitCenters to kmeans, but rather
      # let kmeans figure out 2 new centers
      t <- tryCatch(koutnew <- kmeans(cl, splitCenters), error = function(e) e)
      
      if (is(t, "error")) {
        koutnew <- kmeans(cl, 2)
      }
        
      newCenters <- koutnew$centers

      v <- as.vector(newCenters[1,] - newCenters[2,])
      p <- projection(cl, v)
      
      p <- (p-mean(p))/sd(p)  # Normalize
      
      # Anderson-Darling test on p
      adtest <- nortest::ad.test(p)
      
      if (adtest$p.value <= alpha) {   # Reject H0. cluster is not Gaussian. Accept the split
        if (newCentersToAddInited == FALSE) {
          newCentersToAdd <- as.matrix(newCenters)
          newCentersToAddInited = TRUE
        }
        else {
          newCentersToAdd <- rbind(newCentersToAdd, newCenters)
        }
      }
      else {  # Reject H1. Cluster follows a Gaussian dist. Do not accept split
        centersToKeep <- c(centersToKeep, i)
      }
    } # for end

    # Do a sanity check here.
    # The following equality should be valid:
    # (ncenters - length(centersToKeep)) * 2 = nrow(newCentersToAdd)
    if (newCentersToAddInited == TRUE) {
      if ( 2*(ncenters- length(centersToKeep)) != nrow(newCentersToAdd) ) {
        stop(paste("SANITY check failed. ncenters = ", ncenters, "length(centersToKeep) = ", length(centersToKeep), "nrow(newCentersToAdd) =", nrow(newCentersToAdd)))
      }
    }
    
    GmeansIterationCount = GmeansIterationCount + 1

    # If number of centers is same as before, we are done
    if (length(centersToKeep) == ncenters)
      break   # Break while
    
    # Create the new set of centers to use for the next iteration
    centers <- centers[centersToKeep,]
    centers <- rbind(centers, newCentersToAdd)

    # At this point we may have a new set of centers. Run kmeans again
    kout <- kmeans(M, centers)
    centers <- kout$centers
    clusters <- kout$cluster
    
    ncenters <- dim(centers)[1]
  } # while end
  
  # Below is for increased readibility of output
  # TODO: Should probably create a class here before returning
  clnms <- c()
  rnames <- c()
  for (i in 1:dim(centers)[2]) {
    clnms <- c(clnms, paste("dim", i, sep=""))
  }
  for (i in 1:dim(centers)[1]) {
    rnames <- c(rnames, paste("c", i, sep=""))
  }
  centersDF <- data.frame(centers)
  colnames(centersDF) <- clnms
  rownames(centersDF) <- rnames
  
  return (list(centers=centersDF, numIterations=GmeansIterationCount, clusters=clusters))
}

# Assumes data is N x 2 matrix
plot.clusters.2D <- function(data, clusters) {
  if (nrow(data) != length(clusters))
    stop("Data and clusters vector should be of same length")

  plot(data)
  df <- as.data.frame(cbind(data, clusters))
  names(df) <- c("1", "2", "clusters")

  for (i in unique(df$clusters)) {
    cluster <- subset(df, clusters == i, select = c("1", "2"))
    ellipse(mu=colMeans(cluster), sigma=cov(cluster), alpha=0.05, npoints=nrow(cluster), col="red")
  }
}