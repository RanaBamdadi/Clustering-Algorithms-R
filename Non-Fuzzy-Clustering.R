# Clear all objects from the environment
rm(list = ls(all = TRUE))

# Load necessary libraries
library(mclust)

# Define custom quantile function
computeQuantile <- function(x) {
  quantile(x, probs = (k1 / (k1 + k2)))
}

# Custom function to compute distances based on asymmetric distance measure
computeCustomDistance <- function(X, centers) {
  if (ncol(X) != ncol(centers)) {
    stop("X and centers must have the same number of columns")
  }
  
  dnew <- matrix(0, nrow = nrow(X), ncol = nrow(centers))
  d <- matrix(0, nrow = ncol(X), ncol = nrow(X))
  
  for (k in 1:nrow(centers)) {
    p <- apply(X, 1, '>', centers[k, ])
    
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        if (p[j, i]) {
          d[j, i] <- k1 * (t(X[i, j]) - centers[k, j])
        } else {
          d[j, i] <- k2 * (centers[k, j] - t(X[i, j]))
        }
      }
    }
    dnew[, k] <- apply(d, 2, sum)
  }
  return(dnew)
}

# Custom K-Quantiles Clustering Algorithm with Asymmetric Distance
customKquantiles <- function(X, k, maxIter, totalIterations, tolerance) {
  data.X <- as.matrix(X)
  
  bestCluster <- NULL
  bestCost <- Inf
  bestCenters <- NULL
  
  for (iteration in 1:totalIterations) {
    centers <- inaparc::kmpp(data.X, k = k)$v
    previousCenters <- centers
    
    for (iter in 1:maxIter) {
      # Compute distance and assign clusters
      distances <- computeCustomDistance(data.X, centers)
      clusters <- apply(distances, 1, which.min)
      
      # Update centroids
      if (length(unique(clusters)) != k) {
        break
      }
      
      for (i in 1:k) {
        clusterData <- data.X[clusters == i, ]
        if (nrow(clusterData) == ncol(data.X)) {
          centers[i, ] <- clusterData
        } else {
          centers[i, ] <- apply(clusterData, 2, computeQuantile)
        }
      }
      
      # Check for convergence
      centerShift <- sum(apply(abs(previousCenters - centers), 2, sum))
      if (centerShift < tolerance) {
        break
      }
      previousCenters <- centers
    }
    
    # Compute cost function
    cost <- 0
    for (i in 1:k) {
      p <- matrix(as.numeric(t(apply(X, 1, '>', centers[i, ]))), nrow = nrow(X), ncol = ncol(X))
      weight <- ifelse(p == 1, k1, -k2)
      cost <- cost + sum(rowSums((weight[clusters == i, ] * (X[clusters == i, ] - centers[i, ]))))
    }
    
    # Keep the best configuration based on cost
    if (cost < bestCost) {
      bestCluster <- clusters
      bestCost <- cost
      bestCenters <- centers
    }
  }
  
  # Return the clustering result
  return(list(
    centroids = bestCenters,
    objectiveFunction = bestCost,
    clusters = bestCluster,
    numClusters = k,
    iterations = totalIterations
  ))
}

# Custom Manhattan Distance
computeManhattanDistance <- function(X, centers) {
  if (ncol(X) != ncol(centers)) {
    stop("X and centers must have the same number of columns")
  }
  d <- matrix(0, nrow = nrow(X), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    d[, k] <- colSums(abs(t(X) - centers[k, ]))
  }
  return(d)
}

# Basic K-Medians Clustering with Manhattan Distance
basicKmedians <- function(X, k, maxIter, totalIterations, tolerance) {
  data.X <- as.matrix(X)
  
  bestCluster <- NULL
  bestCost <- Inf
  bestCenters <- NULL
  
  for (iteration in 1:totalIterations) {
    centers <- inaparc::kmpp(data.X, k = k)$v
    previousCenters <- centers
    
    for (iter in 1:maxIter) {
      # Compute Manhattan distance and assign clusters
      distances <- computeManhattanDistance(data.X, centers)
      clusters <- apply(distances, 1, which.min)
      
      # Update centroids
      if (length(unique(clusters)) != k) {
        break
      }
      
      for (i in 1:k) {
        clusterData <- data.X[clusters == i, ]
        centers[i, ] <- apply(clusterData, 2, median)
      }
      
      # Check for convergence
      centerShift <- sum(apply(abs(previousCenters - centers), 2, sum))
      if (centerShift < tolerance) {
        break
      }
      previousCenters <- centers
    }
    
    # Compute cost function
    cost <- 0
    for (i in 1:k) {
      clusterData <- data.X[clusters == i, ]
      cost <- cost + sum(rowSums(abs(clusterData - centers[i, ])))
    }
    
    # Keep the best configuration based on cost
    if (cost < bestCost) {
      bestCluster <- clusters
      bestCost <- cost
      bestCenters <- centers
    }
  }
  
  # Return the clustering result
  return(list(
    centroids = bestCenters,
    objectiveFunction = bestCost,
    clusters = bestCluster,
    numClusters = k,
    iterations = totalIterations
  ))
}

# Custom Euclidean Distance
computeEuclideanDistance <- function(X, centers) {
  if (ncol(X) != ncol(centers)) {
    stop("X and centers must have the same number of columns")
  }
  d <- matrix(0, nrow = nrow(X), ncol = nrow(centers))
  for (k in 1:nrow(centers)) {
    d[, k] <- sqrt(colSums((t(X) - centers[k, ])^2))
  }
  return(d)
}

# Basic K-Means Clustering with Euclidean Distance
basicKmeans <- function(X, k, maxIter, totalIterations, tolerance) {
  data.X <- as.matrix(X)
  
  bestCluster <- NULL
  bestCost <- Inf
  bestCenters <- NULL
  
  for (iteration in 1:totalIterations) {
    centers <- inaparc::kmpp(data.X, k = k)$v
    previousCenters <- centers
    
    for (iter in 1:maxIter) {
      # Compute Euclidean distance and assign clusters
      distances <- computeEuclideanDistance(data.X, centers)
      clusters <- apply(distances, 1, which.min)
      
      # Update centroids
      if (length(unique(clusters)) != k) {
        break
      }
      
      for (i in 1:k) {
        clusterData <- data.X[clusters == i, ]
        centers[i, ] <- apply(clusterData, 2, mean)
      }
      
      # Check for convergence
      centerShift <- sum(apply(abs(previousCenters - centers), 2, sum))
      if (centerShift < tolerance) {
        break
      }
      previousCenters <- centers
    }
    
    # Compute cost function
    cost <- 0
    for (i in 1:k) {
      clusterData <- data.X[clusters == i, ]
      cost <- cost + sum(rowSums(abs(clusterData - centers[i, ])))
    }
    
    # Keep the best configuration based on cost
    if (cost < bestCost) {
      bestCluster <- clusters
      bestCost <- cost
      bestCenters <- centers
    }
  }
  
  # Return the clustering result
  return(list(
    centroids = bestCenters,
    objectiveFunction = bestCost,
    clusters = bestCluster,
    numClusters = k,
    iterations = totalIterations
  ))
}


# Example Usage with Iris Dataset
X <- iris[, -5]
labels <- iris[, 5]

# Define k1 and k2 values
k1 <- 1
k2 <- 1

# Run the custom K-Quantiles clustering with asymmetric distance
result <- customKquantiles(X, k = 3, maxIter = 100, totalIterations = 10, tolerance = 1e-6)

# Compute misclassification and ARI for q=0.5
adjustedRandIndex(as.numeric(label), result$clusters)


# Custom Manhattan Clustering
basicManhattan <- basicKmedians(X, k = 3, maxIter = 100, totalIterations = 10, tolerance = 1e-6)
adjustedRandIndex(as.numeric(label), basicManhattan$clusters)

# Custom Euclidean Clustering 
basicEuclidean <- basicKmeans(X, k = 3, maxIter = 1000, totalIterations = 10, tolerance = 1e-6)
adjustedRandIndex(as.numeric(label), basicEuclidean$clusters)

