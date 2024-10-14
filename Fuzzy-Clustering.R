# Clear all variables from the environment
rm(list = ls(all = TRUE))

# Load required libraries
library(e1071)
library(mclust)
library(fclust)
library(modi)
library(inaparc)
library(bootcluster)

# Custom Manhattan Distance Function
customManhattanDistance <- function(X, centers, k1, k2) {
  if (ncol(X) != ncol(centers)) 
    stop(sQuote("X"), " and ", sQuote("centers"), " must have the same number of columns")
  
  dnew <- matrix(0, nrow = nrow(X), ncol = nrow(centers))
  
  for (k in 1:nrow(centers)) {
    d <- matrix(0, nrow = ncol(X), ncol = nrow(X))
    for (i in 1:nrow(X)) {
      for (j in 1:ncol(X)) {
        # Condition to apply different weights (k1 and k2)
        d[j, i] <- ifelse(X[i, j] > centers[k, j],
                          k1 * (X[i, j] - centers[k, j]),
                          k2 * (centers[k, j] - X[i, j]))
      }
    }
    dnew[, k] <- apply(d, 2, sum)
  }
  return(dnew)
}

# Fuzzy C-Means with custom Manhattan distance
customFuzzyCMeans <- function(X, centers, m, iter.max, iterations, threshold, k1, k2) {
  data.X <- as.matrix(X)
  n <- nrow(data.X)
  p <- ncol(data.X)
  best_func_obj <- Inf
  best_cluster <- best_member <- best_centroid <- NULL
  
  for (iter in 1:iterations) {
    # Initialize centers using k-means++
    V <- inaparc::kmpp(data.X, k = centers)$v
    D <- customManhattanDistance(data.X, V, k1, k2)
    
    U <- matrix(0, nrow = n, ncol = centers)
    for (i in 1:n) {
      U[i, ] <- 1 / (((D[i, ]) ^ (2 / (m - 1))) * sum((1 / D[i, ]) ^ (2 / (m - 1))))
    }
    U[is.nan(U)] <- 1
    U[U < 0] <- 0
    U[U > 1] <- 1
    
    iteration <- 1
    repeat {
      U.old <- U
      V.old <- V
      
      # Update centroids using weighted quantiles
      for (i in 1:centers) {
        for (j in 1:ncol(data.X)) {
          V[i, j] <- weighted.quantile(data.X[, j], U[, i] ^ m, prob = k1 / (k1 + k2), plot = FALSE)
        }
      }
      
      # Update distance matrix
      D <- customManhattanDistance(data.X, V, k1, k2)
      
      # Update membership matrix U
      for (i in 1:n) {
        U[i, ] <- 1 / (((D[i, ]) ^ (2 / (m - 1))) * sum((1 / D[i, ]) ^ (2 / (m - 1))))
      }
      U[is.nan(U)] <- 1
      U[U < 0] <- 0
      U[U > 1] <- 1
      
      # Objective function
      func_obj <- sum(U ^ m * D)
      
      iteration <- iteration + 1
      if ((max(abs(U.old - U)) <= threshold) || (iteration > iter.max)) {
        break
      }
    }
    
    # Store the best clustering result
    label <- apply(U, 1, which.max)
    if (func_obj < best_func_obj) {
      best_func_obj <- func_obj
      best_cluster <- label
      best_member <- U
      best_centroid <- V
    }
  }
  
  return(list(member = best_member, centroid = best_centroid, func_obj = best_func_obj,
              cluster = best_cluster, centers = centers, q = k1 / (k1 + k2), m = m))
}

# Load iris dataset and preprocess
X = iris[,-5] # Exclude label column
label = iris[,5]

# Define center initialization
centers <- inaparc::kmpp(X, 3)$v

# Calculate ARI using Euclidean distance for 1000 iterations
ARI.F.Euclidean <- numeric(1000)
for (i in 1:1000) {
  result.Euclidean <- cmeans(X, iter.max = 100, centers = centers, method = "cmeans", m = 2, dist = "euclidean")
  ARI.F.Euclidean[i] <- ARI.F(as.numeric(label), result.Euclidean$member)
}
ARIF.Euclidean <- max(ARI.F.Euclidean)
FARI.Euclidean <- mean(ARI.F.Euclidean)

# ARI with Manhattan distance
ARI.F.Manhattan <- numeric(1000)
for (i in 1:1000) {
  result.Manhattan <- cmeans(X, iter.max = 100, centers = centers, method = "cmeans", m = 2, dist = "manhattan")
  ARI.F.Manhattan[i] <- ARI.F(as.numeric(label), result.Manhattan$member)
}
ARIF.Manhattan <- max(ARI.F.Manhattan)
FARI.Manhattan <- mean(ARI.F.Manhattan)

# Compute ARI for multiple q values using custom fuzzy C-Means
ARI.F.results <- list()
q_values <- c(0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9)
k1_values <- c(1, 1, 1, 3, 2, 1, 3, 7, 3, 4, 9)
k2_values <- c(9, 4, 3, 7, 3, 1, 2, 3, 1, 1, 1)

for (q in seq_along(q_values)) {
  ARI.F.results[[q]] <- numeric(1000)
  k1 <- k1_values[q]
  k2 <- k2_values[q]
  
  for (i in 1:1000) {
    result <- customFuzzyCMeans(X, centers = 3, m = 2, iter.max = 100, iterations = 1, threshold = 1e-09, k1, k2)
    ARI.F.results[[q]][i] <- ARI.F(as.numeric(label), result$member)
  }
}

# Extract max and mean ARI values for each q
ARIF.q <- sapply(ARI.F.results, max)
FARI.q <- sapply(ARI.F.results, mean)

# Plot results
barplot_values <- c(ARIF.Euclidean, ARIF.Manhattan, ARIF.q)
names <- c("Euclidean", "Manhattan", paste0("q=", q_values))
bar_colors <- c("black", "grey25", rep("white", length(q_values)))

barplot(matrix(barplot_values, nrow = 13, byrow = TRUE), beside = TRUE,  
        main = "FARI of Iris", col = bar_colors, ylim = c(0, 0.95), cex.names = 0.6)
legend("topleft", legend = c("C-means", "C-medians", "C-quantile"), cex = 0.7, fill = bar_colors, bty = "n")