#### a different PCA function ####

normalisation <- function (x) # assuming no missing data!
{
  x <- as.vector(x)
  
  col.mean <- mean(x)
  n <- length(x)  
  p <- (sum(x) + 1) / (2 * (n + 1))  # estimate of the underlying allele freq
  
  y <- (x - col.mean) / sqrt (p * (1 - p)) # normalise the column
  
  return(y)
} 


# now a function that returns PCs

pca.detect_tw <- function(data) # data matrix
{
  n <- nrow(data) # number of samples
  p <- ncol(data) #- 1 # number of snps
  x <- data#[, -(p + 1)]
  Sn <- matrix(data = 0, nrow = n, ncol = n) # initialise matrix
  
  z <- apply(x, 2, normalisation) # normalised matrix
  
  rm(x) # we don't need to remember this
  
  Sn <- Sn + (z %*% t(z)) # compute unscaled matrix
  
  rm(z) # we don't need to remember this
  
  Sn <- Sn / p # finish computing Sn
  
  # now we do pca
  
  eig <- eigen(Sn, only.values = TRUE)  # do eigenanalysis
  
  # return the whole eigen object for now
  
  val <- eig$values
  
  return(val) # i only need the eigenvalues for tracy widom
  
}



#### tracy widom test ####

tracy_widom <- function(eval){ # eigenvalues sorted
  S1 <- rev(cumsum(rev(eval)))
  S2 <- rev(cumsum(rev(eval^2)))
  n_hat <- (length(eval)):1
  p_hat <- (n_hat + 2)*(S1^2) / (n_hat * S2 - S1^2)
  p_tilde <- S1^2 / S2
  
  n_hat_term <- sqrt(n_hat)
  p_hat_term <- sqrt(p_hat - 1)
  p_tilde_term <- sqrt(p_tilde - 1)
  
  L <- n_hat * eval / S1
  
  mu <- (p_hat_term + n_hat_term)^2 / p_hat
  mu_tilde <- (p_tilde_term + n_hat_term)^2 / p_tilde
  
  sigma <- (p_hat_term + n_hat_term) / p_hat * (1 / p_hat_term + 1 / n_hat_term)^(1/3)
  sigma_tilde <- (p_tilde_term + n_hat_term) / p_tilde * (1 / p_tilde_term + 1 / n_hat_term)^(1/3)
  
  tw <- (L - mu) / sigma
  tw_tilde <- (L - mu_tilde) / sigma_tilde
  
  sign <- sum(tw > 2.0234, na.rm = TRUE)
  sign_tilde <- sum(tw_tilde > 2.0234, na.rm = TRUE)
  
  return(list(price = sign, other = sign_tilde))
  
}


