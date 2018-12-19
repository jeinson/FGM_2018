## This script is an R version of Logistic Regression with Black Box Variational Inference
## This was adopted from a blog post by Keyon Vafa

library(e1071)
normalize <- function(x){sqrt(sum(x^2))}

logistic_reg_BBVI <- function(X, y, S = 10, n.iter = 10000, convergence = .01){
  # Add and intercept term
  #featNames <- c("Intercept", colnames(X))
  featNames <- colnames(X)
  #X <- cbind(Intecept  = rep(1, nrow(X)), X)
  # Data Parameters
  N <- nrow(X)
  P <- ncol(X)

  elbo.grad <- function(z.sample, mu, sigma, X, y){
    # Variational updates
    score.mu <- (z.sample - mu)/(sigma)
    score.logsigma <- (-1/(2*sigma) + ((z.sample - mu)^2/(2*sigma^2))) * sigma
    # Objective Functions
    upper <- y * log(sigmoid(X %*% z.sample)) # A hacky way of getting around roundoff issues
    upper[!is.finite(upper)] <- 1
    lower <- (1-y) * log(1-sigmoid(X %*% z.sample))
    lower[!is.finite(lower)] <- 0
    log.p <- sum(upper + lower) + sum(log(dnorm(z.sample)))
    log.q <- sum(log(dnorm(z.sample, mu, sqrt(sigma))))
    c(score.mu, score.logsigma) * (log.p - log.q)
  }
  
  # Model adjustment parameters
  mu <- rnorm(P)
  G <- matrix(0, 2*P, 2*P)
  eta <- 1
  log.sigma <- rnorm(P)
  mus <- matrix(0, n.iter, P)
  delta.lambda <- rep(0, n.iter)
  
  
  print("Beginning to Optimize")
  for(i in 1:n.iter){
    mus[i,] <- mu
    
    if(i %% 500 == 0){
      cat(paste("Iteration:",i, "\n"))
      cat(c("Mu:", round(mu, 3), "\n"))
      cat("Sigma:", round(exp(log.sigma), 3), "\n")
    }
    
    sigma <- exp(log.sigma)
    samples <- sapply(1:S, function(s) rnorm(mu, mu, sqrt(sigma)))
    grad.estimate <- rowMeans(apply(samples, 2, elbo.grad, mu, sigma, X, y))
    G <- G + (grad.estimate %o% grad.estimate)
    mu.new <- mu + (eta * 1/sqrt(diag(G)))[1:P] * grad.estimate[1:P]
    log.sigma.new <- log.sigma + (eta * 1/sqrt(diag(G)))[-(1:P)] * grad.estimate[-(1:P)]
    delta.lambda[i] <- normalize(mu.new - mu)
    
    if(normalize(mu.new - mu) < convergence) break;
    
    mu <- mu.new
    log.sigma <- log.sigma.new
  } 
  
  cat("Optimization Complete\n")
  cat("Final mu: ", mu, "\n")
  cat("Final sigma: ", exp(log.sigma), "\n")
  
  names(mu) <- featNames
  names(sigma) <- featNames
  
  return(list(
    mu = mu,
    sigma = sigma, 
    mus = mus, 
    delta.lambda = delta.lambda
  ))
}

plot.mus <- function(model){
  i <- sum(model$delta.lambda > 0)
  mus <- model$mus

  plot(NULL, xlim = c(0, i), ylim = c(min(mus), max(mus)), 
       ylab = expression(mu), xlab = "Iteration")
  colors <- rainbow(ncol(mus))
  colnum <- 0
  apply(mus[1:i,], 2, function(x){
    colnum <<- colnum + 1
    lines(x, type = "l", col = colors[colnum])
  } )
  
}
