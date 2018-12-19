# Test BBVI on dummy data
source("Documents/FGM/logistic_BBVI_function.R")
setwd("Documents/FGM/")
N <- 100
P <- 10
X <- matrix(rnorm(N*P), nrow = N, ncol = P)
colnames(X) <- paste0("Feature", 1:P)
z.real <- rnorm(ncol(X))
# z.real <- c(-3,-1,1,3)
y <- sapply(sigmoid(X %*% z.real), function(x) rbinom(1, 1, x))

dummy.BBVI.model <- logistic_reg_BBVI(X, y, S = 30, convergence = .01)
results <- rbind(mu = dummy.BBVI.model$mu, sigma = dummy.BBVI.model$sigma)

pdf(file = "figures/BBVI_test_S30_N100.pdf", height=4)
plot(NULL, 
     #xlim = c(min(results)-.2, max(results)+.2), 
     xlim = c(-1.2, 2),
     ylim = c(0,8), 
     xlab = expression(x), 
     ylab = expression('p(y|x)'), 
     sub = "S = 30, N = 100")
colors <- rainbow(ncol(results))
colnum <- 1
apply(results, 2, function(x){
  plot.range <- seq(x[1]-10*x[2],x[1]+10*x[2],.005)
  lines(plot.range, dnorm(plot.range, x[1], x[2]), col = colors[colnum])
  abline(v = z.real[colnum], col = colors[colnum], lty=1)
  colnum <<- colnum + 1
})

dummy.naive.model <- glm(y ~ . - 1, family = "binomial", data = data.frame(X,y))
dummy.naive.model$coefficients
for(i in 1:P){
  abline(v=dummy.naive.model$coefficients[i], col = colors[i], lty=4)
}
dev.off()



# z.real
# plot.mus(dummy.BBVI.model)
# 
# # 

# 
# 
# plot(dummy.BBVI.model$delta.lambda[dummy.BBVI.model$delta.lambda > 0], type = "l", 
#      ylab = expression(paste(Delta, "(", lambda, ")")))


