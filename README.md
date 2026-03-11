# MCMC-method-for-MLE-of-zero-truncated-geometric-ishita

# n=100
rm(list=ls())

# Precompute combination
max_x <- 20
comb_list <- lapply(1:max_x,function(x) choose(x,0:x))

# ZTGID pmf

pmf_ztgid <- function(x,theta){
  
  const <- theta^3/(theta^3+2)
  
  k <- 0:x
  
  num <- sum((-1)^k * comb_list[[x]] *
               (theta/(theta+k+1)+2/(theta+k+1)^3))
  
  den <- 1-const*(theta/(theta+1)+2/(theta+1)^3)
  
  const*num/den
}

# Random sample generator

rztgid <- function(n,theta){
  
  xvals <- 1:max_x
  
  probs <- sapply(xvals,pmf_ztgid,theta=theta)
  probs <- probs/sum(probs)
  
  sample(xvals,n,replace=TRUE,prob=probs)
}

#  log likelihood

loglik <- function(theta,x){
  
  if(theta<=0) return(-Inf)
  
  n <- length(x)
  
  term1 <- n*log(theta^3/(theta^3+2))
  
  inner <- numeric(n)
  
  for(i in 1:n){
    
    xi <- x[i]
    k <- 0:xi
    
    inner[i] <- sum((-1)^k * comb_list[[xi]] *
                      (theta/(theta+k+1)+2/(theta+k+1)^3))
  }
  
  term2 <- sum(log(inner))
  
  term3 <- -n*log(1-(theta^3/(theta^3+2))*
                    (theta/(theta+1)+2/(theta+1)^3))
  
  term1+term2+term3
}

#  MCMC

mcmc_theta <- function(x,iter=800){
  
  theta <- numeric(iter)
  theta[1] <- 2
  
  proposal_sd <- 0.2
  
  for(i in 2:iter){
    
    prop <- rnorm(1,theta[i-1],proposal_sd)
    
    logratio <- loglik(prop,x)-loglik(theta[i-1],x)
    
    if(log(runif(1))<logratio)
      theta[i] <- prop
    else
      theta[i] <- theta[i-1]
  }
  
  burn <- 300
  chain <- theta[(burn+1):iter]
  
  c(mean(chain),sd(chain))
}

# 1000 theta values

theta_values <- runif(1000,2,4)

n <- 100

theta_est <- numeric(1000)
se_est <- numeric(1000)

# Simulation

for(j in 1:1000){
  
  theta_true <- theta_values[j]
  
  x <- rztgid(n,theta_true)
  
  res <- mcmc_theta(x)
  
  theta_est[j] <- res[1]
  se_est[j] <- res[2]
  
  if(j%%50==0) cat("Replication:",j,"\n")
}

# Final result

mean_theta <- mean(theta_est)
mean_se <- mean(se_est)

cat("Mean of 1000 Estimated Theta =",mean_theta,"\n")
cat("Mean of 1000 Standard Error =",mean_se,"\n")

hist(theta_est,
     main="Histogram of 1000 Estimated Values of Theta",
     xlab=expression(hat(theta)),
     col="skyblue",
     border="black",
     breaks=30)
abline(v = mean(theta_est), col = "red", lwd = 2)

hist(se_est,
     main="Histogram of 1000 Estimated Standard Errors",
     xlab="Standard Error",
     col="lightgreen",
     border="black",
     breaks=30)
abline(v = mean(se_est), col = "red", lwd = 2)

#n=50
rm(list=ls())

# Precompute combination
max_x <- 20
comb_list <- lapply(1:max_x,function(x) choose(x,0:x))

# ZTGID pmf

pmf_ztgid <- function(x,theta){
  
  const <- theta^3/(theta^3+2)
  
  k <- 0:x
  
  num <- sum((-1)^k * comb_list[[x]] *
               (theta/(theta+k+1)+2/(theta+k+1)^3))
  
  den <- 1-const*(theta/(theta+1)+2/(theta+1)^3)
  
  const*num/den
}

# Random sample generator

rztgid <- function(n,theta){
  
  xvals <- 1:max_x
  
  probs <- sapply(xvals,pmf_ztgid,theta=theta)
  probs <- probs/sum(probs)
  
  sample(xvals,n,replace=TRUE,prob=probs)
}

#  log likelihood

loglik <- function(theta,x){
  
  if(theta<=0) return(-Inf)
  
  n <- length(x)
  
  term1 <- n*log(theta^3/(theta^3+2))
  
  inner <- numeric(n)
  
  for(i in 1:n){
    
    xi <- x[i]
    k <- 0:xi
    
    inner[i] <- sum((-1)^k * comb_list[[xi]] *
                      (theta/(theta+k+1)+2/(theta+k+1)^3))
  }
  
  term2 <- sum(log(inner))
  
  term3 <- -n*log(1-(theta^3/(theta^3+2))*
                    (theta/(theta+1)+2/(theta+1)^3))
  
  term1+term2+term3
}

#  MCMC

mcmc_theta <- function(x,iter=800){
  
  theta <- numeric(iter)
  theta[1] <- 2
  
  proposal_sd <- 0.2
  
  for(i in 2:iter){
    
    prop <- rnorm(1,theta[i-1],proposal_sd)
    
    logratio <- loglik(prop,x)-loglik(theta[i-1],x)
    
    if(log(runif(1))<logratio)
      theta[i] <- prop
    else
      theta[i] <- theta[i-1]
  }
  
  burn <- 300
  chain <- theta[(burn+1):iter]
  
  c(mean(chain),sd(chain))
}

# 1000 theta values

theta_values <- runif(1000,2,4)

n <- 50

theta_est <- numeric(1000)
se_est <- numeric(1000)

# Simulation

for(j in 1:1000){
  
  theta_true <- theta_values[j]
  
  x <- rztgid(n,theta_true)
  
  res <- mcmc_theta(x)
  
  theta_est[j] <- res[1]
  se_est[j] <- res[2]
  
  if(j%%50==0) cat("Replication:",j,"\n")
}

# Final result

mean_theta <- mean(theta_est)
mean_se <- mean(se_est)

cat("Mean of 1000 Estimated Theta =",mean_theta,"\n")
cat("Mean of 1000 Standard Error =",mean_se,"\n")

hist(theta_est,
     main="Histogram of 1000 Estimated Values of Theta",
     xlab=expression(hat(theta)),
     col="skyblue",
     border="black",
     breaks=30)
abline(v = mean(theta_est), col = "red", lwd = 2)

hist(se_est,
     main="Histogram of 1000 Estimated Standard Errors",
     xlab="Standard Error",
     col="lightgreen",
     border="black",
     breaks=30)
abline(v = mean(se_est), col = "red", lwd = 2)


#n=30
rm(list=ls())

# Precompute combination
max_x <- 20
comb_list <- lapply(1:max_x,function(x) choose(x,0:x))

# ZTGID pmf

pmf_ztgid <- function(x,theta){
  
  const <- theta^3/(theta^3+2)
  
  k <- 0:x
  
  num <- sum((-1)^k * comb_list[[x]] *
               (theta/(theta+k+1)+2/(theta+k+1)^3))
  
  den <- 1-const*(theta/(theta+1)+2/(theta+1)^3)
  
  const*num/den
}

# Random sample generator

rztgid <- function(n,theta){
  
  xvals <- 1:max_x
  
  probs <- sapply(xvals,pmf_ztgid,theta=theta)
  probs <- probs/sum(probs)
  
  sample(xvals,n,replace=TRUE,prob=probs)
}

#  log likelihood

loglik <- function(theta,x){
  
  if(theta<=0) return(-Inf)
  
  n <- length(x)
  
  term1 <- n*log(theta^3/(theta^3+2))
  
  inner <- numeric(n)
  
  for(i in 1:n){
    
    xi <- x[i]
    k <- 0:xi
    
    inner[i] <- sum((-1)^k * comb_list[[xi]] *
                      (theta/(theta+k+1)+2/(theta+k+1)^3))
  }
  
  term2 <- sum(log(inner))
  
  term3 <- -n*log(1-(theta^3/(theta^3+2))*
                    (theta/(theta+1)+2/(theta+1)^3))
  
  term1+term2+term3
}

#  MCMC

mcmc_theta <- function(x,iter=800){
  
  theta <- numeric(iter)
  theta[1] <- 2
  
  proposal_sd <- 0.2
  
  for(i in 2:iter){
    
    prop <- rnorm(1,theta[i-1],proposal_sd)
    
    logratio <- loglik(prop,x)-loglik(theta[i-1],x)
    
    if(log(runif(1))<logratio)
      theta[i] <- prop
    else
      theta[i] <- theta[i-1]
  }
  
  burn <- 300
  chain <- theta[(burn+1):iter]
  
  c(mean(chain),sd(chain))
}

# 1000 theta values

theta_values <- runif(1000,2,4)

n <- 30

theta_est <- numeric(1000)
se_est <- numeric(1000)

# Simulation

for(j in 1:1000){
  
  theta_true <- theta_values[j]
  
  x <- rztgid(n,theta_true)
  
  res <- mcmc_theta(x)
  
  theta_est[j] <- res[1]
  se_est[j] <- res[2]
  
  if(j%%50==0) cat("Replication:",j,"\n")
}

# Final result

mean_theta <- mean(theta_est)
mean_se <- mean(se_est)

cat("Mean of 1000 Estimated Theta =",mean_theta,"\n")
cat("Mean of 1000 Standard Error =",mean_se,"\n")

hist(theta_est,
     main="Histogram of 1000 Estimated Values of Theta",
     xlab=expression(hat(theta)),
     col="skyblue",
     border="black",
     breaks=30)
abline(v = mean(theta_est), col = "red", lwd = 2)

hist(se_est,
     main="Histogram of 1000 Estimated Standard Errors",
     xlab="Standard Error",
     col="lightgreen",
     border="black",
     breaks=30)
abline(v = mean(se_est), col = "red", lwd = 2)

