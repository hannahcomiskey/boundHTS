rm(list = ls())
library(tidyverse)

set.seed(1209)

n_obs <- 400
burn_in <- 150
Time <- 1:n_obs

## --------------------------------------------------
## 1. Bottom nodes AA, AB with correlations
## --------------------------------------------------

# AR parameters
beta0 <- c(0.2, 0.1)
beta1 <- c(0.65, 0.55)

# Stationary mean (for centered parametrization)
mu <- beta0 / (1 - beta1)

# Innovation noise covariance
sd_AA <- sd_AB <- 0.15   # smaller reduces overshadowing
cor_AA_AB <- -0.6
Sigma <- matrix(c(sd_AA^2, sd_AA*sd_AB*cor_AA_AB,
                  sd_AA*sd_AB*cor_AA_AB, sd_AB^2), 2, 2)

# Simulate correlated AR(1)
logit_mat <- matrix(NA, n_obs, 2)
logit_mat[1,] <- rnorm(2, 0, 1)   # starting points

for(i in 2:n_obs) {
  logit_mat[i,] <- mu + beta1 * (logit_mat[i-1,] - mu) + MASS::mvrnorm(1, rep(0,2), Sigma)
}

# Drop burn-in to ensure stationarity
logit_AA <- logit_mat[-(1:burn_in), 1]
logit_AB <- logit_mat[-(1:burn_in), 2]

# Transform to proportions
AA <- plogis(logit_AA)
AB <- plogis(logit_AB)

plot(1:length(AA), AA, type='l')
plot(1:length(AB), AB, type='l')


cor.test(logit_AA, logit_AB, method='kendall')

# ## --------------------------------------------------
# ## 2. A_t = w_A * AA_t + (1-w_A) * AB_t
# ## --------------------------------------------------

A_mean <- 0.75*AA + 0.25*AB

phi_A <- 500  # low precision to minimise nose and keep signal
A <- rbeta(
  n = length(A_mean),
  shape1 = phi_A * A_mean,
  shape2 = phi_A * (1 - A_mean)
)
B <- 1 - A

plot(1:length(A), A -  A_mean, type='l') # +/- 5%

plot(1:length(A), A_mean, type='l', ylim=c(0.5, 0.8))
lines(1:length(A), A , col='red')

## --------------------------------------------------
## 3. Top split B = 1 - A
## --------------------------------------------------

B <- 1 - A

## --------------------------------------------------
## 6. Check final dataset
## --------------------------------------------------
plot(1:length(AA), AA, type='l', ylim=c(0.4,0.8))
lines(1:length(AB), AB, type='l', col='blue')
lines(1:length(A), A, type='l', col='green')
lines(1:length(A), A_mean, type='l', col='red')

## --------------------------------------------------
## 5. Final dataset
## --------------------------------------------------

beta_sim_data <- tibble(
  Time = 1:c(n_obs-burn_in),
  A  = A,
  B  = B,
  AA = AA,
  AB = AB
)

plot(beta_sim_data$Time, beta_sim_data$A, type='l')
plot(beta_sim_data$Time, beta_sim_data$B, type='l')
plot(beta_sim_data$Time, beta_sim_data$AA, type='l')
plot(beta_sim_data$Time, beta_sim_data$AB, type='l')

usethis::use_data(beta_beta_sim_data, overwrite = TRUE)


