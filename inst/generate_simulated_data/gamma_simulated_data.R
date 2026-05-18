# This code generates bottom level time series for the simulation study from a Poisson stationary DGP
require(portes)
require(MASS)
require(tidyr)
require(tsibble)

# Set up variables
set.seed(1209)
N <- 2500
burnin <- 500
M <- 4 # Number of bottom level

# Set up coefficients of GLM ---------------------------------------------------
b0 <- runif(M, min=0.4, max = 0.5)

beta <- 5 # rate parameter of gamma

# In gamma: E(X) = alpha/beta

phi <- 0.7 # ar coefficient
sigma <- 0.1 # noise

log_exp <- matrix(NA, N, M)
log_alpha <- y <- y_star <- matrix(NA, N, M)

for(m in 1:M) {
  log_exp[1,m] <- rnorm(1, 0, 0.1)

  for(t in 2:N) {
    log_exp[t,m] <- b0[m] + phi * log_exp[t-1,m] + rnorm(1, 0, sigma) # log(beta * E(X))
  }

  log_alpha[,m] <- log(beta) + log_exp[,m]

  y[,m] <- stats::rgamma(N, shape = exp(log_alpha[,m]), rate = beta)
}

head(y)
colnames(y) <- c("AA", "AB", "BA", "BB")

# phi_bottom <- matrix(runif(nrow(log_alpha)*M, 0.5, 1), nrow(log_alpha), M)
#
# distort_alpha <- phi_bottom * exp(log_alpha)
#
# for(m in 1:M){
#   for(t in 1:N){
#     y_star[t,m] <- rgamma(1, shape = distort_alpha[t,m], rate = beta)
#   }
# }
# colnames(y_star) <- c("AA", "AB", "BA", "BB")
# plot(density(y[,4]))
# lines(density(y_star[,4]), col='red')
#
# plot(1:N, y[,4] - y_star[,4])

# ## --------------------------------------------------
# ## 2. Mid ~ Gamma(sum(alpha_m), beta)
# ## --------------------------------------------------
A_mid <- B_mid <- vector()
alpha_A <- rowSums(exp(log_alpha[,1:2]))
alpha_B <- rowSums(exp(log_alpha[,3:4]))

for(t in 1:N){
  A_mid[t] <- rgamma(1, shape = alpha_A[t], rate = beta)
  B_mid[t] <- rgamma(1, shape = alpha_B[t], rate = beta)
}


# ## --------------------------------------------------
# ## 3. Total ~ Gamma(sum(alpha_i), beta)
# ## --------------------------------------------------

sum_alpha <- apply(exp(log_alpha), 1, sum)

Tot <- stats::rgamma(n = length(sum_alpha),
                         shape = sum_alpha,
                         rate = beta)

# Put into a wide format (for export)

gamma_sim_data <-tibble(Time = 1:c(N-burnin),
                        Tot = Tot[-c(1:burnin)],
                        A = A_mid[-c(1:burnin)],
                        B = B_mid[-c(1:burnin)],
                        as.data.frame(y[-c(1:burnin),]))

usethis::use_data(gamma_sim_data, overwrite=TRUE)

gamma_sim_data

# Plot the hierarchical time series
library(ggplot2)
hts_data_long <- tidyr::pivot_longer(gamma_sim_data, cols = -Time, names_to = "Level", values_to = "Value")

ggplot(hts_data_long, aes(x = Time, y = Value, color = Level)) +
  geom_line() +
  labs(title = "Simulated Hierarchical Gamma Time Series",
       x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Level)

