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
b1 <- runif(M, min=0.4, max = 0.5)
b2 <- runif(M, min=0.4, max = 0.5)

phi <- 0.7
sigma <- 0.1

log_lambda <- matrix(NA, N, M)
lambda <- y <- matrix(NA, N, M)

for(m in 1:M) {
  log_lambda[1,m] <- rnorm(1, 0, 0.1)

  for(t in 2:N) {
    log_lambda[t,m] <- b0[m] + phi * log_lambda[t-1,m] + rnorm(1, 0, sigma)
  }

  lambda[,m] <- exp(log_lambda[,m])
  y[,m] <- rpois(N, lambda[,m])
}

head(y)

hist(y[,4])

# Simulate epsilon noise
y_star <- eps <- matrix(NA, nrow = N, ncol = M)

eps[,1] <- sample(c(-1, 0, 1), size = N, replace = TRUE, prob = c(0.4, 0.4, 0.2))
eps[,2] <- sample(c(-1, 0, 1), size = N, replace = TRUE, prob = c(0.2, 0.7, 0.1))
eps[,3] <- sample(c(-1, 0, 1), size = N, replace = TRUE, prob = c(0.33, 0.34, 0.33))
eps[,4] <- sample(c(-1, 0, 1), size = N, replace = TRUE, prob = c(0.2, 0.5, 0.3))


# Disturb bottom series
for(m in 1:M) {
  y_star[,m] <- ifelse(y[,m] + eps[,m]<0, 0, y[,m] + eps[,m]) # integer values
}

colnames(y_star) <- c("AA", "AB", "BA", "BB")

Tot = apply(y, 1, sum)
Tot = ifelse(Tot < 0, 0, Tot)

# Put into a wide format (for export)

poisson_sim_data <-tibble(Time = 1:c(N-burnin), Tot = Tot[-c(1:burnin)], as.data.frame(y_star[-c(1:burnin),]))

usethis::use_data(poisson_sim_data, overwrite=TRUE)

# Plot the hierarchical time series
library(ggplot2)
hts_data_long <- tidyr::pivot_longer(wide, cols = -X, names_to = "Level", values_to = "Value")

ggplot(hts_data_long, aes(x = X, y = Value, color = Level)) +
  geom_line() +
  labs(title = "Simulated Hierarchical Integer Time Series",
       x = "Time", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ Level, scales='free')
ggsave(filename = paste0(vispath, "simulated_poisson_glm_data.pdf"))

