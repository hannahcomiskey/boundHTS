# This code generates bottom level time series for the
# simulation study from a Binomial stationary DGP
require(portes)
require(MASS)
require(tidyr)
require(tsibble)

# Set up variables
set.seed(1209)
N <- 2000
M <- 4 # Number of bottom level

# Set up coefficients of GLM ---------------------------------------------------
b0 <- runif(M, min=0.4, max = 0.5)
b1 <- runif(M, min=0.4, max = 0.5)
b2 <- runif(M, min=0.4, max = 0.5)

# Set up x values
x <- runif(N, min = 0, max = 1)

# Calculate estimated lambdas
y <- lambda <- matrix(NA, nrow = N, ncol = M)
for(m in 1:M) {
  lambda[,m] <- exp(b0[m] + b1[m]*x + b2[m]*x^2)
  y[,m] <- rpois(N, lambda = lambda[,m])
}

head(y)

hist(y[,4])

# Simulate epsilon noise
y_star <- eps <- matrix(NA, nrow = N, ncol = M)

eps[,1] <- sample(c(-4:2), size = N, replace = TRUE, prob = c(0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1))
eps[,2] <- sample(c(-2, -1, 0, 1), size = N, replace = TRUE, prob = c(0.2, 0.3, 0.4, 0.1))
eps[,3] <- sample(c(-1, 0, 1, 2, 3), size = N, replace = TRUE, prob = c(0.2, 0.2, 0.2, 0.2, 0.2))
eps[,4] <- sample(c(-3, -2, 0, 1, 2, 3, 4), size = N, replace = TRUE, prob = c(0.126, 0.235, 0.235, 0.031, 0.206, 0.031, 0.136))


# Disturb bottom series
for(m in 1:M) {
  y_star[,m] <- ifelse(y[,m] + eps[,m]<0, 0, y[,m] + eps[,m]) # integer values
}

colnames(y_star) <- c("AA", "AB", "BA", "BB")

Tot = apply(y, 1, sum)
Tot = ifelse(Tot < 0, 0, Tot)

# Put into a wide format (for export to csv)

poisson_sim_data <-tibble(X=x,Tot, as.data.frame(y_star))

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

