#rejection 20x*(1-x)^3

# Define the target density function
f <- function(x) {
  20 * x * (1 - x)^3
}
# Find the maximum of f(x) to determine our bounding constant
# We can find the maximum analytically by taking derivative:
# f'(x) = 20(1-x)^3 - 60x(1-x)^2
# Setting f'(x) = 0 => x = 1/4
x_max <- 1/4
c <- f(x_max)  # c ≈ 2.109375
# Define the proposal density (uniform on [0,1])
# g(x) = 1 for 0 < x < 1
# Rejection sampling algorithm
rejection_sampling <- function(n) {
  samples <- numeric(n)
  accepted <- 0
  
 while (accepted < n) {
    # Step 1: Generate Y from g (uniform)
    Y <- runif(1)
    
   # Step 2: Generate U from uniform(0,1)
    U <- runif(1)
    # Step 3: Accept or reject
    if (U <= f(Y)/(c * 1)) {  # Since g(Y) = 1
      accepted <- accepted + 1
      samples[accepted] <- Y
    }
  }
  
  return(samples)
}

# Generate 10,000 samples
set.seed(123)  # For reproducibility
n_samples <- 10000
samples <- rejection_sampling(n_samples)

# Plot the results to verify
hist(samples, breaks = 50, freq = FALSE, main = "Empirical vs Theoretical Density")
curve(f(x), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Empirical", "Theoretical"), 
       col = c("black", "red"), lty = c(1, 1))

# Print acceptance rate (should be about 1/c)
cat("Theoretical acceptance rate:", round(1/c, 4), "\n")
cat("Empirical acceptance rate:", round(n_samples/(n_samples * c), 4), "\n")
_______________________________________________-
#Rejection (kx^0.5e^-x)

# Define the target density function
target_density <- function(x) {
  if (x > 0) {
    return(x^(1/2) * exp(-x))
  } else {
    return(0)
  }
}

# Define the proposal density function (exponential with rate 1)
proposal_density <- function(x) {
  return(dexp(x, rate = 1))
}

# Define the constant M
M <- optimize(function(x) target_density(x) / proposal_density(x), interval = c(0, 10), maximum = TRUE)$objective

# Rejection sampling function
rejection_sampling <- function(n) {
  samples <- numeric(n)
  for (i in 1:n) {
    accept <- FALSE
    while (!accept) {
      # Sample from the proposal distribution
      x <- rexp(1, rate = 1)
      # Calculate the acceptance probability
      u <- runif(1)
      if (u < target_density(x) / (M * proposal_density(x))) {
        samples[i] <- x
        accept <- TRUE
      }
    }
  }
  return(samples)
}

# Generate samples
set.seed(123)
samples <- rejection_sampling(1000)

# Plot the histogram of the samples
hist(samples, breaks = 30, freq = FALSE, main = "Rejection Sampling for Gamma(3/2, 1)", xlab = "x")
curve(dgamma(x, shape = 3/2, rate = 1), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = "Gamma(3/2, 1)", col = "red", lwd = 2)
__________________________________________________________
#Rejection z(0,1)

# Define the target density function
f <- function(x) {
  (2 / sqrt(2 * pi)) * exp(-x^2 / 2)
}

# Define the proposal density function (e.g., exponential distribution with rate 1)
g <- function(x) {
  dexp(x, rate = 1)
}

# Find the maximum value of f(x) / g(x) to determine M
x_vals <- seq(0, 10, length.out = 1000)
ratio <- f(x_vals) / g(x_vals)
M <- max(ratio)

# Rejection sampling algorithm
n_samples <- 10000
samples <- numeric(n_samples)
accepted <- 0

while (accepted < n_samples) {
  # Sample from the proposal distribution
  x <- rexp(1, rate = 1)
  
  # Sample from a uniform distribution
  u <- runif(1)
  
  # Accept or reject the sample
  if (u < f(x) / (M * g(x))) {
    accepted <- accepted + 1
    samples[accepted] <- x
  }
}

# Since we generated the absolute value, we need to assign random signs
signs <- sample(c(-1, 1), n_samples, replace = TRUE)
normal_samples <- samples * signs

# Check the results
hist(normal_samples, breaks = 50, prob = TRUE, main = "Histogram of Generated Normal Samples")
curve(dnorm(x), add = TRUE, col = "red")
_____________________________________________________________________-
#  find value sqrt(2/pi)
# Set a random seed for reproducibility
set.seed(123)

# Generate 1 million standard normal random variables
n <- 10^6
z <- rnorm(n)  # No libraries needed - rnorm is in base R

# Calculate the mean of absolute values
mean_abs_z <- mean(abs(z))

# Theoretical value
theoretical_value <- sqrt(2/pi)

# Compare results
cat("Simulated E[|Z|]:", mean_abs_z, "\n")
cat("Theoretical value (√(2/π)):", theoretical_value, "\n")
cat("Difference:", abs(mean_abs_z - theoretical_value), "\n")
__________________________________________________________-

##target n(0,1)and proposed n(5,2)
# Set seed for reproducibility
set.seed(123)

# Define the target and proposal distributions
target_dist <- function(x) dnorm(x, mean = 0, sd = 1)
proposal_dist <- function(x) dnorm(x, mean = 5, sd = sqrt(2))

# Generate samples from the proposal distribution
n_samples <- 10000
samples <- rnorm(n_samples, mean = 5, sd = sqrt(2))

# Calculate the importance weights
weights <- target_dist(samples) / proposal_dist(samples)

# Estimate the expected value E(x^2)
expected_value <- mean(samples^2 * weights)

# Print the estimated expected value
cat("Estimated E(x^2):", expected_value, "\n")

# Create a density plot to compare the target and proposal distributions
x_values <- seq(-5, 10, length.out = 1000)
target_density <- target_dist(x_values)
proposal_density <- proposal_dist(x_values)

# Plot the densities using base R
plot(x_values, target_density, type = "l", col = "blue", lwd = 2,
     xlab = "x", ylab = "Density", main = "Comparison of Target and Proposal Distributions")
lines(x_values, proposal_density, col = "red", lwd = 2)
legend("topright", legend = c("Target N(0,1)", "Proposal N(5,2)"),
       col = c("blue", "red"), lwd = 2)
____________________________________________________________________-
##Big question
# Discrete hazard rate method for geometric random variable
simulate_geometric <- function(p) {
  X <- 1
  while (TRUE) {
    U <- runif(1) # Generate a random number U
    if (U < p) {
      return(X) # Stop if U < p
    }
    X <- X + 1 # Increment X
  }
}

# Example usage
set.seed(123) # For reproducibility
p <- 0.3 # Parameter for geometric distribution
n_simulations <- 1000 # Number of simulations
results <- replicate(n_simulations, simulate_geometric(p))

# Compare with theoretical mass function
theoretical_probs <- dgeom(0:max(results), prob = p)
empirical_probs <- table(results) / n_simulations

# Print results
cat("Theoretical probabilities:\n", theoretical_probs, "\n")
cat("Empirical probabilities:\n", empirical_probs, "\n")

# Plot results
plot(0:max(results), theoretical_probs, type = "h", col = "blue", lwd = 2,
     xlab = "X", ylab = "Probability", main = "Geometric Distribution")
lines(as.numeric(names(empirical_probs)), empirical_probs, type = "h", col = "red", lwd = 2)
legend("topright", legend = c("Theoretical", "Empirical"), col = c("blue", "red"), lwd = 2)
___________________________________________________________________-
#Antithetic & montecarlo x/2^x-1
h <- function(x) {
  return(x/2^(x-1))
}
monte_carlo <- function(k) {
  samples <- rexp(k, rate = 1)  #
  estimates <- h(samples)       #
  return(mean(estimates))       #
}
antithetic_variates <- function(k) {
  u <- runif(k / 2)             
  samples1 <- -log(u)           
  samples2 <- -log(1 - u)       #
  estimates1 <- h(samples1)     # 
  estimates2 <- h(samples2)     #
  return(mean(c(estimates1, estimates2)))  
}
# Number of samples
k <- c(1000,2000,60000)
mc_estimates <- replicate(100, monte_carlo(k))  
mc_mean <- mean(mc_estimates) 
mc_var <- var(mc_estimates)    
av_estimates <- replicate(100, antithetic_variates(k)) 
av_mean <- mean(av_estimates)  
av_var <- var(av_estimates)    
# Compare computational efficiency
time_mc <- system.time(replicate(100, monte_carlo(k)))  # Time for Monte Carlo
time_av <- system.time(replicate(100, antithetic_variates(k)))  # Time for Antithetic Variates

# Print results
cat("Monte Carlo Mean:", mc_mean, "\n")
cat("Monte Carlo Variance:", mc_var, "\n")
cat("Antithetic Variates Mean:", av_mean, "\n")
cat("Antithetic Variates Variance:", av_var, "\n")
cat("Monte Carlo Time:", time_mc[3], "seconds\n")
cat("Antithetic Variates Time:", time_av[3], "seconds\n")

# Visualization using base R
# Create a density plot for Monte Carlo and Antithetic Variates estimates
plot(density(mc_estimates), col = "blue", lwd = 2, main = "Comparison of Monte Carlo and Antithetic Variates Estimates",
     xlab = "Estimate", ylab = "Density", xlim = range(c(mc_estimates, av_estimates)))
lines(density(av_estimates), col = "red", lwd = 2)
legend("topright", legend = c("Monte Carlo", "Antithetic Variates"), col = c("blue", "red"), lwd = 2)
________________________________________________________________________

#Suppose Xis a random variable following the Exponential distribution with parameter λ.
Use Monte Carlo simulation to estimate the expectation of the function h(x) = e^-x

set.seed(123)  # For reproducibility

# Parameters
lambda <- 2    # Parameter for exponential distribution
k <- 10000     # Number of samples
true_value <- lambda / (1 + lambda)  # Theoretical value of E[e^-X]

# Basic Monte Carlo estimation
basic_mc <- function(k, lambda) {
  start_time <- Sys.time()
  X <- rexp(k, rate = lambda)
  hX <- exp(-X)
  estimate <- mean(hX)
  variance <- var(hX)/k
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  
  list(estimate = estimate, 
       variance = variance,
       time = as.numeric(time_taken),
       samples = hX)
}

# Control Variates estimation
control_variates <- function(k, lambda) {
  start_time <- Sys.time()
  X <- rexp(k, rate = lambda)
  hX <- exp(-X)
  
  # Use X itself as control variate (since we know E[X] = 1/lambda)
  c <- -cov(hX, X) / var(X)  # Optimal coefficient
  Y <- hX + c * (X - 1/lambda)
  
  estimate <- mean(Y)
  variance <- var(Y)/k
  end_time <- Sys.time()
  time_taken <- end_time - start_time
  
  list(estimate = estimate, 
       variance = variance,
       time = as.numeric(time_taken),
       samples = Y,
       c = c)
}

# Run both methods
basic_result <- basic_mc(k, lambda)
cv_result <- control_variates(k, lambda)

# Print results
cat("True value:", true_value, "\n\n")
cat("Basic Monte Carlo:\n")
cat("  Estimate:", basic_result$estimate, "\n")
cat("  Variance:", basic_result$variance, "\n")
cat("  Time (s):", basic_result$time, "\n\n")

cat("Control Variates:\n")
cat("  Estimate:", cv_result$estimate, "\n")
cat("  Variance:", cv_result$variance, "\n")
cat("  Time (s):", cv_result$time, "\n")
cat("  Optimal c:", cv_result$c, "\n\n")

cat("Variance reduction:", 100*(1 - cv_result$variance/basic_result$variance), "%\n")

# Visualization of variance reduction
par(mfrow = c(1, 2))

# Plot basic MC estimates
plot(cumsum(basic_result$samples) / (1:k), type = "l", 
     ylim = range(c(basic_result$estimate, cv_result$estimate, true_value)),
     main = "Basic Monte Carlo Convergence",
     xlab = "Number of samples", ylab = "Estimate")
abline(h = true_value, col = "red", lty = 2)
legend("topright", legend = c("MC Estimate", "True Value"), 
       col = c("black", "red"), lty = c(1, 2))

# Plot control variates estimates
plot(cumsum(cv_result$samples) / (1:k), type = "l",
     ylim = range(c(basic_result$estimate, cv_result$estimate, true_value)),
     main = "Control Variates Convergence",
     xlab = "Number of samples", ylab = "Estimate")
abline(h = true_value, col = "red", lty = 2)
legend("topright", legend = c("CV Estimate", "True Value"), 
       col = c("black", "red"), lty = c(1, 2))
___________________________________________________________________________-
#AR(1) Also, visualize the simulated
data using a time series plot of one MCMC.
set.seed(123)
n <- 1000
phi <- 0.8
sigma <- 1
y <- numeric(n)
y[1] <- rnorm(1, 0, sigma / sqrt(1 - phi^2))
for (t in 2:n) {
  y[t] <- phi * y[t-1] + rnorm(1, 0, sigma)
}
plot(y,type="l",col="blue",main="simulated AR(1) data")
grid()
acf(y)
____________________________________________________________________
#Run the MCMC simulation using MH algorithm for 10,000 iterations to generate
samples of the parameters.
#Estimate the variance of the posterior mean using the nonoverlapping batch mean
(NBM ) method. Also, provide the necessary plots.

phi_true <- 0.7  # True value for phi
sigma2_true <- 1 # True value for variance of epsilon
n <- 1000        # Number of observations

# Generate synthetic AR(1) data
epsilon <- rnorm(n, mean = 0, sd = sqrt(sigma2_true))
y <- numeric(n)
y[1] <- epsilon[1]  # Set initial value for y_1
for (t in 2:n) {
  y[t] <- phi_true * y[t-1] + epsilon[t]
}
# Metropolis-Hastings MCMC for AR(1) model with estimation of sigma2
mcmc_iterations <- 10000 # Number of MCMC iterations

# Initialize parameters for MCMC
phi_start <- 0.5    # Initial guess for phi
sigma2_start <- 1   # Initial guess for sigma^2
phi_values <- numeric(mcmc_iterations)  # To store the sampled phi values
sigma2_values <- numeric(mcmc_iterations) # To store the sampled sigma^2 values
phi_values[1] <- phi_start
sigma2_values[1] <- sigma2_start

# Metropolis-Hastings Algorithm for parameter estimation
proposal_sd_phi <- 0.05  # Standard deviation for proposal distribution of phi
proposal_sd_sigma2 <- 0.1 # Standard deviation for proposal distribution of sigma^2

for (i in 2:mcmc_iterations) {
  
  # Propose a new value for phi from a normal distribution centered around the current phi
  phi_proposal <- rnorm(1, mean = phi_values[i-1], sd = proposal_sd_phi)
  
  # Propose a new value for sigma^2 (scale parameter) from a normal distribution
  sigma2_proposal <- abs(rnorm(1, mean = sigma2_values[i-1], sd = proposal_sd_sigma2)) # Ensures positivity
  
  # Calculate the likelihood for the current and proposed values of phi and sigma^2
  log_likelihood_current <- sum(dnorm(y[2:n], mean = phi_values[i-1] * y[1:(n-1)], sd = sqrt(sigma2_values[i-1]), log = TRUE))
  log_likelihood_proposal <- sum(dnorm(y[2:n], mean = phi_proposal * y[1:(n-1)], sd = sqrt(sigma2_proposal), log = TRUE))
  
  # Calculate the acceptance ratio
  acceptance_ratio_phi <- min(1, exp(log_likelihood_proposal - log_likelihood_current))
  
  # Accept or reject the proposal for phi
  if (runif(1) < acceptance_ratio_phi) {
    phi_values[i] <- phi_proposal
  } else {
    phi_values[i] <- phi_values[i-1]  # Retain the previous value if the proposal is rejected
  }
  
  # Update the sigma^2 using the residuals from the AR(1) model
  residuals <- y[2:n] - phi_values[i] * y[1:(n-1)]
  log_likelihood_current_sigma2 <- -n / 2 * log(sigma2_values[i-1]) - sum(residuals^2) / (2 * sigma2_values[i-1])
  log_likelihood_proposal_sigma2 <- -n / 2 * log(sigma2_proposal) - sum(residuals^2) / (2 * sigma2_proposal)
  
  # Calculate the acceptance ratio for sigma^2
  acceptance_ratio_sigma2 <- min(1, exp(log_likelihood_proposal_sigma2 - log_likelihood_current_sigma2))
  
  # Accept or reject the proposal for sigma^2
  if (runif(1) < acceptance_ratio_sigma2) {
    sigma2_values[i] <- sigma2_proposal
  } else {
    sigma2_values[i] <- sigma2_values[i-1]  # Retain the previous value if the proposal is rejected
  }
}
mean(phi_values)
mean(sigma2_values)
# Function to compute the NBM variance estimator
nbm_variance <- function(chain, batch_size) {
  n <- length(chain)
  B <- n %/% batch_size  # Number of batches
  
  # Ensure the number of samples is a multiple of batch_size
  chain <- chain[1:(B * batch_size)]
  
  # Compute batch means
  batch_means <- sapply(1:B, function(j) mean(chain[((j-1)*batch_size + 1):(j*batch_size)]))
  
  # Compute overall mean
  overall_mean <- mean(batch_means)
  
  # Compute variance of batch means
  nbm_var <- (batch_size / (B - 1)) * sum((batch_means - overall_mean)^2)
  
  return(nbm_var)
}

# Set batch size
batch_size <- 500

# Compute NBM variance for phi and sigma^2
phi_nbm_var <- nbm_variance(phi_values, batch_size)
sigma2_nbm_var <- nbm_variance(sigma2_values, batch_size)

# Print results
cat("NBM Variance Estimate for Phi:", phi_nbm_var, "\n")
cat("NBM Variance Estimate for Sigma^2:", sigma2_nbm_var, "\n")

# Visualization: Trace plots and histogram
par(mfrow = c(2, 2))

# Trace plot for phi
plot(phi_values, type = "l", col = "blue", main = "Trace Plot of Phi", ylab = "Phi", xlab = "Iteration")

# Histogram of phi
hist(phi_values, col = "lightblue", main = "Histogram of Phi", xlab = "Phi", breaks = 30)

# Trace plot for sigma^2
plot(sigma2_values, type = "l", col = "red", main = "Trace Plot of Sigma^2", ylab = "Sigma^2", xlab = "Iteration")

# Histogram of sigma^2
hist(sigma2_values, col = "pink", main = "Histogram of Sigma^2", xlab = "Sigma^2", breaks = 30)
____________________________________________________________-
# Generate random samples from a Normal distribution using the Metropolis-Hastings
(M-H) algorithm and estimate its parameters. Then, implement a burn-in period to
remove initial bias and compare results using the Mean Squared Error (MSE).
Do it for different sample sizes (n=50,100,250)

# Function for Metropolis-Hastings Algorithm to sample from N(mu, sigma^2)
metropolis_hastings <- function(n, mu, sigma, proposal_sd, burn_in = 0) {
  samples <- numeric(n + burn_in)  # Store samples
  samples[1] <- rnorm(1, mean = mu, sd = sigma)  # Initialize with a random value
  
  for (i in 2:(n + burn_in)) {
    # Propose a new sample from a normal distribution centered at the current value
    proposal <- rnorm(1, mean = samples[i-1], sd = proposal_sd)
    
    # Compute acceptance probability
    acceptance_ratio <- dnorm(proposal, mean = mu, sd = sigma) / 
                        dnorm(samples[i-1], mean = mu, sd = sigma)
    
    # Accept or reject the proposed value
    if (runif(1) < acceptance_ratio) {
      samples[i] <- proposal
    } else {
      samples[i] <- samples[i-1]
    }
  }
  
  # Apply burn-in period (removing first 'burn_in' samples)
  return(samples[(burn_in + 1):(n + burn_in)])
}

# Function to compute Mean Squared Error (MSE)
mse <- function(estimated, true_value) {
  return(mean((estimated - true_value)^2))
}

# True parameters of the normal distribution
mu_true <- 5
sigma_true <- 2

# Proposal standard deviation
proposal_sd <- 1

# Burn-in period
burn_in_period <- 500

# Different sample sizes
sample_sizes <- c(50, 100, 250)

# Store results
results <- data.frame(Sample_Size = integer(), MSE_Without_Burnin = numeric(), MSE_With_Burnin = numeric())

# Run M-H algorithm for different sample sizes and compare MSE
for (n in sample_sizes) {
  # Without burn-in
  samples_no_burn <- metropolis_hastings(n, mu_true, sigma_true, proposal_sd, burn_in = 0)
  mu_est_no_burn <- mean(samples_no_burn)
  mse_no_burn <- mse(mu_est_no_burn, mu_true)
  
  # With burn-in
  samples_with_burn <- metropolis_hastings(n, mu_true, sigma_true, proposal_sd, burn_in = burn_in_period)
  mu_est_with_burn <- mean(samples_with_burn)
  mse_with_burn <- mse(mu_est_with_burn, mu_true)
  
  # Store results
  results <- rbind(results, data.frame(Sample_Size = n, 
                                       MSE_Without_Burnin = mse_no_burn, 
                                       MSE_With_Burnin = mse_with_burn))
}

# Print results
print(results)

# Visualization
par(mfrow = c(1, 2))

# Plot of samples without burn-in
plot(samples_no_burn, type = "l", col = "blue", main = "Trace Plot (No Burn-in)", xlab = "Iteration", ylab = "Sampled Value")

# Plot of samples with burn-in
plot(samples_with_burn, type = "l", col = "red", main = "Trace Plot (With Burn-in)", xlab = "Iteration", ylab = "Sampled Value")





