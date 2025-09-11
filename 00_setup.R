# Purpose -----------------------------------------------------------------
# Run this script to create any custom functions or load data for the main app

# Truncated zero inflated gamma ----
generate_filtered_ZAGA <- function(n, mu, sigma, nu, max_log = 8) {
  samples <- rZAGA(n, mu = mu, sigma = sigma, nu = nu)
  
  while(sum(samples > max_log) > 0) {               # While any values > max_log
    n_replace <- sum(samples > max_log)             # Count how many to replace
    new_vals <- rZAGA(n_replace,                    # Generate replacements
                      mu = mu, 
                      sigma = sigma, 
                      nu = nu)
    samples[samples > max_log] <- new_vals             # Replace the high values
  }
  
  return(samples)
}

# define distribution parameters that aren't user inputted
norm.incom.contam.mu <- 1.631
norm.incom.contam.sigma <- 1.157293

