# Define the number of parameters
num_params = 2

# Define the names of the parameters.
# If you can't be bothered, use names=rep(NA, num_params).
names = c("mu", "sigma")

# Step sizes for proposals
step_sizes = c(1, 1)

# If these are nonzero, they set the scale of a lognormal
# distribution for the step sizes.
sig_step_sizes = c(3, 3)

# A dataset
data = list(x = c(0.044103595329532, -0.922186227194728, 0.0960010565013454, 
-1.12828522333698, 1.36291783960003, -2.23065969145022, -0.93027082294735, 
-2.07175801570109, 1.20069146950599, 3.3937627148158), N = 10)

# A sensible starting point
starting_point = c(0.0, 1.0)
names(starting_point) = names

# Function that takes a vector of parameters and returns the log posterior,
# up to a normalising constant.
log_prob = function(params)
{
    log_prior = 0.0
    
    # Improper flat prior for mu requires no action
    
    # Improper log-uniform for sigma requires this
    if(params["sigma"] <= 0.0)
        return(-Inf)
    log_prior = log_prior - log(params["sigma"])

    log_likelihood = sum(dnorm(data$x,
                         mean=params["mu"], sd=params["sigma"], log=TRUE))

    return(log_prior + log_likelihood)
}

# Function that takes the current parameter values as argument
# and returns a proposed new set of parameter values. Remember
# that in R the original params has been passed by value, i.e., copied,
# so I can modify "params" without destroying the original
generate_proposal = function(params)
{
    # Choose which parameter to change
    k = sample(1:length(params), 1)

    # Change it
    step_size = step_sizes[k]*exp(sig_step_sizes[k]*rnorm(1))
    params[k] = params[k] + step_size*rnorm(1)

    # Return modified vector
    return(params)
}
