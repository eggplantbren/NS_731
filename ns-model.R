# Define the number of parameters
num_params = 2

# Define the names of the parameters.
# If you can't be bothered, use names=rep(NA, num_params).
names = c("mu", "sigma")

# A dataset
data = list(x = c(0.044103595329532, -0.922186227194728, 0.0960010565013454, 
-1.12828522333698, 1.36291783960003, -2.23065969145022, -0.93027082294735, 
-2.07175801570109, 1.20069146950599, 3.3937627148158), N = 10)

# Function that takes a vector of Uniform(0, 1) variables
# and returns a vector of the actual parameters. This
# function implicitly defines the prior for the parameters.
us_to_params = function(us)
{
    # Vector to be returned as the result of the function
    params = rep(NA, num_params)

    # Apply the names
    names(params) = names

    #### You'll only need to edit this function below this line ####

    # A broad normal prior for mu
    params["mu"] = qnorm(us[1], mean=0, sd=1000)

    # A broad lognormal prior for sigma
    params["sigma"] = exp(qnorm(us[2], mean=0, sd=5))

    return(params)
}

# Function that takes a vector of parameters and returns the
# log likelihood.
log_likelihood = function(params)
{
    logL = sum(dnorm(data$x, mean=params["mu"], sd=params["sigma"], log=TRUE))
    return(logL)
}

