library(compiler)

# Define the number of parameters
num_params = 2

# A dataset
# This could easily be loaded from an external file
data = list(t=c(1,   2,  3,  4,  5,  6,  7,  8,  9, 10)-1,
            y=c(21, 19, 23, 40, 27, 22, 31, 39, 28, 42),
            N=10)

# Function that takes a vector of Uniform(0, 1) variables
# and returns a vector of the actual parameters. This
# function implicitly defines the prior for the parameters.
us_to_params = cmpfun(function(us)
{
    #### You'll only need to edit this function below this line ####

    # A lognormal prior for the first parameter
    # and a normal prior for the second parameter
    params = c(lambda0=exp(qnorm(us[1], 0, 10)),
               slope=qnorm(us[2], 0, 1))

    return(params)
})

# Function that takes a vector of parameters and returns the
# log likelihood.
log_likelihood = cmpfun(function(params)
{
    line = params["lambda0"]*exp(params["slope"]*data$t)
    logL = sum(dpois(data$y, line, log=TRUE))
    return(logL)
})

