# Function that does MCMC
do_mcmc = function(particle, log_likelihood, threshold)
{
    accepted = 0

    for(j in 1:mcmc_steps)
    {
        # Generate proposal by perturbing one or several of the underlying
        # coordinates
        new = particle
        num = 1
        if(runif(1) < 0.5)
            num = floor(num_params**runif(1))
        for(k in 1:num)
        {
            which = sample(1:num_params, 1)
            scale = 10.0**(1.0 - 3*abs(rt(1, df=2)))
            new[which] = new[which] + rnorm(1, sd=scale)
            new[which] = new[which] %% 1.0
        }

        # Accept the proposal if the new particle is above the threshold
        logL_new = log_likelihood(us_to_params(new))

        if(logL_new >= threshold)
        {
            particle = new
            log_likelihood = logL_new
            accepted = accepted + 1
        }
    }

    return(list(particle=particle,
                log_likelihood=log_likelihood, accepted=accepted))
}


# Useful function
logsumexp = function(xs)
{
    biggest = max(xs)
    xs = xs - biggest
    result = log(sum(exp(xs))) + biggest
    return(result)
}


