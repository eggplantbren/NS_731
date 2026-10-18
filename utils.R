# Function that does MCMC
do_mcmc = function(particle, logl, threshold)
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
            which = sample.int(num_params, 1)
            scale = 10.0**(1.0 - 3*abs(rt(1, df=2)))
            new[which] = new[which] + rnorm(1, sd=scale)
            new[which] = new[which] %% 1.0
        }

        # Accept the proposal if the new particle is above the threshold
        logl_new = log_likelihood(prior_transform(new))

        if(logl_new >= threshold)
        {
            particle = new
            logl = logl_new
            accepted = accepted + 1
        }
    }

    return(list(particle=particle,
                logl=logl, accepted=accepted))
}


# Useful function
logsumexp = function(xs)
{
    biggest = max(xs)
    xs = xs - biggest
    result = log(sum(exp(xs))) + biggest
    return(result)
}


