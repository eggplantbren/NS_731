# Import the model and some functions
source("ns-model.R")
source("utils.R")

# Number of particles
num_particles = 100

# MCMC steps per NS iteration
mcmc_steps = 1000

# Depth to go to. Set to Inf if you want it to try being automatic.
depth = Inf

# Number of NS iterations to do
steps = floor(num_particles*depth)

# Generate particles from the prior
# and calculate their log likelihoods
particles = array(dim=c(num_particles, num_params))
log_likelihoods = rep(NA, num_particles)
for(i in 1:num_particles)
{
    particles[i, ] = runif(num_params)
    log_likelihoods[i] = log_likelihood(us_to_params(particles[i, ]))
}

# Keep track of best posterior weight seen so far (arbitrary normalisation)
best = -Inf

# Start the output file
header = ""
for(i in 1:num_params)
    header = paste(header, parameter_names[i], ",", sep="")
header = paste(header, "log_likelihood", sep="")
fileConn = file("ns-output.csv", open="w")
writeLines(header, fileConn)
close(fileConn)

# Main NS loop
iteration = 1
while(TRUE)
{
    # Print message
    cat("Iteration ", iteration, ". ", sep="")

    # Find worst particle
    worst = which.min(log_likelihoods)

    # Save its details
    fileConn = file("ns-output.csv", open="a")
    line = ""
    params = us_to_params(particles[worst, ])
    for(i in 1:num_params)
        line = paste(line, params[i], ",", sep="")
    line = paste(line, log_likelihoods[worst], sep="")
    writeLines(line, fileConn)
    close(fileConn)

    threshold = log_likelihoods[worst]
    ln_p = log_likelihoods[worst] - iteration/num_particles
    if(ln_p > best)
        best = ln_p

    # Check for termination
    done = FALSE
    if(iteration == steps)
    {
        done = TRUE
    } else if(is.infinite(steps))
    {
        if(ln_p < best - log(1E6))
            done = TRUE
    }
    if(done)
    {
        cat("done.\n")
        break
    }

    cat("Generating new particle...")

    # Copy a survivor
    if(num_particles > 1)
    {
        while(TRUE)
        {
            which = sample(1:num_particles, 1)
            if(which != worst)
                break
        }
        particles[worst, ] = particles[which, ]
        log_likelihoods[worst] = log_likelihoods[which]
    }

    # Evolve within likelihood constraint using Metropolis
    newpoint = do_mcmc(particles[worst, ], log_likelihoods[worst], threshold)
    particles[worst, ] = newpoint$particle
    log_likelihoods[worst] = newpoint$log_likelihood
    accepted = newpoint$accepted

    cat("done. Accepted ", accepted, "/", mcmc_steps, " steps.\n", sep="")

    # Make a plot
    if(iteration %% num_particles == 0)
    {
        png(filename="ns-progress-plot.png", width=1000, height=1000)

        keep = as.matrix(read.csv("ns-output.csv"))

        logxs = -(1:iteration)/num_particles
        logws = logxs + keep[1:iteration, dim(keep)[2]]
        logws = logws - logsumexp(logws)

        # Smart ylim
        ylim = c()
        temp = sort(keep[1:iteration, dim(keep)[2]])
        if(length(temp) >= 2)
        {
            ylim[1] = temp[floor(0.1*length(temp))]
            ylim[2] = temp[length(temp)]
        }
        if(any(is.infinite(ylim)))
            ylim=c(0, 1)

        # Get plot window ready
        par(mfrow=c(2,1))
        par(mar=c(6,6,4,4))

        plot(logxs, keep[1:iteration, dim(keep)[2]],
             type="b", xlab="ln(X)", ylab="ln(L)", ylim=ylim,
             cex.lab=2, cex.axis=2)
        plot(logxs, exp(logws),
             type="b", xlab="ln(X)", ylab="Posterior weight", cex.lab=2,
             cex.axis=2)

        dev.off()
    }

    iteration = iteration + 1
} # End main NS loop

# Load output
keep = as.matrix(read.csv("ns-output.csv"))

# Prior weights
logws = -(1:iteration)/num_particles
logws = logws - logsumexp(logws)

# Calculate marginal likelihood
logZ = logsumexp(logws + keep[, dim(keep)[2]])

# Normalised posterior weights
post_weights = exp(logws + keep[, dim(keep)[2]] - logZ)

# ESS
ent = -sum(post_weights*log(post_weights + 1E-300))
ess = floor(exp(ent))

# Information
H = sum(post_weights*(keep[, dim(keep)[2]] - logZ))
err = sqrt(H/num_particles)

# Print results
cat("\n")
cat("Marginal likelihood: ln(Z) = ", logZ, " +- ", err, ".", sep="")
cat("\n")
cat("Information: H = ", H, " nats.", sep="")
cat("\n")
cat("Effective posterior sample size = ", ess, ".", sep="")
cat("\n")

# Create posterior samples by resampling
posterior_samples = array(dim=c(ess, dim(keep)[2]))

# Counter
k = 1
top = max(post_weights)
while(TRUE)
{
    # Choose one of the samples
    which = sample(1:iteration, 1)

    # Acceptance probability
    prob = post_weights[which]/top
    if(runif(1) <= prob)
    {
        posterior_samples[k, ] = keep[which, ]
        k = k + 1
    }

    # Check to see if we're done
    if(k == ess + 1)
        break
}

# Name the columns
colnames(posterior_samples) = colnames(keep)

# Save the output. Just posterior
write.table(posterior_samples, file="ns-posterior-samples.csv",
            sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat("Posterior samples saved in ns-posterior-samples.csv.\n")



