# Source the model file
source("metropolis-model.R")

# Number of iterations to do
steps = 100000

# Thinning
thin = 10

# Initialise the particle somewhere
particle = starting_point
logp = log_prob(particle)

# Start the output file
header = ""
for(i in 1:num_params)
    header = paste(header, parameter_names[i], ",", sep="")
header = paste(header, "log_prob", sep="")
fileConn = file("metropolis-output.csv", open="w")
writeLines(header, fileConn)
close(fileConn)

# Storage of results for trace plot
keep = rep(NA, floor(steps/thin))

# Count the number of accepted proposals
accepted = 0

# Main Metropolis loop
for(iteration in 1:steps)
{
    # Generate proposal and measure how good it is
    proposal = generate_proposal(particle)
    proposal_logp = log_prob(proposal)

    # Compute acceptance probability
    log_alpha = proposal_logp - logp
    if(is.nan(log_alpha))
        log_alpha = -Inf
    if(log_alpha >= 0.0)
        alpha = 1.0
    else
        alpha = exp(log_alpha)


    # Make accept/reject decision
    if(runif(1) <= alpha)
    {
        # Accept the proposal
        particle = proposal
        logp = proposal_logp
        accepted = accepted + 1
    }

    if(iteration %% thin == 0)
    {
        # Save particle to disk
        fileConn = file("metropolis-output.csv", open="a")
        line = ""
        for(i in 1:num_params)
            line = paste(line, particle[i], ",", sep="")
        line = paste(line, logp, sep="")
        writeLines(line, fileConn)
        close(fileConn)

        # Print a message to screen.
        accept_frac = accepted/iteration
        cat(paste("Iteration", iteration, ". Acceptance rate =", accept_frac))
        cat("\n")
    }
}

