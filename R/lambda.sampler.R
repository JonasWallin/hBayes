#' Initialize Parameters for Lambda Sampling
#'
#' This function prepares and returns all necessary parameters for sampling lambda in a Gibbs sampling process.
#' It computes sequences based on eigenvalues and specified quantiles for scaling with a chi-squared distribution.
#' Additionally, it precomputes the divisor for the chi-squared density used in the sampling process.
#'
#' @param E A list containing eigenvalues `E$values`.
#' @param n The degrees of freedom for the chi-squared distribution, usually equal to the sample size or number of observations.
#' @param I The number of intervals or segments for the piecewise constant prior model.
#' @param q_low The lower quantile (default is 0.01) used to compute the minimum lambda scaling.
#' @param q_high The upper quantile (default is 0.99) used to compute the maximum lambda scaling.
#'
#' @return A list containing parameters:
#'         - `a_seq` : A sequence of logarithmically scaled lambda values.
#'         - `a_vec` : A sequence of logarithmically scaled lambda values.
#'         - `a_ind`: Indices mapping `a_seq` values to intervals defined by the vector `a_vec`.
#'         - `mcmc_pi`: Initial probabilities for each interval, evenly distributed.
#'         - `exp_a_seq_div_n`: Precomputed divisor for the chi-squared density calculation.
#'
#' @examples
#' E <- list(values = runif(10, 1, 10))  # Simulated eigenvalues
#' init_params <- initLambdaSampling(E, n = 100, I = 64)
#' print(init_params)
initLambdaSampling <- function(E, n, I, q_low = 0.01, q_high = 0.99) {
    a_min <- log(min(E$values) * qchisq(q_low, n) / n)
    a_max <- log(max(E$values) * qchisq(q_high, n) / n)
    a_vec <- seq(a_min, a_max, length = I + 1)
    a_seq <- seq(a_min, a_max, length = I * 50)
    a_ind <- findInterval(a_seq, a_vec)
    a_ind[I * 50] <- I  # Adjust last index
    mcmc_pi <- rep(1 / I, I)

    exp_a_seq <- exp(a_seq)

    return(list(
        a_seq = a_seq,
        a_ind = a_ind,
        a_vec = a_vec,
        mcmc_pi = mcmc_pi,
        exp_a_seq= exp_a_seq  # Add this to the list of outputs
    ))
}

#' Sample Lambda Values for Gibbs Sampling
#'
#' This function samples lambda values using the parameters initialized by `initLambdaSampling`.
#' It employs the chi-squared distribution to compute likelihoods and updates lambda values based on posterior probabilities.
#'
#' @param lambda_hat Numeric vector representing the means of squared transformed data,
#' which is used in the chi-squared likelihood calculations.
#' @param params A list of parameters prepared by `initLambdaSampling`, containing:
#'               - `a_seq`: Sequence for sampling lambda values.
#'               - `a_ind`: Indices for mapping to piecewise constant prior probabilities.
#'               - `mcmc_pi`: Probabilities associated with each interval.
#'               - `exp_a_seq_div_n`: Precomputed values for chi-squared density calculations.
#'
#' @return A numeric vector of sampled lambda values.
#'
#' @examples
#' lambda_hat <- c(1.5, 2.0, 2.5)  # Simulated lambda_hat values
#' sampled_lambdas <- sampleLambda(lambda_hat, init_params)
#' print(sampled_lambdas)
sampleLambda <- function(lambda_hat,n, params) {
    # Extract parameters from the list
    a_seq <- params$a_seq
    a_ind <- params$a_ind
    mcmc_pi <- params$mcmc_pi
    exp_a_seq <- params$exp_a_seq

    # Preallocate the lambda vector
    sampled_lambda <- numeric(length(lambda_hat))

    for (j in seq_along(lambda_hat)) {
        # Compute likelihood for current component
        a_lik_j <- -0.5 * ( lambda_hat[j] / exp_a_seq  + n * a_seq)

        # Compute posterior probabilities
        a_pst_j <- a_lik_j + log(mcmc_pi[a_ind])
        a_pst_j <- a_pst_j - max(a_pst_j)
        a_pst_j <- exp(a_pst_j)
        a_pst_j <- a_pst_j / sum(a_pst_j)

        # Sample lambda based on posterior probabilities
        sampled_lambda[j] <- exp(sample(a_seq, 1, prob = a_pst_j))
    }

    return(sampled_lambda)
}
