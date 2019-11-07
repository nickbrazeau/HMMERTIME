#------------------------------------------------
#' Run polyIBD MCMC using Rcpp functions
#'
#' @param input TODO
#' @param m_max TODO
#' @param k_max TODO
#' @param rho TODO
#' @param burnin TODO
#' @param samples TODO
#' @param e1 TODO
#' @param e2 TODO
#' @param reportIteration TODO
#'
#' @description words todo
#'
#' @return words todo
#'
#' @export


runMCMC <- function(input = NULL,
                    m_max = 5,
                    k_max = 10,
                    rho = 1e-5,
                    burnin = 1e2,
                    samples = 1e3,
                    e1 = 0.05, e2 = 0.05,
                    reportIteration = 1e3) {


  # TODO - input parameter checks
  # note - vcf must have 4 columns, samples in final two columns
  # better management of memory for p and vcf -- if you have it as one object, p continually gets copied and takes up way too much memory
  # assertions



  # ------------------------------
  #     Setup Input for Rcpp
  # ------------------------------

  Rcppcompat <- polyIBDinput_to_Rcppcompat(input)

  # define list of arguments to pass to Rcpp
  args <- list(x = Rcppcompat[["x"]],
               p = unlist(Rcppcompat[["p"]]),
               rho = rho,
               SNP_dist = Rcppcompat[["SNP_dist"]],
               m_max = m_max,
               k_max = k_max,
               burnin = burnin,
               samples = samples,
               e1 = e1,
               e2 = e2,
               reportIteration = reportIteration)

  # R functions to pass to Rcpp
  args_functions <- list(getTransProbs = polyIBD::getTransProbs)


  # ------------------------------
  #            RUN MCMC
  # ------------------------------

  # run efficient Rcpp function
  output_raw <- runMCMC_cpp(args, args_functions)

  # check for convergence
  checkConvergence(output_raw$logLike_burnin, output_raw$logLike)


  # ------------------------------
  #      SAVE MCMC RAW OUTPUT
  # ------------------------------

  # list of raw output
  raw_output <- list(logLike_burnin = coda::mcmc(output_raw$logLike_burnin),
                     logLike = coda::mcmc(output_raw$logLike),
                     m1 = coda::mcmc(output_raw$m1),
                     m2 = coda::mcmc(output_raw$m2),
                     f = coda::mcmc(output_raw$f),
                     f_ind = coda::mcmc(output_raw$f_ind),
                     k = coda::mcmc(output_raw$k)
                     )


  # ------------------------------
  #   SAVE MCMC SUMMARY OUTPUT
  # ------------------------------

  # get marginal IBD matrix
  IBD_marginal <- t(Rcpp_to_mat(output_raw$IBD_marginal))
  colnames(IBD_marginal) <- paste0("z", 0:(ncol(IBD_marginal)-1))
  IBD_marginal <- cbind(input$CHROMPOS, IBD_marginal)

  # get final acceptance rate
  accept_rate <- output_raw$accept_rate/samples

  # get time it took MCMC to run
  runTime <- output_raw$runTime

  # calculate quantiles over parameters
  quants <- t(mapply(function(x){quantile(x, probs=c(0.025, 0.5, 0.975))}, raw_output))
  quants <- quants[rownames(quants) %in% c("m1", "m2", "f", "f_ind", "k"),]

  # list of summary output
  summary_output <- list(IBD_marginal = IBD_marginal,
                         quantiles = quants,
                         accept_rate = accept_rate,
                         runTime = runTime)


  # ------------------------------
  #            RETURN
  # ------------------------------
  # create output class polyIBD
  ret <- list(samples = colnames(input$gtmatrix),
              summary = summary_output,
              iterations = raw_output)
  class(ret) <- "polyIBD"

  # return
  return(ret)
}


#------------------------------------------------
# Trans Probs
#------------------------------------------------
# ------------------------------------------------------------------
#' @title Get transition probabilities
#'
#' @description Takes values of f, rho, and zmax. Produces rate matrix and calculates eigen values and vectors.
#' @param f TODO
#' @param rho TODO
#' @param k TODO
#' @param z_max TODO
#' @export

getTransProbs <- function(f, rho, k, z_max) {

  # generate rate matrix
  z0 <- z_max
  z1 <- z_max + 1
  rateMat <- matrix(0, z1, z1)
  rateMat[cbind(1:z0, 1:z0 + 1)] <- (z0:1)*rho*k*f # fill in matrix with flows up from current state to next state
  rateMat[cbind(1:z0 + 1, 1:z0)] <- (1:z0)*rho*k*(1 - f) # ...
  rateMat[cbind(1:z1, 1:z1)] <- -rowSums(rateMat)

  # obtain Eigen values and vectors
  E <- eigen(t(rateMat))
  Esolve <- solve(E$vectors)

  return(
    list(
      Evalues =  E$values,
      Evectors = mat_to_Rcpp(E$vectors),
      Esolve =   mat_to_Rcpp(Esolve)
    )
  )
}




# -----------------------------------
# checkConvergence
# calculates Geweke statistic from a series of burn-in and sampling draws. Report whether burn-in length was sufficient based on this statistic.
# (not exported)

checkConvergence <- function(burnin, samples) {

  # get number of burnin and sampling iterations
  nburnin <- length(burnin)
  nsamples <- length(samples)

  # calculate Geweke diagnostic on combined chain
  chain <- coda::mcmc(c(burnin, samples))
  geweke_z <- coda::geweke.diag(chain, frac1 = nburnin/(nburnin+nsamples), frac2 = nsamples/(nburnin+nsamples))$z

  if(is.na(geweke_z)){
    stop("NaN p-value was calculated from Geweke statistic")
  }

  # convert to p-value
  geweke_p <- 2*pnorm(abs(geweke_z), lower.tail=FALSE)

  # report convergence
  if (geweke_p > 0.05) {
    cat(paste0("convergence reached within defined burn-in period (Geweke p=", round(geweke_p, 3), ")"))
  } else {
    cat(paste0("WARNING: convergence not reached within defined burn-in period (Geweke p=", round(geweke_p,3), ")"))
  }

}


