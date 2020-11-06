#' @title Get transition probabilities
#' @description Takes values of f, rho, and zmax. Produces rate matrix and calculates eigen values and vectors.
#' @param f TODO
#' @param rho TODO
#' @param k TODO
#' @param z_max TODO
#' @noRd

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



#' @title mat_to_Rcpp
#' @description  takes matrix as input, converts to list format for use within Rcpp code. Needed for trans prob eigen method
#' @noRd

mat_to_Rcpp <- function(x) {
  return(split(x,f=1:nrow(x)))
}

#' @title Rcpp_to_mat
#' @description Takes list format returned from Rcpp and converts to matrix. Needed for trans prob eigen method
#' @noRd

Rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}

#' @title Pairwise Genetic Combinations
#' @description Get genotype pairwise combinations for Cpp
#' @noRd

pair_gen_combns <- function(pairmatrix) {
  # compare two samples and save comparison type in vector x
  # x is an integer vector with values in 0:15. These values indicate genotype combinations that cycle through the four options: {missing, homo REF, het, homo ALT} in the first sample, then the same four options in the second sample, leading to 16 options in total
  # 0 = {NA, NA}
  # 1 = {NA, A}
  # 2 = {NA, Aa}
  # 3 = {NA, a}
  # 4 = {A, NA}
  # 5 = {A,A}
  # 6 = {A,Aa}
  # 7 = {A,a}
  # 8 = {Aa, NA}
  # 9 = {Aa,A}
  # 10 = {Aa,Aa}
  # 11 = {Aa,a}
  # 12 = {a, NA}
  # 13 = {a,A}
  # 14 = {a,Aa}
  # 15 = {a,a}

  ret <- 4*(pairmatrix[,1]+1) + (pairmatrix[,2]+1)
  return(ret)
}


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
#' @description words todo
#' @return words todo
#' @importFrom magrittr %>%
#' @export


runMCMC <- function(vcfRobj = NULL,
                    vcfploid = 2,
                    m_max = 5,
                    k_max = 10,
                    rho = 1e-5,
                    e1 = 0.05,
                    e2 = 0.05,
                    burnin = 1e2,
                    samples = 1e3,
                    verbose = TRUE,
                    reportIteration = 1e3,
                    parallelize = TRUE) {

  #............................................................
  # Assertions
  #...........................................................
  #......................
  # assertions and checks on vcf
  #......................
  # check if input is vcfR
  if(! class(vcfRobj) %in% "vcfR" ){
    stop("vcfR object must be of class vcfR")
  }
  # check if vcf is biallelic snps
  if( ! identical(vcfRobj, vcfR::extract.indels(vcfRobj[vcfR::is.biallelic(vcfRobj)], return.indels = F) ) ){
    stop("Your vcf input contains variants that are not biallelic SNPs. HMMERTIME
          is restricted to use cases with biallelic SNPs. To convert your vcfR
          object to biallelic SNPs, you can use the following command:
          vcfR::extract.indels(vcfRobj[vcfR::is.biallelic(vcfRobj)], return.indels = F)
         ")
  }
  # check ploidy
  assert_in(vcfploid, c(1, 2),
            message = "VCF ploidy %s must be in either 1 or 2")

  ploidy <- as.numeric(sub("ploidy=", "", stringr::str_extract(vcfRobj@meta[grepl("ploidy", vcfRobj@meta)], "ploidy=[0-9]")))
  assert_eq(vcfploid, ploidy,
            message = "You have indcated that your VCF ploidy is %s but the
                        ploidy in the VCF metadata is %s")
  #......................
  # assertions on input parameters
  #......................
  assert_single_pos_int(m_max)
  assert_single_pos_int(k_max)
  assert_single_bounded(rho) # recombination rate bounded by 0 and 1
  assert_single_bounded(e1)
  assert_single_bounded(e2)
  assert_single_pos_int(burnin)
  assert_single_pos_int(samples)
  assert_single_pos_int(reportIteration)



  #............................................................
  # Wrangle data for Cpp input
  #...........................................................
  #......................
  # pieces for input
  #......................
  gtmatrix <- vcfRmanip::gtmat012(vcfRobj)
  # drop any rownames
  rownames(gtmatrix) <- NULL
  CHROM <- vcfR::getCHROM(vcfRobj)
  POS <- vcfR::getPOS(vcfRobj)
  # population loci-specific allele freq
  L <- length(POS)
  p <- rep(NA, length(POS))
  p <- apply(gtmatrix, 1,
             function(x){(2*(length(which(x==0))) + length(which(x==1)))/(2*length(x))}) # since we know homozygous ref is 0, so this counts as 2 As, and then we count hets and then divide by 2*possible alleles. Doing this roundabout way because we aren't in HWE
  # get distances between snps
  SNP_dist <- diff(POS)
  # find "jumps" between chroms
  chromtab <- table(CHROM)
  nc <- length(chromtab)
  n <- as.vector(chromtab)
  SNP_dist[cumsum(n)[1:(nc-1)]] <- -1

  #......................
  # check if pairwise extensions need to be made
  #......................
  # does vcf contain more than two samples
  if (ncol(gtmatrix) > 2) {

    # get combinations
    gtcombs <- combn(colnames(gtmatrix), 2)
    gtcombs_list <- split(gtcombs, 1:ncol(gtcombs))
    # split up pairwise
    pairmatrix_list <- lapply(gtcombs_list, function(x){gtmatrix[, x]})
    # liftover for genotypes
    pairmatrix_list <- lapply(pairmatrix_list, HMMERTIME:::pair_gen_combns)

  } else {

    pairmatrix_list <- HMMERTIME:::pair_gen_combns(gtmatrix)

  }

  # tidy this up for output
  tidyout <- tibble::tibble(smpl1 = gtcombs[1,],
                            smpl2 = gtcombs[2,]) %>%
    dplyr::mutate(pairmat = pairmatrix_list)

  # check
  assert_eq(length(pairmatrix_list), choose(ncol(gtmatrix), 2),
            message = "Issue with splitting your VCF. Check if your VCF
                       has duplicate sample names")

  #............................................................
  # print out message about pairwise user is about to do
  #...........................................................
  if (verbose) {
    cat(crayon::cyan("*~*~*~*~ HMMERTIME Run Summary ~*~*~*~*\n"))
    cat(crayon::green("Samples in VCF: "), ncol(gtmatrix), "\n")
    cat(crayon::magenta("Pairwise Comparisons: "), length(pairmatrix_list), "\n")
    cat(crayon::blue("Burn-in for Each: "), burnin, "\n")
    cat(crayon::blue("Sampling for Each: "), samples, "\n")
  }

  #............................................................
  # run efficient Rcpp function
  #...........................................................
  # internal wrapper
  my_MCMC_wrapper <- function(pairmat) {

    # define list of "set" parameters to pass to Rcpp
    args <- list(x = pairmat,
                 p = p,
                 rho = rho,
                 SNP_dist = SNP_dist,
                 m_max = m_max,
                 k_max = k_max,
                 burnin = burnin,
                 samples = samples,
                 e1 = e1,
                 e2 = e2,
                 reportIteration = reportIteration)


    # R functions to pass to Rcpp
    args_functions <- list(getTransProbs = HMMERTIME:::getTransProbs)

    # MCMC
    output_raw <- runMCMC_cpp(args, args_functions)

    # check for convergence
    checkConvergence(output_raw$logLike_burnin, output_raw$logLike)

    # list of raw output
    raw_output <- list(logLike_burnin = coda::mcmc(output_raw$logLike_burnin),
                       logLike = coda::mcmc(output_raw$logLike),
                       m1 = coda::mcmc(output_raw$m1),
                       m2 = coda::mcmc(output_raw$m2),
                       f = coda::mcmc(output_raw$f),
                       k = coda::mcmc(output_raw$k))

    # get marginal IBD matrix
    IBD_marginal <- t(Rcpp_to_mat(output_raw$IBD_marginal))
    colnames(IBD_marginal) <- paste0("z", 0:(ncol(IBD_marginal)-1))
    # tidy up
    IBD_marginal <- IBD_marginal %>%
      dplyr::mutate(CHROM = CHROM,
                    POS = POS) %>%
      dplyr::select(c("CHROM", "POS", dplyr::everything()))

    # get final acceptance rate
    accept_rate <- output_raw$accept_rate/samples

    # get time it took MCMC to run
    runTime <- output_raw$runTime

    # calculate quantiles over parameters
    quants <- t(mapply(function(x){quantile(x, probs=c(0.025, 0.5, 0.975))}, raw_output))
    quants <- quants[rownames(quants) %in% c("m1", "m2", "f", "k"),]

    # list of summary output
    summary_output <- list(IBD_marginal = IBD_marginal,
                           quantiles = quants,
                           accept_rate = accept_rate,
                           runTime = runTime)

    # create output
    ret <- list(summary = summary_output,
                iterations = raw_output)
    # return
    return(ret)
  }

  # wrapper over pairwise
  if (!parallelize) {
    tidyout <- tidyout %>%
      dplyr::mutate(mcmcout = furrr::future_map(pairmat, my_MCMC_wrapper))
  }


}





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


