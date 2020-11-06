#' @description function for determining if object is of class HMMERTIME
#' @noRd
#' @export

is.HMMERTIME <- function(x) {
  inherits(x, "is.HMMERTIME")
}

#' @description overload print() function to print summary only
#' @noRd
#' @export
print.HMMERTIME <- function(x, ...) {

  # print summary only
  summary(x)

  # return invisibly
  invisible(x)
}

#' @description overload summary() function.
#' @noRd
#' @export
summary.HMMERTIME <- function(x, ...) {

  # print MCMC summary
  cat("# MCMC summary\n")
  cat(paste("burn-in iterations:\t", length(x$iterations$logLike_burnin)) ,"\n")
  cat(paste("sampling iterations:\t", length(x$iterations$logLike)) ,"\n")
  cat(paste("acceptance rate:\t", x$summary$accept_rate) ,"\n")
  cat(paste("run-time (seconds):\t", round(x$summary$runTime, 3)) ,"\n")
  cat("\n")

  # print posterior parameter summary
  cat("# Posterior estimates\n")
  quants <- x$summary$quantiles
  print(quants)
}

#' @description function for determining if object is of class HMMERTIME
#' @noRd
#' @export

is.HMMERTIMEinput <- function(x) {
  inherits(x, "HMMERTIMEinput")
}


#' @description overload print() function to print summary only
#' @noRd
#' @export

print.HMMERTIMEinput <- function(x, ...) {

  # print this output line
  cat("-------------------------------------- \n")
  cat(paste("There are", ncol(x$gtmatrix), "Samples"), "\n")
  cat(paste(nrow(x$gtmatrix), "Biallelic SNPs"), "\n")
  cat("-------------------------------------- \n")
  # return invisibly
  invisible(x)
}

#' @description overload summary() function to print summary only
#' @noRd
#' @export
summary.HMMERTIMEinput <- function(x, ...) {

  # print this output line
  cat("-------------------------------------- \n")
  cat(paste("There are", ncol(x$gtmatrix), "Samples"), "\n")
  cat(paste(nrow(x$gtmatrix), "Biallelic SNPs"), "\n")
  cat("-------------------------------------- \n")
  # return invisibly
  invisible(x)
}


