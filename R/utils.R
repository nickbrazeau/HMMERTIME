#' @title Convert a \code{vcfRobject} to ----TEMP -- polyIBDinput
#' @description This function wraps the \code{vcfR} package to convert a vcffile or \code{vcfR} object to a SNP matrix suitable for polyIBD.
#' @param vcffile is the file path that for the VCF. This VCF will be converted to an object of class \code{vcfR}.
#' @param vcfR is the file path that for the VCF. This VCF will be converted to an object of class \code{vcfR}.
#'
#' @export


# TODO:
#     Change these functions to one call as an S4 object with slots for gtmatrix, GT, PopAf, and d
#     better class definitions
#     snpmatrixlist needs to be in an apply loop not for loop
#     better management of memory with p

vcf2polyIBDinput <- function(vcffile=NULL, vcfR=NULL) {


  # check if input is vcfR
  if(is.null(vcffile)){
    if(! class(vcfR) %in% "vcfR" ){
      stop("vcfR object must be of class vcfR")
    }
    # consistent naming
    vcf <- vcfR

    # check if input is vcf file path
  } else if (!is.null(vcffile)){
    vcf <- vcfR::read.vcfR(file=vcffile, verbose=F) # read vcf
  } else {
    stop("Must specify an input")
  }

  if( ! identical( vcf, vcfR::extract.indels(vcf[vcfR::is.biallelic(vcf)], return.indels = F) ) ){
    stop("Your vcf input contains variants that are not biallelic SNPs. polyIBD
          is restricted to use cases with biallelic SNPs. To convert your vcfR
          object to biallelic SNPs, you can use the following command:
          vcfR::extract.indels(vcf[vcfR::is.biallelic(vcf)], return.indels = F)
         ")
  }


  # -----------------------------------------------------
  # determine ploidy to determine genotype numeric placeholder
  #------------------------------------------------------
  ploidy <- stringr::str_extract(vcf@meta[grepl("ploidy", vcf@meta)], "ploidy=[0-9]")
  ploidy <- stringr::str_split(t(ploidy), "=", simplify = T)
  ploidy <- as.numeric(ploidy[1,2])

  if(ploidy == 2){

    gtmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=F) # numeric as T doesn't parse 0/1 correctly
    gtmatrix[gtmatrix == "0/0"] <- 0
    gtmatrix[gtmatrix == "0/1"] <- 1
    gtmatrix[gtmatrix == "1/1"] <- 2
    gtmatrix[is.na(gtmatrix)] <- -1
    gtmatrix <- apply(gtmatrix, 2, function(x){as.numeric(x)}) # need to convert from char (--dependent on case of "/") to numeric
  } else if(ploidy == 1){
    gtmatrix <- vcfR::extract.gt(vcf, element='GT', as.numeric=T)
    gtmatrix[gtmatrix == 0] <- 0
    gtmatrix[gtmatrix == 1] <- 2
    gtmatrix[is.na(gtmatrix)] <- -1
  } else {
    stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by polyIBD")
  }
  # -----------------------------------------------------
  # Store CHROM and POS from the VCF
  #------------------------------------------------------
  CHROM <- vcfR::getCHROM(vcf)
  POS <- vcfR::getPOS(vcf)

  # -----------------------------------------------------
  # Calculate pop AF, p
  #------------------------------------------------------
  L <- length(POS)
  p <- rep(NA, length(POS))
  p <- apply(gtmatrix, 1,
             function(x){(2*(length(which(x==0))) + length(which(x==1)))/(2*length(x))}) # since we know homozygous ref is 0, so this counts as 2 As, and then we count hets and then divide by 2*possible alleles. Doing this roundabout way because we aren't in HWE


  # -----------------------------------------------------
  # return
  #-----------------------------------------------------
  retlist <- list(CHROMPOS = cbind.data.frame(CHROM, POS),
                  gtmatrix = gtmatrix,
                  p=p)
  class(retlist) <- "polyIBDinput"
  return(retlist)

}



# -----------------------------------
# polyIBDinput_to_Rcppcompat
# Takes the polyIBDinput and makes it easier to be parsed for the Rcpp args.
# (not exported)

polyIBDinput_to_Rcppcompat <- function(polyIBDinput){

  p <- polyIBDinput[["p"]]
  gtmatrix <- polyIBDinput[["gtmatrix"]]

  # extract basic parameters
  CHROMtab <- table(polyIBDinput$CHROMPOS[, 1])
  nc <- length(CHROMtab)
  cnames <- names(CHROMtab)
  n <- as.vector(CHROMtab)

  # get distances between SNPs. Distance=-1 between contigs, indicating infinite distance
  SNP_dist <- diff(polyIBDinput$CHROMPOS[, 2]) # second column in this class is POS
  SNP_dist[cumsum(n)[1:(nc-1)]] <- -1

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

  x <- 4*(gtmatrix[,1]+1) + (gtmatrix[,2]+1)

  ## return
  ret <- list(p=p,
              x=x,
              SNP_dist = SNP_dist)

  return(ret)

}




#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#------------------------------------------------
#' @description function for determining if object is of class polyIBD
#' @noRd
#' @export

is.polyIBD <- function(x) {
  inherits(x, "polyIBD")
}

#------------------------------------------------
#' @description overload print() function to print summary only
#' @noRd
#' @export
print.polyIBD <- function(x, ...) {

  # print summary only
  summary(x)

  # return invisibly
  invisible(x)
}

#------------------------------------------------
#' @description overload summary() function.
#' @noRD
#' @export
summary.polyIBD <- function(x, ...) {

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

#------------------------------------------------
#' @description function for determining if object is of class polyIBD
#' @noRD
#' @export

is.polyIBDinput <- function(x) {
  inherits(x, "polyIBDinput")
}


#------------------------------------------------
#' @description overload print() function to print summary only
#' @noRD
#' @export

print.polyIBDinput <- function(x, ...) {

  # print this output line
  cat("-------------------------------------- \n")
  cat(paste("There are", ncol(x$gtmatrix), "Samples"), "\n")
  cat(paste(nrow(x$gtmatrix), "Biallelic SNPs"), "\n")
  cat("-------------------------------------- \n")
  # return invisibly
  invisible(x)
}


#------------------------------------------------
#' @description overload summary() function to print summary only
#' @noRD
#' @export
summary.polyIBDinput <- function(x, ...) {

  # print this output line
  cat("-------------------------------------- \n")
  cat(paste("There are", ncol(x$gtmatrix), "Samples"), "\n")
  cat(paste(nrow(x$gtmatrix), "Biallelic SNPs"), "\n")
  cat("-------------------------------------- \n")
  # return invisibly
  invisible(x)
}




#------------------------------------------------
# R <> Cpp compatibility for trans probs
#------------------------------------------------

# -----------------------------------
#' @title mat_to_Rcpp
#' @description  takes matrix as input, converts to list format for use within Rcpp code. Needed for trans prob eigen method
#' @noRd

mat_to_Rcpp <- function(x) {
  return(split(x,f=1:nrow(x)))
}

# -----------------------------------
#' @title Rcpp_to_mat
#' @description Takes list format returned from Rcpp and converts to matrix. Needed for trans prob eigen method
#' @noRd

Rcpp_to_mat <- function(x) {
  ret <- matrix(unlist(x), nrow=length(x), byrow=TRUE)
  return(ret)
}
