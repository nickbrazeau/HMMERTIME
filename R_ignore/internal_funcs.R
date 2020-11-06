# This file contains internal functions for the package
#  these are not exported and are for processing data
#  overall structure of data processing is vcfR --> HMMERTIME --> MCMClist

#' @title Convert a \code{vcfRobject} to ----TEMP -- HMMERTIMEinput
#' @description This function wraps the \code{vcfR} package to convert a vcffile or \code{vcfR} object to a SNP matrix suitable for HMMERTIME.
#' @param vcfR An object of class \code{vcfR} for manipulation
#'
#' @noRd


# TODO:
#     Change these functions to one call as an S4 object with slots for gtmatrix, GT, PopAf, and d
#     better class definitions
#     snpmatrixlist needs to be in an apply loop not for loop
#     better management of memory with p

vcf2HMMERTIMEinput <- function(vcfR=NULL) {


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
    stop("Your vcf input contains variants that are not biallelic SNPs. HMMERTIME
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
    stop("You have a ploidy that is less than 1 or greater than 3, which cannot be accomodated by HMMERTIME")
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
  class(retlist) <- "HMMERTIMEinput"
  return(retlist)

}



# -----------------------------------
# HMMERTIMEinput_to_Rcppcompat
# Takes the HMMERTIMEinput and makes it easier to be parsed for the Rcpp args.
# (not exported)

HMMERTIMEinput_to_Rcppcompat <- function(HMMERTIMEinput){

  p <- HMMERTIMEinput[["p"]]
  gtmatrix <- HMMERTIMEinput[["gtmatrix"]]

  # extract basic parameters
  CHROMtab <- table(HMMERTIMEinput$CHROMPOS[, 1])
  nc <- length(CHROMtab)
  cnames <- names(CHROMtab)
  n <- as.vector(CHROMtab)

  # get distances between SNPs. Distance=-1 between contigs, indicating infinite distance
  SNP_dist <- diff(HMMERTIMEinput$CHROMPOS[, 2]) # second column in this class is POS
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
