# ------------------------------------------------------------------
# Contains two sections:
#                       1) Data for the package
#                       2) Data simulation functions
# ------------------------------------------------------------------

#' P. falciparum Cross Project Data
#'
#' A variant call file (vcf) that contains SNPs on chromosome 1 for
#' the 3D7 and HB3 parents and the F1 progeny, C01 and C14.
#'
#' @format A vcf containing 4 samples and 100 biallelic SNPs
#' @section Warning: If you copy and paste the code below, it will write to
#'                   your Desktop/temp directory (if it exists).
#'
#' @section Dependencies: To generate this subsetted VCF, we used the
#'                        \code{NFBtools} package which is freely available
#'                        from GitHub with  \code{devtools::install_github("nickbrazeau/NFBtools")}
#'
#' @section Citation: This VCF has generously been made publicly available by
#' the MalariaGEN Consortium and the P. falciparum Genetic Cross project led by
#' Alistair Miles. The manuscript is available through NCBI (PMID 27531718).
#'
#' @format A \code{vcfRobject} with 3 samples and 100 biallelic SNPs
#' \describe{
#'   This file was generated with the following code:
#'   \dontrun{
#'   url <- "ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/3d7_hb3.combined.final.vcf.gz"
#'   destfile <- "~/Desktop/temp/temp.vcf.gz"
#'   httr::GET(url=url, httr::write_disk(path=destfile, overwrite = F))
#'   vcfRobject <- vcfR::read.vcfR(file=destfile)
#'   vcfRobject <- vcfR::extract.indels(vcfRobject[vcfR::is.biallelic(vcfRobject)], return.indels = F) # subset to SNPs
#'   vcfRobject <- vcfRmanip::vcfR2SubsetChrom(vcfRobject = vcfRobject, chromvect = "Pf3D7_01_v3")
#'   pfcross_subset <- vcfRobject[, c("FORMAT", "3D7/PG0051-C/ERR019061", "HB3/PG0052-C/ERR019054", "C04/PG0061-C/ERR019059", "C02/PG0053-C/ERR019067")]
#'   }
#' }
#'
#' @source \url{ftp://ngs.sanger.ac.uk/production/malaria/pf-crosses/1.0/3d7_hb3.combined.final.vcf.gz}
#'
"pfcross_subset"

