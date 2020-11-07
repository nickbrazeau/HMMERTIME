get_truth_from_arg <- function(swfsim, arg, hosts = NULL){

  # need these details in order to know
  # which hosts you subsetted to to prune the
  # swfsim since the ARG doesn't store this information
  # choose hosts to subset to
  if(is.null(hosts)){
    hosts <- 1:length(swfsim$coi)
  }
  # find which elements in sim2 bvtrees correspond to haplotypes from these hosts
  this_coi <- swfsim$coi[hosts]


  # convert trees into matrix of alleles
  allele_mat <- polySimIBD::get_haplotype_matrix(arg)

  # split the haplotype matrix into individual (host) matrices
  hosts.haplotypes <- NULL
  splitter <- rep(x = 1:length(hosts), times = this_coi)
  for (i in 1:length(unique(splitter))) {
    hosthap <- allele_mat[, c( splitter == i ), drop = F]
    hosts.haplotypes <- c(hosts.haplotypes, list(hosthap))
  }

  #..............................................................
  # Find true IBD between samples
  #..............................................................
  # expand grid for combinations
  paircompar.long <- expand.grid(list( 1:length(hosts), 1:length(hosts) ) ) %>%
    magrittr::set_colnames(c("smpl1", "smpl2")) %>%
    dplyr::filter(smpl1 != smpl2) %>%
    tibble::as_tibble()

  paircompar.long$btwnIBD <- purrr::pmap(paircompar.long[,c("smpl1", "smpl2")],
                                         function(smpl1, smpl2){
                                           # get mat for pairwise
                                           allele_mat_i <- hosts.haplotypes[[smpl1]]
                                           allele_mat_j <- hosts.haplotypes[[smpl2]]

                                           # find number of haplotypes that are IBD between hosts
                                           overlap <- mapply(function(x) {
                                             length(intersect(allele_mat_i[x,], allele_mat_j[x,]))
                                           }, 1:nrow(allele_mat_i))

                                           # just want if any overlap
                                           overlap[overlap >= 1] <- 1

                                           # return
                                           ret <- list(
                                             locioverlap = overlap,
                                             btwn_IBDprop = sum(overlap)/length(overlap)
                                           )
                                           return(ret)

                                         })

  return(paircompar.long)
}

#--------------------------------------------------------------------
# Purpose of this script is to RUN the simulations that will allow us to
#--------------------------------------------------------------------
library(tidyverse)
library(polySimIBD)
#set.seed(1)

#............................................................
# Magic Numbers for Consideration
#...........................................................
# Genomic positions based on setting between chromosomes to a very high number.
# Effectively this breaks apart chromosomes for our forward simulation recombination runs
#
# Miles et al. 2016 (PMC5052046) & Taylor et al. 2019 (PMC6707449) gives us a recombination rate by 7.4e-7 M/bp
# Aimee gets this number by taking the inverse of Mile's estiamte of the CO recombination rate of 13.5 kb/cM
rho <- 7.4e-7

# chrompos
pos <- readRDS("R_ignore/deploy/simparams/sim_POS.rds")
brkpts <- which(diff(pos) > 1e8)
brkpts <- c(0, brkpts, length(pos))
brkpts <- diff(brkpts)
chrmnms <- rplasmodium::chromnames()[1:14]
CHROMPOS <- tibble::tibble(CHROM = rep(chrmnms, brkpts),
                           POS = pos)
# liftover polysim trick
CHROMPOS <- CHROMPOS %>%
  dplyr::left_join(., tibble::tibble(chrmnm = 1:14,
                                     CHROM = rplasmodium::chromnames()[1:14])) %>%
  dplyr::mutate(fixlft = 1e9*chrmnm,
                POS = POS - fixlft) %>%
  dplyr::select(c("CHROM", "POS"))


# tlim for arg
tlim <- 10

# NB, we will consider hosts (N) on a logarithmic-base 10 scale
N <- round(10^seq(1, 3, l = 11))
# Mean COIs for lambda
coilamdas <- readRDS("R_ignore/deploy/simparams/optim_lambda.RDS")

# various levels of M
M <- c(0, 0.25, 0.5, 1)

# expand out combinations
paramsdf <- expand.grid(N, coilamdas, M) %>%
  tibble::as_tibble() %>%
  magrittr::set_colnames(c("N", "mean_coi", "m"))

paramsdf <- paramsdf %>%
  dplyr::mutate(pos = list(pos),
                rho = rho,
                tlim = tlim,
                hosts = list(1:2)
  )

#............................................................
# pick one to run and sim forward
#...........................................................
row <- 60
# run forward
swfsim <- polySimIBD::sim_swf(pos = paramsdf[row, ]$pos[[1]],
                              N =  paramsdf[row, ]$N[[1]],
                              m =  paramsdf[row, ]$m[[1]],
                              rho =  paramsdf[row, ]$rho[[1]],
                              mean_coi =  paramsdf[row, ]$mean_coi[[1]],
                              tlim = 10)

# extract ARG and down sample arg
hosts <- c(1,2)
ARG <- polySimIBD::get_arg(swfsim, host_index = hosts)
ARG

# extract Haplotype Matrix
hapmat <- polySimIBD::get_haplotype_matrix(ARG)

# simulate Reads
this_coi <- swfsim$coi[hosts]
WSAF.list <- polySimIBD::sim_biallelic(COIs = this_coi,
                                       haplotypematrix = hapmat,
                                       shape1 = 1.544,
                                       shape2 = 0.620,
                                       coverage = 100,
                                       alpha = 1,
                                       overdispersion = 0.1,
                                       epsilon = 0.05)

# get True IBD
trueIBD <- get_truth_from_arg(swfsim = swfsim,
                              arg = ARG,
                              hosts = hosts)

trueIBD.btwn <- trueIBD %>%
  dplyr::mutate(btwnIBD = purrr::map(btwnIBD, "btwn_IBDprop")) %>%
  tidyr::unnest(cols = btwnIBD)
trueIBD.btwn
this_coi
#............................................................
# convert to VCF
#...........................................................

meta <- paste("##fileformat=VCFv4.3", "##Simulated with polySimIBD", "##ploidy=2", collapse = "\n")

# get Fix:  must be in this order and only these
fix <- data.frame(CHROM = CHROMPOS$CHROM,
                  POS = CHROMPOS$POS,
                  ID = NA,
                  REF = "A", # arbitrary choice
                  ALT = "T", # arbitrary choice
                  QUAL = NA,
                  FILTER = "PASS",
                  INFO = NA)
# get GT
# here we will need to make a choice on the "threshold" for "calling" a homozyg-ref, heterozyg, or homozyg-alt allele
# selecting 10% as 2x the simulated error rate
gtmatsim <- matrix("0/1", nrow = nrow(WSAF.list$NRWSAcounts), ncol = ncol(WSAF.list$NRWSAcounts))
gtmatsim[WSAF.list$NRWSAF > 0.9] <- "1/1"
gtmatsim[WSAF.list$NRWSAF < 0.1] <- "0/0"
# now that we have genotype calls, we can combine these into depths
gt <- matrix(NA, nrow = nrow(gtmatsim), ncol = ncol(gtmatsim))
# quick for loop to protect against vector setting
for (i in 1:nrow(gtmatsim)) {
  for (j in 1:ncol(gtmatsim)) {
    AD <- paste(WSAF.list$WS.coverage[i,j] - WSAF.list$NRWSAcounts[i,j], WSAF.list$NRWSAcounts[i,j], sep = ",")
    gt[i, j] <- paste(gtmatsim[i,j], AD, WSAF.list$WS.coverage[i,j], sep = ":")
  }
}
gt <- cbind(FORMAT = "GT:AD:DP", gt)

# write out new vcfRobj
require(vcfR) # need this for class
vcfRobj <- new("vcfR", meta = meta, fix = as.matrix(fix), gt = gt)

#............................................................
# run hmmertime
#...........................................................

ret <- HMMERTIME::runMCMC(vcfRobj = vcfRobj, # vcfR object we simulated
                          vcfploid = 2, # ploidy of VCF
                          PLAF = 1 - WSAF.list$rbetaPLAF,
                          m_max = 10, # max COI to consider
                          rho = rho, # recombination rate
                          k_max = 20, # max switch rate to consider
                          e1 = 0.05, # error for going from homozygous to heterozygous
                          e2 = 0.05, # error for going from heterozygous to homozygous
                          burnin = 1e4,
                          samples = 1e4,
                          reportIteration = 1e3,
                          verbose = TRUE,
                          parallelize = TRUE)

ret$mcmcout[[1]]$summary
trueIBD.btwn
this_coi

plot(ret$mcmcout[[1]]$posteriors$f_ind)
plot(ret$mcmcout[[1]]$posteriors$k)


tibble::tibble(iteration = 1:length(ret$mcmcout[[1]]$posteriors$f_ind),
               F_posterior = ret$mcmcout[[1]]$posteriors$f_ind) %>%
  ggplot() +
  geom_line(aes(x = iteration, y = F_posterior)) +
  geom_hline(yintercept = unique(trueIBD.btwn$btwnIBD), color = "red")

