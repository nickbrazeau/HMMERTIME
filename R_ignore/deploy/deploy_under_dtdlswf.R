get_truth_from_pairwise_arg <- function(arg, this_coi){
  # arg doesn't store host information so need to carry COI information
  if (!length(this_coi) == 2) {
    stop("must be a pairwise comparison")
  }
  # get connections
  conn <- purrr::map(arg, "c")
  # get timing of connections
  tm <- purrr::map(arg, "t")
  # find pairwise
  get_pairwise_ibd <- function(conni, tmi, this_coi) {
    smpl1con <- conni[1:this_coi[1]]
    smpl2con <- conni[(this_coi[1]+1):(cumsum(this_coi)[2])]

    #......................
    # get strictly between
    # i.e. ignoring internal IBD connections
    #......................
    eff_pairwiseIBD <- sum(smpl2con %in% 0:(this_coi[1]-1) )

    #......................
    # get true IBD
    #......................
    # connections between 1 and 2
    pwconn <- which(smpl2con %in% 0:(this_coi[1]-1) )
    locimatches <- rep(1, length(pwconn))
    # note we are 0 based in connections
    # note bvtrees always point left
    # catch if there are multiple matches within sample 2 to the pairwise
    # this is a coalescent true that looks like below if host COI is 2,2
    # c: -1 -1 1 2
    # t: -1 -1 5 1
    if (length(pwconn) != 0) {
      for (i in 1:length(pwconn)) {
        haplotypeindex <- this_coi[1] + pwconn[i] - 1 # -1 for 0-based
        internalconn <- which(smpl2con %in% haplotypeindex )
        if (length(internalconn) != 0) {
          for (i in 1:length(internalconn)) {
            internalhaplotypeplace <- this_coi[1] + internalconn[i] # here 1-based in R
            if (tmi[internalhaplotypeplace] < tmi[this_coi[1] + internalconn[i]]) { # here 1-based in R
              locimatches[i] <- locimatches[i] + 1
            }
          }
        }
      }
    }
    # within sample1 is easy since we start at 0
    withinIBD_smpl1 <<- sum(smpl1con %in% (1:this_coi[1]-1))
    # within sample2 adjust slightly for "offset"
    withinIBD_smpl2 <- sum(smpl2con %in%
                             this_coi[1]:(length(conni)-1))

    # return
    out <- list(pairwiseIBD = sum(locimatches),
                withinIBD_smpl1 = withinIBD_smpl1,
                withinIBD_smpl2 = withinIBD_smpl2,
                eff_pairwiseIBD = eff_pairwiseIBD)
    return(out)
  }
  # calculate
  numerator <- purrr::map2(.x = conn, .y = tm,
                           .f = get_pairwise_ibd, this_coi = this_coi)

  #......................
  # true between
  #......................
  pairwiseIBD <- sum(purrr::map_dbl(numerator, "pairwiseIBD"))/(min(this_coi) * length(conn)) # min combn * num Loci

  #......................
  # within
  #......................

  # -1 here for the SELF comparison
  withinIBD_smpl1 <- sum(purrr::map_dbl(numerator, "withinIBD_smpl1")) / ((this_coi[1]-1) * length(conn))
  withinIBD_smpl2 <- sum(purrr::map_dbl(numerator, "withinIBD_smpl2")) / ((this_coi[2]-1) * length(conn))

  # catch when MOI = 1 and no within possible, so set to 0
  withinIBD_smpl1 <- ifelse(is.nan(withinIBD_smpl1), 0, withinIBD_smpl1)
  withinIBD_smpl2 <- ifelse(is.nan(withinIBD_smpl2), 0, withinIBD_smpl2)

  #......................
  # effective between
  #......................
  get_effective_coi <- function(arg, hostcoi_index) {
    # get connections
    conn <- purrr::map(arg, "c")
    # get connections for this specific host
    conn <- lapply(conn, function(x)return(x[hostcoi_index]))
    conn <- unique(conn)
    connmat <- matrix(NA, ncol = length(hostcoi_index), nrow = length(conn))
    for (i in 1:nrow(connmat)) {
      connmat[i,] <- conn[[i]]
    }
    # look to see if all loci are coalesced w/in and only w/in for each strain
    clonecount <- apply(connmat, 2, function(x) {all(x %in% (hostcoi_index-1))})
    return(length(hostcoi_index) - sum(clonecount))
  }
  # run
  effcoi1 <- get_effective_coi(arg = arg, hostcoi_index = 1:this_coi[1])
  effcoi2 <- get_effective_coi(arg = arg, hostcoi_index = (this_coi[1]+1):sum(this_coi))
  eff_coi <- c(effcoi1, effcoi2)
  # now calculate effective IBD between
  effpairwiseIBD <- sum(purrr::map_dbl(numerator, "eff_pairwiseIBD"))/(min(eff_coi) * length(conn)) # min combn * num Loci



  # return
  ret <- list(pairwiseIBD = pairwiseIBD,
              effpairwiseIBD = effpairwiseIBD,
              eff_coi = eff_coi,
              withinIBD_smpl1 = withinIBD_smpl1,
              withinIBD_smpl2 = withinIBD_smpl2)
  return(ret)
}


#--------------------------------------------------------------------
# Purpose of this script is to RUN the simulations that will allow us to
#--------------------------------------------------------------------
devtools::load_all()
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

# # chrompos
# pos <- readRDS("R_ignore/deploy/simparams/sim_POS.rds")
# brkpts <- which(diff(pos) > 1e8)
# brkpts <- c(0, brkpts, length(pos))
# brkpts <- diff(brkpts)
# chrmnms <- rplasmodium::chromnames()[1:14]
# CHROMPOS <- tibble::tibble(CHROM = rep(chrmnms, brkpts),
#                            POS = pos)
# # liftover polysim trick
# CHROMPOS <- CHROMPOS %>%
#   dplyr::left_join(., tibble::tibble(chrmnm = 1:14,
#                                      CHROM = rplasmodium::chromnames()[1:14])) %>%
#   dplyr::mutate(fixlft = 1e9*chrmnm,
#                 POS = POS - fixlft) %>%
#   dplyr::select(c("CHROM", "POS"))

pos <- sort(sample(1.4e6, 1e3))
CHROMPOS <- tibble::tibble(CHROM = "CHROM1",
                           POS = pos)


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
                hosts = list(1:2))

#............................................................
# pick one to run and sim forward
#...........................................................
row <- 95 # 125

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
polySimIBD::plot_coalescence_trees(ARG, loci = 1)
unique(purrr::map(ARG, "c"))
this_coi <- swfsim$coi[hosts]

#......................
# intervals
#......................
intv <- purrr::map(ARG, "c")
intv <- purrr::map_chr(intv, function(x){paste(x, collapse = "")})
intv <- tapply(seq_along(intv), intv, identity)[unique(intv)]
intvdf <- tibble::tibble(nm = rep(names(intv), sapply(intv, length)),
                         loci = unlist(intv)) %>%
  dplyr::arrange(loci)

# extract Haplotype Matrix
hapmat <- polySimIBD::get_haplotype_matrix(ARG)

# simulate Reads
WSAF.list <- polySimIBD::sim_biallelic(COIs = this_coi,
                                       haplotypematrix = hapmat,
                                       shape1 = 1.544,
                                       shape2 = 0.620,
                                       coverage = 100,
                                       alpha = 1,
                                       overdispersion = 0, # it's overdispersion that is causing serious headache in wsaf -- prob fine for MIPs where we were making clonal calls, less here
                                       epsilon = 0.025)

# get True IBD
trueIBD <- get_truth_from_pairwise_arg(arg = ARG,
                                       this_coi = this_coi)

trueIBD$pairwiseIBD
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
gtmatsim <- matrix("0/1", nrow = nrow(WSAF.list$NRWSAcounts), ncol = ncol(WSAF.list$NRWSAcounts))
gtmatsim[WSAF.list$NRWSAF > 0.95] <- "1/1"
gtmatsim[WSAF.list$NRWSAF < 0.05] <- "0/0"
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
                          PLAF = 1-WSAF.list$rbetaPLAF,
                          m_max = 10, # max COI to consider
                          rho = rho, # recombination rate
                          k_max = 10, # max switch rate to consider
                          e1 = 0.05, # error for going from homozygous to heterozygous
                          e2 = 0.05, # error for going from heterozygous to homozygous
                          burnin = 1e4,
                          samples = 1e4,
                          reportIteration = 1e3,
                          verbose = TRUE,
                          parallelize = TRUE)

ret$mcmcout[[1]]$summary$quantiles
trueIBD$pairwiseIBD
this_coi


plot(ret$mcmcout[[1]]$posteriors$f_ind)
plot(ret$mcmcout[[1]]$posteriors$k)



#............................................................
# arg
#...........................................................
plot_coalescence_trees(ARG, loci = 648)



#............................................................
# simple plots
#................................... ........................
get_IBD_climber <- function(ARG, this_coi, vcfRobj) {

  #......................
  # find intervals
  #......................
  intv <- purrr::map(ARG, "c")
  # coerce to char for factor
  intv <- purrr::map_chr(intv, function(x){paste(x, collapse = "")})
  # https://stackoverflow.com/questions/22993637/efficient-r-code-for-finding-indices-associated-with-unique-values-in-vector
  intv <- tapply(seq_along(intv), intv, identity)[unique(intv)]
  intvdf <- tibble::tibble(nm = rep(names(intv), sapply(intv, length)),
                           loci = unlist(intv)) %>%
    dplyr::arrange(loci)

  #......................
  # between relatedness
  #......................
  # get connections
  conn <- purrr::map(ARG, "c")
  # find pairwise
  get_pairwise_ibd <- function(conni, this_coi) {
    smpl1con <- conni[1:this_coi[1]]
    smpl2con <- conni[(this_coi[1]+1):(cumsum(this_coi)[2])]

    # get strictly between
    # i.e. ignoring internal IBD connections
    eff_pairwiseIBD <- sum(smpl2con %in% 0:(this_coi[1]-1) )
    out <- tibble::tibble(nm = paste(conni, collapse = ""),
                          eff_pairwiseIBD = eff_pairwiseIBD)
    return(out)
  }
  conn_lftover <- purrr::map(unique(conn), get_pairwise_ibd, this_coi = this_coi) %>%
    dplyr::bind_rows()
  #......................
  # bring together
  #......................
  ret <- dplyr::left_join(intvdf, conn_lftover, by = "nm")
  ret <- tibble::tibble(CHROM = vcfR::getCHROM(vcfRobj),
                        POS = vcfR::getPOS(vcfRobj),
                        eff_pairwiseIBD = ret$eff_pairwiseIBD)
  return(ret)
}


#......................
# bring together marginal and climber
#......................
x <- ret$mcmcout[[1]]
# get IBD matrix
CHROM <- x$summary$IBD_marginal[,1]
POS <- x$summary$IBD_marginal[,2]
IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])

IBDdf <- cbind.data.frame(CHROM, POS, IBD)
IBDdflong <- tidyr::pivot_longer(data=IBDdf, cols=-c("CHROM", "POS"), names_to="Z", values_to="Prob") %>%
  dplyr::mutate(Znum = as.numeric(gsub("z", "", Z))) %>%
  dplyr::group_by(CHROM, Znum) %>%
  dplyr::mutate(start = dplyr::lag(POS),
                start = ifelse(is.na(start), 0, start),
                end = POS) %>%
  dplyr::ungroup()



# filter unneccessary Znumbers
filtdf <- aggregate(IBDdflong$Prob, list(factor(IBDdflong$Znum)), sum)
if(any(filtdf == 0)){
  filtdf <- filtdf[which(filtdf[,2] == 0), 1]
  IBDdflong <- IBDdflong %>% dplyr::filter(! Znum %in% filtdf )
}

library(ggplot2)
ggplot() +
  geom_rect(data = IBDdflong, mapping = aes(xmin = start, xmax = end,
                                            ymin = Znum - 0.49, ymax = Znum + 0.49,
                                            fill = Prob)) +
  geom_line(data = get_IBD_climber(ARG, this_coi, vcfRobj),
            mapping = aes(x = POS, y = eff_pairwiseIBD),
            color = "#d9d9d9", size = 1.2, linetype = "dashed") +
  viridis::scale_fill_viridis("IBD Probability", option = "plasma", limits = c(0,1)) +
  scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
  xlab("POS") +
  facet_grid(~CHROM) +
  theme_bw()

#......................
# quick plots
#......................
table(extract.gt(vcfRobj))/sum(table(extract.gt(vcfRobj)))
vcfRmanip::plot_vcfRobj_GT(vcfRobj)
tibble::tibble(CHROM = vcfR::getCHROM(vcfRobj),
               POS = vcfR::getPOS(vcfRobj)) %>%
  dplyr::bind_cols(., as.data.frame(extract.gt(vcfRobj))) %>%
  magrittr::set_colnames(c("CHROM", "POS", "Sample1", "Sample2")) %>%
  dplyr::mutate(concord = ifelse(Sample1 == Sample2 & Sample1 != "0/1" & Sample2 != "0/1", 2,
                                 ifelse(Sample1 == "0/1" | Sample2 == "0/1", 1, 0)),
                concord = factor(concord, levels = c(2, 1, 0), labels = c("Con", "Het", "Dis")),
                POSfact = factor(POS, levels = POS)) %>%
  ggplot() +
  geom_tile(aes(x = POSfact, y = concord, fill = concord)) +
  scale_fill_manual("Concord Call", values=c("#AA0A3C", "#313695", "#fee090", "#cccccc"),
                    drop=FALSE) +
  scale_y_discrete(drop = FALSE) +
  ylab("Concordant Status") +
  xlab("Pos Factored") +
  facet_grid(CHROM~.) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold", family = "Arial"),
        axis.text.y = element_text(size=12, face="bold", family = "Arial"))





