n <- 1e3 # number of loci we would like to simulate
rho <- 7.4e-7 # recombination rate
k_true <- 2 # switch-rate that we "integrate over" as a nuisance parameter
f_true <- 0.4 # proportion of genetic relatedness we are trying to infer
m1 <- 1 # multiplicity of infection for sample-1 that we are trying to infer
m2 <- 1 # multiplicity of infection for sample-2 that we are trying to infer
m_true <- c(m1, m2)
pos <- sort(sample(1.4e6, n)) # simulate some positions (genomic coordinates)
hist(rbeta(length(pos), shape1 = 10, shape2 = 10))
PLAF <- rbeta(length(pos), shape1 = 10, shape2 = 10)
sim <- HMMERTIME::simData(pos = list(contig1 = pos),
                          rho = rho,
                          m1 = m_true[1], m2 = m_true[2],
                          k = k_true,
                          f = f_true,
                          p = PLAF,
                          propMissing = 0)

ret <- HMMERTIME::runMCMC(vcfRobj = sim$vcfRobj, # vcfR object we simulated
                          vcfploid = 2, # ploidy of VCF
                          PLAF = PLAF,
                          m_max = 5, # max COI to consider
                          rho = rho, # recombination rate
                          k_max = 10, # max switch rate to consider
                          e1 = 0.05, # error for going from homozygous to heterozygous
                          e2 = 0.05, # error for going from heterozygous to homozygous
                          burnin = 1e4,
                          samples = 1e4,
                          reportIteration = 1e3,
                          verbose = TRUE,
                          parallelize = TRUE)
# store trueIBD from simulation between samples
trueIBD <- data.frame(CHROM = vcfR::getCHROM(sim$vcfRobj),
                      POS =vcfR::getPOS(sim$vcfRobj),
                      z_true = rowSums(sim$IBD[,1:ncol(sim$IBD),drop=FALSE]))

# true Find
mean(trueIBD$z_true)

ret$mcmcout[[1]]$summary$quantiles
plot(ret$mcmcout[[1]]$posteriors$f)
plot(ret$mcmcout[[1]]$posteriors$f_ind)

plotdf <- data.frame(f = unlist(ret$mcmcout[[1]]$iterations$f),
                     k = unlist(ret$mcmcout[[1]]$iterations$k))
plot(plotdf$var1, plotdf$var1.1)


x <- ret$mcmcout[[1]]
# get IBD matrix
CHROM <- x$summary$IBD_marginal[,1]
POS <- x$summary$IBD_marginal[,2]
IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])

IBDdf <- cbind.data.frame(CHROM, POS, IBD)
IBDdflong <- tidyr::gather(data=IBDdf, key="Z", value="Prob", 3:ncol(IBDdf))
IBDdflong$Znum <- as.numeric(gsub("z", "", IBDdflong$Z))

# split by chrom to avoid lag pos from diff chrom
IBDdflonglist <- split(IBDdflong, f=IBDdflong$CHROM)
IBDdflonglist <- lapply(IBDdflonglist, function(df){
  df$start <- dplyr::lag(df$POS)
  df$end <- df$POS
  return(df)
})
IBDdflong <- do.call("rbind", IBDdflonglist)
# filter unneccessary Znumbers
filtdf <- aggregate(IBDdflong$Prob, list(factor(IBDdflong$Znum)), sum)
if(any(filtdf == 0)){
  filtdf <- filtdf[which(filtdf[,2] == 0), 1]
  IBDdflong <- IBDdflong %>% dplyr::filter(! Znum %in% filtdf )
}

library(ggplot2)
ggplot() +
  geom_rect(data = IBDdflong, mapping = aes(xmin = start, xmax = end, ymin = Znum - 0.49, ymax = Znum + 0.49, fill = Prob)) +
  geom_line(data = trueIBD, aes(x = POS, y = z_true),
            colour = "#38A31A", size = 0.75) +
  viridis::scale_fill_viridis("IBD Probability", option = "plasma", limits = c(0,1)) +
  scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
  xlab("POS") +
  facet_grid(~CHROM) +
  theme_bw()

