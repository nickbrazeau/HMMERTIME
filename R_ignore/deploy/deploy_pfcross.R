load("data/pfcross_subset.rda")
set.seed(48)
PLAF <- rbeta(nrow(pfcross_subset@gt), 0.1, 0.1)
# run model
ret <- HMMERTIME::runMCMC(vcfRobj = pfcross_subset,
                          PLAF = PLAF,
                          vcfploid = 1, # ploidy of VCF
                          m_max = 5, # max COI to consider
                          rho = 7.4e-7, # recombination rate
                          k_max = 10, # max switch rate to consider
                          e1 = 0.05, # error for going from homozygous to heterozygous
                          e2 = 0.05, # error for going from heterozygous to homozygous
                          burnin = 1e4,
                          samples = 1e4,
                          reportIteration = 1e3,
                          verbose = TRUE,
                          parallelize = TRUE)

ret[6,]
ret$mcmcout[[6]]$summary$quantiles
vcfRmanip::plot_vcfRobj_GT(pfcross_subset)


#............................................................
# temp
#...........................................................

x <- ret$mcmcout[[6]]
# get IBD matrix
CHROM <- x$summary$IBD_marginal[,1]
POS <- x$summary$IBD_marginal[,2]
IBD <- as.matrix(x$summary$IBD_marginal[,-(1:2)])

IBDdf <- cbind.data.frame(CHROM, POS, IBD)
IBDdflong <- IBDdf %>%
  tidyr::pivot_longer(cols=-c("CHROM", "POS"), names_to="Z", values_to="Prob") %>%
  dplyr::mutate(Znum = as.numeric(gsub("z", "", Z))) %>%
  dplyr::group_by(CHROM, Znum) %>%
  dplyr::mutate(start = dplyr::lag(POS),
                end = POS) %>%
  dplyr::ungroup(.)

# filter unneccessary Znumbers
filtdf <- aggregate(IBDdflong$Prob, list(factor(IBDdflong$Znum)), sum)
if(any(filtdf == 0)){
  filtdf <- filtdf[which(filtdf[,2] == 0), 1]
  IBDdflong <- IBDdflong %>% dplyr::filter(! Znum %in% filtdf )
}

library(ggplot2)
ggplot() +
  geom_rect(data = IBDdflong, mapping = aes(xmin = start, xmax = end, ymin = Znum - 0.49, ymax = Znum + 0.49, fill = Prob)) +
  viridis::scale_fill_viridis("IBD Probability", option = "plasma", limits = c(0,1)) +
  scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
  rplasmodium::scale_x_genome(scale = "Mb") +
  xlab("POS") +
  facet_grid(~CHROM) +
  theme_bw()
