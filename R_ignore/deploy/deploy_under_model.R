devtools::load_all()
n <- 1e3 # number of loci we would like to simulate
rho <- 7.4e-7 # recombination rate
k_true <- 10 # switch-rate that we "integrate over" as a nuisance parameter
f_true <- 0.4 # proportion of genetic relatedness we are trying to infer
m1 <- 3 # multiplicity of infection for sample-1 that we are trying to infer
m2 <- 3 # multiplicity of infection for sample-2 that we are trying to infer
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

# store trueIBD from simulation between samples
trueIBD <- data.frame(CHROM = vcfR::getCHROM(sim$vcfRobj),
                      POS =vcfR::getPOS(sim$vcfRobj),
                      z_true = rowSums(sim$IBD[,1:ncol(sim$IBD),drop=FALSE]))







# run model
ret <- HMMERTIME::runMCMC(vcfRobj = sim$vcfRobj, # vcfR object we simulated
                          vcfploid = 2, # ploidy of VCF
                          PLAF = PLAF,
                          m_max = 5, # max COI to consider
                          rho = rho, # recombination rate
                          k_max = 10, # max switch rate to consider
                          e1 = 0.025, # error for going from homozygous to heterozygous
                          e2 = 0.025, # error for going from heterozygous to homozygous
                          burnin = 1e4,
                          samples = 1e4,
                          reportIteration = 1e3,
                          verbose = TRUE,
                          parallelize = TRUE)

plot(ret$mcmcout[[1]]$posteriors$logLike)
plot(ret$mcmcout[[1]]$posteriors$f); f_true
plot(ret$mcmcout[[1]]$posteriors$k); k_true
plot(ret$mcmcout[[1]]$posteriors$f_ind)
ret$mcmcout[[1]]$summary$quantiles
mean(unlist(sim$IBD[,1:ncol(sim$IBD)]))
plot(ret$mcmcout[[1]]$posteriors$m1)
plot(ret$mcmcout[[1]]$posteriors$m2)


#......................
# plot out
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
  geom_line(data = trueIBD, aes(x = POS, y = z_true),
            colour = "#38A31A", size = 0.75) +
  viridis::scale_fill_viridis("IBD Probability", option = "plasma", limits = c(0,1)) +
  scale_y_continuous("Number of IBD Genotypes", breaks = seq(1:max(IBDdflong$Znum+1))-1) +
  xlab("POS") +
  facet_grid(~CHROM) +
  theme_bw()

#......................
# quick plots
#......................
vcfRmanip::plot_vcfRobj_GT(sim$vcfRobj)
tibble::tibble(CHROM = vcfR::getCHROM(sim$vcfRobj),
               POS = vcfR::getPOS(sim$vcfRobj)) %>%
  dplyr::bind_cols(., as.data.frame(extract.gt(sim$vcfRobj))) %>%
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
  ggtitle(paste("Sim data for MOI:", paste(m_true, collapse = ","))) +
  facet_grid(CHROM~.) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold", family = "Arial"),
        axis.text.y = element_text(size=12, face="bold", family = "Arial"))

#............................................................
# temp
#...........................................................
ret$mcmcout[[1]]$summary$quantiles
f_true
mean(trueIBD$z_true)
mean(unlist(sim$IBD[,1:ncol(sim$IBD)]))

sum(vcfR::extract.gt(sim$vcfRobj)[, "Sample1"] ==
      vcfR::extract.gt(sim$vcfRobj)[, "Sample2"])/nrow(vcfR::extract.gt(sim$vcfRobj))





#............................................................
#
#...........................................................
head(trueIBD)
zeropos <-trueIBD[trueIBD$z_true == 0, ]

zerogt <- vcfR::extract.gt(sim$vcfRobj)[which(trueIBD$z_true == 0),]
readr::write_csv(as.data.frame(sim$haploid$haploid1[which(trueIBD$z_true == 0),]),
                 "~/Desktop/hap1.csv"
                 )
readr::write_csv(as.data.frame(sim$haploid$haploid2[which(trueIBD$z_true == 0),]),
                 "~/Desktop/hap2.csv"
)
