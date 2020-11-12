devtools::load_all()
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
plot(ret$mcmcout[[6]]$posteriors$logLike)
vcfRmanip::plot_vcfRobj_GT(pfcross_subset)


#............................................................
# temp
#...........................................................

x <- ret$mcmcout[[3]]
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



#............................................................
# hmmIBD test
#...........................................................
library(tidyverse)
dir.create("~/Desktop/hmmIBD_test/")
ibd <- tibble::tibble(CHROM = vcfR::getCHROM(pfcross_subset),
                      POS = vcfR::getPOS(pfcross_subset)) %>%
  dplyr::mutate(CHROM = factor(CHROM, levels = rplasmodium::chromnames()),
                CHROM = as.integer(CHROM))
ibdgt <- vcfRmanip::gtmat012(pfcross_subset)
ibdgt[is.na(ibdgt)] <- -1
ibdgt[ibdgt == 2] <- 1
summary(ibdgt)

ibdout <- dplyr::bind_cols(ibd, as.data.frame(ibdgt))

write.table(x = ibdout, file = "~/Desktop/hmmIBD_test/pfcross_ibd.txt")

#plaf
#set.seed(48)
#PLAF <- rbeta(nrow(ibdout), 0.1, 0.1)
plafout <- dplyr::bind_cols(ibd, tibble::tibble(
  a0 = PLAF,
  a1 = 1-PLAF))

write_tsv(x = plafout, file = "~/Desktop/hmmIBD_test/pfcross_plaf.txt")


system("/Users/nbrazeau/Documents/GitHub/hmmIBD/hmmIBD -i /Users/nbrazeau/Desktop/hmmIBD_test/pfcross_ibd.txt -f /Users/nbrazeau/Desktop/hmmIBD_test/pfcross_plaf.txt -o /Users/nbrazeau/Desktop/hmmIBD_test/mytest_pfcross")

#......................
# read in results
#.....................
hmmibdret <- readr::read_tsv("~/Desktop/hmmIBD_test/mytest_pfcross.hmm_fract.txt")
hmmibd_segs <- readr::read_tsv("~/Desktop/hmmIBD_test/mytest_pfcross.hmm.txt")

hmmibd_segs %>%
  dplyr::filter(grepl("HB3", sample1)) %>%
  dplyr::filter(grepl("C02", sample2))

hmmibdret %>%
  dplyr::filter(grepl("HB3", sample1)) %>%
  dplyr::filter(grepl("C02", sample2)) %>%
  dplyr::select(c("sample1", "sample2", "N_generation", "fract_sites_IBD"))


which(diff(ret$mcmcout[[5]]$summary$IBD_marginal$z0) == max(diff(ret$mcmcout[[5]]$summary$IBD_marginal$z0)))
ret$mcmcout[[5]]$summary$IBD_marginal[176:180,]
ret$mcmcout[[5]]$summary$quantiles

sum(vcfR::extract.gt(pfcross_subset)[, "HB3/PG0052-C/ERR019054"] ==
      vcfR::extract.gt(pfcross_subset)[, "C02/PG0053-C/ERR019067"])/nrow(vcfR::extract.gt(pfcross_subset))

1 - 177/nrow(ret$mcmcout[[5]]$summary$IBD_marginal)
colSums(ret$mcmcout[[5]]$summary$IBD_marginal[3:7])/nrow(vcfR::extract.gt(pfcross_subset))
