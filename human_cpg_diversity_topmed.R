##define function that estimates nucleotide diversity
getAverageHeterozygosity <- function(cpg_genomic_ranges_maf, bootstrap = c("yes", "no")){
  n_alleles <- 2*62874 #that's the number of individuals in TOPMed
  
  cpg_maf <- cpg_genomic_ranges_maf
  total_number_of_sites <- length(cpg_maf)
  if (bootstrap == "yes"){
    cpg_maf <- cpg_maf[sample(1:length(cpg_maf), length(cpg_maf), replace = TRUE)]
  }
  
  sum(cpg_maf$AF*(1-cpg_maf$AF)*(n_alleles/(n_alleles-1))/total_number_of_sites)
}


##load(file = "all_snps_in_proms_good_power_granges_df.rda")
####CpG snps
cg_snps <- makeGRangesFromDataFrame(cg_df, keep.extra.columns = TRUE) 
cg_snps$ref_af <- as.numeric(cg_snps$ref_af)
cg_snps$alt_af <- as.numeric(cg_snps$alt_af)
all_cg_snps_hg19_in_proms <- cg_snps[which(nchar(cg_snps$ref) == 1 & nchar(cg_snps$alt) == 1 & cg_snps$ref_af >= 0.5)]

##
load(file = "cpgs_proms_good_power_df.rda")
cpgs_in_proms <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
overlaps <- findOverlaps(cpgs_in_proms, proms_good_power_init)
cpgs_in_proms$tx_id <- proms_good_power_init$tx_id[subjectHits(overlaps)]
genome(cpgs_in_proms) <- "hg19"

cpgs_in_proms_2 <- makeGRangesFromDataFrame(data.frame(chr = seqnames(cpgs_in_proms), 
                                                       start = start(cpgs_in_proms)+1, 
                                                       end = end(cpgs_in_proms)+1, 
                                                       meth_level = cpgs_in_proms$meth_level, 
                                                       coverage = cpgs_in_proms$coverage, 
                                                       tx_id = cpgs_in_proms$tx_id), 
                                            keep.extra.columns = TRUE)


cpgs_in_proms_both_strands <- unlist(GRangesList(cpgs_in_proms, cpgs_in_proms_2))
#

#
overlaps <- findOverlaps(cpgs_in_proms_both_strands, all_cg_snps_hg19_in_proms)
cpgs_in_proms_both_strands$AF <- 0
cpgs_in_proms_both_strands$AF[queryHits(overlaps)] <- all_cg_snps_hg19_in_proms$alt_af[subjectHits(overlaps)]


cpgs_hypometh <- cpgs_in_proms_both_strands[which(cpgs_in_proms_both_strands$meth_level <= 0.2 & 
                                                    cpgs_in_proms_both_strands$coverage >= 10)]
cpgs_meth <- cpgs_in_proms_both_strands[which(cpgs_in_proms_both_strands$meth_level >= 0.8 & 
                                                cpgs_in_proms_both_strands$coverage >= 10)]


###the following two plots serve as QC. methylated CpGs should be more variable in the human population
quartz(file = "cpg_maf_spectrum.pdf", height = 2.2, width = 4.4, pointsize = 8, type = "pdf")
par(mfrow = c(1,2))
hist(log10(cpgs_hypometh$AF), col = alpha("brown", 0.75), lty = 0, 
     breaks = 50, freq = FALSE, xlab = "log10(minor allele freq)", cex.lab = 1, ylim = c(0, 4.2),
     main = "hypomethylated\nCpGs", font.main = 1, cex.main = 1, xaxt = 'n', yaxt = 'n')
axis(2, at = c(0, 2, 4))
axis(1, at = c(-5, -1))
hist(log10(cpgs_meth$AF), col = "cornflowerblue", lty = 0, 
     breaks = 50, freq = FALSE, xlab = "log10(minor allele freq)", cex.lab = 1, ylim = c(0, 4.2),
     yaxt = 'n', main = "methylated\nCpGs", font.main = 1, cex.main = 1, xaxt = 'n')
axis(2, at = c(0, 2, 4))
axis(1, at = c(-5, -1))

dev.off()


quartz(file = "all_cpgs_nucleotide_diversity.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 4, 1, 1) + 0.1)
plot(1, getAverageHeterozygosity(cpgs_hypometh, "no"), 
     pch = 19, col = alpha("brown", 0.75), 
     bty = 'l', xlab = "", cex = 2, ylab = "nucleotide diversity", main = "", ylim = c(0, 0.0012), 
     xlim = c(0.8, 2.2), yaxt = 'n', xaxt= 'n')

points(2, getAverageHeterozygosity(cpgs_meth, "no"), 
       pch = 19, col = "cornflowerblue", cex = 2)

axis(1, at = c(1, 2), labels = c("hypometh\nCpGs", "meth\nCpGs"))
axis(2, at = c(0, 0.0005, 0.001), labels = c("0", "0.0005", "0.001"))
dev.off()






#####now test how strong methylation is versus other forces on CpGs
quartz(file = "cpg_nucleotide_diversity.pdf", width = 1.8, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 4, 1, 1) + 0.1)
plot(1, getAverageHeterozygosity(cpgs_hypometh[which(cpgs_hypometh$tx_id %in% 
                                                       proms_good_power_init$tx_id[which(proms_good_power_init$oe_lof_upper <= 
                                                                                           quantile(proms_good_power_init$oe_lof_upper, 0.25))])], "no"), 
     pch = 19, col = alpha("brown", 0.75), 
     bty = 'l', xlab = "LOEUF quartile", cex = 0.75, ylab = "nucleotide diversity", main = "", ylim = c(0, 0.00125), 
     xlim = c(0.8, 8.2), yaxt = 'n', xaxt= 'n')

points(3, getAverageHeterozygosity(cpgs_hypometh[which(cpgs_hypometh$tx_id %in% 
                                                         proms_good_power$tx_id[which(proms_good_power$oe_lof_upper <= 
                                                                                        quantile(proms_good_power$oe_lof_upper, 0.50) & 
                                                                                        proms_good_power$oe_lof_upper > 
                                                                                        quantile(proms_good_power$oe_lof_upper, 0.25))])], "no"), 
       pch = 19, col = alpha("brown", 0.75), cex = 0.75)

points(5, getAverageHeterozygosity(cpgs_hypometh[which(cpgs_hypometh$tx_id %in% 
                                                         proms_good_power$tx_id[which(proms_good_power$oe_lof_upper <= 
                                                                                        quantile(proms_good_power$oe_lof_upper, 0.75) & 
                                                                                        proms_good_power$oe_lof_upper > 
                                                                                        quantile(proms_good_power$oe_lof_upper, 0.50))])], "no"), 
       pch = 19, col = alpha("brown", 0.75), cex = 0.75)

points(7, getAverageHeterozygosity(cpgs_hypometh[which(cpgs_hypometh$tx_id %in% 
                                                         proms_good_power$tx_id[which(proms_good_power$oe_lof_upper > 
                                                                                        quantile(proms_good_power$oe_lof_upper, 0.75))])], "no"), 
       pch = 19, col = alpha("brown", 0.75), cex = 0.75)





###methylated
points(2, getAverageHeterozygosity(cpgs_meth[which(cpgs_meth$tx_id %in% 
                                                     proms_good_power$tx_id[which(proms_good_power$oe_lof_upper <= 
                                                                                    quantile(proms_good_power$oe_lof_upper, 0.25))])], "no"), 
       pch = 19, col = "cornflowerblue", cex = 0.75)

points(4, getAverageHeterozygosity(cpgs_meth[which(cpgs_meth$tx_id %in% 
                                                     proms_good_power$tx_id[which(proms_good_power$oe_lof_upper <= 
                                                                                    quantile(proms_good_power$oe_lof_upper, 0.50) & 
                                                                                    proms_good_power$oe_lof_upper > 
                                                                                    quantile(proms_good_power$oe_lof_upper, 0.25))])], "no"), 
       pch = 19, col = "cornflowerblue", cex = 0.75)

points(6, getAverageHeterozygosity(cpgs_meth[which(cpgs_meth$tx_id %in% 
                                                     proms_good_power$tx_id[which(proms_good_power$oe_lof_upper <= 
                                                                                    quantile(proms_good_power$oe_lof_upper, 0.70) & 
                                                                                    proms_good_power$oe_lof_upper > 
                                                                                    quantile(proms_good_power$oe_lof_upper, 0.50))])], "no"), 
       pch = 19, col = "cornflowerblue", cex = 0.75)

points(8, getAverageHeterozygosity(cpgs_meth[which(cpgs_meth$tx_id %in% 
                                                     proms_good_power$tx_id[which(proms_good_power$oe_lof_upper > 
                                                                                    quantile(proms_good_power$oe_lof_upper, 0.75))])], "no"), 
       pch = 19, col = "cornflowerblue", cex = 0.75)

axis(1, at = c(1.5, 3.5, 5.5, 7.5), labels = c("1", "2", "3", "4"))
axis(2, at = c(0, 0.0005, 0.001), labels = c("0", "0.0005", "0.001"))
abline(v = c(2.5, 4.5, 6.5), lty = "longdash", col = rgb(0,0,0,0.7))
legend <- legend("bottom", legend = c("hypomethylated\npromoter CpGs", "methylated\npromoter CpGs"), bty = 'n',
                 pch = 19, col = c(alpha("brown", 0.75), "cornflowerblue"), cex = 0.62)
dev.off()








