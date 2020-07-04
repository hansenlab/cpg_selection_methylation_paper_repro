####
load(file = "proms_with_human_methylation_info.rda")
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
cpg_ranges <- cpgs_in_proms_both_strands

proms_no_cgs <- setdiff(proms_good_power_init, cpg_ranges, ignore.strand = TRUE) #this excludes CpG sites
overlaps <- findOverlaps(proms_no_cgs, proms_good_power_init)
proms_no_cgs$tx_id <- proms_good_power_init$tx_id[subjectHits(overlaps)]

###add PhastCons info
library(phastCons100way.UCSC.hg19)
gsco <- phastCons100way.UCSC.hg19
###
proms_no_cgs$score <- gscores(gsco, proms_no_cgs)$default

regions_by_tx <- split(proms_no_cgs, proms_no_cgs$tx_id)
regions_by_tx <- regions_by_tx[proms_good_power_init$tx_id]

proms_good_power_init$score_no_cpgs <- as.numeric(sapply(regions_by_tx, function(xx) {
  gr <- xx
  widths <- width(gr)
  weighted_average_numerator <- sum(gr$score*widths)
  weighted_average_denominator <- sum(widths)
  weighted_average_numerator/weighted_average_denominator
}))

cpgs_in_proms_both_strands$score <- gscores(gsco, cpgs_in_proms_both_strands)$default

cpgs_hypometh <- cpgs_in_proms_both_strands[which(cpgs_in_proms_both_strands$meth_level <= 0.2 & 
                                                    cpgs_in_proms_both_strands$coverage >= 10)]
cpgs_meth <- cpgs_in_proms_both_strands[which(cpgs_in_proms_both_strands$meth_level >= 0.8 & 
                                                cpgs_in_proms_both_strands$coverage >= 10)]

cpgs_hypometh_by_tx <- split(cpgs_hypometh, cpgs_hypometh$tx_id)
cpgs_hypometh_by_tx <- cpgs_hypometh_by_tx[proms_good_power_init$tx_id]

proms_good_power_init$score_cpgs <- as.numeric(sapply(cpgs_hypometh_by_tx, function(xx) {
  mean(xx$score)
}))


quartz(file = "CpGs_vs_non_CpGs_phastcons.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(6, 4, 1, 1)+0.2)
proms_good_power_init$n_cpgs <- lengths(cpgs_hypometh_by_tx)
proms_good_power_init2 <- proms_good_power_init[which(proms_good_power_init$n_cpgs >= 30)]
boxplot(proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper <= 0.35)], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n',
        xlim = c(0.8, 8.2), medlty = 1, medlwd = 0.8, ylim = c(0, 0.75), at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "PhastCons")
boxplot(proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper <= 0.35)], 
        frame = FALSE, lty = "solid", col = alpha("deep pink", 0.5), 
        medlty = 1, medlwd = 0.8, at = 2, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
#
boxplot(proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.5) & 
                                                  proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.25))], 
        frame = FALSE, lty = "solid", col = alpha("brown", 0.75), 
        medlty = 1, medlwd = 0.8, at = 3, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
boxplot(proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.5) & 
                                                     proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.25))], 
        frame = FALSE, lty = "solid", col = alpha("deep pink", 0.5), 
        medlty = 1, medlwd = 0.8, at = 4, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
#
boxplot(proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.75) & 
                                                  proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.5))], 
        frame = FALSE, lty = "solid", col = alpha("brown", 0.75), 
        medlty = 1, medlwd = 0.8, at = 5, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
boxplot(proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.75) & 
                                                     proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.5))], 
        frame = FALSE, lty = "solid", col = alpha("deep pink", 0.5), 
        medlty = 1, medlwd = 0.8, at = 6, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
#
boxplot(proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper >= 1)], 
        frame = FALSE, lty = "solid", col = alpha("brown", 0.75), 
        medlty = 1, medlwd = 0.8, at = 7, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)

boxplot(proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper >= 1)], 
        frame = FALSE, lty = "solid", col = alpha("deep pink", 0.5), 
        medlty = 1, medlwd = 0.8, at = 8, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
axis(2, at = c(0, 0.3, 0.6))
axis(1, at = c(1.5, 3.5, 5.5, 7.5), labels = c("1", "2", "3", "4"), las = 2)
legend <- legend("topright", legend = c("promoter\nhypomethylated CpGs", "promoter non-CpGs"), 
                 fill = c(alpha("brown", 0.75), alpha("deep pink", 0.5)), 
                 border = "white", bty = 'n', cex = 0.75)
dev.off()


quartz(file = "CpGs_vs_non_CpGs_phastcons_median_diff.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
diff_1 <- proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.25))]-
  proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.25))]

diff_2 <- proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.5) & 
                                                    proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.25))] - 
  proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.5) & 
                                               proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.25))]

diff_3 <- proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.75) & 
                                                    proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.5))] - 
  proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper <= quantile(proms_good_power_init2$oe_lof_upper, 0.75) & 
                                               proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.5))]

diff_4 <- proms_good_power_init2$score_cpgs[which(proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.75))]-
  proms_good_power_init2$score_no_cpgs[which(proms_good_power_init2$oe_lof_upper > quantile(proms_good_power_init2$oe_lof_upper, 0.75))]

par(mar = c(4, 5, 1, 1)+0.2)
boxplot(diff_1, frame = FALSE, 
        lty = "solid", col = alpha("red", 0.62), yaxt = 'n',
        xlim = c(0.8, 4.2), medlty = 1, medlwd = 0.8, ylim = c(-0.2, 0.42), at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "average PhastCons difference\n(CpGs vs non-CpGs)")
boxplot(diff_2, 
        frame = FALSE, lty = "solid", col = alpha("red", 0.62), 
        medlty = 1, medlwd = 0.8, at = 2, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
boxplot(diff_3, 
        frame = FALSE, lty = "solid", col = alpha("red", 0.62), 
        medlty = 1, medlwd = 0.8, at = 3, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)
boxplot(diff_4, 
        frame = FALSE, lty = "solid", col = alpha("red", 0.62), 
        medlty = 1, medlwd = 0.8, at = 4, yaxt = 'n',
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE)


axis(2, at = c(-0.2, 0, 0.2, 0.4))
axis(1, at = c(1, 2, 3, 4))

abline(h = 0, lty = "dashed", col = rgb(0,0,0,0.7))
dev.off()
