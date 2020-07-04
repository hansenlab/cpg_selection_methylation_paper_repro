load(file = "proms_with_human_methylation_info.rda")
load(file = file = "tx_stop_sites_with_human_methylation_info.rda")

###define function used for plotting
getPercentOfMethPromoters <- function(promoter_granges, loeuf_quantile, what = c("less", "greater", "in between")){
  if (what == "less"){
    length(which(promoter_granges$meth_cg_percentage[which(promoter_granges$oe_lof_upper <= 
                                                             quantile(promoter_granges$oe_lof_upper, loeuf_quantile[1]))] >= 80)) / 
      length(promoter_granges$meth_cg_percentage[which(promoter_granges$oe_lof_upper <= 
                                                         quantile(promoter_granges$oe_lof_upper, loeuf_quantile[1]))])
  } else if (what == "in between"){
    length(which(promoter_granges$meth_cg_percentage[which(promoter_granges$oe_lof_upper <= 
                                                             quantile(promoter_granges$oe_lof_upper, loeuf_quantile[2]) & 
                                                             promoter_granges$oe_lof_upper > 
                                                             quantile(promoter_granges$oe_lof_upper, loeuf_quantile[1]))] >= 80)) / 
      length(promoter_granges$meth_cg_percentage[which(promoter_granges$oe_lof_upper <= 
                                                         quantile(promoter_granges$oe_lof_upper, loeuf_quantile[2]) & 
                                                         promoter_granges$oe_lof_upper > 
                                                         quantile(promoter_granges$oe_lof_upper, loeuf_quantile[1]))])
  } else {
    length(which(promoter_granges$meth_cg_percentage[which(promoter_granges$oe_lof_upper > 
                                                             quantile(promoter_granges$oe_lof_upper, loeuf_quantile[1]))] >= 80)) / 
      length(promoter_granges$meth_cg_percentage[which(promoter_granges$oe_lof_upper > 
                                                         quantile(promoter_granges$oe_lof_upper, loeuf_quantile[1]))])
  }
}

###obtain null regions
null_permutation_meth_promoter <- replicate(10000, {
  proms_good_power_boot <- proms_good_power[sample(1:length(proms_good_power), length(which(proms_good_power$oe_lof_upper <= 
                                                                                              quantile(proms_good_power$oe_lof_upper, 0.1))))]
  length(which(proms_good_power_boot$meth_cg_percentage >= 80))/length(proms_good_power_boot$meth_cg_percentage)
})

null_permutation_meth_3_prime <- replicate(10000, {
  tx_boot <- tx_all[sample(1:length(tx_all), length(which(tx_all$oe_lof_upper <= quantile(tx_all$oe_lof_upper, 0.1))))]
  length(which(tx_boot$meth_cg_percentage >= 80))/length(tx_boot$meth_cg_percentage)
})

###make plots
quartz(file = "methylation_loeuf_distribution.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 1) + 0.1)
plot(1, getPercentOfMethPromoters(proms_good_power, 0.1, "less"), pch = 19, col = alpha("red", 0.59), 
     bty = 'l', xlab = "LOEUF decile", cex = 1, ylab = "% of genes with\nmethylated promoter", main = "", xlim = c(0.8, 10.2), 
     ylim = c(0, 0.3), yaxt = 'n', xaxt= 'n')
rect(0, quantile(null_permutation_meth_promoter, 0.0025), 12, quantile(null_permutation_meth_promoter, 1 - 0.0025), 
     col = alpha("gray", 0.4), lty = 0)
points(2, getPercentOfMethPromoters(proms_good_power, c(0.1, 0.2), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(3, getPercentOfMethPromoters(proms_good_power, c(0.2, 0.3), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(4, getPercentOfMethPromoters(proms_good_power, c(0.3, 0.4), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(5, getPercentOfMethPromoters(proms_good_power, c(0.4, 0.5), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(6, getPercentOfMethPromoters(proms_good_power, c(0.5, 0.6), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(7, getPercentOfMethPromoters(proms_good_power, c(0.6, 0.7), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(8, getPercentOfMethPromoters(proms_good_power, c(0.7, 0.8), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(9, getPercentOfMethPromoters(proms_good_power, c(0.8, 0.9), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(10, getPercentOfMethPromoters(proms_good_power, 0.9, "greater"), pch = 19, col = alpha("red", 0.59), cex = 1)
axis(1, at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), cex.axis = 0.72)
axis(2, at = c(0, 0.15, 0.30), labels = c("0", "15", "30"))
dev.off()

#
quartz(file = "methylation_loeuf_distribution_3_prime_end.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 1) + 0.1)
plot(1, getPercentOfMethPromoters(tx_all, 0.1, "less"), pch = 19, col = alpha("red", 0.59), 
     bty = 'l', xlab = "LOEUF decile", cex = 1, ylab = "% of genes with\nmethylated 3' end", main = "", xlim = c(0.8, 10.2), 
     ylim = c(0.8, 1), yaxt = 'n', xaxt= 'n')
rect(0, quantile(null_permutation_meth_3_prime, 0.0025), 12, quantile(null_permutation_meth_3_prime, 1 - 0.0025), 
     col = alpha("gray", 0.4), lty = 0)
points(2, getPercentOfMethPromoters(tx_all, c(0.1, 0.2), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(3, getPercentOfMethPromoters(tx_all, c(0.2, 0.3), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(4, getPercentOfMethPromoters(tx_all, c(0.3, 0.4), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(5, getPercentOfMethPromoters(tx_all, c(0.4, 0.5), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(6, getPercentOfMethPromoters(tx_all, c(0.5, 0.6), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(7, getPercentOfMethPromoters(tx_all, c(0.6, 0.7), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(8, getPercentOfMethPromoters(tx_all, c(0.7, 0.8), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(9, getPercentOfMethPromoters(tx_all, c(0.8, 0.9), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(10, getPercentOfMethPromoters(tx_all, 0.9, "greater"), pch = 19, col = alpha("red", 0.59), cex = 1)
axis(1, at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), cex.axis = 0.72)
axis(2, at = c(0.8, 1), labels = c("80", "100"))
dev.off()


###make same plot for promoter using s_het instead of loeuf
library(readr)
s_het_df <- read_csv('s_het_estimates.csv') #from Cassa et al., 2017
proms_good_power2 <- proms_good_power
proms_good_power2$oe_lof_upper_true <- proms_good_power2$oe_lof_upper
proms_good_power2$s_het <- as.numeric(sapply(proms_good_power2$gene_name, function(xx) 
  s_het_df$s_het[which(s_het_df$gene_symbol == xx)]))
proms_good_power2$oe_lof_upper <- proms_good_power2$s_het #this is just for convenience so that I don't have to change the code I use for the plots
proms_good_power2 <- proms_good_power2[-which(is.na(proms_good_power2$s_het))]


null_permutation_meth_promoter_s_het <- replicate(10000, {
  proms_good_power_boot <- proms_good_power2[sample(1:length(proms_good_power2), length(which(proms_good_power2$oe_lof_upper <= 
                                                                                                quantile(proms_good_power2$oe_lof_upper, 0.1))))]
  length(which(proms_good_power_boot$meth_cg_percentage >= 80))/length(proms_good_power_boot$meth_cg_percentage)
})

quartz(file = "methylation_s_het_distribution.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
par(mar = c(4, 5, 1, 1) + 0.1)
plot(1, getPercentOfMethPromoters(proms_good_power2, 0.1, "less"), pch = 19, col = alpha("red", 0.59), 
     bty = 'l', xlab = "s_het decile", cex = 1, ylab = "% of genes with\nmethylated promoter", main = "", xlim = c(0.8, 10.2), 
     ylim = c(0, 0.3), yaxt = 'n', xaxt= 'n')
rect(0, quantile(null_permutation_meth_promoter_s_het, 0.0025), 12, quantile(null_permutation_meth_promoter_s_het, 1 - 0.0025), 
     col = alpha("gray", 0.4), lty = 0)
points(2, getPercentOfMethPromoters(proms_good_power2, c(0.1, 0.2), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(3, getPercentOfMethPromoters(proms_good_power2, c(0.2, 0.3), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(4, getPercentOfMethPromoters(proms_good_power2, c(0.3, 0.4), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(5, getPercentOfMethPromoters(proms_good_power2, c(0.4, 0.5), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(6, getPercentOfMethPromoters(proms_good_power2, c(0.5, 0.6), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(7, getPercentOfMethPromoters(proms_good_power2, c(0.6, 0.7), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(8, getPercentOfMethPromoters(proms_good_power2, c(0.7, 0.8), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(9, getPercentOfMethPromoters(proms_good_power2, c(0.8, 0.9), "in between"), pch = 19, col = alpha("red", 0.59), cex = 1)
points(10, getPercentOfMethPromoters(proms_good_power2, 0.9, "greater"), pch = 19, col = alpha("red", 0.59), cex = 1)
axis(1, at = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), cex.axis = 0.72)
axis(2, at = c(0, 0.15, 0.30), labels = c("0", "15", "30"))
dev.off()






###finally, assess the rank correlation between loeuf and s_het
quartz(file = "loeuf_s_het_scatterplot.pdf", height = 2.2, width = 2.2, pointsize = 8, type = "pdf")
plot(proms_good_power2$oe_lof_upper_true, proms_good_power2$s_het, 
     col = rgb(1,0,0,0.035), pch = 19, cex = 0.5, bty = 'l', 
     main = paste0("rank correlation = ", 
                   round(cor(proms_good_power2$oe_lof_upper_true, proms_good_power2$s_het, method = "spearman"), 2)), 
     xaxt = 'n', yaxt = 'n', ylab = "s_het", xlab = "LOEUF", cex.main = 1.35, cex.lab = 1.35, font.main = 1)
axis(1, at = c(0, 1, 2), cex.axis = 1.2)
axis(2, at = c(0, 0.35, 0.7), labels = c(0, 0.35, 0.7), cex.axis = 1.2)
dev.off()


