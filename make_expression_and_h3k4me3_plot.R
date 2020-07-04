#

library(matrixStats)
library(tximport)
library(EnsDb.Hsapiens.v75)

####
load(file = "median_expr_testis.rda") #this is median expression across individuals in gtex
expr_germline <- log2(median_expr_testis[which(names(median_expr_testis) %in% proms_good_power$gene_id)] + 1)
expr_germline <- sapply(proms_good_power$gene_id, function(xx) expr_germline[which(names(expr_germline) == xx)])
names(expr_germline) <- gsub("[.].*", "", names(expr_germline))

####
files <- paste0("quants/", 
                list.files("quants"), "/quant.sf") #expression from Swanson et al., 2020

names(files) <- paste0("healthy", 1:8)

edb <- EnsDb.Hsapiens.v75
k <- keys(edb, keytype = "GENEID")
df <- select(edb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]
txi <- tximport(files, type = "salmon", tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)

expr_germline2 <- log2(txi$counts + 1)
no_expr_germline2 <- rownames(expr_germline2)[which(rowMedians(expr_germline2) == 0 & rowMaxs(expr_germline2) <= 0.25)]
high_expr_germline2 <- rownames(expr_germline2)[which(rowMedians(expr_germline2) >= 7.5 & rowMins(expr_germline2) >= 5)]

###assess concordance between testis and sperm
quartz(file = "germline_expression_concordance.pdf", width = 4.4, height = 2.2, pointsize = 8, type = "pdf")
par(mfrow = c(1, 2))
plot(density(expr_germline[which(names(expr_germline) %in% no_expr_germline2)], from = 0), lwd = 2.5, 
     col = "forest green", xlab = "expression in testis", main = "", xlim = c(0, 12), bty = 'l', xaxt = 'n', yaxt = 'n')
lines(density(expr_germline[which(names(expr_germline) %in% high_expr_germline2)], from = 0), lwd = 2.5, 
      col = "orange")
legend <- legend("topright", cex = 0.72, bty = 'n', 
                 legend = c("not expressed in sperm", "highly expressed in sperm"), col = c("forest green", "orange"), 
                 lty = "solid", lwd = 2.5)
axis(1, at = c(0, 5, 10), cex.axis = 1.2)
axis(2, at = c(0, 0.5), cex.axis = 1.2)

plot(expr_germline_matched, expr_germline2_matched, 
     col = alpha("forest green", 0.02), pch = 19, cex = 0.5, bty = 'l', 
     main = "", 
     xaxt = 'n', yaxt = 'n', xlab = "testis expression (GTEx)", ylab = "sperm expression", cex.main = 1.5, cex.lab = 1.35)
axis(1, at = c(0, 5, 10), cex.axis = 1.2)
axis(2, at = c(0, 6, 12), cex.axis = 1.2)

dev.off()



#####the following is ran after read mapping to hg19 with bowtie2, deduplication with picard, and peak calling with macs2
mypeaks <- read.delim("h3k4me3_peaks.narrowPeak", header=F)
colnames(mypeaks) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "fold.enrichment",
                       "log10.pval", "log10.qval", "peak")

h3k4me3_chip_granges <- GRanges(seqnames = mypeaks$chrom, IRanges(mypeaks$chromStart+1, mypeaks$chromEnd), 
                                fold.enrichment = mypeaks$fold.enrichment, log10.pval = mypeaks$log10.pval, 
                                log10.qval = mypeaks$log10.qval)


proms_good_power$has_h3k4me3_peak <- "no"
proms_good_power$has_h3k4me3_peak[unique(queryHits(findOverlaps(proms_good_power, h3k4me3_chip_granges)))] <- "yes"


###make plots
quartz(file = "loeuf_meth_expression_germline.pdf", width = 4.4, height = 2.2, pointsize = 8, type = "pdf")
#all genes
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40)], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        xlim = c(0.8, 14.2), ylim = c(0, 2.25), at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "downstream gene LOEUF", bowwex = 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80)], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 2, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

#not expressed in testis
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              expr_germline <= quantile(expr_germline, 0.025))], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        at = 3, xaxt = 'n', main = "", font.main = 1,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              expr_germline <= quantile(expr_germline, 0.025))], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 4, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
#not expressed in sperm
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              proms_good_power$gene_id %in% no_expr_germline2)], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        at = 5, xaxt = 'n', main = "", font.main = 1,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              proms_good_power$gene_id %in% no_expr_germline2)], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 6, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

#highly expressed in testis
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              expr_germline >= quantile(expr_germline, 0.75))], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        at = 7, xaxt = 'n', main = "", font.main = 1,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              expr_germline >= quantile(expr_germline, 0.75))], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 8, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

#highly expressed in sperm
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              proms_good_power$gene_id %in% high_expr_germline2)], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        at = 9, xaxt = 'n', main = "", font.main = 1,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              proms_good_power$gene_id %in% high_expr_germline2)], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 10, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

#now stratify according the presence of H3K4me3 at the promoter in sperm
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              proms_good_power$has_h3k4me3_peak == "no")], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        at = 11, xaxt = 'n', main = "", font.main = 1,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              proms_good_power$has_h3k4me3_peak == "no")], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 12, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              proms_good_power$has_h3k4me3_peak == "yes")], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        at = 13, xaxt = 'n', main = "", font.main = 1,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              proms_good_power$has_h3k4me3_peak == "yes")], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 14, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)
axis(1, at = c(1.5, 4.5, 8.5, 11.5, 13.5), labels = c("all", "not\nexpressed", "highly\nexpressed", "with\nH3K4me3", "w/out\nH3K4me3"), 
     cex.axis = 0.8)
axis(2, at = c(0, 1, 2), las = 2)
legend <- legend("top", legend = c("hypomethylated promoters", "methylated promoters"), 
                 fill = c(alpha("brown", 0.75), "cornflowerblue"), bty = 'n', cex = 0.62, 
                 border = "white", horiz = TRUE)

dev.off()





#####supplemental plot for H3K4me3

###run this in cluster
library(Rsubread)
combined_granges <- proms_good_power
combined_granges$id <- paste0(seqnames(combined_granges), "_", 
                              as.character(start(combined_granges)), "_", as.character(end(combined_granges)))
combined_peaks_ann <- data.frame(GeneID = combined_granges$id, 
                                 Chr = seqnames(combined_granges), 
                                 Start = start(combined_granges), 
                                 End = end(combined_granges), 
                                 strand = strand(combined_granges), 
                                 stringsAsFactors = FALSE)
bamFiles_filepaths <- list.files(pattern = "dedup_merged.bam")
count_reads_in_peaks <- featureCounts(bamFiles_filepaths, annot.ext = combined_peaks_ann, 
                                      isPairedEnd = FALSE, countChimericFragments = FALSE, 
                                      nthreads = 30, countMultiMappingReads = FALSE)


save(count_reads_in_peaks, file = "count_reads_in_proms_h3k4me3.rda")

###and run this locally
load(file = "count_reads_in_proms_h3k4me3.rda")

proms_good_power$h3k4me3_antib_input_ratio <- count_reads_in_peaks$counts[which(
  proms_good_power_init$total_cg_count_with_good_coverage >= 10), 1] / count_reads_in_peaks$counts[
    which(proms_good_power_init$total_cg_count_with_good_coverage >= 10), 2]



quartz(file = "loeuf_meth_h3k4me3_germline_antib_control_ratio.pdf", width = 2.2, height = 2.2, pointsize = 8, type = "pdf")
#par(mar = c(5, 5.4, 0.75, 0.5) + 0.1)
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              proms_good_power$h3k4me3_antib_input_ratio >= quantile(proms_good_power$h3k4me3_antib_input_ratio, 0.9))], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), yaxt = 'n', medlwd = 0.75,
        xlim = c(0.8, 4.2), ylim = c(0, 2.25), at = 1, xaxt = 'n', 
        boxlty = 0, staplelwd = 0, outline = FALSE, ylab = "downstream gene LOEUF", bowwex = 0.1)
#points(jitter(rep(1, length(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40)])), 5), 
#       proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40)], pch = 19, cex= 0.25, 
#       col = rgb(0,0,0,0.01))
boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              proms_good_power$h3k4me3_antib_input_ratio >= quantile(proms_good_power$h3k4me3_antib_input_ratio, 0.9))], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 2, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage <= 40 & 
                                              proms_good_power$h3k4me3_antib_input_ratio <= quantile(proms_good_power$h3k4me3_antib_input_ratio, 0.1))], frame = FALSE, 
        lty = "solid", col = alpha("brown", 0.75), at = 3, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

boxplot(proms_good_power$oe_lof_upper[which(proms_good_power$meth_cg_percentage >= 80 & 
                                              proms_good_power$h3k4me3_antib_input_ratio <= quantile(proms_good_power$h3k4me3_antib_input_ratio, 0.1))], frame = FALSE, 
        lty = "solid", col = "cornflowerblue", at = 4, yaxt = 'n', medlwd = 0.75,
        boxlty = 0, staplelwd = 0, outline = FALSE, add = TRUE, bowwex = 0.1)

legend <- legend("top", legend = c("hypomethylated promoters", "methylated promoters"), 
                 fill = c(alpha("brown", 0.75), "cornflowerblue"), bty = 'n', cex = 0.57, 
                 border = "white")

axis(2, at = c(0, 1, 2), las = 2)
axis(1, at = c(1.5, 3.5), labels = c("low h3k4me3/\ninput ratio", "high h3k4me3/\ninput ratio"), cex.axis = 0.75)

dev.off()





