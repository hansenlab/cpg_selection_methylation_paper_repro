########
###the following are ran on the cluster
load(file = "proms_good_power_exp_lof_more_than_10.rda")

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

tx_all <- transcripts(edb, columns = c("tx_id"))
tx_all <- tx_all[proms_good_power$tx_id]
tx_all <- tx_all+2000
tx_all <- resize(tx_all, 4000, fix = "end")
tx_all$oe_lof_upper <- proms_good_power$oe_lof_upper
seqlevelsStyle(tx_all) <- "ucsc"
genome(tx_all) <- "hg19"
seqlevels(tx_all) <- seqlevels(tx_all)[1:22]
overlaps <- findOverlaps(cpgs_with_meth_level, tx_all)
cpgs_tx_ends_good_power <- cpgs_with_meth_level[unique(queryHits(overlaps))]


df <- data.frame(chr = seqnames(cpgs_tx_ends_good_power), 
                 start = start(cpgs_tx_ends_good_power), 
                 end = end(cpgs_tx_ends_good_power), 
                 meth_level = cpgs_tx_ends_good_power$meth_level, 
                 coverage = cpgs_tx_ends_good_power$coverage)

save(df, file = "/dcl01/hansen/data/bsseq_germline_molaro_etal/cpg_level_meth_track_from_ucsc/cpgs_tx_ends_good_power_df.rda")

###now these are ran locally
load(file = "cpgs_tx_ends_good_power_df.rda")

cpgs_in_tx_ends <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
overlaps <- findOverlaps(cpgs_in_tx_ends, tx_all)

#next line removes CpGs overlapping more than 1 transcript
cpgs_in_tx_ends <- cpgs_in_tx_ends[-queryHits(overlaps)[which(queryHits(overlaps) %in% 
                                                                queryHits(overlaps)[which(duplicated(queryHits(overlaps)))])]]

#add tx_id info
overlaps <- findOverlaps(cpgs_in_tx_ends, tx_all)
cpgs_in_tx_ends$tx_id <- tx_all$tx_id[subjectHits(overlaps)]
genome(cpgs_in_tx_ends) <- "hg19"

cpgs_in_tx_ends_list <- split(cpgs_in_tx_ends, cpgs_in_tx_ends$tx_id)
tx_all <- tx_all[names(cpgs_in_tx_ends_list)]

meth_percentage <- sapply(cpgs_in_tx_ends_list, function(xx) {
  meth_cg_count <- length(which(xx$meth_level >= 0.8 & xx$coverage >= 10))
  total_cg_count <- length(which(xx$meth_level >= 0.8 & xx$coverage >= 10)) + length(which(xx$meth_level <= 0.2 & xx$coverage >= 10))
  100*(meth_cg_count/total_cg_count)
})

total_cg_count_with_good_coverage <- sapply(cpgs_in_tx_ends_list, function(xx) {
  length(which(xx$meth_level >= 0.8 & xx$coverage >= 10)) + length(which(xx$meth_level <= 0.2 & xx$coverage >= 10))
})

tx_all$total_cg_count_with_good_coverage <- sapply(tx_all$tx_id, function(xx) 
  total_cg_count_with_good_coverage[which(names(total_cg_count_with_good_coverage) == xx)])

tx_all$meth_cg_percentage <- sapply(tx_all$tx_id, function(xx) 
  meth_percentage[which(names(meth_percentage) == xx)])

tx_all_init <- tx_all
tx_all <- tx_all[which(tx_all$total_cg_count_with_good_coverage >= 10)]
tx_all$meth_cg_percentage <- as.numeric(tx_all$meth_cg_percentage)
save(tx_all_init, tx_all, file = "tx_stop_sites_with_human_methylation_info.rda")



