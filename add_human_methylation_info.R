###the following are ran on the cluster
library(rtracklayer)

#load meth level and coverage for human CpGs genome-wide. Both are downloaded from ucsc
cpgs_with_meth_level <- import("/dcl01/hansen/data/bsseq_germline_molaro_etal/cpg_level_meth_track_from_ucsc/Human_Sperm.meth.bw", 
                               format = "BigWig")
cpgs_with_coverage <- import("/dcl01/hansen/data/bsseq_germline_molaro_etal/cpg_level_meth_track_from_ucsc/Human_Sperm.read.bw", 
                             format = "BigWig")

cpgs_with_meth_level$meth_level <- cpgs_with_meth_level$score
cpgs_with_meth_level$coverage <- cpgs_with_coverage$score

load(file = "proms_good_power_exp_lof_more_than_10.rda")
overlaps <- findOverlaps(cpgs_with_meth_level, proms_good_power)
cpgs_proms_good_power <- cpgs_with_meth_level[unique(queryHits(overlaps))]

df <- data.frame(chr = seqnames(cpgs_proms_good_power), 
                 start = start(cpgs_proms_good_power), 
                 end = end(cpgs_proms_good_power), 
                 meth_level = cpgs_proms_good_power$meth_level, 
                 coverage = cpgs_proms_good_power$coverage)

save(df, file = "/dcl01/hansen/data/bsseq_germline_molaro_etal/cpg_level_meth_track_from_ucsc/cpgs_proms_good_power_df.rda")


###from now on everything is ran locally
load(file = "cpgs_proms_good_power_df.rda")

cpgs_in_proms <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
overlaps <- findOverlaps(cpgs_in_proms, proms_good_power)
cpgs_in_proms$tx_id <- proms_good_power$tx_id[subjectHits(overlaps)]
genome(cpgs_in_proms) <- "hg19"

cpgs_in_proms_list <- split(cpgs_in_proms, cpgs_in_proms$tx_id)

meth_percentage <- sapply(cpgs_in_proms_list, function(xx) {
  meth_cg_count <- length(which(xx$meth_level >= 0.8 & xx$coverage >= 10))
  total_cg_count <- length(which(xx$meth_level >= 0.8 & xx$coverage >= 10)) + length(which(xx$meth_level <= 0.2 & xx$coverage >= 10))
  100*(meth_cg_count/total_cg_count)
})

total_cg_count_with_good_coverage <- sapply(cpgs_in_proms_list, function(xx) {
  length(which(xx$meth_level >= 0.8 & xx$coverage >= 10)) + length(which(xx$meth_level <= 0.2 & xx$coverage >= 10))
})

proms_good_power$total_cg_count_with_good_coverage <- sapply(proms_good_power$tx_id, function(xx) 
  total_cg_count_with_good_coverage[which(names(total_cg_count_with_good_coverage) == xx)])

proms_good_power$meth_cg_percentage <- sapply(proms_good_power$tx_id, function(xx) 
  meth_percentage[which(names(meth_percentage) == xx)])

proms_good_power_init <- proms_good_power
proms_good_power <- proms_good_power[which(proms_good_power$total_cg_count_with_good_coverage >= 10)]
proms_good_power$meth_cg_percentage <- as.numeric(proms_good_power$meth_cg_percentage)

save(proms_good_power_init, proms_good_power, file = "proms_with_human_methylation_info.rda")




