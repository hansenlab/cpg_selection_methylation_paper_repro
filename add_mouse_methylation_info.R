load(file = "proms_with_human_methylation_info.rda")

#load meth level and coverage for mouse CpGs mapped to hg19 genome-wide. Both are downloaded from ucsc
cpgs_with_meth_level_mouse <- import("~/Downloads/mouse_germline_meth.bw", 
                                     format = "BigWig")
cpgs_with_coverage_mouse <- import("~/Downloads/mouse_germline_meth_coverage.bw", 
                                   format = "BigWig")

cpgs_with_meth_level_mouse$meth_level <- cpgs_with_meth_level_mouse$score
cpgs_with_meth_level_mouse$coverage <- cpgs_with_coverage_mouse$score

overlaps <- findOverlaps(cpgs_with_meth_level_mouse, proms_good_power)
cpgs_mouse_human_proms <- cpgs_with_meth_level_mouse[queryHits(overlaps)]
cpgs_mouse_human_proms$tx_id <- proms_good_power$tx_id[subjectHits(overlaps)]
cpgs_mouse_human_proms <- cpgs_mouse_human_proms[which(cpgs_mouse_human_proms$meth_level >= 0.8 | cpgs_mouse_human_proms$meth_level <= 0.2)]
cpgs_mouse_human_proms <- cpgs_mouse_human_proms[which(cpgs_mouse_human_proms$coverage >= 10)]

cpgs_mouse_human_proms_list <- split(cpgs_mouse_human_proms, cpgs_mouse_human_proms$tx_id)


meth_percentage_mouse <- sapply(cpgs_mouse_human_proms_list, function(xx) {
  meth_cg_count <- length(which(xx$meth_level >= 0.8 & xx$coverage >= 10))
  total_cg_count <- length(which(xx$meth_level >= 0.8 & xx$coverage >= 10)) + length(which(xx$meth_level <= 0.2 & xx$coverage >= 10))
  100*(meth_cg_count/total_cg_count)
})

total_cg_count_with_good_coverage_mouse <- sapply(cpgs_mouse_human_proms_list, function(xx) {
  length(which(xx$meth_level >= 0.8 & xx$coverage >= 10)) + length(which(xx$meth_level <= 0.2 & xx$coverage >= 10))
})

proms_good_power$total_cg_count_with_good_coverage_mouse <- as.numeric(sapply(proms_good_power$tx_id, function(xx) 
  total_cg_count_with_good_coverage_mouse[which(names(total_cg_count_with_good_coverage_mouse) == xx)]))


proms_good_power$meth_cg_percentage_mouse <- as.numeric(sapply(proms_good_power$tx_id, function(xx) 
  meth_percentage_mouse[which(names(meth_percentage_mouse) == xx)]))

proms_good_power_for_comparison <- proms_good_power[which(proms_good_power$total_cg_count_with_good_coverage_mouse >= 10)]

save(proms_good_power_for_comparison, file = "proms_with_human_and_mouse_methylation_info.rda")

