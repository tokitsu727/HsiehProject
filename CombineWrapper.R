library('VariantCombiner')
args <- commandArgs(trailingOnly = TRUE)

writeVcf(SomaticCombiner(args[1], args[2], c(args[3], args[4])), args[5])
