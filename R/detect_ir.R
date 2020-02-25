irDetect <- function(genome){
  genome_rc <- Biostrings::reverseComplement(genome)
  ir <- genome_rc[8000:12000]
  t <- matchPattern(ir, genome, max.mismatch = 100, with.indels = TRUE)
}
