gc_count <- function(genome, view.width = 10){
  g <- paste(genome, collapse = "")
  g <- Biostrings::DNAString(g)
  n <- length(g)
  g <- c(g, g[1:view.width - 1])
  count <- Biostrings::letterFrequencyInSlidingView(g, letters = c("C", "G"),
                                        view.width = view.width)
  count <- count[, 1] + count[, 2]
  total <- Biostrings::letterFrequency(g, letters = c("C", "G"))
  total <- (total[1] + total[2])/n
  res <- data.frame(position = seq(1, n),
                    gc_count = count/view.width,
                    stringsAsFactors = FALSE)
  res <- res[seq(1, n, by = view.width), ]
  res <- list(gc_count_table = res, total_count = total)
  return(res)
}

