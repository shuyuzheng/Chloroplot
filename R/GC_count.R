gc_count <- function(genome, view.width = 10){
  n <- Biostrings::nchar(genome)
  genome <- c(genome, genome[1:view.width - 1])
  count <- Biostrings::letterFrequencyInSlidingView(genome,
                                                    letters = c("C", "G"),
                                                    view.width = view.width)
  count <- count[, 1] + count[, 2]
  total <- Biostrings::letterFrequency(genome, letters = c("C", "G"))
  total <- (total[1] + total[2])/n
  res <- data.frame(position = seq(1, n),
                    gc_count = count/view.width,
                    stringsAsFactors = FALSE)
  res <- res[seq(1, n, by = view.width), ]
  res <- list(gc_count_table = res, total_count = total)
  return(res)
}

gc_count_gene <- function(genome, gene_table){
  genes <- Biostrings::DNAStringSet(genome, start = gene_table$start,
                                    end = gene_table$end)
  gc_genes <- Biostrings::alphabetFrequency(genes)
  gene_table$gc_count <- unlist(apply(gc_genes, 1, function(x){
    return(round(sum(x[2:3])/sum(x), 2))
  }))
  return(gene_table)
}

gc_count_ir <- function(genome, ir_table){
  genes <- Biostrings::DNAStringSet(genome, start = ir_table$start + 1,
                                    end = ir_table$end)
  gc_genes <- Biostrings::alphabetFrequency(genes)
  if (nrow(ir_table) > 4) {
    gc_genes <- as.data.frame(gc_genes)
    gc_genes[1, ] <- gc_genes[5, ] <- gc_genes[1, ] + gc_genes[5, ]
    ir_table$gc_count <- unlist(apply(gc_genes, 1, function(x){
      return(round(sum(x[2:3])/sum(x), 2))
    }))
  } else {
    ir_table$gc_count <- unlist(apply(gc_genes, 1, function(x){
      return(round(sum(x[2:3])/sum(x), 2))
    }))
  }
  return(ir_table)
}
