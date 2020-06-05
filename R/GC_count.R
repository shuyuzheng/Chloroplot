gc_count <- function(genome, view.width = 10){
  n <- Biostrings::nchar(genome)
  genome <- c(genome, genome[1:view.width])
  count <- as.vector(Biostrings::letterFrequencyInSlidingView(genome,
                                                    letters = "GC",
                                                    view.width = view.width))
  total <- as.vector(Biostrings::letterFrequency(genome, letters = "GC"))/n
  res <- data.frame(position = seq(1, n+1),
                    gc_count = count/view.width,
                    stringsAsFactors = FALSE)
  res <- res[c(seq(1, n, by = view.width), n), ]
  res <- list(gc_count_table = res, total_count = total)
  return(res)
}

gc_count_gene <- function(genome, gene_table){
  genes <- Biostrings::DNAStringSet(genome, start = gene_table$start,
                                    end = gene_table$end)
  gc_genes <- as.vector(Biostrings::letterFrequency(genes, letters = "GC"))
  len_genes <- Biostrings::nchar(genes)
  gene_table$gc_count <- gc_genes/len_genes
  return(gene_table)
}

gc_count_ir <- function(genome, ir_table){
  genes <- Biostrings::DNAStringSet(genome, start = ir_table$start + 1,
                                    end = ir_table$end)
  gc_genes <- as.vector(Biostrings::letterFrequency(genes, letters = "GC"))
  len_genes <- Biostrings::nchar(genes)
  if (nrow(ir_table) == 5) {
    gc_genes[1] <- gc_genes[5] <- gc_genes[1] + gc_genes[5]
    len_genes[1] <- len_genes[5] <- len_genes[1] + len_genes[5]
    ir_table$gc_count <- gc_genes/len_genes
  } else {
    ir_table$gc_count <- gc_genes/len_genes
  }
  return(ir_table)
}
