#' Generate gene table
#'
#' @param gb parsed gbfile
#' @param genome a DNAstring object. It contains the genome sequence.
#'
#' @return a data frame.
#' @importFrom magrittr %>%
#' @import dplyr
#' @export

geneTable <- function(gb, genome){
  type <- lapply(gb$FEATURES, "[[", "type")

  cols <- c("start", "end", "strand",
            "type", "gene", "pseudo", "product")
  info <- NULL
  for(i in 1:length(gb$FEATURES)){
    tmp <- gb$FEATURES[[i]]
    miscol <- cols[!cols %in% colnames(tmp)]
    df <- data.frame(matrix(rep(NA, length(miscol) * nrow(tmp)),
                            nrow = nrow(tmp), ncol = length(miscol)),
                     stringsAsFactors = FALSE)
    colnames(df) <- miscol
    tmp <- cbind.data.frame(tmp, df)
    tmp <- tmp[, which(colnames(tmp) %in% cols)]
    info <- rbind.data.frame(info, tmp)
  }


  # gene
  info$gene[is.na(info$gene)] <- info$product[is.na(info$gene)]
  info$pseudo[is.na(info$pseudo)] <- FALSE

  info$gene[grepl(".*([0-9\\.]+)S.*", info$gene)] <-
    rrnFixer(info$gene[grepl(".*([0-9\\.]+)S.*", info$gene)])
  info$gene[grepl("^trn.*", info$gene, ignore.case=TRUE)] <-
    trnFixer(info$gene[grepl("^trn.*", info$gene, ignore.case=TRUE)])
  # remove duplicated tRNA and rRNA

  gene_table <- info %>%
    dplyr::filter(type %in% c("gene", "tRNA", "rRNA")) %>%
    dplyr::select(start, end, strand, gene) %>% #, pseudo) %>%
    stats::na.omit() %>%
    unique() %>%
    dplyr::mutate(chr = rep("chr1", n()))

  gene_table <- gene_table[order(gene_table[, "start"], -gene_table[, "end"]), ]
  gene_table <- gene_table[!duplicated(gene_table[, c("start", "strand", "gene")]),]
  gene_table <- gene_table[!duplicated(gene_table[, c("end", "strand", "gene")]),]
  # codon usage
  # cds <- info[which(info$type == "CDS"),]
  # cds_f <- cds[which(cds$strand == "+"),]
  # cds_seq_f <- DNAStringSet(genome, start = cds_f$start,
  #                           end = cds_f$end)
  # names(cds_seq_f) <- cds_f$gene
  # cds_cu_f <- coRdon::codonTable(cds_seq_f)
  # cds_cu_f <- coRdon::codonCounts(cds_cu_f)
  # cds_r <- cds[which(cds$strand == "+"),]
  # cds_seq_f <- DNAStringSet(genome, start = cds_r$start,
  #                           end = cds_r$end)
  return(gene_table)
}

rrnFixer <- function(rRNA){
  rRNA <- sub("[a-zA-Z]*([0-9\\.]*)S.*", "\\1", rRNA)
  rRNA <- paste("rrn", rRNA, sep = "")
}


trnFixer <- function(tRNA) {
  #tRNA <- gene_table$gene[grepl("^trn.*", gene_table$gene, ignore.case=TRUE)]
  tRNA <- sub("-", "", tRNA)
  tRNA <- sub("^tRNA", "trn", tRNA)
  aa_table <- rbind(c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu",
                      "Gln", "Gly", "His", "He", "Leu", "Lys",
                      "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
                      "Tyr", "Val"),
                    c("A", "R", "N", "D", "C", "E", "Q", "G", "H",
                      "I", "L", "K", "M", "F", "P", "S", "T", "W",
                      "Y", "V"))

  for (i in 1:ncol(aa_table)){
    tRNA <- sub(aa_table[1, i], aa_table[2, i], tRNA)
  }
  tRNA <- sub("(trnf*[A-Z]).*", "\\1", tRNA)
  #gene_table$gene[grepl("^trn.*", gene_table$gene, ignore.case=TRUE)] <- tRNA
  #return(gene_table)
  return(tRNA)
}

