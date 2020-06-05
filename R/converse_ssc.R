
convert_gene<- function(gene_table, region_start, region_end, l){
  res <- gene_table
  append_gene <- NULL
  gene_out_region_end <- gene_table$gene[which(gene_table$start >= region_end)[1]]
  if (is.na(gene_out_region_end)){
    gene_out_region_end <- gene_table$gene[1]
  }

  gene_out_region_start <- gene_table$gene[which(gene_table$end <= region_start)]
  gene_out_region_start <- tail(gene_out_region_start, 1)
  if (is.na(gene_out_region_start)){
    gene_out_region_start <- tail(gene_table$gene)
  }

  gene_cut_by_start <- gene_table$gene[gene_table$start < region_start & gene_table$end > region_start]
  if (length(gene_cut_by_start) == 0){
    gene_cut_by_start <- "NO_CUT_GENE"
  }
  gene_cut_by_end <- gene_table$gene[gene_table$start < region_end & gene_table$end > region_end]
  if (length(gene_cut_by_end) == 0){
    gene_cut_by_end <- "NO_CUT_GENE"
  }
  for (i in 1:nrow(gene_table)) {
    tmp <- res[i, c("start", "end")]
    if (region_start < region_end){
      if (all(tmp >= region_start) & all(tmp <= region_end)){
        tmp <- region_start + region_end - tmp
        tmp <- tmp[c(2, 1)]
        res[i, c("start", "end")] <- tmp
        res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
      } else if (tmp[1] < region_start & tmp[2] > region_start){
        if (res$gene[i] != gene_out_region_end){
          if (res$gene[i] != gene_cut_by_end){
            tmp2 <- res[i, ]
            tmp2$end <- region_start
            tmp2$pseudo <- TRUE
            append_gene <- rbind.data.frame(append_gene, tmp2)

            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            res[i, c("start", "end")] <- tmp
            res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            res[i, "end"] <- region_end
            res[i, "pseudo"] <- TRUE
          } else {
            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            res[i, c("start", "end")] <- tmp
            res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
          }
        }
      } else if (tmp[1] < region_end & tmp[2] > region_end){
        if (res$gene[i] != gene_out_region_end){
          if (res$gene[i] != gene_cut_by_start) {
            tmp2 <- res[i, ]
            tmp2$start<- region_end
            tmp2$pseudo <- TRUE
            append_gene <- rbind.data.frame(append_gene, tmp2)

            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            res[i, c("start", "end")] <- tmp
            res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            res[i, "start"] <- region_start
            res[i, "pseudo"] <- TRUE
          } else {
            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            res[i, c("start", "end")] <- tmp
            res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
          }
        }
      }
    } else {
      if (all(tmp <= region_end) | all(tmp >= region_start)){
        tmp <- region_start + region_end - tmp
        tmp <- tmp[c(2, 1)]
        tmp[tmp > l] <- tmp[tmp > l] - l
        tmp[tmp < 0] <- l - tmp[tmp < 0]
        if (tmp[1]>tmp[2]){
          res[i, c("start", "end")] <- c(tmp[1], l + tmp[2])
          res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
        } else{
          res[i, c("start", "end")] <- tmp
          res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
        }
      } else if (tmp[1] < region_start & tmp[2] > region_start){
        if (res$gene[i] != gene_out_region_end ){
          if (res$gene[i] != gene_cut_by_end){
            tmp2 <- res[i, ]
            tmp2$end <- region_start
            tmp2$pseudo <- TRUE
            append_gene <- rbind.data.frame(append_gene, tmp2)

            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            tmp[tmp > l] <- tmp[tmp > l] - l
            tmp[tmp < 0] <- l - tmp[tmp < 0]
            if (tmp[1]>tmp[2]){
              res[i, c("start", "end")] <- c(tmp[1], l + tmp[2])
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            } else{
              res[i, c("start", "end")] <- tmp
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            }
            res[i, "end"] <- region_end
            res[i, "pseudo"] <- TRUE
          } else {
            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            tmp[tmp > l] <- tmp[tmp > l] - l
            tmp[tmp < 0] <- l - tmp[tmp < 0]
            if (tmp[1]>tmp[2]){
              res[i, c("start", "end")] <- c(tmp[1], l + tmp[2])
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            } else{
              res[i, c("start", "end")] <- tmp
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            }
          }

        }
      } else if (tmp[1] < region_end & tmp[2] > region_end){
        if (res$gene[i] != gene_out_region_end){
          if (res$gene[i] != gene_cut_by_start){
            tmp2 <- res[i, ]
            tmp2$start<- region_end
            tmp2$pseudo <- TRUE
            append_gene <- rbind.data.frame(append_gene, tmp2)

            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            tmp[tmp > l] <- tmp[tmp > l] - l
            tmp[tmp < 0] <- l - tmp[tmp < 0]
            if (tmp[1]>tmp[2]){
              res[i, c("start", "end")] <- c(tmp[1], l + tmp[2])
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            } else{
              res[i, c("start", "end")] <- tmp
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            }
            res[i, "start"] <- region_start
            res[i, "pseudo"] <- TRUE
          } else {
            tmp <- region_start + region_end - tmp
            tmp <- tmp[c(2, 1)]
            tmp[tmp > l] <- tmp[tmp > l] - l
            tmp[tmp < 0] <- l - tmp[tmp < 0]
            if (tmp[1]>tmp[2]){
              res[i, c("start", "end")] <- c(tmp[1], l + tmp[2])
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            } else{
              res[i, c("start", "end")] <- tmp
              res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
            }
          }
        }
      }
    }
  }

  res <- rbind.data.frame(res, append_gene)
  return(res)
}

convert_small_region <- function(gene_table, region_start, region_end, l){
  res <- gene_table

  for (i in 1:nrow(gene_table)) {
    tmp <- res[i, c("start", "end")]
    if (region_start < region_end){
      if (!(all(tmp < region_start) | all(tmp > region_end))){
        tmp <- region_start + region_end - tmp
        tmp <- tmp[c(2, 1)]
        res[i, c("start", "end")] <- tmp
      }
    } else {
      if (!(all(tmp > region_end) & all(tmp < region_start))){
        tmp <- region_start + region_end - tmp
      tmp <- tmp[c(2, 1)]
      tmp[tmp > l] <- tmp[tmp > l] - l
      tmp[tmp < 0] <- l - tmp[tmp < 0]
        if (tmp[1]>tmp[2]){
          res[i, c("start", "end")] <- c(tmp[1], l + tmp[2])
        } else{
        res[i, c("start", "end")] <- tmp
        }
      }
    }
  }
  return(res)
}

convert_point <- function(position, region_start, region_end){
    if (region_start < region_end){
      position[which(position >= region_start & position <= region_end)] <- region_start + region_end -
        position[which(position >= region_start & position <= region_end)]
    } else {
      position[which(position >= region_start | position <= region_end)] <- region_start + region_end -
        position[which(position >= region_start | position <= region_end)]
    }
  return(position)
}

convert_region <- function(ir_table, l, region = "SSC", genome, gene_table,
                           customize.ring2, customize.ring3, customize.ring1,
                           indel_table = NULL) {
  if (sum(ir_table$name == region) == 2){
    region_start <- ir_table$start[ir_table$name == region][2]
    region_end <- ir_table$end[ir_table$name == region][1]
    genome_ssc_convert <- c(genome[l:region_start], genome[region_end:region_start], genome[region_end:1])
  } else {
    region_start <- ir_table$start[ir_table$name == region]
    region_end <- ir_table$end[ir_table$name == region]
    genome_ssc_convert <- c(genome[1:region_start], genome[region_end:region_start],
                            genome[region_end:l])
  }
  gene_table <- convert_gene(gene_table, region_start = region_start, region_end = region_end, l)

  gc.window <- 100
  if (l > 500000){
    gc.window <- 200
  } else if (l < 100000){
    gc.window <- 50
  }
  gc_count_list <- gc_count(genome_ssc_convert, view.width = gc.window)
  gc_count <- gc_count_list[[1]]
  gc_total <- gc_count_list[[2]]
  gc_count$chr <- rep("chr1", nrow(gc_count))
  gc_count <- select(gc_count, chr, position, gc_count)

  if (!is.null(customize.ring3)){
    customize.ring3 <- convert_small_region(customize.ring3, region_start = region_start, region_end = region_end, l)
  }
  if (!is.null(customize.ring1)){
    customize.ring1$position <- convert_point(customize.ring1$position, region_start = region_start, region_end = region_end)
  }
  if (!is.null(customize.ring2)){
    customize.ring2$position <- convert_point(customize.ring2$position, region_start = region_start, region_end = region_end)
  }
  if (!is.null(indel_table)){
    indel_table$position <- convert_point(indel_table$position, region_start = region_start, region_end = region_end)
  }
  res <- list(gene_table = gene_table, gc_count = gc_count, gc_total = gc_total,
              customize.ring1 = customize.ring1,
              customize.ring2 = customize.ring2,
              customize.ring3 = customize.ring3, indel_table = indel_table)
  return(res)
}
