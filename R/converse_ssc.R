convert_gene<- function(gene_table, SSCs, SSCe, l){
  res <- gene_table
  for (i in 1:nrow(gene_table)) {
    tmp <- res[i, c("start", "end")]
    if (SSCs < SSCe){
      if (!(all(tmp < SSCs) | all(tmp > SSCe))){
        tmp <- SSCs + SSCe - tmp
        tmp <- tmp[c(2, 1)]
        res[i, c("start", "end")] <- tmp
        res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
      }
    } else {
      if (!(all(tmp > SSCe) & all(tmp < SSCs))){
        tmp <- SSCs + SSCe - tmp
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
  return(res)
}

convert_point <- function(position, SSCs, SSCe){
    if (SSCs < SSCe){
      position[which(position >= SSCs & position <= SSCe)] <- SSCs + SSCe -
        position[which(position >= SSCs & position <= SSCe)]
    } else {
      position[which(position >= SSCs | position <= SSCe)] <- SSCs + SSCe -
        position[which(position >= SSCs | position <= SSCe)]
    }
  return(position)
}

convert_region <- function(ir_table, gene_table, l, region = "SSC") {
  if (nrow(ir_table) == 1){
    warning("Didn't get IR region from thid species. It's impossible to convert ",
            region, " region.")
    gene_table <- gene_table
  } else {
    if (sum(ir_table$name == region) == 2){
      SSCs <- ir_table$start[ir_table$name == region][2]
      SSCe <- ir_table$end[ir_table$name == region][1]
      genome_ssc_convert <- c(genome[L:SSCs], genome[SSCe:SSCs], genome[SSCe:1])
    } else {
      SSCs <- plot.tables$ir_table$start[plot.tables$ir_table$name == region]
      SSCe <- plot.tables$ir_table$end[plot.tables$ir_table$name == region]
      genome_ssc_convert <- c(plot.tables$genome[1:SSCs], plot.tables$genome[SSCe:SSCs],
                              plot.tables$genome[SSCe:L])
    }
    gene_table <- SSCrev(plot.tables$gene_table, SSCs = SSCs, SSCe = SSCe, L)
    gc.window <- 100
    if (L > 500000){
      gc.window <- 200
    } else if (L < 100000){
      gc.window <- 50
    }
    gc_count_list <- gc_count(genome_ssc_convert, view.width = gc.window)
    gc_count <- gc_count_list[[1]]
    gc_total <- gc_count_list[[2]]
    gc_count$chr <- rep("chr1", nrow(gc_count))
    gc_count <- select(gc_count, chr, position, gc_count)

    if (!is.null(customize.ring3)){
      customize.ring3 <- SSCrev(customize.ring3, SSCs = SSCs, SSCe = SSCe, L)
    }
    if (!is.null(customize.ring1)){
      customize.ring1$position <- SSCrev_point(customize.ring1$position)
    }
    if (!is.null(customize.ring2)){
      customize.ring2$position <- SSCrev_point(customize.ring2$position)
    }
  }
}
