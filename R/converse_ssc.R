
SSCrev<- function(gene_table, SSCs, SSCe, l){
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
