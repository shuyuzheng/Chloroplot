
SSCrev<- function(gene_table, SSCs, SSCe){
  res <- gene_table
  for (i in 1:nrow(gene_table)) {
    tmp <- res[i, c("start", "end")]
    if (!(all(tmp < SSCs) | all(tmp > SSCe))){
      tmp <- SSCs + SSCe - tmp
      tmp <- tmp[c(2, 1)]
      res[i, c("start", "end")] <- tmp
      res[i, "strand"] <- ifelse(res[i, "strand"]=="+", "-", "+")
    }
  }
  return(res)
}
