# SSCrev<- function(genecord, SSCs, SSCe){
#   g<- genecord
#   a<- SSCs
#   b<- SSCe
#   g[, 2]<- as.numeric(g[,2])
#   g[, 3]<- as.numeric(g[,3])
#   g[, 4]<- as.numeric(g[,4])
#   g[, 5]<- as.numeric(g[,5])
#   for (i in 1:length(g[,1])){
#     #if the two values of the second and third column are in the SSC, then reverse them
#     if(max(g[i, 2], g[i, 3]) <= b  && min(g[i, 2], g[i, 3]) >= a){
#       g[i, 2]<- paste(b-(as.numeric(g[i, 2])-a))
#       g[i, 3]<- b-(as.numeric(g[i, 3])-a)
#     }
#     #If the gene is starting from SSC then it should also be fixed
#     #on the positive strand
#     if(g[i, 2] < a & g[i, 3] > a){
#       g[i, 2]<- b+(a-as.numeric(g[i,2]))
#       g[i, 3]<- b-(as.numeric(g[i,3])-a)
#     }
#     #reverseve of that
#     else if(g[i, 2] > b & g[i, 3] < b){
#       g[i, 2]<- a-(as.numeric(g[i,2])-b)
#       g[i, 3]<- a+(b-as.numeric(g[i,3]))
#     }
#     #on the negative strand
#     else if(g[i, 2] < b & g[i, 3] > b){
#       g[i, 2]<- a+(b-as.numeric(g[i,2]))
#       g[i, 3]<- a-(as.numeric(g[i,3])-b)
#     }
#     #reverseve of that
#     else if(g[i, 2] > a & g[i, 3] < a){
#       g[i, 2]<- b-(as.numeric(g[i,2])-a)
#       g[i, 3]<- b+(a-as.numeric(g[i,3]))
#     }
#     ###for the 4th and 5th columns
#     #if the two values of the second and third column are in the SSC, then reverse them
#     if(max(g[i, 4], g[i, 5]) <= b  && min(g[i, 4], g[i, 5]) >= a){
#       g[i, 4]<- paste(b-(as.numeric(g[i, 4])-a))
#       g[i, 5]<- b-(as.numeric(g[i, 5])-a)
#     }
#     #If the gene is starting from SSC then it should also be fixed
#     #on the positive strand
#     if(g[i, 4] < a & g[i, 5] > a){
#       g[i, 4]<- b+(a-as.numeric(g[i,4]))
#       g[i, 5]<- b-(as.numeric(g[i,5])-a)
#     }
#     #reverseve of that
#     else if(g[i, 4] > b & g[i, 5] < b){
#       g[i, 4]<- a-(as.numeric(g[i,4])-b)
#       g[i, 5]<- a+(b-as.numeric(g[i,5]))
#     }
#     #on the negative strand
#     else if(g[i, 4] < b & g[i, 5] > b){
#       g[i, 4]<- a+(b-as.numeric(g[i,4]))
#       g[i, 4]<- a-(as.numeric(g[i,5])-b)
#     }
#     #reverseve of that
#     else if(g[i, 4] > a & g[i, 5] < a){
#       g[i, 4]<- b-(as.numeric(g[i,4])-a)
#       g[i, 5]<- b+(a-as.numeric(g[i,5]))
#     }
#   }
#   g[, 2]<- as.character(g[,2])
#   g[, 3]<- as.character(g[,3])
#   g[, 4]<- as.character(g[,4])
#   g[, 5]<- as.character(g[,5])
#   return(g)
# }

SSCrev<- function(gene_table, SSCs, SSCe){
  res <- gene_table
  for (i in 1:nrow(gene_table)) {
    tmp <- res[i, c("start", "end")]
    if (!(all(tmp < SSCs) | all(tmp > SSCe))){
      tmp <- SSCs + SSCe - tmp
      res[i, c("start", "end")] <- tmp
    }
  }
  return(res)
}
