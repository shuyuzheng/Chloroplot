# intronCordinates<- function(gb){#The Gene bank file as the input an
#   t<- gb[grep("  intron ", gb)]
#   m<- matrix(0, length(t), 2)
#   for (i in 1:length(t)){
#     crude<-strsplit(t[i], " ")[[1]]
#     cord<-crude[which(crude!=rep("", length(crude)))]
#     if(length(cord)!= 2){
#       stop("Check the GeneBank file: error with their intron row(s)","\n")
#     }
#     cord<-cord[2]
#     if(length(strsplit(cord, "complement")[[1]])==2){
#       cord<-strsplit(cord, "complement")[[1]][2]
#       m[i,]<- strsplit(cord, "\\..")[[1]][c(2,1)]
#     }
#     else if (length(strsplit(cord, "complement")[[1]])==1){
#       m[i,]<- strsplit(cord, "\\..")[[1]]
#     }
#     else {
#       stop("Check the Genebank file: error with cordinates of the intron row(s)")
#     }
#   }
#   m<- gsub(")", "", m)
#   m<- gsub("\\(", "", m)
#   return(m)#giving the table format of the cordinates of the genes with introns (need imporvment to return the gene name as well and their exoin cordinates too)
# }

gene.name<-function(gb, type){
  if(type=="gene"){
    t <- gb[grep("  gene  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  gene  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  gene  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  gene  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="tRNA"){
    t <- gb[grep(" tRNA  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){
        t[i] <- gb[grep("tRNA  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  tRNA  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  tRNA  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("tRNA ", gb)+1]), ""), ""), "\n")
      }
    }
  }
  else if(type=="rRNA"){
    t <- gb[grep("rRNA  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){
        t[i] <- gb[grep("rRNA  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  rRNA  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  rRNA  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else if(type=="mRNA"){
    t <- gb[grep("mRNA  ", gb)+1]
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      if(length(grep("=", crude))==0){
        t[i] <- gb[grep("mRNA  ", gb)+2][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  mRNA  ", gb)+3][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      if(length(grep("=", crude))==0){#if it cannot find the gene name in the following line go to the next
        t[i] <- gb[grep("  mRNA  ", gb)+4][i]
        crude<- strsplit(t[i], " ")[[1]]
      }
      crude.name<-crude[which(crude!=rep("", length(crude)))]
      t[i]<-gsub("\"", "", strsplit(crude.name, "=")[[1]][2])
    }
    na<-which(is.na(t))
    if(length(na) > 0){
      for (i in 1:length(na)){
        warning(paste(paste(paste("Gene No.", na[i], " "), paste("is not properly named and is deleted from the list."), ""), paste("Check the gb file on line", which(gb==gb[grep("gene ", gb)+1][na]), ""), ""))
      }
    }
  }
  else {
    stop("The type should be defined as either gene or tRNA")
  }
  t<-t[!is.na(t)]
  return(t)
  #intermediate gene name function to substract the gene names of either gene or tRNA, mRNA or rRNA from their second line information(or third), the input is the gb file.
}

gb.gene.cor<- function(gb){
  gene<- gb[grep("  gene  ", gb)]
  # trna<- gb[grep("  tRNA ", gb)]
  rrna<- gb[grep("  rRNA ", gb)]
  # t<-c(gene, trna, rrna)
  t <- c(gene, rrna)
  m<- matrix(0, length(t), 2)
  for (i in 1:length(t)){
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      s<-t[i]
      t[i]<-paste(t[i], gsub(" ", "", gb[which(t[i]==gb)+1]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]==","){#check again if the last characted if its ',' and add the next one from gb file to that
      #previous code--> if(strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])]=="," && strsplit(t[i], "")[[1]][length(strsplit(t[i],"")[[1]])-1]==")")
      t[i]<-paste(t[i], gsub(" ", "", gb[which(s==gb)+2]), "")
    }
    t[i]<- gsub(", ", ",", t[i])
    t[i]<- gsub(") ", ")", t[i])
    m[i,2]<-strsplit(t[i], " ")[[1]][length(strsplit(t[i], " ")[[1]])]
  }
  # names<-  c(gene.name(gb, "gene"), gene.name(gb, "tRNA"), gene.name(gb, "rRNA") )
  names<-  c(gene.name(gb, "gene"), gene.name(gb, "rRNA") )
  if(length(t)!=length(names)){
    stop("Error: Check the gene names of the gb file")
  }
  m[,1]<-names
  m #gives the crude table extracted from gb file for genes
}

gene.cordinates<- function(gb){
  gb.gene.cor.out<-gb.gene.cor(gb)
  l<- length(gb.gene.cor.out[,1])
  m<- matrix(0, l, 5)
  m[,1]<- gb.gene.cor.out[,1]
  for (i in 1:l){
    if(length(strsplit(gb.gene.cor.out[i,2], ",")[[1]])==2){
      m[i,2]<-sub("join\\(", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][1])
      m[i,2]<-sub("order\\(", "", m[i,2])
      m[i,4]<-sub(")", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][2])
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
      cord<-m[i,4]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
    else {
      m[i,2]<-gb.gene.cor.out[i,2]
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
  }
  m<- gsub(")", "", m)
  m<- gsub("\\(", "", m)
  m<- gsub("<", "", m)
  m<- gsub(">", "", m)
  return(m)#gives a polished table of genes and their cordinates, as the gene names their start and end of first part and the second parts as the first to fifth columns. The orders are reflected with the swapping of values
}

gene.cordinates<- function(gb){
  gb.gene.cor.out<-gb.gene.cor(gb)
  l<- length(gb.gene.cor.out[,1])
  m<- matrix(0, l, 5)
  m[,1]<- gb.gene.cor.out[,1]
  for (i in 1:l){
    if(length(strsplit(gb.gene.cor.out[i,2], ",")[[1]])==2){
      m[i,2]<-sub("join\\(", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][1])
      m[i,2]<-sub("order\\(", "", m[i,2])
      m[i,4]<-sub(")", "", strsplit(gb.gene.cor.out[i,2], ",")[[1]][2])
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
      cord<-m[i,4]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(4,5)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
    else {
      m[i,2]<-gb.gene.cor.out[i,2]
      cord<-m[i,2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,c(2,3)]<- strsplit(cord, "\\..")[[1]]
      }
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
  }
  m<- gsub(")", "", m)
  m<- gsub("\\(", "", m)
  m<- gsub("<", "", m)
  m<- gsub(">", "", m)
  return(m)#gives a polished table of genes and their cordinates, as the gene names their start and end of first part and the second parts as the first to fifth columns. The orders are reflected with the swapping of values
}
#
# trnNfixer<- function(Genelist){
#   for (i in 1:length(Genelist[,1])){
#     if (Genelist[i,1] %in% c("trnN-GUU", "tRNA-Asn")){
#       Genelist[i,1]<- "trnN"
#     }
#   }
#   return(Genelist)
# }


chr.count<- function(word){
  if (length(word)==1){
    return(length(strsplit(word, "")[[1]]))
  }
  else {
    t<- numeric(length(word))
    for (i in 1:length(word)){
      t[i]<- length(strsplit(word[i], "")[[1]])
    }
    return(t)
  }
}
