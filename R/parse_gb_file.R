# Questions:
# Is it necessary we check column number of ORIGIN field in function "FasExtract"?

#' Extract specie's name from GB file
#'
#' @param gb a large character vector returned by function
#' \code{\link{fetch.gb}} or \code{\link{read.gb}}
#'
#' @return a character contains the specie's name.
#' @export
#'
sp.name<- function(definition, text = FALSE){
  if (text){
    sp <- gsub("(DEFINITION\\ \\ )|(\\,.*)", "", gb[2], perl = TRUE)
  } else {
    sp <- gsub("(\\,.*)", "", definition, perl = TRUE)
  }
  return(sp)
}

#' Genome extracter
#'
#' Extracting the fasta format chloroplast genome from the GeneBank File
#' @param gb a large character vector returned by function
#' \code{\link{fetch.gb}} or \code{\link{read.gb}}
#'
#' @return The genome of the sequence that is deposited at the end of the
#' GeneBank file in fasta format
#' @export
FasExtract<- function(gb){

  fasta<-gb[(grep("ORIGIN", gb)+1):length(gb)]

  # Remove empty lines
  fasta <- fasta[-which(fasta == "" | fasta == "//")]
  # while(fasta[length(fasta)]=="") {
  #   fasta<- fasta[1:length(fasta)-1]
  # }
  #
  # while(fasta[length(fasta)]=="//") {
  #   fasta<- fasta[1:length(fasta)-1]
  # }

  # Remove position indexes
  fasta <- substring(fasta, 11)
  # Extract all letters from fasta
  fasta <- gsub("[^a-zA-Z\\-]", "", fasta)
  fasta <- Reduce(c, strsplit(fasta, ""))
  fasta <- paste(fasta, collapse = "")
  fasta <- Biostrings::DNAString(fasta)
  fasta <- rdnFixer(fasta)
  return(fasta)
}

#' Randomly fill nucleotides
#'
#' This fuction change liters "u", "r", "y", "s", "w", "k", "m", "b", "d", "h",
#' "v", "n" or "-" in genome sequence (randomly) into corresponding nucleotide
#' "a", "t", "c" or "g".
#'
#' @param gb A large character vector. It contains the sequence to fix. Each
#' element of the vector represents a nucleotide.
#'
#' @return A large character vector. It countains the nucleotide sequence with
#' all letters except "a", "t", "c", "g" changed.
#' @export
#'
rdnFixer<- function(genome){
  seq <- Biostrings::toString(genome)
  seq <- unlist(strsplit(seq, ""))
  seq[which(seq=="u")]<-sample(c("t"), length(which(seq=="u")), TRUE)
  seq[which(seq=="r")]<-sample(c("a", "g"), length(which(seq=="r")), TRUE)
  seq[which(seq=="y")]<-sample(c("c", "t"), length(which(seq=="y")), TRUE)
  seq[which(seq=="s")]<-sample(c("c", "g"), length(which(seq=="s")), TRUE)
  seq[which(seq=="w")]<-sample(c("a", "t"), length(which(seq=="w")), TRUE)
  seq[which(seq=="k")]<-sample(c("g", "t"), length(which(seq=="k")), TRUE)
  seq[which(seq=="m")]<-sample(c("c", "a"), length(which(seq=="m")), TRUE)
  seq[which(seq=="b")]<-sample(c("c", "g", "t"), length(which(seq=="b")), TRUE)
  seq[which(seq=="d")]<-sample(c("a", "g", "t"), length(which(seq=="d")), TRUE)
  seq[which(seq=="h")]<-sample(c("c", "a", "t"), length(which(seq=="h")), TRUE)
  seq[which(seq=="v")]<-sample(c("c", "a", "g"), length(which(seq=="v")), TRUE)
  seq[which(seq=="n")]<-sample(c("c", "g", "t", "a"), length(which(seq=="n")), TRUE)
  seq[which(seq=="-")]<-sample(c("c", "g", "t", "a"), length(which(seq=="-")), TRUE)
  seq <- paste(seq, collapse = "")
  seq <- Biostrings::DNAString(seq)
  return(seq)
}

