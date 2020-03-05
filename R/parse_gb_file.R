# Questions:
# Is it necessary we check column number of ORIGIN field in function "FasExtract"?

#' Extract specie's name from GB file
#'
#' @param definition a charactor. It contain the definition field from the
#' GenBank file.
#'
#' @return a character contains the specie's name.
#' @export
#'
sp.name<- function(definition){
  # if (text){
    # sp <- gsub("(DEFINITION\\ \\ )", "", definition[2], perl = TRUE)
    # sp <- sub("(\\w+\\s+\\w+).*", "\\1", sp, perl = TRUE)
  # } else {
    sp <- sub("(\\w+\\s+\\w+).*", "\\1", definition, perl = TRUE)
  # }
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
#' @param genome A DNAString object. It is a class defined in \code{Biostrings}
#' package. It contains the whole genome sequence of a specie.
#'
#' @return A large character vector. It countains the nucleotide sequence with
#' all letters except "a", "t", "c", "g" changed.
#' @export
#'
rdnFixer<- function(genome){
  seq <- Biostrings::toString(genome)
  seq <- unlist(strsplit(seq, ""))
  seq[which(seq=="U")]<-sample(c("T"), length(which(seq=="U")), TRUE)
  seq[which(seq=="R")]<-sample(c("A", "G"), length(which(seq=="R")), TRUE)
  seq[which(seq=="Y")]<-sample(c("C", "T"), length(which(seq=="Y")), TRUE)
  seq[which(seq=="S")]<-sample(c("C", "G"), length(which(seq=="S")), TRUE)
  seq[which(seq=="W")]<-sample(c("A", "T"), length(which(seq=="W")), TRUE)
  seq[which(seq=="K")]<-sample(c("G", "T"), length(which(seq=="K")), TRUE)
  seq[which(seq=="M")]<-sample(c("C", "A"), length(which(seq=="M")), TRUE)
  seq[which(seq=="B")]<-sample(c("C", "G", "T"), length(which(seq=="B")), TRUE)
  seq[which(seq=="D")]<-sample(c("A", "G", "T"), length(which(seq=="D")), TRUE)
  seq[which(seq=="H")]<-sample(c("C", "A", "T"), length(which(seq=="H")), TRUE)
  seq[which(seq=="V")]<-sample(c("C", "A", "G"), length(which(seq=="V")), TRUE)
  seq[which(seq=="N")]<-sample(c("C", "G", "T", "A"), length(which(seq=="N")), TRUE)
  seq[which(seq=="-")]<-sample(c("C", "G", "T", "A"), length(which(seq=="-")), TRUE)
  seq <- paste(seq, collapse = "")
  seq <- Biostrings::DNAString(seq)
  return(seq)
}

