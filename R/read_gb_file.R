
#' Fetch GB files online
#'
#' Fetching GB file from GenBank database by GenInfo Identifier (GI) or
#' accession number and writeing it to a file named "GI.gb" in current work
#' directory.
#' @param GI GenInfo Identifier(GI) for interested genome file.
#' @param read A logical value. If it is \code{TRUE}, the function will read
#' the fetched GB file into current R session environment.
#' @return If \code{read} is set as \code{TRUE}. The function return a large
#' character vector which contains the content of GB file for input \code{GI}
#' from GeneBank.
#' @export
fetch.gb<- function(GI, read=TRUE){
  #the GI can be also accession number
  p<- reutils::efetch(GI, "nucleotide", "gb", complexity = 1)
  #write(reutils::content(p, "text"), file = paste(GI, ".gb", sep=""))
  p <- strsplit(reutils::content(p), "\n")[[1]]
  return(p)
}

#' Read the local GB file
#'
#' Reading the GB file in the working directory producing the RGB format file
#' needed for the functions of the package
#' @param file The path to the GB file which is in GeneBank format
#' @return a large character vector returned by function
#' \code{\link{fetch.gb}} or \code{\link{read.gb}}
#' @export
read.gb <- function(file){
  return(readLines(file))
}
