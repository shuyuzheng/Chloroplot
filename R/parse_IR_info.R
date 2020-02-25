# Guestions:
# What should be the correct setting for shifter in function "phase.detector"?
# Why change lowercase to uppercase and then convert all uppercases to lower cases
# in function "genome.comp.rev"

#circular rotative function
cirtick <- function(tick, vector){
  l <- Biostrings::nchar(vector)
  if(tick > l-1 || tick < 1){
    return(vector)
  }
  else {
    return(c(vector[(tick+1):l], vector[1:tick]))
  }
}

#reverse complement function
# genome.comp.rev <- function(genome){
#   gcr<-genome[length(genome):1]
#   gcr<-gsub("a", "T", gcr)
#   gcr<-gsub("t", "A", gcr)
#   gcr<-gsub("g", "C", gcr)
#   gcr<-gsub("c", "G", gcr)
#   return(tolower(gcr))
# }

#' Auxiliary function for parallely detecting phase difference
#'
#' the function to be passed to the slaves for the parallel computing, the
#' output is with max 10 second either NA or the phase.detector
#' @param shifter a integer. It indicates the number of nucleotides which will
#' be shifted to the end of the genome sequence.
#'
#' @param genome a large character vector. It
#'
#' @return
#' @export
fun <- function(shifter, genome){
  gcr<- Biostrings::reverseComplement(genome)
  genome.tick<-cirtick(shifter, genome)
  l<-Biostrings::nchar(genome)
  stop <- round(l/4)
  for (i in 1:l){
    tmp <- sum(compareDNA(cirtick(i, genome.tick), gcr))
    if ((tmp - stop)/l > 0.07) {
      break
    }
  }
  # if (no.value) {
  #   a <- NA
  #   return(a)
  # } else {
    a <- i +shifter
    if ( a > l){
      return(a - l)
    }
    else {
      return(a)
    # }
  }
}

#' Detect the phase difference of the two inverted genomes
#'
#' @param genome A large character vector. It is the plastid genome of a species
#' as a simple vector of nucleotides.
#' @param nCPU a integer. It represent the number of CPU will be used to
#' run the function.
#'
#' @return A integer. It indicates the start position of the first IR region.
#' @export
#' @import snowfall
#'
phase.detector<- function(genome, nCPU = 1){
  #detecting the phase difference of the two inverted genomes
  if (nCPU == 1){
    shifter=84000 #this shifter is faster mode, to be sure set the value to 80000,
    gcr<- genome.comp.rev(genome)
    genome<-cirtick(shifter, genome)
    l<-length(genome)
    track<- numeric(l+1)
    track[l+1]<- round(l/4)
    for (i in 1:l){
      track[i]<- sum(cirtick(i, genome)==gcr)
      if ((track[i] - track[l+1])/l > 0.07) {
        break
      }
    }
    a<- which(track==max(track))+shifter
    if ( a > l){
      return(a - l)
    }
    else {
      return(a)
    }
  } else {
    # Parallel version
    genome <- genome
    ini.forw <- 84000
    ini.back <- ini.forw

    mm <- rep(NA, nCPU)
    snowfall::sfStop()
    snowfall::sfInit(parallel=TRUE, cpus=nCPU)
    while(sum(is.na(mm))==length(mm)){
      snowfall::sfExport("genome", "cirtick", "genome.comp.rev", "fun", "mm",
               "ini.forw", "ini.back", namespace = "circhler")
      mm <-unlist(snowfall::sfLapply(c(seq(ini.forw, ini.forw+nCPU/2*1000-1, 1000),
                            seq(ini.back, ini.back- nCPU/2*1000, -1000)[-1]),
                          fun, genome = genome))
      ini.forw <<- ini.forw + nCPU / 2 * 1000
      ini.back <<-ini.back - nCPU / 2 * 1000
    }
    sfStop()
    return(unique(mm)[which((is.na(unique(mm))==FALSE))])
  }
}

#finding the cordinate of the IR region
True.sequence.finder<- function(genome, phase.difference){
  #phase.difference<-phase.detector(genome)
  true.search<-compareDNA(cirtick(phase.difference, genome),
                          reverseComplement(genome))
  count <- rle(true.search)
}

IR1start<-function(phase.difference,True.sequence.finder){
  return(phase.difference+True.sequence.finder)
}

IR.length<- function(IR1start, True.sequence.finder, genome){
  s<-IR1start
  t<-True.sequence.finder
  r<- genome.comp.rev(genome)
  T<-cirtick(s, genome)==cirtick(t, r)
  count<-1
  while(T[count]){
    count<- count+1
  }
  if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
    count<- count+10
    while(T[count]){
      count<- count+1
    }
  }
  if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
    count<- count+10
    while(T[count]){
      count<- count+1
    }
  }
  if (sum(T[(count+10):(count+25)]==rep(TRUE, 16))==16){
    count<- count+10
    while(T[count]){
      count<- count+1
    }
  }
  count
}

IR2start<- function(True.sequence.finder, IR.length, genome){
  return(length(genome)-(True.sequence.finder+IR.length-2))
}

#' IR information
#'
#' Detecting the starts and the length of the IR regions on the chloroplast genome
#'
#' @param genome The plastid genome of a species as a simple vector of nucleotides.
#' @param nCPU number of CPU used for detecting phase difference of the two
#' inverted genomes.
#' @return a vector of four elements as the start of the first and second IR region,
#' their length and the total lenght of the genome, respectively
#' @export
#'
IRinfo<- function(genome, nCPU=4){
  phase.difference<- phase.detector(genome, nCPU)
  Trsf<- True.sequence.finder(genome, phase.difference)
  IR1s<- IR1start(phase.difference, Trsf)
  IR.l<- IR.length(IR1s, Trsf, genome)
  IR2s<- IR2start(Trsf, IR.l, genome)

  #calculation
  return(c(IR1s, IR2s, IR.l, length(genome)))#returning the start of IR one and
  #two follow by their lenght and the lenght of genome gives a vector with four
  #elemens as the start of the IRb and IRa and their length and the lenght of
  #the genome sequence
}

#' Make the table of the IR information
#'
#' This function generates a table for information about IR region from multiple
#' GenBank files.
#'
#' @param GBFiles A vector of characters. It contains the GIs or GenBank
#' accessions for the interested species' chloroplast genome files.
#'
#' @return A matrix with the IR information for each input species
#' @export
#'
IRtab<- function(GBFiles){#making the table of the IR information
  l<-length(GBFiles)
  IRtable<- matrix(0, l, 5)
  for (i in 1:l){
    IRtable[i,1]<-GBFiles[i]
    IRtable[i, 2:5]<- IRinfo(FasExtract(read.gb(GBFiles[i])))
  }
  return(IRtable)
}


LSC<- function(IRinfo){#a vecotr as a list
  c<-as.vector(IRinfo)
  IR<- c[3]
  SSC<- c[2]-(c[1]+IR)
  return(c[4]-(SSC+2*IR))
}

SSC<- function(IRinfo){
  c<-as.vector(IRinfo)
  return(c[2]-(c[1]+c[3]))
}
