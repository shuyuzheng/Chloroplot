


#colouring

#photosystem I {psa}
#photosystem II{psb}
#cytochrome b/f complex {pet}
#ATP synthesis {atp}
#NADH dehydrogenase {ndh}
#RubisCO larg subunit {rbc}
#RNA polymerase {rpo}
#small ribosomal protein {rps}
#large ribosomal protein {rpl}
#clpP, matK, infA {clp, mat, inf}
#hypothetical reading frame {ycf}
#transfer RNA {trn}
#ribosomal RNA {rrn}


gbFile <- "Solanum lycopersicum copy.gb"


#
  
  #calling libraries
  
  
  library(plotrix)
  library(snow)
  library(snowfall)
  
  Arc.plotter <- function (theta0, theta1, radius.in, radius.out, center.x = 0, 
                          center.y = 0, edges = 10, col = "black", border = NA, 
                          ...) {
    if (length(edges) == 1) 
      edges <- rep(edges, length = length(theta0))
    col <- rep(col, length = length(theta0))
    ok <- !(is.na(theta0) | is.na(theta1))
    for (i in seq(along = theta0[ok])) {
      theta.seq <- seq(theta0[ok][i], theta1[ok][i], length = edges[ok][i])
      x <- c(radius.in * cos(theta.seq), radius.out * cos(rev(theta.seq))) + 
        center.x
      y <- c(radius.in * sin(theta.seq), radius.out * sin(rev(theta.seq))) + 
        center.y
      polygon(x, y, col = col[ok][i], border = border, ...)
    }
  }
  
  FasExtract<- function(gb){
    #' Genome extracter
    #' 
    #' Extracting the fasta format chloroplast genome from the GeneBank File
    #' @param file Name of the GeneBank file of a sequence on the working reposotory
    #' @return The genome of the sequence that is deposited at the end of the GeneBank file in fasta format
    #' @export
    fasta<-gb[(grep("ORIGIN", gb)+1):length(gb)]
    while(fasta[length(fasta)]=="") {
      fasta<- fasta[1:length(fasta)-1]
    }
    while(fasta[length(fasta)]=="//") {
      fasta<- fasta[1:length(fasta)-1]
    }
    fas<-""
    for (i in 1:length(fasta)){
      sort.let <- sort(unique(c(grep("c", strsplit(fasta[i], " ")[[1]]),
                               grep("a", strsplit(fasta[i], " ")[[1]]), 
                               grep("t", strsplit(fasta[i], " ")[[1]]), 
                               grep("g", strsplit(fasta[i], " ")[[1]]))))
      try(if(length(sort.let)==!6) stop("Check the gb file; the columns of ORIGIN should be 6"))
      fasta[i] <- paste(strsplit(fasta[i], " ")[[1]][sort.let[1]],
                       strsplit(fasta[i], " ")[[1]][sort.let[2]],
                       strsplit(fasta[i], " ")[[1]][sort.let[3]],
                       strsplit(fasta[i], " ")[[1]][sort.let[4]],
                       strsplit(fasta[i], " ")[[1]][sort.let[5]],
                       strsplit(fasta[i], " ")[[1]][sort.let[6]], sep="")
      fas <- paste(fas, fasta[i], sep="")
    }
    #fasta[length(fasta)]<- gsub("NA", "", fasta[length(fasta)])
    fas <- gsub("NA", "", fas)
    strsplit(fas, "")[[1]]
  } 
  
  read.gb <- function(file){
    #' R version of the GB file
    #'
    #' Reading the GB file in the working directory producing the RGB format file needed for the functions of the package
    #'@param file The name of the file in the directory in GeneBank format
    #'@return An object of class Rgb
    #'@export
    return(readLines(file))
  }
  
  fetch.gb<- function(GI, read=TRUE){
    #the GI can be also accession number
    p<- efetch(GI, "nucleotide", "gb")
    write(content(p, "text"), file = paste(GI, ".gb", sep=""))
    if (read){
      read.gb(paste(GI, ".gb", sep=""))
    }
  }
  
  rdnFixer<- function(gb){
    seq<- FasExtract(gb)
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
    return(seq)
  } 
  
  sp.name<- function(gb){
    paste(strsplit(gb[2], " ")[[1]][3], strsplit(gb[2], " ")[[1]][4], sep=" ")
  }
  
  IRinfo<- function(genome, parallel=TRUE){
    
    #' IR information
    #' 
    #' Detecting the starts and the length of the IR regions on the chloroplast genome
    #' 
    #' @param genome The plastid genome of a species as a simple vector of nucleotides 
    #' @return a vector of four elements as the start of the first and second IR region, 
    #' their length and the total lenght of the genome, respectively
    #' @export 
    
    ###Preliminary functions 
    
    cirtick <<- function(tick, vector){#circular rotative function
      if(tick > length(vector)-1 || tick < 1){
        return(vector)
      }
      else {
        return(c(vector[(tick+1):length(vector)], vector[1:tick]))
      }
    }
    
    genome.comp.rev<<- function(genome){#reverse complement function
      gcr<-genome[length(genome):1]
      gcr<-gsub("a", "T", gcr)
      gcr<-gsub("t", "A", gcr)
      gcr<-gsub("g", "C", gcr)
      gcr<-gsub("c", "G", gcr)
      return(tolower(gcr))
    }
    
    #Checked
    phase.detector<- function(genome){
      #detecting the phase difference of the two inverted genomes
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
    }
    #the parallel version of the phase.detector function. Set with default 4
    p.d <<- function(genome, nCPU=4){#tested in the GlobEnv, it might fail when embedded in the bigger function
      genome<<-genome
      fun<<- function(shifter){#the function to be passed to the slaves for the parallel computing, the out put is with max 10 second either NA or the phase.detector
        gcr<- genome.comp.rev(genome)
        cir.genome<<-cirtick(shifter, genome)
        l<-length(genome)
        track<- numeric(l+1)
        track[l+1]<- round(l/4)
        s.time<- Sys.time()
        no.value<- FALSE
        for (i in 1:l){
          track[i]<- sum(cirtick(i, cir.genome)==gcr)
          i.time<- Sys.time()
          if ((track[i] - track[l+1])/l > 0.07) {
            break
          }
          if (i.time - s.time > 11){
            no.value<- TRUE
            break
          }
        }
        if (no.value) {
          a<- NA
          return(a)
        }
        else {
          a<- which(track==max(track))+shifter
          if ( a > l){
            return(a - l)
          }
          else {
            return(a)
          } 
        }
      }
      ini.forw<<- 84000
      ini.back<<- ini.forw
      mm<<- rep(NA, nCPU)
      sfStop()
      sfInit(parallel=TRUE, cpus=nCPU)
      while(sum(is.na(mm))==length(mm)){
        sfExport("genome", "cirtick", "genome.comp.rev", "fun", "mm", "ini.forw", "ini.back")
        mm<-unlist(sfLapply(c(seq(ini.forw, ini.forw+nCPU/2*1000-1, 1000), seq(ini.back, ini.back- nCPU/2*1000, -1000)[-1]),  fun))
        ini.forw<<- ini.forw+nCPU/2*1000
        ini.back<<-ini.back-nCPU/2*1000
      }
      sfStop()
      return(unique(mm)[which((is.na(unique(mm))==FALSE))])
    }
    
    #Checked
    True.sequence.finder<- function(genome, phase.difference){#finding the cordinate of the IR region
      #phase.difference<-phase.detector(genome)
      true.search<-cirtick(phase.difference, genome)==genome.comp.rev(genome)
      true.arm<- round(length(genome)/100)
      for (i in 1:length(genome)){
        if (sum(true.search[i:(true.arm+i-1)])==true.arm) {
          return(i); break
        }
      }
    }
    
    #Checked
    IR1start<-function(phase.difference,True.sequence.finder){
      return(phase.difference+True.sequence.finder)
    }
    
    #Checked
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
    
    #declaration
    if (parallel){
      phase.difference<- p.d(genome)
    }
    else{
      phase.difference<- phase.detector(genome)
    }
    Trsf<- True.sequence.finder(genome, phase.difference)
    IR1s<- IR1start(phase.difference, Trsf)
    IR.l<- IR.length(IR1s, Trsf, genome)
    IR2s<- IR2start(Trsf, IR.l, genome)
    
    #calculation
    return(c(IR1s, IR2s, IR.l, length(genome)))#returning the start of IR one and two follow by their lenght and the lenght of genome
    #gives a vector with four elemens as the start of the IRb and IRa and their length and the lenght of the genome sequence
  }
  
  IRtab<- function(GBFiles){#making the table of the IR information
    l<-length(GBFiles)
    IRtable<- matrix(0, l, 5)
    for (i in 1:l){
      IRtable[i,1]<-GBFiles[i]
      IRtable[i, 2:5]<- IRinfo(FasExtract(read.gb(GBFiles[i])))
    }
    return(IRtable)
  }
  
  intronCordinates<- function(gb){#The Gene bank file as the input an
    t<- gb[grep("  intron ", gb)]
    m<- matrix(0, length(t), 2)
    for (i in 1:length(t)){
      crude<-strsplit(t[i], " ")[[1]]
      cord<-crude[which(crude!=rep("", length(crude)))]
      if(length(cord)!= 2){
        stop("Check the GeneBank file: error with their intron row(s)","\n")
      }
      cord<-cord[2]
      if(length(strsplit(cord, "complement")[[1]])==2){
        cord<-strsplit(cord, "complement")[[1]][2]
        m[i,]<- strsplit(cord, "\\..")[[1]][c(2,1)]
      }
      else if (length(strsplit(cord, "complement")[[1]])==1){
        m[i,]<- strsplit(cord, "\\..")[[1]]
      } 
      else {
        stop("Check the Genebank file: error with cordinates of the intron row(s)")
      }
    }
    m<- gsub(")", "", m)
    m<- gsub("\\(", "", m)
    return(m)#giving the table format of the cordinates of the genes with introns (need imporvment to return the gene name as well and their exoin cordinates too)
  }  
  
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
  
  JunctRadiusFinder<- function(gene.cordinates, IRinfo, J.pos, radius, silence=TRUE){
    #function to find the genes near the junction cordinate based on the output of the gene.cordinate function
    if(J.pos==1){
      J<- IRinfo[1]
    }
    else if(J.pos==2){
      J<- IRinfo[1]+IRinfo[3]
    }
    else if(J.pos==3){
      J<- IRinfo[2]
    }
    else if(J.pos==4){
      J<- IRinfo[2]+IRinfo[3]
    }
    else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
    g<-gene.cordinates
    l<- ncol(g)
    r<- radius
    gs<- IRinfo[4]
    t1<- subset(g, J-r<as.numeric(g[,2]) & as.numeric(g[,2])<J+r)
    t2<- subset(g, J-r<as.numeric(g[,3]) & as.numeric(g[,3])<J+r)
    t3<- subset(g, J-r<as.numeric(g[,4]) & as.numeric(g[,4])<J+r)
    t4<- subset(g, J-r<as.numeric(g[,5]) & as.numeric(g[,5])<J+r)
    for (i in 1:nrow(g)){
      if(sum(g[i, ]=="0")==2){
        g[i,4:5]<-c(NA, NA)
      }
    }
    t5<- subset(g, J-r<(as.numeric(g[,2])+gs) & (as.numeric(g[,2])+gs)<J+r)
    t6<- subset(g, J-r<(as.numeric(g[,3])+gs) & (as.numeric(g[,3])+gs)<J+r)
    t7<- subset(g, J-r<(as.numeric(g[,4])+gs) & (as.numeric(g[,4])+gs)<J+r)
    t8<- subset(g, J-r<(as.numeric(g[,5])+gs) & (as.numeric(g[,5])+gs)<J+r)
    t<-unique(do.call("rbind", list(t1, t2, t3, t4, t5, t6, t7, t8)))
    if((length(t)+silence)==0){
      warning("There is no gene on the specified junction with the given radius, increasing the radius might help")
    }
    return(t)
  }
  
  JunctNoFinder<- function(gene.cordinates, IRinfo, J.pos ,Number){#the same function like it precusor but returns the number of genes found near the junction site
    g<-gene.cordinates
    l<- ncol(g)
    n<-Number
    counter=0
    count=0
    while (nrow(JunctRadiusFinder(g, IRinfo, J.pos ,count)) <= n-1){
      count<-count+1
    }
    t<-JunctRadiusFinder(g, IRinfo, J.pos ,count)
    t[is.na(t)]<- "0"
    tup<- matrix(0, n, 3)
    for (i in 1:n){
      dist<-abs(as.numeric(t[i,][2:5])-J)
      if (which(dist==min(dist))>3){
        tup[i,]<-t[i, c(1,4,5)] 
      }
      else {
        tup[i,]<-t[i, c(1,2,3)]
      }
    }
    return(tup)
  }
  
  thFinder<- function(gene.cordinates, IRinfo, J.pos, n){#finding the nth closest gene to the junctions site
    g<-gene.cordinates
    l<- ncol(g)
    if(n==1){
      return(JunctNoFinder(g, IRinfo, J.pos, n))
    }
    else{
      previous<-JunctNoFinder(g, IRinfo, J.pos, n-1)
      current<- JunctNoFinder(g, IRinfo, J.pos, n  )
      return(current[!current[,2] %in% previous[,2],])
    }
  }
  
  RadiusNoFinder<- function(gene.cordinates, IRinfo, J.pos, Number){
    g<-gene.cordinates
    l<- ncol(g)
    n<-Number
    counter=0
    count=0
    while (nrow(JunctRadiusFinder(g, IRinfo, J.pos, count)) <= n-1){
      count<-count+1
    }
    return(count)
  }
  
  GinRadius<- function(Radius, J.pos, track){#thosse who will be plotted
    if(J.pos==1){
      pc<-15
      J<- IRList[[track]][1]
    }
    else if(J.pos==2){
      pc<-40
      J<- IRList[[track]][1]+IRList[[track]][3]
    }
    else if(J.pos==3){
      pc<-70
      J<- IRList[[track]][2]
    }
    else if(J.pos==4){
      pc<-95
      J<- IRList[[track]][2]+IRList[[track]][3]
    }
    else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    t[is.na(t)]<- "0"
    n<- length(t[,1])
    tup<- matrix(0, n, 3)
    for (i in 1:n){
      dist<-abs(as.numeric(t[i,][2:5])-J)
      if (which(dist==min(dist))>3){
        tup[i,]<-t[i, c(1,4,5)] 
      }
      else {
        tup[i,]<-t[i, c(1,2,3)]
      }
    }               
    return(tup)
  }
  
  JG.plotter<- function(Radius, J.pos, track){#JunctionGene plotter
    #' Junction site gene plotter
    #' 
    #' Plotting the genes in the vicinity of the junction site of the chloroplast. 
    #' @param gene.cordinate The output of the gene.cordinate function with the name of the genes and their cordinates on genome
    #' @param junction.cordinate The cordinate of the junction site for which the genes are to be plotted
    #' @param n The number of the genes that is needed to be plotted
    #' @param J.pos The position of the J site: JLB, JSB, JSA, and JLA for 1,2,3, and 4 respectively
    #' @param margin The margin value to be added to the interval as the proportion of the minimum lenght of the nth gene in vicinity
    #' @param track The track on which the genes are to be plotted, strating from the bottom to up as integers 1,2,...
    #' @export
    if(J.pos==1){
      pc<-15
      J<- IRList[[track]][1]
    }
    else if(J.pos==2){
      pc<-40
      J<- IRList[[track]][1]+IRList[[track]][3]
    }
    else if(J.pos==3){
      pc<-70
      J<- IRList[[track]][2]
    }
    else if(J.pos==4){
      pc<-95
      J<- IRList[[track]][2]+IRList[[track]][3]
    }
    else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)    
    t[is.na(t)]<- "0"
    n<- length(t[,1])
    tup<- matrix(0, n, 3)
    for (i in 1:n){
      dist<-abs(as.numeric(t[i,][2:5])-J)
      if (which(dist==min(dist))>3){
        tup[i,]<-t[i, c(1,4,5)] 
      }
      else {
        tup[i,]<-t[i, c(1,2,3)]
      }
    }               
    bw<- 10/Radius
    Rcord<- matrix(0, n, 2)
    Rcord[,1]<-as.numeric(tup[,2])
    Rcord[,2]<-as.numeric(tup[,3])
    Pcord<- (Rcord-J)*bw+pc
    gcol<- function(tup){
      l<-length(tup[,1])
      col<-numeric(l)
      for (i in 1:l){
        if(tup[i,1]=="trnH"){col[i]<- "burlywood4"}
        else if(tup[i,1]=="trnN"){col[i]<- "violet"}
        else if(tup[i,1]=="psbA"){col[i]<- "purple2"}
        else if(tup[i,1]=="rps19"){col[i]<- "firebrick1"}
        else if(tup[i,1]=="ycf1"){col[i]<- "dodgerblue3"}
        else if(tup[i,1]=="ndhF"){col[i]<- "darkred"}
        else if(tup[i,1]=="rpl22"){col[i]<- "navy"}
        else if(tup[i,1]=="rpl2"){col[i]<- "forestgreen"}
        else if(tup[i,1]=="rpl23"){col[i]<- "blue3"}
        else {col[i]<- "black"}
      }
      return(col)
    }
    if (J.pos==4){
      Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
    }
    for (i in 1:n){
      if(Rcord[i,1]>Rcord[i,2]){
        x1<-min(Pcord[i,1], pc+10)
        x2<-max(Pcord[i,2], pc-10)
        for (j in seq(0.10, 0.70, 0.05)){
          segments(x1, track*5+j+5, x2, track*5+j+5, lwd=1, col=paste(gcol(tup)[i]))
        }
      }
      else {
        for (j in seq(1.1, 1.7, 0.05)){
          segments(max(Pcord[i,1], pc-10), track*5-j+5, min(Pcord[i,2], pc+10), track*5-j+5, lwd=1, col=paste(gcol(tup)[i]))
        }
      }
    }
  }
  
  GN.plotter<- function(Radius, J.pos, track){#GeneName plotter
    #' Gene Name plotter
    #' 
    #' Plotting the gene names on a given gene which is already plotted on the tracks of the IR plot
    if(J.pos==1){
      pc<-15
      J<- IRList[[track]][1]
    }
    else if(J.pos==2){
      pc<-40
      J<- IRList[[track]][1]+IRList[[track]][3]
    }
    else if(J.pos==3){
      pc<-70
      J<- IRList[[track]][2]
    }
    else if(J.pos==4){
      pc<-95
      J<- IRList[[track]][2]+IRList[[track]][3]
    }
    else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
    txtincol<- "white"
    txtoutcol="black"
    txtfont<- 4
    numfont<- 1
    numcex<- 0.44
    txtcex<- 0.46
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    t[is.na(t)]<- "0"
    n<- length(t[,1])
    tup<- matrix(0, n, 3)
    for (i in 1:n){
      dist<-abs(as.numeric(t[i,][2:5])-J)
      if (which(dist==min(dist))>3){
        tup[i,]<-t[i, c(1,4,5)] 
      }
      else {
        tup[i,]<-t[i, c(1,2,3)]
      }
    }               
    bw<- 10/Radius
    Rcord<- matrix(0, n, 2)
    Rcord[,1]<-as.numeric(tup[,2])
    Rcord[,2]<-as.numeric(tup[,3])
    Pcord<- (Rcord-J)*bw+pc
    if (J.pos==4){
      Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
    }
    for (i in 1:n){
      if(Rcord[i,1]>Rcord[i,2]){
        x1<-min(Pcord[i,1], pc+10)
        x2<-max(Pcord[i,2], pc-10)
        if (min(x1, x2) >= 104){
          text(103.5, track*5+1.05+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
        }
        else if (abs(x1-x2) >= 9.7){
          text(min(x1, x2)+1.8, track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
          text(max(x1, x2)+1, track*5+0.3+5, paste(Rcord[i,1]-Rcord[i,2], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
        }
        else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
          text(((min(x1, x2)+max(x1, x2))/2), track*5+0.35+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        }
        else if (abs(x1-x2) < 3){
          text(((min(x1, x2)+max(x1, x2))/2), track*5+1.15+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
        }
      }
      else {
        x1<-max(Pcord[i,1], pc-10)
        x2<-min(Pcord[i,2], pc+10)
        if (abs(x1-x2) >= 9.7){
          text(min(x1, x2)+1.8, track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
          text(max(x1, x2)+1, track*5-1.45+5, paste(Rcord[i,2]-Rcord[i,1], "bp", " "), cex=numcex, col=txtincol, font=numfont, pos=2)
        }
        else if (abs(x1-x2) >= 3 & abs(x1-x2) < 9.7){
          text(((min(x1, x2)+max(x1, x2))/2), track*5-1.40+5, tup[i,1], cex=txtcex, col=txtincol, font=txtfont)
        }
        else if (abs(x1-x2) < 3){
          text(((min(x1, x2)+max(x1, x2))/2), track*5-2.25+5, tup[i,1], cex=txtcex, col=txtoutcol, font=txtfont)
        }
      }
    }
  }
  
  OJ.plotter<- function(Radius, J.pos, track){#GeneName plotter
    #' On Junction plotter
    #' 
    #' Plotting the fine tuned narrow lines showing the limits of the genes which are passing through the junction sites with their bp
    if(J.pos==1){
      pc<-15
      J<- IRList[[track]][1]
    }
    else if(J.pos==2){
      pc<-40
      J<- IRList[[track]][1]+IRList[[track]][3]
    }
    else if(J.pos==3){
      pc<-70
      J<- IRList[[track]][2]
    }
    else if(J.pos==4){
      pc<-95
      J<- IRList[[track]][2]+IRList[[track]][3]
    }
    else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
    txtincol<- "white"
    txtoutcol="black"
    txtfont<- 4
    numfont<- 1
    numcex<- 0.44
    txtcex<- 0.46
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    t[is.na(t)]<- "0"
    n<- length(t[,1])
    tup<- matrix(0, n, 3)
    for (i in 1:n){
      dist<-abs(as.numeric(t[i,][2:5])-J)
      if (which(dist==min(dist))>3){
        tup[i,]<-t[i, c(1,4,5)] 
      }
      else {
        tup[i,]<-t[i, c(1,2,3)]
      }
    }               
    bw<- 10/Radius
    Rcord<- matrix(0, n, 2)
    Rcord[,1]<-as.numeric(tup[,2])
    Rcord[,2]<-as.numeric(tup[,3])
    Pcord<- (Rcord-J)*bw+pc
    if (J.pos==4){
      Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
    }
    for (i in 1:n){
      if(Rcord[i,1]>Rcord[i,2]){
        x1<-min(Pcord[i,1], pc+10)
        x2<-max(Pcord[i,2], pc-10)
        if (Rcord[i,1] > J & Rcord[i,2] <J ){
          Arrows(min(x1, x2), track*5+1.1+5, pc-0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          Arrows(pc-0.15, track*5+1.1+5, min(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          Arrows(pc+0.15, track*5+1.1+5, max(x1, x2), track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          Arrows(max(x1, x2), track*5+1.1+5, pc+0.15, track*5+1.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          text(pc-1.4, track*5+1.4+5, paste(Rcord[i, 1]-J+1, "bp", sep=" "), cex=0.4, pos=4)#up-right
          text(pc+1.4, track*5+1.4+5, paste(J-Rcord[i, 2], "bp", sep=" "), cex=0.4, pos=2)#up-left
        }
      }
      else {
        x1<-max(Pcord[i,1], pc-10)
        x2<-min(Pcord[i,2], pc+10)
        if (Rcord[i,2] > J & Rcord[i,1] <J ){
          Arrows(max(x1, x2), track*5-2.1+5, pc+0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          Arrows(pc+0.15, track*5-2.1+5, max(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          Arrows(pc-0.15, track*5-2.1+5, min(x1, x2), track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          Arrows(min(x1, x2), track*5-2.1+5, pc-0.15, track*5-2.1+5, arr.type = "T", cex=0.5, arr.length = 0.12, lwd=0.5, arr.width = 0.4)
          text(pc-1.4, track*5-2.6+5, paste(Rcord[i, 2]-J, "bp", sep=" "), cex=0.4, pos=4)#low-right
          text(pc+1.4, track*5-2.6+5, paste(J-Rcord[i, 1], "bp", sep=" "), cex=0.4, pos=2)#low-left
        }
      }
    }
  }
  
  JD.plotter<- function(Radius, J.pos, track){#GeneName plotter
    #' Junction Distance plotter
    #' 
    #' plotting the narrow lines of the distance of the genes for the junction sites which are not passing through any gene and their bp
    if(J.pos==1){
      pc<-15
      J<- IRList[[track]][1]
    }
    else if(J.pos==2){
      pc<-40
      J<- IRList[[track]][1]+IRList[[track]][3]
    }
    else if(J.pos==3){
      pc<-70
      J<- IRList[[track]][2]
    }
    else if(J.pos==4){
      pc<-95
      J<- IRList[[track]][2]+IRList[[track]][3]
    }
    else {stop("J.pos missing or out bound. It should be either 1,2,3, or 4 for JLB, JSB, JSA, and JLA, recpectively ")}
    txtincol<- "white"
    txtoutcol="black"
    txtfont<- 4
    numfont<- 1
    numcex<- 0.44
    txtcex<- 0.46
    t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
    t[is.na(t)]<- "0"
    n<- length(t[,1])
    tup<- matrix(0, n, 3)
    for (i in 1:n){
      dist<-abs(as.numeric(t[i,][2:5])-J)
      if (which(dist==min(dist))>3){
        tup[i,]<-t[i, c(1,4,5)] 
      }
      else {
        tup[i,]<-t[i, c(1,2,3)]
      }
    }               
    bw<- 10/Radius
    Rcord<- matrix(0, n, 2)
    Rcord[,1]<-as.numeric(tup[,2])
    Rcord[,2]<-as.numeric(tup[,3])
    Pcord<- (Rcord-J)*bw+pc
    if (J.pos==4){
      ind<-Pcord < 0
      Pcord[Pcord < 0]<- ((Rcord[which(Pcord< 0)] + IRList[[track]][4])-J)*bw+pc
      Rcord[ind]<- Rcord[ind] + IRList[[track]][4]
    }
    counter<- 0
    for (i in 1:n){
      if(Rcord[i,1]>Rcord[i,2]){
        x1<-min(Pcord[i,1], pc+10)
        x2<-max(Pcord[i,2], pc-10)
        if (Rcord[i,1] > J & Rcord[i,2] <J ){
          counter<- counter+1
        }
      }
      else {
        x1<-max(Pcord[i,1], pc-10)
        x2<-min(Pcord[i,2], pc+10)
        if (Rcord[i,2] > J & Rcord[i,1] <J ){
          counter<- counter+1
        }
      }
    }
    if (counter==0){#find the closest gene to the junction site
      nearest<- which(abs(Pcord-pc)==min(abs(Pcord-pc)))
      col.cor<- floor((nearest-0.01)/n)+1
      row.cor<- nearest-n*floor((nearest-0.01)/n)###Now we have the row of the nearest gene
      #with this setting the position of the zero (genes tangant to the junction site) will not be plotted. If interested either put "=" for the middle condition of the left or right binary operator or better develop zero only handling if function
      if     (Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top,right, big
        curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5+1.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, right, small
        arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(Rcord[row.cor, 2]-J, "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, right, big
        curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 + 3, track*5-2.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE) 
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4, track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] < 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, right, small
        arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 -3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 -4.5 , track*5-2.3-0.2+5, paste(Rcord[row.cor, 1]-J, "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) > 3 ){#low, left, big
        curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5-2.3+5), curve = -0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE) 
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4, track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] < Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & abs(max(Pcord[row.cor, 1], pc-10) - min(Pcord[row.cor, 2], pc+10)) <= 3 ){#low, left, small
        arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5-1.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 3 , track*5-2.3+5, angle = 15, length = 0.05, lwd=0.6)
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4.5 , track*5-2.3-0.2+5, paste(J-Rcord[row.cor, 2], "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) > 3 ){#top, left, big
        curvedarrow (from=c(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5), to=c(pc+(Pcord[row.cor, col.cor]-pc)/2 - 3, track*5+1.3+5), curve = 0.21, lwd=0.6, arr.type = "curved", arr.col = "white", arr.length=0.08, arr.lwd=0.4, arr.pos=0.69, endhead=TRUE)
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 - 4.5 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
      }
      else if(Rcord[row.cor, 1] > Rcord[row.cor, 2] & pc-Pcord[row.cor, col.cor] > 0 & min(Pcord[row.cor, 1], pc+10) - max(Pcord[row.cor, 2], pc-10) <= 3 ){#top, left, small
        arrows(pc+(Pcord[row.cor, col.cor]-pc)/2, track*5+0.3+5 , pc+(Pcord[row.cor, col.cor]-pc)/2 + 2 , track*5+1.3+5, angle = 15, length = 0.05, lwd=0.6)
        text(pc+(Pcord[row.cor, col.cor]-pc)/2 + 4 , track*5+1.3+0.2+5, paste(J-Rcord[row.cor, 1], "bp", sep=" "), cex=0.4)
      }
    }
  }
  
  Max.Radius<-function(J.pos){
    if(J.pos==1){
      Radius<-680
    }
    if(J.pos==2){
      Radius<-100
    }
    if(J.pos==3){
      Radius<-1800
    }
    if(J.pos==4){
      Radius<-1000
    }
    R<- numeric(l)
    for (track in 1:l){
      t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
      while(nrow(t)==0){
        Radius<- Radius+(0.2)*Radius
        t<- JunctRadiusFinder(GeneList[[track]], IRList[[track]], J.pos , Radius)
      }
      R[track]<-Radius
    }
    return(round(max(R)+1))
  }
  
  trnHfixer<- function(Genelist){
    for (i in 1:length(Genelist[,1])){
      if (Genelist[i,1] %in% c("trnH-GUG", "tRNA-His")){
        Genelist[i,1]<- "trnH"
      }
    }
    return(Genelist)
  }
  
  trnNfixer<- function(Genelist){
    for (i in 1:length(Genelist[,1])){
      if (Genelist[i,1] %in% c("trnN-GUU", "tRNA-Asn")){
        Genelist[i,1]<- "trnN"
      }
    }
    return(Genelist)
  }
  
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

mboard <- function(gbFile, rad=100, cexlab=0.4, gene.width=5, offset=1, gay=TRUE, intron.width=3){
  gb<- read.gb(gbFile)
  genome<- rdnFixer(gb)
  gene.tab<- gene.cordinates(gb)
  IR<- IRinfo(genome)
  L<- length(genome)
  Nrg<- length(gene.tab[,1])
  RmidIR<- IR[1]+IR[3]+(IR[2]-(IR[1]+IR[3]))/2#it can be either full or half number  #mid SSC cordinate should move to hL
  hL<- L/2#can be round or half length
  shifter<- RmidIR-hL-0.1
  #shifter<- RmidIR+hL-0.1#shifter is making the genes exactly on the line of it to the first trigonometric area!
  IR[1]<- IR[1]-shifter
  IR[2]<- IR[2]-shifter
  if (IR[1]<=0){IR[1]<- IR[1]+L}
  if (IR[2]<=0){IR[2]<- IR[2]+L}
  nc<- ncol(gene.tab)-1
  m<- matrix(0, Nrg, nc)
  m[, 1:nc]<- as.numeric(gene.tab[,2:(nc+1)])
  m<- m-shifter
  m[m[,]==-shifter]<- 0
  m[m[,]<0]<- m[m[,]<0]+L
  t<- m
  t[t>0]<- pi/2+2*pi*m[m>0]/L  #theta convertor of the cordinate  pi/2+2*pi*cordinate/total length
  t[t==0]<- NA
  for (i in 1:Nrg){
    if (abs(t[i, 1]-t[i,2])> pi ){
      ph<-min(t[i,2], t[i,1])
      t[i,1]<- min(t[i,2], t[i,1])+(2*pi-abs(t[i,2]- t[i,1]))
      t[i,2]<-ph
    }
  }
  x<- cos(t)*rad
  y<- sin(t)*rad
  op=par(mai=c(0,0,0,0),bg = "white")
  if (gay){pal<- sample(1:999, Nrg)} else {pal<- rep(Nrg, "black")}
  pal<- sample(1:999, Nrg)
  #Draw a circle
  #dev.new(8.3, 8.3)
  jpeg("Circ.jpg",  width=8.3, height=8.9, units="in", res=1000)
  plot(NULL,xlim=c(-115,115),ylim=c(-115,115),axes=FALSE,ann=FALSE)
  radi<- rad-1
  rado<- rad+1
  draw.circle(0,0, 100, nv=70000)
  draw.circle(0, 0, 45, nv=70000)
  draw.circle(0, 0, 55, nv=70000)
  draw.circle(0,0, 0.01, lwd = 2)
  (t[,1]+t[,2])/2 -> m1t
  (t[,3]+t[,4])/2 -> m2t##and so on for other columns of the genetable if they existed 
  int.table<- intronCordinates(gb)
  gw<- gene.width
  iw<- intron.width
  for (i in 1:Nrg){#gene plotting
    if (m[i,1]<m[i,2]){Arc.plotter(t[i,1], t[i,2], 100, 100+gw, col=pal[i], border=T, lwd=0.4)} #forward is out circle
    else if (m[i,3]<m[i,4]){Arc.plotter(t[i,3], t[i,4], 100, 100+gw, col=pal[i], border=T, lwd=0.4)}
    else if(m[i,1]>m[i,2]){Arc.plotter(t[i,1], t[i,2], 100-gw, 100, col=pal[i], border=T, lwd=0.4)}
    else if (m[i,3]>m[i,4]){Arc.plotter(t[i,3], t[i,4], 100-gw, 100, col=pal[i], border=T, lwd=0.4)}
  }
  #for (i in 1:length(int.table[,1])){
   # if (int.table[i,1]<int.table[i,2]){Arc.plotter(int.table[i,1], int.table[i,2], 100, 100+iw, col="white", border = F)}
  #  else                              {Arc.plotter(int.table[i,1], int.table[i,2], 100, 100-iw, col="white", border = F)}  
  #}
  IRtheta<- pi/2+2*pi*IR/L #pi/2+2*pi*cordinate/total
  draw.arc(x=0, y=0, radius = 100, angle1= IRtheta[1], angle2= IRtheta[1]+IRtheta[3]-pi/2, col = "orange", lwd = 3)
  draw.arc(x=0, y=0, radius = 100, angle1= IRtheta[2], angle2= IRtheta[2]+IRtheta[3]-pi/2, col = "orange", lwd = 3)
  Arc.plotter( 0, 2*pi, 53, 57, col="lightblue", border=F, edges = 400, lwd=0.4)
  Arc.plotter( IRtheta[1], IRtheta[1]+IRtheta[3]-pi/2, 53, 57, col="blue", border=F, edges = 100, lwd=0.4)
  Arc.plotter( IRtheta[2], IRtheta[2]+IRtheta[3]-pi/2, 53, 57, col="blue", border=F, edges = 100, lwd=0.4)
  Arc.plotter( IRtheta[1]+IRtheta[3]-pi/2, IRtheta[2], 53, 57, col="lightgreen", border=F, edges = 100, lwd=0.4)
    for (i in 1:Nrg){#annotating
    c<-(100-gw)/100-(offset)/100
    if (m[i,1]>m[i,2]){
      theta<- m1t[i]
      if ((m1t[i]>=pi/2)&(m1t[i]<3*pi/2)){#inner left
        theta<- m1t[i]+0.007
        ttext<-pi+theta
        text(c*radi*cos(theta),
             c*radi*sin(theta),
             gene.tab[i,1],
             col="black",
             cex=cexlab,
             srt=360*ttext/(2*pi),
             adj=c(0,0))
      }
      else{                               #inner right
        theta<- m1t[i]-0.007
        ttext=theta
        dist=c*radi-strwidth(as.character(gene.tab[i,1]))*cexlab
        text(dist*cos(theta),dist*sin(theta),gene.tab[i,1],col="black",cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
    }
    if (m[i,3]>m[i,4]){                   #inner left
      theta<- m2t[i]
      if ((m2t[i]>=pi/2)&(m2t[i]<3*pi/2)){
        theta<- m2t[i]+0.007
        ttext<-pi+theta
        text(c*radi*cos(theta),c*radi*sin(theta),gene.tab[i,1],col="black",cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
      else{                               #inner right
        theta<- m2t[i]-0.007
        ttext=theta
        dist=c*radi-strwidth(as.character(gene.tab[i,1]))*cexlab
        text(dist*cos(theta),dist*sin(theta),gene.tab[i,1],col="black",cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
    }
    c<- (100+gw)/100+offset/100
    if (m[i,1]<m[i,2]){
      theta<- m1t[i]
      if ((m1t[i]>=pi/2)&(m1t[i]<3*pi/2)){#outter left
        theta<- m1t[i]+0.007
        ttext<-pi+theta
        dist=c*radi+(strwidth(as.character(gene.tab[i,1])))*cexlab+(strwidth("a")+1)*cexlab
        text(dist*cos(theta),dist*sin(theta),gene.tab[i,1],col="black",cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
      else{                               #outter right
        theta<- m1t[i]-0.007
        ttext=theta
        dist=c*radi+(strwidth("a")+1)*cexlab
        text(dist*cos(theta),dist*sin(theta),gene.tab[i,1],col="black",cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
    }
    if (m[i,3]<m[i,4]){                   #outter left
      if ((m2t[i]>=pi/2)&(m2t[i]<3*pi/2)){
        theta<- m2t[i]+0.007
        ttext<-pi+theta
        dist=c*radi+(strwidth(as.character(gene.tab[i,1])))*cexlab+(strwidth("a")+1)*cexlab
        text(dist*cos(theta),dist*sin(theta),gene.tab[i,1],col="black",cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
      else{                               #outter right
        theta<- m2t[i]-0.007
        ttext=theta
        dist=c*radi+(strwidth("a")+1)*cexlab
        text(dist*cos(theta),dist*sin(theta),gene.tab[i,1],col=222,cex=cexlab,srt=360*ttext/(2*pi),adj=c(0,0))
      }
    }
  }
  text(0,8, sp.name(gb), font=4, cex = 1)
  text(0,-3, paste(prettyNum(L, big.mark = ","), "bp", " "),font=4, cex = 1)
  dev.off()
  dev.off()
}
#}

