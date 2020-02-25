compareDNA <- function(x,y){
  as.integer(x) == as.integer(y)
}

irDetect <- function(genome) {
  # detect shifter
  genome_rc <- Biostrings::reverseComplement(genome)
  l <- Biostrings::nchar(genome_rc)
  ir <- genome_rc[1:1000]
  m <- Biostrings::matchPattern(ir,genome,max.mismatch = 50,with.indels = TRUE)
  if (length(m) < 1){
    ir_table <- data.frame(chr = "chr1",
                           start = 0,
                           end = l,
                           center = round(l/2),
                           name = "LSC",
                           text = paste("LSC:", l),
                           # y_t = 0.25,
                           # y_b = -0.25,
                           stringsAsFactors = FALSE)
    return(ir_table)
  }
  shifter <- m@ranges@start - 1

  # IRA start, end and lenght
  for (i in 1:l){
    true.search<-compareDNA(cirtick((shifter - 50 + i), genome),
                            genome_rc)
    count <- rle(true.search)
    count_T <- count$lengths[count$values]
    if ((sum(count_T) - l/4)/l > 0.07) { #length of IR region
      break()
    }
  }



  count_F <- count$lengths[!count$values]
  count_len <- length(count_T)
  ira_s <- shifter - 50 + i + 1
  ira_len <- count$lengths[1]
  if (tail(count$values, 1)){
    ira_len <- ira_len + tail(count$lengths, 1)
    ira_s <- ira_s - tail(count$lengths, 1)
    count_T <- count_T
  }


  for (i in 2:count_len) {
    if (count_T[i] >= 16 & count_F[i] <= 3){
      ira_len <- ira_len + count_T[i]
    } else {
      break()
    }
  }

  ira_e <- ira_s + ira_len -1

  # IRB start, end and length

  irb_len <- sort(count_T, decreasing = TRUE)[1:2]

  if (min(irb_len) > 20000){
    pos <- which(count_T %in% irb_len)
    irb_len <- count_T[max(pos)]
    if (length(pos) %in% c(1,2)){
      irb_s <- sum(count$lengths[1:(max(pos) - 1)]) + shifter + 1
    }
  } else {
    warning("Didn't find IRB")
  }

  for (i in (max(pos) + 1):count_len) {
    if (count_T[i] >= 16 & count_F[i] <= 3){
      irb_len <- irb_len + count_T[i]
    } else {
      break()
    }
  }
  irb_e <- irb_s + irb_len - 1

  # IR information table
  if (irb_e == l){
    ir_table <- data.frame(chr = rep("chr1", 4),
                          start = c(0, ira_s, ira_e + 1, irb_s),
                          end = c(ira_s - 1, ira_e, irb_s -1, irb_e),
                          name = c("LSC", "IRA", "SSC", "IRB"),
                          text = c(paste("LSC:", ira_s -1),
                                   paste("IRA:", ira_len),
                                   paste("SSC:", irb_s - 1 - ira_e),
                                   paste("IRB:", irb_len)),
                          # y_b = rep(-0.25, 4),
                          # y_t = rep(0.25, 4),
                          stringsAsFactors = FALSE)
    ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
  } else if (irb_e < l) {
    ir_table <- data.frame(chr = rep("chr1", 5),
                           start = c(0, ira_s, ira_e + 1, irb_s, irb_e + 1),
                           end = c(ira_s - 1, ira_e, irb_s -1, irb_e, l),
                           name = c("LSC", "IRA", "SSC", "IRB", "LSC"),
                           text = c(paste("LSC:", ira_s -1 + l - irb_e),
                                    paste("IRA:", ira_len),
                                    paste("SSC:", irb_s - 1 - ira_e),
                                    paste("IRB:", irb_len),
                                    ""),
                           # y_b = rep(-0.25, 5),
                           # y_t = rep(0.25, 5),
                           stringsAsFactors = FALSE)
    ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
    ir_table$center[1] <- ir_table$center[1] - l + ir_table$center[5]
  } else if (irb_e > l) {
    irb_e <- irb_e - l
    ir_table <- data.frame(chr = rep("chr1", 5),
                           start = c(0, irb_e, ira_s, ira_e + 1, irb_s),
                           end = c(irb_e, ira_s - 1, ira_e, irb_s -1, l),
                           name = c("IRB", "LSC", "IRA", "SSC", "IRB"),
                           text = c("",
                                    paste("LSC:", ira_s -1),
                                    paste("IRA:", ira_len),
                                    paste("SSC:", irb_s - 1 - ira_e),
                                    paste("IRB:", irb_len)),
                           # y_b = rep(-0.25, 5),
                           # y_t = rep(0.25, 5),
                           stringsAsFactors = FALSE)
    ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
    ir_table$center[1] <- ir_table$center[1] + ir_table$center[5]
  }

  return(ir_table)
}
