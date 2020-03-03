compareDNA <- function(x,y){
  as.integer(x) == as.integer(y)
}

irDetect <- function(genome) {
  # detect shifter
  genome_rc <- Biostrings::reverseComplement(genome)
  l <- Biostrings::nchar(genome_rc)

  s <- seq(1, (l - 1000), 100)
  seeds <- Biostrings::DNAStringSet(genome, start = s, width = 1000)
  seeds <- PDict(seeds)
  m <- Biostrings::matchPDict(seeds, genome_rc)
  m <- as.data.frame(m)
  m <- m[m$width > 900, ]
  if (nrow(m) == 0){
  ir_table <- data.frame(chr = "chr1",
                         start = 0,
                         end = l,
                         center = round(l/2),
                         name = "LSC",
                         text = paste("LSC:", l),
                         stringsAsFactors = FALSE)
    return(ir_table)
  }

  # IRA start, end and lenght
  for (i in 1:nrow(m)){
    shifter <- m$start[i] - s[m$group[i]]
    if (shifter < 0){
      shifter <- l + shifter
    }
    true.search<-compareDNA(cirtick((shifter), genome_rc), genome)
    count <- rle(true.search)
    count_T <- count$lengths[count$values]
    ir_len <- sort(count_T, decreasing = TRUE)[1:2] #length of IR region
    if (all(ir_len > 1000)) {
      break()
    }
  }

  if (!all(ir_len > 1000)){
    ir_table <- data.frame(chr = "chr1",
                           start = 0,
                           end = l,
                           center = round(l/2),
                           name = "LSC",
                           text = paste("LSC:", l),
                           stringsAsFactors = FALSE)
    return(ir_table)
  }

  count_len <- length(count$lengths)
  pos <- which(count$lengths %in% ir_len)

  # IRA start, end and length (0-base coordinate system)

  ira_s2 <- NA
  ira_len <- count$lengths[min(pos)]
  if (min(pos) == 1){
    ira_s <- 0
  } else {
    ira_s <- sum(count$lengths[1:(min(pos) - 1)])
  }

  # Enable SNP in IRA
  # trace forward
  for (i in seq((min(pos) + 2), count_len, by = 2)) {
    if (count$lengths[i] >= 16 & count$lengths[i - 1] <= 3){
      ira_len <- ira_len + count$lengths[i] + count$lengths[i - 1]
    } else {
      break()
    }
  }

  # trace backward
  if (min(pos) < 10){
    count_lengths <- c(count$lengths[11:count_len], count$lengths[1:10])
    count_values <- c(count$values[11:count_len], count$values[1:10])
    pos_new <- min(pos) + count_len - 10
    if (count_values[pos_new - 1]){
      ira_len <- ira_len + count_lengths[pos_new - 1]
      ira_s <- ira_s - count_lengths[pos_new - 1]
      for (i in seq(pos_new - 3, 10, by = - 2)) {
        if (count_lengths[i] >= 16 & count_lengths[i + 1] <= 3){
          ira_len <- ira_len + count_lengths[i] + count_lengths[i + 1]
          ira_s <- ira_s - (count_lengths[i] + count_lengths[i + 1])
        } else {
          break()
        }
      }
    } else {
      for (i in seq(pos_new - 2, 10, by = - 2)) {
        if (count_lengths[i] >= 16 & count_lengths[i + 1] <= 3){
          ira_len <- ira_len + count_lengths[i] + count_lengths[i + 1]
          ira_s <- ira_s - (count_lengths[i] + count_lengths[i + 1])
        } else {
          break()
        }
      }
    }
  } else{
    for (i in seq(min(pos) - 2, 0, by = - 2)) {
      if (count$lengths[i] >= 16 & count$lengths[i + 1] <= 3){
        ira_len <- ira_len + count$lengths[i] + count$lengths[i + 1]
        ira_s <- ira_s - (count$lengths[i] + count$lengths[i + 1])
      } else {
        break()
      }
    }
  }

  ira_e <- ira_s + ira_len

  if(ira_s < 0) {
    ira_s2 <-  l + ira_s
    ira_e <- ira_len + ira_s
    ira_s <- 0
  }

  # IRB start, end and length (0-base coordinate system)

  irb_e2 <- NA
  irb_len <- count$lengths[max(pos)]
  irb_s <- sum(count$lengths[1:(max(pos) - 1)])

  # Enable SNP in IRB
  # trace forward
  if ((count_len - max(pos)) < 10){
    count_lengths <- c(count$lengths[11:count_len], count$lengths[1:10])
    count_values <- c(count$values[11:count_len], count$values[1:10])
    pos_new <- max(pos) - 10
    if (count_values[pos_new + 1]){
      irb_len <- irb_len + count_lengths[pos_new + 1]
      irb_s <- irb_s - count_lengths[pos_new + 1]
      for (i in seq(pos_new + 3, count_len, by = 2)) {
        if (count_lengths[i] >= 16 & count_lengths[i - 1] <= 3){
          irb_len <- irb_len + count_lengths[i] + count_lengths[i - 1]
        } else {
          break()
        }
      }
    } else {
      for (i in seq(pos_new + 2, count_len, by = 2)) {
        if (count_lengths[i] >= 16 & count_lengths[i - 1] <= 3){
          irb_len <- irb_len + count_lengths[i] + count_lengths[i - 1]
        } else {
          break()
        }
      }
    }
  } else{
    for (i in seq(max(pos) + 2, count_len, by = 2)) {
      if (count$lengths[i] >= 16 & count$lengths[i - 1] <= 3){
        irb_len <- irb_len + count$lengths[i] + count$lengths[i - 1]
      } else {
        break()
      }
    }
  }

  # trace backward
  for (i in seq((max(pos) - 2), 0, by = - 2)) {
    if (count$lengths[i] >= 16 & count$lengths[i + 1] <= 3){
      irb_len <- irb_len + count$lengths[i] + count$lengths[i + 1]
      irb_s <-  irb_s - (count$lengths[i] + count$lengths[i + 1])
    } else {
      break()
    }
  }

  irb_e <- irb_s + irb_len

  if (irb_e > l){
    irb_e2 <- irb_e - l
    irb_e <- l
  }

  # IR information table
  ir_table <- data.frame(chr = rep("chr1", 5),
                         start = c(0, ira_s, ira_e + 1, irb_s, irb_e + 1),
                         end = c(ira_s - 1, ira_e, irb_s -1, irb_e, l),
                         name = c("LSC", "IRA", "SSC", "IRB", "LSC"),
                         text = c(paste("LSC:", ira_s -1 + l - irb_e),
                                  paste("IRA:", ira_len),
                                  paste("SSC:", irb_s - 1 - ira_e),
                                  paste("IRB:", irb_len),
                                  ""),
                         stringsAsFactors = FALSE)
  ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
  ir_table$center[1] <- ir_table$center[1] + round(l/2) + round(ir_table$start[5]/2)
  if(ir_table$center[1] > l){
    ir_table$center[1] <- ir_table$center[1] - l
  }

  if (ira_s == 0) {
    ir_table$text[5] <- ir_table$text[1]
    ir_table <- ir_table[-1, ]
  }
  if (irb_e == l) {
    ir_table <- ir_table[-5, ]
  }

  if (!is.na(irb_e2)){
    df <- data.frame(chr = "chr1", start = 0, end = irb_e2, name = "IRB",
                     text = "",center = 0, stringsAsFactors = FALSE)
    ir_table <- rbind.data.frame(df, ir_table)
    ir_table$start[2] <- irb_e2 + 1
    ir_table$center[2] <- ir_table$center[2] - round(irb_e2/2)
    ir_table$center[5] <- ir_table$center[5] + round(irb_e2/2)
    if (ir_table$center[5] > l) {
      ir_table$center[5] <- ir_table$center[5] - l
    }
  }
  if (!is.na(ira_s2)){
    df <- data.frame(chr = "chr1", start = ira_s2, end =  l, name = "IRA",
                     text = "",center = 0, stringsAsFactors = FALSE)
    ir_table <- rbind.data.frame(ir_table, df)
    ir_table$end[4] <- ira_s2 - 1
    ir_table$center[1] <- ir_table$center[1] - round((l-ira_s2)/2)
    if (ir_table$center[1] < 0) {
      ir_table$center[1] <- l + ir_table$center[1]  - 1
    }
  }

  # if (irb_e = l){
  #   ir_table <- data.frame(chr = rep("chr1", 4),
  #                         start = c(0, ira_s, ira_e + 1, irb_s),
  #                         end = c(ira_s - 1, ira_e, irb_s -1, irb_e),
  #                         name = c("LSC", "IRA", "SSC", "IRB"),
  #                         text = c(paste("LSC:", ira_s -1),
  #                                  paste("IRA:", ira_len),
  #                                  paste("SSC:", irb_s - 1 - ira_e),
  #                                  paste("IRB:", irb_len)),
  #                         stringsAsFactors = FALSE)
  #   ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
  # } else if (irb_e < l) {
  #   ir_table <- data.frame(chr = rep("chr1", 5),
  #                          start = c(0, ira_s, ira_e + 1, irb_s, irb_e + 1),
  #                          end = c(ira_s - 1, ira_e, irb_s -1, irb_e, l),
  #                          name = c("LSC", "IRA", "SSC", "IRB", "LSC"),
  #                          text = c(paste("LSC:", ira_s -1 + l - irb_e),
  #                                   paste("IRA:", ira_len),
  #                                   paste("SSC:", irb_s - 1 - ira_e),
  #                                   paste("IRB:", irb_len),
  #                                   ""),
  #                          stringsAsFactors = FALSE)
  #   ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
  #   ir_table$center[1] <- ir_table$center[1] - l + ir_table$center[5]
  # } else if (irb_e > l) {
  #   irb_e <- irb_e - l
  #   ir_table <- data.frame(chr = rep("chr1", 5),
  #                          start = c(0, irb_e, ira_s, ira_e + 1, irb_s),
  #                          end = c(irb_e, ira_s - 1, ira_e, irb_s -1, l),
  #                          name = c("IRB", "LSC", "IRA", "SSC", "IRB"),
  #                          text = c("",
  #                                   paste("LSC:", ira_s -1),
  #                                   paste("IRA:", ira_len),
  #                                   paste("SSC:", irb_s - 1 - ira_e),
  #                                   paste("IRB:", irb_len)),
  #                          stringsAsFactors = FALSE)
  #   ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
  #   ir_table$center[1] <- ir_table$center[1] + ir_table$center[5]
  # }

  return(ir_table)
}
