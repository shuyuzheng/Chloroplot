compareDNA <- function(x,y){
  as.integer(x) == as.integer(y)
}
#circular rotative function
cirtick <- function(tick, vector){
  if(tick > length(vector)-1 || tick < 1){
    return(vector)
  }
  else {
    return(c(vector[(tick+1):length(vector)], vector[1:tick]))
  }
}
irDetect <- function(genome) {
  # detect shifter
  genome_rc <- Biostrings::reverseComplement(genome)
  l <- Biostrings::nchar(genome_rc)

  s <- seq(1, (l - 1000), 100)
  seeds <- Biostrings::DNAStringSet(genome, start = s, width = 1000)
  seeds <- Biostrings::PDict(seeds)
  m <- Biostrings::matchPDict(seeds, genome_rc)
  m <- as.data.frame(m)
  m <- m[m$width > 900, ]
  if (nrow(m) == 0){
  gc_count <- Biostrings::letterFrequency(genome, letters = c("C", "G"))
  gc_count <- round((gc_count[1] + gc_count[2])/l)
  ir_table <- data.frame(chr = "chr1",
                         start = 0,
                         end = l,
                         center = round(l/2),
                         name = "LSC",
                         text = paste("LSC:", l),
                         gc_count = gc_count,
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
    gc_count <- Biostrings::letterFrequency(genome, letters = c("C", "G"))
    gc_count <- round((gc_count[1] + gc_count[2])/l)
    ir_table <- data.frame(chr = "chr1",
                           start = 0,
                           end = l,
                           center = round(l/2),
                           name = "LSC",
                           text = paste("LSC:", l),
                           gc_count = gc_count,
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
  return(ir_table)
}

irBounder <- function(genome, shifter = 0){

  l <- Biostrings::nchar(genome)
  if (shifter == 0){
    g <- genome
  } else if(shifter > 0) {
    g <- c(genome[(shifter + 1):l], genome[1:shifter])
  } else if (shifter < 0) {
    g <- c(genome[(l + shifter + 1):l], genome[1:(l + shifter)])
  }

  g_rc<- Biostrings::reverseComplement(g)
  # split exanded genome to 1000 bp small pieces and match with reverse
  # complemented expanded genome
  s <- seq(1, (l - 999))
  seeds <- Biostrings::DNAStringSet(g, start = s, width = 1000)
  seeds <- Biostrings::PDict(seeds)
  m <- Biostrings::matchPDict(seeds, g_rc)
  m <- as.data.frame(m)
  m <- m[m$width > 0, ]

  if (nrow(m) == 0){
    gc_count <- Biostrings::letterFrequency(genome, letters = c("C", "G"))
    gc_count <- round((gc_count[1] + gc_count[2])/l)
    ir_table <- data.frame(chr = "chr1",
                           start = 0,
                           end = l,
                           center = round(l/2),
                           name = "LSC",
                           text = paste("LSC:", l),
                           gc_count = gc_count,
                           stringsAsFactors = FALSE)
    return(ir_table)
  }

  # IR start, end and length

  m <- m[!duplicated(m$group), ]
  m_l <- nrow(m)
  step <- m$group[-1] - m$group[-m_l]
  count <- rle(step)
  index <- which(count$lengths %in% sort(count$lengths, decreasing = TRUE)[1:2])
  if (min(index) == 1) {
    ira_s <- m$group[1] - 1
    ira_e <- m$group[count$lengths[1] + 1] + 999
  } else {
    ira_s_index <- sum(count$lengths[1:(min(index) - 1)]) + 1
    ira_e_index <- sum(count$lengths[1:min(index)])
    for (i in seq(min(index) - 2, 1, -2)){
      if (count$values[i + 1] < 1500){
        ira_s_index <- ira_s_index - 1 - count$lengths[i]
      } else {
        break()
      }
    }
    ira_s <- m$group[ira_s_index] - 1
    ira_e <- m$group[ira_e_index + 1] + 999
  }

  if (max(index) == length(count$lengths)) {
    irb_e <- m$group[m_l] + 999
    irb_s <- m$group[sum(count$lengths[1:(max(index) - 1)]) + 1] - 1
  } else {
    irb_e_index <- sum(count$lengths[1:max(index)])
    for (i in seq(max(index) + 2, length(count$lengths), 2)){
      if (count$values[i - 1] < 1500){
        irb_e_index <- irb_e_index + 1 + count$lengths[i]
      } else {
        break()
      }
    }
    irb_e <- m$group[irb_e_index + 1] + 999
    irb_s <- m$group[sum(count$lengths[1:(max(index) - 1)]) + 1] - 1
  }
  indel_table <- NULL

  # mismatch in IRA
  ira_step <- step[which(m$group > ira_s & m$group < (ira_e - 999))]
  ira_mismatch <- which(ira_step > 1)
  if (length(ira_mismatch) > 0){
    ira_mismatch <- ira_mismatch + length(step[which(m$group <= ira_s)])
    start <- m$group[ira_mismatch] + 999
    end <- m$group[ira_mismatch + 1] - 1
    len_jump <- (m$group[ira_mismatch + 1] - m$group[ira_mismatch])
    len_jump_rc <- (m$start[ira_mismatch  + 1] - m$start[ira_mismatch])
    type <- rep("replace", length(ira_mismatch))
    type[which(len_jump < len_jump_rc)] <- "delet"
    type[which(len_jump > len_jump_rc)] <- "insert"
    indel_table_a <- data.frame(name = type, start = start, end = end,
                                stringsAsFactors = FALSE)
    indel_table_a$text <- rep("", nrow(indel_table_a))
    indel_table_a$center <- round((indel_table_a$start + indel_table_a$end)/2)
    indel_table_a$gc_count <- rep("", nrow(indel_table_a))
    indel_table_a$chr <- rep("chr1", nrow(indel_table_a))
    indel_table <- rbind.data.frame(indel_table, indel_table_a)
  }

  # mismatch in IRB
  irb_step <- step[which(m$group > irb_s & m$group < (irb_e - 999))]
  irb_mismatch <- which(irb_step > 1)
  if (length(irb_mismatch) > 0){
    irb_mismatch <- irb_mismatch + length(step[which(m$group <= irb_s)])
    start <- m$group[irb_mismatch] + 999
    end <- m$group[irb_mismatch + 1] - 1
    len_jump <- (m$group[irb_mismatch + 1] - m$group[irb_mismatch])
    len_jump_rc <- (m$start[irb_mismatch + 1] - m$start[irb_mismatch])
    type <- rep("replace", length(irb_mismatch))
    type[which(len_jump < len_jump_rc)] <- "delet"
    type[which(len_jump > len_jump_rc)] <- "insert"
    indel_table_b <- data.frame(name = type, start = start, end = end,
                                stringsAsFactors = FALSE)
    indel_table_b$text <- rep("", nrow(indel_table_b))
    indel_table_b$center <- round((indel_table_b$start + indel_table_b$end)/2)
    indel_table_b$gc_count <- rep("", nrow(indel_table_b))
    indel_table_b$chr <- rep("chr1", nrow(indel_table_b))
    indel_table <- rbind.data.frame(indel_table, indel_table_b)
  }

  if (!is.null(indel_table)){
    for (i in 1:nrow(indel_table)){
      if (indel_table$start[i] > indel_table$end[i]){
        s <- g[indel_table$end[i]:indel_table$start[i]]
        if (indel_table$name[i] == "delet" & (length(Biostrings::uniqueLetters(s[-1]))==1)){
          indel_table$end[i] <- indel_table$start[i]
          indel_table$center[i] <- indel_table$start[i]
        } else if (indel_table$name[i] == "insert" & (length(Biostrings::uniqueLetters(s))==1)){
          indel_table$end[i] <- indel_table$start[i] + 1
          indel_table$center[i] <- indel_table$start[i]
        }
      }
    }
    indel_table$start <- indel_table$start + shifter
    indel_table$end <- indel_table$end + shifter
    indel_table$center <- indel_table$center + shifter
  }

  res <- c(ira_s, ira_e, irb_s, irb_e)
  names(res) <- c("ira_s", "ira_e", "irb_s", "irb_e")
  res <- res + shifter

  result <- list(ir_bounder = res, indel_table = indel_table)
  return(result)
}

irDetect_indel <- function(genome) {
  # expand genome with 1000 bp at the head and tail
  l <- Biostrings::nchar(genome)
  ir_bounder <- irBounder(genome, shifter = 0)

  if(is.data.frame(ir_bounder)) {
    return(ir_bounder)
  }

  # In case origin cut IRB
  if (ir_bounder$ir_bounder['irb_e'] == l) {
    # shift 1000 bp from head to tail and retest ir bounders
    ir_bounder <- irBounder(genome, shifter = 1000)

  }
  # In case origin cut IRA
  if (ir_bounder$ir_bounder['ira_s'] == 0) {
    # shift 1000 bp from head to tail and retest ir bounders
    ir_bounder <- irBounder(genome, shifter = -1000)
  }

  # Deal with IR region cutted by origins
  ira_s2 <- NA
  irb_e2 <- NA
  ira_s <- ir_bounder$ir_bounder['ira_s']
  ira_e <- ir_bounder$ir_bounder['ira_e']
  irb_s <- ir_bounder$ir_bounder['irb_s']
  irb_e <- ir_bounder$ir_bounder['irb_e']
  indel_table <- ir_bounder$indel_table
  ira_len <- ira_e - ira_s
  irb_len <- irb_e - irb_s

  if (ira_s < 0) {
    ira_s2 <- l + ira_s
    ira_s <- 0
  }
  if (irb_e > l) {
    irb_e2 <- irb_e - l
    irb_e <- l
  }

  # IR information table
  lsc_len <- ira_s + l - irb_e
  ssc_len <- irb_s - ira_e
  ir_table <- data.frame(chr = rep("chr1", 5),
                         start = c(0, ira_s, ira_e, irb_s, irb_e),
                         end = c(ira_s, ira_e, irb_s, irb_e, l),
                         name = c("LSC", "IRA", "SSC", "IRB", "LSC"),
                         text = c(paste("LSC:", lsc_len),
                                  paste("IRA:", ira_len),
                                  paste("SSC:", ssc_len),
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
    lsc_len <- ira_s - irb_e2
    ir_table$text[2] <- paste("LSC:", lsc_len)
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
    lsc_len <- ira_s2 - irb_e
    ir_table$text[4] <- paste("LSC:", lsc_len)
  }
  ir_table <- gc_count_ir(genome, ir_table)

  if (lsc_len < ssc_len){
    ir_table$name <- sub("LSC", "tmp", ir_table$name, fixed = TRUE)
    ir_table$name <- sub("SSC", "LSC", ir_table$name, fixed = TRUE)
    ir_table$name <- sub("tmp", "SSC", ir_table$name, fixed = TRUE)
    ir_table$text <- sub("LSC", "tmp", ir_table$text, fixed = TRUE)
    ir_table$text <- sub("SSC", "LSC", ir_table$text, fixed = TRUE)
    ir_table$text <- sub("tmp", "SSC", ir_table$text, fixed = TRUE)
  }
  if (is.null(indel_table)){
    return(ir_table)
  } else {
    return(list(ir_table = ir_table, indel_table = indel_table))
  }

}
