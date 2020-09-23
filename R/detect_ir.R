irDetect <- function(genome, seed.size = 1000) {
  tick = 0

  # Matching seeds to genome ------------------------------------------------
  # Get the reverse complement version of genome
  genome_rc <- Biostrings::reverseComplement(genome)

  # Get the length of genome
  l <- Biostrings::nchar(genome_rc)

  # set seeds start points
  seed_starts <- seq(1, (l - seed.size + 1))

  other_letter <- Biostrings::letterFrequency(genome, letters = "ATCG") != l

  m <- map_genome(genome, seed_starts = seed_starts, seed.size = seed.size,
                  other_letter = other_letter)
  if (nrow(m) == 0){ # if there is no matches return empty table
    ir_table <- data.frame(chr = "chr1",
                           start = 0,
                           end = l,
                           center = round(l/2),
                           name = "Genome",
                           text = paste("Genome:", l),
                           stringsAsFactors = FALSE)
    ir_table <- gc_count_ir(genome, ir_table)
    return(list(ir_table = ir_table, indel_table = NULL))
  }

  # if IRA cover genome start point, shift the genome sequence backward with seed.size
  while(m$group[1] == 1){
    tick <- tick + seed.size
    genome <- c(genome[(l - tick + 1):l], genome[1:(l - tick)])
    m <- map_genome(genome, seed_starts = seed_starts, seed.size = seed.size,
                    other_letter = other_letter)
  }

  # if IRB cover genome end point, shift the genome sequence forward with seed.size
  while(m$group[nrow(m)] + seed.size - 1 == l){
    tick <- tick - seed.size
    genome <- c(genome[(- tick + 1):l], genome[1: -tick])
    m <- map_genome(genome, seed_starts = seed_starts, seed.size = seed.size,
                    other_letter = other_letter)
  }

  m <- m[!duplicated(m$group),]
  # IRA start, end and lenght
  df <- data.frame(group_before = m$group[-nrow(m)],
                   group_after = m$group[-1],
                   group_diff = m$group[-1]- m$group[-nrow(m)] - 1,
                   start_diff = m$start[-1]- m$start[-nrow(m)] - 1,
                   start_before = m$start[-nrow(m)],
                   start_after = m$start[-1],
                   end_diff = m$end[-1]- m$end[-nrow(m)] - 1,
                   end_before = m$end[-nrow(m)],
                   end_after = m$end[-1],
                   stringsAsFactors = FALSE)

  # extract the first and last rows and the rows where group diff != 0
  pos <- df[c(1, nrow(df)), ] %>%
    rbind.data.frame(dplyr::filter(df, group_diff != 0)) %>%
    dplyr::arrange(group_before)

  # get start and end points for IRA and IRB (1-base)
  if (sum(pos$start_diff < 0) > 0){
    tmp <- pos[pos$start_diff < 0, ]
    pos <- rbind.data.frame(pos, df[which(df$group_diff %in% tmp$group_diff) - 1, ],
                            df[which(df$group_diff %in% tmp$group_diff) + 1, ]) %>%
      dplyr::filter(start_diff >= 0) %>%
      arrange(group_before)

    mismatch_group <- which(pos$group_diff > 0)
    ira_s <- pos$group_before[mismatch_group[1] - 1]
    ira_e <- pos$group_before[which.max(pos$group_diff)] + seed.size - 1
    irb_s <- l - pos$start_before[which.max(pos$group_diff)] - seed.size + 2
    irb_e <- pos$group_after[mismatch_group[length(mismatch_group)] + 1] + seed.size -1
  } else {
    ira_s <- pos$group_before[1]
    ira_e <- pos$group_before[(nrow(pos)+1)/2] + seed.size - 1
    irb_s <- pos$group_after[(nrow(pos)+1)/2]
    irb_e <- pos$group_after[nrow(pos)] + seed.size - 1
  }

  # Calculate lengths of each regions (1-base)
  ira_len <- ira_e - ira_s + 1
  irb_len <- irb_e - irb_s + 1
  lsc_len <- ira_s - 1 + l - irb_e
  ssc_len <- irb_s - ira_e - 1

  # Detecting indels and replaces in IR-----------------------------------------

  if (nrow(pos) > 3){
    ira_seq <- genome[ira_s:ira_e]
    irb_seq <- genome[irb_s:irb_e]
    other_letter <- (Biostrings::letterFrequency(ira_seq, letters = "ATCG") !=
                       Biostrings::nchar(ira_seq)) |
      (Biostrings::letterFrequency(irb_seq, letters = "ATCG") !=
                                            Biostrings::nchar(irb_seq))
    indel_table <- detect_mismatch(ira_seq = ira_seq, irb_seq = irb_seq,
                                   ira_s = ira_s, irb_s = irb_s,
                                   other_letter = other_letter)
  } else {
    indel_table <- NULL
  }

  # If the difference between IRA length and IRB length is not equal to the
  # number of the "insert" base pairs, the result should be wrong
  # ir_dff <- indel_table[indel_table$mismatch_type == c("insert", "delete"),]
  # ir_dff$ira <- ir_dff$position < ira_e & ir_dff$position > ira_s
  # ir_dff$dff <- (ir_dff$mismatch_type == "insert") == ir_dff$ira
  # ir_dff <- abs(sum(ifelse(ir_dff$dff, 1, -1)))
  if (ira_len != irb_len){
    if (is.null(indel_table) |
        sum(indel_table$mismatch_type %in% c("insert", "delete")) == 0 |
        abs(ira_len - irb_len) > 100){
      stop("The IR regions are not in similar length and no indel was detected")
    }
  }

  # recover original coordinates (0-base)
  if (tick == 0) {
    ir_table <- data.frame(chr = rep("chr1", 5),
                           start = c(0, ira_s - 1, ira_e, irb_s - 1, irb_e),
                           end = c(ira_s - 1, ira_e, irb_s - 1, irb_e, l),
                           name = c("LSC", "IRA", "SSC", "IRB", "LSC"),
                           text = c(paste("LSC:", lsc_len),
                                    paste("IRA:", ira_len),
                                    paste("SSC:", ssc_len),
                                    paste("IRB:", irb_len),
                                    ""),
                           stringsAsFactors = FALSE)
  } else if (tick > 0) {
    if ((tick - ira_s) == 0){
      ir_table <- data.frame(chr = rep("chr1", 4),
                             start = c(0, ira_e - tick, irb_s - tick,
                                       irb_e - tick),
                             end = c(ira_e - tick, irb_s - tick,
                                     irb_e - tick, l),
                             name = c("IRA", "SSC", "IRB", "LSC"),
                             text = c(paste("IRA:", ira_len),
                                      paste("SSC:", ssc_len),
                                      paste("IRB:", irb_len),
                                      paste("LSC:", lsc_len)),
                             stringsAsFactors = FALSE)
    } else {
      ir_table <- data.frame(chr = rep ("chr1", 5),
                             start = c(0, ira_e - tick, irb_s - tick - 1,
                                       irb_e - tick, l-(tick - ira_s) - 1),
                             end = c(ira_e - tick, irb_s - tick - 1,
                                     irb_e - tick, l-(tick - ira_s) - 1, l),
                             name = c("IRA", "SSC", "IRB", "LSC", "IRA"),
                             text = c(paste("IRA:", ira_len),
                                      paste("SSC:", ssc_len),
                                      paste("IRB:", irb_len),
                                      paste("LSC:", lsc_len),
                                      ""),
                             stringsAsFactors = FALSE)
    }
  } else if (tick < 0) {
    if ((irb_e-tick-l) == 0) {
      ir_table <- data.frame(chr = rep ("chr1", 4),
                             start = c(0, ira_s - tick - 1,
                                       ira_e - tick, irb_s - tick - 1),
                             end = c(ira_s - tick - 1, ira_e - tick,
                                     irb_s - tick - 1, l),
                             name = c("LSC", "IRA", "SSC", "IRB"),
                             text = c(paste("LSC:", lsc_len),
                                      paste("IRA:", irb_len),
                                      paste("SSC:", ssc_len),
                                      paste("IRB:", irb_len)),
                             stringsAsFactors = FALSE)
    } else {
      ir_table <- data.frame(chr = rep ("chr1", 5),
                             start = c(0, (irb_e - tick - l), ira_s - tick - 1,
                                       ira_e - tick, irb_s - tick - 1),
                             end = c((irb_e - tick - l), ira_s - tick - 1,
                                     ira_e - tick, irb_s - tick - 1, l),
                             name = c("IRB", "LSC", "IRA", "SSC", "IRB"),
                             text = c(paste("IRB:", irb_len),
                                      paste("LSC:", lsc_len),
                                      paste("IRA:", irb_len),
                                      paste("SSC:", ssc_len),
                                      ""),
                             stringsAsFactors = FALSE)
    }
  }

  # Add text coordinates
  ir_table$center <- round((ir_table$start + ir_table$end)/2, 0)
  if (nrow(ir_table) == 5){
    ir_table$center[1] <- ir_table$center[1] + round(l/2) + round(ir_table$start[5]/2)
  }

  if(ir_table$center[1] > l){
    ir_table$center[1] <- ir_table$center[1] - l
  }

  # Add gc count
  ir_table <- gc_count_ir(genome, ir_table)

  if (lsc_len < ssc_len){
    ir_table$name <- sub("LSC", "tmp", ir_table$name, fixed = TRUE)
    ir_table$name <- sub("SSC", "LSC", ir_table$name, fixed = TRUE)
    ir_table$name <- sub("tmp", "SSC", ir_table$name, fixed = TRUE)
    ir_table$text <- sub("LSC", "tmp", ir_table$text, fixed = TRUE)
    ir_table$text <- sub("SSC", "LSC", ir_table$text, fixed = TRUE)
    ir_table$text <- sub("tmp", "SSC", ir_table$text, fixed = TRUE)
  }

  return(list(ir_table = ir_table, indel_table = indel_table))

}

map_genome <- function(genome, seed_starts, seed.size = 1000,
                       other_letter = FALSE){

  l <- Biostrings::nchar(genome)
  genome_rc <- Biostrings::reverseComplement(genome)
  # Breack the whole genome into small pieces (seeds) with size [seed.size]bp

  seeds <- Biostrings::DNAStringSet(genome, start = seed_starts, width = seed.size,)
  if (other_letter){
    names(seeds) <- 1:length(seeds)
    seeds <- seeds[grepl("^[ATCG]", seeds)] # remove the reads start with letters rather than "ATCG"
    seeds <- Biostrings::PDict(seeds, tb.start = 1, tb.end = 1)
    m <- Biostrings::matchPDict(seeds, genome_rc, fixed = FALSE)
    deleted_group <- setdiff(1:(length(seeds) + 1), names(m))
    m <- as.data.frame(m)
    for (i in 1:length(deleted_group)){
      m$group[which(m$group > deleted_group[i])] <- m$group[which(m$group > deleted_group[i])] + 1
    }
  } else{
    seeds <- Biostrings::PDict(seeds)
    m <- Biostrings::matchPDict(seeds, genome_rc, max.mismatch = )
    m <- as.data.frame(m)
  }

  # Mapping seeds to the reverse conplemented genome
  m <- m[m$width > (seed.size * 0.9), ] # Keeping the seads with 90% matching to the genome

  return(m)
}

detect_mismatch <- function(ira_seq, irb_seq, ira_s, irb_s, other_letter){

  if (other_letter){
    ir_map_a <- Biostrings::pairwiseAlignment(pattern = ira_seq,
                                              subject = Biostrings::reverseComplement(irb_seq),
                                              fuzzyMatrix = Biostrings::nucleotideSubstitutionMatrix(),
                                              substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix())

    ir_map_b <- Biostrings::pairwiseAlignment(pattern = irb_seq,
                                              subject = Biostrings::reverseComplement(ira_seq),
                                              fuzzyMatrix = Biostrings::nucleotideSubstitutionMatrix(),
                                              substitutionMatrix = Biostrings::nucleotideSubstitutionMatrix())
  } else {
    ir_map_a <- Biostrings::pairwiseAlignment(pattern = ira_seq,
                                              subject = Biostrings::reverseComplement(irb_seq))

    ir_map_b <- Biostrings::pairwiseAlignment(pattern = irb_seq,
                                              subject = Biostrings::reverseComplement(ira_seq))
  }

  ## insert table for IRA and delete table for IRB

  insert_table_a <- suppressWarnings(data.frame(Biostrings::insertion(ir_map_a)))

  if (nrow(insert_table_a) != 0){

    insert_table_a$string <- as.character(Biostrings::DNAStringSet(ir_map_a@pattern,
                      start = insert_table_a$start, end = insert_table_a$end))
    n_skip <- Biostrings::letterFrequency(Biostrings::DNAStringSet(ir_map_a@pattern,
                      start = rep(1, nrow(insert_table_a)), end = insert_table_a$start),
                      letters = "-")
    insert_table_a <- insert_table_a %>%
      dplyr::mutate(start = start - n_skip + ira_s - 1) %>%
      dplyr::mutate(end = end - n_skip + ira_s - 1)

    insert_a <- NULL
    for (i in 1:nrow(insert_table_a)){
      tmp <- data.frame(mismatch_type = rep("insert", insert_table_a$width[i]),
                        position = seq(insert_table_a$start[i],
                                       insert_table_a$end[i]),
                        string = strsplit(insert_table_a$string[i], "")[[1]],
                        col = rep("green", insert_table_a$width[i]),
                        stringsAsFactors = FALSE)
      insert_a <- rbind.data.frame(insert_a, tmp)
    }

    delete_table_b <- suppressWarnings(data.frame(Biostrings::deletion(ir_map_b))) %>%
      dplyr::mutate(position = start + irb_s - 1) %>%
      dplyr::mutate(string = rep("D", n())) %>%
      dplyr::mutate(mismatch_type = rep("delete", n())) %>%
      dplyr::mutate(col = rep("yellow", n())) %>%
      dplyr::select(mismatch_type, position, string, col)
  } else {
    insert_a <- NULL
    delete_table_b <- NULL
  }

  ## insert table for IRB and delete table for IRA
  insert_table_b <- suppressWarnings(data.frame(Biostrings::insertion(ir_map_b)))

  if (nrow(insert_table_b) != 0) {

    insert_table_b$string <- as.character(Biostrings::DNAStringSet(ir_map_b@pattern,
                      start = insert_table_b$start, end = insert_table_b$end))
    n_skip <- Biostrings::letterFrequency(Biostrings::DNAStringSet(ir_map_b@pattern,
                              start = rep(1, nrow(insert_table_b)),
                              end = insert_table_b$start),
                              letters = "-")
    insert_table_b <- insert_table_b %>%
      dplyr::mutate(start = start - n_skip + irb_s - 1) %>%
      dplyr::mutate(end = end - n_skip + irb_s - 1)

    insert_b <- NULL
    for (i in 1:nrow(insert_table_b)){
      tmp <- data.frame(mismatch_type = rep("insert", insert_table_b$width[i]),
                        position = seq(insert_table_b$start[i],
                                       insert_table_b$end[i]),
                        string = strsplit(insert_table_b$string[i], "")[[1]],
                        col = rep("green", insert_table_b$width[i]),
                        stringsAsFactors = FALSE)
      insert_b <- rbind.data.frame(insert_b, tmp)
    }


    delete_table_a <- suppressWarnings(data.frame(Biostrings::deletion(ir_map_a))) %>%
      dplyr::mutate(position = start + ira_s - 1) %>%
      dplyr::mutate(string = rep("D", n())) %>%
      dplyr::mutate(mismatch_type = rep("delete", n())) %>%
      dplyr::mutate(col = rep("yellow", n())) %>%
      dplyr::select(mismatch_type, position, string, col)
  } else{
    insert_b <- NULL
    delete_table_a <- NULL
  }

  ## replace table
  replace_table_a <- Biostrings::mismatchTable(ir_map_a)

  if (nrow(replace_table_a) != 0){
    tmp <- NULL
    for (i in 1:nrow(replace_table_a)){
      if (Biostrings::nucleotideSubstitutionMatrix()[as.character(replace_table_a$PatternSubstring[i]),
                as.character(replace_table_a$SubjectSubstring[i])] != 0){
        tmp <- c(tmp, i)
      }
    }

    if (!is.null(tmp)){
      replace_table_a <- replace_table_a[!(1:nrow(replace_table_a) %in% tmp), ]
    }
  }

  if (nrow(replace_table_a) != 0) {
    n_skip <- as.vector(Biostrings::letterFrequency(Biostrings::DNAStringSet(ir_map_a@pattern,
                                           start = rep(1, nrow(replace_table_a)),
                                           end = replace_table_a$PatternStart),
                              letters = "-"))
    replace_table_ira <- data.frame(position = replace_table_a$PatternStart - n_skip + ira_s - 1,
                                    string = replace_table_a$PatternSubstring,
                                    mismatch_type = rep("replace",
                                                        nrow(replace_table_a)),
                                    col = rep("red", nrow(replace_table_a)),
                                    stringsAsFactors = FALSE)

    replace_table_b <- Biostrings::mismatchTable(ir_map_b)


    tmp <- NULL
    for (i in 1:nrow(replace_table_b)){
      if (Biostrings::nucleotideSubstitutionMatrix()[as.character(replace_table_b$PatternSubstring[i]),
                as.character(replace_table_b$SubjectSubstring[i])] != 0){
        tmp <- c(tmp, i)
      }
    }

    if (!is.null(tmp)){
      replace_table_b <- replace_table_b[!(1:nrow(replace_table_b) %in% tmp), ]
    }

    n_skip <- as.vector(Biostrings::letterFrequency(Biostrings::DNAStringSet(ir_map_b@pattern,
                                                       start = rep(1, nrow(replace_table_b)),
                                                       end = replace_table_b$PatternStart),
                              letters = "-"))

    replace_table_irb <- data.frame(position = replace_table_b$PatternStart - n_skip + irb_s - 1,
                                    string = replace_table_b$PatternSubstring,
                                    mismatch_type = rep("replace",
                                                        nrow(replace_table_b)),
                                    col = rep("red", nrow(replace_table_b)),
                                    stringsAsFactors = FALSE)


  } else {
    replace_table_ira <- NULL
    replace_table_irb <- NULL
  }
  indel_table <- Reduce(rbind.data.frame,
                           list(insert_a, insert_b, delete_table_a,
                                delete_table_b, replace_table_ira,
                                replace_table_irb))
  if (nrow(indel_table) == 0){
    indel_table <- NULL
  }
  return(indel_table)
}


