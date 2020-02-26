#' Generate table for PlotGenome function
#'
#' @param gbfile A character. It is the GI or GenBank accession for the
#' interested specie's chloroplast genome files.
#' @param local.file A logical value. If it is \code{TRUE}, the value passed to
#' \code{gbfile} should be a local GB file.
#' @param gc.window A integer. It indicates the window sized for plot
#' GC count.
#' @return A list. It contains 3 tables:
#' \itemize{
#'   \item ir_table: a data frame contains information for IR region.
#'   \item gene_table: a data frame contains information of genes.
#'   \item gc_count: a data frame contains GC count information.
#' }
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
PlotTab <- function(gbfile, local.file = FALSE, gc.window = 100){
  if (local.file){
    tryCatch({
      gb <- genbankr::readGenBank(gbfile)
      genome <- genbankr::getSeq(gb)
      genome <- rdnFixer(genome)
      L<- Biostrings::nchar(genome)
      sp_name <- sp.name(gb@definition)
      genes <- as.data.frame(genbankr::genes(gb))

      for (i in 1:nrow(genes)){
        if (is.na(genes$gene[i])){
          genes$gene[i] <- genes$gene_id[i]
        }
      }

      features <- as.data.frame(genbankr::otherFeatures(gb))

      features <- features %>%
        dplyr::filter(!gene %in% genes$gene)

      for (i in 1:nrow(features)){
        if (is.na(features$gene[i])){
          features$gene[i] <- features$product[i]
        }
      }
      gene_table <- genes %>%
        select(start, end, gene, strand) %>%
        rbind.data.frame(select(features, start, end, gene, strand)) %>%
        mutate(chr = rep("chr1", n()))

    }, error = function(e){
      gb <<- readLines(gbfile)
      genome <<- FasExtract(gb)
      L <<- Biostrings::nchar(genome)
      sp_name <<- sp.name(gb, text = TRUE)
      gene_table <<- data.frame(chr = rep("chr1", nrow(gene.tab)),
                               start = as.numeric(gene.tab[, 3]),
                               end = as.numeric(gene.tab[, 2]),
                               gene = gene.tab[,1],
                               stringsAsFactors = FALSE) %>%
        rbind.data.frame(data.frame(chr = rep("chr1", nrow(gene.tab)),
                                    start = as.numeric(gene.tab[, 5]),
                                    end = as.numeric(gene.tab[, 4]),
                                    gene = gene.tab[,1])) %>%
        filter(end != 0) %>%
        unique()
    })

  } else {
    tryCatch({
      gb <- genbankr::readGenBank(text = fetch.gb(gbfile))
      genome <- genbankr::getSeq(gb)
      genome <- rdnFixer(genome)
      L<- Biostrings::nchar(genome)
      sp_name <- sp.name(gb@definition)
      genes <- as.data.frame(genbankr::genes(gb))

      for (i in 1:nrow(genes)){
        if (is.na(genes$gene[i])){
          genes$gene[i] <- genes$gene_id[i]
        }
      }

      features <- as.data.frame(genbankr::otherFeatures(gb))

      features <- features %>%
        dplyr::filter(!gene %in% genes$gene)

      for (i in 1:nrow(features)){
        if (is.na(features$gene[i])){
          features$gene[i] <- features$product[i]
        }
      }
      gene_table <- genes %>%
        select(start, end, gene, strand) %>%
        rbind.data.frame(select(features, start, end, gene, strand)) %>%
        mutate(chr = rep("chr1", n()))

    }, error = function(e){
      gb <<- fetch.gb(gbfile)
      genome <<- FasExtract(gb)
      L <<- Biostrings::nchar(genome)
      sp_name <<- sp.name(gb, text = TRUE)
      gene_table <<- data.frame(chr = rep("chr1", nrow(gene.tab)),
                                start = as.numeric(gene.tab[, 3]),
                                end = as.numeric(gene.tab[, 2]),
                                gene = gene.tab[,1],
                                stringsAsFactors = FALSE) %>%
        rbind.data.frame(data.frame(chr = rep("chr1", nrow(gene.tab)),
                                    start = as.numeric(gene.tab[, 5]),
                                    end = as.numeric(gene.tab[, 4]),
                                    gene = gene.tab[,1])) %>%
        filter(end != 0) %>%
        unique()
    })
  }


  # 1. IR LSC SSC -----------------------------------------------------------
  ils <- irDetect(genome)

  # 2. Gene table -----------------------------------------------------------

  #colouring
  color_table <- geneColor(-10, -10)
  gene_table$col <- rep("grey", nrow(gene_table))
  for (i in 1:nrow(color_table)){
    gene_table$col[which(grepl(as.character(color_table$acronym[i]),
                               gene_table$gene, perl = TRUE))] <-
      color_table$col[i]
  }

  # 4. GC count -------------------------------------------------------------

  gc_count_list <- gc_count(genome, view.width = gc.window)
  gc_count <- gc_count_list[[1]]
  gc_total <- gc_count_list[[2]]
  gc_count$chr <- rep("chr1", nrow(gc_count))
  gc_count <- select(gc_count, chr, position, gc_count)

  tables <- list(ir_table = ils, gc_count = gc_count, gc_total = gc_total,
                 sp_name = sp_name, genome_len = L, gene_table = gene_table,
                 gene_color = color_table)
  return(tables)
}

#' Generate Genome plot
#'
#' @param plot_tables a list contains information of IR region, gene, and gc
#' count of the genome. It can be generated by function \code{\link{PlotTab}}.
#' @param save A logical value. If it is \code{TRUE}, the plot will be saved in
#' work directory.
#' @param file.type A charactor. It indicates the format of the file in which
#' the plot will be saved. Options are: "pdf", "png", "jpeg","bmp", "tiff".
#' By default, it is set as "pdf".
#' @param file.name A charactor. It indicates the name of the file in which
#' the plot will be saved. By default, it is set as the specie's name.
#' @param shadow A logical value. If it is \code{TRUE}, the section of IR
#' regions will be highlighted by shadows.
#' @param legend A logical value. If it is \code{TRUE}, the legend for gene
#' colors will be shown.
#' @param ssc.converse A logical value. If it is \code{TRUE}, the SSC region
#' @param background A character. It indicates the color for the background of
#' entire plot area.
#' @param gc.color A character. It indicates the color for the lines in gc count
#' plot.
#' @param gc.background A character. It indicates the color for the background
#' of gc count plot area.
#' @param info.background A character. It indicates the color for the background
#' of central area where species' information was shown.
#' @param ir.color A character. It indicates the color for the background
#' of IR sectors.
#' @param shadow.color A character. It indicates the color for the shadow
#' casted from IR sectors.
#' @param ssc.color A character. It indicates the color for the background
#' of SSC sectors.
#' @param lsc.color A character. It indicates the color for the background
#' of LSC sectors.
#' will be converted to its reverse complementary version
#' @return A plot for chloroplast genome.
#' @importFrom magrittr %>%
#' @import dplyr
#' @export
#'
PlotGenome <- function(plot_tables, save = TRUE, file.type = "pdf",
                       file.name = NULL, shadow = TRUE, legend = TRUE,
                       ssc.converse = FALSE,
                       background = "grey90", gc.color = "grey30",
                       gc.background = "grey70", info.background = "black",
                       ir.color = "#2F3941", shadow.color = "#0000FF20",
                       ssc.color = "#82B6E2", lsc.color = "#299E96"){

  # Modify gene table
  if (ssc.converse){
    if (nrow(plot_tables$ir_table) < 4){
      warning("Didn't get IR region from '", plot_tables$sp_name,
              "genome'.It's impossible to convert SSC region.")
      gene_table <- plot_tables$gene_table
    } else {
      gene_table <- SSCrev(plot_tables$gene_table,
                           SSCs = plot_tables$ir_table$start[plot_tables$ir_table$name == "SSC"],
                           SSCe = plot_tables$ir_table$end[plot_tables$ir_table$name == "SSC"])
    }

  } else {
    gene_table <- plot_tables$gene_table
  }

  if ("strand" %in% colnames(gene_table)){
    gene_table_f <- filter(gene_table, strand == "+") %>%
      select(chr, start, end, gene, col) %>%
      arrange(start)
    gene_table_r <- filter(gene_table, strand == "-") %>%
      select(chr, start, end, gene, col) %>%
      arrange(start)
  } else {
    gene_table_f <- filter(gene_table, start < end) %>%
      select(chr, start, end, gene, col) %>%
      arrange(start)
    gene_table_r <- filter(gene_table, start > end) %>%
      rename(start = "end", end = "start") %>%
      select(chr, start, end, gene, col) %>%
      arrange(start)
  }


  # Set colors
  plot_tables$ir_table$bg_col <- rep(NA, nrow(plot_tables$ir_table))
  plot_tables$ir_table$bg_col[plot_tables$ir_table$name == "LSC"] <- lsc.color
  plot_tables$ir_table$bg_col[plot_tables$ir_table$name == "SSC"] <- ssc.color
  plot_tables$ir_table$bg_col[grepl("IR", plot_tables$ir_table$name)] <- ir.color

  # Automatically adjust colors
  info.color <- CompColor(info.background)

  plot_tables$ir_table$inf_col <- rep(NA, nrow(plot_tables$ir_table))
  plot_tables$ir_table$inf_col[plot_tables$ir_table$name=="LSC"] <- CompColor(lsc.color)
  plot_tables$ir_table$inf_col[plot_tables$ir_table$name=="SSC"] <- CompColor(ssc.color)
  plot_tables$ir_table$inf_col[grepl("IR", plot_tables$ir_table$name)] <- CompColor(ir.color)

  L <- plot_tables$genome_len

  # Initialize the plot device
  if (save){
    if (is.null(file.name)){
      file <- plot_tables$sp_name
    } else {
      file <- file.name
    }
    if (!file.type %in% c("jpeg", "bmp", "png", "tiff", "pdf")){
      warning("Can not save plot in ", file.type, " format. Avaliable formats
                are 'jpeg', 'bmp', 'png', 'tiff',and 'pdf'.")
    } else if (file.type  == "pdf") {
      grDevices::pdf(paste(file, file.type, sep="."), width=10.3, height=8.9)
    } else {
      do.call(file.type, args=list(filename=paste(file, file.type, sep="."),
                                   width = 10.3, height = 8.9, units = "in",
                                   res = 1000))
    }
  }

  # Initialize the layout
  circlize::circos.clear()
  circlize::circos.par(start.degree = 180,
                       gap.after = 0,
                       track.margin = c(0, 0), cell.padding = c(0, 0, 0, 0))
  circlize::circos.genomicInitialize(data=data.frame(chr="chr1", start=0, end=L,
                                                     stringsAsFactors = FALSE),
                                     plotType = NULL)


  # 1.2. gene label outside

  circlize::circos.genomicLabels(gene_table_r, labels.column = 4,
                                 side = "outside",cex = 0.5,
                                 connection_height = circlize::convert_height(3, "mm"),
                                 labels_height = circlize::convert_height(3, "mm"))

  circlize::draw.sector(0, 360,
                        # rou1 = circlize::get.cell.meta.data("cell.top.radius",
                        #                                     track.index = 3),
                        rou1 = circlize::get.cell.meta.data("cell.bottom.radius",
                                                            track.index = 2),
                        rou2 = 0,
                        col = "grey90", border = NA)
  # circlize::circos.genomicTrackPlotRegion(gene_table_r, ylim = c(0, 1),
  #                                         #bg.col = bg.col,
  #                               track.margin = c(convert_height(0.5, "mm"),
  #                                                circos.par("track.margin")[2]),
  #                               cell.padding = c(0, 0, 0, 0),
  #                               track.height = connection_height,
  #                               bg.border = NA)
  # 3. Gene rectangles
  circlize::circos.genomicTrack(gene_table_f, factors = as.factor("chr1"),
                                ylim = c(-5, 5), bg.border = NA,
                                track.height = circlize::convert_height(10, "mm"),
                                panel.fun = function(region, value, ...) {
                                  # genes lay on forward chain
                                  circlize::circos.genomicRect(
                                    gene_table_f[, c("start", "end")],
                                    value = gene_table_f$gene,
                                    ybottom = -5,
                                    ytop = 0,
                                    col = gene_table_f$col)
                                  # genes lay on reverse chain
                                  circlize::circos.genomicRect(
                                    gene_table_r[, c("start", "end")],
                                    value = gene_table_r$gene,
                                    ybottom = 0,
                                    ytop = 5,
                                    col = gene_table_r$col)
                                  # the bars in the middle indicate the IR, SSR, LSR
                                  circlize::circos.rect(xleft = plot_tables$ir_table$start,
                                                        ybottom = - 0.25,
                                                        xright = plot_tables$ir_table$end,
                                                        ytop = 0.25,
                                                        col = plot_tables$ir_table$bg_col,
                                                        border = NA)
                                })

  # 4.5. gene labels inside

  circlize::circos.genomicLabels(gene_table_f, labels.column = 4,
                                 side = "inside", cex = 0.5,
                                 connection_height = circlize::convert_height(3, "mm"),
                                 labels_height = circlize::convert_height(3, "mm"))


  # 6. Arrow outside gnome axis
  circlize::circos.track(ylim = c(0, 1), track.height = 0.1, bg.border = NA,
                         panel.fun = function(x, y) {
                           # circos.arrow(x1 = 0, x2 = L %/% 10000, y = -1.5,
                           #              arrow.position = "start",
                           #              col = "grey", border = NA, )
                           circlize::circos.arrow(x1 = 0,
                                                  x2 = L %/% 40,
                                                  y = 0.5,
                                                  col = MildColor(background), border = NA,
                                                  arrow.head.length = circlize::convert_x(3, "mm"),
                                                  width = circlize::convert_y(1, "mm"))
                           # circlize::circos.genomicAxis(h = "bottom")
                         })

  # 7. GC count
  circlize::circos.track(factors =as.factor(plot_tables$gc_count$chr),
                         x=plot_tables$gc_count$position,
                         y = plot_tables$gc_count$gc_count,
                         ylim = c(0, 1),
                         track.height = 0.2,
                         bg.border = NA, bg.col = gc.background,
                         panel.fun = function(x, y) {
                           circlize::circos.arrow(x1=L - L %/% 35, x2=L, y=0.9,
                                                  col=MildColor(gc.background), border=NA,
                                                  arrow.position="start",
                                                  arrow.head.length=circlize::convert_x(3, "mm"),
                                                  width=circlize::convert_y(1, "mm"))
                           # circlize::circos.lines(x, y, type = "h", area=TRUE,
                           #                        col = "grey", lwd = 0.3,
                           #                        pt.col=plot_tables$gc_count$peak_valley,
                           #                        cex = 0.2, pch = 16)
                           circlize::circos.lines(x, y, type = "h",
                                                  col = gc.color, lwd = 0.6)
                           circlize::circos.segments(x0=0, x1=L, y0=0.25,
                                                     y1=0.25, lwd=0.5,
                                                     lty="16", col="grey30")
                           circlize::circos.segments(x0=0, x1=L, y0=0.5,
                                                     y1=0.5, lwd=0.6,
                                                     lty="11", col="grey40")
                           circlize::circos.segments(x0=0, x1=L, y0=0.75,
                                                     y1=0.75, lwd=0.5,
                                                     lty="16", col="grey30")
                           circlize::circos.genomicAxis(h = "top")
                         })
  #circos.yaxis(at=c(0.25,0.5,0.75),labels.cex=0.25,lwd=0,tick.length=0,
  #             labels.col="grey40",col="#FFFFFF")

  # 8. inner ring for IR, SSR, LSR
  circlize::circos.genomicTrack(plot_tables$ir_table, bg.border = NA,
                                ylim = c(0, 3),track.height = 0.1,
                                panel.fun = function(region, value, ...) {
                                  circlize::circos.rect(xleft = plot_tables$ir_table$start,
                                                        ybottom = 0,#circlize::convert_y(0, "mm"),
                                                        xright = plot_tables$ir_table$end,
                                                        ytop = 3,#circlize::convert_y(3, "mm"),
                                                        col = plot_tables$ir_table$bg_col,
                                                        border = NA)
                                  circlize::circos.text(x = plot_tables$ir_table$center,
                                                        y = 1.5, cex = 0.7,
                                                        labels = plot_tables$ir_table$text,
                                                        facing = "bending.inside",
                                                        col = plot_tables$ir_table$inf_col,
                                                        niceFacing = TRUE)
                                })

  # Specie's information inthe central
  circlize::draw.sector(0, 360,
                        # rou1 = circlize::get.cell.meta.data("cell.top.radius",
                        #                                     track.index = 3),
                        rou1 = circlize::get.cell.meta.data("cell.bottom.radius",
                                                            track.index = 8),
                        rou2 = 0,
                        col = info.background, border = NA)
  text(0,0.05, plot_tables$sp_name, font=4, cex = 0.9, col = info.color)
  text(0,0, paste(prettyNum(L, big.mark = ","), "bp", " "),
       font=4, cex = 0.9, col = info.color)
  text(0,-0.05,paste("GC count:",round(plot_tables$gc_total, 4)*100,"%"),
       font=4, cex = 0.9, col = info.color)

  # Highlight IR

  if (shadow & (nrow(plot_tables$ir_table) >= 4)){
    pos = circlize::circlize(c(plot_tables$ir_table$start[2],
                               c(plot_tables$ir_table$end[2]),
                               c(plot_tables$ir_table$start[4]),
                               c(plot_tables$ir_table$end[4])),
                             c(0, 1), sector.index = "chr1", track.index = 1)
    circlize::draw.sector(pos[1, "theta"], pos[2, "theta"],
                          #start.degree = 0, end.degree = 360 - 10,
                          rou1 = circlize::get.cell.meta.data("cell.top.radius",
                                                              track.index = 3),
                          rou2 = circlize::get.cell.meta.data("cell.bottom.radius",
                                                              track.index = 7),
                          col = shadow.color, border = NA)
    circlize::draw.sector(pos[3, "theta"], pos[4, "theta"],
                          #start.degree = 0, end.degree = 360 - 10,
                          rou1 = circlize::get.cell.meta.data("cell.top.radius",
                                                              track.index = 3),
                          rou2 = circlize::get.cell.meta.data("cell.bottom.radius",
                                                              track.index = 7),
                          col = shadow.color, border = NA)
  }


  # Add legend
  if (legend){
    legend(x = -1.25, y = -0.6, legend = plot_tables$gene_color$label, cex = 0.6,
           fill = plot_tables$gene_color$col, box.col = "white")
  }

  if (save){
    dev.off()
  }

}
