library(circlize)
library(dplyr)

gbFile <- "Nicotiana tabacum copy.gb"
gb<- read.gb(gbFile)
genome<- rdnFixer(gb)
Gene_table <-  gb.gene.cor(gb)
genome <- data.frame(name = c("genome"),
                     start = 1,
                     end = L)

gene_c <- Gene_table[grepl(".*complement.*", Gene_table[, 2]), ]


gb<- read.gb(gbFile)
genome<- rdnFixer(gb)
L<- length(genome)

# 1. IR LSC SSC -----------------------------------------------------------
IR<- IRinfo(genome)
ils <- data.frame(chr = rep("chr1", 5),
                  start = c(0, IR[1], IR[1] + IR[3] - 1, IR[2], 
                            IR[2] + IR[3] - 1),
                  end = c(IR[1], IR[1] + IR[3] - 1, IR[2], IR[2] + IR[3] - 1, L),
                  name = c("LSC", "IRA", "SSC", "IRB", ""),
                  y_t = c(0.25, 0.5, 0.25, 0.5, 0.25),
                  y_b = (0 - c(0.25, 0.5, 0.25, 0.5, 0.25))) %>% 
  mutate(center = (start + end)/2)

# 2. Gene table -----------------------------------------------------------
gene.tab<- gene.cordinates(gb)
gene_table <- data.frame(chr = rep("chr1", nrow(gene.tab)),
                         start = as.numeric(gene.tab[, 3]),
                         end = as.numeric(gene.tab[, 2]),
                         gene = gene.tab[,1]) %>%
  rbind.data.frame(data.frame(chr = rep("chr1", nrow(gene.tab)),
                        start = as.numeric(gene.tab[, 5]),
                        end = as.numeric(gene.tab[, 4]),
                        gene = gene.tab[,1])) %>%
  filter(end != 0) %>%
  unique() %>% 
  mutate(col = as.numeric(gene))
  
gene_table$chr <- as.character(gene_table$chr)
gene_table_f <- filter(gene_table, start < end) %>% 
  arrange(start)
gene_table_r <- filter(gene_table, start > end) %>%
  rename(start = "end", end = "start") %>%
  select(chr, start, end, gene, col) %>% 
  arrange(start)

# 4. Draw the plot --------------------------------------------------------

# Initialize the layout
jpeg("Circ2.jpg",  width=8.3, height=8.9, units="in", res=1000)
circos.clear()
circos.par(start.degree = 190, gap.after = 0, track.margin = c(0, 0), 
           cell.padding = c(0, 0, 0, 0))
circos.genomicInitialize(data = data.frame(chr = "chr1", start = 0, end = L,
                                           stringsAsFactors = FALSE),
                         plotType = NULL)


# 1.2. gene label outside

circos.genomicLabels(gene_table_r, labels.column = 4, side = "outside", 
                     cex = 0.5, connection_height = convert_height(3, "mm"),
                     labels_height = convert_height(3, "mm"))

# 3. Gene rectangles
circos.genomicTrack(gene_table_f, factors = as.factor("chr1"), ylim = c(-5, 5), 
                    track.height = convert_height(10, "mm"), 
                    bg.border = NA, 
                    panel.fun = function(region, value, ...) {
                      # genes lay on forward chain
                       circos.genomicRect(gene_table_f[, c("start", "end")], 
                                          value = gene_table_f$gene,
                                          ybottom = -5, 
                                          ytop = 0, 
                                          col = gene_table_f$col)
                      # genes lay on reverse chain
                       circos.genomicRect(gene_table_r[, c("start", "end")], 
                                          value = gene_table_r$gene,
                                          ybottom = 0, 
                                          ytop = 5,
                                          col = gene_table_r$col)
                       # the bars in the middle indicate the IR, SSR, LSR
                       circos.rect(xleft = ils$start, ybottom = ils$y_b, 
                                   xright = ils$end, ytop = ils$y_t,
                                   col = c("black", "orange", "black", 
                                           "orange", "black"),
                                   border = NA)
             })

# 4.5. gene labels inside

circos.genomicLabels(gene_table_f, labels.column = 4, side = "inside", cex = 0.5,
                     connection_height = convert_height(3, "mm"), 
                     labels_height = convert_height(3, "mm"))



# 6. Genome axis
circos.track(ylim = c(-1, 1), track.height = 0.15, bg.border = NA, 
             panel.fun = function(x, y) {
               # circos.arrow(x1 = 0, x2 = L %/% 10000, y = -1.5,
               #              arrow.position = "start", 
               #              col = "grey", border = NA, )
               circos.arrow(x1 = 0, x2 = L %/% 40, y = -0.3,
                            col = "grey", border = NA, 
                            arrow.head.length = convert_x(3, "mm"),
                            width = convert_y(1, "mm"))
               circos.genomicAxis(h = "bottom")
               })

# 7. arrow inside
circos.track(ylim = c(-1, 1), track.height = 0.1, bg.border = NA,
             panel.fun = function(x, y){
               circos.arrow(x1 = L - L %/% 35, x2 = L, y = 0.3,
                            col = "grey", border = NA, 
                            arrow.position = "start", 
                            arrow.head.length = convert_x(3, "mm"),
                            width = convert_y(1, "mm"))
             })

# 8. inner ring for IR, SSR, LSR
circos.genomicTrack(ils, bg.border = NA, 
                    ylim = c(-1, 10), 
                    panel.fun = function(region, value, ...) {
                      circos.rect(xleft = ils$start, 
                                  ybottom = convert_y(-1, "mm"), 
                                  xright = ils$end, ytop = convert_y(0.5, "mm"),
                                  col = c("skyblue", "blue", "green", 
                                          "blue", "skyblue"),
                                  border = NA)
                      circos.text(x = ils$center, y = convert_y(1.5, "mm"),
                                  labels = ils$name, cex = 0.5,
                                  facing = "inside",
                                  niceFacing = TRUE)
                    })

text(0,0.05, sp.name(gb), font=4, cex = 1)
text(0,-0.05, paste(prettyNum(L, big.mark = ","), "bp", " "),font=4, cex = 1)

dev.off()
