geneColor <- function(x, y){
  color.table <- data.frame(
    acronym = c("psa","psb","pet","atp","ndh","rbc","rpo","rps","rpl",
                "clp|mat|inf","ycf","trn","rrn", "OTHER"),
    col = c(grDevices::rgb(42, 99, 50, maxColorValue = 255), # {psa}
            grDevices::rgb(76, 136, 5, maxColorValue = 255),# {psb}
            grDevices::rgb(127, 153, 44, maxColorValue = 255),# {pet}
            grDevices::rgb(159, 187, 61, maxColorValue = 255),# {atp}
            grDevices::rgb(254, 238, 80, maxColorValue = 255),# {ndh}
            grDevices::rgb(77, 158, 63, maxColorValue = 255),# {rbc}
            grDevices::rgb(174, 45, 41, maxColorValue = 255),# {rpo}
            grDevices::rgb(214, 173, 124, maxColorValue = 255),# {rps}
            grDevices::rgb(156, 122, 75, maxColorValue = 255),# {rpl}
            grDevices::rgb(217, 102, 45, maxColorValue = 255),# {clp, mat, inf}
            grDevices::rgb(113, 184, 169, maxColorValue = 255),# {ycf}
            grDevices::rgb(23, 44, 127, maxColorValue = 255),# {trn}
            grDevices::rgb(209, 56, 42, maxColorValue = 255),# {rrn}
            grDevices::rgb(125, 125, 125, maxColorValue = 255)# {OTHER}
          ),
    label = c("photosystem I", "photosystem II", "cytochrome b/f complex",
              "ATP synthesis", "NADH dehydrogenase", "RubisCO larg subunit",
              "RNA polymerase", "small ribosomal protein",
              "large ribosomal protein", "clpP, matK, infA",
              "hypothetical reading frame", "transfer RNA", "ribosomal RNA",
              "other"),
    stringsAsFactors = FALSE
  )
}
