geneColor <- function(x, y, transparent = 180){
  color.table <- data.frame(
    acronym = c("psa","psb","pet","atp","ndh","rbc","rpo","rps","rpl",
                "clp|mat|inf","ycf","trn","rrn", "OTHER"),
    col = c(rgb(42, 99, 50, alpha = transparent, maxColorValue = 255), # {psa}
            rgb(76, 136, 5, alpha = transparent, maxColorValue = 255),# {psb}
            rgb(127, 153, 44, alpha = transparent, maxColorValue = 255),# {pet}
            rgb(159, 187, 61, alpha = transparent, maxColorValue = 255),# {atp}
            rgb(254, 238, 80, alpha = transparent, maxColorValue = 255),# {ndh}
            rgb(77, 158, 63, alpha = transparent, maxColorValue = 255),# {rbc}
            rgb(174, 45, 41, alpha = transparent, maxColorValue = 255),# {rpo}
            rgb(214, 173, 124, alpha = transparent, maxColorValue = 255),# {rps}
            rgb(156, 122, 75, alpha = transparent, maxColorValue = 255),# {rpl}
            rgb(217, 102, 45, alpha = transparent, maxColorValue = 255),# {clp, mat, inf}
            rgb(113, 184, 169, alpha = transparent, maxColorValue = 255),# {ycf}
            rgb(23, 44, 127, alpha = transparent, maxColorValue = 255),# {trn}
            rgb(209, 56, 42, alpha = transparent, maxColorValue = 255),# {rrn}
            rgb(125, 125, 125, alpha = transparent, maxColorValue = 255)# {OTHER}
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
