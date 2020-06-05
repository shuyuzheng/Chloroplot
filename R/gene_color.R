plastidGeneColor <- function(psa.color = "#2A6332", psb.color = "#4C8805",
                      pet.color = "#7F992C", atp.color = "#9FBB3D",
                      ndh.color = "#FEEE50", rbc.color = "#4D9E3F",
                      rpo.color = "#AE2D29", rsp.color = "#D6AD7C",
                      rpl.color = "#9C7A4B", clp_mat_inf.color = "#D9662D",
                      ycf.color = "#71B8A9", trn.color = "#172C7F",
                      rrn.color = "#D1382A", other_gene.color = "#7D7D7D"){
  color.table <- data.frame(
    acronym = c("psa","psb","pet","atp","ndh","rbc","rpo","rps","rpl",
                "clp|mat|inf","ycf","trn","rrn", "OTHER"),
    col = c(psa.color, psb.color, pet.color , atp.color, ndh.color, rbc.color,
            rpo.color, rsp.color, rpl.color, clp_mat_inf.color, ycf.color,
            trn.color,rrn.color, other_gene.color),
    label = c("photosystem I", "photosystem II", "cytochrome b/f complex",
              "ATP synthesis", "NADH dehydrogenase", "RubisCO larg subunit",
              "RNA polymerase", "small ribosomal protein",
              "large ribosomal protein", "clpP, matK, infA",
              "hypothetical reading frame", "transfer RNA", "ribosomal RNA",
              "other"),
    stringsAsFactors = FALSE
  )
}

mitGeneColor <- function(nad.color = "#2A6332", sdh.color = "#4C8805",
                         cob.color = "#7F992C", cox.color = "#FEEE50",
                         atp.color = "#9FBB3D", ccmF.color = "#4D9E3F",
                         mmt.color = "#AE2D29", rps.color = "#D6AD7C",
                         rpl.color = "#9C7A4B", mat.color = "#D9662D",
                         orf.color = "#71B8A9", trn.color = "#172C7F",
                         rrn.color = "#D1382A", other_gene.color = "#7D7D7D"){
  color.table <- data.frame(
    acronym = c("nad|nd","sdh","cob","cox|cytb","atp","ccmF","mtt","rps","rpl",
                "mat","orf","trn","rrn", "OTHER"),
    col = c(nad.color, sdh.color, cob.color, cox.color, atp.color, ccmF.color,
            mmt.color, rps.color, rpl.color, mat.color, orf.color, trn.color,
            rrn.color, other_gene.color),
    label = c("complex I (NADH dehydrogenase)", "complex II (succinate dehydrogenase)",
              "complex III (ubichinol cytochrome reductase)",
              "complex IV (cytochrome c oxidase)", "ATP synthase",
              "cytochrome c biogenesis",
              "transport membrane protein", "ribosomal proteins (SSU)",
              "ribosomal proteins (LSU)", "maturases",
              "ORFs", "transfer RNA", "ribosomal RNA",
              "other"),
    stringsAsFactors = FALSE
  )
}
