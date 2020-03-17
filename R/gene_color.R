geneColor <- function(psa.color = "#2A6332", psb.color = "#4C8805",
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
