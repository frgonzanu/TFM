contraste_pb <- function(pb,
                         target,
                         ref = "Tumor Primario",
                         lfc = log2(1.2),
                         robust = FALSE) {
  
  library(edgeR)
  library(limma)
  
  # 1. Subconjunto del contraste concreto
  pb_sub <- pb[, pb$region_long %in% c(ref, target)]
  
  # 2. Rehacer factores después de filtrar
  pb_sub$region_long <- droplevels(factor(pb_sub$region_long))
  pb_sub$region_long <- relevel(pb_sub$region_long, ref = ref)
  pb_sub$patient <- droplevels(factor(pb_sub$patient))
  
  # 4. Diseño estadístico
  design <- model.matrix(
    ~ patient + region_long,
    data = colData(pb_sub)
    # data = as.data.frame(colData(pb_sub))
  )
  
  # 5. Coeficiente del contraste
  coef_name <- paste0("region_long", target)
  
  
  # 6. Objeto edgeR
  y <- DGEList(
    counts = counts(pb_sub)
  )
  
  # 7. Filtrado de genes poco expresados
  keep <- filterByExpr(y, 
                       design = design,
                       min.count = 10)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # 8. Normalización y estimación de dispersión
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design = design)
  
  # 9. Ajuste quasi-likelihood
  
  fit <- glmQLFit(
    y,
    design = design,
    robust = robust
  )
  
  # 10. Test TREAT frente a umbral mínimo de logFC
  res <- glmTreat(
    fit,
    coef = coef_name,
    lfc = lfc
  )
  
  # 11. Tabla de resultados
  tab <- topTags(res, 
                 n = Inf)$table
  
  tab$gene_id <- rownames(tab)
  
  tab <- tab[, c(
    "gene_id",
    base::setdiff(colnames(tab), 
            "gene_id")
  )]

  
  return(tab)
}