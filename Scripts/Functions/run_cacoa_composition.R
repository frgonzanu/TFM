run_cacoa_composition <- function(
    sce,
    ref,
    target,
    celltype_col = "cell_type",
    region_col = "region_long",
    patient_col = "patient",
    assay_name = "logcounts",
    n.cores = 4,
    estimate_loadings = TRUE
) {
  
  regions_use <- c(ref, target)
  sce_sub <- sce[, colData(sce)[[region_col]] %in% regions_use]
  
  keep <- !is.na(colData(sce_sub)[[celltype_col]]) &
    !is.na(colData(sce_sub)[[region_col]]) &
    !is.na(colData(sce_sub)[[patient_col]])
  
  sce_sub <- sce_sub[, keep]
  
  sce_sub$sample_id <- paste(
    colData(sce_sub)[[patient_col]],
    colData(sce_sub)[[region_col]],
    sep = "_"
  )
  
  cm <- assay(sce_sub, assay_name)
  
  cell.groups <- colData(sce_sub)[[celltype_col]]
  names(cell.groups) <- colnames(sce_sub)
  
  sample.per.cell <- sce_sub$sample_id
  names(sample.per.cell) <- colnames(sce_sub)
  
  sample.meta <- unique(
    as.data.frame(colData(sce_sub))[, c("sample_id", region_col)]
  )
  
  sample.groups <- sample.meta[[region_col]]
  names(sample.groups) <- sample.meta$sample_id
  
  cao <- cacoa::Cacoa$new(
    cm,
    sample.groups = sample.groups,
    cell.groups = cell.groups,
    sample.per.cell = sample.per.cell,
    ref.level = ref,
    target.level = target,
    n.cores = n.cores,
    verbose = TRUE
  )
  if(estimate_loadings){
    cao$estimateCellLoadings()
  } 
  
  return(cao)
}
