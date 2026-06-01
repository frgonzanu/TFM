extract_selectVar <- function(model, comp) {
  
  sel <- selectVar(model, comp = comp)
  
  val <- sel$value
  
  loading <- as.numeric(val$value.var)
  
  tibble(
    gene = sel$name,
    component = paste0("comp", comp),
    loading = loading
  ) %>%
    arrange(desc(abs(loading)))
}