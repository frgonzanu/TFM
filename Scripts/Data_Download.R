library(future.apply)
plan(multisession, workers = 6)

safe_download <- function(url, destfile) {
  if (file.exists(destfile)) {
    message("Ya existe: ", destfile, " — se omite.")
    return(destfile)
  }
  download.file(url, destfile, mode = "wb")
  return(destfile)
}

future_mapply(
  FUN = safe_download,
  url = sub_data$supplementary_file_1,
  destfile = file.path("GSE97693/TPMs", files$name)
)



library(R.utils)

input_dir  <- "../Datos/GSE97693/TPMs"
output_dir <- "GSE97693/TPMs/TXTs"

# Crear la subcarpeta si no existe
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

files <- list.files(input_dir, pattern = "\\.gz$", full.names = TRUE)

for (f in files) {
  # Nombre del archivo sin .gz
  out_file <- file.path(output_dir, sub("\\.gz$", "", basename(f)))
  
  gunzip(f, destname = out_file, remove = FALSE)
}
